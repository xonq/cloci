import os
import sys
import pickle
import subprocess
import multiprocessing as mp
from tqdm import tqdm
from itertools import chain
from collections import defaultdict
from mycotools.acc2fa import dbmain as acc2fa
from mycotools.lib.biotools import gff2list
from mycotools.lib.kontools import multisub, tardir, collect_files, checkdir
from orthocluster.orthocluster.lib.input_parsing import compileCDS2, hg_fa_mngr


def retroactive_grab_hgx_genes(
    gff_path, ome_loc, gene2hg, clusplusminus
    ):                  
                    
    gff_list = gff2list(gff_path)
    clus_hgs = defaultdict(list)
    cds_dict, cds_dict2 = compileCDS2(
        gff_list, os.path.basename(gff_path).replace('.gff3','')
        )               
                            
    # the goal here is to take each HGx and compile and OG by OG list of each
    # gene that is associated with it - I think I am going to need some sort of
    # while loop to accomodate searches outside of the plusminus range.
    # Essentially this will operate by first checking if the entirety of the
    # HGx is within the plusminus range. The next step is that, if not, go
    # ahead and use a while loop to simultaneously parse up and down beyond the
    # plusminus to grab the rest of the HGx and stop when acquired. This can be
    # accomplished by taking the set of HGx and one by one removing the hits
    # found. The while loop would then operate on what is contained within the
    # set. 
    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]):
            if seq0 in ome_loc: # if the og is part of a significant seed 
            # locus
                for hgx, clanI in ome_loc[seq0]:
                    hgs = set(hgx)
                    clus_hgs[hgx].append((clanI, defaultdict(list),))
                    start, end = None, None
                    if i0 < clusplusminus: # is i0 - clusplusminus < 0 ?
                        locus = cds_dict[scaf][:i0+clusplusminus+1] # then gather
                        # all the beginning
                    else:
                        locus = cds_dict[scaf][i0-clusplusminus:i0+clusplusminus+1]
                        # instead get the +/- and adjust for python
                    for i1, seq1 in enumerate(locus): # for each index and sequence
                    # in the locus
                        try:
                            og = gene2hg[seq1]
                        except KeyError:
                            continue
                        if og in hgs:
                            clus_hgs[hgx][-1][1][og].append(seq1)
                        # I'm envisioning a scenario where a P450 from outside
                        # the cluster joins in, even though a P450 OG was an
                        # original part of the HGx; the newly joining P450
                        # would not hit another organism with the cluster and
                        # thus impact the score. it is an assumption that
                        # this is minor.
    return clus_hgs
    # clus_hgs = {hgx: [{og: []}]} list is per locus

def run_make_dmnddb(db, diamond, hgx_dir, og, genes):
    fa_dict = acc2fa(db, genes)
    makeDBcmd = subprocess.Popen([
        diamond, 'makedb', '--db',
        hgx_dir + str(og) + '.dmnd'
        ], stdin = subprocess.PIPE, stdout = subprocess.DEVNULL,
        stderr = subprocess.DEVNULL
        )
    makeDBcmd.communicate(input=dict2fa(fa_dict).encode())[0]
    makeDBcmd.stdin.close()
    makeDBcmd.wait()


def hgx2omes2gbc_calc(
    hgx, omesI, omes, hgxDict, hgx_dir
    ):
    # hgxDict = {og: set(gene1, gene2)}

    # res = {og: {}}
    res = {}
    for og, qs in hgxDict.items():
        hgx_genes, hgx_omes = set(qs), set([q[:q.find('_')] for q in qs])
        hgx_gene_len = len(hgx_genes)
        if hgx_gene_len == 1: # no homologs
            continue
        res[og] = {}
        with open(hgx_dir + str(og) + '.out', 'r') as raw:
        # open the blast results
            geneInfo = defaultdict(list)
            for line in raw:
                d = line.rstrip().split('\t')
                if len(d) == 4:
                    q,s,i,p = d
                elif len(d) == 5:
                    q,s,e,i,p = d
                if q in hgx_genes:
                    geneInfo[q].append((s, float(i), float(p)))
                    # dict(geneInfo) = {query: [(sbj, evalue, id, pos)]}
        geneInfo = {
            k: sorted(v, key = lambda x: x[1], reverse = True) \
            for k,v in geneInfo.items()
            } # sort each query by percent ID (should maybe have a coverage filter)

        # what about paralogs within the same group?
        for gene in qs: # for each gene
            ome = gene[:gene.find('_')] # identify the ome
            if ome not in res[og]:
                res[og][ome] = [0, {}]

#            if sbjct_ome in omes: # SHARED OMES PASS, a conservative approach
            # to identify if the subject is in the family of omes
#            if geneInfo[gene][0][0] in geneInfo: # ONLY SHARED GENES PASS
            hits = {gene}
            # while all cluster homologs aren't accounted for and there remain hits
            if gene not in geneInfo:
                continue
            while geneInfo[gene]:
                sbj_gene = geneInfo[gene][0][0]
                sbj_ome = sbj_gene[:sbj_gene.find('_')]
 #               if sbj_ome == ome: # same ome, could be recent paralog
                if sbj_gene in hgx_genes and sbj_ome != ome:
                    # grab the hit positives and identity
                    res[og][ome][1][sbj_gene] = tuple(geneInfo[gene][0][1:])
                    hits.add(sbj_gene)
                    if not hgx_genes.difference(hits):
                        break
#                elif sbj_ome in hgx_omes: # ome in subjects, could be recent paralog
                else:
                    res[og][ome][0] += 1 # add one to the failed count
                try:
                    geneInfo[gene].pop(0)
                except IndexError:
                    break

            # if there are missing hits
#            print(hgx_gene_len -  len(res[og][ome][1]) + 1)
            if hgx_gene_len - len(res[og][ome][1]) - 1 > 0:
                res[og][ome][0] = 0
            else:
                res[og][ome][0] = hgx_gene_len / (hgx_gene_len + res[og][ome][0])

    # populate a binary response dictionary for each ome and its genes that are
    # in the shared HGx; 0 = the OG's best blast hit in this ome is not another
    # ome code that shares the HGx, 1 = the OG's best blast hit is another ome
    # code that share the HGx
    omeScores = defaultdict(dict)
    ids_y_pos = defaultdict(dict)
    for og in res:
        for ome, data in res[og].items():
            d = data[1]
            omeScores[ome][og] = data[0]
            ids, pos = [x[0] for x in d.values()], [x[1] for x in d.values()]
            if ids:
                ids_y_pos[og][ome] = (min(ids), min(pos))
            else:
                ids, pos = 0, 0
    try:
        # average score of each og (perhaps needs to be weighted by presence)
        # sum of percent top hits that're in GCF adjusted number of ogs appear
        # in) averaged over number of omes
        gbcScore = 0
        ome_av = [sum(v.values())/len(v) for ome, v in omeScores.items()]
        gbcScore = sum(ome_av)/len(omeScores)
    except ZeroDivisionError:
        print(f'\t{hgx} {ome_av}')
        gbcScore = 0

    # average minimum identity / minimum positive; for each gene homolog group
    gcf_min_id, gcf_min_pos, total = 0, 0, 0
    for og, ome_dict in ids_y_pos.items():
        gcf_min_id += sum([x[0] for x in ome_dict.values()])
        gcf_min_pos += sum([x[1] for x in ome_dict.values()])
        total += len(ome_dict)
    try:
        gcf_min_id /= total
        gcf_min_pos /= total
    except ZeroDivisionError:
        gcf_min_id, gcf_min_pos = 0, 0

    return hgx, omesI, gbcScore, gcf_min_id, gcf_min_pos


def gbc_mngr(
    hgs, omes, hgx_dir, hgx2loc,
    db, gene2hg, clusplusminus, hg2gene, modules,
    moduleOmes, moduleHGxs, ome2i, cpus = 1
    ):

    i2ome = {v: k for k, v in ome2i.items()}
    if not os.path.isfile(hgx_dir + 'clusOGs.pickle'):
        hgxGene_cmds = []
        ome_locs = {ome: defaultdict(list) for ome in omes}

        for i, hgx2omes in enumerate(modules):
            modHGx = moduleHGxs[i]
            for hgx, omes in hgx2omes.items():
                omes_set = set(omes)
                for seq in hgx2loc[hgx]:
                    ome = seq[:seq.find('_')]
                    if ome2i[ome] in omes_set:
    #                    if seq not in ome_locs[ome]:
     #                       ome_locs[ome][seq] = []
                        ome_locs[ome][seq].append((modHGx, i,))

        for ome, ome_loc in ome_locs.items():
            gff = db[ome]['gff3']
            hgxGene_cmds.append([gff, ome_loc, gene2hg, clusplusminus])

        print('\tAssimilating GCF loci', flush = True)
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            clus_hgs_prep = pool.starmap(retroactive_grab_hgx_genes, 
                                         tqdm(hgxGene_cmds, total = len(hgxGene_cmds)))
        with open(hgx_dir + 'clusOGs.pickle', 'wb') as out:
            pickle.dump(clus_hgs_prep, out)
    else:
        with open(hgx_dir + 'clusOGs.pickle', 'rb') as raw:
            clus_hgs_prep = pickle.load(raw)

    clus_hgs = {i: (modHGx, defaultdict(list),) \
                for i, modHGx in enumerate(moduleHGxs)}

    for res in clus_hgs_prep:
        # clus_hgs = {hgx: [{og: []}]} list is per locus
        for hgx, loci in res.items():
            for clanI, locus in loci:
                for og in locus:
                    clus_hgs[clanI][1][og].extend(locus[og])
    clus_hgs = {
        clanI: [d[0], {og: set(v) for og, v in d[1].items()}] \
        for clanI, d in clus_hgs.items()
        } # make sets from it
    # {clanI: [hgx, {og: set(seqs)}}

    print('\tCalculating gene blast congruence (GBC) scores', flush = True)
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        gbc_res = pool.starmap(
            hgx2omes2gbc_calc,
            ([moduleHGxs[clanI], moduleOmes[clanI],
            [i2ome[i] for i in moduleOmes[clanI]],
            d[1], hgx_dir] \
             for clanI, d in clus_hgs.items())
            )

    hgx2omes2gbc = defaultdict(dict)
    hgx2omes2id = defaultdict(dict)
    hgx2omes2pos = defaultdict(dict)
    for hgx, omes, score, id_, pos in gbc_res:
        hgx2omes2gbc[hgx][omes] = score
        hgx2omes2id[hgx][omes] = id_
        hgx2omes2pos[hgx][omes] = pos

    return hgx2omes2gbc, hgx2omes2id, hgx2omes2pos


def prep_blast_cmds(db, hgs, hg_dir, hgx_dir, minid = 30, diamond = 'diamond'):

    alnd_hgs = set([int(os.path.basename(x[:-4])) \
                    for x in collect_files(hgx_dir, 'out')])
    dmnd_dbs = set([int(os.path.basename(x[:-5])) \
                    for x in collect_files(hgx_dir, 'dmnd')])
    missing_alns = set(hgs).difference(alnd_hgs)
    missing_dmnds = set(missing_alns).difference(dmnd_dbs)

    cmds1 = [(diamond, 'makedb', '--db', f'{hgx_dir}{hg}.dmnd',
              '--in', f'{hg_dir}{hg}.faa', '--threads', '2') \
              for hg in missing_dmnds]
    cmds2 = [((diamond, 'blastp', '--query', f'{hg_dir}{hg}.faa',
              '--db', f'{hgx_dir}{hg}.dmnd', '-o', f'{hgx_dir}{hg}.out.tmp',
              '--outfmt', '6', 'qseqid', 'sseqid', 'evalue', 'pident',
              'ppos', '--threads', '2', '--id', str(minid),
              '--no-self-hits', '&&'), ('mv', f'{hgx_dir}{hg}.out.tmp',
              f'{hgx_dir}{hg}.out')) for hg in missing_alns]

    return cmds1, cmds2


def run_blast(hgs, db, hg_dir, hgx_dir, 
              diamond = 'diamond', printexit = False, cpus = 1):
#    hgs = sorted(set(chain(*list(hgx2loc.keys()))))
    db_cmds, dmnd_cmds = prep_blast_cmds(db, hgs, hg_dir, hgx_dir, diamond = 'diamond')
    if printexit and dmnd_cmds:
        print('\tPreparing diamond commands and exiting', flush = True)
        with open(hgx_dir + '../../makedb.sh', 'w') as out:
            out.write('\n'.join([' '.join(x) for x in db_cmds]))
        with open(hgx_dir + '../../srch.sh', 'w') as out:
            out.write('\n'.join([' '.join(x) for x in dmnd_cmds]))
        print('\nGBC commands outputted to `<OUTPUT>/*.sh` ' \
            + 'Run `makedb.sh` first', flush = True)
        sys.exit(0)
    elif dmnd_cmds:
        print(f'\tBuilding {len(db_cmds)} diamond dbs', flush = True)
        multisub(db_cmds, verbose = 2, processes = cpus)
        print(f'\tAligning {len(dmnd_cmds)} HGs', flush = True)
        multisub(dmnd_cmds, verbose = 2, processes = cpus,
                 injectable = True)


def gbc_main(
    hgx2loc, wrk_dir, ome2i, hg_dir, hgx_dir,
    blastp, db, gene2hg, plusminus, hg2gene,
    old_path = 'hgx2gbc.pickle',
    modules = None, moduleHGxs = None,
    moduleOmes = None, cpus = 1, printexit = False
    ):

    if not os.path.isfile(wrk_dir + old_path):
        if not checkdir(hgx_dir, unzip = True, rm = True):
            os.mkdir(hgx_dir)

        
        hgs = list(chain(*moduleHGxs))
        hg_dir = hg_fa_mngr(wrk_dir, hg_dir, hgs, 
                            db, hg2gene, cpus = cpus,
                            low_mem = True)
        run_blast(hgs, db, hg_dir, 
                  hgx_dir, cpus = cpus,
                  diamond = 'diamond', printexit = printexit)
    
        d2gbc, d2id_, d2pos = gbc_mngr(
            list(hgs), list(ome2i.keys()), hgx_dir, hgx2loc,
            db, gene2hg, plusminus, hg2gene, modules,
            moduleOmes, moduleHGxs, ome2i, cpus = cpus
            )
        hgx_dirTar = mp.Process(target=tardir, args=(hgx_dir, True))
        hgx_dirTar.start() # when to join ...
        with open(wrk_dir + 'hgx2omes2gbc.full.pickle', 'wb') as pickout:
            pickle.dump(d2gbc, pickout)
        with open(wrk_dir + 'hgx2omes2id.full.pickle', 'wb') as pickout:
            pickle.dump(d2id_, pickout)
        with open(wrk_dir + 'hgx2omes2pos.full.pickle', 'wb') as pickout:
            pickle.dump(d2pos, pickout)
    else:
        print('\tLoading previous coevolution results', flush = True)
        with open(wrk_dir + 'hgx2omes2gbc.full.pickle', 'rb') as pickin:
            d2gbc = pickle.load(pickin)
        with open(wrk_dir + 'hgx2omes2id.full.pickle', 'rb') as pickin:
            d2id_ = pickle.load(pickin)
        with open(wrk_dir + 'hgx2omes2pos.full.pickle', 'rb') as pickin:
            d2pos = pickle.load(pickin)

    return d2gbc, d2id_, d2pos
