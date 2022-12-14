import os
import sys
import pickle
import subprocess
import multiprocessing as mp
from collections import defaultdict
from mycotools.lib.kontools import multisub, tardir
from mycotools.lib.biotools import gff2list
from mycotools.acc2fa import dbmain as acc2fa
from orthocluster.orthocluster.lib.input_parsing import compileCDS2


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


def blast_homolog(db, hgs, hg_dir, hgx_dir,
                  diamond, hg2gene, cpus = 1,
                  printexit = False):

    cmds1, cmds2 = [], []
    for og in hgs:
        if not os.path.isfile(hgx_dir + str(og) + '.out'):
            hg_file = hg_dir + str(og) + '.faa'
            out_file = hgx_dir + str(og) + '.out'
            if not os.path.isfile(hg_file): # try OrthoFinder check
                digits = len(str(og))
                zeros = 7 - digits
                hg_file = hg_dir + 'OG' + '0' * zeros + str(og) + '.fa'
#            if len(hg2gene[og]) > 100: # key here is efficiency
        if not os.path.isfile(hgx_dir + str(og) + '.dmnd'):
#            cmds1.append([
 #               db, blastp, hgx_dir, og, hg2gene[og]
  #              ])
            cmds1.append([diamond, 'makedb', '--db',
                          hgx_dir + str(og) + '.dmnd',
                          '--in', hg_file])

            cmds2.append([
                diamond, 'blastp', '--query', hg_file,
                '--db', hgx_dir + str(og) + '.dmnd',
                '-o', out_file, '-outfmt',
                '-6', 'qseqid', 'sseqid', 'evalue', 'pident', 'ppos'
                ])
 #               search_cmd = subprocess.call(('mmseqs easy-search ' + hg_file + ' ' \
  #                          + hg_file + ' ' + out_file + '.tmp ' \
   #                         + hgx_dir + 'tmp/ --threads ' + str(cpus) + ' ' \
    #                        + '--format-output query,target,evalue').split(' '),
     #                       stdout = subprocess.DEVNULL,
      #                      stderr = subprocess.DEVNULL)
       #     else:
        #        search_cmd = subprocess.call('blastp -query {hg_file} -subject {hg_file} \
         #                    -out {out_file}.tmp -num_threads {cpus} -outfmt \
          #                   "6 qseqid sseqid evalue"' %
           #                  {'hg_file': hg_file,
            #                 'out_file': out_file,
             #                'cpus': str(cpus)},
              #               shell = True, stdout = subprocess.DEVNULL,
               #              stderr = subprocess.DEVNULL)
   #         if os.path.isfile(out_file + '.tmp'):
  #              os.rename(out_file + '.tmp', out_file)
 #           else:
#                raise FileNotFoundError('BLASTp failed: ' + str(search_cmd) \
    #                                  + ' ' + out_file)

    if printexit:
        with open(hgx_dir + '../../gbc_makedb.sh', 'w') as out:
            out.write('\n'.join([' '.join(x) for x in cmds1]))
        with open(hgx_dir + '../../gbc_srch.sh', 'w') as out:
            out.write('\n'.join([' '.join(x) for x in cmds2]))
        print('\nGBC commands outputted to `<OUTPUT>/gbc*.sh` \
               Run makedb first', flush = True)
        sys.exit(0)

    print('\tMaking ' + str(len(cmds1)) + ' diamond databases', flush = True)
#    with mp.get_context('fork').Pool(processes = cpus) as pool:
 #       pool.starmap(run_make_dmnddb, cmds1)
    multisub(cmds1, processes = cpus)

    print('\tRunning ' + str(len(cmds2)) + ' diamonds', flush = True)
    multisub(cmds2, processes = cpus)


def gbc_mngr_2( 
    hgs, omes, hg_dir, hgx_dir, diamond, hgx2loc, 
    db, gene2hg, clusplusminus, hg2gene, cpus = 1,
    printexit = False
    ):
                    
    blast_homolog(db, hgs, hg_dir, hgx_dir, diamond, hg2gene, cpus = 1,
                  printexit = printexit) 
    hgxGene_cmds = []
    ome_locs = {ome: defaultdict(list) for ome in omes}
    for hgx in hgx2loc:
        for seq in hgx2loc[hgx]:
            ome = seq[:seq.find('_')]
            ome_locs[ome][seq].append((hgx, None,))
    
    for ome in ome_locs:
        gff = db[ome]['gff3']
        hgxGene_cmds.append([gff, ome_locs[ome], gene2hg, clusplusminus])
  
    print('\tAssimilating loci with significant HGxs', flush = True)
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        clus_hgs_prep = pool.starmap(retroactive_grab_hgx_genes, hgxGene_cmds)
#    {hgx:[{og:[seq]}]}
        
    clus_hgs = {} 
    for res in clus_hgs_prep:
        for hgx, loci in res.items():
            if hgx not in clus_hgs:
                clus_hgs[hgx] = {x: [] for x in hgx}
            for null, locus in loci:
                for og in locus:
                    clus_hgs[hgx][og].extend(locus[og])
    clus_hgs = {    
        hgx: {og: set(v) for og, v in hgs.items()} for hgx, hgs in clus_hgs.items()
        } # make sets from it
    # {hgx: {og: set(seqs)}}
            
    print('\tCalculating gene blast congruence (GBC) scores', flush = True)
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        gbc_res = pool.starmap(
            calc_gbc, [[hgx, clus_hgs[hgx], hgx_dir] for hgx in clus_hgs]
            )
                
    gcb_scores = dict(gbc_res)
     
    return gcb_scores 


def hgx2omes2gbc_calc(
    hgx, omesI, omes, hgxDict, hgx_dir
    ):
    # hgxDict = {og: set(gene1, gene2)}

    # res = {og: {}}
    res = {}
    for og, qs in hgxDict.items():
        hgx_genes, hgx_omes = set(qs), set([q[:q.find('_')] for q in qs])
        res[og] = {}
        with open(hgx_dir + str(og) + '.out', 'r') as raw:
        # open the blast results
            geneInfo = defaultdict(list)
            for line in raw:
                q,s,e,i,p = line.rstrip().split('\t')
                if q in qs:
                    geneInfo[q].append((s, float(e), float(i), float(p)))
                    # dict(geneInfo) = {query: [(sbj, evalue, id, pos)]}
        geneInfo = {
            k: sorted(v, key = lambda x: x[2], reverse = True) \
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
            while hgx_genes.difference(hits) and geneInfo[gene]:
                sbj_gene = geneInfo[gene][0][0]
                sbj_ome = sbj_gene[:sbj_gene.find('_')]
                if sbj_ome == ome: # same ome, could be recent paralog
                    geneInfo[gene].pop(0)
                elif sbj_gene in hgx_genes:
                    # grab the hit positives and identity
                    res[og][ome][1][sbj_gene] = tuple(geneInfo[gene][0][2:])
                    hits.add(sbj_gene)
                elif sbj_ome in hgx_omes: # ome in subjects, could be recent paralog
                    continue
                else:
                    res[og][ome][0] += 1 # add one to the failed count
            # conservative penality
                # penalize the cluster family by 1 for each hit that's missing
#            res[og][ome][0] += len(hgx_genes) - len(res[og][ome]) - 1
 #           res[og][ome][0] = len(hgx_genes) / res[og][ome][0]
            # normalize to 0.5-1 because 0.5 is lowest possible completely penalized
  #          if res[og][ome][0] < 0.5:
   #             res[og][ome][0] = 0
    #        else:
     #           res[og][ome][0] = (res[og][ome][0] - 0.5) / 0.5
            if len(hgx_genes) - len(res[og][ome]) + 1:
                res[og][ome][0] = 0
            else:
                res[og][ome][0] = len(hgx_genes) / res[og][ome][0]


    # populate a binary response dictionary for each ome and its genes that are
    # in the shared HGx; 0 = the OG's best blast hit in this ome is not another
    # ome code that shares the HGx, 1 = the OG's best blast hit is another ome
    # code that share the HGx
    omeScores = defaultdict(dict)
    ids_y_pos = defaultdict(dict)
    for og in res:
        for ome, d in res[og].items():
            omeScores[ome][og] = d[0]
            ids, pos = [x[0] for x in d.values()], [x[1] for x in d.values()]
            ids_y_pos[og][ome] = (min(ids), min(pos))
    try:
        # average score of each og (perhaps needs to be weighted by presence)
        # sum of percent top hits that're in GCF adjusted number of ogs appear
        # in) averaged over number of omes
        gbcScore = \
            sum([sum(omeScores[ome].values()) \
            / len(omeScores[ome]) for ome in omeScores]) \
            / len(omeScores)
    except ZeroDivisionError:
        gbcScore = 0

    # average minimum identity / minimum positive; for each gene homolog group
    gcf_min_id, gcf_min_pos, total = 0, 0, 0
    for og, ome_dict in ids_y_pos.items():
        gcf_min_id += sum([x[0] for x in ome_dict.values()]) / len(ome_dict)
        gcf_min_pos += sum([x[1] for x in ome_dict.values()]) / len(ome_dict)
        total += len(ome_dict)
    gcf_min_id /= total
    gcf_min_pos /= total

    return hgx, omesI, gbcScore, gcf_min_id, gcf_min_pos


def gbc_mngr_3(
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
            clus_hgs_prep = pool.starmap(retroactive_grab_hgx_genes, hgxGene_cmds)
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


def calc_gbc(
    hgx, hgxDict, hgx_dir
    ):
    # hgxDict = {og: set(gene1, gene2)}

    res= {}
    for og in hgxDict:
        res[og] = {}
        with open(hgx_dir + str(og) + '.out', 'r') as raw:
            geneInfo = defaultdict(list)
            for line in raw:
                # query,subject,...,bit
                d = line.rstrip().split('\t')
                q = d[0]
                s = d[1]
                evalue = float(d[-1])
                if q in hgxDict[og]:
                    geneInfo[q].append((s, evalue))
        geneInfo = {
            k: sorted(v, key = lambda x: x[1]) for k,v in
            geneInfo.items()
            } # would need to reverse if using bitscore

        for gene in list(hgxDict[og]):
            ome = gene[:gene.find('_')]
            if ome not in res[og]:
                res[og][ome] = False
            if res[og][ome]:
                continue
            try:
                while geneInfo[gene][0][0][:geneInfo[gene][0][0].find('_')] == ome:
                    geneInfo[gene].pop(0)
            except IndexError:
                continue
            except KeyError: # the gene is not in the blast results (too short?)
                continue # would be nice to quantitate the percent of this
            if geneInfo[gene][0][0] in geneInfo:
                res[og][ome] = geneInfo[gene][0][0] # this isn't a reciprocal
                # analysis as it currently stands, its just a check for
                # correlation

    omeScores = defaultdict(dict)
    for og in res:
        for ome in res[og]:
            omeScores[ome][og] = 0
            if res[og][ome]:
                omeScores[ome][og] = 1

    try:
        gbcScore = \
            sum([sum(omeScores[ome].values()) / len(omeScores[ome]) for ome in omeScores]) \
            / len(omeScores) # overall average of average binary positive per ome 
    except ZeroDivisionError:
        gbcScore = 0

    return hgx, gbcScore


def gbc_main_0(
    hgx2loc, wrk_dir, ome2i, hg_dir, hgx_dir,
    blastp, db, gene2hg, plusminus, hg2gene, 
    old_path = 'hgx2gbc.pickle',
    cpus = 1, printexit = False
    ):
    
    hgs_list = []
    for hgx in hgx2loc:
        hgs_list.extend(list(hgx))
    hgs = set(hgs_list) 
    if not os.path.isfile(wrk_dir + old_path):
        d2gbc = gbc_mngr_2(
            list(hgs), list(ome2i.keys()), hg_dir, hgx_dir, blastp, hgx2loc,
            db, gene2hg, plusminus, hg2gene, cpus = cpus, 
            printexit = printexit
            )
    else:
        print('\tLoading previous coevolution results', flush = True)
        with open(wrk_dir + old_path, 'rb') as pickin: 
            d2gbc = pickle.load(pickin)
    return d2gbc


def gbc_main_1(
    hgx2loc, wrk_dir, ome2i, hg_dir, hgx_dir,
    blastp, db, gene2hg, plusminus, hg2gene,
    old_path = 'hgx2gbc.pickle',
    modules = None, moduleHGxs = None,
    moduleOmes = None, cpus = 1, printexit = False
    ):

    hgs_list = []
    for hgx in hgx2loc:
        hgs_list.extend(list(hgx))
    hgs = set(hgs_list)
    if not os.path.isfile(wrk_dir + old_path):
        d2gbc, d2id_, d2pos = gbc_mngr_3(
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
