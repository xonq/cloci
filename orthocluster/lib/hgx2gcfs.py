import os
import shutil
import pickle
import subprocess
import networkx as nx
import multiprocessing as mp
from datetime import datetime
from itertools import combinations
from scipy.sparse import lil_matrix
from collections import defaultdict, Counter
from mycotools.lib.biotools import fa2dict, dict2fa, gff2list
from mycotools.lib.kontools import write_json, read_json, collect_files
from orthocluster.orthocluster.lib import input_parsing, phylocalcs

def hash_hgx(gff_path, ome, hgx_genes, gene2hg, clusplusminus):
    """sliding window of size clusplusminus, identify sequences that may have
    hgs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus"""
        
    gff_list, gene_dict = gff2list(gff_path), {}
    cds_dict = input_parsing.compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))
                
    hgx_dict = {}   
    for scaf in cds_dict: 
        for i0, seq0 in enumerate(cds_dict[scaf]):
            if seq0 in hgx_genes: # if the og is part of a significant seed
            # locus 
                if i0 < clusplusminus: # is i0 - clusplusminus < 0 ?
                    locus = cds_dict[scaf][:i0+clusplusminus+1] # then gather
                    # all the beginning
                else:
                    locus = cds_dict[scaf][i0-clusplusminus:i0+clusplusminus+1]
                    # instead get the +/- and adjust for python
                og0 = gene2hg[seq0]
                for hgx in hgx_genes[seq0]:
                    if hgx not in hgx_dict:
                        hgx_dict[hgx] = {og: [] for og in hgx}
                    hgx_dict[hgx][og0].append(seq0) # add the sequence to the hgx_dict
                    start, end = None, None
                    for i1, seq1 in enumerate(locus): # for each index and sequence
                    # in the locus
                        try:
                            og1 = gene2hg[seq1]
                        except KeyError: # missing entry
                            continue
                        if og1 in set(hgx) and start is None: # if the og is
                        # in the sig clus and we haven't started
                            start = i1 # this is the start index
                            hgx_dict[hgx][og1].append(seq1)

                        elif og1 in set(hgx): # otherwise if it is in the
                        # sig clus label this the other border unless another is 
                        # found
                            end = i1 + 1
                            hgx_dict[hgx][og1].append(seq1)
                    for gene in locus[start:end]:
                        if gene not in gene_dict:
                            gene_dict[gene] = []
                        gene_dict[gene].append(hgx)

    gene_tup = tuple([
        tuple([gene, tuple(sorted(vals))]) \
        for gene, vals in gene_dict.items()
        ])
    # ((gene, {gene: [HGx...]}))

    hgx_tup = tuple([
        tuple([
            hgx,
            tuple([tuple([og, tuple(sorted(set(seqs)))]) \
                for og, seqs in hgs.items()]) \
            ])
            for hgx, hgs in hgx_dict.items()
        ])

    return ome, gene_tup, hgx_tup


def find_hgx_pairs(gene2hgx, ome, hgx2i, pairsDict, minimum = 2):
    print(ome, flush = True)
    # is this going to merge hgxs that overlap by one gene in minimum loci?
    hgx_pairs_raw = []
    # check for gene overlap in hgxs
    for hgxs in list(gene2hgx.values()):
        # add each possible hgxpair
        hgx_pairs_raw.extend([x for x in combinations(sorted(hgxs), 2)])

    hgx_pairs_count = [x for x, g in Counter(hgx_pairs_raw).items() \
                       if g > minimum]
    null = [pairsDict[hgxpair].append(ome) for hgxpair in hgx_pairs_count]

    return pairsDict

def MCL(adj_path, clusFile, inflation = 1.0, threads = 1):

    subprocess.call([
        'mcl', adj_path, '--abc', '-I', str(inflation),
        '-o', clusFile, '-te', str(threads)
        ], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
        )
#    mclFile = prep_mcl(adj_path + '.tmp', sym = True, output = None) 
    # symmetrical matrix, no verbose, mcxload must be in path
 #   run_mcl(mclFile, clusFile, inflation, cpus = cpus)
    # parallelization could be more efficient overall

    return clusFile


def acquire_clus_gcf_sim_noq(
    i0, i1, index, loc0, loc1, ogL0, ogL1, set0, set1, blast_ids,
    mingcfid = 0.15
    ):
    ogDict = {}
    for i, og in enumerate(ogL0):
        if og not in ogDict:
            ogDict[og] = [[], []]
        ogDict[og][0].append(loc0[i])
    for i, og in enumerate(ogL1):
        if og not in ogDict:
            ogDict[og] = [[], []]
        ogDict[og][1].append(loc1[i])

    intersection = set0.intersection(set1)
    #jaccard = len(intersection)/len(set0.union(set1))
    # currently this is biased against contig edges and will also
    # force fragments with a much smaller subset OG # then the
    # primary HGx into separate gcfs

    # could alternatively weight just by incorporating 0 values
    # for OGs that don't overlap

    overlap_coef = len(intersection)/min([len(set0), len(set1)])
    scores = []
    for og in list(intersection):
        if ogDict[og][0] and ogDict[og][1]:
            scores.append([])
            for gene0 in ogDict[og][0]:
                for gene1 in ogDict[og][1]:
                    try:
                        scores[-1].append(blast_ids[og][gene0][gene1])
                    except KeyError: # missing gene
                        scores[-1].append(0)
        else:
            scores[-1].append(0)
    maxScores = [max(i) for i in scores]
    try:
        total = (sum(maxScores)/len(maxScores)) * overlap_coef # * jaccard
        if total > mingcfid:
            return str(i0 + index) + '\t' + str(i1 + index) + '\t' + str(total) + '\n'
    except ZeroDivisionError: # no overlapping OGs
#        print(set0,set1, flush = True)
        return


def acquire_clus_gcf_sim(
    i0, i1, index, loc0, loc1, ogL0, ogL1, set0, set1, blast_ids, Q
    ):
    ogDict = {}
    for i, og in enumerate(ogL0):
        if og not in ogDict:
            ogDict[og] = [[], []]
        ogDict[og][0].append(loc0[i])
    for i, og in enumerate(ogL1):
        if og not in ogDict:
            ogDict[og] = [[], []]
        ogDict[og][1].append(loc1[i])

    intersection = set0.intersection(set1)
    #jaccard = len(intersection)/len(set0.union(set1))
    # currently this is biased against contig edges and will also
    # force fragments with a much smaller subset OG # then the
    # primary HGx into separate gcfs

    # could alternatively weight just by incorporating 0 values
    # for OGs that don't overlap

    overlap_coef = len(intersection)/min([len(set0), len(set1)])
    scores = []
    for og in list(intersection):
        if ogDict[og][0] and ogDict[og][1]:
            scores.append([])
            for gene0 in ogDict[og][0]:
                for gene1 in ogDict[og][1]:
                    try:
                        scores[-1].append(blast_ids[og][gene0][gene1])
                    except KeyError: # missing gene
                        scores[-1].append(0)
        else:
            scores[-1].append(0)
    maxScores = [max(i) for i in scores]
    try:
        total = (sum(maxScores)/len(maxScores)) * overlap_coef # * jaccard
        if total > 0.0:
            Q.put(str(i0 + index) + '\t' + str(i1 + index) + '\t' + str(total) + '\n')
    except ZeroDivisionError: # no overlapping OGs
#        print(set0,set1, flush = True)
        return


def clan_to_gcf_loci_resume(
    db, clanI, loci, hgLoci, hg_dir, hgx_dir, gcf_dir, index,
    minid = 30, mingcfid = 0.15
    ):

    finished_file = f'{gcf_dir}{clanI}.sadj.tmp'
    # is this clan already completed?
    if os.path.isfile(finished_file):
        return
    # was this clan partially written?
    elif os.path.isfile(f'{finished_file}.w'):
        # identify the last definitively known clan
        with open(f'{finished_file}.w1', 'w') as out:
            with open(f'{finished_file}.w', 'r') as raw:
                data, max_complete = [], 0
                for line in raw:
                    d = line.rstrip().split()
                    try:
                        if d[0] != data[-1][0]:
                            out.write('\n'.([' '.join(x) for x in data]) \
                                    + '\n')
                            out.flush()
                            max_complete = int(data[-1][0]) - index
                            data = []
                    except IndexError: # first line
                        pass
                    data.append(d)
        shutil.move(f'{finished_file}.w1', f'{finished_file}.w')
    # nothing completed
    else:
        max_complete = 0

    blast_hash = defaultdict(list)
    # only examine loci that have not been completed
    for i in range(max_complete, len(loci)):
        locus = loci[i]
        hgs = hgLoci[i]
        for i1, hg in enumerate(hgs):
            if hg is not None:
                blast_hash[hg].append(locus[i1])

    blast_ids = defaultdict(dict)
    for hg, genes in blast_hash.items():
        if len(genes) > 1:
            gene_set = set(genes)
            algn_base = f'{hgx_dir}{hg}.out'
            gene_len = len(genes)
            try:
                with open(algn_base, 'r') as raw:
                    for line in raw:
                        d = line.rstrip().split()
                        q, s, pident = d[0], d[1], d[-2]
                        if {q, s}.issubset(gene_set):
                            if q not in blast_ids[hg]:
                                blast_ids[hg][q] = {}
                            blast_ids[hg][q][s] = float(pident)/100 # adjust diamond to decimal
                        if len(blast_ids[hg]) == gene_len:
                            break
            except FileNotFoundError:
                 blast_ids = run_gcf_blast(blast_ids, hg_dir, hg, genes, 
                              gcf_dir, clanI, minid = minid,
                              diamond = 'diamond')  

    if blast_ids:
        # begin appending to temporary writing output
        with open(f'{finished_file}.w', 'a') as out:
            # for all unique, unexamined pairwise loci comparisons
            for i0, i1 in combinations(range(max_complete, len(loci)), 2):
                loc0, loc1 = loci[i0], loci[i1]
                hgL0, hgL1 = hgLoci[i0], hgLoci[i1]
                sHGl0 = set([x for x in hgL0 if x is not None])
                sHGl1 = set([x for x in hgL1 if x is not None])
                # if there is an intersection between HGs, calculate the gcf sim
                if not sHGl0.isdisjoint(sHGl1):
                    data = acquire_clus_gcf_sim_noq(i0, i1, index, loc0, 
                                                    loc1, hgL0, hgL1, 
                                                    sHGl0, sHGl1, blast_ids,
                                                    mingcfid = mingcfid)
                    if data:
                        out.write(data)
                        out.flush()
        # the file is now complete
        shutil.move(f'{gcf_dir}{clanI}.sadj.tmp.w', f'{gcf_dir}{clanI}.sadj.tmp')
    

def run_gcf_blast(blast_ids, hg_dir, hg, genes, 
                  gcf_dir, clanI, minid = 30,
                  diamond = 'diamond'):
    fileBase = gcf_dir + str(clanI) + '.' + str(hg)
    if not os.path.isfile(fileBase + '.out'):
        f_fa_dict = fa2dict(f'{hg_dir}{hg}.faa')
        fa_dict = {g: f_fa_dict[g] for g in genes}
        with open(fileBase + '.faa', 'w') as out:
            out.write(dict2fa(fa_dict))
        makeDBcmd = subprocess.call([
            diamond, 'makedb', '--in', fileBase + '.faa', '--db',
            fileBase + '.dmnd', '--threads', str(2)
            ], stdout = subprocess.DEVNULL,
            stderr = subprocess.DEVNULL
            )
        dmndBlast = subprocess.call([diamond, 'blastp', '--query', 
            fileBase + '.faa', '--db', fileBase, '--threads', str(2), 
            '--id', str(minid), '--no-self-hits', '-o', 
            fileBase + '.out.tmp', '--outfmt', '6', 'qseqid', 'sseqid', 
            'pident', 'ppos'
            ], stdin = subprocess.PIPE, stdout = subprocess.DEVNULL, 
             stderr = subprocess.DEVNULL)
        shutil.move(fileBase + '.out.tmp', fileBase + '.out')
    
    with open(fileBase + '.out', 'r') as raw:
        for line in raw:
            d = line.rstrip().split()
            q, s, pident = d[0], d[1], d[-2]
            if q not in blast_ids[hg]:
                blast_ids[hg][q] = {}
            blast_ids[hg][q][s] = float(pident)/100 # adjust diamond to decimal

    return blast_ids


def clan_to_gcf_loci_sensitive(
    db, clanI, loci, hgLoci, hg_dir, gcf_dir, Q, index,
    diamond = 'diamond', minid = 30
    ):

    blast_hash = defaultdict(list)
    for i, locus in enumerate(loci):
        hgs = hgLoci[i]
        for i1, hg in enumerate(hgs):
            if hg is not None:
                blast_hash[hg].append(locus[i1])

    blast_ids = defaultdict(dict)

    for hg, genes in blast_hash.items():
        if len(genes) > 1:
             blast_ids = run_gcf_blast(blast_ids, hg_dir, hg, genes, 
                              gcf_dir, clanI, minid = 30,
                              diamond = 'diamond')  
            
    blast_ids = dict(blast_ids)
    if blast_ids:
        for i0, i1 in combinations(range(len(loci)), 2):
            loc0, loc1 = loci[i0], loci[i1]
            hgL0, hgL1 = hgLoci[i0], hgLoci[i1]
            sHGl0 = set([x for x in hgL0 if x is not None])
            sHGl1 = set([x for x in hgL1 if x is not None])
            if not sHGl0.isdisjoint(sHGl1): # if there is an intersection 
#                print('\t', i0, i1, flush = True)
                acquire_clus_gcf_sim(i0, i1, index, loc0, loc1, hgL0, hgL1, sHGl0,
                                 sHGl1, blast_ids, Q)



def hash_clan_loci(ome, gff_path, ome_sig_clus, gene2hg, clusplusminus):
    """sliding window of size clusplusminus, identify sequences that may have
    hgs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus"""

    gff_list, protoclus, clus_out = gff2list(gff_path), {}, []
    cds_dict = input_parsing.compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))

    preclanLoci = defaultdict(list)
    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]):
            if seq0 in ome_sig_clus: # if the og is part of a significant seed
            # locus
                for sig_clus in ome_sig_clus[seq0]:
                    clan = sig_clus[1]
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
                        if og in sig_clus[0] and start is None: # if the og is
                        # in the sig clus and we haven't started
                            start = i1 # this is the start index
                        elif og in sig_clus[0]: # otherwise if it is in the
                        # sig clus label this the other border unless another is
                        # found
                            end = i1 + 1
                    hgx = tuple(sorted(sig_clus[0]))
                    preclanLoci[clan].append([
                         locus[start:end], hgx, [hgx]
                         ])
                    # clus_out = [(og0, og1, ... ogn), [prot0, prot1, ... protn])]

    outclanLoci, outclanHGx = defaultdict(list), defaultdict(list)
    for clan, loci in preclanLoci.items():
        while loci: # exhaustively merge overlapping loci
            loc0, hgxs0, allHGxs = set(loci[0][0]), loci[0][1], loci[0][2]
            locIntersect = None
            for i1, loc1d in enumerate(loci[1:]):
                loc1, hgxs1, allHGxs1 = set(loc1d[0]), loc1d[1], loc1d[2]
                locIntersect = loc0.intersection(loc1)
                if locIntersect: # check for overlap
                    newHGx = list(hgxs0)
                    newHGx.extend(list(hgxs1))
                    allHGxs.extend(allHGxs1)
                    loci.append([list(loc1.union(loc0)), newHGx, allHGxs])
                    break
            if locIntersect: # if there was overlap, delete the overlappers
                del loci[0]
                del loci[i1 + 1]
            else: # no more overlap for this locus, add to the final output
                outclanLoci[clan].append(sorted(loc0))
                outclanHGx[clan].append(tuple(sorted(set(allHGxs))))
#                outclanHGx[clan].append(tuple(sorted(set(hgxs0))))
                del loci[0]

# blast each locus OG against itself
    outogLoci = {}
    for clan, loci in outclanLoci.items():
        outogLoci[clan] = []
        for locus in loci:
            outogLoci[clan].append([])
            for gene in locus:
                try:
                    outogLoci[clan][-1].append(gene2hg[gene])
                except KeyError:
                    outogLoci[clan][-1].append(None)

    return ome, outclanLoci, outclanHGx, outogLoci

def write_adj_matrix_noq(out_file, gcf_dir):
    with open(out_file, 'w') as out:
        init_files = collect_files(gcf_dir, 'sadj.tmp.r')
        for f in init_files:
            with open(f, 'r') as raw:
                for line in raw:
                    out.write(line.rstrip() + '\n')
            out.flush()

        files = collect_files(gcf_dir, 'sadj.tmp')
        f_set = set(files)
        while f'{gcf_dir}done.sadj.tmp' not in f_set:
            for f in files:
                print(f, flush = True)
                with open(f, 'r') as raw:
                    for line in raw:
                        out.write(line.rstrip() + '\n')
                shutil.move(f, f'{f}.r')
                out.flush()
            files = collect_files(gcf_dir, 'sadj.tmp')
        for f in files:
            if f == f'{gcf_dir}done.sadj.tmp':
                os.remove(f'{gcf_dir}done.sadj.tmp')
            else:
                with open(f, 'r') as raw:
                    for line in raw:
                        out.write(line.rstrip() + '\n')
                shutil.move(f, f'{f}.r')

    for f in collect_files(gcf_dir, 'sadj.tmp.r'):
        os.remove(f)

def write_adj_matrix(Q, out_file):
    with open(out_file, 'w') as out: # might be nice to compress this/write to MCL binary matrix
        x = True
        while x:
            x = Q.get()
            if x:
                out.write(x)
#                out.flush() # shouldn't be a bottleneck


def classify_gcfs(
    hgx2loc, db, gene2hg, i2hgx, hgx2i,
    phylo, bordScores, ome2i, hgx2omes, hg_dir, hgx_dir,
    wrk_dir, ome2partition, omes2dist = {}, clusplusminus = 3,
    inflation = 1.5, minimum = 2, min_omes = 2, 
    minid = 30, cpus = 1, diamond = 'diamond',
    sensitive = True
    ):

    groupI = wrk_dir + 'group.I.pickle'
    groupII = wrk_dir + 'group.II.pickle'
    if not os.path.isfile(groupI):
        print('\tAggregating HGxs with overlapping loci', flush = True)

        hgx_genes = {}
        for hgx in hgx2loc:
            if hgx in hgx2i: # if it passed the border threshold
                for gene in hgx2loc[hgx]:
                    ome = gene[:gene.find('_')]
                    if ome not in hgx_genes:
                        hgx_genes[ome] = {}
                    if gene not in hgx_genes[ome]:
                        hgx_genes[ome][gene] = []
                    hgx_genes[ome][gene].append(hgx)

        hashOgx_cmds = []
        start = datetime.now()
        for ome in hgx_genes:
            hashOgx_cmds.append([
                db[ome]['gff3'],
                ome, hgx_genes[ome], gene2hg, clusplusminus
                ])
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            gene_tups = pool.starmap(hash_hgx, hashOgx_cmds)

        gene2hgx, hgx2genes = {}, {}
        for res in gene_tups:
            gene2hgx[ome2i[res[0]]] = {x[0]: set(x[1]) for x in res[1]}
            for hgx, hgs in res[2]:
                if hgx not in hgx2genes:
                    hgx2genes[hgx] = {og: [] for og in hgx}
                for og, seqs in hgs:
                    hgx2genes[hgx][og].extend(seqs)
            # hgx2genes = {hgx: og: [gene, gene]}
            # gene2hgx = {omeI: gene: {HGx}}
        with open(groupI, 'wb') as out:
            pickle.dump([gene2hgx, hgx2genes], out)
        print('\t\t\t' + str(datetime.now() - start))
    elif not os.path.isfile(groupII):
        print('\tLoading aggregated HGxs', flush = True)
        with open(groupI, 'rb') as in_:
            gene2hgx, hgx2genes = pickle.load(in_)

    groupIII = wrk_dir + 'group.III.pickle'
    if not os.path.isfile(groupIII):
        print('\tIdentifying significant HGx aggregations', flush = True)
        tpairDict = defaultdict(list)
        for ome, gene2hgx_ome in gene2hgx.items():
            tpairDict = find_hgx_pairs(
                gene2hgx_ome, ome, hgx2i, tpairDict,
                minimum = minimum
                )

        pairDict = {}
        for id_, omes in tpairDict.items():
            omes_set = set(omes)
            if len(omes_set) > 1: # need at least one ome to calc a branch length
                pairDict[id_] = tuple(sorted(omes_set))

        omes2dist = phylocalcs.update_dists(phylo, pairDict, cpus = cpus, omes2dist = omes2dist)

        print('\tClassifying HGx clans and gene cluster families (GCFs)', flush = True)
        # populate a lil_matrix here, then use that network to identify modules
        print('\t\tBuilding binary HGx-HGx network', flush = True)
        matrix = lil_matrix((len(i2hgx), len(i2hgx)), dtype=bool)
        for idPair, omes in pairDict.items():
            i0, i1 = idPair[0], idPair[1]
            nullSize = max([len(i2hgx[x]) for x in idPair])
            parts = set(ome2partition[x] for x in omes)
            if None in parts:
                parts = parts.remove(None)
                if not parts:
                    continue
            bord_score = min([bordScores[i][nullSize] for i in list(parts)])
            if omes2dist[omes] >= bord_score:
                matrix[i0, i1] = True
                matrix[i1, i0] = True
        print('\t\tIsolating subgraphs (HGx clans)', flush = True)
        network = nx.from_scipy_sparse_matrix(matrix)
        subgraphs = [network.subgraph(c) for c in nx.connected_components(network)]
        # use .copy() after iterator to copy
        clans, clanOmes, clanHGxs = [], [], []
        for sg in subgraphs:
            clans.append({})
            clanOmes.append([])
            clanHGxs.append([])
            for id_ in list(sg.nodes):
                hgx = i2hgx[id_]
                omes = hgx2omes[hgx]
                clans[-1][hgx] = omes # storing in a hash to allow future
                # manipulation
                clanHGxs[-1].extend(hgx)
                clanOmes[-1].extend(omes)
            clanOmes[-1] = tuple(sorted(set(clanOmes[-1])))
            clanHGxs[-1] = tuple(sorted(set(clanHGxs[-1])))

        omes2dist = phylocalcs.update_dists(
            phylo, {clanHGxs[i]: omes for i, omes in enumerate(clanOmes)},
            cpus = cpus, omes2dist = omes2dist
            )
        with open(groupIII, 'wb') as out:
            pickle.dump([clans, clanOmes, clanHGxs], out)
    else:
        print('\tLoading HGx clans', flush = True)
        with open(groupIII, 'rb') as in_:
            clans, clanOmes, clanHGxs = pickle.load(in_)

    print('\t\t' + str(len(clans)) + ' clans', flush = True)

    print('\tClassifying HGx clans into GCFs', flush = True)
    clanFile = wrk_dir + 'clan2loci.'
    if not os.path.isfile(clanFile + 'hg.json.gz'):
        print('\t\tPreparing loci extraction', flush = True)
        ome2clus2extract = defaultdict(dict)
        for i, clan in enumerate(clans):
            for hgx, omes in clan.items():
                for gene in hgx2loc[hgx]:
                    ome = gene[:gene.find('_')]
                    if gene not in ome2clus2extract[ome]:
                        ome2clus2extract[ome][gene] = [[set(hgx), i]]
                    else:
                        ome2clus2extract[ome][gene].append([set(hgx), i])
                    # {ome: {gene: [[set(hgx), clanI]]}}

        print('\t\tExtracting clan loci', flush = True)
        cmds = []
        for ome, clus2extract in ome2clus2extract.items():
            # prepare commands for locus extraction
            cmds.append([
                ome, db[ome]['gff3'],
                clus2extract, gene2hg, clusplusminus
                ])
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            hash_res = pool.starmap(hash_clan_loci, cmds)

        clanLoci, clanHGx4gcfs, clanHGloci = defaultdict(list), defaultdict(list), defaultdict(list)
        for ome, outclanLoci, outclanHGx, outclanHGloci in hash_res:
            for clan in outclanLoci:
                clanLoci[clan].extend(outclanLoci[clan])
                clanHGx4gcfs[clan].extend(outclanHGx[clan])
                clanHGloci[clan].extend(outclanHGloci[clan])
        write_json(clanLoci, clanFile + 'json.gz')
        write_json(clanHGx4gcfs, clanFile + 'hgx.json.gz')
        write_json(clanHGloci, clanFile + 'hg.json.gz')
    else:
        print('\t\tLoading clan loci', flush = True)
        clanLoci = {int(k): tuple(v) for k,v in sorted(
                read_json(clanFile + 'json.gz').items(), 
                key = lambda x: len(x[1]), reverse = True)}
        clanHGx4gcfs = {int(k): tuple([tuple(i) for i in v]) for k,v in read_json(clanFile + 'hgx.json.gz').items()}
        clanHGloci = {int(k): tuple(v) for k, v in read_json(clanFile + 'hg.json.gz').items()}

    # dict(clanLoci) = {int(clanI): ((ome0_gene0, ome0_gene1,)...,)}
    # dict(clanHGx4gcfs) = {int(clanI): ((og0, og1, ogn,)...,)}
    # dict(clanHGloci) = {int(clanI): ((ome0_gene0_og, ome0_gene1_og,)...,)}

    gcf_dir = wrk_dir + 'gcf/'
    if not os.path.isdir(gcf_dir):
        os.mkdir(gcf_dir)

    print('\t\tCalling GCFs', flush = True)
    print('\t\t\t' + str(sum([len(v) for v in list(clanLoci.values())])) \
        + ' loci', flush = True)

    m, adj_mtr = mp.Manager(), gcf_dir + 'loci.adj.tmp'
    if cpus < 2:
        cpus = 2
    index, cmds = 0, []

    loci, hgxXloci = {}, {}
    if sensitive:
        Q = m.Queue()
        W = mp.get_context('spawn').Process(target=write_adj_matrix, args=(Q, adj_mtr))
        for clanI in clanLoci:
            cmds.append([
                db, clanI, clanLoci[clanI], clanHGloci[clanI],
                hg_dir, gcf_dir, Q, index, diamond, minid
                ])
            for i, locus in enumerate(clanLoci[clanI]):
                loci[index] = clanLoci[clanI][i]
                hgxXloci[index] = clanHGx4gcfs[clanI][i]
                index += 1
    else:
        W = mp.Process(target = write_adj_matrix_noq, args = (adj_mtr, gcf_dir))
        for clanI in clanLoci:
            cmds.append([
                db, clanI, clanLoci[clanI], clanHGloci[clanI],
                hg_dir, hgx_dir, gcf_dir, index, minid
                ])
            for i, locus in enumerate(clanLoci[clanI]):
                loci[index] = clanLoci[clanI][i]
                hgxXloci[index] = clanHGx4gcfs[clanI][i]
                index += 1

    W.start()
    
 #   bigClan = cmds[0][1]
  #  del cmds[0] # remove the first one because it is huge enough to mp on its own
    if not os.path.isfile(gcf_dir + 'loci.adj'):
        if sensitive:
            with mp.get_context('spawn').Pool(processes = cpus - 1) as pool:
                pool.starmap(
                    clan_to_gcf_loci_sensitive, cmds
                    )
            Q.put(None)
        else:
            with mp.get_context('spawn').Pool(processes = cpus - 1) as pool:
                pool.starmap(clan_to_gcf_loci_resume, cmds)
            with open(f'{gcf_dir}done.adj.tmp', 'w') as out:
                pass
        shutil.move(adj_mtr, gcf_dir + 'loci.adj')

    W.join()

    if not os.path.isfile(gcf_dir + 'loci.clus'):
        print('\t\t\tRunning MCL', flush = True)
        if cpus > 5:
            mcl_threads = 10 # cap at 5 processors, ten threads
        else:
            mcl_threads = (cpus * 2) - 1
        MCL(gcf_dir + 'loci.adj', gcf_dir + 'loci.clus',
            inflation = 1.5, threads = mcl_threads)
        # could add iterative subsample MCL option here

    t_gcfs = []
    with open(gcf_dir + 'loci.clus', 'r') as raw:
        for line in raw: # loci indices
            d = line.rstrip().split('\t')
            if len(d) > 1:
                indices = [int(x) for x in d]
            else: # singletons won't pass thresholds, disregard them
                continue
            t_gcfs.append(indices)

    # list(gcf) = [{hgx: (omes,)}]
    gcfs, gcf_hgxs, gcf_omes = [], [], []
    for gcf, locIs in enumerate(t_gcfs):
        gcfs.append(defaultdict(list))
        gcfHGx, gcfOme_list = [], []
        for locI in locIs:
            loc = loci[locI]
            hgxs = tuple([tuple(hgx) for hgx in hgxXloci[locI]])
            # really should be done above to make hgxs formatted right
            omeI = ome2i[loc[0][:loc[0].find('_')]]
            [gcfs[-1][hgx].append(omeI) for hgx in hgxs]
            [gcfHGx.extend(hgx) for hgx in hgxs]
            gcfOme_list.append(omeI)
        if len(set(gcfOme_list)) > min_omes: # need more than 1
            gcfs[-1] = {k: sorted(set(v)) for k,v in gcfs[-1].items()}
            gcf_hgxs.append(tuple(sorted(set(gcfHGx))))
            gcf_omes.append(tuple(sorted(set(gcfOme_list))))
        else:
            del gcfs[-1]

    print('\t\t\t' + str(len(gcfs)) + ' GCFs w/' \
        + str(sum([len(x) for x in gcfs])) + ' loci', flush = True)
    omes2dist = phylocalcs.update_dists(
        phylo, {gcf_hgxs[i]: omes for i, omes in enumerate(gcf_omes)},
        cpus = cpus, omes2dist = omes2dist
        )

    return gcfs, gcf_hgxs, gcf_omes, omes2dist
