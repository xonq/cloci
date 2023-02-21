import os
import sys
import copy
import shutil
import pickle
import subprocess
import networkx as nx
import numpy as np
import multiprocessing as mp
from tqdm import tqdm
from datetime import datetime
from itertools import combinations, chain
from scipy.sparse import lil_matrix, csr_matrix
from collections import defaultdict, Counter
from mycotools.lib.biotools import gff2list
from mycotools.lib.kontools import write_json, read_json, collect_files, \
                                   checkdir, eprint
from orthocluster.orthocluster.lib import input_parsing, treecalcs, evo_conco

def hash_hgx(gff_path, ome, hgx_genes, gene2hg, clusplusminus):
    """sliding window of size clusplusminus, identify sequences that may have
    hgs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus"""
        
    gff_list, gene_dict = gff2list(gff_path), defaultdict(list)
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
                    try:
                        hgx_dict[hgx][og0].append(seq0) # add the sequence to the hgx_dict
                    except KeyError:
                        pass
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


def find_hgx_pairs(gene2hgx, ome, hgx2i, pairsDict, min_hgx_overlap = 2):
    # is this going to merge hgxs that overlap by one gene in min_hgx_overlap loci?
    hgx_pairs_raw = []
    # check for gene overlap in hgxs
    for hgxs in list(gene2hgx.values()):
        # add each possible hgxpair
        hgx_pairs_raw.extend(list(combinations(sorted(hgxs), 2)))

    hgx_pairs_count = [tuple((hgx2i[y] for y in x)) \
                       for x, g in Counter(hgx_pairs_raw).items() \
                       if g >= min_hgx_overlap]
    for hgxpair in hgx_pairs_count:
        pairsDict[hgxpair].append(ome)
    return pairsDict

def MCL(adj_path, clusFile, inflation = 1.5, threads = 1):

    subprocess.call([
        'mcl', adj_path, '--abc', '-I', str(inflation),
        '-o', clusFile, '-te', str(threads)
        ], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
        )

    return clusFile

def sorensen(interlen, len0, len1):
    return 2*interlen/(len0 + len1)
def jaccard(interlen, len0, len1):
    return interlen/(len0 + len1 - interlen)
def overlap(interlen, len0, len1):
    return interlen/min([len0, len1])

def acquire_clus_gcf_sim_noq(
    i0, i1, index, loc0, loc1, hgL0, hgL1, set0, set1, blast_ids,
    simfun = overlap, mingcfid = 0.15
    ):

    intersection = set0.intersection(set1)
    inter_list = list(intersection)
    hg_dict = {k: [[], []] for k in inter_list}
    for i, hg in enumerate(hgL0):
        if hg in intersection:
            hg_dict[hg][0].append(loc0[i])
    for i, hg in enumerate(hgL1):
        if hg in intersection:
            hg_dict[hg][1].append(loc1[i])

    #jaccard = len(intersection)/len(set0.union(set1))
    # currently this is biased against contig edges and will also
    # force fragments with a much smaller subset OG # then the
    # primary HGx into separate gcfs

    # could alternatively weight just by incorporating 0 values
    # for OGs that don't overlap

    scores = []
    for hg in inter_list:
#        if hg_dict[hg][0] and hg_dict[hg][1]: #unnec
        scores.append([])
        for gene0 in hg_dict[hg][0]:
            for gene1 in hg_dict[hg][1]:
                try:
                    scores[-1].append(blast_ids[hg][gene0][gene1])
                except KeyError: # missing gene
                    #adjust overlap coef, not 0
   #                 scores[-1].append(0)
                    pass
    #    else:
    #        scores[-1].append(0)
    scores = [x for x in scores if x]
    if scores:
        maxScores = [max(i) for i in scores]
        coef = simfun(len(scores), len(set0), len(set1))
        gcf_id = (sum(maxScores)/len(maxScores)) * coef
        if gcf_id > mingcfid:
            return str(i0 + index) + '\t' + str(i1 + index) + '\t' + str(gcf_id) + '\n'

    return


def clan_to_gcf_sim(
    db, clanI, loci, hgLoci, hg_dir, hgx_dir, gcf_dir, index,
    minid = 30, mingcfid = 0.15, simfun = overlap, min_overlap = 2
    ):

    finished_file = f'{gcf_dir}{clanI}.sadj.tmp'
    # is this clan already completed?
    if os.path.isfile(finished_file) or os.path.isfile(finished_file + '.r'):
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
                            out.write('\n'.join([' '.join(x) for x in data]) \
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
    hg2i = {}
    for i in range(max_complete, len(loci)):
        locus = loci[i]
        hgs = hgLoci[i]
        for i1, hg in enumerate(hgs):
            if hg is not None:
                blast_hash[hg].append(locus[i1])
                if hg not in hg2i:
                    hg2i[hg] = len(hg2i)
    clan_arr = lil_matrix((len(loci), len(blast_hash)), dtype=bool)
    for i in range(max_complete, len(loci)):
        hgs = hgLoci[i]
        for i1, hg in enumerate(hgs):
            if hg is not None:
                clan_arr[i, hg2i[hg]] = True

    blast_ids = defaultdict(dict)
    for hg, genes in blast_hash.items():
        gene_len = len(genes)
        if gene_len > 1:
            gene_set = set(genes)
            algn_base = f'{hgx_dir}{hg}.out'
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
                 blast_ids = run_hgx_blast(blast_ids, hg_dir, hg, genes, 
                              hgx_dir, clanI, minid = minid,
                              diamond = 'diamond')  

    # could be expedited by multiprocessing a writer and numba njit the math
    if blast_ids:
        clan_arr = clan_arr.tocsr()
        adj_arr = clan_arr @ clan_arr.transpose()
        hits = np.split(adj_arr.indices, adj_arr.indptr)[1:-1]
        with open(f'{finished_file}.w', 'a') as out:
            for ti0, overlap_loci in enumerate(hits):
                i0 = ti0 + max_complete
                try:
                    loc0, hgL0 = loci[i0], hgLoci[i0]
                except IndexError:
                    continue
                sHGl0 = set([x for x in hgL0 if x is not None])
                for ti1 in overlap_loci[:]:
                    if ti1 > ti0:
                        i1 = ti1 + max_complete
                        loc1, hgL1 = loci[i1], hgLoci[i1]
                        sHGl1 = set([x for x in hgL1 if x is not None])
                        data = acquire_clus_gcf_sim_noq(i0, i1, index, loc0,
                                                        loc1, hgL0, hgL1,
                                                        sHGl0, sHGl1, blast_ids,
                                                        mingcfid = mingcfid,
                                                        simfun = simfun)
                   
                        if data:
                            out.write(data)
                    out.flush()
        # begin appending to temporary writing output
        # the file is now complete
        shutil.move(f'{gcf_dir}{clanI}.sadj.tmp.w', f'{gcf_dir}{clanI}.sadj.tmp')

 
def run_hgx_blast(blast_ids, hg_dir, hg, genes, 
                  hgx_dir, clanI, minid = 30,
                  diamond = 'diamond'):
    hg_file = f'{hg_dir}{hg}.faa'
    fileBase = hgx_dir + str(hg)
    if not os.path.isfile(fileBase + '.out'):
#        if not os.path.isfile(fileBase + '.dmnd'):
        makeDBcmd = subprocess.call([
            diamond, 'makedb', '--in', hg_file, 
            '--db', fileBase + '.dmnd', '--threads', str(2)
            ], stdout = subprocess.DEVNULL)# just remake the database
        if not os.path.isfile(fileBase + '.out'):
            dmndBlast = subprocess.call([diamond, 'blastp', '--query', 
                hg_file, '--db', fileBase, '--threads', str(2), 
                '--id', str(minid), '--no-self-hits', '-o', 
                fileBase + '.out.tmp', '--outfmt', '6', 'qseqid', 'sseqid', 
                'pident', 'ppos'
                ], stdin = subprocess.PIPE, stdout = subprocess.DEVNULL)#, 
#                 stderr = subprocess.DEVNULL)
            shutil.move(fileBase + '.out.tmp', fileBase + '.out')
    
    with open(fileBase + '.out', 'r') as raw:
        for line in raw:
            d = line.rstrip().split()
            q, s, pident = d[0], d[1], d[-2]
            if q not in blast_ids[hg]:
                blast_ids[hg][q] = {}
            blast_ids[hg][q][s] = float(pident)/100 # adjust diamond to decimal

    return blast_ids


def merge_clan_loci(preclan_loci, gene2hg):
    out_clan_loci, outclanHGx = defaultdict(list), defaultdict(list)
    for clan, loci in preclan_loci.items():
        while loci: # exhaustively merge overlapping loci
            # the locus for pairwise comparisons is the first
            loc0, hgxs0, allHGxs = set(loci[0][0]), loci[0][1], loci[0][2]
            locIntersect = None
            # for all other loci
            for i1, loc1d in enumerate(loci[1:]):
                loc1, hgxs1, allHGxs1 = set(loc1d[0]), loc1d[1], loc1d[2]
                locIntersect = loc0.intersection(loc1)
                # if any genes overlap between these loci
                if locIntersect:
                    # grab the overlap and add it to the end of the loci list
                    newHGx = list(hgxs0)
                    newHGx.extend(list(hgxs1))
                    allHGxs.extend(allHGxs1)
                    loci.append([list(loc1.union(loc0)), newHGx, allHGxs])
                    break
            if locIntersect: # if there was overlap, delete the overlappers
                del loci[0]
                del loci[i1] # don't add 1 because we removed the first and i1 relative to 1
            else: # no more overlap for this locus, add to the final output
                out_clan_loci[clan].append(sorted(loc0))
                outclanHGx[clan].append(tuple(sorted(set(allHGxs))))
#                outclanHGx[clan].append(tuple(sorted(set(hgxs0))))
                del loci[0]

    out_hg_loci = {}
    for clan, loci in out_clan_loci.items():
        out_hg_loci[clan] = []
        for locus in loci:
            out_hg_loci[clan].append([])
            for gene in locus:
                try:
                    out_hg_loci[clan][-1].append(gene2hg[gene])
                except KeyError:
                    out_hg_loci[clan][-1].append(None)

    return out_clan_loci, outclanHGx, out_hg_loci


def merge_clan_loci_heat(preclan_loci, gene2hg, min_sim, hgx2dist_path):
    with open(hgx2dist_path, 'rb') as raw:
        hgx2dist = pickle.load(raw)

    gene2hgx = defaultdict(list)
    out_clan_loci, outclanHGx = defaultdict(list), defaultdict(list)
    for clan, loci in preclan_loci.items():
        for i, loc_d in enumerate(loci):
            loc, hgx, null = loc_d
            loci[i].append([hgx2dist[hgx]])
            [gene2hgx[x].append(hgx) for x in loc]
        merge, check_loci = True, loci
        while merge: # need to rerun to verify
            loci = copy.deepcopy(check_loci)
            check_loci, merge = [], False
            while loci: # exhaustively merge overlapping loci
                # the locus for pairwise comparisons is the first
                loc0, hgxs0, allHGxs, dists0 = set(loci[0][0]), loci[0][1], loci[0][2], loci[0][3]
                locIntersect = None
                # for all other loci
                for i1, loc1d in enumerate(loci[1:]):
                    loc1, hgxs1, allHGxs1, dists1 = set(loc1d[0]), loc1d[1], loc1d[2], loc1d[3]
                    locIntersect = loc0.intersection(loc1)
                    # if any genes overlap between these loci
                    if locIntersect:
                        # compile the distances of the overlapping HGx genes
                        inter_hgxs = []
                        [inter_hgxs.extend(gene2hgx[gene]) for gene in list(locIntersect)]
                        inter_hgx0 = [x for x in inter_hgxs if x in set(allHGxs)]
                        inter_hgx1 = [x for x in inter_hgxs if x in set(allHGxs1)]
                        inter_dists0 = [hgx2dist[x] for x in inter_hgx0]
                        inter_dists1 = [hgx2dist[x] for x in inter_hgx1]
                        min0, max0 = min(inter_dists0), max(inter_dists0)
                        min1, max1 = min(inter_dists1), max(inter_dists1)
                        # do the ranges of observed distances overlap?
                        if min1 <= min0 <= max1 or min1 <= max0 <= max1 \
                            or min0 <= min1 <= max0 or min0 <= max1 <= max0:  
                            # grab the overlap and add it to the end of the loci list
                            newHGx = list(hgxs0)
                            newHGx.extend(list(hgxs1))
                            allHGxs.extend(allHGxs1)
                            newdists = dists0 + dists1
                            loci.append([list(loc1.union(loc0)), newHGx, allHGxs, newdists])
                            break
                        else:
                            # determine the closest bounds of the two distance lists
                            min_list, max_list = [min0, min1], [max0, max1]
                            sml_up_bound = min(max_list)
                            lrg_low_bound = max(min_list)
                            # are they similar enough?
                            if sml_up_bound/lrg_low_bound >= min_sim:
                                # grab the overlap and add it to the end of the loci list
                                newHGx = list(hgxs0)
                                newHGx.extend(list(hgxs1))
                                allHGxs.extend(allHGxs1)
                                newdists = dists0 + dists1
                                loci.append([list(loc1.union(loc0)), newHGx, allHGxs, newdists])
                                break
                            else:
                                locIntersect = None
                if locIntersect: # if there was overlap, delete the overlappers
                    merge = True
                    del loci[0]
                    del loci[i1] # don't add 1 because we removed the first and i1 relative to 1
                else:
                    check_loci.append(loci[0])
                    del loci[0]

        loci = check_loci
        # give the remaining overlapping genes to the highest distance group
        # will need to add the new HGxs in the end
        for i0, locd0 in enumerate(loci):
            loc0, hgx0, ahgx0, dist0 = set(locd0[0]), locd0[1], locd0[2], max(locd0[-1])
            for i1, locd1 in enumerate(loci[i0+1:]):
                loc1, hgx1, ahgx1, dist1 = set(locd1[0]), locd1[1], locd1[2], max(locd1[-1])
                locIntersection = loc0.intersection(loc1)
                if locIntersection:
# currently this doesnt check if there are others of the same hg that are removed in the
# overlap - ideally, this should check if there are other members of the same hg before
# removing them IF they are not part of the overlap
                    hg_intersection = set(gene2hg[gene] for gene in list(locIntersection) \
                                          if gene in gene2hg)
                    if dist0 > dist1:
                        loci[i1 + 1 + i0][0] = list(loc1.difference(loc0))
 #                       loci[i1 + 1 + i0][1] = list(set(hgx1).difference(hg_intersection))
                        loci[i1 + 1 + i0][2] = [tuple(sorted(set(x).difference(hg_intersection))) \
                                                for x in ahgx1]
                    else:
                        new_loc = list(loc0.difference(loc1))
#                        hgx0 = list(set(hgx0).difference(hg_intersection))
                        ahgx0 = [tuple(sorted(set(x).difference(hg_intersection))) \
                                 for x in ahgx0]
                        loci[i0] = [new_loc, hgx0, ahgx0, [dist0]]
                        loc0 = set(new_loc)

        for loc, hgx, ahgx, null in loci:
            out_loc = [loc, [x for x in ahgx if len(x) > 2]]
            if out_loc[1] and len(out_loc[0]) > 2:
                out_clan_loci[clan].append(sorted(out_loc[0]))
                outclanHGx[clan].append(tuple(sorted(set(out_loc[1]))))
    #                outclanHGx[clan].append(tuple(sorted(set(hgxs0))))

    out_hg_loci = {}
    for clan, loci in out_clan_loci.items():
        out_hg_loci[clan] = []
        for locus in loci:
            out_hg_loci[clan].append([])
            for gene in locus:
                try:
                    out_hg_loci[clan][-1].append(gene2hg[gene])
                except KeyError:
                    out_hg_loci[clan][-1].append(None)

    return out_clan_loci, outclanHGx, out_hg_loci



def hash_clan_loci(ome, gff_path, ome_sig_clus, gene2hg, clusplusminus,
                   min_merge_perc = 0, hgx2dist_path = None):
    """sliding window of size clusplusminus, identify sequences that may have
    hgs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus"""

    gff_list, protoclus, clus_out = gff2list(gff_path), {}, []
    cds_dict = input_parsing.compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))

    preclan_loci = defaultdict(list)
    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]):
            if seq0 in ome_sig_clus: # if the og is part of a significant seed
            # locus [[set(hgx)], i]
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
                    preclan_loci[clan].append([
                         locus[start:end], hgx, [hgx]
                         ])

    if not min_merge_perc:
        out_clan_loci, outclanHGx, out_hg_loci = merge_clan_loci(preclan_loci, gene2hg)
    else:
        out_clan_loci, outclanHGx, out_hg_loci = merge_clan_loci_heat(preclan_loci, gene2hg, 
                                                                      min_merge_perc, hgx2dist_path)
    return ome, out_clan_loci, outclanHGx, out_hg_loci


def read_to_write(out_h, in_adj, min_gcf):
    with open(in_adj, 'r') as raw:
        for line in raw:
            d = line.rstrip().split()
            gcf_id = round(float(d[2]) * 100)
            if gcf_id >= min_gcf:
                out_h.write(line.rstrip() + '\n')
    out_h.flush()


def write_adj_matrix(out_file, gcf_dir, min_gcf_id):
    comp_gcf_id = round(float(min_gcf_id) * 100)
    with open(out_file + '.tmp', 'w') as out:
        # get previously read adjacency matrices
        init_files = collect_files(gcf_dir, 'sadj.tmp.r')
        for f in init_files:
            read_to_write(out, f, comp_gcf_id)

        # get ready to read adjacency matrices
        files = collect_files(gcf_dir, 'sadj.tmp')
        f_set = set(files)

        # while the completion signal has not been received
        while f'{gcf_dir}done.sadj.tmp' not in f_set:
            for f in files:
                read_to_write(out, f, comp_gcf_id)
                shutil.move(f, f'{f}.r')
            files = collect_files(gcf_dir, 'sadj.tmp')
            f_set = set(files)

        # get final files
        for f in files:
            if f == f'{gcf_dir}done.sadj.tmp':
                os.remove(f'{gcf_dir}done.sadj.tmp')
            else:
                read_to_write(out, f, comp_gcf_id)
                shutil.move(f, f'{f}.r')
    shutil.move(out_file + '.tmp', out_file)
    for f in collect_files(gcf_dir, 'sadj.tmp.r'):
        os.remove(f)


def check_tune(tune, gcf_hgxs, gcf_omes):
    failed = []
    satisfied = False
    for name, data in tune.items():
        hgx, omes, false_omes = data[0], data[1], data[2]
        set_hgx = set(hgx)
        set_omes = set(omes)
        set_false = set(false_omes)
        check = [v for i, v in enumerate(gcf_hgxs) \
                 if set_hgx.issubset(set(v)) \
                 and set_omes.issubset(gcf_omes[i]) \
                 and not set_false.intersection(gcf_omes[i])]
        if check:
            continue
        else:
            failed.append(name)
    if failed:
        for name in failed:
            print(f'\t\t\t{name} was not recovered', flush = True)
    else:
        print(f'\t\t\tTuned successfully', flush = True)
        satisfied = True
    return satisfied


def mcl2gcfs(gcf_dir, loci, hgxXloci, ome2i, inflation, 
             loc2clan, tune = False, min_omes = 1, cpus = 1):
    satisfied = False
    if not os.path.isfile(gcf_dir + 'loci.clus') or tune:
        print(f'\t\t\tRunning MCL - inflation: {inflation}', flush = True)
        if cpus > 5:
            mcl_threads = 10 # cap at 5 processors, ten threads
        else:
            mcl_threads = (cpus * 2) - 1
        MCL(gcf_dir + 'loci.adj', gcf_dir + 'loci.clus',
            inflation = inflation, threads = mcl_threads)
        # could add iterative subsample MCL option here

    t_gcfs = []
    with open(gcf_dir + 'loci.clus', 'r') as raw:
        for line in raw: # loci indices
            indices = [int(x) for x in line.rstrip().split()]
            # singletons won't pass thresholds, disregard them
            if len(indices) > 1:
                t_gcfs.append(indices)

    # list(gcf) = [{hgx: (omes,)}]
    gcfs, gcf_hgxs, gcf_omes, gcf2clan = [], [], [], {}
    with open(f'{gcf_dir}loci.txt', 'w') as out:
        out.write('#omeI locI clanI pregcf locus\n')
        for gcf, locIs in enumerate(t_gcfs):
#            gcfs.append(defaultdict(list))
            gcfs.append([])
            gcf_hgxs.append([])
            gcfOme_list = [] 
            locs = []
            for locI in locIs:
                loc = loci[locI]
                locs.append((loc, locI))
                hgxs = tuple([tuple(hgx) for hgx in hgxXloci[locI]])
                # really should be done above to make hgxs formatted right
                omeI = ome2i[loc[0][:loc[0].find('_')]]
#                [gcfs[-1][hgx].append(omeI) for hgx in hgxs]
                gcfs[-1].append(loc)
                gcfOme_list.append(omeI)
                gcf_hgxs[-1].extend(hgxs)
            if len(set(gcfOme_list)) > min_omes: # need more than 1
#                gcfs[-1] = {k: sorted(set(v)) for k,v in gcfs[-1].items()}
                if gcfs[-1]:
                    gcfs[-1] = tuple([tuple(x) for x in gcfs[-1]])
                    gcf_hgxs[-1] = tuple(sorted(set(chain(*gcf_hgxs[-1]))))
                    gcf_omes.append(tuple(sorted(set(gcfOme_list))))
#                    gcf_hgxs.append(tuple(sorted(set(chain(*list(gcfs[-1].keys()))))))
  #                  gcf_omes.append(tuple(sorted(set(chain(*list(gcfs[-1].values()))))))
                    pregcf = len(gcfs) - 1
                    gcf2clan[pregcf] = loc2clan[locs[0][1]]
                    for i, omeI in enumerate(gcfOme_list):
                        locus, locI = locs[i]
                        out.write(f'{omeI} {locI} {loc2clan[locI]} {pregcf} {",".join(locus)}\n')
                else:
                    del gcfs[-1]
                    del gcf_hgxs[-1]
            else:
                del gcfs[-1]
                del gcf_hgxs[-1]

    failed = []
    if tune:
        satisfied = check_tune(tune, gcf_hgxs, gcf_omes)
    else:
        satisfied = True

    list_gcf2clan = []
    for i in range(len(gcfs)):
        list_gcf2clan.append(gcf2clan[i])

    return satisfied, {i: v for i, v in enumerate(gcfs)}, \
           {i: v for i, v in enumerate(gcf_hgxs)}, \
           {i: v for i, v in enumerate(gcf_omes)}, \
           {i: v for i, v in enumerate(list_gcf2clan)}


def classify_gcfs(
    hgx2loc, db, gene2hg, i2hgx, hgx2i,
    phylo, ome2i, hgx2omes, hg_dir, hgx_dir,
    wrk_dir, ome2partition, bord_scores_list,
    hg2gene, omes2dist = {}, clusplusminus = 3,
    inflation = 1.5, min_hgx_overlap = 1, min_omes = 2, 
    minid = 30, cpus = 1, diamond = 'diamond',
    min_gcf_id = 0.3, simfun = overlap, printexit = False,
    tune = False, dist_func = treecalcs.calc_tmd,
    uniq_sp = False, min_merge_perc = 0
    ):

    i2ome = {v: k for k, v in ome2i.items()}

    if not checkdir(hgx_dir, unzip = True, rm = True):
        os.mkdir(hgx_dir)

    groupI = wrk_dir + 'group.I.pickle'
    groupII = wrk_dir + 'group.II.pickle'

    if not os.path.isfile(groupI):
        print('\tAggregating HGxs with overlapping loci', flush = True)
        hgx_genes = {}
        # compile a hash that shows the HGx each gene corresponds to
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
        # expand the hgx 2 gene hash to include all genes in the locus corresponding
        # to the hgx
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            gene_tups = pool.starmap(hash_hgx, hashOgx_cmds)
            pool.close()
            pool.join()    

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
        print('\t\t' + str(datetime.now() - start))
    else:
        print('\tLoading aggregated HGxs', flush = True)
        with open(groupI, 'rb') as in_:
            gene2hgx, hgx2genes = pickle.load(in_)


    if not os.path.isfile(groupII):
        print(f'\tIdentifying HGxs with at least {min_hgx_overlap} gene overlap', 
              flush = True)
        tpairDict = defaultdict(list)
        for ome, gene2hgx_ome in tqdm(gene2hgx.items(),
                                      total = len(gene2hgx)):
            tpairDict = find_hgx_pairs(
                gene2hgx_ome, ome, hgx2i, tpairDict,
                min_hgx_overlap = min_hgx_overlap
                )

        pairDict = {}
        for id_, omes in tpairDict.items():
            omes_set = set(omes)
   #         if len(omes_set) > 1: # need at least one ome to calc a branch length
            pairDict[id_] = tuple(sorted(omes_set))
#        omes2dist = treecalcs.update_dists(phylo, pairDict, cpus = cpus, omes2dist = omes2dist,
 #                                          func = dist_func, uniq_sp = uniq_sp, i2ome = i2ome)


        print('\tClassifying HGx clans and gene cluster families (GCFs)', flush = True)
        # populate a lil_matrix here, then use that network to identify modules
        print('\t\tBuilding binary HGx-HGx network', flush = True)
        matrix = lil_matrix((len(i2hgx), len(i2hgx)), dtype=bool)
        for idPair, omes in tqdm(pairDict.items(), total = len(pairDict)):
            i0, i1 = idPair[0], idPair[1]
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

        omes2dist = treecalcs.update_dists(
            phylo, {clanHGxs[i]: omes for i, omes in enumerate(clanOmes)},
            cpus = cpus, omes2dist = omes2dist, func = dist_func, uniq_sp = uniq_sp,
            i2ome = i2ome
            )
        with open(groupII, 'wb') as out:
            pickle.dump([clans, clanOmes, clanHGxs], out)
    else:
        print('\tLoading HGx clans', flush = True)
        with open(groupII, 'rb') as in_:
            clans, clanOmes, clanHGxs = pickle.load(in_)

    print('\t\t' + str(len(clans)) + ' clans', flush = True)
    print('\tClassifying HGx clans into GCFs', flush = True)
    clanFile = wrk_dir + 'clan2loci.'
    if not os.path.isfile(clanFile + 'hg.json.gz'):
        print('\t\tPreparing loci extraction', flush = True)
        ome2clus2extract = defaultdict(dict)
        for i, clan in enumerate(clans):
            for hgx, omes in clan.items():
                # create a seed hash of the hgx and its seed gene
                for gene in hgx2loc[hgx]:
                    ome = gene[:gene.find('_')]
                    if gene not in ome2clus2extract[ome]:
                        ome2clus2extract[ome][gene] = [[set(hgx), i]]
                    else:
                        ome2clus2extract[ome][gene].append([set(hgx), i])
                    # {ome: {gene: [[set(hgx), clanI]]}}

        print('\t\tExtracting clan loci', flush = True)
        if not min_merge_perc:
            with mp.get_context('fork').Pool(processes = cpus) as pool:
                hash_res = pool.starmap(hash_clan_loci, 
                                        tqdm(((ome, db[ome]['gff3'], clus2extract,
                                              gene2hg, clusplusminus) \
                                              for ome, clus2extract \
                                              in ome2clus2extract.items()),
                                             total = len(ome2clus2extract)))
                pool.close()
                pool.join()
        else:
             with mp.get_context('fork').Pool(processes = cpus) as pool:
                hash_res = pool.starmap(hash_clan_loci, 
                                        tqdm(((ome, db[ome]['gff3'], clus2extract,
                                              gene2hg, clusplusminus, min_merge_perc,
                                              f'{wrk_dir}hgx2dist.pickle') \
                                              for ome, clus2extract \
                                              in ome2clus2extract.items()),
                                             total = len(ome2clus2extract)))
                pool.close()
                pool.join()
    
        clanLoci = defaultdict(list)
        clanHGx4gcfs = defaultdict(list)
        clanHGloci = defaultdict(list)
        for ome, out_clan_loci, outclanHGx, outclanHGloci in hash_res:
            for clan in out_clan_loci:
                clanLoci[clan].extend(out_clan_loci[clan])
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
        clanHGx4gcfs = {int(k): tuple([tuple(i) for i in v]) \
                        for k,v in read_json(clanFile + 'hgx.json.gz').items()}
        clanHGloci = {int(k): tuple(v) \
                      for k, v in read_json(clanFile + 'hg.json.gz').items()}

    # dict(clanLoci) = {int(clanI): ((ome0_gene0, ome0_gene1,)...,)}
    # dict(clanHGx4gcfs) = {int(clanI): ((og0, og1, ogn,)...,)}
    # dict(clanHGloci) = {int(clanI): ((ome0_gene0_og, ome0_gene1_og,)...,)}

    # need to move up to before writing
    clanLoci = {k: v for k,v in sorted(clanLoci.items(), key = lambda x: len(x[1]), reverse = True)}
    clanHGx4gcfs = [clanHGx4gcfs[k] for k in clanLoci]
    clanHGloci = [clanHGloci[k] for k in clanLoci]
    clanLoci = [v for v in clanLoci.values()]

    gcf_dir = wrk_dir + 'gcf/' 
    if not os.path.isdir(gcf_dir):
        os.mkdir(gcf_dir)

    print('\t\tOutputting HG fastas', flush = True)
#    hgs_prep = set(chain(*list(chain(*list(clanHGloci.values())))))
    hgs_prep = set(chain(*list(chain(*clanHGloci))))
    hgs = sorted([x for x in hgs_prep if x is not None])
    hg_dir = input_parsing.hg_fa_mngr(wrk_dir, hg_dir, 
                             hgs, db, hg2gene, cpus = cpus)
    evo_conco.run_blast(hgs, db, hg_dir, hgx_dir, cpus = cpus,
                        diamond = diamond, printexit = printexit)


    print('\t\tCalling GCFs', flush = True)
    clan_lens = [len(v) for v in clanLoci]
#    print('\t\t\t' + str(sum([len(v) for v in list(clanLoci.values())])) \
    print(f'\t\t\t{sum(clan_lens)}' \
        + ' loci', flush = True)
    print(f'\t\t\t{max(clan_lens)} loci in largest clan', flush = True)
    print('\t\t\tNOTE: for large datasets, clan 0 will hang progress bar', flush = True)
    m, adj_mtr = mp.Manager(), gcf_dir + 'loci.adj'
    if cpus < 2:
        cpus = 2
    index, cmds = 0, []


    # if the min_gcf_id is greater than 33% then 
    # at least HGx overlaps < 2 can be eliminated
    min_overlap = round((min_gcf_id * 3) - 0.5)

    loci, hgxXloci, loc2clan = {}, {}, {}
    for clanI, locs in enumerate(clanLoci):
        cmds.append([
            db, clanI, locs, clanHGloci[clanI],
            hg_dir, hgx_dir, gcf_dir, index, minid, min_gcf_id,
            simfun, min_overlap
            ])
        for i, locus in enumerate(locs):
            loci[index] = clanLoci[clanI][i]
            hgxXloci[index] = clanHGx4gcfs[clanI][i]
            loc2clan[index] = clanI
            index += 1


    if not os.path.isfile(adj_mtr):
        W = mp.Process(target = write_adj_matrix, 
                       args = (adj_mtr, gcf_dir, min_gcf_id))
        W.start()
        with mp.get_context('forkserver').Pool(processes = cpus - 1) as pool:
            pool.starmap(clan_to_gcf_sim, tqdm(cmds, total = len(cmds)))
            pool.close()
            pool.join()    
        with open(f'{gcf_dir}done.sadj.tmp', 'w') as out:
            pass
        W.join()


    satisfied = False
    if tune:
        inflation = 1.1
    while not satisfied:
        if tune:
            inflation += 0.1
            if inflation > 3:
                break
        satisfied, gcfs, gcf_hgxs, gcf_omes, gcf2clan = mcl2gcfs(gcf_dir, loci,
                                                           hgxXloci, ome2i, 
                                                           inflation, loc2clan,
                                                           tune, min_omes,
                                                           cpus)
    if tune:
        if not satisfied:
            eprint('\nERROR: tuning failed to recover clusters', flush = True)
            sys.exit(135)
        inf_strp = str(inflation * 10)[:2]
        inf_str = inf_strp[0] + '.' + inf_strp[1]
        with open(f'{wrk_dir}../log.txt', 'r') as raw:
            d = raw.read().replace('inflation\tNone', 'inflation\t' + inf_str)
        with open(f'{wrk_dir}../log.txt', 'w') as out:
            out.write(d)

    print('\t\t\t' + str(len(gcfs)) + ' GCFs', flush = True)
    omes2dist = treecalcs.update_dists(
        phylo, {gcf_hgxs[i]: omes for i, omes in gcf_omes.items()},
        cpus = cpus, omes2dist = omes2dist, func = dist_func, 
        uniq_sp = uniq_sp, i2ome = i2ome
        )

    return gcfs, gcf_hgxs, gcf_omes, gcf2clan, omes2dist
