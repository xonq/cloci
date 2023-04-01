import os
import sys
import copy
import shutil
import pickle
import subprocess
import numpy as np
import multiprocessing as mp
from tqdm import tqdm
from graph_tool.all import *
from datetime import datetime
from itertools import combinations, chain
from scipy.sparse import lil_matrix, csr_matrix, save_npz, load_npz, triu
from collections import defaultdict, Counter
from mycotools.lib.biotools import gff2list, dict2fa
from mycotools.lib.kontools import write_json, read_json, collect_files, \
                                   checkdir, eprint
from orthocluster.orthocluster.lib import input_parsing, treecalcs, evo_conco, output_data


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


def read_mci_rows(rows_file):
    mci2pre = {}
    with open(rows_file, 'r') as raw:
        for line in raw:
            mci, pre = [int(x) for x in line.rstrip().split()]
            mci2pre[mci] = pre
    return mci2pre


def MCL(adj_path, clus_file, inflation = 1.5, threads = 1):

    mci_file = adj_path[:-4] + '.mci'
    rows_file = os.path.dirname(adj_path) + '/mcl_rows.tsv'
    if not os.path.isfile(mci_file):
        subprocess.call(['mcxload', '-abc', adj_path, '-o',
                         mci_file + '.tmp', '--write-binary',
                         '-write-tab', rows_file, '--stream-mirror'],
                         stdout = subprocess.DEVNULL, 
                         stderr = subprocess.DEVNULL)
        shutil.move(mci_file + '.tmp', mci_file)
    subprocess.call([
        'mcl', mci_file, '-I', str(inflation),
        '-o', clus_file + '.tmp', '-te', str(threads)
        ], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
        )
    subprocess.call(['mcxdump', '-icl', clus_file + '.tmp',
                     '-o', clus_file, '--dump-pairs'],
                    stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    os.remove(clus_file + '.tmp')


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

    # Identify the overlapping HGs
    intersection = set0.intersection(set1)
    inter_list = list(intersection)
    hg_dict = {k: [[], []] for k in inter_list}

    # populate an hg2gene dict for ome0 and ome1
    for i, hg in enumerate(hgL0):
        if hg in intersection:
            hg_dict[hg][0].append(loc0[i])
    for i, hg in enumerate(hgL1):
        if hg in intersection:
            hg_dict[hg][1].append(loc1[i])

    scores = []
    for hg in inter_list:
        # check to verify the genes are recovered
        if hg in hg_dict:
            scores.append([])
            # for each pairwise gene comparison
            for gene0 in hg_dict[hg][0]:
                for gene1 in hg_dict[hg][1]:
                    try:
                        scores[-1].append(blast_ids[hg][gene0][gene1])
                    # if the genes are missing its a 0
                    # the alignment algorithm MUST output enough alignments
                    # for this is to reliably work on widely dispersed clusters
                    except KeyError: 
                        scores[-1].append(0)
    scores = [x for x in scores if x]
    if scores:
        maxScores = [max(i) for i in scores]
        coef = simfun(len(scores), len(set0), len(set1))
        gcf_id = (sum(maxScores)/len(maxScores)) * coef
        if gcf_id > mingcfid:
            return f'{i0 + index}\t{i1 + index}\t{gcf_id}\n'

    return


def rnd2_gen_blastids_mp(hg, group_dict, hgx_dir, grp_dir, minid):#, blast_ids_mngr):
    gene_set = set(chain(*list(group_dict.values())))
    gene_len = len(gene_set)
    algn_base = f'{hgx_dir}{hg}.out'
    json_out = f'{grp_dir}/algn/{hg}.json'
    if not os.path.isfile(json_out):
        try:
            hg_dict = parse_blast(algn_base, gene_len, gene_set, minid)
            new_hg_dict = defaultdict(dict)
            for g, genes in group_dict.items():
                for g0, g1 in combinations(sorted(genes), 2):
                    try:
                        new_hg_dict[g0][g1] = hg_dict[g0][g1]
                    except KeyError:
                        pass
                    try:
                        new_hg_dict[g1][g0] = hg_dict[g1][g0] 
                    except KeyError:
                        pass
            write_json(new_hg_dict, json_out + '.tmp')
            os.rename(json_out + '.tmp', json_out)
        except FileNotFoundError:
            pass



def gen_blastids_mp(hg, genes, hgx_dir, grp_dir, minid):#, blast_ids_mngr):
    # CURRENTLY NOT SETUP FOR ANY OTHER CLANS THAN 0
    # would need to incorporate clanI in the output json
    # or else it runs the risk of using a failed runs' old jsons
    # cannot multiprocess the runs either
    gene_set = set(genes)
    gene_len = len(gene_set)
    algn_base = f'{hgx_dir}{hg}.out'
    json_out = f'{grp_dir}/algn/{hg}.json'
    if not os.path.isfile(json_out):
        try:
            hg_dict = parse_blast(algn_base, gene_len, gene_set, minid)
            if hg_dict is not None:
                write_json(hg_dict, json_out + '.tmp')
                os.rename(json_out + '.tmp', json_out)
    #        blast_ids_mngr[hg].update(hg_dict)
        except FileNotFoundError:
            pass
 #       print(f'\t\t\t\t\tAttempting to BLAST {hg}', flush = True)
#        blast = run_hgx_blast(hg_dir, hg, genes, hgx_dir, 
  #                    blastp = 'blastp', cpus = 1)
   #     if blast:
    #        blast_ids_mngr[hg].update(parse_blast(algn_base, 
     #                                             gene_len, gene_set))


def read_gcf_sim_jsons(gcf_dir):
    jsons = collect_files(gcf_dir, 'json')
    blast_ids = {}
    for j_f in jsons:
        hg = os.path.basename(j_f).replace('.json','')
        blast_ids[int(hg)] = read_json(j_f)
    return blast_ids


def parse_blast(algn_base, gene_len, gene_set, minid = 30):
    hg_dict = defaultdict(dict)
    with open(algn_base, 'r') as raw:
        for line in raw:
            d = line.rstrip().split()
            # adjust to just reference d when safe
            try:
                q, s, pident = d[0], d[1], d[-1]
            except IndexError:
                eprint(f'ERROR: malformatted alignment: {os.path.basename(algn_base)}', 
                       flush = True) 
                return False
            if {q, s}.issubset(gene_set):
                # NEED a no self hits option in mmseqs or blast to avoid
                if q != s:
                    pid = float(pident)
                    if pid >= minid:
                        hg_dict[q][s] = pid/100 # adjust diamond to decimal
            if len(hg_dict) == gene_len:
                break
    return hg_dict


def gcf_sim_write_mngr(Q, sim_file):
    x = Q.get()
    with open(sim_file, 'a') as out:
        while x:
            out.write(x)
            x = Q.get()
            out.flush()


def calc_sim_by_row(adj_mat, hgLoci, loci, blast_ids, 
                    index, Q, max_complete, mingcfid, simfun):
 #   adj_mat = np.load(f'{mat_dir}{rowI}.npz')
    for ti0, overlap_loci in enumerate(adj_mat):
        i0 = ti0 + max_complete
        print(hgLoci, overlap_loci)
        loc0, hgL0 = loci[i0], hgLoci[i0]
        sHGl0 = set([x for x in hgL0 if x is not None])
        for ti1 in overlap_loci:
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
                    Q.put(data)



def gcf_sim_mp_mngr(adj_mat, sim_file, max_complete, loci, hgLoci, index, 
                    blast_ids, gcf_dir, mat_dir, mingcfid, simfun, 
                    chunkscale = 100, cpus = 1):
    if not os.path.isdir(mat_dir):
        os.mkdir(mat_dir)
    
    # break up the adj_mat by chunksize and offload to file while preparing cmds
    row_len = len(adj_mat)
    chunksize = round(row_len/chunkscale)
    cmds = []
    Q = mp.Manager().Queue()
    for i, row in enumerate(range(0, row_len, chunksize)):
        if len(adj_mat) >= row + chunksize:
            sub_mat = [list(adj_mat[x]) for x in range(row, row + chunksize)]
        else:
            sub_mat = [list(adj_mat[x]) for x in range(row, len(adj_mat))]
#        np.save(f'{mat_dir}{i}.npz', sub_mat)
        locIs = [x for x in list(set(chain(*sub_mat))) if x >= row]
        t_hg_loc = {locI: hgLoci[locI + max_complete] for locI in locIs}
        t_loc = {locI: loci[locI + max_complete] for locI in locIs}
        t_hgs = set(chain(*list(t_hg_loc.values())))
        if None in t_hgs:
            t_hgs.remove(None)
        t_ids = {hg: blast_ids[hg] for hg in list(t_hgs) if hg in blast_ids}
        cmds.append((sub_mat, t_hg_loc, t_loc, t_ids, 
                     index, Q, max_complete, mingcfid, simfun))

    write_proc = mp.Process(target=gcf_sim_write_mngr, args=(Q, sim_file))
    write_proc.start()

    with mp.get_context('forkserver').Pool(processes = cpus - 1) as pool:
        pool.starmap(calc_sim_by_row, tqdm(cmds, total = len(cmds)))
        pool.close()
        pool.join()

    Q.put(None)
    write_proc.close()
    write_proc.join()


def loc2loc_sim_mngr(clan_arr, finished_file, max_complete, loci, hgLoci,
                     index, blast_ids, mingcfid, simfun):
    clan_arr = clan_arr.tocsr()
    # a matrix representing the cardinality of the intersection of hg_locs
    adj_arr = csr_matrix(triu(clan_arr @ clan_arr.transpose()))
#    adj_arr = clan_arr @ clan_arr.transpose()
    # less than two gene overlap is set to 0
    adj_arr.data[adj_arr.data < 2] = 0
    hits = np.split(adj_arr.indices, adj_arr.indptr)[1:-1]
    with open(f'{finished_file}.w', 'a') as out:
        for ti0, overlap_loci in enumerate(hits):
            i0 = ti0 + max_complete
#            try:
            loc0, hgL0 = loci[i0], hgLoci[i0]
#            except IndexError:
 #               continue
            sHGl0 = set([x for x in hgL0 if x is not None])
            for ti1 in overlap_loci[:]:
  #              if ti1 <= ti0: 
 #                   continue
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

def open_prev_res(finished_file, index):
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
                        max_complete = int(data[0][0]) - index
                        data = []
                except IndexError: # first line
                    pass
                data.append(d)
    shutil.move(f'{finished_file}.w1', f'{finished_file}.w')
    return max_complete

def rnd2_loc2loc_mngr(gcf_dir, loci, hg_loci,
                     hgx_dir, groups, minid, mingcfid, simfun, cpus = 1):

    finished_file = f'{gcf_dir}loci.adj.tmp'
    # is this clan already completed?
    if os.path.isfile(finished_file) or os.path.isfile(finished_file + '.r'):
        return
    # was this clan partially written?
    elif os.path.isfile(f'{finished_file}.w'):
        os.remove(f'{finished_file}.w')

    # only examine loci that have not been completed
    ghg2i = defaultdict(dict)
    blast_hash = defaultdict(lambda: defaultdict(list))
    for g,v in enumerate(groups):
        i0, i1 = v
        for i2 in range(i0, i1):
            locus = loci[i2]
            hgs = hg_loci[i2]
            for i4, hg in enumerate(hgs):
                if hg is not None:
                    blast_hash[hg][g].append(locus[i4])
                    if hg not in ghg2i[g]:
                        ghg2i[g][hg] = len(ghg2i[g])

    clan_arrs = []
    for g, v in enumerate(groups):
        i0, i1 = v
        clan_arr = lil_matrix((i1 - i0, len(ghg2i[g])), dtype=bool)
        for i2 in range(i0, i1):
            hgs = hg_loci[i2]
            for hg in hgs:
                if hg is not None:
                    clan_arr[i2 - i0, ghg2i[g][hg]] = True
        clan_arrs.append(clan_arr)


    blast_ids = defaultdict(dict)
    print('\t\t\t\tLoading alignment results', flush = True)
    if not os.path.isdir(f'{gcf_dir}algn/'):
        os.mkdir(f'{gcf_dir}algn/')
    jsons = collect_files(f'{gcf_dir}algn/', 'json')
    finished_hgs = set(int(os.path.basename(x).replace('.json', '')) \
                    for x in jsons)
    missing_hgs = set(blast_hash.keys()).difference(finished_hgs)
    if missing_hgs:
        blast_hash = {k: blast_hash[k] for k in list(missing_hgs)}
        with mp.get_context('forkserver').Pool(processes = cpus) as pool:
            pool.starmap(rnd2_gen_blastids_mp, 
                         ((hg, group_dict, hgx_dir, gcf_dir, minid) \
                           for hg, group_dict in blast_hash.items()))
            pool.close()
            pool.join()
    blast_ids = read_gcf_sim_jsons(f'{gcf_dir}algn/')

    # could be expedited by multiprocessing a writer and numba njit the math
    print('\t\t\t\tCalculating similarity', flush = True)
    if blast_ids:
        for g, v in enumerate(groups):
            i0, i1 = v
            clan_arr = clan_arrs[g]
            loc2loc_sim_mngr(clan_arr, finished_file, i0, loci, hg_loci,
                             0, blast_ids, mingcfid, simfun)
        # the file is now complete
        shutil.move(f'{gcf_dir}loci.adj.tmp.w', f'{gcf_dir}loci.adj.tmp')


def rnd1_loc2loc_mngr(
    db, clanI, loci, hgLoci, hg_dir, hgx_dir, gcf_dir, index,
    minid = 30, mingcfid = 0.15, simfun = overlap, min_overlap = 2,
    cpus = 1, Q = None
    ):

    finished_file = f'{gcf_dir}{clanI}.sadj.tmp'
    # is this clan already completed?
    if os.path.isfile(finished_file) or os.path.isfile(finished_file + '.r'):
        return
    # was this clan partially written?
    elif os.path.isfile(f'{finished_file}.w'):
        max_complete = open_prev_res(finished_file, index)
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
    clan_arr = lil_matrix((len(loci) - max_complete, len(blast_hash)), dtype=bool)
    for i, v in enumerate(loci[max_complete:]):
        hgs = hgLoci[i + max_complete]
        for i1, hg in enumerate(hgs):
            if hg is not None:
                clan_arr[i, hg2i[hg]] = True


    blast_ids = defaultdict(dict)
    if cpus < 2:
        for hg, genes in blast_hash.items():
            gene_len = len(genes)
            if gene_len > 1:
                gene_set = set(genes)
                algn_base = f'{hgx_dir}{hg}.out'
                try:
                    blast_result = parse_blast(algn_base, gene_len, gene_set,
                                                minid)
                    if blast_result is not None:
                        blast_ids[hg] = blast_result                
                except FileNotFoundError:
                    pass

    else:
        print('\t\t\t\t\tLoading alignment results', flush = True)
        if not os.path.isdir(f'{gcf_dir}algn/'):
            os.mkdir(f'{gcf_dir}algn/')
        jsons = collect_files(f'{gcf_dir}algn/', 'json')
        finished_hgs = set(int(os.path.basename(x).replace('.json', '')) \
                        for x in jsons)
        missing_hgs = set(blast_hash.keys()).difference(finished_hgs)
        if missing_hgs:
            blast_hash = {k: blast_hash[k] for k in list(missing_hgs)}
            with mp.get_context('forkserver').Pool(processes = cpus) as pool:
                pool.starmap(gen_blastids_mp, 
                             ((hg, genes, hgx_dir, gcf_dir, minid) \
                               for hg, genes in blast_hash.items() \
                               if len(genes) > 1))
                pool.close()
                pool.join()
        blast_ids = read_gcf_sim_jsons(f'{gcf_dir}algn/')

    # could be expedited by multiprocessing a writer and numba njit the math
    if blast_ids:
        if cpus > 1:
            print('\t\t\t\t\tCalculating locus similarity', flush = True)
            Q.put(None)
        loc2loc_sim_mngr(clan_arr, finished_file, max_complete, loci, hgLoci,
                         index, blast_ids, mingcfid, simfun)

#        else:
 #           gcf_sim_mp_mngr(hits, f'{finished_file}.w', 
  #                          max_complete, loci, hgLoci, index, 
   #                         blast_ids, gcf_dir, f'{gcf_dir}mtrx/', mingcfid, 
    #                        simfun, chunkscale = 100, cpus = cpus)
        # begin appending to temporary writing output
        # the file is now complete
        shutil.move(f'{gcf_dir}{clanI}.sadj.tmp.w', f'{gcf_dir}{clanI}.sadj.tmp')
   #     if cpus > 1: # leave up to zipping 
#            jsons = collect_files(f'{gcf_dir}algn/', 'json')
 #           for j_f in jsons:
  #              os.remove(j_f)


def merge_clan_loci(preclan_loci, gene2hg, min_sim, hgx2omes_path, hgx2dist_path, 
                    scaf2gene2i, min_heat_len = 25, min_branch_sim = 0.5, phylo = None):

    if min_heat_len > 3:
        vclan_loci = copy.deepcopy(preclan_loci)
        out_clan_loci = defaultdict(lambda: defaultdict(list))
        outclanI = defaultdict(lambda: defaultdict(list))
        for clan, scaf_dict in preclan_loci.items():
            for scaf, loci in scaf_dict.items():
                while loci: # exhaustively merge overlapping loci
                    # the locus for pairwise comparisons is the first
                    loc0, hgxs0, allHGxs, locI0 = set(loci[0][0]), loci[0][1], loci[0][2], loci[0][3]
                    locIntersect = None
                    # for all other loci
                    for i1, loc1d in enumerate(loci[1:]):
                        loc1, hgxs1, allHGxs1, locI1 = set(loc1d[0]), loc1d[1], loc1d[2], loc1d[3]
                        locIntersect = loc0.intersection(loc1)
                        # if any genes overlap between these loci
                        if locIntersect:
                            # grab the overlap and add it to the end of the loci list
                            newHGx = list(hgxs0)
                            newHGx.extend(list(hgxs1))
                            allHGxs.extend(allHGxs1)
                            loci.append([list(loc1.union(loc0)), newHGx, allHGxs, locI0 + locI1])
                            break
                    if locIntersect: # if there was overlap, delete the overlappers
                        del loci[0]
                        del loci[i1] # don't add 1 because we removed the first and i1 relative to 1
                    else: # no more overlap for this locus, add to the final output
                        out_clan_loci[scaf][clan].append(sorted(loc0))
                        outclanI[scaf][clan].append(locI0)
        #                out_clan_hgx[clan].append(tuple(sorted(set(hgxs0))))
                        del loci[0]
    
        # initiate the data structure for heat-based locus merging
        todel = defaultdict(list)
        n_loci = defaultdict(list)
        for clan, scaf_dict in out_clan_loci.items():
            for scaf, loci in scaf_dict.items():
                for i, loc in enumerate(loci):
                    if len(loc) >= min_heat_len:
                        for locI in outclanI[clan][scaf][i]:
                            n_loci[clan].append(vclan_loci[clan][scaf][locI])
                        todel[clan].append((i, scaf))
        for clan, td in todel.items():
            for i, scaf in reversed(td):
                del out_clan_loci[clan][scaf][i]
    
        if n_loci:
            heat_clan_loci = \
                merge_clan_loci_heat(n_loci, gene2hg, min_sim, hgx2omes_path,
                                     hgx2dist_path, scaf2gene2i, min_branch_sim,
                                     phylo)
            for clan, scaf_dict in heat_clan_loci.items():
                for scaf, loci in scaf_dict.items():
                    out_clan_loci[clan][scaf].extend(loci)
    else:
        out_clan_loci = merge_clan_loci_heat(preclan_loci, gene2hg, 
                                             min_sim, hgx2omes_path,
                                             hgx2dist_path, scaf2gene2i,
                                             min_branch_sim, phylo)

    clan_loci = {clan: list(chain(*list(scaf_dict.values()))) \
                 for clan, scaf_dict in out_clan_loci.items()}

    hg_loci = {}
    for clan, loci in clan_loci.items():
        hg_loci[clan] = []
        for locus in loci:
            hg_loci[clan].append([])
            for gene in locus:
                try:
                    hg_loci[clan][-1].append(gene2hg[gene])
                except KeyError:
                    hg_loci[clan][-1].append(None)

    return clan_loci, hg_loci


def merge_clan_loci_heat(preclan_loci, gene2hg, min_sim, 
                         hgx2omes_path, hgx2dist_path, scaf2gene2i,
                         min_branch_sim = 0.5, phylo = None):
    with open(hgx2dist_path, 'rb') as raw:
        hgx2dist = pickle.load(raw)
    with open(hgx2omes_path, 'rb') as raw:
        hgx2omes = pickle.load(raw)

    gene2hgx = defaultdict(list)
    out_clan_loci = defaultdict(lambda: defaultdict(list))
    for clan, scaf_dict in preclan_loci.items():
        for scaf, loci in scaf_dict.items():
            for i, loc_d in enumerate(loci):
                loc, hgx = loc_d[0], loc_d[1]
                loci[i].append(hgx2dist[hgx])
                loci[i].append(hgx2omes[hgx])
 #               loci[i].append([])
                [gene2hgx[x].append(hgx) for x in loc]
            # sort the scaffolds loci from highest to lowest TMD
            loci = sorted(loci, key = lambda x: x[-2], reverse = True)
  #          check_loci = loci
    #        merge, check_loci = True, loci
            # placeholder, currently we don't want to merge this way
            # in the future, this should be based on the sorensen/jaccard
            # similarity of the phylogenetic distributions of the traits
            # e.g. overlapping ome TMD / overall ome TMD
        #    while merge: # need to rerun to verify
         #       loci = copy.deepcopy(check_loci)
          #      check_loci, merge = [], False
           #     while loci: # exhaustively merge overlapping loci
            #        # the locus for pairwise comparisons is the first
             #       loc0, hgxs0, allHGxs, dist0 = set(loci[0][0]), loci[0][1], loci[0][2], loci[0][-1]
              #      locIntersect = None
               #     # for all other loci
                #    for i1, loc1d in enumerate(loci[1:]):
                 #       loc1, hgxs1, allHGxs1, dist1 = set(loc1d[0]), loc1d[1], loc1d[2], loc1d[-1]
                  #      locIntersect = loc0.intersection(loc1)
                        # if any genes overlap between these loci
                   #     if locIntersect:
                            # are they similar enough?
                    #        if dist1/dist0 >= min_sim:
    #                            if sml_up_bound/lrg_low_bound >= min_sim:
                                    # grab the overlap and add it to the end of the loci list
                     #           newHGx = list(hgxs0)
                      #          newHGx.extend(list(hgxs1))
                       #         allHGxs.extend(allHGxs1)
     #                           newdists = dists0 + dists1
                        #        loci.insert(0, [list(loc1.union(loc0)), newHGx, allHGxs, dist0])
                         #       break
                          #  else:
                           #     locIntersect = None
#                    if locIntersect: # if there was overlap, delete the overlappers
 #                       merge = True
  #                      del loci[1]
   #                     del loci[i1 + 1] # don't add 1 because we removed the first and i1 relative to 1
    #                else:
     #                   check_loci.append(loci[0])
      #                  del loci[0]
    
            # add an empty list to each locus to populate with merge indices
#            loci = [x + [[]] for x in \
 #                   sorted(check_loci, key = lambda x: x[-1], reverse = True)]
            # give the remaining overlapping genes to the highest distance group
            # will need to add the new HGxs in the end
            final_loci = []
            while loci:
                locd0 = loci[0]
                loc0 = set(locd0[0])
                if len(loc0) == 1:
                    break
#                overi, check = 1, True
#                while check:
 #                   check = False
                overi = 1
                for i1, locd1 in enumerate(loci[1:]):
                    loc1 = set(locd1[0])
                    locIntersection = loc0.intersection(loc1)
                    if locIntersection:
                        # index from the back because some may have hgx still
                        dist1, omes1 = locd1[-2], locd1[-1]
#                        check = True
 #                       overi += i1
    # currently this doesnt check if there are others of the same hg that are removed in the
    # overlap - ideally, this should check if there are other members of the same hg before
    # removing them IF they are not part of the overlap
                        # the loci are sorted from largest to lowest distribution, so
                        # we award the overlap to the largest TMD
                        new_locs = [sorted(loc1.difference(loc0), 
                                    key = lambda x: scaf2gene2i[scaf][x])]
                        if not new_locs[0]:
                            new_locs = []
                        else:
                            break_i = True
                            while break_i:
                                new_loc = new_locs[0]
                                break_i = False
                                i_loc = [scaf2gene2i[scaf][x] for x in new_loc]
                                for i2, v in enumerate(i_loc[:-1]):
                                    if i2:
                                        if not i_loc[i2+1] - v == 1:
                                            break_i = i2 + 1
                                            break
                                if break_i:
                                    n_loc0 = new_loc[:break_i]
                                    n_loc1 = new_loc[break_i:]
                                    if len(n_loc0) > 1:
                                        new_locs.append(n_loc0)
                                    else:
                                        loci.append([n_loc0, dist1, omes1])
                                    if len(n_loc1) > 1:
                                        new_locs.append(n_loc1)
                                    else:
                                        loci.append([n_loc1, dist1, omes1])
                                    new_locs.pop(0)
                        # delete the initial version with overlap
                        loci.pop(i1+overi)
                        for new_loc in new_locs:
                            # insert the new loci
                            loci.insert(overi+i1, [new_loc, dist1, omes1])
                        # add one because these aren't overlapping
                        overi += len(new_locs) - 1

#                        break
                final_loci.append(locd0)
                loci.pop(0)
    

#            for i0, locd0 in enumerate(loci):
 #               loc0, hgx0, ahgx0 = set(locd0[0]), locd0[1], locd0[2]
  #              if len(loc0) == 1:
   #                 break
                # preserve a slice to modify
    #            t_loci = loci[i0+1:] 
     #           todel = []
      #          for i1, locd1 in enumerate(t_loci):
       #             loc1, hgx1, ahgx1 = set(locd1[0]), locd1[1], locd1[2]
        #            dist1, omes1 = locd1[-2], locd1[-1]
         #           locIntersection = loc0.intersection(loc1)
          #          if locIntersection:
    # currently this doesnt check if there are others of the same hg that are removed in the
    # overlap - ideally, this should check if there are other members of the same hg before
    # removing them IF they are not part of the overlap
           #             hg_intersection = set(gene2hg[gene] for gene in list(locIntersection) \
            #                                  if gene in gene2hg)
                        # the loci are sorted from largest to lowest distribution, so
                        # we award the overlap to the largest TMD
             #           new_locs = [sorted(loc1.difference(loc0), 
              #                      key = lambda x: scaf2gene2i[scaf][x])]
               #         if not new_locs[0]:
                #            new_locs = []
                 #       else:
                  #          break_i = True
                   #         while break_i:
                    #            new_loc = new_locs[0]
                     #           break_i = None
                      #          i_loc = [scaf2gene2i[scaf][x] for x in new_loc]
                       #         for i, v in enumerate(i_loc[:-1]):
                        #            if i:
                         #               if not i_loc[i+1] - v == 1:
                          #                  break_i = i + 1
                           #                 break
                            #    if break_i:
                             #       n_loc0 = new_loc[:break_i]
                              #      n_loc1 = new_loc[break_i:]
                               #     if len(n_loc0) > 1:
                                #        new_locs.append(n_loc0)
                                 #   else:
                                  #      t_loci.append([[n_loc0[0]], None, None, dist1, omes1])
                                   # if len(n_loc1) > 1:
                                    #    new_locs.append(n_loc1)
                                   # else:
                                    #    t_loci.append([[n_loc1[0]], None, None, dist1, omes1])
                         #           new_locs.pop(0)

  #                      for new_loc in new_locs:
 #                           t_loci.insert(i1+1, [new_loc, hgx1, 
   #                                               [tuple(sorted(set(x).difference(hg_intersection))) \
    #                                               for x in ahgx1],
     #                                             dist1, omes1])
                        # remove the original locus
      #                  todel.append(i1)
       #         for i in reversed(todel):
        #            t_loci.pop(i)
         #       loci[i0+1:] = t_loci
                            # collect the index of the merge for later
#                            loci[i1 + 1 + i0][-1].append(i0)
            # sort the singletons by TMD, highest to lowest
            singletons_prep = sorted([x for x in loci],
                                key = lambda y: y[-2], reverse = True)
            singletons, s_set = [], set()
            for s in singletons_prep:
                if s[0][0] not in s_set:
                    singletons.append(s)
                    s_set.add(s[0][0])
            singletons = sorted(singletons, key = lambda y: scaf2gene2i[scaf][y[0][0]])

            singleton_iz = [scaf2gene2i[scaf][x[0][0]] for x in singletons]
            continuities = []
            for i, v in enumerate(singleton_iz[:-1]):
                if singleton_iz[i + 1] - v == 1:
                    continuities.append(i)

            merges = []
            for i in continuities:
                s0, s1 = singletons[i], singletons[i+1]
                omes0 = s0[-1]
                omes1 = s1[-1]
                branch_sim = treecalcs.calc_branch_sim(phylo, omes0, omes1)
                if branch_sim > min_branch_sim:
                    merges.append(set(s0[0] + s1[0]))

            for i0, loc0 in enumerate(merges[:-1]):
                loc1 = merges[i0+1]
                locIntersection = loc0.intersection(loc1)
                if locIntersection:
                    merges[i0 + 1] = list(set(loc0).union(set(loc1)))
                merges[i0] = None

            for i, loc in enumerate(merges):
                if loc:
                    final_loci.append([loc])

#            final_loci = [x for x in loci if len(x[0]) > 1]
            for locI, locx in enumerate(final_loci):   
#            for locI in range(len(loci) - 1, -1, -1):
#                locx = loci[locI]
                loc = locx[0]
                loc_len = len(loc)
            #    if loc_len > 1:
                out_clan_loci[clan][scaf].append(sorted(loc))
                # award single gene locus to what it is derived from merging into
            #    elif loc_len == 1:
             #       mergeI = locx[-1][0]
                    # if this isnt a nested merge
              #      if out_clan_loci[clan][scaf][mergeI]:
               #         out_clan_loci[clan][scaf][mergeI].extend(loc)
                #        out_clan_hgx[clan][scaf][mergeI] = \
                  #          tuple(sorted(set(list(out_clan_hgx[clan][scaf][mergeI]) \
                 #                              + ahgx)))
                   # out_clan_loci[clan][scaf].append(None)
                    #out_clan_hgx[clan][scaf].append(None)
              #  else:
               #     out_clan_loci[clan][scaf].append(None)
                #    out_clan_hgx[clan][scaf].append(None)

    clan_loci = {} 
    for clan, scaf_dict in out_clan_loci.items():
        clan_loci[clan] = {}
        for scaf, loci in scaf_dict.items():
            clan_loci[clan][scaf] = []
            for i, loc in enumerate(loci):
                if loc:
                    clan_loci[clan][scaf].append(loc)

    return clan_loci #, out_hg_loci



def hash_clan_loci(ome, gff_path, ome_sig_clus, gene2hg, clusplusminus,
                   min_merge_perc = 0, hgx2dist_path = None, hgx2omes_path = None,
                   min_heat_len = 25, min_branch_sim = 0.5, phylo = None):
    """sliding window of size clusplusminus, identify sequences that may have
    hgs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus"""

    gff_list, protoclus, clus_out = gff2list(gff_path), {}, []
    cds_dict = input_parsing.compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))
    scaf2gene2i = {}
    for scaf, prots in cds_dict.items():
        scaf2gene2i[scaf] = {}
        for i, gene in enumerate(prots):
            scaf2gene2i[scaf][gene] = i

    preclan_loci, loc_index = defaultdict(lambda: defaultdict(list)), defaultdict(int)
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
                    preclan_loci[clan][scaf].append([
                         locus[start:end], hgx, [hgx], [loc_index[clan]]
                         ])
                    loc_index[clan] += 1

    out_clan_loci, out_hg_loci = merge_clan_loci(preclan_loci, gene2hg,
                                                             min_merge_perc, hgx2omes_path,
                                                             hgx2dist_path,
                                                             scaf2gene2i, min_heat_len,
                                                             min_branch_sim, phylo)
    return ome, out_clan_loci, out_hg_loci


def read_to_write(out_h, in_adj, min_gcf):
    with open(in_adj, 'r') as raw:
        for line in raw:
            d = line.rstrip().split()
            gcf_id = round(float(d[2]) * 100)
            if gcf_id >= min_gcf:
                out_h.write(line.rstrip() + '\n')
    out_h.flush()


def write_adj_matrix(out_file, gcf_dir, min_loc_id):
    comp_gcf_id = round(float(min_loc_id) * 100)
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


def read_clus(gcf_dir, mci2pre):
    t_gcfs = defaultdict(list)
    with open(gcf_dir + 'loci.clus', 'r') as raw:
        for line in raw: # loci indices
            gcf, mci = [int(x) for x in line.rstrip().split()]
            t_gcfs[gcf].append(mci2pre[mci])

    return [v for k, v in sorted(t_gcfs.items(), key = lambda x: x[0]) \
            if len(v) > 1]
  #          indices = [int(x) for x in line.rstrip().split()]
            # singletons won't pass thresholds, disregard them
    #        if len(indices) > 1:
     #           t_gcfs.append(indices)



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

    mci2pre = read_mci_rows(gcf_dir + 'mcl_rows.tsv')
    t_gcfs = read_clus(gcf_dir, mci2pre)

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
                hgxs = hgxXloci[locI]
                # really should be done above to make hgxs formatted right
                omeI = ome2i[loc[0][:loc[0].find('_')]]
#                [gcfs[-1][hgx].append(omeI) for hgx in hgxs]
                gcfs[-1].append((loc, locI))
                gcfOme_list.append(omeI)
                gcf_hgxs[-1].extend(hgxs)
            if len(set(gcfOme_list)) > min_omes: # need more than 1
#                gcfs[-1] = {k: sorted(set(v)) for k,v in gcfs[-1].items()}
                if gcfs[-1]:
                    gcfs[-1] = [tuple(x) for x in gcfs[-1]]
                    gcf_hgxs[-1] = tuple(sorted(set(gcf_hgxs[-1])))
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

def parse_and_sort_gff(gff_path):
    # order the loci in each scaffold
    cds_dict = input_parsing.compileCDS(gff2list(gff_path),
                                        os.path.basename(gff_path).replace('.gff3',''))
    scaf2gene2i = {}
    for scaf, prots in cds_dict.items():
        scaf2gene2i[scaf] = {}
        for i, gene in enumerate(prots):
            scaf2gene2i[scaf][gene] = i
   
    return scaf2gene2i, cds_dict

def hash_ome_lg_loc(gff_path, locs_y_lg):
    scaf2gene2i, cds_dict = parse_and_sort_gff(gff_path)
    scaf2i2gene = {}
    for scaf, gene2i in scaf2gene2i.items():
        scaf2i2gene[scaf] = {v: k for k, v in gene2i.items()}
    scaf2locs = output_data.clus2scaf(cds_dict, locs_y_lg)

    gcf2scaf2locs = defaultdict(lambda: defaultdict(list))
    for scaf, locs in scaf2locs.items():
        for loc, gcf in locs:
            gcf2scaf2locs[gcf][scaf].append(loc)
    
    gcf2locs = {}
    for gcf, scaf_dict in gcf2scaf2locs.items():
        gcf2locs[gcf] = []
        for scaf, locs in scaf_dict.items():
            while len(locs) > 1:
                loc0 = locs[0]
                loc1 = locs[1]
                top_i = scaf2gene2i[scaf][loc0[-1]]
                bot_i = scaf2gene2i[scaf][loc1[0]]
                if bot_i - top_i == 1:
                    locs[0].extend(locs[1])
                    locs.pop(1)
                elif bot_i - top_i == 2: # one intervening gene, so award it
                    locs[0].extend([scaf2i2gene[scaf][top_i + 1]] + locs[1])
                    locs.pop(1)
                else:
                    gcf2locs[gcf].append(loc0)
                    locs.pop(0)
            gcf2locs[gcf].append(locs[0])

    return gcf2locs


def lg_loc_mngr(db, lgs, gene2hg, cpus = 1):
    omes2loci = defaultdict(list)
    for lg, locs in lgs.items():
        for loc, locI in locs:
            ome = loc[0][:loc[0].find('_')]
            omes2loci[ome].append((loc, lg))

    with mp.get_context('forkserver').Pool(processes = cpus) as pool:
        lg_loci = pool.starmap(hash_ome_lg_loc,
                               tqdm(((db[ome]['gff3'], omes2loci[ome]) \
                                     for ome in db.keys()), 
                                    total = len(db)))

    lg2locs = defaultdict(list)
    lg2hg_locs = defaultdict(list)
    lg2hgxs_locs = defaultdict(list)
    for res in lg_loci:
        for lg, locs in res.items():
            for loc in locs:
                hg_loc = []
                for gene in loc:
                    try:
                        hg_loc.append(gene2hg[gene])
                    except KeyError:
                        hg_loc.append(None)
                lg2locs[lg].append(loc)
                lg2hg_locs[lg].append(hg_loc)

    return lg2locs, lg2hg_locs    
#            db, clanI, locs, clan_hg_loci[clanI],
 #           hg_dir, hgx_dir, gcf_dir, index, minid, min_loc_id,
  #          simfun, min_overlap

def refine_group(db, ome2i, group_hg_loci, group_loci,
                 hg_dir, hgx_dir, wrk_dir, grp_dir, minid, min_loc_id,
                 simfun, inflation, tune, min_omes, cpus = 1, rnd = 1):

    m, adj_mtr = mp.Manager(), grp_dir + 'loci.adj'
    if cpus < 2:
        cpus = 2

#    min_overlap = round((min_loc_id * 2) - 0.5)
    min_overlap = 2 # need at least two genes of overlap to call it a cluster
    index, cmds = 0, []

    protog2h = {gI: Counter(chain(*hg_locs)) \
                for gI, hg_locs in enumerate(group_hg_loci)}
    group2hgx = {gI: set([k for k, v in h.items() if v > 1]) \
                 for gI, h in protog2h.items()}

    if rnd != 2:
        loci, hgxXloci, loc2group = {}, {}, {}
        for gI, locs in enumerate(group_loci):
            if not os.path.isfile(f'{grp_dir}{gI}.sadj.tmp.r'):
                cmds.append([
                    db, gI, locs, group_hg_loci[gI],
                    hg_dir, hgx_dir, grp_dir, index, minid, min_loc_id,
                    simfun, min_overlap
                    ])
            for i, locus in enumerate(locs):
                hg_loc = set(group_hg_loci[gI][i])
                loci[index] = group_loci[gI][i]
                hgxXloci[index] = tuple(sorted([x for x in \
                                       list(hg_loc.intersection(group2hgx[gI])) \
                                       if x is not None]))
                loc2group[index] = gI
                index += 1
        if not os.path.isfile(adj_mtr):
            W = mp.Process(target = write_adj_matrix, 
                           args = (adj_mtr, grp_dir, min_loc_id))
            W.start()
    #        print(f'\t\t\t{max([len(x) for x in )} loci in clan 0', flush = True)
            procs = []

            if cmds:
                if len(cmds[0][2]) > 30000 and cmds[0][1] == 0:
                    print(f'\t\t\t\tRunning group 0', flush = True)
                    Q = mp.Manager().Queue()
                    procs = [mp.Process(target=rnd1_loc2loc_mngr,
                                    args=(*cmds[0], cpus - 1, Q))]
                    procs[0].start()
                    Q.get()
                    print('\t\t\t\tRunning all groups', flush = True)
                    del cmds[0]
                with mp.get_context('forkserver').Pool(processes = cpus - 1) as pool:
                    pool.starmap(rnd1_loc2loc_mngr, tqdm(cmds, total = len(cmds), miniters = 1))
                    pool.close()
                    pool.join()
                if procs:
                    procs[0].join()
                    procs[0].close()
            with open(f'{grp_dir}done.sadj.tmp', 'w') as out:
                pass
            W.join()
    else:
        index, hgxXloci, loc2group, loci, hg_loci, groups = 0, {}, {}, {}, {}, []
        for gI, locs in enumerate(group_loci):
            groups.append([index, -1])
            for i, locus in enumerate(locs):
                hgs = group_hg_loci[gI][i]
                hg_loci[index] = hgs
                hg_loc = set(hgs)
                loci[index] = group_loci[gI][i]
                hgxXloci[index] = tuple(sorted([x for x in \
                                       list(hg_loc.intersection(group2hgx[gI])) \
                                       if x is not None]))
                loc2group[index] = gI
                index += 1
            groups[-1][1] = index
        if not os.path.isfile(f'{grp_dir}loci.adj.tmp') \
            and not os.path.isfile(f'{grp_dir}loci.adj'):
            rnd2_loc2loc_mngr(grp_dir, loci, hg_loci, hgx_dir, 
                              groups, minid, min_loc_id, simfun, cpus)
            os.rename(f'{grp_dir}loci.adj.tmp', f'{grp_dir}loci.adj')

    satisfied = False
    if tune:
        inflation = 1.1
    while not satisfied:
        if tune:
            inflation += 0.1
            if inflation > 3:
                break
        satisfied, grps, grp_hgxs, grp_omes, lilg2bigg = mcl2gcfs(grp_dir, loci,
                                                           hgxXloci, ome2i, 
                                                           inflation, loc2group,
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

    return grps, grp_hgxs, grp_omes, lilg2bigg, hgxXloci


def comp_scaf2hgs(cds_dict, gene2hg):
    scaf2hgs = {}
    for scaf, prots in cds_dict.items():
        scaf2hgs[scaf] = []
        for prot in prots:
            try:
                scaf2hgs[scaf].append(gene2hg[prot])
            except KeyError:
                scaf2hgs[scaf].append(None)
    return scaf2hgs


def extend_loci(ome, gff_path, loc2need, gene2hg, forbid_genes, 
                ext_len = 1, min_new_cont = 2, announcements = {87332}):
    forbidden_genes = set(forbid_genes)
    scaf2gene2i, cds_dict = parse_and_sort_gff(gff_path)
    
    gene2scaf, index_map = {}, {}
    for scaf, genes in cds_dict.items():
        index_map[scaf] = {v: i for i, v in enumerate(genes)}
        for gene in genes:
            gene2scaf[gene] = scaf

    # currently only built for an extension length of 1, will need modification
    # to accomodate greater lengths, e.g. determining how to merge an extension
    # if only the nearest-to-the-cluster entry of the extension meets missing_hg
    mod_locs, edge_locs = {}, {}
    for locI, d in loc2need.items():
        # opt to declare contig edges based on the starting locus, not the end
        # this is a decision that has arguments for/against both mechanisms
        # I'm setting it on the start to be conservative about this module
        contig_edge = False
        # keep track of modified loci, or contig edge loci for future mod
        modified = False
        loc, gcf, miss_hg_tup = d

        # announce this GCF for auditing
#        if gcf in announcements:
 #           announce = True
  #      else:
   #         announce = False

        miss_hg = set(miss_hg_tup)
        gene0 = loc[0]
        scaf = gene2scaf[gene0]
        scaf_len = len(cds_dict[scaf])
        ord_loc = sorted(loc, key = lambda pair: index_map[scaf][pair])
        # identify whether the top/bottom can be extended
        bot_i, top_i = scaf2gene2i[scaf][ord_loc[0]], scaf2gene2i[scaf][ord_loc[-1]]
        if bot_i:
            bot_i -= ext_len
            bot_ex = cds_dict[scaf][bot_i]
            # do NOT interfere with existing clusters
            if bot_ex in forbidden_genes:
                bot_ex = False
        else:
            contig_edge = True
            bot_ex = False
        if top_i < scaf_len - 1:
            top_i += ext_len
            top_ex = cds_dict[scaf][top_i]
            if top_ex in forbidden_genes:
                top_ex = False
        else:
            contig_edge = True
            top_ex = False

        # while there is a true extension possibility
        while top_ex or bot_ex:
            todel = set()
            if top_ex:
                try:
                    top_hg = gene2hg[top_ex]
                except KeyError:
                    top_hg = None
                if top_hg in miss_hg:
                    modified = True
                    ord_loc.append(top_ex)
                    todel.add(top_hg)
                    if top_i < scaf_len - 1:
                        top_i += ext_len
                        top_ex = cds_dict[scaf][top_i]
                        if top_ex in forbidden_genes:
                            top_ex = False
                else:
                    top_ex = False
            if bot_ex:
                try:
                    bot_hg = gene2hg[bot_ex]
                except KeyError:
                    bot_hg = None
                if bot_hg in miss_hg:
                    modified = True
                    ord_loc.insert(0, bot_ex)
                    todel.add(bot_hg)
                    if bot_i:
                        bot_i -= ext_len
                        bot_ex = cds_dict[scaf][bot_i]
                        if bot_ex in forbidden_genes:
                            bot_ex = False
                else:
                    bot_ex = False
            # dont edit the miss_hg until the end to not bias for/against top/bot    
            miss_hg = miss_hg.difference(todel)
        if modified:
            mod_locs[locI] = (ord_loc, gcf)
        
        # its risky to extend to different contigs based on a single HG, so we
        # should only do it if there's evidence that exceeds min_new_cont (2)
        if contig_edge and len(miss_hg) >= min_new_cont:
            edge_locs[locI] = (ord_loc, gcf, miss_hg)
        
    # if there are loci at the edge of a contig
    to_add_locs = []
    if edge_locs:
        scaf2hgs = comp_scaf2hgs(cds_dict, gene2hg)
        # we are only looking for extensions of min_new_cont length (2)
        # anything greater than this should have been picked up by a different
        # mechanism
        for locI, d in edge_locs.items():
            loc, gcf, miss_hg = d

            gene0 = loc[0]
            loc_scaf = gene2scaf[gene0]
            found = False
            # it will extend based on the first scaffold that meets the criteria
            for scaf, hg_list in scaf2hgs.items():
                # ignore self-scaffold check
                if loc_scaf == scaf:
                    continue
                bot_prot = cds_dict[scaf][:min_new_cont]
                top_prot = cds_dict[scaf][-min_new_cont:]
                # ignore other cluster entries
                if all(x not in forbidden_genes for x in bot_prot):
                    # this approach allows the rare situation for both sides of 
                    # the contig to meet the criteria before breaking
                    bot_set = hg_list[:min_new_cont]
                    if len(set(bot_set).intersection(miss_hg)) >= min_new_cont:
                        to_add_locs.append((bot_prot, gcf))
                        found = True
                if all(x not in forbidden_genes for x in top_prot):
                    top_set = hg_list[-min_new_cont:]
                    if len(set(top_set).intersection(miss_hg)) >= min_new_cont:
                        if top_set != bot_set:
                            to_add_locs.append((top_prot, gcf))
                            found = True
                if found:
                    break

    # ideally would like to propagate contig edge to downstream gbk because
    # antiSMASH does it and software like BiG-SCAPE use that data
    return ome, mod_locs, to_add_locs


def extend_loc_mngr(db, gcfs, gcf_hgxs, ome2gene2hg, hgxXloci, gcf_dir,
                    ext_len = 1, min_new_cont = 2, cpus = 1):
    ome2genes = defaultdict(list)
    ome2loc2need = defaultdict(dict)
    for gcf, locs in gcfs.items():
        gcf_hg_set = set(gcf_hgxs[gcf])
        for loc, locI in locs:
            ome = loc[0][:loc[0].find('_')]
            hg_loc = hgxXloci[locI]
            miss_hg = gcf_hg_set.difference(set(hg_loc))
            if miss_hg:
                ome2loc2need[ome][locI] = (loc, gcf, tuple(miss_hg))
            ome2genes[ome].extend(loc)
   
    with mp.Pool(processes = cpus) as pool:
        modifications = pool.starmap(extend_loci, 
                                     ((ome, db[ome]['gff3'], loc2need,
                                       ome2gene2hg[ome], ome2genes[ome],
                                       ext_len, min_new_cont) \
                                      for ome, loc2need in ome2loc2need.items()))

    with open(f'{gcf_dir}edge_clusters.tsv', 'w') as edge_out:
        edge_out.write('#gcf\tlocus\n')
        with open(f'{gcf_dir}extension_clusters.tsv', 'w') as extension_out:
            extension_out.write('#gcf\toriginal_locus\tupdate_locus\n')
            mod_count, cont_count = 0, 0
            for ome, mod_locs, to_add_locs in modifications:
                mod_count += len(mod_locs)
                for locI, d in mod_locs.items():
                    loc, gcf = d
                    ome2genes[ome].extend(loc)
                    for i, v in enumerate(gcfs[gcf]):
                        t_loc, t_locI = v
                        if t_locI == locI:
                            extension_out.write(f'{gcf}\t{",".join(t_loc)}\t{",".join(loc)}\n')
                            gcfs[gcf][i] = (loc, None)
                ome2genes[ome] = set(ome2genes[ome])
                for loc, gcf in to_add_locs:
                    # this awards the new loci to the first entry, which doesn't
                    # have an order that is necessarily meaningful
                    if all(x not in ome2genes[ome] for x in loc):
                        cont_count += 1
                        gcfs[gcf].append((loc, None))
                        edge_out.write(f'{gcf}\t{",".join(loc)}\n')
                        ome2genes[ome] = ome2genes[ome].union(set(loc))

        
    # we aren't modifying the hgs of each gcf so can ignore them
    return {gcf: tuple([loc[0] for loc in locs]) for gcf, locs in gcfs.items()}, \
           mod_count, cont_count


def classify_gcfs(
    hgx2loc, db, gene2hg, i2hgx, hgx2i,
    phylo, ome2i, hgx2omes, hg_dir, hgx_dir,
    wrk_dir, ome2partition, bord_scores_list,
    hg2gene, omes2dist = {}, clusplusminus = 3,
    inflation = 1.5, min_hgx_overlap = 1, min_omes = 2, 
    minid = 30, cpus = 1, algorithm = 'diamond',
    min_loc_id = 0.3, simfun = overlap, printexit = False,
    tune = False, dist_func = treecalcs.calc_tmd,
    uniq_sp = False, min_merge_perc = 0, algn_sens = '',
    skipalgn = False, fallback = False, min_heat_len = 25,
    min_branch_sim = 0.8
    ):

    i2ome = {v: k for k, v in ome2i.items()}

    if not checkdir(hgx_dir, unzip = True, rm = True):
        os.mkdir(hgx_dir)

    groupI = wrk_dir + 'gene2hgx.pickle'
    groupII = wrk_dir + 'hgcs.pickle'

    # make an ome-centric hash for less MP overhead
    ome2gene2hg = defaultdict(dict)
    for gene, hg in gene2hg.items():
        ome = gene[:gene.find('_')]
        ome2gene2hg[ome][gene] = hg

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
                ome, hgx_genes[ome], ome2gene2hg[ome], clusplusminus
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
        g = Graph(directed = False)
        idx = matrix.nonzero()
        g.add_edge_list(np.transpose(matrix.nonzero()))
        components = label_components(g)[0]
        subgraph_is = sorted(set(components.a))
#        network = nx.from_scipy_sparse_matrix(matrix)
 #       subgraphs = [network.subgraph(c) for c in nx.connected_components(network)]
        # use .copy() after iterator to copy
        clans, clanOmes, clanHGxs = [], [], []
#        for sg in subgraphs:
        for sgi in subgraph_is:
            sg = GraphView(g, vfilt = components.a == sgi)
            clans.append({})
            clanOmes.append([])
            clanHGxs.append([])
#            for id_ in list(sg.nodes):
            for id_ in sg.vertices():
                hgx = i2hgx[int(id_)]
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
        print('\t\tPreparing locus extraction', flush = True)
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
        with mp.get_context('forkserver').Pool(processes = cpus) as pool:
            hash_res = pool.starmap(hash_clan_loci, 
                                    tqdm(((ome, db[ome]['gff3'], clus2extract,
                                          ome2gene2hg[ome], clusplusminus, min_merge_perc,
                                          f'{wrk_dir}hgx2dist.pickle', 
                                          f'{wrk_dir}hgx2omes.pickle', min_heat_len, min_branch_sim,
                                          phylo) \
                                          for ome, clus2extract \
                                          in ome2clus2extract.items()),
                                         total = len(ome2clus2extract)))
            pool.close()
            pool.join()

        print('\t\tWriting clan loci', flush = True)
        clan_loci = defaultdict(list)
#        clan_hgx4gcfs = defaultdict(list)
        clan_hg_loci = defaultdict(list)
        for ome, out_clan_loci, out_clan_hg_loci in hash_res:
            for clan in out_clan_loci:
                clan_loci[clan].extend(out_clan_loci[clan])
 #               clan_hgx4gcfs[clan].extend(out_clan_hgx[clan])
                clan_hg_loci[clan].extend(out_clan_hg_loci[clan])

        write_json(clan_loci, clanFile + 'json.gz')
  #      write_json(clan_hgx4gcfs, clanFile + 'hgx.json.gz')
        write_json(clan_hg_loci, clanFile + 'hg.json.gz')
    else:
        print('\t\tLoading clan loci', flush = True)
        clan_loci = {int(k): tuple(v) for k,v in sorted(
                read_json(clanFile + 'json.gz').items(), 
                key = lambda x: len(x[1]), reverse = True)}
   #     clan_hgx4gcfs = {int(k): tuple([tuple(i) for i in v]) \
    #                    for k,v in read_json(clanFile + 'hgx.json.gz').items()}
        clan_hg_loci = {int(k): tuple(v) \
                      for k, v in read_json(clanFile + 'hg.json.gz').items()}

    # dict(clan_loci) = {int(clanI): ((ome0_gene0, ome0_gene1,)...,)}
    # dict(clan_hgx4gcfs) = {int(clanI): ((og0, og1, ogn,)...,)}
    # dict(clan_hg_loci) = {int(clanI): ((ome0_gene0_og, ome0_gene1_og,)...,)}

    # need to move up to before writing
    clan_loci = {k: v for k,v in sorted(clan_loci.items(), key = lambda x: len(x[1]), reverse = True)}
    #clan_hgx4gcfs = [clan_hgx4gcfs[k] for k in clan_loci]
    clan_hg_loci = [clan_hg_loci[k] for k in clan_loci]
    clan_loci = [v for v in clan_loci.values()]

    gcf_dir = wrk_dir + 'gcf/' 
    if not os.path.isdir(gcf_dir):
        os.mkdir(gcf_dir)
    lg_dir = gcf_dir + 'lg/'
    if not os.path.isdir(lg_dir):
        os.mkdir(lg_dir)

    print('\t\tOutputting HG fastas', flush = True)
#    hgs_prep = set(chain(*list(chain(*list(clan_hg_loci.values())))))
    hgs_prep = list(chain(*list(chain(*clan_hg_loci))))
    hg_count = Counter(hgs_prep)
    if None in hg_count:
        del hg_count[None]
    hgs = sorted([x for x, v in hg_count.items() \
                  if v > 1])
    hg_dir = input_parsing.hg_fa_mngr(wrk_dir, None, 
                             hgs, db, hg2gene, cpus = cpus)
    evo_conco.run_blast(hgs, db, hg_dir, hgx_dir,
                        algorithm = algorithm, printexit = printexit,
                        sensitivity = algn_sens, hg2gene = hg2gene,
                        skipalgn = skipalgn, #minid = minid, let minid be at locus step
                        fallback = fallback, cpus = cpus)

    print('\t\tCalling LGs', flush = True)
    clan_lens = [len(v) for v in clan_loci]
#    print('\t\t\t' + str(sum([len(v) for v in list(clan_loci.values())])) \
    print(f'\t\t\t{sum(clan_lens)}' \
        + ' loci', flush = True)
    print(f'\t\t\t{max(clan_lens)}' \
        + ' loci in clan 0', flush = True)
    lgs, lg_hgxs, lg_omes, lg2clan, hgxXloci = refine_group(
                 db, ome2i, clan_hg_loci, clan_loci,
                 hg_dir, hgx_dir, wrk_dir, lg_dir, minid, min_loc_id,
                 simfun, inflation, tune, min_omes, cpus = cpus)
    print('\t\t\t' + str(len(lgs)) + ' LGs', flush = True)

    print('\t\tMerging loci', flush = True)
    lg2locs, lg2hg_locs = lg_loc_mngr(db, lgs, gene2hg, cpus)
    lg_locs = {k: v for k,v in sorted(lg2locs.items(), key = lambda x: len(x[1]), reverse = True)}
    lg_hg_loci = [lg2hg_locs[k] for k in lg_locs]
    lg_locs = [v for v in lg_locs.values()]

    print('\t\tCalling GCFs', flush = True)
    lg_lens = [len(v) for v in lg_locs]
#    print('\t\t\t' + str(sum([len(v) for v in list(clan_loci.values())])) \
    print(f'\t\t\t{sum(lg_lens)}' \
        + ' loci', flush = True)
    print(f'\t\t\t{max(lg_lens)}' \
        + ' loci in LG 0', flush = True)
    gcfs, gcf_hgxs, gcf_omes, gcf2lg, hgxXloci = refine_group(
                 db, ome2i, lg_hg_loci, lg_locs,
                 hg_dir, hgx_dir, wrk_dir, gcf_dir, minid, min_loc_id,
                 sorensen, inflation, tune, min_omes, cpus = cpus, rnd = 2)
    print('\t\t\t' + str(len(gcfs)) + ' GCFs', flush = True)

    print('\tExtending incomplete clusters and contig edge clusters', flush = True)
    gcfs, mod_count, cont_count = extend_loc_mngr(db, gcfs, gcf_hgxs, ome2gene2hg, 
                                                  hgxXloci, gcf_dir, ext_len = 1,
                                                  min_new_cont = 2, cpus = cpus)
    print(f'\t\t{mod_count} clusters extended', flush = True)
    print(f'\t\t{cont_count} clusters merged from contig edges', flush = True)

    omes2dist = treecalcs.update_dists(
        phylo, {gcf_hgxs[i]: omes for i, omes in gcf_omes.items()},
        cpus = cpus, omes2dist = omes2dist, func = dist_func, 
        uniq_sp = uniq_sp, i2ome = i2ome
        )

    return gcfs, gcf_hgxs, gcf_omes, gcf2lg, omes2dist
