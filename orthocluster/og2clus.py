#! /usr/bin/env python3

#NEED to make gbc_mngr_3 split different families with the same OCG
#NEED to model rearrangement for lineages instead of unexpected microsynteny
#NEED to change border percentile entries to clus percentile
#NEED to allow hard threshold for microsynteny distance
#NEED locus-based OCG hypergeometric average
#NEED to post filter clusters after grouping
#NEED to accept an alternative OG input
#NEED to update logging and intelligent resuming
    # md5sum for database integrity
#NEED TO OUTPUT FAMILY TO info.out
#NEED TO ADD PERCENTILE TO OUTPUT
#NEED to implement gff2svg
#NEED to delete excess when checkpoints are reached
    # implement tarring effectively

import os
import re
import sys
import copy
import gzip
import time
import shutil
import random
import pickle
import argparse
import subprocess
import numpy as np
import multiprocessing as mp
from datetime import datetime
from itertools import combinations
from collections import defaultdict, Counter
from mycotools.lib.kontools import \
    intro, outro, format_path, collect_files, \
    dictSplit, file2list, fmt_float, findExecs, \
    multisub, tardir, untardir, checkdir, eprint, \
    write_json, read_json
from mycotools.lib.biotools import \
    gff2list, list2gff, fa2dict, dict2fa
from mycotools.lib.dbtools import mtdb, masterDB
from mycotools.gff2seq import ntmain
from mycotools.acc2gff import grabGffAcc
from mycotools.gff2svg import main as gff2svg
from mycotools.acc2fa import dbmain as acc2fa
from mycotools.db2files import soft_main as symlink_files


def calc_branch_len(phylo, omes):
    """calculate descending branch length from cogent3 tree"""
    # need to verify wtf descending branch length v total of supplied nodes is
    omes = [str(x) for x in omes]
    try:
        return addPatch(phylo.lowest_common_ancestor(
            omes
            ), set(omes)), tuple([int(i) for i in omes])
    except ValueError:
        print(phylo)
        eprint('\t\t' + ','.join([str(x) for x in omes]) + ' raised a tip not found error', flush = True)
        return 0, tuple([int(i) for i in omes])
    except AttributeError:
        print(omes, '\n', phylo)

def update_dists(phylo, cooccur_dict, cpus = 1, omes2dist = {}):
    """update the omes2dist with a new set of data"""
    results = calc_dists(phylo, cooccur_dict, cpus, omes2dist = omes2dist)
    omes2dist = {**omes2dist, **{x[1]: x[0] for x in results}} 
    return omes2dist

def parse_1to1(hg_file, useableOmes = set()):
    derivations = defaultdict(list)
    with open(hg_file, 'r') as raw:
        for line in raw:
            k, v = line.rstrip().split() # all white space
            derivations[k].append(v)

    hgs_list = []
    for k, v in derivations.items():
        v.append(k)
        hgs_list.append(sorted(set(v)))

    hg2gene = {i: v for i, v in enumerate(
                    sorted(hgs_list, key = lambda x: len(x), 
                    reverse = True)
                    )}

    todel = []
    for og, genes in hg2gene.items():
        genes = [x for x in genes if x[:x.find('_')] in useableOmes]
        omes = set([x[:x.find('_')] for x in genes])
        if len(omes) < 2: # skip singleton omes
            todel.append(og)
    for og in todel:
        del hg2gene[og]
    hg2gene = {k: v for k, v in sorted(hg2gene.items(), key = lambda x:
                                       len(x[1]), reverse = True)}

    gene2hg = {}
    for i, og_list in hg2gene.items():
        for gene in og_list:
            gene2hg[gene] = i
       
    i2ome = sorted(set([x[:x.find('_')] for x in list(gene2hg.keys())]))
    ome_num = {v: i for i, v in enumerate(i2ome)}

    return ome_num, gene2hg, i2ome, hg2gene


def parse_orthofinder(hg_file, useableOmes = set()):
    """
    imports orthofinder Orthogroups.txt "hg_file". outputs several data structures:
    ome_num = {ome: number}, gene2hg = {gene: og}, i2ome = [ome0, ome1, ome2]
    """

    gene2hg, ome_num, i2ome, hg2gene = \
        {}, {}, [], {}
    with open(hg_file, 'r') as raw:
        for line in raw:
            data = line.rstrip().split(' ') # OrthoFinder input
            og = int(data[0].replace(':','').replace('OG',''))
            hits = [x for x in data[1:] if x[:x.find('_')] in useableOmes]
            omes = [x[:x.find('_')] for x in hits]
            if len(set(omes)) < 2: # skip singleton organisms
                continue
            hg2gene[og] = hits
            for i, gene in enumerate(hits):
                ome = omes[i]
#                ome = gene[:gene.find('_')] # wtf, parsing omes raises and index error
                if ome not in ome_num:
                    i2ome.append(ome)
                    ome_num[ome] = len(i2ome) - 1
                gene2hg[gene] = og

    return ome_num, gene2hg, i2ome, hg2gene

def run_mmseqs(db, wrk_dir, algorithm = 'mmseqs easy-cluster',
                 min_id = 0.3, min_cov = 0.5, cpus = 1):
    symlink_files(['faa'], db, wrk_dir, verbose = False) # symlink proteomes
    cluster_res_file = wrk_dir + 'homolog_groups.tsv'
    if not os.path.isfile(cluster_res_file): # NEED to add to log removal
        # be extra cautious about shell injection because we need to glob
        int(cpus)
        float(min_id)
        float(min_cov)
        if not os.path.isdir(wrk_dir):
            raise OSError('invalid working directory')
        elif not algorithm in {'mmseqs easy-linclust', 'mmseqs easy-cluster'}:
            raise OSError('invalid mmseqs binary')
        mmseqs_cmd = subprocess.call(algorithm + ' ' + ' '.join([wrk_dir + 'faa/*faa',
                                      wrk_dir + 'cluster', wrk_dir + 'tmp/',
                                      '--min-seq-id', str(min_id), '--threads',
                                      str(cpus), '--compressed', '1',
                                      '--cov-mode', '0', '-c', str(min_cov)]),
                                      shell = True)
 #                                     stdout = subprocess.DEVNULL,
#                                      stderr = subprocess.DEVNULL)
        shutil.move(wrk_dir + 'cluster_cluster.tsv', cluster_res_file)
    elif os.path.getsize(cluster_res_file):
        mmseqs_cmd = 0
    else:
        mmseqs_cmd = 1
    if mmseqs_cmd:
        eprint('\tERROR: cluster failed')
        sys.exit(1)
    if os.path.isfile(wrk_dir + 'cluster_all_seqs.fasta'):
        os.remove(wrk_dir + 'cluster_all_seqs.fasta')
    if os.path.isfile(wrk_dir + 'cluster_rep_seq.fasta'):
        os.remove(wrk_dir + 'cluster_rep_seq.fasta')
    if os.path.isdir(wrk_dir + 'tmp/'):
        shutil.rmtree(wrk_dir + 'tmp/')
    return cluster_res_file

def compile_homolog_groups(hg_file, wrk_dir, 
                           method = 'mmseqs easy-cluster', 
                           useableOmes = set()):
    if method in {'mmseqs easy-cluster', 'mmseqs easy-linclust'}:
        ogInfo = parse_1to1(hg_file, useableOmes)
        hg2genes = ogInfo[-1]
        with open(wrk_dir + 'homolog_groups.txt', 'w') as out:
            for hg, genes in hg2genes.items():
                out.write(str(hg) + '\t' + ' '.join(genes) + '\n')
    elif method == 'orthofinder':
        ogInfo = parse_orthofinder(hg_file, useableOmes)

    with open(wrk_dir + 'ome2i.tsv', 'w') as out:
        out.write(
            '\n'.join([x + '\t' + str(ogInfo[0][x]) for x in ogInfo[0].keys()])
            )
    max_ome =  max([int(x) for x in ogInfo[1].values()])
    return ogInfo


def compile_loci(
    db, ome2i, gene2hg, plusminus, cpus = 1
    ):

    loci_hash_cmds = [
        [x, ome2i[db['ome'][i]], 
        gene2hg, plusminus]
        for i, x in enumerate(db['gff3']) \
        if db['ome'][i] in ome2i
        ]
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        loci_hashes = pool.starmap(parseLoci, loci_hash_cmds)
    pool.join()
    pairs = {x[0]: x[1] for x in loci_hashes}

    return pairs 


def genDistArr(phylo):
    t_dist_arr = phylo.tip_tip_distances()
    dist_arr = np.zeros([len(t_dist_arr.ids), len(t_dist_arr.ids)], dtype=float)
    for i0, v0 in enumerate(t_dist_arr.ids):
        ome0 = int(v0)
        for i1, v1 in enumerate(t_dist_arr.ids):
            ome1 = int(v1)
            dist_arr[ome0, ome1] = float(t_dist_arr[i0, i1])
            dist_arr[ome1, ome0] = float(t_dist_arr[i1, i0])
    return dist_arr

def compileTree(i2ome, tree_path, root = []):
#    if microsynt_dict: # implies the tree needs to be made from scratch
 #       pre_phylo = nj.nj(microsynt_dict, show_progress = False)
  #      # change to ome code for exporting phylogeny
   #     pre_phylo.reassign_names({str(i): v for i, v in enumerate(i2ome)})

    #    if not root or root == 'midpoint':
     #       root_pre = pre_phylo.root_at_midpoint()
      #  elif len(root) == 1:
       #     root_pre = pre_phylo.rooted_with_tip(root[0])
        #else:
         #   root_pre = pre_phylo.rooted_at(
          #      pre_phylo.get_edge_names(root[0], root[1])[0]
           #     )
       # root_pre.write(tree_path, with_distances = True)

    phylo = load_tree(tree_path) 
    if root:
        if len(root) > 1:
            phylo = phylo.rooted_at(
                phylo.get_edge_names(root[0], root[1])[0]
                    )
        else:
            phylo = phylo.rooted_with_tip(root[0])
        phylo.write(tree_path, with_distances = True)

    phylo.reassign_names({v: str(i) for i, v in enumerate(i2ome)})
    return phylo


def compileCDS(gff_list, ome):
    """
    Inputs the gff_list and organism ome code. Compiles CDS entries for loci parsing and outputs:
    cds_dict = {contig: {protein: [[CDSstart, CDSstop]] } }
    Sorts the cds_dict for each protein from smallest to largest
    """

    cds_dict, fail = defaultdict(dict), False
    for entry in gff_list:
        if entry['type'] == 'CDS':
            prot_prep_i0 = entry['attributes'].index(';Alias=') # grab the
            # mycotools accession index
            try: # grab the end of the Alias by grabbing the index for a ';',
            # that is after the accession tag's occurrence.
            # then add that to the start of the accession tag and add one to
            # remove that semicolon. If there is a value error, the accession 
            # is at the end of the attributes (gff col 8 (?))
                prot_prep_i1 = entry['attributes'][prot_prep_i0+1:].index(';') + prot_prep_i0+1
            except ValueError:
                prot_prep_i1 = len(entry['attributes'])
            prot = entry['attributes'][prot_prep_i0:prot_prep_i1].replace(';Alias=','')
            # obtain the protein accession
            if prot not in cds_dict[entry['seqid']] and prot: # add the protein
            # to the cds_dict entry
                cds_dict[entry['seqid']][prot] = [] 
            elif not prot: # if there isn't a valid accession it may mean the
            # mycotools curation did not work or the user did not curate
            # correctly
                print(entry['attributes'], prot_prep_i0, prot_prep_i1)
                if not fail:
                    print('\tWARNING: ' + ome + ' has proteins in gff with no Alias', flush = True)
                fail = True
                continue
            cds_dict[entry['seqid']][prot].extend([int(entry['start']),
                int(entry['end'])]) # add the coordinates

    for contig in cds_dict:
        for prot in cds_dict[contig]:
            cds_dict[contig][prot].sort() # sort the coordinates of the proteins
            # lowest to highest
        cds_dict[contig] = list(sorted(cds_dict[contig].keys(), key = lambda k:
            cds_dict[contig][k][0])) # sort the proteins in the contigs lowest
            # to highest coordinates

    return dict(cds_dict)


def compileCDS2(gff_list, ome):
    """
    Inputs the gff_list and organism ome code. Compiles CDS entries for loci parsing and outputs:
    cds_dict = {contig: {protein: [[CDSstart, CDSstop]] } }
    Sorts the cds_dict for each protein from smallest to largest
    """

    cds_dict, cds_dict2, fail = defaultdict(dict), {}, False
    for entry in gff_list:
        if entry['type'] == 'CDS':
            prot_prep_i0 = entry['attributes'].index(';Alias=') # grab the
            # mycotools accession index
            try: # grab the end of the Alias by grabbing the index for a ';',
            # that is after the accession tag's occurrence.
            # then add that to the start of the accession tag and add one to
            # remove that semicolon. If there is a value error, the accession 
            # is at the end of the attributes (gff col 8 (?))
                prot_prep_i1 = entry['attributes'][prot_prep_i0+1:].index(';') + prot_prep_i0+1
            except ValueError:
                prot_prep_i1 = len(entry['attributes'])
            prot = entry['attributes'][prot_prep_i0:prot_prep_i1].replace(';Alias=','')
            # obtain the protein accession
            if prot not in cds_dict[entry['seqid']] and prot: # add the protein
            # to the cds_dict entry
                cds_dict[entry['seqid']][prot] = [] 
                cds_dict2[prot] = []
            elif not prot: # if there isn't a valid accession it may mean the
            # mycotools curation did not work or the user did not curate
            # correctly
                print(entry['attributes'], prot_prep_i0, prot_prep_i1)
                if not fail:
                    print('\tWARNING: ' + ome + ' has proteins in gff with no Alias', flush = True)
                fail = True
                continue
            cds_dict[entry['seqid']][prot].extend([int(entry['start']),
                int(entry['end'])]) # add the coordinates
            cds_dict2[prot].append(entry)

    for contig in cds_dict:
        for prot in cds_dict[contig]:
            cds_dict[contig][prot].sort() # sort the coordinates of the proteins
            # lowest to highest
        cds_dict[contig] = list(sorted(cds_dict[contig].keys(), key = lambda k:
            cds_dict[contig][k][0])) # sort the proteins in the contigs lowest
            # to highest coordinates

    return dict(cds_dict), cds_dict2


def parseLoci(
    gff_path, ome_num, gene2hg, plusminus = 6
    ):
    """obtain a set of tuples of HG pairs {(OG0, OG1)...}"""

    gff_list = gff2list(gff_path) # open here to improve pickling
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))
    pairs = []
    for scaf in cds_dict: # for each contig
        for i, seq in enumerate(cds_dict[scaf]): # for each gene
            if i < plusminus: # does the locus start at the contig border?
                locus = cds_dict[scaf][:i+plusminus+1]
            else: # or I extend it to the +/- or opposite border
                locus = cds_dict[scaf][i-plusminus:i+plusminus+1]
            loc_og = []
            for gene in locus:
                try: # failure should be relatively rare, so faster than `if`
                    loc_og.append(gene2hg[gene])
                except KeyError: # it's either a singleton or not in an OG
                    pass
            pairs.extend([tuple(sorted(x)) for x in combinations(set(loc_og), 2)]) # don't deal with
            # same OGs - tandem duplications will just naturally cooccur more often

    out_pairs = set(pairs) # unique pairs of OGs

    return ome_num, out_pairs

def form_cooccur_dict(cooccur_dict):

    # sort the values by the length of the combination
    cooccur_dict = {k: tuple(sorted(v)) for k, v in sorted(cooccur_dict.items(), key = lambda x: len(x[0]))}

    hgx2i, i2hgx = {}, {}
    for i, hgx in enumerate(list(cooccur_dict.keys())):
        i2hgx[i] = hgx
        hgx2i[hgx] = i

    return i2hgx, hgx2i, cooccur_dict


def form_cooccur_array(cooccur_dict, ome_len):

    count, hgx2i, size_dict, cooccur_arrays, i2hgx = 0, {}, {}, {}, {}
    cooccur_dict = {k: tuple(sorted(v)) for k, v in sorted(cooccur_dict.items(), key = lambda x: len(x[0]))}
    cooccur_arrays = np.zeros([ome_len, len(cooccur_dict)], dtype = np.int32)

    old_len = len(list(cooccur_dict.keys())[0])
    for i, hgx in enumerate(list(cooccur_dict.keys())):
        i2hgx[i] = hgx
        hgx2i[hgx] = i
        for ome in cooccur_dict[hgx]:
#            cooccur_arrays[ome, i] = hgx_dict[ome][hgx]
            cooccur_arrays[ome, i] = 1

    return i2hgx, hgx2i, cooccur_dict, cooccur_arrays


def prep_microsynt_dict(cooccur_arr):

 #   microsynt_arr = np.zeros((cooccur_arr.shape[0], cooccur_arr.shape[0]))
    microsynt_dict = {}
    i0 = 0
    for row0 in cooccur_arr[:, ]:
        microsynt_dict[(i0, i0)] = 0
        i1 = i0 + 1
        for row1 in cooccur_arr[i1:, ]:
            row_sum = row0 + row1
            try:
                dist = 1 - (len(row_sum[row_sum > 1]) / len(row0[row0 > 0]))
                microsynt_dict[(i0, i1)] = dist
            except ZeroDivisionError: # no overlap
                microsynt_dict[(i0, i1)] = 0
#            microsynt_arr[(i1, i0)] = dist
            i1 += 1
        i0 += 1

    return microsynt_dict

try:
    from numba import njit, jit
    @njit
    def est_conds(ome_arr, cooccur_array):
        p = 1
        for i in range(len(ome_arr) - 1): # for all but the last organisms hgx
            others, ome0 = ome_arr[i + 1:], cooccur_array[ome_arr[i]] #all others and i's hits
            other_m = cooccur_array[others, :] #grab others hits
            other_sum = np.sum(other_m, axis = 0) #sum the others' hits
            other_sum[other_sum > 0] = 1 #convert to binary
            sucs_prep = ome0 + other_sum 
            sucs = len(sucs_prep[sucs_prep > 1]) #find overlapping clusters
            tot = np.sum(ome0)
            p *= sucs/tot # conditional p = overlap/total i
        return p


except ModuleNotFoundError:
    print('Package numba: error: the following arguments are required: filename not detected. May increase throughput', flush = True)

    def est_conds(ome_arr, cooccur_array):
        p = 1
        for i in range(len(ome_arr) - 1):
            others, ome0 = ome_arr[i + 1:], cooccur_array[ome_arr[i]]
            other_m = cooccur_array[others, :]
            other_sum = np.sum(other_m, axis = 0)
            other_sum[other_sum > 0] = 1
            sucs_prep = ome0 + other_sum 
            sucs = len(sucs_prep[sucs_prep > 1])
            tot = np.sum(ome0)
            p *= sucs/tot
        return p

def est_combo_probs(hgx, omes, genesInOGinOme, genesInOme, window,
                    cooccur_arr):
#    win_size = [win_size[x] for x in win_size]
#    p_coeff = calc_coeff(og0, og1, tot_genes, win_size)
    p_conds = est_conds(np.array(omes), cooccur_arr)
    p_coef = 1
    for og in hgx:
        p_coef *= est_hypergeo(omes, genesInOGinOme[og],
                               genesInOme, window)
    p = p_conds * p_coef
    return hgx, p


def est_hypergeo(
    omes, genesInOGinOme, genesInOme, window
    ):
    # NEED to modify to account for sliding window

    # should this be for all omes, or each og in each ome and then multiplied?
    pval = hypergeom.sf(1, genesInOme, genesInOGinOme, window)
    return pval


def combo_prob_mngr(
    hgx2omes, omes2hg2genes, omes2genes, window, cooccur_array, cpus = 1
    ):

    cmds = []
    for hgx, omes in hgx2omes.items(): # take each hgx
        ome = omes[-1] # coefficient ome
        genesInOGinOme = {og: len(omes2hg2genes[ome][og]) for og in hgx}
        genesInOme = len(omes2genes[ome])
        cmds.append([
            hgx, omes, genesInOGinOme, genesInOme, window,
            cooccur_array
            ])

    # run og-by-og, then accumulate via hgx at the end
    with mp.get_context('fork').Pool(processes = cpus) as pool: # will fail on Windows
        hypergeoRes = pool.starmap(est_combo_probs, cmds)

    comparisons = len(hgx2omes) # for Bonferroni correction
    hgx2pval = {}
    for hgx, pval in hypergeoRes:
        hgx2pval[hgx] = comparisons * pval

    return hgx2pval


def gen_null_dict(combo_dict, sample = 10000):
    """combo_dict = {ome: [(OG0...OGn)]}"""

    nulls = []
    for combos in list(combo_dict.values()):
        nulls.extend(list(combos)) # extend all combos to the null

    if not isinstance(combo_dict[list(combo_dict.keys())[0]], set):
        combo_dict = {k: set(v) for k, v in combo_dict.items()}


    null_set_list = list(set(nulls)) # acquire unique combinations
    try: # sample
        null_list = random.sample(null_set_list, sample)
    except ValueError:
        print('\t\t\t\tWARNING: requested null sample greater than HGx size', flush = True)
        null_list = null_set_list

    cooccur_dict = {}
    for null in null_list: # for combo
        cooccur_dict[null] = []
        for ome, combos in combo_dict.items():
            if null in combos:
                cooccur_dict[null].append(ome)

    cooccur_dict = {x: cooccur_dict[x] for x in cooccur_dict}
    # {(OG0,OGn): [omei0...]}

    hgx2i, i2hgx, cooccur_dict = form_cooccur_dict(cooccur_dict)
    return hgx2i, i2hgx, cooccur_dict


def calc_dists(phylo, cooccur_dict, cpus = 1, omes2dist = {}):
    # multiprocessing calculating only new distances for omes2dist
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        results = pool.starmap(
            calc_branch_len, 
            [(phylo, x,) for x in list(set(cooccur_dict.values())) \
            if x not in omes2dist]
            )

    return results


def gen_nulls(pairs, phylo, samples = 10000, cpus = 1):
    """pairs = {ome: [(OG0,OG1)...]}"""

    hgpair2i, i2hgpair, null_dict = gen_null_dict(pairs, samples)
    oldLen = len(null_dict)
    null_dict = {k: v for k, v in null_dict.items() if len(v) > 1}
    results = calc_dists(phylo, null_dict, cpus = cpus)
    omes2dist, pair_scores = {x[1]: x[0] for x in results}, []
    for k in null_dict:
        pair_scores.append(omes2dist[null_dict[k]])
    for i in range(oldLen - len(null_dict)):
        pair_scores.append(0)

    return omes2dist, sorted(pair_scores)


def form_cooccur_structures(pairs, min_omes, ome_len, cc_arr_path = None):
    """
    Imports the out_dicts from parseLoci and creates index vectors that bear
    the ome_num's with a given  cooccurence. e.g. cooccur_dict[(og1, og2)] = [1, 2, 3]
    """

    cooccur_dict = defaultdict(list)
    for ome in pairs:
        for key in pairs[ome]:
            new_key = tuple(sorted(key))
            cooccur_dict[new_key].append(ome)

    cooccur_dict = {
        x: tuple(sorted(cooccur_dict[x])) for x in cooccur_dict \
        if len(cooccur_dict[x]) > min_omes
        }

    i2hgpair, hgpair2i, cooccur_dict, cooccur_array = form_cooccur_array(
        cooccur_dict, ome_len
        )

    return cooccur_array, dict(cooccur_dict), hgpair2i, i2hgpair


def id_near_schgs(hg2gene, omes, max_hgs = 10000, max_median = 2, max_mean = 2):
    """Identify near single copy homology groups. Criteria:
       1) has all omes; 2) minimal overall size of homology group;
       3) lowest median"""
    near_schgs = []
    min_hg2gene = {k: v for k, v in sorted(hg2gene.items(), 
                    key = lambda x: len(x[1])) \
                   if len(v) >= len(omes)}

    max_len = max_mean * len(omes) # dont need to compute means iteratively
    median_i = round((len(omes) / 2) - 0.5)

    for hg, genes in min_hg2gene.items():
        hg_omes = [x[:x.find('_')] for x in genes]
        if not omes.difference(hg_omes): # all omes are present
            count_omes = list(Counter(hg_omes).values()) # prepare for median calculation
            count_omes.sort()
            median = count_omes[median_i]
            print(median, max_median, len(genes), max_len)
            if median <= max_median and len(genes) <= max_len:
                print('\t', median, len(genes))
                near_schgs.append(hg)
            elif len(genes) >= max_len: # everything else will be higher
                break
            if len(near_schgs) == max_hgs:
                break

    return near_schgs


def extract_nschg_pairs(nschgs, hgpair2i, m_arr):
    nschg_set = set(nschgs)
    nschg_pairs = [i for hgpair, i in hgpair2i.items() \
                   if any(x in nschg_set for x in hgpair)]
    nschgs_arr = np.array(nschg_pairs)
    valid_arr = m_arr[:, nschgs_arr]
    return valid_arr
    

def trim_microsynt_arr(m_arr):#, perc = 0):
    ome_num, cols = m_arr.shape
    sum_arr = m_arr.sum(axis = 0) #/ ome_num # for percentage
#    tru_arr = np.where(sum_arr <= perc) # trim via percentage
    tru_arr = np.where(sum_arr <= 1)
    trm_arr = np.delete(m_arr, tru_arr, axis = 1)
    return trm_arr

def align_microsynt_np(m_arr, i2ome, hg2gene, hgpair2i, wrk_dir):
#    trm_arr = trim_microsynt_arr(m_arr)
    nschgs = id_near_schgs(hg2gene, set(i2ome), max_hgs = 100, 
                           max_median = 4, max_mean = 3)
    trm_arr = extract_nschg_pairs(nschgs, hgpair2i, m_arr)
    with open(wrk_dir + 'microsynt.align.phy', 'w') as out:
        out.write(f'{trm_arr.shape[0]} {trm_arr.shape[1]}\n')
        [out.write(f'{i2ome[i]} {"".join([str(x) for x in trm_arr[i]])}\n') \
                   for i in range(trm_arr.shape[0])]
    return wrk_dir + 'microsynt.align.phy'


def run_tree(alignment, wrk_dir, constraint = False, iqtree = 'iqtree',
             model = 'GTR2+FO+ASC+R5', verbose = False, cpus = 1):

    tree_dir = wrk_dir + 'tree/'
    if not os.path.isdir(tree_dir):
        os.mkdir(tree_dir)

    prefix = tree_dir + 'microsynt'
    tree_cmd = [iqtree, '-s', alignment, '-m', model,
                '-p', prefix]
    if constraint:
        tree_cmd.extend(['-te', constraint])
        comp_file = prefix + '.treefile'
    else: # use bootstrap
        tree_cmd.extend(['-B', '1000'])
        comp_file = prefix + '.contree'
    
    if verbose:
        v = None
    else:
        v = subprocess.DEVNULL
    tree_res = subprocess.call(tree_cmd, stdout = v, stderr = v)
    shutil.copy(comp_file, wrk_dir + '../microsynt.newick')

    return tree_res


            

def remove_nulls(cc_arr):
    sum_arr = np.sum(cc_arr, axis = 1) # sum all ome rows
    null_i_list = list(np.where(sum_arr == 0)[0])
    del_list = sorted(null_i_list, reverse = True)
    for i in sorted(null_i_list, reverse = True):
        cc_arr = np.delete(cc_arr, i, axis = 0)
    return cc_arr, del_list


def load_seedScores(file_, seed_thresh):#, seed_thresh):

    out_hgs = []
    with gzip.open(file_, 'rt') as raw:
        for line in raw:
            if not line.startswith('#'):
                data = [x.rstrip() for x in line.split('\t')]
                if float(data[3]) > seed_thresh:
                    out_hgs.append(line.split('\t'))

    return out_hgs


def form_hgpairDict(out_hgs):

    hgs = set([int(x[0]) for x in out_hgs]) # make a set of the first og in an
    # og-pair
    hgs = hgs.union(set([int(x[1]) for x in out_hgs])) # add in the second og
    hgpair_dict = {x: set() for x in hgs} # prepare the keys to bypass `if`'s
    for i in out_hgs:
        hgpair_dict[int(i[0])].add(int(i[1])) # {og0: set(og1, ..., ogn), } or 
        # {og0: set(og0, og1, ..., ogn)}
        hgpair_dict[int(i[1])].add(int(i[0]))

    return hgpair_dict


def hash_protoclusters(gff_path, hgpair_dict, gene2hg, clusplusminus = 10):
    """parse a gff and compile its organized CDS_dict. identify where og-pairs
    co-occur and retrieve the set of hgs for the locus with the locus seed
    protein"""

    gff_list, protoclus = gff2list(gff_path), defaultdict(list)
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))

    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]): # for each index and protein
            try:
                og0 = gene2hg[seq0] 
            except KeyError: # gene not in OG / is a singleton
                continue
            if og0 in hgpair_dict: # if this og is part of a significant seed
                if i0 < clusplusminus: # check if index - clusplusminus < 0
                    locus = cds_dict[scaf][:i0+clusplusminus+1] # obtain locus
                    # as all the sequences from the beginning of the contig to
                    # i0 plus the locus size tolerance (add 1 for python)
                else:
                    locus = cds_dict[scaf][i0-clusplusminus:i0+clusplusminus+1]
                    # otherwise obtain the locus from the +/-
                for seq1 in locus:
                    if seq1 != seq0: # if the protein isn't the first queried
                        try:
                            og1 = gene2hg[seq1]
                        except KeyError:
                            continue
                        if og1 in hgpair_dict[og0]: # if it is in the set of
                        # hgs for the first queried og
                            og_loc_list = []
                            for gene in locus:
                                try:
                                    og_loc_list.append(gene2hg[gene])
                                except KeyError: # singleton or missing
                                    pass
                            og_loc = set(og_loc_list)
                            # the set of hgs in this locus
                            hgpair = tuple(sorted([og0, og1])) # sort the hgs
                            protoclus[hgpair].append([og_loc, seq0, seq1])
                            # {(og0, og1}: [[set(hgx, .. ogy), seq0, seq1], ...]} 
                            
    return dict(protoclus)


def merge_protos(protoclus_res):
    """combine protoclus results from each organism"""

    protohgx2omes = defaultdict(list)
    for ome_protoclus in protoclus_res: # for each protoclus from the mp
    # results
        for hgpair in ome_protoclus: 
            protohgx2omes[hgpair].extend(ome_protoclus[hgpair])
            # {(og0, og1)}: [[set(hgx, ogy), seq], ...]

    return dict(protohgx2omes)


def gen_clusters(loc0, loc1, hgpair_set):

    ome0 = loc0[1][:loc0[1].find('_')] # extract ome from gene acc
    ome1 = loc1[1][:loc1[1].find('_')]
    if ome0 != ome1:
        loc0vloc1 = loc0[0].intersection(loc1[0])
        if loc0vloc1 != hgpair_set: #if there are more overlapping hgs
            hgx = tuple(sorted(list(loc0vloc1))) # sort the loci (not by
            return hgx, ome0, ome1, loc0[1], loc1[1]


def FormulateClusDicts(gen_clus_cmds, hgx2omes, hgx2loc, ome2i, cpus):

    with mp.get_context('fork').Pool(processes = cpus) as pool:
        clus_res = pool.starmap(gen_clusters, gen_clus_cmds)
    for i, res in enumerate(clus_res):
        if res:
            for comb_len in range(3, len(res[0]) + 1):
                hgx = res[0]
                if hgx not in hgx2omes:
                    hgx2omes[hgx] = set()
                    hgx2loc[hgx] = set()
                hgx2omes[hgx] = hgx2omes[hgx].union({
                     ome2i[res[1]], ome2i[res[2]]
                     })
                hgx2loc[hgx] = hgx2loc[hgx].union({res[3], res[4]})

    return hgx2omes, hgx2loc


def par_rm(hgx, hgx2omes):

    t_comb_sets, t_combs, todel = [set(hgx)], [hgx], []
    for comb_len in reversed(range(3, len(hgx))):
        for comb in combinations(hgx, comb_len):
            set_comb = set(comb)
            for i, set_hgx in enumerate(t_comb_sets):
                if set_comb < set_hgx:
                    try:
                        if hgx2omes[comb] == hgx2omes[t_combs[i]]:
                            todel.append(comb)
                            break
                    except KeyError:
                        break
            t_combs.append(set_comb)

    return todel


def rm_subsets(hgx2omes, hgx2loc, cpus = 1):

    hgx2omes = {
        k: hgx2omes[k] for k in sorted(
            list(hgx2omes.keys()), 
            key = lambda x: len(x), reverse = True
            )
        }
    max_len = len(list(hgx2omes.keys())[0])
    for hgx_len in reversed(range(4, max_len + 1)):
        i2hgx = list(hgx2omes.keys())
        rm_cmds = [[hgx, hgx2omes] for hgx in i2hgx if len(hgx) == hgx_len]
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            todels = pool.starmap(par_rm, rm_cmds)
        for todel in todels:
            for hgx in todel:
                try:
                    del hgx2omes[hgx]
                    del hgx2loc[hgx]
                except KeyError:
                    continue

    return hgx2omes, hgx2loc
                            

def hgpair2hgx(db, hgpair_dict, gene2hg, ome2i, cpus, clusplusminus = 10):

    hash_protoclus_cmds = []
    gffs = [v['gff3'] for k, v in db.items() if k in ome2i]
    print('\t\tForming protocluster hashes', flush = True)
    for gff in gffs: # prepare for protocluster hashing by organism
        hash_protoclus_cmds.append((gff, hgpair_dict, gene2hg, clusplusminus,))
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        protoclus_res = pool.starmap(hash_protoclusters, hash_protoclus_cmds)
    pool.join()

    print('\t\tForming proto-HGxs', flush = True)
    protohgx2omes = merge_protos(protoclus_res) # merge mp results

    print('\t\tGenerating high order HGxs', flush = True)
    print('\t\t\t' + str(sum([
           len(protohgx2omes[x])**2 - len(protohgx2omes[x]) \
           for x in protohgx2omes
        ])) + ' to run', flush = True)
    count = 0

    hgx2omes, hgx2loc, gen_clus_cmds = {}, {}, []
    for hgpair in protohgx2omes:
        hgpair_set = set(hgpair)
        for loci in combinations(protohgx2omes[hgpair], 2):
            loc0, loc1 = loci[0], loci[1]
            gen_clus_cmds.append([loc0, loc1, hgpair_set])
        if len(gen_clus_cmds) > 5000000:
            count += len(gen_clus_cmds)
            print('\t\t\t' + str(count), flush = True)
            hgx2omes, hgx2loc = FormulateClusDicts(
                gen_clus_cmds, hgx2omes, hgx2loc, ome2i, cpus
                )
            gen_clus_cmds = []

    hgx2omes, hgx2loc = FormulateClusDicts(
        gen_clus_cmds, hgx2omes, hgx2loc, ome2i, cpus
        )

    print('\t\tRemoving subset HGxs', flush = True)
    hgx2omes, hgx2loc = rm_subsets(hgx2omes, hgx2loc, cpus)
    hgx2omes = {x: tuple(sorted(list(hgx2omes[x]))) for x in hgx2omes}

    return hgx2omes, hgx2loc


def Hash4nulls(
    gff_path, gene2hg, max_size, plusminus
    ):

    gff_list = gff2list(gff_path) # open here for decreased serialization
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))
    sml_sizes = {i: [] for i in range(3, (plusminus * 2) + 2)} # for HGx sizes
    if max_size > (plusminus*2) + 1: # if biggest HGx is beyond the window size
        lrg_sizes = {i: [] for i in range((plusminus * 2) + 1, max_size + 1)}
        # create a separate dictionary to handle these
    else:
        lrg_sizes = {}

    for scaf, seqs in cds_dict.items():
        for i, seq in enumerate(seqs):
            if i < plusminus: # compile the locus
                locus = seqs[:i+plusminus+1]
            else:
                locus = seqs[i-plusminus:i+plusminus+1]
            loc_og = set([gene2hg[x] for x in locus if x in gene2hg]) # translate to OGs
            for size, size_list in sml_sizes.items():
                size_list.extend((
                    [tuple(sorted(x)) for x in combinations(loc_og, size)]
                    ))
        for size, size_list in lrg_sizes.items():
            size_list.extend([
                tuple(sorted(seqs[x:x+size])) for x in range(len(scaf)) \
                if len(seqs[x:x+size+1]) == size
                ])

    nulls = {
        **{k: set(v) for k,v in sml_sizes.items()}, 
        **{k: set(v) for k,v in lrg_sizes.items()}
        } # make them unique sets of combinations
    return nulls


def genHGxNulls(
    gffs, gene2hg, max_clus_size, plusminus, phylo,
    hgx_perc, clus_perc, wrk_dir,
    omes2dist = {}, samples = 10000, cpus = 1
    ): # NEED to adjust; this currently uses way too much memory
    """Generates a null distribution of randomly sampled HGxs for each 
    # of orthogroups observed in HGxs. Applies the percentile for each
    size as the minimum value for significance.
    Outputs a dictionary of the minimum values for each size HGx
    based on the inputted percentiles. {# of OGs: minimum value}"""

    print('\t\tParsing for random samples', flush = True)
    hash_null_cmds = [
        (x, gene2hg, max_clus_size, plusminus,) \
        for x in gffs
        ]
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        hashRes = pool.starmap(Hash4nulls, hash_null_cmds)
        # hashRes = [({size: [HGx]})...] by ome

    print('\t\tCalculating microsyteny distances of random samples\n\t\tHGx size:', flush = True)
    hgxBordPercs, hgxClusPercs = {}, {}
    for size in range(3, max_clus_size + 1): # for all observed sizes
        print('\t\t\t' + str(size), flush = True)
        size_dict = {i: v[size] for i, v in enumerate(hashRes)}
        # spoof i keys, size_dict values
        hgx2i, i2hgx, null_dict = gen_null_dict(size_dict, samples)
        oldLen = len(null_dict)
        null_dict = {k: v for k, v in null_dict.items() if len(v) > 1} # 0 will n t compute
        distRes = calc_dists(phylo, null_dict, omes2dist = omes2dist, cpus = cpus)
        omes2dist, scores = {
            **omes2dist, **{x[1]: x[0] for x in distRes}
            }, []

        nullSizes = []
        with open(wrk_dir + str(size) + '.null.txt', 'w') as out:
            for omes in list(null_dict.values()):
                nullSizes.append(omes2dist[omes])
                out.write(str(omes2dist[omes]) + '\n')
            out.write('\n'.join([str(0) for x in range(oldLen-len(null_dict))])) 
            # account for missing values

        # sort to apply percentile threshold for each size HGx and for both the
        # border percentile and cluster percentile
        nullSizes.sort()
        hgxBordPercs[size] = nullSizes[round(hgx_perc * len(nullSizes) + .5)]
        hgxClusPercs[size] = nullSizes[round(clus_perc * len(nullSizes) + .5)]

    return hgxBordPercs, hgxClusPercs


def loadNulls(max_clus_size, wrk_dir, hgx_perc, clus_perc):

    hgxBordPercs, hgxClusPercs = {}, {}
    for size in range(3, max_clus_size + 1):
        nullSizes = []
        with open(wrk_dir + str(size) + '.null.txt', 'r') as raw:
            for line in raw:
                nullSizes.append(float(line.rstrip()))
        nullSizes.sort()
        hgxBordPercs[size] = nullSizes[round(hgx_perc * len(nullSizes) + .5)]
        hgxClusPercs[size] = nullSizes[round(clus_perc * len(nullSizes) + .5)]

    return hgxBordPercs, hgxClusPercs


def hash_hgx(gff_path, ome, hgx_genes, gene2hg, clusplusminus):
    """sliding window of size clusplusminus, identify sequences that may have
    hgs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus"""

    gff_list, gene_dict = gff2list(gff_path), {}
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))

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

def findHGxPairSingle(gene2hgx, ome, hgx2i, pairsDict, minimum = 2):
    omePairs_dict = defaultdict(list)
    for gene, hgxs in gene2hgx.items():
        for pairTup in combinations(hgxs, 2):
            omePairs_dict[pairTup].append(gene) # gene centric 

    for pairTup, genes in omePairs_dict.items():
        if len(genes) >= minimum:
            hgxpair = tuple(sorted([hgx2i[pairTup[0]], hgx2i[pairTup[1]]]))
            pairsDict[hgxpair].append(ome)
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

def WriteAdjMatrix(Q, out_file):
    with open(out_file, 'w') as out: # might be nice to compress this/write to MCL binary matrix
        x = True
        while x:
            x = Q.get()
            if x:
                out.write(x)
#                out.flush() # shouldn't be a bottleneck
    
def BLASTclanOG(db, og, fileBase, genes, minid, blastp = 'blastp'):
    blast_ids = defaultdict(dict)
    if not os.path.isfile(fileBase + '.out'):
        fa_dict = acc2fa(db, genes)
        with open(fileBase + '.fa', 'w') as out:
            out.write(dict2fa(fa_dict))
#        makeDBcmd = subprocess.call([
 #           diamond, 'makedb', '--in', fileBase + '.fa', '--db',
  #          fileBase + '.dmnd', '--threads', str(2)
   #         ], stdout = subprocess.DEVNULL,
    #        stderr = subprocess.DEVNULL
     #       )
#        dmndBlast = subprocess.call([diamond, 'blastp', '--query', fileBase + '.fa', 
 #           '--db', fileBase, '--threads', str(2), '--id', str(minid), '--no-self-hits',
  #          '-o', fileBase + '.out.tmp', '--outfmt', '6', 'qseqid', 'sseqid', 'pident'
   #         ], stdin = subprocess.PIPE, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
        blastp = subprocess.call(['blastp', '-query', fileBase + '.fa', 
            '-subject', fileBase + '.fa', '-num_threads', str(2),
            '-out', fileBase + '.out.tmp', '-outfmt', '6 qseqid sseqid pident'
            ], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
        shutil.move(fileBase + '.out.tmp', fileBase + '.out')

    with open(fileBase + '.out', 'r') as raw:
        for line in raw:
            q, s, pident = line.split()
            if q not in blast_ids[og]:
                blast_ids[og][q] = {}
            blast_ids[og][q][s] = float(pident)/100 # adjust to decimal

    return dict(blast_ids)

def mpacquire_clus_ocg_sim(
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
    # primary HGx into separate ocgs

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


def acquire_clus_ocg_sim(
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
    # primary HGx into separate ocgs

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



def mpclan_to_ocg_loci(
    db, clanI, loci, ogLoci, ocg_dir, Q, index, 
    diamond = 'diamond', minid = 30, cpus = 1
    ):

    blast_hash = defaultdict(list)
    for i, locus in enumerate(loci):
        hgs = ogLoci[i]
        for i1, og in enumerate(hgs):
            if og is not None:
                blast_hash[og].append(locus[i1])

    blast_ids = defaultdict(dict)
    cmds = []
    for og, genes in blast_hash.items():
        if len(genes) > 1:
            fileBase = ocg_dir + str(clanI) + '.' + str(og)
            cmds.append([
                db, og, fileBase, genes, minid, diamond
                ])

    with mp.get_context('fork').Pool(processes = cpus) as pool:
        tblast_res = pool.starmap(BLASTclanOG, cmds)

    for tblast_ids in tblast_res:
        for og, qdict in tblast_ids.items():
            for q, sdict in qdict.items():
                if q in blast_ids[og]:
                    blast_ids[og][q] = {**blast_ids[og][q], **sdict}
                else:
                    blast_ids[og][q] = sdict

    blast_ids = dict(blast_ids)
    setOGloci = [set(v) for v in ogLoci] # prestage set making
    [v.remove(None) for v in setOGloci if None in v] # prestage removal
    if blast_ids:
        procs = []
        for i0, ogLoc0 in enumerate(ogLoci):
            for i1, ogLoc1 in enumerate(ogLoci[i0+1:]):
                SogLoc0, SogLoc1 = setOGloci[i0], setOGloci[i1]
                if not SogLoc0.isdisjoint(SogLoc1):
                    print('\t',i0,i1)
                    while len(procs) >= cpus:
                        todel = None
                        for i2, proc in enumerate(procs):
                            if not proc.is_alive():
                                todel = i2
                                proc.join()
                                break
                        if todel is not None:
                            del procs[i2]
                    Thgs = copy.copy(ogLoc0)
                    Thgs.extend(ogLoc1)
                    Tids = {
                        og: blast_ids[og] for og in Thgs if og in blast_ids
                        }
                    procs.append(mp.Process(
                        target = mpacquire_clus_ocg_sim, 
                        args = (
                            i0, i1, index, loci[i0], loci[i1], 
                            ogLoci[i0], ogLoci[i1], SogLoc0, SogLoc1,
                            Tids, Q
                            )
                        ))
                    procs[-1].start()
        while procs:
            todel = None
            for i, proc in enumerate(procs):
                if not proc.is_alive():
                    todel = i
                    proc.join()
                    break
            if todel is not None:
                del procs[i2]
                
#        with mp.get_context('fork').Pool(processes = cpus) as pool:
 #           pool.starmap(
  # ome0_gene1_og,)...,)}             acquire_clus_ocg_sim,
   #             ([i0, i1, index, loci[i0], loci[i1], ogLoci[i0], ogLoci[i1], blast_ids, Q] \
    #            for i0, i1 in combinations(range(len(loci)), 2))
     #           )

def clan_to_ocg_loci(
    db, clanI, loci, ogLoci, ocg_dir, Q, index, 
    diamond = 'diamond', minid = 30
    ):
#    print(clanI, flush = True)
    blast_hash = defaultdict(list)
    for i, locus in enumerate(loci):
        hgs = ogLoci[i]
        for i1, og in enumerate(hgs):
            if og is not None:
                blast_hash[og].append(locus[i1])
    
    blast_ids = defaultdict(dict)
    for og, genes in blast_hash.items():
        if len(genes) > 1:
            fileBase = ocg_dir + str(clanI) + '.' + str(og)
            if not os.path.isfile(fileBase + '.out'):
                fa_dict = acc2fa(db, genes)
                with open(fileBase + '.fa', 'w') as out:
                    out.write(dict2fa(fa_dict))
                makeDBcmd = subprocess.call([
                    diamond, 'makedb', '--in', fileBase + '.fa', '--db',
                    fileBase + '.dmnd', '--threads', str(2)
                    ], stdout = subprocess.DEVNULL,
                    stderr = subprocess.DEVNULL
                    )
                dmndBlast = subprocess.call([diamond, 'blastp', '--query', fileBase + '.fa', 
                    '--db', fileBase, '--threads', str(2), '--id', str(minid), '--no-self-hits',
                    '-o', fileBase + '.out.tmp', '--outfmt', '6', 'qseqid', 'sseqid', 'pident'
                    ], stdin = subprocess.PIPE, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
                shutil.move(fileBase + '.out.tmp', fileBase + '.out')
        
            with open(fileBase + '.out', 'r') as raw:
                for line in raw:
                    q, s, pident = line.split()
                    if q not in blast_ids[og]:
                        blast_ids[og][q] = {}
                    blast_ids[og][q][s] = float(pident)/100 # adjust diamond to decimal
        
    blast_ids = dict(blast_ids)    
    if blast_ids:
        for i0, i1 in combinations(range(len(loci)), 2):
            loc0, loc1 = loci[i0], loci[i1]
            ogL0, ogL1 = ogLoci[i0], ogLoci[i1]
            sOGl0 = set([x for x in ogL0 if x is not None])
            sOGl1 = set([x for x in ogL1 if x is not None])
            if not sOGl0.isdisjoint(sOGl1): # if there is an intersection 
#                print('\t', i0, i1, flush = True)
                acquire_clus_ocg_sim(i0, i1, index, loc0, loc1, ogL0, ogL1, sOGl0,
	                             sOGl1, blast_ids, Q)
    

def hash_clan_loci(ome, gff_path, ome_sig_clus, gene2hg, clusplusminus):
    """sliding window of size clusplusminus, identify sequences that may have
    hgs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus"""

    gff_list, protoclus, clus_out = gff2list(gff_path), {}, []
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))

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

def classify_ocgs(
    hgx2loc, db, gene2hg, i2hgx, hgx2i,
    phylo, clusScores, bordScores, ome2i, hgx2omes,
    wrk_dir, omes2dist = {}, clusplusminus = 3, 
    inflation = 1.5, minimum = 2, min_omes = 2, cpus = 1
    ):

    groupI = wrk_dir + 'group.I.pickle'
    groupII = wrk_dir + 'group.II.pickle'
    if not os.path.isfile(groupI):
        print('\tDiscovering HGxs with overlapping loci', flush = True)

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
        print('\tLoading HGxs with overlapping loci', flush = True)
        with open(groupI, 'rb') as in_:
            gene2hgx, hgx2genes = pickle.load(in_)

    groupIII = wrk_dir + 'group.III.pickle'
    if not os.path.isfile(groupIII):
        print('\tIdentifying significantly overlapping HGxs', flush = True)
        tpairDict = defaultdict(list)
        for ome, gene2hgx_ome in gene2hgx.items():
            tpairDict = findHGxPairSingle(
                gene2hgx_ome, ome, hgx2i, tpairDict,
                minimum = minimum
                )

        pairDict = {}
        for id_, omes in tpairDict.items():
            omes_set = set(omes)
            if len(omes_set) > 1: # need at least one ome to calc a branch length
                pairDict[id_] = tuple(sorted(omes_set)) 


        omes2dist = update_dists(phylo, pairDict, cpus = cpus, omes2dist = omes2dist)

        print('\tClassifying HGx clans and Orthologous Cluster Groups (OCGs)', flush = True)
        # populate a lil_matrix here, then use that network to identify modules
        print('\t\tBuilding binary HGx-HGx network', flush = True)
        matrix = lil_matrix((len(i2hgx), len(i2hgx)), dtype=bool)
        for idPair, omes in pairDict.items():
            i0, i1 = idPair[0], idPair[1]
            nullSize = max([len(i2hgx[x]) for x in idPair])
            if omes2dist[omes] >= bordScores[nullSize]:
                matrix[i0, i1] = True
                matrix[i1, i0] = True
    
        print('\t\tIsolating subgraphs', flush = True)
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

        omes2dist = update_dists(
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

    print('\tClassifying HGx clans into OCGs', flush = True)
    clanFile = wrk_dir + 'clan2loci.'
    if not os.path.isfile(clanFile + 'og.json.gz'):
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
    
        clanLoci, clanHGx4ocgs, clanOGloci = defaultdict(list), defaultdict(list), defaultdict(list)
        for ome, outclanLoci, outclanHGx, outclanOGloci in hash_res:
            for clan in outclanLoci:
                clanLoci[clan].extend(outclanLoci[clan])
                clanHGx4ocgs[clan].extend(outclanHGx[clan])
                clanOGloci[clan].extend(outclanOGloci[clan])
        write_json(clanLoci, clanFile + 'json.gz')
        write_json(clanHGx4ocgs, clanFile + 'hgx.json.gz')
        write_json(clanOGloci, clanFile + 'og.json.gz')
    else:
        print('\t\tLoading clan loci', flush = True)
        clanLoci = {
            int(k): tuple(v) for k,v in sorted(
                read_json(clanFile + 'json.gz').items(), key = lambda x: len(x[1]), reverse = True
            )}
        clanHGx4ocgs = {int(k): tuple([tuple(i) for i in v]) for k,v in read_json(clanFile + 'hgx.json.gz').items()}
        clanOGloci = {int(k): tuple(v) for k, v in read_json(clanFile + 'og.json.gz').items()}

    # dict(clanLoci) = {int(clanI): ((ome0_gene0, ome0_gene1,)...,)}
    # dict(clanHGx4ocgs) = {int(clanI): ((og0, og1, ogn,)...,)}
    # dict(clanOGloci) = {int(clanI): ((ome0_gene0_og, ome0_gene1_og,)...,)}

    ocg_dir = wrk_dir + 'ocg/'
    if not os.path.isdir(ocg_dir):
        os.mkdir(ocg_dir)
   
    print('\t\tCalling OCGs', flush = True)
    print('\t\t\t' + str(sum([len(v) for v in list(clanLoci.values())])) \
        + ' loci considered', flush = True)

    m, adj_mtr = mp.Manager(), ocg_dir + 'loci.adj.tmp'
    Q = m.Queue()
    W = mp.Process(target=WriteAdjMatrix, args=(Q, adj_mtr))
    W.start()

    if cpus < 2:
        cpus = 2
    index, cmds = 0, []

    loci, hgxXloci = {}, {}
    for clanI in clanLoci:
        cmds.append([
            db, clanI, clanLoci[clanI], clanOGloci[clanI],
            ocg_dir, Q, index
            ])
        for i, locus in enumerate(clanLoci[clanI]):
            loci[index] = clanLoci[clanI][i]
            hgxXloci[index] = clanHGx4ocgs[clanI][i]
            index += 1

 #   bigClan = cmds[0][1]
  #  del cmds[0] # remove the first one because it is huge enough to mp on its own
    if not os.path.isfile(ocg_dir + 'loci.adj'):
        with mp.get_context('fork').Pool(processes = cpus - 1) as pool:
            pool.starmap(
                clan_to_ocg_loci, cmds
                )
    #    mpclan_to_ocg_loci(
     #       db, bigClan, clanLoci[bigClan], clanOGloci[bigClan], 
      #      clanHGx4ocgs[bigClan], ocg_dir, Q, 0,
       #     diamond = 'diamond', minid = 30, cpus = cpus - 1
        #    ) # now do the biggest OCG
   
        shutil.move(adj_mtr, ocg_dir + 'loci.adj')

    Q.put(None)
    W.join()

    if not os.path.isfile(ocg_dir + 'loci.clus'):
        print('\t\t\tRunning MCL', flush = True)
        if cpus > 5:
            mcl_threads = 10 # cap at 5 processors, ten threads
        else:
            mcl_threads = (cpus * 2) - 1
        MCL(ocg_dir + 'loci.adj', ocg_dir + 'loci.clus', 
            inflation = 1.5, threads = mcl_threads)
        # could add iterative subsample MCL option here
    
    t_ocgs = []
    with open(ocg_dir + 'loci.clus', 'r') as raw:
        for line in raw: # loci indices
            d = line.rstrip().split('\t')
            if len(d) > 1:
                indices = [int(x) for x in d]
            else: # singletons won't pass thresholds, disregard them
                continue
            t_ocgs.append(indices)

    # list(ocg) = [{hgx: (omes,)}]
    ocgs, ocg_hgxs, ocg_omes = [], [], []
    for ocg, locIs in enumerate(t_ocgs):
        ocgs.append(defaultdict(list))
        ocgHGx, ocgOme_list = [], []
        for locI in locIs:
            loc = loci[locI]
            hgxs = tuple([tuple(hgx) for hgx in hgxXloci[locI]])
            # really should be done above to make hgxs formatted right
            omeI = ome2i[loc[0][:loc[0].find('_')]]
            [ocgs[-1][hgx].append(omeI) for hgx in hgxs]
            [ocgHGx.extend(hgx) for hgx in hgxs]
            ocgOme_list.append(omeI)
        if len(set(ocgOme_list)) > min_omes: # need more than 1
            ocgs[-1] = {k: sorted(set(v)) for k,v in ocgs[-1].items()}
            ocg_hgxs.append(tuple(sorted(set(ocgHGx))))
            ocg_omes.append(tuple(sorted(set(ocgOme_list))))
        else:
            del ocgs[-1]

    print('\t\t\t' + str(len(ocgs)) + ' OCGs w/' \
        + str(sum([len(x) for x in ocgs])) + ' loci', flush = True)

    omes2dist = update_dists(
        phylo, {ocg_hgxs[i]: omes for i, omes in enumerate(ocg_omes)},
        cpus = cpus, omes2dist = omes2dist
        )

    return ocgs, ocg_hgxs, ocg_omes, omes2dist


def extrct_sig_clus(
    clus_scores_dict, hgx2loc, top_hgxs, ome2i
    ):
    """Create a hash to seed retrieving clusters based on
    their entry in hgx2loc."""
 
    sig_clus = defaultdict(dict)
    for top in top_hgxs:
        hgx, clan, omeIs = top[0], top[1], top[2]
        for gene in hgx2loc[hgx]:
            ome = gene[:gene.find('_')]
            if ome2i[ome] in omeIs:
                if gene not in sig_clus[ome]:            
                    sig_clus[ome][gene] = [[set(hgx), clan]]
                else:
                    sig_clus[ome][gene].append([set(hgx), clan])

    return dict(sig_clus)


def calc_counts(phylo, omes, set_omes):

    count = 1
    for sub in iter(phylo):
        if any(str(x)[:str(x).find(':')] in set_omes for x in sub.tips(True)):
            count += calc_counts(sub, omes, set_omes)
        else:
            count += 1

    return count

def addPatch(phylo, omes_set):
    """recursively add branches without the trait to the missing total"""

    total = 0
    for sphylo in phylo:
        subOmes = [str(x)[:str(x).find(':')] for x in sphylo.tips(True)]
        if all(
            x not in omes_set \
            for x in subOmes
            ): # if all descendents aren't in the set, add the total length
            total += sphylo.total_descending_branch_length()
        elif not all(
            x in omes_set \
            for x in subOmes
            ): # if some descendents aren't in the set, search this branch
            total += addPatch(sphylo, omes_set)
    return total

def calc_patchiness(phylo, omes):
    """calculate the percent branch length that a trait is missing over the
    MRCA of where the trait is present"""
    omes_set = set(omes)
    try:
        mrca = phylo.lowest_common_ancestor(omes)
    except ValueError:
        print(phylo, omes)
    except AttributeError: # 1 ome
        eprint('\t\t' + ','.join([str(x) for x in omes]) + ' raised a tip not found error', flush = True)
        print(phylo)
    totalDist = mrca.total_descending_branch_length()
    subDist = addPatch(mrca, omes_set)
    return tuple([int(x) for x in omes]), subDist/totalDist


def makeBlastDB(hg_file, out_file, makeblastdb):
#    og_int = int(re.sub(r'\.fa$', '', os.path.basename(hg_file).replace('OG','')))
 #   out_file = blast_dir + str(og_int) + '.db'
    cmd = [
        makeblastdb, '-in', hg_file, '-dbtype', 'prot',
        '-parse_contigs', '-out', out_file
        ]
    makeblastExec = subprocess.call(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    if makeblastExec:
        print(
            '\t\tERROR: code: ' + str(makeblastExec) + ' ' + \
            os.path.basename(outfile), flush = True
            )


def runBlast(queryfile, blastdb, blastp, outfile, max_seqs = 25):
    
    cmd = [
        blastp, '-query', queryfile, '-db', blastdb,
        '-out', outfile, '-max_hsps', '1', '-max_target_seqs',
        str(max_seqs), '-outfmt', '6 qcontig scontig bitscore pident'
        ]
    blastExec = subprocess.call(
        cmd, stdout = subprocess.DEVNULL,
        stderr = subprocess.DEVNULL
        )
    if blastExec:
        print(
            '\t\tERROR: code: ' + str(blastExec) + ' ' + \
            os.path.basename(outfile), flush = True
            )
#    for seq in queryfa:
 #       faStr = dict2fa({seq: queryfa[seq]}).rstrip() + '\n'
  #      blastExec.stdin.write(bytes(faStr, encoding = 'utf-8'))
   # blastExec.stdin.flush()
    #blastExec.stdin.close()
    #blastExec.wait()


def parseRes(res_file, ome_set, hgx, og):

    hits = defaultdict(list)
    with open(res_file, 'r') as raw:
        for line in raw:
            data = line.rstrip().split('\t')
            query, subject, bit = data[0], data[1], float(data[2])
            s_ome = subject[:subject.find('_')]
            hits[query].append((s_ome, bit))

    hits = {q: sorted(hits[q], key = lambda x: x[1], reverse = True) for q in hits}
    results = {q: 0 for q in hits}
    for q in hits:
        q_ome = q[:q.find('_')]
        for hit in hits[q]:
            if hit[0] != q_ome:
                if hit[0] in ome_set:
                    results[q] = 1
                break

    return results, hgx, og


def calcGBCcorr(hgx_res, hgx):

    omes, total = set(), 0
    for og in hgx_res:
        t_dict = {}
        for q in hgx_res[og]:
            ome = q[:q.find('_')]
            if ome not in t_dict:
                t_dict[ome] = 0
            if hgx_res[og][q] == 1: # if there are any hits for an ome w/multiple queries
                t_dict[ome] = 1
        try:
            total += sum(t_dict.values()) / len(t_dict)
        except ZeroDivisionError:
            print('\t\t' + str(og) + ' no hits', flush = True)
    try:
        total /= len(hgx_res)
    except ZeroDivisionError:
        print('\t\t' + str(hgx) + ' no hits', flush = True)

    return total, hgx


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
                '-6', 'qseqid', 'sseqid', 'evalue'
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
    multisub(cmds1, process = cpus)


    print('\tRunning ' + str(len(cmds2)) + ' diamonds', flush = True)
    multisub(cmds2, processes = cpus)


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

def hgx2omes2gbc_calc(
    hgx, omesI, omes, hgxDict, hgx_dir
    ):
    # hgxDict = {og: set(gene1, gene2)}

#    if hgx == (1225, 1407, 1753, 2875, 3830, 5417, 6211, 6997, 8298, 8898,
 #   11992, 12077, 16076, 33272):
  #      print(omesI, omes)
    res = {}
    for og, qs in hgxDict.items():
        res[og] = {}
        with open(hgx_dir + str(og) + '.out', 'r') as raw:
        # open the blast results
            geneInfo = defaultdict(list) 
            for line in raw:
                d = line.rstrip().split('\t')
                q = d[0] # qryid
                s = d[1] # sbjid
                e = float(d[-1])
                if q in qs:
                    geneInfo[q].append((s, e))
                    # dict(geneInfo) = {query: [(sbj, evalue)]}
        geneInfo = {
            k: sorted(v, key = lambda x: x[1]) for k,v in
            geneInfo.items() 
            } # sort each query by evalue (reverse = True for bit)

        for gene in list(hgxDict[og]): # for each gene
            ome = gene[:gene.find('_')] # identify the ome
            if ome not in res[og]:
                res[og][ome] = False
            if res[og][ome]:
                continue
            try:
                sbjct_ome = geneInfo[gene][0][0][:geneInfo[gene][0][0].find('_')]
                while sbjct_ome == ome:
                # while the first subject of the query is the same ome as the
                # query 
                    geneInfo[gene].pop(0) # delete the first subject
            except IndexError: # if there are no more subjects, just continue
                continue
            except KeyError: # the gene is not in the blast results
                continue # would be nice to quantitate the percent of this

#            if sbjct_ome in omes: # SHARED OMES PASS, a conservative approach
            # to identify if the subject is in the family of omes
#            if geneInfo[gene][0][0] in geneInfo: # ONLY SHARED GENES PASS
            res[og][ome] = geneInfo[gene][0][0]
                # dict(res) = {og: {ome: highest bit score subject}}

    # populate a binary response dictionary for each ome and its genes that are
    # in the shared HGx; 0 = the OG's best blast hit in this ome is not another 
    # ome code that shares the HGx, 1 = the OG's best blast hit is another ome 
    # code that share the HGx
    omeScores = defaultdict(dict)
    for og in res:
        for ome, seq in res[og].items():
            if seq:
                omeScores[ome][og] = 0
                if seq[:seq.find('_')] in omes:
                    omeScores[ome][og] = 1
                # dict(omeScores) = {ome: {og: [0,1]}}
    try:
        gbcScore = \
            sum([sum(omeScores[ome].values()) \
            / len(omeScores[ome]) for ome in omeScores]) \
            / len(omeScores) 
            # total average of ome-by-ome averages of binary 
            # positives
    except ZeroDivisionError:
        gbcScore = 0

    return hgx, omesI, gbcScore



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
    
        print('\tAssimilating OCG loci', flush = True)
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

#    sys.exit()
    hgx2omes2gbc = defaultdict(dict)
    for hgx, omes, score in gbc_res:
        hgx2omes2gbc[hgx][omes] = score

    return hgx2omes2gbc    


def dndsGeneGrab(
    gff_path, assembly_path, proteome_path,
    ome_sig_clus, gene2hg, clusplusminus
    ):

    gff_list, prot_dict = gff2list(gff_path), fa2dict(proteome_path)
    assem_dict, clus_out = fa2dict(assembly_path), {}
    cds_dict, cds_dict2 = compileCDS2(
        gff_list, os.path.basename(gff_path).replace('.gff3','')
        )

    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]):
            if seq0 in ome_sig_clus: # if the og is part of a significant seed
            # locus
                try:
                    og0 = gene2hg[seq0]
                except KeyError:
                    continue
                for hit in ome_sig_clus[seq0]:
                    hgx = tuple(sorted(list(hit)))
                    if hgx not in clus_out:
                        clus_out[hgx] = [{og: [{}, {}] for og in hgx}]
                    clus_out[hgx][-1][og0][0][seq0] = prot_dict[seq0]
                    clus_out[hgx][-1][og0][1] = {
                        **clus_out[hgx][-1][og0][1], 
                        **ntmain(cds_dict2[seq0], assem_dict)
                        }
                    
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
                        if og in hit: # if the og is
                        # in the sig clus
                            clus_out[hgx][-1][og][0][seq1] = prot_dict[seq1]
                            clus_out[hgx][-1][og][1] = {
                                **clus_out[hgx][-1][og][1],
                                **ntmain(cds_dict2[seq1], assem_dict)
                                }
    return clus_out


def dnds_preparation(db, omes4hgs, hgx2omes, gene2hg, plusminus, hgx_dir, i2ome):

    gene_checks, hgx_check = {}, {}
    
    for ome in omes4hgs:
        gff_path = db[i2ome[ome]]['gff3'] 
        proteome_path = db[i2ome[ome]]['faa']
        assembly_path = db[i2ome[ome]]['fna']
        res = dndsGeneGrab(
            gff_path, assembly_path, proteome_path,
            omes4hgs[ome], gene2hg, plusminus
            )
        for hgx in res:
            if hgx not in hgx_check:
                hgx_check[hgx] = set()
            if hgx not in gene_checks:
                gene_checks[hgx] = {og: [{}, {}] for og in hgx}
            for loc in res[hgx]:
                for og in loc:
                    gene_checks[hgx][og][0] = {
                        **gene_checks[hgx][og][0], **loc[og][0]
                        }
                    gene_checks[hgx][og][1] = {
                        **gene_checks[hgx][og][1], **loc[og][1]
                        }
            hgx_check[hgx].add(ome)
            if hgx_check[hgx] == set(hgx2omes[hgx]):
                for og in hgx:
                    out_base = hgx_dir + '-'.join([str(x) for x in hgx]) + \
                        '.' + str(og)
                    with open(out_base + '.aa.fa', 'w') as out:
                        out.write(dict2fa(gene_checks[hgx][og][0]))
                    with open(out_base + '.nt.fa', 'w') as out:
                        out.write(dict2fa(gene_checks[hgx][og][1]))
                del gene_checks[hgx]

    if gene_checks:
        print('\t\tERROR: discrepancies with previous run', flush = True)
        for hgx in gene_checks:
            print('\t\t\t' + str(hgx), flush = True)
            for og in hgx:
                out_base = hgx_dir + '-'.join([str(x) for x in hgx]) + \
                    '.' + str(og)
                with open(out_base + '.aa.fa', 'w') as out:
                    out.write(gene_checks[hgx][og][0])
                with open(out_base + '.nt.fa', 'w') as out:
                    out.write(gene_checks[hgx][og][1])
        

def calc_dnds(mafft):

    og_catch = re.search(r'(^[^\.]+)\.(\d+)\.mafft$', os.path.basename(mafft))
    hgx = tuple([int(x) for x in og_catch[1].split('-')])
    og = int(og_catch[2])

    err = os.path.basename(mafft).replace('.mafft','')
    nt_path = mafft.replace('.mafft', '.nt.fa')
    try:
        prot_align = AlignIO.read(mafft, 'fasta')
    except ValueError:
        print('\t\t' + err + ' empty alignment?', flush = True)
        return None

    nt_fa = SeqIO.parse(nt_path, 'fasta')

    try:
        aln = codonalign.build(prot_align, nt_fa) #, ids)
    except KeyError:
        print('\t\t' + err + ' ambiguous AAs', flush = True)
        return
    except RuntimeError:
        print('\t\t' + err + ' BioPython error', flush = True)
        return
    except IndexError:
        print('\t\t' + err + ' unknown error', flush = True)
        return
    try:
        dn_matr, ds_matr = aln.get_dn_ds_matrix()
        np_dn, np_ds = np.array(dn_matr), np.array(ds_matr)
        dnds_matrNulls = np.divide(np_dn, np_ds)
        dnds_matr = np.ma.array(dnds_matrNulls, mask=np.isnan(dnds_matrNulls))
        dnds_matr = np.abs(dnds_matr)
        dnds = np.mean(dnds_matr)
        return hgx, og, dnds
    except ZeroDivisionError:
        print('\t\t' + err + ' unknown error')
    except KeyError:
        print('\t\t' + err + ' ambiguous NTs')
    return None


def parse_dnds(mpRes):

    hgxdNdS_dict = {}
    for res in mpRes:
        if res:
            hgx, og, dnds = res[0], res[1], res[2]
            try:
                float(dnds)
            except ValueError:
                continue
            if hgx not in hgxdNdS_dict:
                hgxdNdS_dict[hgx] = {
                    'selection_total': 0,
                    'og_total': 0,
                    'hgs': {}
                    }

            if dnds < 1: # if purifying selection
                hgxdNdS_dict[hgx]['selection_total'] += dnds
            else: # if positive selection use reciprocal
                hgxdNdS_dict[hgx]['selection_total'] += 1/dnds
            hgxdNdS_dict[hgx]['og_total'] += 1
            hgxdNdS_dict[hgx]['hgs'][og] = dnds
    for hgx in hgxdNdS_dict:
        selection_coefficient = hgxdNdS_dict[hgx]['selection_total']/hgxdNdS_dict[hgx]['og_total']
        mean_dnds = sum(hgxdNdS_dict[hgx]['hgs'].values())/len(hgxdNdS_dict[hgx]['hgs'])
        hgxdNdS_dict[hgx] = [selection_coefficient, mean_dnds, hgxdNdS_dict[hgx]['hgs']]

    return hgxdNdS_dict


def hash_clusters(ome, gff_path, ome_sig_clus, gene2hg, clusplusminus):
    """sliding window of size clusplusminus, identify sequences that may have
    hgs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus"""

    gff_list, protoclus, clus_out = gff2list(gff_path), {}, []
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))

    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]):
            if seq0 in ome_sig_clus: # if the og is part of a significant seed
            # locus
                for sig_clus in ome_sig_clus[seq0]:
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
                    hgx = tuple(sorted(list(sig_clus[0])))
                    clus_out.append([
                         hgx, locus[start:end], {sig_clus[1]}, scaf
                         ])
                    # clus_out = [(og0, og1, ... ogn), [prot0, prot1, ... protn])]
    
    return ome, clus_out, cds_dict


def merge_clusters(clus_out):
    """because we have significant clusters with borders, we can merge each
    cluster and assume that at a certain number of hgs in a combo will never
    happen via random sampling; therefore, their significance in the
    microsynteny adjustment is sufficient for justification for merging"""

    comp_clus, change = [], False
    while clus_out: # while there are clusters
        any_intersect, toDel = False, [0] # we haven't found intersection
        clus0 = clus_out[0] # grab the first tuple
        loc0 = set(clus0[1]) # grab the set of proteins in the cluster
        for i, clus1 in enumerate(clus_out[1:]): # look at all other clusters
            loc1 = set(clus1[1]) # grab their loci
            intersect = loc0.intersection(loc1) # check for overlap between
            # proteins
            if intersect:
                any_intersect, change = True, True # we found overlap in
                # clus_out and in this specific cluster comparison
                hgx = tuple(sorted(list(set(clus0[0]).union(set(clus1[0])))))
                # obtain the higher order og combination as the union of the
                # two cluster sets, then sort it and store as a tuple
                comp_clus.append([
                    hgx, list(loc0.union(loc1)), 
                    clus0[2].union(clus1[2]), clus0[3] #, merge_x
                    ])
                # append the og combo and the union of the two loci's proteins
                # protein order is random
                toDel.append(i) # we can delete the intersected locus as well
                # as the first

        if not any_intersect: # if there is not any intersection in the locus
            comp_clus.append(clus0) #, clus0[2]]) # simply append the original
            # result
        toDel.sort(reverse = True) # sort the todel from highest to lowest to
        # not mess up order when indexing
        for i in toDel:
            clus_out.pop(i)

    return comp_clus, change


def write_clusters(ome_sig_clus, ome, out_file, gff_path, gene2hg, clusplusminus = 10):
    # is multiprocessing this messing something and memory and that's why
    # output is wonky? (sometimes fewer proteins than should be in locus) 
    # or is output consistent fucking up elsewhere

    ome, clus_out, cds_dict = hash_clusters(ome, gff_path, ome_sig_clus, gene2hg, clusplusminus)
    change = True
    while change: # so long as there is a single merge
        clus_out, change = merge_clusters(clus_out) # attempt to merge clusters
    clus_out = sorted(clus_out, key = lambda x: len(x[0]), reverse = True)
    # sort the clusters by the size of the OG combination tuple highest to
    # lowest

    scafs = set()
    for clus in clus_out:
        clus[1] = sorted(clus[1], key = cds_dict[clus[3]].index)
        clusOGs = []
        for gene in clus[1]:
            try:
                clusOGs.append(gene2hg[gene])
            except KeyError: # gene not in an OG/is a singleton
                clusOGs.append('')
        clus[0] = clusOGs 
        scaf = clus[3][:20]
        if clus[3] != scaf:
            scaf = clus[3][:20] + '..'
        count = 0
        while scaf + '_' + str(count) in scafs:
            count += 1
        name = scaf + '_' + str(count)
        scafs.add(name)
        clus.append(name)
 

    out_genes = []
    with open(out_file, 'w') as out:
        out.write('#name\thgs\tgenes\tclans\n')
        for clus in clus_out:
            clus[2] = sorted(clus[2])
            hgs = ','.join([str(x) for x in clus[0]]) # og1,og2,...ogn
            genes = ','.join(clus[1]) # prot1,prot2,prot3,...protn
            clans= ';'.join([str(x) for x in list(clus[2])])
            out.write(clus[4] + '\t' + hgs + '\t' + genes + '\t' + clans + '\n') 
            out_genes.append(clus[1])

    return ome, out_genes


def write_cluster_scores(out_dir, clus_probs, i2ome, hgx2omes): #, rand_clus_probs):
    
    out_list, clus_len_dict = [], {}
    for hgx in clus_probs:
        out_list.append([
            hgx, clus_probs[hgx], 
            [i2ome[z] for z in hgx2omes[hgx]] #, rand_clus_probs[hgx]
            ])
#        clus_len_dict[hgx] = len(hgx2omes[hgx])
         
    out_list = sorted(out_list, key = lambda x: len(x[0]))

    with gzip.open(out_dir + 'clusters.probs.tsv.gz', 'wt') as out:
        out.write('#orthogroups\tmicrosynt_prob\tomes\n') #\trandom_prob\n')
        for i in out_list:
            out.write(','.join([str(x) for x in i[0]]) + '\t')
            out.write(str(i[1]) + '\t' + ','.join(i[2]) + '\n')#'\t' + str(i[3]) + '\n')


def setupHmmsearch(genes_in, prot_paths, hmm_dir):

    genes_out, todel = {}, []
    for ome in genes_in:
        prot_dict = fa2dict(prot_paths[ome])
        try:
            for sub_list in genes_in[ome]:
                for gene in sub_list:
                    genes_out[gene] = prot_dict[gene]
        except KeyError:
            eprint('\t\tERROR: ' + ome + ' accessions modified from input', flush = True)
            todel.append(ome)
            continue

    genes = list(genes_out.keys())
    for ome in todel:
        omeKeys = [x for x in genes if x.startswith(ome + '_')]
        for key in omeKeys:
            del genes_out[key]

    with open(hmm_dir + 'complete.fa', 'w') as out:
        out.write(dict2fa(genes_out))

    return list(genes_out.keys()), hmm_dir + 'complete.fa', todel


def runHmmsearch(pfam, fa, hmm_dir, cpus):
    hmmsearch = subprocess.call([
        'hmmsearch', '--cpu', str(cpus), '--domtblout', hmm_dir + 'pfam.tmp', pfam, fa
        ], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
        )
    if hmmsearch:
        print('\t\tERROR: hmmsearch failed: ' + str(hmmsearch), flush = True)
    else:
       shutil.move(hmm_dir + 'pfam.tmp', hmm_dir + 'pfam.out')
    

def parseHmmRes(hmm_out, evalue = 0.01, threshold = 0.5):

    lineComp, res = re.compile(r'([^ ]+)'), defaultdict(list)
    with open(hmm_out, 'r') as raw:
        for line in raw:
            if not line.startswith('#'):
                data = lineComp.findall(line.rstrip())
                hit_perc = (int(data[15]) - int(float(data[14])))/int(data[5])
                if hit_perc > threshold:
                    if float(data[6]) < evalue: # and float(data[7]) > score:
                        # pfam, name, cov_perc, evalue, bit, algn_start, algn_end
                        res[data[0]].append((
                            data[4], data[3], hit_perc, float(data[6]), 
                            float(data[7]), int(data[16]), int(data[17])
                            ))

    ome_res = defaultdict(dict)
    for gene in res:
        res[gene] = sorted(res[gene], key = lambda x: x[3], reverse = True)
        todel, hits = [], set()
        for i, hit in enumerate(res[gene]):
            if hit[0] in hits:
                todel.append(i)
            hits.add(hit)
        for i in reversed(todel):
            del res[gene][i]
        ome = gene[:gene.find('_')]
        ome_res[ome][gene] = res[gene]

    return ome_res


def pfamMngr(genes_list, prot_paths, wrk_dir, pfam, evalue = 0.01, threshold = 0.5, cpus = 1):

    hmm_dir = wrk_dir + 'hmm/'
    if not os.path.isdir(hmm_dir):
        os.mkdir(hmm_dir)

    genes, hmm_fa, failedOmes = setupHmmsearch(genes_list, prot_paths, hmm_dir)
    if not os.path.isfile(hmm_dir + 'pfam.out'):
        print("\tHmmsearch'ing Pfam database", flush = True)
        runHmmsearch(pfam, hmm_fa, hmm_dir, cpus)
    print('\tParsing hmmsearch output', flush = True)
    hmm_res = parseHmmRes(hmm_dir + 'pfam.out', evalue, threshold)
    
    return hmm_res, set(failedOmes)


def grabClus(genes_list, gff_path, prot_path, ome, ome_dir, gene2hg, pfamRes = {}):

    gff_dir, fa_dir = ome_dir + ome + '/gff/', ome_dir + ome + '/fa/'
    if not os.path.isdir(gff_dir):
        os.mkdir(gff_dir)
    if not os.path.isdir(fa_dir):
        os.mkdir(fa_dir)

    clus_gffs, clus_fas, clus_anns, scafs, svg_dict = [], [], {}, set(), {}
    gff_list, prot_dict = gff2list(gff_path), fa2dict(prot_path)
    print('\t' + ome, flush = True)
    for genes in genes_list:
        clus_gffs.append([[]])
        clus_fas.append({})
        clus_ann = ''
        for gene in genes:
            geneCoord = None
            geneGff = grabGffAcc(gff_list, gene)
            for entry in geneGff:
                if entry['type'].lower() == 'gene':
                    geneCoord = (int(entry['start']), int(entry['end']), entry['strand'])
            geneFa = prot_dict[gene]
            if gene in pfamRes:
                pfamStr = ';Pfam=' + '|'.join([
                    (hit[0] + '-' + hit[1]).replace('|','&') for hit in pfamRes[gene]
                    ])
                clus_ann += pfamStr[6:]
            else:
                pfamStr = ''
            clus_ann += ','
            try:
                ogStr = ';OG=' + str(gene2hg[gene])
            except KeyError: # gene not in an OG or is a singleton
                ogStr = ';OG='
            for entry in geneGff:
                if entry['type'] == 'gene' or entry['type'].lower() == 'mrna':
                    entry['attributes'] += pfamStr + ogStr
            geneFa['description'] = pfamStr[6:]
            clus_gffs[-1][0].extend(geneGff)
            clus_fas[-1][gene] = geneFa
        scaf_prep = clus_gffs[-1][0][0]['seqid']
        scaf = scaf_prep[:20]
        if scaf_prep != scaf:
            scaf = scaf_prep[:20] + '..'
        count = 0
        while scaf + '_' + str(count) in scafs:
            count += 1
        name = scaf + '_' + str(count)
        scafs.add(name)
        clus_anns[name] = clus_ann[:-1]
        clus_gffs[-1].append(name)
#        svg_dict[name] = copy.deepcopy(t_svg_dict)
        with open(gff_dir + name + '.gff3', 'w') as out:
            out.write(list2gff(clus_gffs[-1][0]))
        with open(fa_dir + name + '.fa', 'w') as out:
            out.write(dict2fa(clus_fas[-1]))

    with open(ome_dir + ome + '/info.out.tmp', 'w') as out:
        out.write('#name\thgs\tgenes\thgx_clans\tpfams\n')
        with open(ome_dir + ome + '/info.out', 'r') as raw:
            for line in raw:
                if not line.startswith('#'):
                    clus = line.split('\t')[0]
                    out.write(line.rstrip() + '\t' + clus_anns[clus] + '\n')
    shutil.move(ome_dir + ome + '/info.out.tmp', ome_dir + ome + '/info.out')

    return ome#, svg_dict
#    return ome, clus_gffs, clus_fas, svg_dict


def mk_3d(labels, axes, axes_labels, alpha = 0.6):

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        mode = 'markers',
        x = axes[0], y = axes[1], z = axes[2],
        marker = dict(
            color = 'rgba(163, 25, 169,' + str(alpha) + ')',
            size = 5, opacity = alpha
            ),
        text = labels, showlegend = False
        )
    )

    return fig


def mk_subplots(labels, axes, axes_labels, alpha = 0.6):

    colors = [
        [240,163,255], [0,117,220], [153,63,0], [76,0,92], [25,25,25], 
        [0,92,49], [43,206,72], [255,204,153], [128,128,128], [148,255,181], 
        [143,124,0], [157,204,0], [194,0,136], [0,51,128], [255,164,5], 
        [255,168,187], [66,102,0], [255,0,16], [94,241,242], [0,153,143], 
        [224,255,102], [116,10,255], [153,0,0], [255,255,128], [255,255,0], [255,80,5]
        ]
#    random.shuffle(colors)
    pl_clr = [
        'rgba(' + ', '.join([str(y) for y in x]) + ', ' + str(alpha) + ')' \
        for x in colors
        ]
    fig = make_subplots(
        rows=len(axes) - 1, cols=len(axes) - 1,
        start_cell = "bottom-left"
        )

    clr_count = 0
    for i0, axis0 in enumerate(axes):
        for i1, axis1 in enumerate(axes[i0+1:]):
            fig.add_trace(
                go.Scatter(
                    mode = 'markers', x=axis0, y=axis1,
                    marker = dict(
                        color = pl_clr[clr_count],
                        size = 5,
                        opacity = alpha
                        ),
                    text = labels
                    ), row = i0+1, col = i0+i1+1,
                )
            fig.update_xaxes(
                title_text=axes_labels[i0], row = i0+1, col =
                i0+i1+1
                )
            fig.update_yaxes(
                title_text=axes_labels[i0+i1+1], row = i0+1, col =
                i0+i1+1
                )
            clr_count += 1

    return fig


def outputSVG(clus, svg_dict, svg_dir, width):

    colors = [
        '#000000', '#010067', '#d5ff00', '#ff0056', '#9e008e', '#0e4ca1', 
        '#ffe502', '#005f39', '#00ff00', '#95003a', '#ff937e', '#a42400', 
        '#001544', '#91d0cb', '#620e00', '#6b6882', '#0000ff', '#007db5', 
        '#6a826c', '#00ae7e', '#c28c9f', '#be9970', '#008f9c', '#5fad4e', 
        '#ff0000', '#ff00f6', '#ff029d', '#683d3b', '#ff74a3', '#968ae8', 
        '#98ff52', '#a75740', '#01fffe', '#ffeee8', '#fe8900', '#bdc6ff', 
        '#01d0ff', '#bb8800', '#7544b1', '#a5ffd2', '#ffa6fe', '#774d00', 
        '#7a4782', '#263400', '#004754', '#43002c', '#b500ff', '#ffb167', 
        '#ffdb66', '#90fb92', '#7e2dd2', '#bdd393', '#e56ffe', '#deff74', 
        '#00ff78', '#009bff', '#006401', '#0076ff', '#85a900', '#00b917', 
        '#788231', '#00ffc6', '#ff6e41', '#e85ebe'
        ]

    features, color_dict = [], {}
    geneMin = min([svg_dict[x]['gene'][0] for x in svg_dict])
    geneMax = max([svg_dict[x]['gene'][1] for x in svg_dict])
    for gene in svg_dict:
        feature = GraphicFeature(
            start=svg_dict[gene]['gene'][0] - geneMin, 
            end=svg_dict[gene]['gene'][1] - geneMin,
            color='#ffffff', strand=int(svg_dict[gene]['gene'][2] + '1'),
            label=gene
            )
        features.append(feature)
    colorI = 0
    for gene in svg_dict:
        for pfam in svg_dict[gene]['pfam']:
            if pfam[0] not in color_dict:
                color_dict[pfam[0]] = colors[colorI]
                colorI += 1
                if colorI >= len(colors):
                    colorI = 0
            feature = GraphicFeature(
                start=pfam[2] - geneMin, end=pfam[3] - geneMin,
                label=pfam[0], strand=int(svg_dict[gene]['gene'][2] + '1'),
                color=color_dict[pfam[0]]
                )
            features.append(feature)

    record = GraphicRecord(
        sequence_length = geneMax - geneMin, features = features
        )
    ax, _ = record.plot(figure_width = width)
    ax.figure.savefig(svg_dir + clus + '.svg')


def runGFF2SVG(ome_dir, regex = r'Pfam=[^;]+'):
    if not os.path.isdir(ome_dir + 'svg/'):
        os.mkdir(ome_dir + 'svg/')
    gffs = collect_files(ome_dir + 'gff/', 'gff3')
    for gff in gffs:
        svg_path = ome_dir + 'svg/' + re.sub(r'\.gff3$', '.svg', os.path.basename(gff))
        gff2svg(gff2list(gff), svg_path, prod_comp = regex, types = types)


def initLog(
    log_file, log_dict
    ):
    with open(log_file, 'w') as out:
        out.write(
            'hg_file\t' + str(log_dict['hg_file']) + '\n' + \
            'plusminus\t' + str(log_dict['plusminus']) + '\n' + \
            'pair_percentile\t' + str(log_dict['pair_percentile']) + '\n' + \
            'hgx_percentile\t' + str(log_dict['hgx_percentile']) + '\n' + \
            'border_percentile\t' + str(log_dict['border_percentile']) + '\n' + \
            'null_samples\t' + str(log_dict['null_samples']) + '\n' + \
            'n50\t' + str(log_dict['n50'])
            )

def readLog(
    log_file, log_dict
    ):
    log_res = {}
    with open(log_file, 'r') as raw:
        for line in raw:
            key = line[:line.find('\t')]
            res = line[line.find('\t') + 1:].rstrip()
            if res != str(log_dict[key]):
                log_res[key] = False
            else:
                log_res[key] = True
    try:
        if not log_res['n50']:
            log_res['plusminus'] = False
        if not log_res['plusminus']:
            log_res['null_samples'] = False
        if not log_res['null_samples']:
            log_res['pair_percentile'] = False
        if not log_res['pair_percentile']:
            log_res['border_percentile'] = False
        if not log_res['border_percentile']:
            log_res['hgx_percentile'] = False
    except KeyError:
        print('\nERROR: corrupted log.txt.' + \
            '\nIf not rectified, future runs will completely overwrite the current\n')
        sys.exit(149)

    return log_res


def rmOldData(
    log_res, out_dir, wrk_dir
    ):
    if not log_res['null_samples']:
        nulls = collect_files(wrk_dir, 'null.txt')
        for null in nulls:
            os.remove(null)
    if not log_res['pair_percentile']:
        seed_file = out_dir + 'seed_scores.tsv.gz'
        seed_arr = wrk_dir + '.arr.npy'
        if os.path.isfile(seed_file):
            os.remove(seed_file)
        if os.path.isfile(seed_arr):
            os.remove(seed_arr)
        clus_pickle = wrk_dir + 'hgx2loc.pickle'
        hgx_pickle = wrk_dir + 'hgx_scores.pickle'
        ome_pickle = wrk_dir + 'hgx_omes.pickle'
        if os.path.isfile(clus_pickle):
            os.remove(clus_pickle)
        if os.path.isfile(hgx_pickle):
            os.remove(hgx_pickle)
        if os.path.isfile(ome_pickle):
            os.remove(ome_pickle)
    if not log_res['border_percentile']:
        row_file = wrk_dir + 'mtx/mcl.prep.rows'
        prep_file = wrk_dir + 'mtx/mcl.prep.gz'
        if os.path.isfile(row_file):
            os.remove(row_file)    
        if os.path.isfile(prep_file):
            os.remove(prep_file)
    if not log_res['hgx_percentile']:
        kern_file = out_dir + 'hgx_clans.tsv.gz'
        clus_file = out_dir + 'hgxs.tsv.gz'
        patch_pickle = wrk_dir + 'patchiness.scores.pickle'
        ome_dir = wrk_dir + 'ome/'
        hgx_dir = wrk_dir + 'hgx/'
        hmm_dir = wrk_dir + 'hmm/'
        if os.path.isfile(kern_file):
            os.remove(kern_file)
        if os.path.isfile(clus_file):
            os.remove(clus_file)
        if os.path.isdir(ome_dir):
            shutil.rmtree(ome_dir)
        if os.path.isdir(hgx_dir):
            shutil.rmtree(hgx_dir)
        if os.path.isdir(hmm_dir):
            shutil.rmtree(hmm_dir)
        if os.path.isfile(wrk_dir + 'hgx.tar.gz'):
            os.remove(wrk_dir + 'hgx.tar.gz')
        if os.path.isfile(patch_pickle):
            os.remove(patch_pickle)
    if not log_res['hg_file']:
        shutil.rmtree(wrk_dir)
        os.mkdir(wrk_dir)
        for key in log_res:
            log_res[key] = False


def logCheck(log_dict, log_path, out_dir, wrk_dir, flag = True):

    if not os.path.isfile(log_path):
        log_res = {x: False for x in log_dict}
        rmOldData(log_res, out_dir, wrk_dir)
        initLog(log_path, log_dict)
    log_res = readLog(log_path, log_dict)
    if any(not log_res[x] for x in log_res):
        if not flag:
            print('\nInitializing new run', flush = True)
            rmOldData(log_res, out_dir, wrk_dir)
            initLog(log_path, log_dict)
        else:
            eprint('\nERROR: -n not specified and run parameters changed \
                    \nSee ' + log_path, flush = True)
            sys.exit(15)

    return log_res


def patch_main(
    phylo, hgx2omes, hgxs, wrk_dir, 
    old_path = 'patchiness.scores.pickle', cpus = 1
    ):

    if not os.path.isfile(wrk_dir + old_path):
        if isinstance(hgx2omes, dict): # round 1 patchiness
            clusOmes = set([tuple([str(x) for x in hgx2omes[y]]) for y in hgxs])
        else: # round 2 is a list
            clusOmes = set([
                tuple([str(x) for x in y]) for y in hgxs
                ])
#        more = set([tuple(omes) for omes in moduleOmes])
#        allHGxs = list(clusOgxs.union(more))    
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            patch_res = pool.starmap(
                calc_patchiness, [(phylo, x) for x in clusOmes]
                )
        pool.join()
        omes2patch = {ome_tup: patchiness for ome_tup, patchiness in patch_res}

        with open(wrk_dir + old_path, 'wb') as out:
            pickle.dump(omes2patch, out)

    else:
        print('\tLoading previous patchiness results', flush = True)
        with open(wrk_dir + old_path, 'rb') as in_pick:
            omes2patch = pickle.load(in_pick)
 
    return omes2patch

def output_og_fas(db, genes, hg_file):
    fa_str = dict2fa(acc2fa(
                db, genes, error = False, spacer = '\t\t',
                coord_check = False
                ))
    with open(hg_file, 'w') as out:
        out.write(fa_str)


def gbc_main(
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
        if not modules: # for hgxs
            d2gbc = gbc_mngr_2(
                list(hgs), list(ome2i.keys()), hg_dir, hgx_dir, blastp, hgx2loc,
                db, gene2hg, plusminus, hg2gene, cpus = cpus, 
                printexit = printexit
                )
        else: # run the kernel detection version
            d2gbc = gbc_mngr_3(
                list(hgs), list(ome2i.keys()), hgx_dir, hgx2loc, 
                db, gene2hg, plusminus, hg2gene, modules,
                moduleOmes, moduleHGxs, ome2i, cpus = cpus
                )
            hgx_dirTar = mp.Process(target=tardir, args=(hgx_dir, True))
            hgx_dirTar.start() # when to join ...
        with open(wrk_dir + old_path, 'wb') as pickout:
            pickle.dump(d2gbc, pickout)
    else:
        print('\tLoading previous coevolution results', flush = True)
        with open(wrk_dir + old_path, 'rb') as pickin:
            d2gbc = pickle.load(pickin)

    return d2gbc


def main(
    db, hg_file, out_dir, plusminus = 1, hg_dir = None,
    seed_perc = 0.2, clus_perc = 0.7, hgx_perc = 0.7,
    minimum_omes = 2, samples = 10000, pfam = None,
    constraint_path = None, blastp = 'blastp',
    run_dnds = False, cpus = 1, n50thresh = None, 
    root = None, coevo_thresh = 0, patch_thresh = 0,
    microsyn_thresh = 0, method = 'mmseqs easy-cluster',
    printexit = False, flag = True
    ):
    """
    The general workflow:
    log management -> input data parsing -> orthogroup pair identification ->
    microsynteny distance thresholding -> HGx formation ->
    microsynteny distance border thresholding -> HGx grouping ->
    microsynteny distance cluster thresholding -> patchiness calculation ->
    coevolution calculation -> optional dN/dS calculations ->
    HGx data output -> cluster retrieving -> data output
    """

    # initialize the log and working directory
    wrk_dir = out_dir + 'working/'
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)
    log_path = out_dir + 'log.txt'
    log_dict = {
        'hg_file': hg_file, 'plusminus': plusminus,
        'pair_percentile': seed_perc, 'hgx_percentile': clus_perc,
        'border_percentile': hgx_perc,
        'null_samples': samples, 'n50': n50thresh
        }
    log_res = logCheck(log_dict, log_path, out_dir, wrk_dir, flag)

    # if not using an external tree
#    if constraint_path:
    tree_path = out_dir + 'microsynt.newick'

    # obtain useable omes
    useableOmes, dbOmes = set(), set(db['ome'])
    print('\nI. Inputting data', flush = True)
    if n50thresh: # optional n50 threshold via mycotoolsDB
        assemblyPath = format_path('$MYCODB/../data/assemblyStats.tsv')
        if os.path.isfile(assemblyPath):
            with open(assemblyPath, 'r') as raw:
                for line in raw:
                    d = line.rstrip().split('\t')
                    ome, omeN50 = d[1], float(d[2])
                    if omeN50 > n50thresh and ome in dbOmes:
                        useableOmes.add(ome)
        else:
            raise FileNotFoundError(assemblyPath + ' not found. N50 threshold not applied')
    else:
        useableOmes = dbOmes

    # initialize orthogroup data structures            
    if not orthogroups and not os.path.isfile(wrk_dir + 'homology_groups.tsv'):
        hg_file = run_mmseqs(db, wrk_dir, algorithm = method,
                               min_id = 0.3, cpus = cpus)
    print('\tParsing orthogroups', flush = True)
    ome2i, gene2hg, i2ome, hg2gene = compile_homolog_groups(hg_file, wrk_dir, 
                                                            method, useableOmes)
    print('\t\tOmes:', len(ome2i), flush = True)
    print('\t\tHomology groups:', len(hg2gene), flush = True)
    print('\t\tGenes:', len(gene2hg), flush = True)

    with open(wrk_dir + 'og_info.txt', 'w') as out:
        out.write('\n'.join([str(k) + '\t' + ' '.join(v) \
                  for k, v in hg2gene.items()]))

    cc_arr_path = wrk_dir + ''
    if not os.path.isfile(out_dir + 'seed_scores.tsv.gz'):
        # compile cooccuring pairs of orthogroups in each genome
        print('\tCompiling all loci', flush = True)
        pairs = compile_loci(
            db, ome2i, gene2hg, plusminus, 
            cpus = cpus
            )
    
        # assimilate cooccurrences across omes
        print('\tIdentifying cooccurences', flush = True) 
    
        seed_len = sum([len(pairs[x]) for x in pairs])
        print('\t\t' + str(seed_len) + ' initial OG-pairs', flush = True)
        cooccur_array, cooccur_dict, hgpair2i, i2hgpair = \
            form_cooccur_structures(pairs, 2, len(ome2i), cc_arr_path)
        max_ome = max([len(cooccur_dict[x]) for x in cooccur_dict])
        print('\t\t' + str(max_ome) + ' maximum organisms with OG-pair', flush = True)
        cooccur_array[cooccur_array > 0] = 1
        cooccur_array.astype(np.int8)
        print('\t\t' + str(sys.getsizeof(cooccur_array)/1000000) + ' MB', flush = True)
        cooccur_array, del_omes = remove_nulls(cooccur_array)
        for i in del_omes:
            del i2ome[i]

        ome2i = {v: i for i, v in enumerate(i2ome)}
        with open(wrk_dir + 'ome2i.tsv', 'w') as out:
            out.write(
                '\n'.join([k + '\t' + str(v) for k, v in ome2i.items()])
                )
        np.save(cc_arr_path, cooccur_array)
#    elif not os.path.isfile(tree_path):
    else:
        cooccur_array = np.load(cc_arr_path + '.npy')

    # reload ome2i to make sure it is up-to-date with lost data
    ome2i, i2ome = {}, []
    with open(wrk_dir + 'ome2i.tsv', 'r') as raw:
        for line in raw:
            ome, i = line.rstrip().split()
            ome2i[ome] = int(i)
            i2ome.append(ome) # dependent on ome2i being sorted in output
    print('\t' + str(len(ome2i)) + ' omes with shared homolog combinations',
          flush = True)


    db = db.set_index('ome')
    microsynt_dict = {}
    if not os.path.isfile(tree_path):
        # create microsynteny distance matrix and make tree
        print('\tPreparing microsynteny alignment', flush = True)
#        microsynt_dict = prep_microsynt_dict(cooccur_array)
        align_file = align_microsynt_np(cooccur_array, i2ome, hg2gene, 
                                        hgpair2i, wrk_dir)
        print('\tBuilding microsynteny tree', flush = True)
        run_tree(align_file, wrk_dir, constraint = constraint_path, iqtree = 'iqtree',
                 model = 'GTR2+FO+ASC+R5', verbose = False, cpus = cpus)   

    print('\tReading microsynteny tree', flush = True)
    phylo = compileTree(
        i2ome, tree_path, root = root
        )

    # create null distribution for orthogroup pairs
    if not os.path.isfile(wrk_dir + '2.null.txt'):
        print('\tGenerating null distributions', flush = True)
        omes2dist, pair_null = gen_nulls(
            pairs, phylo, samples = samples, cpus = cpus
            )
        with open(wrk_dir + '2.null.txt', 'w') as out:
            out.write('\n'.join([str(i) for i in pair_null]))
        with open(wrk_dir + 'ome_scores.pickle', 'wb') as out:
            pickle.dump(omes2dist, out)
    else: # or just load what is available
        with open(wrk_dir + 'ome_scores.pickle', 'rb') as raw:
            omes2dist = pickle.load(raw)
        with open(wrk_dir + '2.null.txt', 'r') as raw:
            pair_null = [float(x.rstrip()) for x in raw]
    scores_i = round(seed_perc * len(pair_null) + .5)
    min_score = pair_null[scores_i]


    # seed clusters by calcuating total microsynteny distance for 
    # orthogroup pairs
    print('\nII. Seeding HG pairs', flush = True) 
    if not os.path.isfile(out_dir + 'seed_scores.tsv.gz'):
        print('\tCalculating seed OG-pair scores', flush = True)
        seed_score_start = datetime.now()

        results = calc_dists(phylo, cooccur_dict, cpus)
        omes2dist, top_hgs = {x[1]: x[0] for x in results}, []
        with open(wrk_dir + 'ome_scores.pickle', 'wb') as out:
            pickle.dump(omes2dist, out)
        for hgpair in cooccur_dict:
            score = omes2dist[cooccur_dict[hgpair]]
            if score > min_score:
                i = revs[hgpair]
                top_hgs.append([
                    [0], hgpair[1], len(cooccur_dict[hgpair]), score#, score/ome_dist
                    ])

        # write scores            
        print('\t\t' + str(datetime.now() - seed_score_start), flush = True)
        top_hgs = sorted(top_hgs, key = lambda x: x[3])
        print('\tWriting seed scores', flush = True)
        with gzip.open(out_dir + 'seed_scores.tsv.gz', 'wt') as out:
            out.write('#og0\tog1\tcooccurrences\tscore\tadj_score\n')
            for line in top_hgs:
                out.write('\t'.join([str(x) for x in line]) + '\n')
        print('\t\t' + str(len(top_hgs)) + ' significant seeds', flush = True)
    
    elif not os.path.isfile(wrk_dir + 'hgx_scores.pickle'): # load previous og pairs
        print('\tLoading previous seed HG pairs', flush = True)
        top_hgs = load_seedScores(out_dir + 'seed_scores.tsv.gz', min_score)
        print('\t\t' + str(len(top_hgs)) + ' significant seeds', flush = True)


    # begin sifting for HGxs using pairs as seeds for HGx detection
    print('\nIII. Sprouting high order OG combinations (HGx)', flush = True)
    if not os.path.isfile(wrk_dir + 'hgx_omes.pickle'):
        hgpair_dict = form_hgpairDict(top_hgs)
        print('\tForming HGxs', flush = True)
        form_clus_start = datetime.now()
        hgx2omes, hgx2loc = hgpair2hgx(
            db, hgpair_dict, gene2hg, 
            ome2i, cpus, clusplusminus = plusminus
            )
        with open(wrk_dir + 'hgx2loc.pickle', 'wb') as pickout:
            pickle.dump(hgx2loc, pickout)
        with open(wrk_dir + 'hgx_omes.pickle', 'wb') as pickout:
            pickle.dump(hgx2omes, pickout)
        formHGxTime = datetime.now() - form_clus_start
        print('\t\t' + str(formHGxTime), flush = True)

    else: # or just load available structures
        with open(wrk_dir + 'hgx2loc.pickle', 'rb') as pickin:
            hgx2loc = pickle.load(pickin)
        with open(wrk_dir + 'hgx_omes.pickle', 'rb') as pickin:
            hgx2omes = pickle.load(pickin)

    ome_combos = set([tuple(sorted(list(x))) for x in list(hgx2omes.values())])
    if not os.path.isfile(wrk_dir + 'hgx_scores.pickle'):
        hgx2i, i2hgx, hgx_cooccur_dict = form_cooccur_dict(
            hgx2omes
            ) # create hgx data structures

        print('\tCalculating HGx microsynteny distances', flush = True)
        clus_score_start = datetime.now()
        clus_obs = len(hgx2omes)
        print('\t\t' + str(clus_obs) + ' observed HGx', flush = True)

        # calculate HGx microsynteny distances
        results = calc_dists(phylo, hgx_cooccur_dict, cpus, omes2dist = omes2dist)
        hgx2dist = {}
        omes2dist, top_hgs = {**omes2dist, **{x[1]: x[0] for x in results}}, []
        with open(wrk_dir + 'ome_scores.pickle', 'wb') as out:
            pickle.dump(omes2dist, out)
        for hgx in hgx_cooccur_dict:
            score = omes2dist[tuple(hgx_cooccur_dict[hgx])]
            hgx2dist[hgx] = score
    
        print('\t\t' + str(datetime.now() - clus_score_start), flush = True)
        with open(wrk_dir + 'hgx_scores.pickle', 'wb') as pickout:
            pickle.dump(hgx2dist, pickout)

    else: # just load previous hgx scores
        print('\tLoading previous HGxs', flush = True)
        with open(wrk_dir + 'hgx_scores.pickle', 'rb') as pickin:
            hgx2dist = pickle.load(pickin) 

    calc_hgx_p = False
    if not os.path.isfile(wrk_dir + 'hgx2pval.pickle'):
        if calc_hgx_p:
            print('\tEstimating HGx probability', flush = True)
            omes2hg2genes, omes2genes = {}, defaultdict(list)
            for gene, og in gene2hg.items():
                ome = gene[:gene.find('_')]
                omeI = ome2i[ome]
                if omeI not in omes2hg2genes:
                    omes2hg2genes[omeI] = defaultdict(list)
                omes2hg2genes[omeI][og].append(gene)
                omes2genes[omeI].append(gene)
            unadjHGx2pval = combo_prob_mngr(
                hgx2omes, omes2hg2genes, omes2genes, (plusminus*2)-1, 
                cooccur_array, cpus = cpus
                )
   #     comparisons = len(hgx2omes)
  #      hgx2pval = {
 #           k: v * comparisons for k, v in unadjHGx2pval.items()
#            } # apply bonferroni correction
            hgx2pval = unadjHGx2pval
            with open(wrk_dir + 'hgx2pval.pickle', 'wb') as out:
                pickle.dump(unadjHGx2pval, out)
    else:
        with open(wrk_dir + 'hgx2pval.pickle', 'rb') as raw:
            hgx2pval = pickle.load(raw)
    

    # prepare null distributions for each size (# of OGs) observed
    # in HGxs    
    max_hgx_size = max([len(x) for x in hgx2omes])
    if not os.path.isfile(wrk_dir + str(max_hgx_size) + '.null.txt'):
        print('\tPreparing HGx nulls', flush = True)
        bordScores, clusScores = genHGxNulls(
            [v['gff3'] for k, v in db.items() if k in ome2i], 
            gene2hg, max_hgx_size, plusminus, phylo,
            hgx_perc, clus_perc, wrk_dir,
            omes2dist, samples = samples, cpus = cpus
            )
    else: # or just load the old ones
        print('\tLoading HGx nulls', flush = True)
        bordScores, clusScores = loadNulls(max_hgx_size, wrk_dir, hgx_perc, clus_perc)

    # collect for normalizing relative to the whole dataset later
    # should I move into absolute space? it may just be better for comparing
    # datasets in the future and some algorithm is going to pick up that gene clusters
    # are in some absolute microsynteny distance, so long as its not normalized. 
#    max_dist = max(hgx2dist.values())
 #   min_dist = min(hgx2dist.values())

    i2hgx, hgx2i, thgx2dist, count = {}, {}, {}, 0
    for hgx in list(hgx2dist.keys()):
        if hgx2dist[hgx] >= bordScores[len(hgx)]:
             thgx2dist[hgx] = hgx2dist[hgx]
             i2hgx[count], hgx2i[hgx] = hgx, count
             count += 1
             # apply the border threshold for ALL pre-grouped hgxs
    print(
        '\t\t' + str(len(hgx2dist)) + ' HGx pass border threshold', 
        flush = True
        )
    hgx2dist = thgx2dist
    
    hgx2gbc, omes2patch = {}, {}
#    print('\nIV. Quantifying HGx patchiness', flush = True)
 #   omes2patch = patch_main(
  #      phylo, hgx2omes, hgx2dist, wrk_dir,
   #     old_path = 'patchiness.scores.pickle', cpus = cpus
    #    )

    hgx_dir = wrk_dir + 'hgx/'
    if not checkdir(hgx_dir, unzip = True, rm = True):
        os.mkdir(hgx_dir)

    print('\nIV. Quantifying HGx coevolution', flush = True)
    if not hg_dir:
        og_fa_cmds = []
        hg_dir = wrk_dir + 'og/'
        if not os.path.isdir(hg_dir):
            os.mkdir(hg_dir)
        for og, genes in hg2gene.items():
            hg_file = hg_dir + str(og) + '.faa'
            if not os.path.isfile(hg_file):
                og_fa_cmds.append([db, genes, hg_file])
        with mp.Pool(processes = cpus) as pool:
            pool.starmap(output_og_fas, og_fa_cmds)
#                fa_str = dict2fa(acc2fa(
 #                           db, genes, error = False, spacer = '\t\t'
  #                          ))
   #             with open(hg_file, 'w') as out:
    #                out.write(fa_str)

    hgx2gbc = gbc_main(
        hgx2loc, wrk_dir, ome2i, hg_dir, hgx_dir,
        'diamond', db, gene2hg, plusminus, hg2gene, 
        old_path = 'hgx2gbc.pickle', cpus = cpus, printexit = printexit
        ) # should be able to skip this

    # Group hgxs
    print('\nV. Inferring HGx clans', flush = True)
    if not os.path.isfile(wrk_dir + 'ocgs.pickle'): # need to add this to
    # log parsing
        ocgs, ocg_hgxs, ocg_omes, omes2dist = classify_ocgs(
            hgx2loc, db, gene2hg, i2hgx, hgx2i,
            phylo, clusScores, bordScores, ome2i,
            hgx2omes, wrk_dir, #hgx2dist,
            omes2dist = omes2dist, clusplusminus = plusminus, 
            min_omes = 2, cpus = cpus
            )
        with open(wrk_dir + 'ocgs.pickle', 'wb') as pickout:
            pickle.dump([ocgs, ocg_omes, ocg_hgxs], pickout)
        with open(wrk_dir + 'ome_scores.pickle', 'wb') as pickout:
            pickle.dump(omes2dist, pickout)
    else: # or just load old data
        with open(wrk_dir + 'ocgs.pickle', 'rb') as raw:
            ocgs, ocg_omes, ocg_hgxs = pickle.load(raw)

    runOmes = [
        omes for omes in ocg_omes \
        if omes not in omes2patch
        ] # omes without patchiness scores
    runHGxs = [
        hgx for i, hgx in enumerate(ocg_hgxs) \
        if ocg_omes[i] not in omes2patch
        ] # hgxs without patchiness scores
   
    print('\nVI. Quantifying OCG patchiness', flush = True)
    omes2patch = {**patch_main(
        phylo, runHGxs, runOmes, wrk_dir,
        old_path = 'patchiness.full.pickle', cpus = cpus
        ), **omes2patch} # could make more efficient by skipping redos

    hgx_dir = wrk_dir + 'hgx/'
    if not checkdir(hgx_dir, unzip = True, rm = True):
        os.mkdir(hgx_dir)

    print('\nVII. Quantifying OCG gene evolution congruence', flush = True)
    hgx2omes2gbc = gbc_main(
        hgx2loc, wrk_dir, ome2i, hg_dir, hgx_dir,
        blastp, db, gene2hg, plusminus, hg2gene, 
        old_path = 'hgx2omes2gbc.full.pickle',
        modules = ocgs, moduleHGxs = ocg_hgxs,
        moduleOmes = ocg_omes, cpus = cpus
        )

    if coevo_thresh > 0 or patch_thresh > 0 or microsyn_thresh > 0:
        print('\tApplying thresholds', flush = True)
        print('\t\t' + str(len(ocgs)) + ' ocgs before', flush = True)
        # edit ocgs, hgx2dist
 #       newOgx2dist, newOCGs, newOCGOmes, newOCGHGxs = {}, [], [], []
        newOCGs, newOCGOmes, newOCGHGxs = [], [], []
        for i0, ocg in enumerate(ocgs):
            check = False # have we added a new list
            ogc = ocg_hgxs[i0]
            omesc = ocg_omes[i0]
            for x, omes in ocg.items():
                if hgx2omes2gbc[ogc][omesc] >= coevo_thresh \
                    and omes2patch[omesc] >= patch_thresh:
                    if check:
                        newOCGs[-1][x] = omes
                    else:
                        newOCGs.append({x: omes})
                        check = True
            if check:
                newOCGOmes.append(ocg_omes[i0])
                newOCGHGxs.append(ogc)
        ocgs, ocg_omes, ocg_hgxs = \
            newOCGs, newOCGOmes, newOCGHGxs
        print('\t\t' + str(len(ocgs)) + ' ocgs after', flush = True)

    if run_dnds: # need to bring file naming to speed
        print('\nIIX. Quantifying OCG dn/ds', flush = True)
        omes4hgs, hgx_files = {x: defaultdict(list) for x in range(len(ome2i))}, []
        for i, ocg in enumerate(ocgs):
            for hgx in ocg:
#                hgx = hgx2i[x]
                for gene in hgx2loc[hgx]:
                    omeI = ome2i[gene[:gene.find('_')]]
                    if omeI in set(ocg[i][hgx]):
                        omes4hgs[omeI][gene].append(ocg_hgxs[i])

        hgx_files = []
        for hgx in ocg_hgxs:
            for og in hgx:
                out_base = hgx_dir + '-'.join([str(x) for x in hgx]) + '.' + str(og)
                hgx_files.append(out_base)

        maffts = [file + '.mafft' for file in hgx_files]
        if not all(os.path.isfile(x + '.aa.fa') for x in hgx_files):
            print("\tCompiling HGx kernels' genes", flush = True)
            dnds_preparation(db, omes4hgs, hgx2omes, gene2hg, plusminus, hgx_dir, i2ome)
        if not all(os.path.isfile(x) for x in maffts):
            print('\tAligning proteins', flush = True)
            maffts = [x for x in maffts if not os.path.isfile(x)]
            mafft_cmds = [
                [['mafft', '--auto', '--thread', '2', x.replace('.mafft', '.aa.fa')], x, 
                x.replace('.mafft', '.aa.fa')] \
                for x in maffts
                ]
            mafft_res = multisub(mafft_cmds, processes = cpus, stdout = True, rm = True)
            maffts = [x['output'] for x in mafft_res]

        print('\tCalculating dn/ds', flush = True)
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            dnds_res = pool.map(calc_dnds, maffts)
        pool.join()
    
        hgx2dnds = parse_dnds(dnds_res)
    else:
        print('\nIIX. Skipping quantifying selection', flush = True)
        hgx2dnds = {}

    print('\nIX. Writing and annotating clusters', flush = True)
    print('\tWriting cluster scores', flush = True)
    ocg_output, done, maxval = [], set(), max(hgx2dist.values())
    for ocg, modHGx2omes in enumerate(ocgs):
        for hgx, omes in modHGx2omes.items():
            hgx_id = hgx2i[hgx]
            if hgx_id not in done:
                done.add(hgx_id)
            else:
                continue
            ocg_output.append([
                ','.join([str(x) for x in hgx]), hgx_id,
                ocg, hgx2dist[hgx]/maxval, hgx2gbc[hgx], 
#                omes2patch[tuple(hgx2omes[hgx])],
                ','.join([i2ome[x] for x in hgx2omes[hgx]])
                ]) # HGxs at this stage are not segregated into groups
    ocg_output = sorted(ocg_output, key = lambda x: x[3], reverse = True)
    with gzip.open(out_dir + 'hgxs.tsv.gz', 'wt') as out:
        out.write('#hgs\thgx_id\tocg\tdistance\tcoevolution\tomes')#\tpatchiness\tomes')
        for entry in ocg_output:
            out.write('\n' + '\t'.join([str(x) for x in entry]))

    kern_output, top_hgxs = [], []
    for i, ocg in enumerate(ocg_hgxs):
        omesc = ocg_omes[i]
        try:
            if ocg not in dnds_dict:
                dnds_dict[ocg] = ['na' for x in range(3)]
        except NameError:
            dnds_dict = {ocg: ['na' for x in range(3)]}
        if ocg in hgx2omes2gbc:
            kern_output.append([
                ','.join([str(x) for x in ocg]), i, 
                omes2dist[omesc]/maxval, hgx2omes2gbc[ocg][omesc], omes2patch[omesc], 
                ','.join([str(i2ome[x]) for x in omesc])
#                dnds_dict[hgx][0], dnds_dict[hgx][1], str(dnds_dict[hgx][2]),
#                omes2dist[omesc]
                ])
        else:
            kern_output.append([
                ','.join([str(x) for x in ogc]), i,
                omes2dist[omesc]/maxval, hgx2omes2gbc[ogc][omesc], omes2patch[omesc], 
                ','.join([str(i2ome[x]) for x in omesc]) 
#                dnds_dict[hgx][0], dnds_dict[hgx][1], str(dnds_dict[hgx][2]),
               # omes2dist[omesc]
                ])
        for shgx, omes in ocgs[i].items():
            top_hgxs.append([shgx, i, set(omes)])

    kern_output = sorted(kern_output, key = lambda x: x[4], reverse = True)


    with gzip.open(out_dir + 'ocgs.tsv.gz', 'wt') as out:
        out.write(
            '#hgs\tocg\tdistance\t' + \
            'coevolution\tpatchiness\tomes' #+ \
            #'selection_coef\tmean_dnds\tog_dnds\t' + \
         #   'total_dist'
            )
        for entry in kern_output:
            out.write('\n' + '\t'.join([str(x) for x in entry]))

    if run_dnds:
        axes = [[],[],[],[],[]]
    else:
        axes = [[],[],[]]
    labels = []

    for entry in kern_output:
        labels.append(
            ['Clan:', str(entry[1]) + ' | Omes: ' + entry[-1] + ' | OGs: ' + str(entry[0])]
            )
        axes[0].append(entry[4])
        axes[1].append(entry[3])
        axes[2].append(entry[2])
        if run_dnds:
            axes[3].append(entry[6])
            axes[4].append(entry[7])
    print('\tOutputting scatter plots', flush = True)
    axes_labels = [
        'Distribution Patchiness', 'Gene coevolution', 'Microsynteny Distance' #,
#        'Selection Coefficient', 'Mean dN/dS'
        ]
    if run_dnds:
        axes_labels.extend(['Selection Coefficient', 'Mean dN/dS'])
    fig = mk_subplots(labels, axes, axes_labels, alpha = 0.6)
    fig.write_html(out_dir + 'scores.scatters.html')
    fig = mk_3d(labels, axes[:3], axes_labels[:3], alpha = 0.7)
    fig.write_html(out_dir + 'kernel.scatter.html')


    print('\tCompiling clusters from annotations', flush = True)
    sig_clus = extrct_sig_clus(
        hgx2dist, hgx2loc, top_hgxs, ome2i
        )
    write_clus_cmds = []
    ome_dir = out_dir + 'ome/'
    if not os.path.isdir(ome_dir):
        os.mkdir(ome_dir)
    for ome in sig_clus:
        if not os.path.isdir(ome_dir + ome):
            os.mkdir(ome_dir + ome)
        gff = db[ome]['gff3']
        out_file = ome_dir + ome + '/info.out'
        write_clus_cmds.append([sig_clus[ome], ome, out_file, gff, gene2hg, plusminus])
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        out_genes = pool.starmap(write_clusters, write_clus_cmds)
    pool.join()
    genes = {x[0]: x[1] for x in out_genes if x[1]}

#    hgx_dirTar.join()

    print('\tApplying Pfam annotations', flush = True)
    prot_paths = {}
    for ome in genes:
        prot_paths[ome] = db[ome]['faa']
    pfamRes, failedOmes = pfamMngr(
        genes, prot_paths, wrk_dir, pfam, evalue = 0.01, threshold = 0.5, cpus = cpus
        )


    print('\nX. Outputting clusters', flush = True)
    print('\tCluster biofiles', flush = True)
    grabClus_cmds = []
    for ome in genes:
        if ome in failedOmes:
            continue
        gff_path = db[ome]['gff3']
        pro_path = db[ome]['faa']
        try:
            grabClus_cmds.append([genes[ome], gff_path, pro_path, ome, ome_dir, gene2hg, pfamRes[ome]])
        except KeyError:
            grabClus_cmds.append([genes[ome], gff_path, pro_path, ome, ome_dir, gene2hg])


    with mp.get_context('fork').Pool(processes = cpus) as pool:
        clus_info = pool.starmap(grabClus, grabClus_cmds)
    pool.join()

#    print('\tCluster SVGs', flush = True)
 #   ome_svg_cmds = []
  #  svg_width = 10
   # for ome in clus_info:
    #    ome_svg_cmds.append([ome_dir + ome + '/'])
#    with mp.get_context('fork').Pool(processes = cpus) as pool:
 #       pool.starmap(runGFF2SVG, ome_svg_cmds)


if __name__ == '__main__':
    # need these here because spawn mp context forces reimport
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from cogent3.phylo import nj
    from cogent3 import PhyloNode, load_tree
    from io import StringIO
    from Bio import BiopythonWarning
    import networkx as nx
    from scipy.sparse import lil_matrix
    from scipy.stats import hypergeom
    import warnings
    warnings.simplefilter('ignore', BiopythonWarning)
    from Bio import SeqIO, AlignIO, codonalign
    from dna_features_viewer import GraphicFeature, GraphicRecord
    

    description = \
    """Pharmaceuticals are primarily derived from biologically-produced compounds, their
    derivatives, and synthetic compounds inspired by natural pharmacophores.
    Many natural product specialized metabolites are produced by gene clusters, or 
    regions of genomes that contain colocalized genes with concerted 
    biosynthetic function. 
    Genomes rearrange and recombine through macroevolutionary timescales, so the 
    persistence of colocalized gene families within clusters is indicative of 
    selection. 
    og2clus detects these regions by identifying unexpected conserved microsynteny 
    of homology groups throughout a sample of organisms. 
    Some unexpectedly conserved regions are not directly involved in specialized 
    metabolite biosynthesis - og2clus filters results by implementing a
    multivariate model derived from 1) total distance between organisms with a 
    particular homolog combination, 2) phylogenetic topological correlation between the
    constituent genes of the homolog combination, and 3) phylogenetic patchiness of
    the homolog combination."""
    parser = argparse.ArgumentParser(description = description)
    i_opt = parser.add_argument_group('Inputs')
    i_opt.add_argument('-d', '--database', required = True, default = masterDB(), 
        help = 'MycotoolsDB. DEFAULT: masterdb')
 #   parser.add_argument('-i', '--input', 
#        help = 'Precomputed whitespace delimitted file of homologous sequences')
    i_opt.add_argument('-of', '--orthofinder',
        help = 'Precomputed OrthoFinder output directory')
    i_opt.add_argument('-i', '--input',
        help = 'Precomputed cluster results file')
    i_opt.add_argument('-c', '--constraint',
        help = 'Constrain microsynteny topology to species tree')
    i_opt.add_argument('--pfam', help = 'Pfam-A.hmm')

    det_opt = parser.add_argument_group('Detection parameters')
    det_opt.add_argument('-l', '--linclust', action = 'store_true',
        help = 'Homology inference via linclust; DEFAULT: mmseqs cluster')
    det_opt.add_argument('-w', '--window', default = 5, type = int,
        help = 'Max genes +/- for each locus window. DEFAULT: 5 (11 gene window)')
    det_opt.add_argument('-r', '--root', 
        help = 'Ome code or 2 ome codes to root upon, e.g. psicya1, ' + \
        '"psicya1 psicub1"; DEFAULT: midpoint')

    thr_opt = parser.add_argument_group('Thresholding')
    thr_opt.add_argument('-sp', '--seed_percentile', type = int, default = 75,
        help = 'Percentile of HG pair distances. DEFAULT: 75')
    thr_opt.add_argument('-op', '--hgx_percentile', type = int,
        help = 'Percentile of HGx microsynteny distances. ' + \
        'Must be less than -cp. DEFAULT: -cp')
    thr_opt.add_argument('-cp', '--clus_percentile', type = int, default = 80, 
        help = 'Percentile of HGx microsynteny distances . DEFAULT: 80')
    thr_opt.add_argument('-pt', '--patch_threshold', default = 0, type = float,
        help = "Threshold [0 < value < 1] of high order OG combinations' patchiness scores.")
    thr_opt.add_argument('-ct', '--coevo_threshold', default = 0, type = float, 
        help = "Threshold [0 < value < 1] of high order OG combinations' coevolution scores.")
    thr_opt.add_argument('--null_sample', type = int, default = 10000,
        help = 'Samples for null distributions. DEFAULT: 10,000')

    run_opt = parser.add_argument_group('Runtime options')
    run_opt.add_argument('-s', '--dnds', action = 'store_true', help = 'Run dN/dS calculations')
    run_opt.add_argument('--n50', help = 'Minimum assembly N50')
    run_opt.add_argument('--stop', action = 'store_true', 
                        help = 'Print GBC commands and exit')
    run_opt.add_argument('-o', '--output_dir')
    run_opt.add_argument('-n', '--new', action = 'store_true', 
        help = 'Overwrite old run')
    run_opt.add_argument('--cpus', default = 1, type = int)
    args = parser.parse_args()

    if args.root:
        if '"' in args.root or "'" in args.root:
            args.root = args.root.replace('"','').replace("'",'')
        root = args.root.split()
        root_txt = ','.join(root)
    else:
        root = []
        root_txt = 'midpoint'

    if not args.hgx_percentile:
        args.hgx_percentile = args.clus_percentile # set the default

    execs = ['mafft', 'hmmsearch', 'blastp', 'mcxload', 'mcxdump', 'mcl',
             'diamond']
    if args.orthofinder:
        of_out = format_path(args.orthofinder)
        if os.path.isdir(of_out):
            orthogroups = of_out + '/Orthogroups/Orthogroups.txt'
            hg_dir = of_out + '/Orthogroup_Sequences/'
        else:
            orthogroups = of_out
            hg_dir = os.path.dirname(of_out) + '../Orthogroup_Sequences/'
        if not os.path.isfile(hg_dir + 'OG0000000.fa'):
            hg_dir = None
        method = 'orthofinder'
    elif args.input:
        orthogroups = format_path(args.input)
        hg_dir = None
        method = 'mmseqs easy-cluster'
    elif args.linclust:
        method = 'mmseqs easy-linclust'
        hg_dir = None
        orthogroups = None
        execs.append('mmseqs')
    else:
        method = 'mmseqs easy-cluster'
        hg_dir = None
        orthogroups = None
        execs.append('mmseqs')

    args_dict = {
        'Homology groups': orthogroups, 'Sequence clusters': method, 'MycotoolsDB': args.database, 
        'Root': root_txt, 'Pfam DB': args.pfam, 'Window': args.window*2+1,
        'Seed Percentile': args.seed_percentile, #'Precluster threshold': args.clus_threshold,
        'Cluster percentile': args.clus_percentile, 'HGx percentile': args.hgx_percentile,
        'Patchiness threshold': args.patch_threshold, 'Coevolution threshold': args.coevo_threshold,
        'Null samples': args.null_sample, 'Calculate dN/dS': args.dnds, 'Minimum N50': args.n50,
        'Processors': args.cpus, 'Output directory': args.output_dir,
        'Overwrite': bool(args.new)
        }

    pfam = format_path(args.pfam)
    if not os.path.isfile(pfam):
        print('\nERROR: invalid Pfam-A.hmm path', flush = True)
        sys.exit(4)
    

    findExecs(execs, exit = set(execs))
    if args.clus_percentile < args.hgx_percentile:
        print('\nERROR: hgx percentile is greater than cluster percentile',
            flush = True)
        sys.exit(3)

    start_time = intro('orthogroups to clusters - og2clus', args_dict, 'Zachary Konkel, Jason Slot')
    date = datetime.strftime(start_time, '%Y%m%d')

    if not args.output_dir:
        out_dir = os.getcwd() + '/og2clus_' + date + '/'
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
    else:
        out_dir = format_path(args.output_dir)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
            out_dir += '/'

#    if args.constraint:
 #       constraint_path = format_path(args.constraint)
  #  else:
   #     constraint_path = None
    if len(str(args.seed_percentile)) > 1:
        seed_perc = float('.' + str(args.seed_percentile))
    else:
        seed_perc = float('.0' + str(args.seed_percentile))
    if len(str(args.clus_percentile)) > 1:
        clus_perc = float('.' + str(args.clus_percentile))
    else:
        clus_perc = float('.0' + str(args.clus_percentile))   
    if len(str(args.hgx_percentile)) > 1:
        hgx_perc = float('.' + str(args.hgx_percentile))
    else:
        hgx_perc = float('.0' + str(args.hgx_percentile))   


    db = mtdb(format_path(args.database))
    main(
        db, orthogroups, out_dir, plusminus = args.window,
        cpus = args.cpus, hg_dir = hg_dir, 
        seed_perc = seed_perc, #clus_thresh = args.clus_threshold,
        clus_perc = clus_perc, blastp= 'blastp',#seed_thresh = args.seed_threshold,
        hgx_perc = hgx_perc, pfam = pfam, samples = args.null_sample,
        run_dnds = args.dnds, n50thresh = args.n50,
        root = root, coevo_thresh = args.coevo_threshold, 
        patch_thresh = args.patch_threshold, method = method,
        printexit = args.stop, flag = bool(not args.new)
        )

    outro(start_time)
