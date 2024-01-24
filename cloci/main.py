#! /usr/bin/env python3
"""
Co-occurrence Locus and Orthologous Cluster Identifier (CLOCI)
(C) Zachary Konkel, Jason Slot 2023
BSD 3-clause License
"""

#NEED to make similarity coefficient adjustable for domains and HLGs
#NEED to detect when a bypass failed
#NEED to uniformize percent outputs
#NEED filter output option
#NEED to implement genome-based null
#NEED call OrthoFinder option
#NEED use hg_dir from OrthoFinder in input_parsing/hgx2gcf, currently doesnt work
#NEED to make gcl_mngr split different families with the same GCF
#NEED locus-based GCF hypergeometric average
#NEED TO ADD PERCENTILE TO OUTPUT
#NEED to implement gff2svg
#NEED to delete excess when checkpoints are reached
    # implement tarring effectively and --compress
#NEED to rerun evo_conco if gcfs.pickle doesn't exist
   # update to remove gcf dir if logs are violated

import os
import re
import sys
import copy
import gzip
import shutil
import pickle
import hashlib
import argparse
import numpy as np
import multiprocessing as mp
from math import log
from tqdm import tqdm
from itertools import chain
from datetime import datetime
from collections import defaultdict
from mycotools.lib.kontools import \
    intro, outro, format_path, collect_files, \
    findExecs, eprint, tardir, write_json, \
    mkOutput
from mycotools.lib.biotools import \
    gff2list
from mycotools.lib.dbtools import mtdb, primaryDB
#from mycotools.gff2svg import main as gff2svg
from mycotools import db2microsyntree
from cloci.lib import treecalcs, evo_conco, \
     hgx2hlgs, input_parsing, hgp2hgx, generate_nulls, output_data

# if numba is available, let's import it
try:
    # this is deprecated
    from numba import njit, jit
    @njit
    def est_conds(ome_arr, cooccur_array):
        """Estimate the conditional probability a given HG combination would occur based
        on the hypergeometric distribution of each HGs' presence-absence"""
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
    def est_conds(ome_arr, cooccur_array):
        """Estimate the conditional probability a given HG combination would occur based
        on the hypergeometric distribution of each HGs' presence-absence"""
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

# deprecated
def est_combo_probs(hgx, omes, genes_in_hg_in_ome, genes_in_ome, window,
                    cooccur_arr):
    """Estimate the probability a given HGx would exist given the abundance of
    each HG and the hypergeometric probability of each HG co-occurrence"""
#    win_size = [win_size[x] for x in win_size]
#    p_coeff = calc_coeff(og0, og1, tot_genes, win_size)
    p_conds = est_conds(np.array(omes), cooccur_arr)
    p_coef = 1
    for og in hgx:
        p_coef *= est_hypergeo(omes, genes_in_hg_in_ome[og],
                               genes_in_ome, window)
    p = p_conds * p_coef
    return hgx, p

# deprecated
def est_hypergeo(
    omes, genes_in_hg_in_ome, genes_in_ome, window
    ):
    """Estimate the hypergeometric-derived probability of a particular
    HG being randomly sampled"""
    # NEED to modify to account for sliding window

    # should this be for all omes, or each og in each ome and then multiplied?
    pval = hypergeom.sf(1, genes_in_ome, genes_in_hg_in_ome, window)
    return pval

# deprecated
def combo_prob_mngr(
    hgx2omes, omes2hg2genes, omes2genes, window, cooccur_array, cpus = 1
    ):
    """Manage the estimation of the probability of randomly sampling each HG
    co-occurrence combination"""

    cmds = []
    # prepare commands for each hgx
    for hgx, omes in hgx2omes.items():
        ome = omes[-1] # coefficient ome
        genes_in_hg_in_ome = {og: len(omes2hg2genes[ome][og]) for og in hgx}
        genes_in_ome = len(omes2genes[ome])
        cmds.append([
            hgx, omes, genes_in_hg_in_ome, genes_in_ome, window,
            cooccur_array
            ])

    # run hg-by-hg, then accumulate via hgx at the end
    with mp.get_context('fork').Pool(processes = cpus) as pool: # will fail on Windows
        hypergeoRes = pool.starmap(est_combo_probs, cmds)
        pool.close()
        pool.join()

    # prepare a Bonferroni correction
    comparisons = len(hgx2omes)
    hgx2pval = {}
    for hgx, pval in hypergeoRes:
        hgx2pval[hgx] = comparisons * pval

    return hgx2pval


def threshold_hgx(hgx2omes, hgx2loc, omes2dist, 
                  ome2partition, bord_scores_list):
    """Acquire HGxs that meet minimum microsynteny tree branch length
    distributions"""
    i2hgx, hgx2i, hgx2dist, count = {}, {}, {}, 0
    for hgx, omes in hgx2omes.items():
        dist = omes2dist[omes]
        parts = set(ome2partition[x] for x in omes)
        if None in parts:
            parts = parts.remove(None)
            if not parts:
                continue
        # choose the minimum threshold [liberal]
        bord_score = min([bord_scores_list[i][len(hgx)] for i in list(parts)])
        if dist > bord_score:
             hgx2dist[hgx] = dist
             i2hgx[count], hgx2i[hgx] = hgx, count
             count += 1
             # apply the border threshold for ALL pre-grouped hgxs

    print(
        '\t\t' + str(len(i2hgx)) + ' HGx pass border threshold', 
        flush = True
        )
    hgx2omes = {v: hgx2omes[v] for v in i2hgx.values()}
    hgx2loc = {v: hgx2loc[v] for v in i2hgx.values()}
    return hgx2omes, hgx2loc, hgx2dist, i2hgx, hgx2i


def read_tune_file(tune_file, gene2hg, ome2i):
    """Read a file to tune GCF circumscription based on the genomes that should
    be included in a particular GCF"""
    with open(tune_file, 'r') as raw:
        tune_data = [x.rstrip() for x in raw if not x.startswith('#')]
    tune = {}
    for cluster in tune_data:
        try:
            name, rawgenes, rawomes, false_omes = cluster.rstrip().split('\t')
        except ValueError:
            name, rawgenes, rawomes = cluster.split('\t')
            false_omes = []
        if name in tune:
            eprint(f'\nERROR: duplicate entry for {name} in tune file',
                   flush = True)
        if ',' in rawgenes:
            genes = [x.rstrip().lstrip() for x in rawgenes.split(',')]
        else:
            genes = [x.rstrip().lstrip() for x in rawgenes.split()]
        if ',' in rawomes:
            omes = [x.rstrip().lstrip() for x in rawomes.split(',')]
        else:
            omes = [x.rstrip().lstrip() for x in rawomes.split()]
        if false_omes:
            if ',' in false_omes:
                false_omes = [x.rstrip().lstrip() for x in false_omes.split(',')]
            else:
                false_omes = [x.rstrip().lstrip() for x in false_omes.split()]
        try:
            hgs = [int(i) for i in genes]
        except ValueError:
            hgs = []
            for gene in genes:
                if gene in gene2hg:
                    hgs.append(gene)
        if not hgs:
            eprint(f'\t\tWARNING: no homology groups in {name}, skipping',
                   flush = True)
            continue
        omes = [ome2i[x] for x in omes]
        false_omes = [ome2i[x] for x in false_omes]
        if not omes:
            eprint(f'\t\tWARNING: no omes in {name}, skipping', flush = True)
        tune[name] = [tuple(sorted(set(hgs))), 
                      tuple(sorted(set(omes))),
                      tuple(sorted(set(false_omes)))]
    return tune


def main(
    db, hg_file, out_dir, plusminus = 1, hg_dir = None, hgx_dir = None,
    hgp_perc = 0.2, clus_perc = 0, hgx_perc = 0.7,
    id_perc = 30, inflation_1 = 1.1, inflation_2 = 1.3, pos_perc = 0,
    csb_thresh = 0, minimum_omes = 2, samples = 10000, pfam = None,
    min_hlg_id = 0.3, simfun = hgx2hlgs.overlap,
    constraint_path = None, aligner = 'diamond',
    run_dnds = False, cpus = 1, n50thresh = None, 
    root = None, gcl_thresh = 0, patch_thresh = 0,
    method = 'mmseqs easy-cluster', dist_thresh = 0,
    printexit = False, skipalgn = False, flag = True, partition_file = None,
    near_single_copy_genes = [], tree_path = None, 
    verbose = False, sim = 'overlap', tune_file = None,
    dist_type = 'tmd', uniq_sp = False, partition_rank = None,
    min_branch_sim= 0, algn_sens = '', min_gene_id = 30,
    fallback = False, merge_via_sim = False, ipr_path = None, force = False
    ):
    """
    The general workflow:
    log management -> input data parsing -> homology group pair (HGp) 
    identification -> microsynteny distance thresholding -> homology 
    group combination (HGx) HGx formation -> microsynteny distance 
    border thresholding -> HGx-GCF grouping ->
    microsynteny distance cluster thresholding -> PDS calculation ->
    GCL calculation -> optional dN/dS calculations (deprecated) ->
    HGx data output -> cluster retrieving -> data output
    """
    db = db.set_index()
    if uniq_sp:
        uniq_sp = db
    if not tree_path:
        tree_path = f'{out_dir}microsynt.newick'

    # determine the distance function: TMD = total microsynteny distance, MMD =
    # maximum microsynteny distance
    if dist_type == 'tmd':
        dist_func = treecalcs.calc_tmd
    else:
        dist_func = treecalcs.calc_mmd

    # import the null distribution partition data, use an established rank, or
    # generate a single null distribution for all genomes
    if partition_file:
        partition = partition_file
    elif partition_rank:
        partition = partition_rank
    else:
        partition = None

    # initialize the run by parsing and updating the log
    wrk_dir, nul_dir, inflation_1, inflation_2, gcf_ready = \
                       input_parsing.init_run(db, out_dir, 
                                              near_single_copy_genes, constraint_path,
                                              tree_path, hg_file, plusminus, 
                                              hgp_perc, #clus_perc,
                                              hgx_perc, f'{aligner}${algn_sens}', 
                                              id_perc, pos_perc, csb_thresh,
                                              patch_thresh, gcl_thresh, dist_thresh,
                                              samples, n50thresh, flag, min_gene_id,
                                              min_hlg_id, inflation_1, inflation_2, sim, 
                                              tune_file, dist_type, uniq_sp, partition,
                                              min_branch_sim, merge_via_sim, hg_dir, hgx_dir,
                                              ipr_path, pfam)

    # we are going to rerun, regardless what the log parser claims
    if force and gcf_ready:
        print('\tJust kidding', flush = True)
        gcf_ready = False

    # !! sometimes will not enter on benign changes to log
    # Skip straight to output thresholding to avoid rerunning and regenerating
    # unnecessary data structures
    if gcf_ready:
        print('\nThresholding and outputting GCFs', flush = True)
        try:
            ome_dir = out_dir + 'ome/'
            annotate = False
            if ipr_path or pfam and annotate: # check if all are annotated
                hlg_files_p = collect_files(ome_dir, 'tsv', recursive = True)
                hlg_files = [x for x in hlg_files_p \
                            if os.path.basename(x) == 'hlg.tsv']
                for hlg_f in hlg_files:
                    with open(hlg_f, 'r') as raw:
                        for line in raw:
                            if len(line.split()) < 5: # annotations missing
                                annotate = True
                            break
                    # a long filter is necessary to annotate missing
                    # annotations
                    if annotate:
                        print('\t\tAnnotations missing and requested. Long filter',
                              flush = True)
                        break
            # proceed to the long process
            if annotate:         
                ome2i, gene2hg, i2ome, hg2gene, ome2pairs, cooccur_dict = \
                    db2microsyntree.main(db, hg_file, out_dir, wrk_dir,
                                    method, tree_path, plusminus = plusminus,
                                    min_cov = 0, min_id = 0.3, n50thresh = n50thresh,
                                    near_single_copy_genes = near_single_copy_genes,
                                    constraint = constraint_path, verbose = verbose,
                                    return_post_compile = gcf_ready, cpus = cpus)
                   
                output_data.threshold_gcf_bypass(db, out_dir, wrk_dir, i2ome, gene2hg, 
                                                 dist_thresh, gcl_thresh, patch_thresh,
                                                 id_perc, pos_perc, csb_thresh, ipr_path,
                                                 pfam, cpus)
            # proceed to the short thresholding process
            else:
                output_data.threshold_gcf_quick(db, out_dir, ome_dir, dist_thresh, gcl_thresh,
                                    patch_thresh, id_perc, pos_perc, csb_thresh, cpus)
            print('\nSUCCESS!', flush = True)
            sys.exit(0)

        # there are missing files that are necessary for one of the above
        # thresholding processes, so we have to proceed through the standard
        # script
        except FileNotFoundError:
            eprint('\t\tWARNING: necessary files missing; must resume run', flush = True)
            raise FileNotFoundError
            ome2i, gene2hg, i2ome, hg2gene, ome2pairs, cooccur_dict = \
                db2microsyntree.main(db, hg_file, out_dir, wrk_dir,
                                    method, tree_path, plusminus = plusminus,
                                    min_cov = 0, min_id = 0.3, n50thresh = n50thresh,
                                    near_single_copy_genes = near_single_copy_genes,
                                    constraint = constraint_path, verbose = verbose,
                                    return_post_compile = gcf_ready, cpus = cpus)

    # build the microsynteny tree, and toss-out any genomes that
    # failed to have any overlapping HGps with others
    else:
        ome2i, gene2hg, i2ome, hg2gene, ome2pairs, cooccur_dict = \
            db2microsyntree.main(db, hg_file, out_dir, wrk_dir,
                                method, tree_path, plusminus = plusminus,
                                min_cov = 0, min_id = 0.3, n50thresh = n50thresh,
                                near_single_copy_genes = near_single_copy_genes,
                                constraint = constraint_path, verbose = verbose,
                                return_post_compile = gcf_ready, cpus = cpus)
    
    # read a tune file to tune GCF formation downstream
    if tune_file:
        print('\tReading tune file', flush = True)
        tune = read_tune_file(tune_file, gene2hg, ome2i)
    else:
        tune = None

    # compile the microsynteny tree
    print('\tReading microsynteny tree', flush = True)
    phylo = input_parsing.compile_tree(
        i2ome, tree_path, root = root
        )
    # generate null distributions for HGps
    partition_omes, ome2partition, omes2dist, min_pair_scores = \
                                generate_nulls.gen_pair_nulls(
                                            db, phylo, ome2i, wrk_dir,
                                            nul_dir, hgp_perc, ome2pairs,
                                            i2ome, samples = samples,
                                            partition_file = partition_file,
                                            partition_rank = partition_rank,
                                            uniq_sp = uniq_sp, dist_func = dist_func,
                                            cpus = cpus
                                            )
    # seed clusters by calcuating total microsynteny distance for HGps
    print('\nIII. Seeding HG pairs (HGp)', flush = True) 
    if not os.path.isfile(out_dir + 'hgps.tsv.gz'):
        print('\tCalculating seed HG-pair scores', flush = True)
        seed_score_start = datetime.now()
        # we are looking at all genomes
        if not uniq_sp:
            omes2dist = treecalcs.update_dists(phylo, cooccur_dict, cpus, omes2dist,
                                               func = dist_func)
        # we are only looking at genomes of unique species
        else:
            omes2dist = treecalcs.update_dists(phylo, cooccur_dict, cpus, omes2dist,
                                               func = treecalcs.calc_tmd,
                                               uniq_sp = db, i2ome = i2ome)
        # output the updated dictionary that relates a sorted tuple of ome
        # indices to their microsynteny distance (TMD/MMD) between them
        with open(wrk_dir + 'omes2dist.pickle', 'wb') as out:
            pickle.dump(omes2dist, out)

        # identify the HGps that pass thresholding
        top_hgs = []
        for hgp, omes in cooccur_dict.items():
            score = omes2dist[omes]
            parts = set(ome2partition[x] for x in omes)
            if None in parts:
                parts = parts.remove(None)
                if not parts:
                    continue
            # choose the minimum threshold [liberal]
            min_pair_score = min([min_pair_scores[i] for i in list(parts)])
            if score >= min_pair_score:
                top_hgs.append([
                    hgp[0], hgp[1], len(omes), score#, score/ome_dist
                    ])

        # write scores            
        print('\t\t' + str(datetime.now() - seed_score_start), flush = True)
        top_hgs = sorted(top_hgs, key = lambda x: x[3])
        print('\tWriting seed scores', flush = True)
        with gzip.open(out_dir + 'hgps.tsv.gz', 'wt') as out:
            out.write('#hg0\thg1\tcooccurrences\tscore\tadj_score\n')
            for line in top_hgs:
                out.write('\t'.join([str(x) for x in line]) + '\n')
        print('\t\t' + str(len(top_hgs)) + ' significant HG pairs', flush = True)
        top_hgs = [(x[0], x[1]) for x in top_hgs]
    # load previous HGps
    elif not os.path.isfile(wrk_dir + 'hgx2loc.pickle'):
        print('\tLoading previous seed HG pairs', flush = True)
        top_hgs = input_parsing.load_seedScores(out_dir + 'hgps.tsv.gz')
        print('\t\t' + str(len(top_hgs)) + ' significant HG pairs', flush = True)
    else:
        top_hgs = None

    # begin sifting for HGxs using pairs as seeds for HGx detection
    print('\nIV. Sprouting high order HG combinations (HGx)', flush = True)
    hgx2omes, hgx2loc = hgp2hgx.hgp2hgx(db, wrk_dir, top_hgs,
                                        gene2hg, ome2i, phylo, 
                                        plusminus, cpus) 

    print('\tCalculating HGx microsynteny distances', flush = True)
    hgx_start = datetime.now()
    hgx_obs = len(hgx2omes)
    print('\t\t' + str(hgx_obs) + ' observed HGx', flush = True)

    # calculate HGx microsynteny distances
    if not uniq_sp:
        omes2dist = treecalcs.update_dists(phylo, hgx2omes, cpus, omes2dist = omes2dist,
                                           func = dist_func)
    else:
        omes2dist = treecalcs.update_dists(phylo, hgx2omes, cpus, omes2dist = omes2dist,
                                           func = treecalcs.calc_tmd, uniq_sp = db,
                                           i2ome = i2ome)
    with open(wrk_dir + 'omes2dist.pickle', 'wb') as out:
        pickle.dump(omes2dist, out)
    print('\t\t' + str(datetime.now() - hgx_start), flush = True)

    # deprecated mechanism of estimating the probability of randomly sampling
    # an HGx from the genomes provided - needs conceptual work
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
    bord_scores_list, clus_scores_list = generate_nulls.partitions2hgx_nulls(
                                            db, partition_omes, ome2i, i2ome, gene2hg,
                                            max_hgx_size, plusminus, hgx_perc, 
                                            clus_perc, nul_dir, omes2dist, 
                                            phylo, samples, dist_func = dist_func,
                                            uniq_sp = uniq_sp, cpus = cpus)


    # collect for normalizing relative to the whole dataset later
    # should I move into absolute space? it may just be better for comparing
    # datasets in the future and some algorithm is going to pick up that gene clusters
    # are in some absolute microsynteny distance, so long as its not normalized. 

    hgx2omes, hgx2loc, hgx2dist, i2hgx, hgx2i = threshold_hgx(hgx2omes, hgx2loc,
                                                    omes2dist, ome2partition, 
                                                    bord_scores_list)
    with open(wrk_dir + 'hgx2dist.pickle', 'wb') as out:
        pickle.dump(hgx2dist, out)
    output_data.output_hgxs(hgx2dist, hgx2omes, hgx2i, i2ome, out_dir)

    hgx_dir = f'{wrk_dir}hgx/'

    # placeholder for legacy runs
    if not os.path.isfile(f'{wrk_dir}hgx2omes.pickle'):
        with open(wrk_dir + 'hgx2omes.pickle', 'wb') as out:
            pickle.dump(hgx2omes, out)

    print('\nV. Circumscribing homologous locus groups (HLGs) from HGxs', flush = True)
    if not os.path.isfile(wrk_dir + 'hlgs.pickle'): # need to add this to
    # log parsing
        # Group hgxs
        hlgs, hlg_hgxs, hlg_omes, hlg2clan, omes2dist = hgx2hlgs.classify_hlgs(
            hgx2loc, db, gene2hg, i2hgx, hgx2i, phylo, 
            ome2i, hgx2omes, hg_dir, hgx_dir, wrk_dir, ome2partition, 
            bord_scores_list, hg2gene, tune = tune, algorithm = aligner,
            omes2dist = omes2dist, clusplusminus = plusminus, 
            inflation_1 = inflation_1, inflation_2 = inflation_2,
            min_loc_id = min_hlg_id, algn_sens = algn_sens,
            min_omes = 2, cpus = cpus, simfun = simfun, printexit = printexit,
            dist_func = dist_func, uniq_sp = uniq_sp, min_branch_sim = min_branch_sim,
            skipalgn = skipalgn, minid = min_gene_id, fallback = fallback,
            merge_via_sim = merge_via_sim
            )
        with open(wrk_dir + 'hlgs.pickle', 'wb') as pickout:
            pickle.dump([hlgs, hlg_omes, hlg_hgxs, hlg2clan], pickout)
        with open(wrk_dir + 'omes2dist.pickle', 'wb') as pickout:
            pickle.dump(omes2dist, pickout)
    else: # or just load old data
        with open(wrk_dir + 'hlgs.pickle', 'rb') as raw:
            hlgs, hlg_omes, hlg_hgxs, hlg2clan = pickle.load(raw)
        # placeholder for legacy data
#        if not isinstance(hlgs, dict):
 #           hlgs = {i: v for i, v in enumerate(hlgs)}
  #          hlg_omes = {i: v for i, v in enumerate(hlg_omes)}
   #         hlg_hgxs = {i: v for i, v in enumerate(hlg_hgxs)}
    #        hlg2clan = {i: v for i, v in enumerate(hlg2clan)}

    if not uniq_sp:
        omes2dist = treecalcs.update_dists(phylo, {i: omes for i, omes in hlg_omes.items()}, 
                                        cpus, omes2dist, func = dist_func)
    else:
        omes2dist = treecalcs.update_dists(phylo, {i: omes for i, omes in hlg_omes.items()},
                                        cpus, omes2dist, func = treecalcs.calc_tmd,
                                        uniq_sp = db, i2ome = i2ome)
    with open(wrk_dir + 'omes2dist.pickle', 'wb') as out:
        pickle.dump(omes2dist, out)


#    if dist_thresh:
 #       print(f'\t{len(hlgs)} HLGs pass threshold', flush = True)
    # this is a crude pseudo cluster family expectedness check; its pseudo because the HGx
    # is now the conglomerate of the extracted loci's, and if that locus is derived from
    # multiple HGxs or the set of omes with that locus is truncated relative to the original
    # HGx then it is not obtained in the same way the null distribution was constructed;
    # nevertheless, this is a liberal mechanism that will remove edge cases to reduce computation
    # If the user doesn't want this, then remove the family percentile argument.
 #   todel = []
#    for i, omes in hlg_omes.items():
  #      dist = omes2dist[omes]
   #     parts = set(ome2partition[x] for x in omes)
    #    if None in parts:
     #       parts = parts.remove(None)
      #      if not parts:
       #         continue
        # get the threshold corresponding to the maximum sized HGx 
        # in the family, and the minimum lineages' value [liberal]
#        max_len = max([len(x) for x in hlgs[i]])
   #     max_len = len(hlg_hgxs[i]) # get the len of the complete HGx
        # if its too big, we cant evaluate its pseudo-expectedness, so pass
    #    if max_len <= plusminus * 2 + 1:
     #       clus_score = min([clus_scores_list[v][max_len] \
      #                        for v in list(parts)])
       #     if not dist > clus_score:
        #        todel.append(i)

#    todel.reverse()
 #   for i in todel:
  #      del hlg_omes[i]
   #     del hlgs[i]
    #    del hlg_hgxs[i]
     #   del hlg2clan[i]

    print('\nVI. Quantifying GCL and intra-HLG alignment similarity', flush = True)
    hgx2omes2gcl, hgx2omes2id, hgx2omes2pos, hlgs, hlg_omes = evo_conco.gcl_main( 
                            hgx2loc, wrk_dir, ome2i, hg_dir, hgx_dir,
                            aligner, db, gene2hg, plusminus, hg2gene, phylo,
                            old_path = 'mmi.pickle',
                            hlgs = hlgs, hlg_hgxs = hlg_hgxs, hlg2clan = hlg2clan,
                            minid = min_gene_id, fallback = fallback,
                            hlg_omes = hlg_omes, printexit = printexit,
                            skipalgn = skipalgn, cpus = cpus
                            )

    runOmes = [
        omes for omes in hlg_omes.values()
        ] # omes without pds scores
    runHGxs = [
        hgx for i, hgx in hlg_hgxs.items()
        ] # hgxs without pds scores

    print('\nVII. Quantifying PDS', flush = True)
    omes2patch = treecalcs.patch_main(phylo, runOmes, wrk_dir,
                            old_path = 'pds.pickle', 
                            cpus = cpus) # could make more efficient by skipping redos

    # this is the output directory for hgx2hlg and evo_conco
    hg_dir = f'{wrk_dir}hg/'

    # deprecated placeholder for dnds methodology
    hgx2dnds = {}
    if not uniq_sp:
        omes2dist = treecalcs.update_dists(phylo, 
                                        {i: v for i, v in hlg_omes.items()}, 
                                        cpus, omes2dist, func = dist_func) 
    else:
        omes2dist = treecalcs.update_dists(phylo,
                                        {i: v for i, v in hlg_omes.items()},
                                        cpus, omes2dist, func = treecalcs.calc_tmd,
                                        uniq_sp = db, i2ome = i2ome)
    with open(wrk_dir + 'omes2dist.pickle', 'wb') as out:
        pickle.dump(omes2dist, out)


    print('\nIIX. Writing and annotating clusters', flush = True)
    output_data.output_hlgs(db, wrk_dir, hlgs, hlg_omes, 
         i2ome, out_dir, hlg_hgxs,
         omes2dist, omes2patch, hgx2omes2gcl, hgx2omes2id,
         hgx2omes2pos, gene2hg, plusminus, ome2i,
         hlg2clan, dist_thresh, gcl_thresh, patch_thresh, id_perc,
         pos_perc, csb_thresh, ipr_path = ipr_path, pfam_path = pfam, dnds_dict = {}, 
         cpus = cpus)


def cli():
    # need these here because spawn mp context forces reimport
#    from scipy.stats import hypergeom
 #   from dna_features_viewer import GraphicFeature, GraphicRecord
    null_ranks = ['kingdom', 'phylum', 'subphylum', 'class', 
                  'order', 'family', 'genus', 'species']
    aligners = ['mmseqs', 'diamond', 'blastp']
    description = \
           """CLOCI elucidates gene clusters de novo by first
           classifying microsynteny-informed sub-cluster domains and
           homologous locus groups. CLOCI calls gene cluster families 
           by filtering homologous locus groups using measurements
           of gene cluster selection. CLOCI requires evidence of 
           unexpected microsynteny to call loci, and thus requires an
           adequate sample. The recommended sample will vary depending
           on your lineage's rate of microsynteny decay, horizontal transfer,
           and other factors. A subphylum-level analysis is a generally good
           for fungi. If you are filtering homologous locus groups for
           a particular type of syntenic locus, such as gene clusters, then
           run your dataset, and tune your minimum values to a reference set
           of gene clusters. These values can vary depending on the sample,
           so it is ideal to retune filtering parameters for each sample."""
    parser = argparse.ArgumentParser(description = description)
    i_opt = parser.add_argument_group('Input parameters')
    i_opt.add_argument('-d', '--database', required = True, default = primaryDB(), 
        help = 'MycotoolsDB. DEFAULT: masterdb')
 #   parser.add_argument('-i', '--input', 
#        help = 'Precomputed whitespace delimitted file of homologous sequences')
    i_opt.add_argument('--pfam', help = 'Pfam-A.hmm for Pfam annotations')
    i_opt.add_argument('--interpro', help = 'interproscan.sh for IPR annotations')

    mt_opt = parser.add_argument_group('Microsynteny tree parameters')
    mt_opt.add_argument('-f', '--focal_genes',
        help = 'File of genes for neighborhood extraction of microsynteny ' \
             + 'tree')
    mt_opt.add_argument('-c', '--constraint',
        help = 'Constrain microsynteny topology to species tree w/ome code tips')
    mt_opt.add_argument('-r', '--root', 
        help = 'Ome(s) to root upon; DEFAULT: midpoint')
    mt_opt.add_argument('-t', '--tree',
        help = 'Precomputed microsynteny tree path')

    det_opt = parser.add_argument_group('Detection parameters')
    det_opt.add_argument('-of', '--orthofinder',
        help = 'Precomputed OrthoFinder output directory')
    det_opt.add_argument('-i', '--input',
        help = 'Precomputed cluster results file')
    det_opt.add_argument('-l', '--linclust', action = 'store_true',
        help = 'Less sensitive Homology inference via linclust; DEFAULT: mmseqs cluster')
    det_opt.add_argument('-w', '--window', default = 2, type = int,
        help = 'Max genes +/- for each locus window. DEFAULT: 2 (5 gene window)')
    det_opt.add_argument('-mmd', '--maximum_dist', action = 'store_true',
        help = 'Calculated maximum microsynteny distance; DEFAULT: total microsynteny distance')
    det_opt.add_argument('-us', '--unique_sp', action = 'store_true',
        help = 'Only consider one genome for replicate species in microsynteny distance' \
             + ' calculations', default = False)

    nul_opt = parser.add_argument_group('Null parameters')
    nul_opt.add_argument('-nr', '--null_rank', default = 'species',
        help = f'Taxonomic rank for local null models {null_ranks}; DEFAULT: species')
    nul_opt.add_argument('-np', '--null_partitions', 
        help = 'Tab-delimited file with omes for each null sample group on ' \
             + 'separate lines "#<LINEAGE>\\n<OME1> <OME2> <OMEn>"')
    nul_opt.add_argument('-ns', '--null_sample', type = int, default = 10000,
        help = 'Samples for null distributions; DEFAULT: 10,000')


    fam_opt = parser.add_argument_group('Locus aggregation parameters')
    fam_opt.add_argument('-a', '--aligner', default = 'diamond',
        help = f'Alignment algorithm: {aligners}; DEFAULT: diamond')
    fam_opt.add_argument('-sa', '--sensitive_align',
        help = f'[diamond ultra-sensitive] or [mmseqs -s 7.5 --num_iterations 3]')
    fam_opt.add_argument('-s', '--similarity', default = 'sorensen',
        help = 'HLG similarity metric: [J]accard, [O]verlap, [S]orensen; ' \
             + 'DEFAULT: Sorensen')
    fam_opt.add_argument('-mg', '--minimum_gene_id', type = float, default = 45,
        help = 'Percent [30 < value < 100] ID minimum between gene for loci; ' \
             + 'DEFAULT: 45')
    fam_opt.add_argument('-ml', '--minimum_loc_id', type = float, default = 30,
        help = 'Percent [0 < value < 100] ID minimum between loci for HLG ' \
             + 'DEFAULT: 30')
#    fam_opt.add_argument('-tm', '--topology_merge', action = 'store_true', default = False,
 #       help = 'Merge loci with similar topologies prior to 1st MCL round via -ts')
    fam_opt.add_argument('-ts', '--min_topology_sim', type = float, default = 25,
        help = 'Percent [0 < value < 100] topology similarity (Jaccard) minimum ' \
             + 'for singleton merging AND merging loci prior to HLG aggregation; ' \
             + 'DEFAULT: 25')
#    fam_opt.add_argument('-I', '--inflation', default = 1.3, type = float,
 #       help = 'MCL inflation during family detection; DEFAULT: 1.3')
    fam_opt.add_argument('-I1', '--inflation_rnd1', default = 1.1, type = float,
        help = 'MCL inflation 1: affects domain/merging granularity; DEFAULT: 1.1')
    fam_opt.add_argument('-I2', '--inflation_rnd2', default = 1.3, type = float,
        help = 'MCL inflation 2: affects HLG/GCF granularity; DEFAULT: 1.3')
    fam_opt.add_argument('-T', '--tune', 
        help = 'Tune inflation to subset data of clusters represented ' \
             + 'in a tab-delimited ' \
             + 'file "<CLUSTER>\\t<CONSERVED_HGS/GENES>\\t<OMES>\\t<EXCLUDED_OMES>"')
   

    thr_opt = parser.add_argument_group('GCF filtering parameters')
    thr_opt.add_argument('-hp', '--hgp_percentile', type = int, default = 20,
        help = 'Null percentile of HG pair distances; DEFAULT: 20')
    thr_opt.add_argument('-xp', '--hgx_percentile', type = int,
        help = 'Null percentile [0 < value < 100] of HGx microsynteny distances. ' \
             + 'Must be less than -fp; DEFAULT: 60', default = 60)
#    thr_opt.add_argument('-fp', '--gcf_percentile', type = int, default = 0, 
 #       help = 'Pseudonull percentile [0 < value < 100] of GCF microsynteny distances')
    thr_opt.add_argument('-ip', '--id_percent', default = 0, type = float,
        help = 'Percent [0 < value < 100] identity minimum for gene cluster family')
    thr_opt.add_argument('-pp', '--pos_percent', default = 0, type = float,
        help = 'Percent [0 < value < 100] positive minimum for gene cluster family')
    thr_opt.add_argument('-ct', '--csb_threshold', default = 0, type = float,
        help = "Threshold [0 < value < 1] conservative substitution bias minimum for " \
             + ' gene cluster family')
    thr_opt.add_argument('-pt', '--pds_threshold', default = 0, type = float,
        help = "Threshold [0 < value < 1] of gene cluster family " \
             + " phylogenetic distribution sparsity scores")
    thr_opt.add_argument('-gt', '--gcl_threshold', default = 0, type = float,
        help = "Threshold [0 < value < 1] of gene cluster committment scores")
    thr_opt.add_argument('-tt', '--md_threshold', default = 0, type = float,
        help = "Threshold [0 < value < 1] of log normalized GCF TMDs")

    run_opt = parser.add_argument_group('Runtime parameters')
#    run_opt.add_argument('-s', '--dnds', action = 'store_true', help = 'Run dN/dS calculations')
    run_opt.add_argument('--n50', help = 'Minimum assembly N50')
    run_opt.add_argument('--stop', action = 'store_true', 
                        help = 'Export HG alignment commands ' \
                             + 'for parallelization and stop')
    run_opt.add_argument('--skip', action = 'store_true', 
                        help = 'Ignore missing HG alignments as assumed failures')
    run_opt.add_argument('--fallback', action = 'store_true',
                        help = 'Fallback to diamond from failed alignments')
    run_opt.add_argument('-n', '--new', action = 'store_true', 
        help = 'Rerun with new parameters and overwrite incompatible data')
    run_opt.add_argument('--force', action = 'store_true',
        help = 'Force rerun over bypassing')
    run_opt.add_argument('--compress', action = 'store_true', 
        help = 'Compress run; SEMI-FUNCTIONAL')
    run_opt.add_argument('--cpus', default = mp.cpu_count(), type = int,
                         help = 'DEFAULT: all')

    dir_opt = parser.add_argument_group('Alternative directories')
    dir_opt.add_argument('-o', '--output_dir', 
        help = 'Output/resume directory; DEFAULT: cloci_YYYYmmdd')
    dir_opt.add_argument('-hg', '--hg_dir', help = 'HG faa dir, format <HG>.faa')
    dir_opt.add_argument('-hgx', '--hgx_dir', 
        help = 'HGx alignment DB and results dir, format <HG>.out and <HG>.dmnd/<HG>.mmseqs')
    args = parser.parse_args()

    # NEED to reinstate, currently does not function correctly
    if args.compress:
        if not args.output_dir:
            eprint('\nERROR: compression requires -o', flush = True)
        elif not os.path.isdir(args.output_dir):
            eprint('\nERROR: -o directory does not exist', flush = True)
        dirs = [format_path(args.output_dir) + 'working/hgx/']
        for d in dirs:
            tardir(d, True)
        sys.exit(4)

    # grab the proposed root for the microsynteny tree
    if args.root:
        if '"' in args.root or "'" in args.root:
            args.root = args.root.replace('"','').replace("'",'')
        root = args.root.split()
        root_txt = ','.join(root)
    else:
        root = []
        root_txt = 'midpoint'

 #   if not args.hgx_percentile:
#        args.hgx_percentile = args.gcf_percentile # set the default

    min_topology_sim = args.min_topology_sim/100

    # parse the similarity coefficient used in HLG circumscription
    if not args.similarity:
        simfun = hgx2hlgs.overlap # set default function
        sim = 'overlap'
    elif args.similarity.lower() in {'o', 'overlap', 'oc'}:
        simfun = hgx2hlgs.overlap
        sim = 'overlap'
    elif args.similarity.lower() in {'j', 'jaccard', 'jacard', 'jac'}:
        simfun = hgx2hlgs.jaccard
        sim = 'Jaccard'
    elif args.similarity.lower() in {'s', 'sorensen', 'sorenson', 'sor'}:
        simfun = hgx2hlgs.sorensen
        sim = 'Sorensen'
    else:
        eprint(f'\nERROR: invalid -s: {args.similarity}',
               flush = True)
        sys.exit(43)

    # check for dependency presence
    execs = ['mafft', 'mcl',
             'mcxdump', 'mcxload', 'iqtree']
    if args.aligner.lower() not in aligners:
        eprint(f'\nERROR: invalid -a: {args.aligner}',
               flush = True)
        sys.exit(47)
    else:
        aligner = args.aligner.lower()
        if args.sensitive_align:
            algn_sens = 'sensitive'
        else:
            algn_sens = ''
        execs.append(aligner)
        if aligner == 'mmseqs' and args.pos_percent:
            eprint(f'\nERROR: -a mmseqs and -pp are incompatible',
                   flush = True)
            sys.exit(613)

    # choose the annotation software
    # NEED to implement interpro
    if args.pfam and args.interpro:
        eprint('\nERROR: --pfam and --interpro are incompatible',
               flush = True)
        sys.exit(921)
    elif args.pfam:
        execs.append('hmmsearch')
    if args.interpro:
        execs.append(format_path(args.interpro))

    # if the proposed alignment methodology fails, implement Diamond as a
    # back-up
    # NEED to remove MMseqs as an option
    if args.fallback:
        if 'diamond' not in set(execs):
            execs.append('diamond')

    # parse OrthoFinder data for the homology groups
    if args.orthofinder:
        of_out = format_path(args.orthofinder)
        if os.path.isdir(of_out):
            homogroups = of_out + '/Orthogroups/Orthogroups.txt'
            hg_dir = of_out + '/Orthogroup_Sequences/'
        else:
            homogroups = of_out
            hg_dir = os.path.dirname(of_out) + '../Orthogroup_Sequences/'
        if not os.path.isfile(hg_dir + 'OG0000000.fa'):
            hg_dir = format_path(args.hg_dir)
        method = 'orthofinder'
    # parse an inputted set of homology groups
    elif args.input:
        homogroups = format_path(args.input)
        hg_dir = format_path(args.hg_dir)
        method = 'mmseqs easy-cluster'
    # circumscribe homology groups using linclust
    elif args.linclust:
        method = 'mmseqs easy-linclust'
        hg_dir = format_path(args.hg_dir)
        homogroups = None
        execs.append('mmseqs')
    # circumscribe homology groups using easy-cluster
    else:
        method = 'mmseqs easy-cluster'
        hg_dir = format_path(args.hg_dir)
        homogroups = None
        execs.append('mmseqs')

    hgx_dir = format_path(args.hgx_dir)

    # set the minimum ID for considering something a homolog
    if args.minimum_gene_id < 30:
        eprint('\nERROR: -mg must be > 30', flush = True)
        sys.exit(107)

    # set the tune arguments if there is a tune file
    if args.tune:
        tune_file = format_path(args.tune)
        args.inflation_rnd1 = None
        args.inflation_rnd2 = None
        if not os.path.isfile(tune_file):
            eprint('\nERROR: -T does not exist', flush = True)
            sys.exit(103)
    else:
        tune_file = None

    # set the distance type used for microsynteny distance qunatitation
    if args.maximum_dist:
        dist = 'mmd'
        if args.unique_sp:
            eprint('\nERROR: -mmd is not compatible with -us', flush = True)
            sys.exit(104)
        uniq_sp = False
    else:
        dist = 'tmd'
        uniq_sp = args.unique_sp

    # set the method of generating null distributions
    if args.null_partitions and args.null_rank:
        eprint('\nERROR: -np and -nr are incompatible', flush = True)
        sys.exit(105)
    elif args.null_partitions:
        partition_file = format_path(args.null_partitions)
        partition_rank = None
        partition = partition_file
    elif args.null_rank:
        if args.null_rank.lower() not in set(null_ranks):
            eprint(f'\nERROR: invalid rank {args.null_rank.lower()}', flush = True)
            sys.exit(106)
        partition_rank = args.null_rank.lower()
        partition = partition_rank
        partition_file = None
    else:
        partition_rank, partition_file = None, None
        partition = None

#    if not os.path.isfile(pfam):
 #       print('\nERROR: invalid Pfam-A.hmm path', flush = True)
  #      sys.exit(4)

    # check the dependencies
    findExecs(execs, exit = set(execs))
#    if args.gcf_percentile < args.hgx_percentile and args.gcf_percentile:
 #       print('\nERROR: HGx percentile is greater than GCF percentile',
  #          flush = True)
   #     sys.exit(3)
    #elif not args.gcf_percentile:
     #   args.gcf_percentile = args.hgx_percentile

    # create/check the output directory
    if not args.output_dir:
        args.output_dir = format_path('./')
        out_dir = mkOutput(format_path(args.output_dir), 'cloci')
    elif os.path.isdir(args.output_dir):
        out_dir = format_path(args.output_dir)
    else:
        os.mkdir(format_path(args.output_dir + '/'))
        out_dir = format_path(args.output_dir)

    # set the topological constraint for microsynteny tree reconstruction
    if args.constraint:
        constraint_path = format_path(args.constraint)
    else:
        constraint_path = None

    # set the percentiles for HGps and HGxs
    if len(str(args.hgp_percentile)) > 1:
        hgp_perc = float('.' + str(args.hgp_percentile))
    else:
        hgp_perc = float('.0' + str(args.hgp_percentile))
#    if len(str(args.gcf_percentile)) > 1:
 #       clus_perc = float('.' + str(args.gcf_percentile))
  #  else:
   #     clus_perc = float('.0' + str(args.gcf_percentile))   
    if len(str(args.hgx_percentile)) > 1:
        hgx_perc = float('.' + str(args.hgx_percentile))
    else:
        hgx_perc = float('.0' + str(args.hgx_percentile))   

    # prepare the introduction dictionary
    args_dict = {
        'Homology groups': homogroups, 'Sequence clusters': method, 'MycotoolsDB': args.database, 
        'Focal genes': args.focal_genes, 
        'Microsynteny tree': format_path(args.tree), 'Toplogy constraint': constraint_path,
        'Root': root_txt, 
        'Pfam DB': format_path(args.pfam), 'InterProScan': format_path(args.interpro), 
        'Window': args.window*2+1,
        'Distance type': dist, 'Unique species': uniq_sp,
        'HGp Percentile': args.hgp_percentile, #'Precluster threshold': args.clus_threshold,
        'HGx percentile': args.hgx_percentile, 'HLG aligner': f'{aligner} {algn_sens}',
        'Similarity index': sim, 'Minimum gene id': args.minimum_gene_id,
        'Minimum locus id': args.minimum_loc_id,
        'Inflation 1': args.inflation_rnd1, 'Inflation 2': args.inflation_rnd2,
        'Tune clusters': tune_file, 'Minimum CSB': args.csb_threshold,
#        'GCF percentile': args.gcf_percentile, 
        'Minimum GCF %id': args.id_percent, 'Minimum GCF %pos': args.pos_percent,
        'PDS threshold': args.pds_threshold, 'GCL threshold': args.gcl_threshold,
        'Log Normal TMD threshold': args.md_threshold, 'Null partitions': partition,
        'Null samples': args.null_sample, #'Calculate dN/dS': args.dnds, 
        'Minimum N50': args.n50, # 'Family sensitivity': args.sensitivity,
        'Processors': args.cpus, 'Output directory': args.output_dir,
        'HG directory': hg_dir, 'HGx directory': hgx_dir,
        'Overwrite': bool(args.new)
        }
    start_time = intro('CLOCI', args_dict, 
                       'Zachary Konkel, Laura Kubatko, Jason Slot')
    date = datetime.strftime(start_time, '%Y%m%d')

    # choose focal genes for microsynteny tree reconstruction
    if args.focal_genes:
        with open(format_path(args.focal_genes), 'r') as raw:
            focal_genes = [x.rstrip() for x in raw.read().split() \
                           if x.rstrip()]
    else:
        focal_genes = []

    db = mtdb(format_path(args.database))
    main(
        db, homogroups, out_dir, plusminus = args.window,
        cpus = args.cpus, hg_dir = hg_dir, inflation_1 = args.inflation_rnd1,
        inflation_2 = args.inflation_rnd2,
        hgp_perc = hgp_perc, hgx_dir = hgx_dir, #clus_thresh = args.clus_threshold,
  #      clus_perc = clus_perc, aligner = aligner,#seed_thresh = args.seed_threshold,
        aligner = aligner, ipr_path = None, csb_thresh = args.csb_threshold,
        id_perc = args.id_percent, pos_perc = args.pos_percent,
        hgx_perc = hgx_perc, pfam = format_path(args.pfam), samples = args.null_sample,
        partition_file = args.null_partitions, #        run_dnds = args.dnds, 
        n50thresh = args.n50, near_single_copy_genes = focal_genes,
        root = root, gcl_thresh = args.gcl_threshold, 
        constraint_path = constraint_path, simfun = simfun,
        patch_thresh = args.pds_threshold, method = method,
        printexit = args.stop, skipalgn = args.skip, flag = bool(not args.new),
        min_hlg_id = args.minimum_loc_id / 100, sim = sim, 
        tree_path = format_path(args.tree), tune_file = tune_file,
        dist_thresh = args.md_threshold, uniq_sp = uniq_sp,
        dist_type = dist, partition_rank = partition_rank, 
        min_branch_sim = min_topology_sim, algn_sens = algn_sens,
        min_gene_id = args.minimum_gene_id, fallback = args.fallback,
        merge_via_sim = False, force = args.force #args.topology_merge
        )

    outro(start_time)


if __name__ == '__main__':
    cli()
