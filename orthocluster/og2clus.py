#! /usr/bin/env python3

#NEED thresholding for percent pos and percent id
#NEED a coverage filter for the diamond searches
#NEED to make gbc_mngr_3 split different families with the same GCF
#NEED to change border percentile entries to clus percentile
#NEED to allow hard threshold for microsynteny distance
#NEED locus-based GCF hypergeometric average
#NEED to accept an alternative OG input
#NEED to update logging and intelligent resuming
    # md5sum for database integrity
#NEED TO OUTPUT FAMILY TO info.out
#NEED TO ADD PERCENTILE TO OUTPUT
#NEED to implement gff2svg
#NEED to delete excess when checkpoints are reached
    # implement tarring effectively
#NEED to use conservative thresholds depending on partition

import os
import re
import sys
import copy
import gzip
import shutil
import pickle
import argparse
import numpy as np
import multiprocessing as mp
from datetime import datetime
from collections import defaultdict
from mycotools.lib.kontools import \
    intro, outro, format_path, collect_files, \
    findExecs, multisub, checkdir, eprint, \
    write_json, read_json
from mycotools.lib.biotools import \
    gff2list, dict2fa
from mycotools.lib.dbtools import mtdb, masterDB
from mycotools.gff2svg import main as gff2svg
from mycotools.acc2fa import dbmain as acc2fa
from orthocluster.orthocluster.tools import db2microsynt
from orthocluster.orthocluster.lib.generate_nulls import \
    gen_pair_nulls, partitions2hgx_nulls
from orthocluster.orthocluster.lib import phylocalcs, evo_conco, \
    hgx2gcfs
from orthocluster.orthocluster.lib.output_data import output_res


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


def load_seedScores(file_, seed_thresh):#, seed_thresh):

    out_hgs = []
    with gzip.open(file_, 'rt') as raw:
        for line in raw:
            if not line.startswith('#'):
                data = [x.rstrip() for x in line.split('\t')]
                if float(data[3]) > seed_thresh:
                    out_hgs.append(line.split('\t'))

    return out_hgs


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
        nulls = collect_files(wrk_dir + 'null/', 'null.txt')
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
                phylocalcs.calc_patchiness, [(phylo, x) for x in clusOmes]
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


def main(
    db, hg_file, out_dir, plusminus = 1, hg_dir = None,
    seed_perc = 0.2, clus_perc = 0.7, hgx_perc = 0.7,
    minimum_omes = 2, samples = 10000, pfam = None,
    constraint_path = None, blastp = 'blastp',
    run_dnds = False, cpus = 1, n50thresh = None, 
    root = None, coevo_thresh = 0, patch_thresh = 0,
    microsyn_thresh = 0, method = 'mmseqs easy-cluster',
    printexit = False, flag = True, partition_file = None,
    near_single_copy_genes = [], verbose = False
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
    db = db.set_index('ome')
    wrk_dir = out_dir + 'working/'
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)

    nul_dir = wrk_dir + 'null/'
    if not os.path.isdir(nul_dir):
        os.mkdir(nul_dir)

    log_path = out_dir + 'log.txt'
    log_dict = {
        'hg_file': hg_file, 'plusminus': plusminus,
        'pair_percentile': seed_perc, 'hgx_percentile': clus_perc,
        'border_percentile': hgx_perc,
        'null_samples': samples, 'n50': n50thresh
        }
    log_res = logCheck(log_dict, log_path, out_dir, wrk_dir, flag)
    tree_path = out_dir + 'microsynt.newick'

    ome2i, gene2hg, i2ome, hg2gene, ome2pairs, cooccur_array =  \
        db2microsynt.main(db, hg_file, out_dir, wrk_dir,
                            method, tree_path, plusminus = plusminus,
                            min_cov = 0, min_id = 0.3, n50thresh = n50thresh,
                            near_single_copy_genes = near_single_copy_genes,
                            constraint = constraint_path, verbose = verbose,
                            cpus = cpus)
    print('\tReading microsynteny tree', flush = True)
    phylo = phylocalcs.compileTree(
        i2ome, tree_path, root = root
        )

    partition_phylos, partition_omes, ome2partition, \
    omes2dist, min_pair_scores = generate_nulls.gen_pair_nulls(
                                            phylo, i2ome, wrk_dir,
                                            nul_dir, seed_perc,
                                            samples = samples,
                                            partition_file = partition_file,
                                            cpus = cpus
                                            )

    # seed clusters by calcuating total microsynteny distance for 
    # orthogroup pairs
    print('\nIII. Seeding HG pairs', flush = True) 
    if not os.path.isfile(out_dir + 'seed_scores.tsv.gz'):
        print('\tCalculating seed HG-pair scores', flush = True)
        seed_score_start = datetime.now()
        results = phylocalcs.calc_dists(phylo, cooccur_dict, cpus)
        omes2dist, top_hgs = {x[1]: x[0] for x in results}, []
        with open(wrk_dir + 'ome_scores.pickle', 'wb') as out:
            pickle.dump(omes2dist, out)
        for hgpair, omes in cooccur_dict.items():
            score = omes2dist[omes]
            parts = set(ome2partition[x] for x in omes)
            if None in parts:
                parts = parts.remove(None)
                if not parts:
                    continue
            min_pair_score = min([min_pair_scores[i] for i in list(parts)])
            if score > min_pair_score:
                i = revs[hgpair]
                top_hgs.append([
                    [0], hgpair[1], len(omes), score#, score/ome_dist
                    ])

        # write scores            
        print('\t\t' + str(datetime.now() - seed_score_start), flush = True)
        top_hgs = sorted(top_hgs, key = lambda x: x[3])
        print('\tWriting seed scores', flush = True)
        with gzip.open(out_dir + 'seed_scores.tsv.gz', 'wt') as out:
            out.write('#hg0\thg1\tcooccurrences\tscore\tadj_score\n')
            for line in top_hgs:
                out.write('\t'.join([str(x) for x in line]) + '\n')
        print('\t\t' + str(len(top_hgs)) + ' significant seeds', flush = True)
    
    elif not os.path.isfile(wrk_dir + 'hgx_scores.pickle'): # load previous og pairs
        print('\tLoading previous seed HG pairs', flush = True)
        top_hgs = load_seedScores(out_dir + 'seed_scores.tsv.gz', min_score)
        print('\t\t' + str(len(top_hgs)) + ' significant seeds', flush = True)


    # begin sifting for HGxs using pairs as seeds for HGx detection
    print('\nIV. Sprouting high order HG combinations (HGx)', flush = True)
    hgx2omes, hgx2loc = hgpairs2hgx.hgpairs2hgx(db, wrk_dir, top_hgs,
                                                hgpair_dict, gene2hg, 
                                                ome2i, omes2dist, phylo, 
                                                plusmimus, cpus) 
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
        results = phylocalcs.calc_dists(phylo, hgx_cooccur_dict, cpus, omes2dist = omes2dist)
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
    bordScores_list, clusScores_list = generate_nulls.partitions2hgx_nulls(
                                            db, partition_omes, ome2i, gene2hg,
                                            max_hgx_size, plusminus, hgx_perc, 
                                            clus_perc, nul_dir, omes2dist, samples,
                                            cpus)


    # collect for normalizing relative to the whole dataset later
    # should I move into absolute space? it may just be better for comparing
    # datasets in the future and some algorithm is going to pick up that gene clusters
    # are in some absolute microsynteny distance, so long as its not normalized. 
#    max_dist = max(hgx2dist.values())
 #   min_dist = min(hgx2dist.values())

    i2hgx, hgx2i, thgx2dist, count = {}, {}, {}, 0
    for hgx, dist in hgx2dist.items():
        omes = hgx2omes[hgx]
        parts = set(ome2partition[x] for x in omes)
        if None in parts:
            parts = parts.remove(None)
            if not parts:
                continue
        bord_score = min([bordScores[i][len(hgx)] for i in list(parts)])
        if dist >= bord_score:
             thgx2dist[hgx] = dist
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

    print('\nV. Quantifying HGx coevolution', flush = True)
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

    hgx2gbc = evo_conco.gbc_main_0(
        hgx2loc, wrk_dir, ome2i, hg_dir, hgx_dir,
        'diamond', db, gene2hg, plusminus, hg2gene, 
        old_path = 'hgx2gbc.pickle', cpus = cpus, printexit = printexit
        ) # should be able to skip this

    # Group hgxs
    print('\nVI. Inferring HGx clans', flush = True)
    if not os.path.isfile(wrk_dir + 'gcfs.pickle'): # need to add this to
    # log parsing
        gcfs, gcf_hgxs, gcf_omes, omes2dist = hgx2gcfs.classify_gcfs(
            hgx2loc, db, gene2hg, i2hgx, hgx2i,
            phylo, bordScores, ome2i,
            hgx2omes, wrk_dir, ome2partion, #hgx2dist,
            omes2dist = omes2dist, clusplusminus = plusminus, 
            min_omes = 2, cpus = cpus
            )
        with open(wrk_dir + 'gcfs.pickle', 'wb') as pickout:
            pickle.dump([gcfs, gcf_omes, gcf_hgxs], pickout)
        with open(wrk_dir + 'ome_scores.pickle', 'wb') as pickout:
            pickle.dump(omes2dist, pickout)
    else: # or just load old data
        with open(wrk_dir + 'gcfs.pickle', 'rb') as raw:
            gcfs, gcf_omes, gcf_hgxs = pickle.load(raw)

    runOmes = [
        omes for omes in gcf_omes \
        if omes not in omes2patch
        ] # omes without patchiness scores
    runHGxs = [
        hgx for i, hgx in enumerate(gcf_hgxs) \
        if gcf_omes[i] not in omes2patch
        ] # hgxs without patchiness scores
   
    print('\nVII. Quantifying GCF patchiness', flush = True)
    omes2patch = {**patch_main(
        phylo, runHGxs, runOmes, wrk_dir,
        old_path = 'patchiness.full.pickle', cpus = cpus
        ), **omes2patch} # could make more efficient by skipping redos

    hgx_dir = wrk_dir + 'hgx/'
    if not checkdir(hgx_dir, unzip = True, rm = True):
        os.mkdir(hgx_dir)

    print('\nIIX. Quantifying GCF gene evolution concordancy', flush = True)
    hgx2omes2gbc, hgx2omes2id, hgx2omes2pos = evo_conco.gbc_main_1(
                            hgx2loc, wrk_dir, ome2i, hg_dir, hgx_dir,
                            blastp, db, gene2hg, plusminus, hg2gene, 
                            old_path = 'hgx2omes2pos.full.pickle',
                            modules = gcfs, moduleHGxs = gcf_hgxs,
                            moduleOmes = gcf_omes, cpus = cpus
                            )

    if coevo_thresh > 0 or patch_thresh > 0 or microsyn_thresh > 0:
        print('\tApplying thresholds', flush = True)
        print('\t\t' + str(len(gcfs)) + ' gcfs before', flush = True)
        # edit gcfs, hgx2dist
 #       newOgx2dist, newGCFs, newGCFOmes, newGCFHGxs = {}, [], [], []
        newGCFs, newGCFOmes, newGCFHGxs = [], [], []
        for i0, gcf in enumerate(gcfs):
            check = False # have we added a new list
            ogc = gcf_hgxs[i0]
            omesc = gcf_omes[i0]
            for x, omes in gcf.items():
                if hgx2omes2gbc[ogc][omesc] >= coevo_thresh \
                    and omes2patch[omesc] >= patch_thresh:
                    if check:
                        newGCFs[-1][x] = omes
                    else:
                        newGCFs.append({x: omes})
                        check = True
            if check:
                newGCFOmes.append(gcf_omes[i0])
                newGCFHGxs.append(ogc)
        gcfs, gcf_omes, gcf_hgxs = \
            newGCFs, newGCFOmes, newGCFHGxs
        print('\t\t' + str(len(gcfs)) + ' GCFs after', flush = True)

#    if run_dnds: # need to bring file naming to speed
#    else:
 #       print('\nIIX. Skipping quantifying selection', flush = True)
    hgx2dnds = {}

    print('\nIX. Writing and annotating clusters', flush = True)
    output_res(hgx2dist, gcfs, i2ome, hgx2omes, out_dir, gcf_hgxs,
         omes2dist, hgx2omes2gbc, omes2patch, hgx2omes2id,
         hgx2omes2pos, hgx2loc, gene2hg, plusminus, ome2i,
         hgx2i, hgx2gbc, hgx2omes, dnds_dict = {}, cpus = cpus)


if __name__ == '__main__':
    # need these here because spawn mp context forces reimport
    from scipy.stats import hypergeom
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
    i_opt.add_argument('--pfam', help = 'Pfam-A.hmm')

    mt_opt = parser.add_argument_group('Tree')
    mt_opt.add_argument('-f', '--focal_genes',
        help = 'File of genes for neighborhood extraction of microsynteny ' \
             + 'tree')
    mt_opt.add_argument('-c', '--constraint',
        help = 'Constrain microsynteny topology to species tree w/ome code tips')
    mt_opt.add_argument('-r', '--root', 
        help = 'Ome code or 2 ome codes to root upon, e.g. psicya1, ' + \
        '"psicya1 psicub1"; DEFAULT: midpoint')

    det_opt = parser.add_argument_group('Detection parameters')
    det_opt.add_argument('-of', '--orthofinder',
        help = 'Precomputed OrthoFinder output directory')
    det_opt.add_argument('-i', '--input',
        help = 'Precomputed cluster results file')
    det_opt.add_argument('-l', '--linclust', action = 'store_true',
        help = 'Less sensitive Homology inference via linclust; DEFAULT: mmseqs cluster')
    det_opt.add_argument('-w', '--window', default = 5, type = int,
        help = 'Max genes +/- for each locus window. DEFAULT: 5 (11 gene window)')
    det_opt.add_argument('-np', '--null_partitions', 
        help = 'Tab-delimited file with omes for each null sample group on ' \
             + 'separate lines')
    det_opt.add_argument('-ns', '--null_sample', type = int, default = 10000,
        help = 'Samples for null distributions; DEFAULT: 10,000')

    thr_opt = parser.add_argument_group('Thresholding')
    thr_opt.add_argument('-sp', '--seed_percentile', type = int, default = 75,
        help = 'Percentile of HG pair distances; DEFAULT: 75')
    thr_opt.add_argument('-op', '--hgx_percentile', type = int,
        help = 'Percentile of HGx microsynteny distances. ' + \
        'Must be less than -cp; DEFAULT: -cp')
    thr_opt.add_argument('-cp', '--clus_percentile', type = int, default = 80, 
        help = 'Percentile of HGx microsynteny distances; DEFAULT: 80')
    thr_opt.add_argument('-pt', '--patch_threshold', default = 0, type = float,
        help = "Threshold [0 < value < 1] of high order OG combinations' patchiness scores")
    thr_opt.add_argument('-ct', '--coevo_threshold', default = 0, type = float, 
        help = "Threshold [0 < value < 1] of high order OG combinations' coevolution scores")

    run_opt = parser.add_argument_group('Runtime options')
#    run_opt.add_argument('-s', '--dnds', action = 'store_true', help = 'Run dN/dS calculations')
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
             'diamond', 'iqtree']
    if args.orthofinder:
        of_out = format_path(args.orthofinder)
        if os.path.isdir(of_out):
            homogroups = of_out + '/Orthogroups/Orthogroups.txt'
            hg_dir = of_out + '/Orthogroup_Sequences/'
        else:
            homogroups = of_out
            hg_dir = os.path.dirname(of_out) + '../Orthogroup_Sequences/'
        if not os.path.isfile(hg_dir + 'OG0000000.fa'):
            hg_dir = None
        method = 'orthofinder'
    elif args.input:
        homogroups = format_path(args.input)
        hg_dir = None
        method = 'mmseqs easy-cluster'
    elif args.linclust:
        method = 'mmseqs easy-linclust'
        hg_dir = None
        homogroups = None
        execs.append('mmseqs')
    else:
        method = 'mmseqs easy-cluster'
        hg_dir = None
        homogroups = None
        execs.append('mmseqs')

    args_dict = {
        'Homology groups': homogroups, 'Sequence clusters': method, 'MycotoolsDB': args.database, 
        'Root': root_txt, 'Pfam DB': args.pfam, 'Window': args.window*2+1,
        'Seed Percentile': args.seed_percentile, #'Precluster threshold': args.clus_threshold,
        'Cluster percentile': args.clus_percentile, 'HGx percentile': args.hgx_percentile,
        'Patchiness threshold': args.patch_threshold, 'Coevolution threshold': args.coevo_threshold,
        'Null samples': args.null_sample, #'Calculate dN/dS': args.dnds, 
        'Minimum N50': args.n50,
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

    start_time = intro('homogroups to clusters - og2clus', args_dict, 
                       'Zachary Konkel, Laura Kubatko, Jason Slot')
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

    if args.focal_genes:
        with open(format_path(args.focal_genes), 'r') as raw:
            focal_genes = [x.rstrip() for x in raw.read().split() \
                           if x.rstrip()]

    db = mtdb(format_path(args.database))
    main(
        db, homogroups, out_dir, plusminus = args.window,
        cpus = args.cpus, hg_dir = hg_dir, 
        seed_perc = seed_perc, #clus_thresh = args.clus_threshold,
        clus_perc = clus_perc, blastp= 'blastp',#seed_thresh = args.seed_threshold,
        hgx_perc = hgx_perc, pfam = pfam, samples = args.null_sample,
#        run_dnds = args.dnds, 
        n50thresh = args.n50, near_single_copy_genes = focal_genes,
        root = root, coevo_thresh = args.coevo_threshold, 
        patch_thresh = args.patch_threshold, method = method,
        printexit = args.stop, flag = bool(not args.new)
        )

    outro(start_time)
