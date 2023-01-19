#! /usr/bin/env python3

#NEED to implement GCC network analysis from loci.adj
#NEED GBC threshold
#NEED a coverage filter for the diamond searches
#NEED to make gbc_mngr split different families with the same GCF
#NEED to allow hard threshold for microsynteny distance
#NEED locus-based GCF hypergeometric average
#NEED to accept an alternative OG input
#NEED to update logging and intelligent resuming
    # md5sum for database integrity
    # inflation parameter
    # minimum gcf id
#NEED TO OUTPUT FAMILY TO info.out
#NEED TO ADD PERCENTILE TO OUTPUT
#NEED to implement gff2svg
#NEED to output to gbk
#NEED to delete excess when checkpoints are reached
    # implement tarring effectively

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
from tqdm import tqdm
from itertools import chain
from datetime import datetime
from collections import defaultdict
from mycotools.lib.kontools import \
    intro, outro, format_path, collect_files, \
    findExecs, checkdir, eprint
from mycotools.lib.biotools import \
    gff2list
from mycotools.lib.dbtools import mtdb, masterDB
from mycotools.gff2svg import main as gff2svg
from mycotools import db2microsyntree
from orthocluster.orthocluster.lib import phylocalcs, evo_conco, \
     hgx2gcfs, input_parsing, hgp2hgx, generate_nulls
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


def patch_main(
    phylo, hgxs, wrk_dir, 
    old_path = 'patchiness.scores.pickle', cpus = 1
    ):

    if not os.path.isfile(wrk_dir + old_path):
        clusOmes = set([
            tuple([str(x) for x in y]) for y in hgxs
            ])
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            patch_res = pool.starmap(
                phylocalcs.calc_patchiness, tqdm([(phylo, x) for x in clusOmes],
                                                 total = len(clusOmes))
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


def main(
    db, hg_file, out_dir, plusminus = 1, hg_dir = None,
    hgp_perc = 0.2, clus_perc = 0.7, hgx_perc = 0.7,
    id_perc = 30, pos_perc = 30, inflation = 1.5,
    minimum_omes = 2, samples = 10000, pfam = None,
    min_gcf_id = 0.3, simfun = hgx2gcfs.overlap,
    constraint_path = None, blastp = 'blastp',
    run_dnds = False, cpus = 1, n50thresh = None, 
    root = None, coevo_thresh = 0, patch_thresh = 0,
    method = 'mmseqs easy-cluster',
    printexit = False, flag = True, partition_file = None,
    near_single_copy_genes = [], tree_path = None, 
    verbose = False, sim = 'overlap'
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
    db = db.set_index()
    if not tree_path:
        tree_path = f'{out_dir}microsynt.newick'

    wrk_dir, nul_dir = input_parsing.init_run(db, out_dir, 
                                              near_single_copy_genes, constraint_path,
                                              tree_path, hg_file, plusminus, 
                                              hgp_perc, clus_perc,
                                              hgx_perc, id_perc, pos_perc, 
                                              patch_thresh, coevo_thresh,
                                              samples, n50thresh, flag,
                                              min_gcf_id, inflation, sim)


    ome2i, gene2hg, i2ome, hg2gene, ome2pairs, cooccur_dict = \
        db2microsyntree.main(db, hg_file, out_dir, wrk_dir,
                            method, tree_path, plusminus = plusminus,
                            min_cov = 0, min_id = 0.3, n50thresh = n50thresh,
                            near_single_copy_genes = near_single_copy_genes,
                            constraint = constraint_path, verbose = verbose,
                            cpus = cpus)
    print('\tReading microsynteny tree', flush = True)
    phylo = input_parsing.compile_tree(
        i2ome, tree_path, root = root
        )
    partition_omes, ome2partition, omes2dist, min_pair_scores = \
                                generate_nulls.gen_pair_nulls(
                                            phylo, ome2i, wrk_dir,
                                            nul_dir, hgp_perc, ome2pairs,
                                            i2ome, samples = samples,
                                            partition_file = partition_file,
                                            cpus = cpus
                                            )

    # seed clusters by calcuating total microsynteny distance for 
    # orthogroup pairs
    print('\nIII. Seeding HG pairs (HGp)', flush = True) 
    if not os.path.isfile(out_dir + 'hgps.tsv.gz'):
        print('\tCalculating seed HG-pair scores', flush = True)
        seed_score_start = datetime.now()
        omes2dist = phylocalcs.update_dists(phylo, cooccur_dict, cpus, omes2dist)
        with open(wrk_dir + 'omes2tmd.pickle', 'wb') as out:
            pickle.dump(omes2dist, out)

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
            if score > min_pair_score:
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
    elif not os.path.isfile(wrk_dir + 'hgx2omes.pickle'): # load previous hg pairs
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
    ome_combos = set([tuple(sorted(list(x))) for x in list(hgx2omes.values())])

    print('\tCalculating HGx microsynteny distances', flush = True)
    hgx_start = datetime.now()
    hgx_obs = len(hgx2omes)
    print('\t\t' + str(hgx_obs) + ' observed HGx', flush = True)

    # calculate HGx microsynteny distances
    omes2dist = phylocalcs.update_dists(phylo, hgx2omes, cpus, omes2dist = omes2dist)
    with open(wrk_dir + 'omes2tmd.pickle', 'wb') as out:
        pickle.dump(omes2dist, out)
    print('\t\t' + str(datetime.now() - hgx_start), flush = True)


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
                                            db, partition_omes, ome2i, gene2hg,
                                            max_hgx_size, plusminus, hgx_perc, 
                                            clus_perc, nul_dir, omes2dist, 
                                            phylo, samples,
                                            cpus)


    # collect for normalizing relative to the whole dataset later
    # should I move into absolute space? it may just be better for comparing
    # datasets in the future and some algorithm is going to pick up that gene clusters
    # are in some absolute microsynteny distance, so long as its not normalized. 

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
        if dist >= bord_score:
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

    hgx2gbc, omes2patch = {}, {}
    hgx_dir = wrk_dir + 'hgx/'
    print('\nV. Calling gene clusters from HGxs', flush = True)
    if not os.path.isfile(wrk_dir + 'gcfs.pickle'): # need to add this to
    # log parsing
        # Group hgxs
        if not checkdir(hgx_dir, unzip = True, rm = True):
            os.mkdir(hgx_dir)
    
        gcfs, gcf_hgxs, gcf_omes, omes2dist = hgx2gcfs.classify_gcfs(
            hgx2loc, db, gene2hg, i2hgx, hgx2i, phylo, bord_scores_list, 
            ome2i, hgx2omes, hg_dir, hgx_dir, wrk_dir, ome2partition, hg2gene, 
            omes2dist = omes2dist, clusplusminus = plusminus, 
            inflation = inflation, min_gcf_id = min_gcf_id,
            min_omes = 2, cpus = cpus, simfun = simfun, printexit = printexit
            )
        with open(wrk_dir + 'gcfs.pickle', 'wb') as pickout:
            pickle.dump([gcfs, gcf_omes, gcf_hgxs], pickout)
        with open(wrk_dir + 'omes2tmd.pickle', 'wb') as pickout:
            pickle.dump(omes2dist, pickout)
    else: # or just load old data
        with open(wrk_dir + 'gcfs.pickle', 'rb') as raw:
            gcfs, gcf_omes, gcf_hgxs = pickle.load(raw)

    omes2dist = phylocalcs.update_dists(phylo, {i: omes for i, omes in enumerate(gcf_omes)}, 
                                        cpus, omes2dist)

    todel = []
    for i, omes in enumerate(gcf_omes):
        dist = omes2dist[omes]
        parts = set(ome2partition[x] for x in omes)
        if None in parts:
            parts = parts.remove(None)
            if not parts:
                continue
        # get the most threshold corresponding to the minimum sized HGx 
        # in the family, and the minimum lineages' value [liberal]
        max_len = max([len(x) for x in gcfs[i]])
        clus_score = min([clus_scores_list[v][max_len] \
                          for v in list(parts)])
        if not dist >= clus_score:
            todel.append(i)

    todel.reverse()
    for i in todel:
        del gcf_omes[i]
        del gcfs[i]
        del gcf_hgxs[i]
    print(f'\t{len(gcfs)} GCFs pass threshold', flush = True)


    runOmes = [
        omes for omes in gcf_omes \
        if omes not in omes2patch
        ] # omes without patchiness scores
    runHGxs = [
        hgx for i, hgx in enumerate(gcf_hgxs) \
        if gcf_omes[i] not in omes2patch
        ] # hgxs without patchiness scores
    print('\nVI. Quantifying GCF patchiness', flush = True)
    omes2patch = {**patch_main(
        phylo, runOmes, wrk_dir,
        old_path = 'patchiness.full.pickle', cpus = cpus
        ), **omes2patch} # could make more efficient by skipping redos


    print('\nVII. Quantifying GCF gene BLAST concordance (GBC)', flush = True)
    hgx2omes2gbc, hgx2omes2id, hgx2omes2pos = evo_conco.gbc_main(
                            hgx2loc, wrk_dir, ome2i, hg_dir, hgx_dir,
                            blastp, db, gene2hg, plusminus, hg2gene, 
                            old_path = 'hgx2omes2pos.full.pickle',
                            modules = gcfs, moduleHGxs = gcf_hgxs,
                            moduleOmes = gcf_omes, printexit = printexit,
                            cpus = cpus
                            )
    # this is the output directory for hgx2gcf and evo_conco
    hg_dir = f'{wrk_dir}hg/'


    if any(x > 0 for x in [coevo_thresh, patch_thresh, id_perc, pos_perc]):
        print('\tApplying thresholds', flush = True)
        print('\t\t' + str(len(gcfs)) + ' gcfs before', flush = True)
        newGCFs, newGCFOmes, newGCFHGxs = [], [], []
        for i0, gcf in enumerate(gcfs):
            check = False # have we added a new list
            ogc = gcf_hgxs[i0]
            omesc = gcf_omes[i0]
            for x, omes in gcf.items():
                if hgx2omes2gbc[ogc][omesc] >= coevo_thresh \
                    and omes2patch[omesc] >= patch_thresh \
                    and hgx2omes2id[ogc][omesc] >= id_perc \
                    and hgx2omes2id[ogc][omesc] >= pos_perc:
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


    hgx2dnds = {}
    gcf_dict = {i: v for i,v in enumerate(gcf_omes)}

    omes2dist = phylocalcs.update_dists(phylo, 
                                        {i: v for i, v in enumerate(gcf_omes)}, 
                                        cpus, omes2dist) 

    print('\nIIX. Writing and annotating clusters', flush = True)
    output_res(db, wrk_dir, hgx2dist, gcfs, gcf_omes, 
         i2ome, hgx2omes, out_dir, gcf_hgxs,
         omes2dist, hgx2omes2gbc, omes2patch, hgx2omes2id,
         hgx2omes2pos, hgx2loc, gene2hg, plusminus, ome2i,
         hgx2i, pfam_path = pfam, dnds_dict = {}, cpus = cpus)


if __name__ == '__main__':
    # need these here because spawn mp context forces reimport
#    from scipy.stats import hypergeom
 #   from dna_features_viewer import GraphicFeature, GraphicRecord
    
    description = \
    """Pharmaceuticals are primarily derived from biologically-produced compounds, their
    derivatives, and synthetic compounds inspired by natural pharmacophores.
    Many natural product specialized metabolites are produced by gene clusters, or 
    regions of genomes that contain colocalized genes with concerted 
    biosynthetic function. 
    Genomes rearrange and recombine through macroevolutionary timescales, so the 
    persistence of colocalized gene families within clusters is indicative of 
    selection. 
    hg2clus detects these regions by identifying unexpected conserved microsynteny 
    of homology groups throughout a sample of organisms. 
    Some unexpectedly conserved regions are not directly involved in specialized 
    metabolite biosynthesis - hg2clus filters results by implementing a
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

    mt_opt = parser.add_argument_group('Microsynteny tree')
    mt_opt.add_argument('-f', '--focal_genes',
        help = 'File of genes for neighborhood extraction of microsynteny ' \
             + 'tree')
    mt_opt.add_argument('-c', '--constraint',
        help = 'Constrain microsynteny topology to species tree w/ome code tips')
    mt_opt.add_argument('-r', '--root', 
        help = 'Ome code or 2 ome codes to root upon, e.g. psicya1, ' + \
        '"psicya1 psicub1"; DEFAULT: midpoint')
    mt_opt.add_argument('-t', '--tree',
        help = 'Precomputed microsynteny tree path')

    det_opt = parser.add_argument_group('Detection parameters')
    det_opt.add_argument('-of', '--orthofinder',
        help = 'Precomputed OrthoFinder output directory')
    det_opt.add_argument('-i', '--input',
        help = 'Precomputed cluster results file')
    det_opt.add_argument('-l', '--linclust', action = 'store_true',
        help = 'Less sensitive Homology inference via linclust; DEFAULT: mmseqs cluster')
    det_opt.add_argument('-w', '--window', default = 3, type = int,
        help = 'Max genes +/- for each locus window. DEFAULT: 3 (7 gene window)')
    det_opt.add_argument('-np', '--null_partitions', 
        help = 'Tab-delimited file with omes for each null sample group on ' \
             + 'separate lines')
    det_opt.add_argument('-ns', '--null_sample', type = int, default = 10000,
        help = 'Samples for null distributions; DEFAULT: 10,000')

    fam_opt = parser.add_argument_group('Gene cluster family parameters')
    fam_opt.add_argument('-s', '--similarity', 
        help = 'GCF similarity metric: [J]accard, [O]verlap, [S]orensen. ' \
             + 'DEFAULT: Overlap coefficient')
    fam_opt.add_argument('-m', '--minimum_id', type = float, default = 30,
        help = 'Percent [0 < value < 100] ID minimum between loci for GCF. ' \
             + 'DEFAULT: 30')
    fam_opt.add_argument('-I', '--inflation', default = 1.5, type = float,
        help = 'MCL inflation during family detection. DEFAULT: 1.5')
   

    thr_opt = parser.add_argument_group('Thresholding')
    thr_opt.add_argument('-hp', '--hgp_percentile', type = int, default = 75,
        help = 'Percentile of HG pair distances; DEFAULT: 75')
    thr_opt.add_argument('-xp', '--hgx_percentile', type = int,
        help = 'Percentile [0 < value < 100] of HGx microsynteny distances. ' \
             + 'Must be less than -cp; DEFAULT: -cp')
    thr_opt.add_argument('-gp', '--gcf_percentile', type = int, default = 70, 
        help = 'Percentile [0 < value < 100] of GCF microsynteny distances. ' \
             + 'DEFAULT: 70')
    thr_opt.add_argument('-ip', '--id_percent', default = 0, type = float,
        help = 'Percent [0 < value < 100] identity minimum for gene cluster family')
    thr_opt.add_argument('-pp', '--pos_percent', default = 0, type = float,
        help = 'Percent [0 < value < 100] positive minimum for gene cluster family')
    thr_opt.add_argument('-pt', '--patch_threshold', default = 0, type = float,
        help = "Threshold [0 < value < 1] of gene cluster family " \
             + " patchiness scores")
    thr_opt.add_argument('-ct', '--coevo_threshold', default = 0, type = float, 
        help = "Threshold [0 < value < 1] of high order HG combinations' coevolution scores")

    run_opt = parser.add_argument_group('Runtime options')
#    run_opt.add_argument('-s', '--dnds', action = 'store_true', help = 'Run dN/dS calculations')
    run_opt.add_argument('--n50', help = 'Minimum assembly N50')
    run_opt.add_argument('--stop', action = 'store_true', 
                        help = 'Export homology group diamond commands ' \
                             + 'for massive parallelization and stop')
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
        args.hgx_percentile = args.gcf_percentile # set the default

    if not args.similarity:
        simfun = hgx2gcfs.overlap # set default function
        sim = 'overlap'
    elif args.similarity.lower() in {'o', 'overlap', 'oc'}:
        simfun = hgx2gcfs.overlap
        sim = 'overlap'
    elif args.similarity.lower() in {'j', 'jaccard', 'jacard', 'jac'}:
        simfun = hgx2gcfs.jaccard
        sim = 'Jaccard'
    elif args.similarity.lower() in {'s', 'sorensen', 'sorenson', 'sor'}:
        simfun = hgx2gcfs.sorensen
        sim = 'Sorensen'
    else:
        eprint(f'\nERROR: invalid -s: {args.similarity}',
               flush = True)
        sys.exit(43)

    execs = ['mafft', 'hmmsearch', 'blastp', 'mcl',
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
        'Microsynteny tree': format_path(args.tree), 'Root': root_txt, 'Pfam DB': args.pfam, 'Window': args.window*2+1,
        'HGp Percentile': args.hgp_percentile, #'Precluster threshold': args.clus_threshold,
        'HGx percentile': args.hgx_percentile, 'Similarity index': sim,
        'Inflation': args.inflation,
        'GCF percentile': args.gcf_percentile, 
        'Minimum family %id': args.id_percent, 'Minimum family %pos': args.pos_percent,
        'Patchiness threshold': args.patch_threshold, 'Coevolution threshold': args.coevo_threshold,
        'Null samples': args.null_sample, #'Calculate dN/dS': args.dnds, 
        'Minimum N50': args.n50, # 'Family sensitivity': args.sensitivity,
        'Processors': args.cpus, 'Output directory': args.output_dir,
        'Overwrite': bool(args.new)
        }

    pfam = format_path(args.pfam)
#    if not os.path.isfile(pfam):
 #       print('\nERROR: invalid Pfam-A.hmm path', flush = True)
  #      sys.exit(4)

    findExecs(execs, exit = set(execs))
    if args.gcf_percentile < args.hgx_percentile:
        print('\nERROR: HGx percentile is greater than GCF percentile',
            flush = True)
        sys.exit(3)

    start_time = intro('homologs to clusters - hg2clus', args_dict, 
                       'Zachary Konkel, Laura Kubatko, Jason Slot')
    date = datetime.strftime(start_time, '%Y%m%d')

    if not args.output_dir:
        out_dir = os.getcwd() + '/hg2clus_' + date + '/'
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
    else:
        out_dir = format_path(args.output_dir)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
            out_dir += '/'

    if args.constraint:
        constraint_path = format_path(args.constraint)
    else:
        constraint_path = None
    if len(str(args.hgp_percentile)) > 1:
        hgp_perc = float('.' + str(args.hgp_percentile))
    else:
        hgp_perc = float('.0' + str(args.hgp_percentile))
    if len(str(args.gcf_percentile)) > 1:
        clus_perc = float('.' + str(args.gcf_percentile))
    else:
        clus_perc = float('.0' + str(args.gcf_percentile))   
    if len(str(args.hgx_percentile)) > 1:
        hgx_perc = float('.' + str(args.hgx_percentile))
    else:
        hgx_perc = float('.0' + str(args.hgx_percentile))   

    if args.focal_genes:
        with open(format_path(args.focal_genes), 'r') as raw:
            focal_genes = [x.rstrip() for x in raw.read().split() \
                           if x.rstrip()]

    sensitive = False

    db = mtdb(format_path(args.database))
    main(
        db, homogroups, out_dir, plusminus = args.window,
        cpus = args.cpus, hg_dir = hg_dir, inflation = args.inflation,
        hgp_perc = hgp_perc, #clus_thresh = args.clus_threshold,
        clus_perc = clus_perc, blastp= 'blastp',#seed_thresh = args.seed_threshold,
        id_perc = args.id_percent, pos_perc = args.pos_percent,
        hgx_perc = hgx_perc, pfam = pfam, samples = args.null_sample,
        partition_file = args.null_partitions, #        run_dnds = args.dnds, 
        n50thresh = args.n50, near_single_copy_genes = focal_genes,
        root = root, coevo_thresh = args.coevo_threshold, 
        constraint_path = constraint_path, simfun = simfun,
        patch_thresh = args.patch_threshold, method = method,
        printexit = args.stop, flag = bool(not args.new),
        min_gcf_id = args.minimum_id / 100, sim = sim, 
        tree_path = format_path(args.tree)
        )

    outro(start_time)
