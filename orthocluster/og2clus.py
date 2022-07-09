#! /usr/bin/env python3

#NEED to make patchiness/coevolution based on OCGs
#NEED to model rearrangement for lineages instead of unexpected microsynteny
#NEED to change border percentile entries to clus percentile
#NEED to allow hard threshold for microsynteny distance
#NEED corrected alpha threshold for hypergeometric
#NEED locus-based family hypergeometric average
#NEED to post filter clusters after grouping
#NEED to accept an alternative OG input
#NEED an intelligent resume for patch and coevo percentiles
#NEED to update logging
#NEED TO OUTPUT FAMILY TO info.out
#NEED TO ADD PERCENTILE TO OUTPUT
#NEED to implement gff2svg
#NEED to delete excess when checkpoints are reached


import plotly.graph_objects as go
from plotly.subplots import make_subplots
import gzip, os, shutil, random, argparse, sys, subprocess, time, re, pickle, copy
import multiprocessing as mp, numpy as np
from itertools import combinations
from mycotools.lib.kontools import \
    intro, outro, formatPath, collect_files, \
    dictSplit, file2list, fmt_float, findExecs, \
    multisub, tardir, untardir, checkdir, eprint, \
    writeJson, readJson
from mycotools.lib.biotools import \
    gff2list, list2gff, fa2dict, dict2fa
from mycotools.lib.dbtools import mtdb, masterDB
from mycotools.gff2seq import ntmain
from mycotools.acc2gff import grabGffAcc
from mycotools.gff2svg import main as gff2svg
from mycotools.acc2fa import dbmain as acc2fa
from datetime import datetime
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
from collections import defaultdict

def calcBranchLen(phylo, omes):
    '''calculate descending branch length from cogent3 tree'''
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

def updateScores(phylo, cooccur_dict, cpus = 1, omes2dist = {}):
    '''update the omes2dist with a new set of data'''
    results = calcDists(phylo, cooccur_dict, cpus, omes2dist = omes2dist)
    omes2dist = {**omes2dist, **{x[1]: x[0] for x in results}} 
    return omes2dist


def parse_ogs(og_file, useableOmes = set()):
    '''
    imports orthofinder Orthogroups.txt "og_file". outputs several data structures:
    ome_num = {ome: number}, gene2og = {gene: og}, i2ome = [ome0, ome1, ome2]
    '''

    gene2og, ome_num, i2ome, og2gene = \
        {}, {}, [], {}
    with open(og_file, 'r') as raw:
        for line in raw:
            data = line.rstrip().split(' ') # OrthoFinder input
            og = int(data[0].replace(':','').replace('OG',''))
            hits = [x for x in data[1:] if x[:x.find('_')] in useableOmes]
            omes = [x[:x.find('_')] for x in hits]
            if len(set(omes)) < 2: # skip singleton organisms
                continue
            og2gene[og] = hits
            for i, gene in enumerate(hits):
                ome = omes[i]
#                ome = gene[:gene.find('_')] # wtf, parsing omes raises and index error
                if ome not in ome_num:
                    i2ome.append(ome)
                    ome_num[ome] = len(i2ome) - 1
                gene2og[gene] = og

    return ome_num, gene2og, i2ome, og2gene


def compileOGs(og_file, wrk_dir, useableOmes = set()):
    ogInfo = parse_ogs(og_file, useableOmes)
    with open(wrk_dir + 'ome2i.tsv', 'w') as out:
        out.write(
            '\n'.join([x + '\t' + str(ogInfo[0][x]) for x in ogInfo[0].keys()])
            )
    max_ome =  max([int(x) for x in ogInfo[1].values()])
    return ogInfo


def compileLoci(
    db, ome2i, gene2og, plusminus, cpus = 1
    ):

    loci_hash_cmds = [
        [x, ome2i[x.replace('.gff3','')], 
        gene2og, plusminus]
        for x in list(db['gff3']) if x.replace('.gff3','') in ome2i
        ]
    with mp.Pool(processes = cpus) as pool:
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

def compileTree(microsynt_dict, i2ome, tree_path, root = []):
    if microsynt_dict: # implies the tree needs to be made from scratch
        pre_phylo = nj.nj(microsynt_dict, show_progress = False)
        # change to ome code for exporting phylogeny
        pre_phylo.reassign_names({str(i): v for i, v in enumerate(i2ome)})

        if not root or root == 'midpoint':
            root_pre = pre_phylo.root_at_midpoint()
        elif len(root) == 1:
            root_pre = pre_phylo.rooted_with_tip(root[0])
        else:
            root_pre = pre_phylo.rooted_at(
                pre_phylo.get_edge_names(root[0], root[1])[0]
                )
        root_pre.write(tree_path, with_distances = True)

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
    '''
    Inputs the gff_list and organism ome code. Compiles CDS entries for loci parsing and outputs:
    cds_dict = {contig: {protein: [[CDSstart, CDSstop]] } }
    Sorts the cds_dict for each protein from smallest to largest
    '''

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
    '''
    Inputs the gff_list and organism ome code. Compiles CDS entries for loci parsing and outputs:
    cds_dict = {contig: {protein: [[CDSstart, CDSstop]] } }
    Sorts the cds_dict for each protein from smallest to largest
    '''

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
    gff_path, ome_num, gene2og, plusminus = 6
    ):
    '''obtain a set of tuples of OG pairs {(OG0, OG1)...}'''

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
                    loc_og.append(gene2og[gene])
                except KeyError: # it's either a singleton or not in an OG
                    pass
            pairs.extend([tuple(sorted(x)) for x in combinations(set(loc_og), 2)]) # don't deal with
            # same OGs - tandem duplications will just naturally cooccur more often

    out_pairs = set(pairs) # unique pairs of OGs

    return ome_num, out_pairs

def FormCooccurDict(cooccur_dict):

    # sort the values by the length of the combination
    cooccur_dict = {k: tuple(sorted(v)) for k, v in sorted(cooccur_dict.items(), key = lambda x: len(x[0]))}

    ogx2i, i2ogx = {}, {}
    for i, ogX in enumerate(list(cooccur_dict.keys())):
        i2ogx[i] = ogX
        ogx2i[ogX] = i

    return i2ogx, ogx2i, cooccur_dict


def FormCooccurArray(cooccur_dict, ome_len):

    count, ogx2i, size_dict, cooccur_arrays, i2ogx = 0, {}, {}, {}, {}
    cooccur_dict = {k: tuple(sorted(v)) for k, v in sorted(cooccur_dict.items(), key = lambda x: len(x[0]))}
    cooccur_arrays = np.zeros([ome_len, len(cooccur_dict)], dtype = np.int32)

    old_len = len(list(cooccur_dict.keys())[0])
    for i, ogX in enumerate(list(cooccur_dict.keys())):
        i2ogx[i] = ogX
        ogx2i[ogX] = i
        for ome in cooccur_dict[ogX]:
#            cooccur_arrays[ome, i] = ogX_dict[ome][ogX]
            cooccur_arrays[ome, i] = 1

    return i2ogx, ogx2i, cooccur_dict, cooccur_arrays


def PrepMicrosyntDict(cooccur_arr):

 #   microsynt_arr = np.zeros((cooccur_arr.shape[0], cooccur_arr.shape[0]))
    microsynt_dict = {}
    i0 = 0
    for row0 in cooccur_arr[:, ]:
        microsynt_dict[(i0, i0)] = 0
        i1 = i0 + 1
        for row1 in cooccur_arr[i1:, ]:
            row_sum = row0 + row1
            dist = 1 - (len(row_sum[row_sum > 1]) / len(row0[row0 > 0]))
            microsynt_dict[(i0, i1)] = dist
#            microsynt_arr[(i1, i0)] = dist
            i1 += 1
        i0 += 1

    return microsynt_dict


def CalcHypergeo(
    ogx, omes, genesInOGinOmes, genesInOmes, window
    ):
    # NEED to modify to account for sliding window

    # should this be for all omes, or each og in each ome and then multiplied?
    pval = hypergeom.sf(len(omes), genesInOmes, genesInOGinOmes, window)
    return ogx, pval


def HypergeoMNGR(
    ogx2omes, omes2og2genes, omes2genes, window, cpus = 1
    ):

    cmds = []
    for ogx, omes in ogx2omes.items(): # take each ogx
        for og in ogx:
            genesInOGinOmes = sum([
                len(omes2og2genes[ome][og]) for ome in omes
                ])
            genesInOmes = sum([
                len(omes2genes[ome]) for ome in omes
                ])
            cmds.append([
                ogx, omes, genesInOGinOmes, genesInOmes, window
                ])

    # run og-by-og, then accumulate via ogx at the end
    with mp.Pool(processes = cpus) as pool:
        hypergeoRes = pool.starmap(CalcHypergeo, cmds)

    comparisons = len(ogx2omes) # for Bonferroni correction
    ogx2pval = {}
    for ogx, pval in hypergeoRes:
        if ogx not in ogx2pval:
            ogx2pval[ogx] = comparisons
        ogx2pval[ogx] *= pval

    return ogx2pval


def gen_null_dict(combo_dict, sample = 10000):
    '''combo_dict = {ome: [(OG0...OGn)]}'''

    nulls = []
    for combos in list(combo_dict.values()):
        nulls.extend(list(combos)) # extend all combos to the null

    if not isinstance(combo_dict[list(combo_dict.keys())[0]], set):
        combo_dict = {k: set(v) for k, v in combo_dict.items()}


    null_set_list = list(set(nulls)) # acquire unique combinations
    try: # sample
        null_list = random.sample(null_set_list, sample)
    except ValueError:
        print('\t\tWARNING: requested null sample greater than pairs/triads', flush = True)
        null_list = null_set_list

    cooccur_dict = {}
    for null in null_list: # for combo
        cooccur_dict[null] = []
        for ome, combos in combo_dict.items():
            if null in combos:
                cooccur_dict[null].append(ome)

    cooccur_dict = {x: cooccur_dict[x] for x in cooccur_dict}
    # {(OG0,OGn): [omei0...]}

    ogx2i, i2ogx, cooccur_dict = FormCooccurDict(cooccur_dict)
    return ogx2i, i2ogx, cooccur_dict


def calcDists(phylo, cooccur_dict, cpus = 1, omes2dist = {}):
    # multiprocessing calculating only new distances for omes2dist
    with mp.Pool(processes = cpus) as pool:
        results = pool.starmap(
            calcBranchLen, 
            [(phylo, x) for x in list(set(cooccur_dict.values())) \
            if x not in omes2dist]
            )

    return results


def gen_nulls(pairs, phylo, samples = 10000, cpus = 1):
    '''pairs = {ome: [(OG0,OG1)...]}'''

    ogpairs, revogpairs, null_dict = gen_null_dict(pairs, samples)
    oldLen = len(null_dict)
    null_dict = {k: v for k, v in null_dict.items() if len(v) > 1}
    results = calcDists(phylo, null_dict, cpus = cpus)
    omes2dist, pair_scores = {x[1]: x[0] for x in results}, []
    for ogpair in null_dict:
        pair_scores.append(omes2dist[null_dict[ogpair]])
    for i in range(oldLen - len(null_dict)):
        pair_scores.append(0)

    return omes2dist, sorted(pair_scores)


def formCooccurStructures(pairs, min_omes, cc_arr_path = None):
    '''
    Imports the out_dicts from parseLoci and creates index vectors that bear
    the ome_num's with a given ogpair cooccurence. e.g. cooccur_dict[(og1, og2)] = [1, 2, 3]
    '''

    cooccur_dict = defaultdict(list)
    for ome in pairs:
        for key in pairs[ome]:
            new_key = tuple(sorted(key))
            cooccur_dict[new_key].append(ome)

    cooccur_dict = {
        x: tuple(sorted(cooccur_dict[x])) for x in cooccur_dict \
        if len(cooccur_dict[x]) > min_omes
        }

    ogpairs, revogpairs, cooccur_dict, cooccur_array = FormCooccurArray(
        cooccur_dict, len(pairs)
        )

    return cooccur_array, dict(cooccur_dict), ogpairs, revogpairs


def load_seedScores(file_, seed_thresh):#, seed_thresh):

    out_ogs = []
    with gzip.open(file_, 'rt') as raw:
        for line in raw:
            if not line.startswith('#'):
                data = [x.rstrip() for x in line.split('\t')]
                if float(data[3]) > seed_thresh:
                    out_ogs.append(line.split('\t'))

    return out_ogs


def form_ogPairDict(out_ogs):

    ogs = set([int(x[0]) for x in out_ogs]) # make a set of the first og in an
    # og-pair
    ogs = ogs.union(set([int(x[1]) for x in out_ogs])) # add in the second og
    ogPair_dict = {x: set() for x in ogs} # prepare the keys to bypass `if`'s
    for i in out_ogs:
        ogPair_dict[int(i[0])].add(int(i[1])) # {og0: set(og1, ..., ogn), } or 
        # {og0: set(og0, og1, ..., ogn)}
        ogPair_dict[int(i[1])].add(int(i[0]))

    return ogPair_dict


def hash_protoclusters(gff_path, ogPair_dict, gene2og, clusplusminus = 10):
    '''parse a gff and compile its organized CDS_dict. identify where og-pairs
    co-occur and retrieve the set of ogs for the locus with the locus seed
    protein'''

    gff_list, protoclus = gff2list(gff_path), defaultdict(list)
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))

    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]): # for each index and protein
            try:
                og0 = gene2og[seq0] 
            except KeyError: # gene not in OG / is a singleton
                continue
            if og0 in ogPair_dict: # if this og is part of a significant seed
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
                            og1 = gene2og[seq1]
                        except KeyError:
                            continue
                        if og1 in ogPair_dict[og0]: # if it is in the set of
                        # ogs for the first queried og
                            og_loc_list = []
                            for gene in locus:
                                try:
                                    og_loc_list.append(gene2og[gene])
                                except KeyError: # singleton or missing
                                    pass
                            og_loc = set(og_loc_list)
                            # the set of ogs in this locus
                            ogPair = tuple(sorted([og0, og1])) # sort the ogs
                            protoclus[ogPair].append([og_loc, seq0, seq1])
                            # {(og0, og1}: [[set(ogx, .. ogy), seq0, seq1], ...]} 
                            
    return list(protoclus)


def merge_protos(protoclus_res):
    '''combine protoclus results from each organism'''

    protoogx2omes = defaultdict(list)
    for ome_protoclus in protoclus_res: # for each protoclus from the mp
    # results
        for ogPair in ome_protoclus: 
            protoogx2omes[ogPair].extend(ome_protoclus[ogPair])
            # {(og0, og1)}: [[set(ogx, ogy), seq], ...]

    return dict(protoogx2omes)


def gen_clusters(loc0, loc1, ogPair_set):

    ome0 = loc0[1][:loc0[1].find('_')] # extract ome from gene acc
    ome1 = loc1[1][:loc1[1].find('_')]
    if ome0 != ome1:
        loc0vloc1 = loc0[0].intersection(loc1[0])
        if loc0vloc1 != ogPair_set: #if there are more overlapping ogs
            ogX = tuple(sorted(list(loc0vloc1))) # sort the loci (not by
            return ogX, ome0, ome1, loc0[1], loc1[1]


def FormulateClusDicts(gen_clus_cmds, ogx2omes, ogx2loc, ome2i, cpus):

    with mp.Pool(processes = cpus) as pool:
        clus_res = pool.starmap(gen_clusters, gen_clus_cmds)
    for i, res in enumerate(clus_res):
        if res:
            for comb_len in range(3, len(res[0]) + 1):
                ogX = res[0]
                if ogX not in ogx2omes:
                    ogx2omes[ogX] = set()
                    ogx2loc[ogX] = set()
                ogx2omes[ogX] = ogx2omes[ogX].union({
                     ome2i[res[1]], ome2i[res[2]]
                     })
                ogx2loc[ogX] = ogx2loc[ogX].union({res[3], res[4]})

    return ogx2omes, ogx2loc


def par_rm(ogX, ogx2omes):

    t_comb_sets, t_combs, todel = [set(ogX)], [ogX], []
    for comb_len in reversed(range(3, len(ogX))):
        for comb in combinations(ogX, comb_len):
            set_comb = set(comb)
            for i, set_ogX in enumerate(t_comb_sets):
                if set_comb < set_ogX:
                    try:
                        if ogx2omes[comb] == ogx2omes[t_combs[i]]:
                            todel.append(comb)
                            break
                    except KeyError:
                        break
            t_combs.append(set_comb)

    return todel


def rm_subsets(ogx2omes, ogx2loc, cpus = 1):

    ogx2omes = {
        k: ogx2omes[k] for k in sorted(
            list(ogx2omes.keys()), 
            key = lambda x: len(x), reverse = True
            )
        }
    max_len = len(list(ogx2omes.keys())[0])
    for ogX_len in reversed(range(4, max_len + 1)):
        i2ogx = list(ogx2omes.keys())
        rm_cmds = [[ogX, ogx2omes] for ogX in i2ogx if len(ogX) == ogX_len]
        with mp.Pool(processes = cpus) as pool:
            todels = pool.starmap(par_rm, rm_cmds)
        for todel in todels:
            for ogX in todel:
                try:
                    del ogx2omes[ogX]
                    del ogx2loc[ogX]
                except KeyError:
                    continue

    return ogx2omes, ogx2loc
                            

def ogPair2ogX(db, ogPair_dict, gene2og, ome2i, cpus, clusplusminus = 10):

    hash_protoclus_cmds = []
    gffs = [x for x in list(db['gff3']) if x.replace('.gff3','') in ome2i]
    print('\t\tForming protocluster hashes', flush = True)
    for gff in gffs: # prepare for protocluster hashing by organism
        hash_protoclus_cmds.append((formatPath(gff), ogPair_dict, gene2og, clusplusminus))
    with mp.Pool(processes = cpus) as pool:
        protoclus_res = pool.starmap(hash_protoclusters, hash_protoclus_cmds)
    pool.join()

    print('\t\tForming proto-OGxs', flush = True)
    protoogx2omes = merge_protos(protoclus_res) # merge mp results

    print('\t\tGenerating high order OGxs', flush = True)
    print('\t\t\t' + str(sum([
           len(protoogx2omes[x])**2 - len(protoogx2omes[x]) \
           for x in protoogx2omes
        ])) + ' to run', flush = True)
    count = 0

    ogx2omes, ogx2loc, gen_clus_cmds = {}, {}, []
    for ogPair in protoogx2omes:
        ogPair_set = set(ogPair)
        for loci in combinations(protoogx2omes[ogPair], 2):
            loc0, loc1 = loci[0], loci[1]
            gen_clus_cmds.append([loc0, loc1, ogPair_set])
        if len(gen_clus_cmds) > 5000000:
            count += len(gen_clus_cmds)
            print('\t\t\t' + str(count), flush = True)
            ogx2omes, ogx2loc = FormulateClusDicts(
                gen_clus_cmds, ogx2omes, ogx2loc, ome2i, cpus
                )
            gen_clus_cmds = []

    ogx2omes, ogx2loc = FormulateClusDicts(
        gen_clus_cmds, ogx2omes, ogx2loc, ome2i, cpus
        )

    print('\t\tRemoving subset OGxs', flush = True)
    ogx2omes, ogx2loc = rm_subsets(ogx2omes, ogx2loc, cpus)
    ogx2omes = {x: tuple(sorted(list(ogx2omes[x]))) for x in ogx2omes}

    return ogx2omes, ogx2loc


def Hash4nulls(
    gff_path, gene2og, max_size, plusminus
    ):

    gff_list = gff2list(gff_path) # open here for decreased serialization
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))
    sml_sizes = {i: [] for i in range(3, (plusminus * 2) + 2)} # for OGx sizes
    if max_size > (plusminus*2) + 1: # if biggest OGx is beyond the window size
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
            loc_og = set([gene2og[x] for x in locus if x in gene2og]) # translate to OGs
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


def genOGxNulls(
    gffs, gene2og, max_clus_size, plusminus, phylo,
    ogx_perc, clus_perc, wrk_dir,
    omes2dist = {}, samples = 10000, cpus = 1
    ): # NEED to adjust; this currently uses way too much memory
    '''Generates a null distribution of randomly sampled OGxs for each 
    # of orthogroups observed in OGxs. Applies the percentile for each
    size as the minimum value for significance.
    Outputs a dictionary of the minimum values for each size OGx
    based on the inputted percentiles. {# of OGs: minimum value}'''

    print('\t\tParsing for random samples', flush = True)
    hash_null_cmds = [
        (x, gene2og, max_clus_size, plusminus,) \
        for x in gffs
        ]
    with mp.Pool(processes = cpus) as pool:
        hashRes = pool.starmap(Hash4nulls, hash_null_cmds)
        # hashRes = [({size: [OGx]})...] by ome

    print('\t\tCalculating microsyteny distances of random samples\n\t\tOGx size:', flush = True)
    ogxBordPercs, ogxClusPercs = {}, {}
    for size in range(3, max_clus_size + 1): # for all observed sizes
        print('\t\t\t' + str(size), flush = True)
        size_dict = {i: v[size] for i, v in enumerate(hashRes)}
        # spoof i keys, size_dict values
        ogx2i, i2ogx, null_dict = gen_null_dict(size_dict, samples)
        oldLen = len(null_dict)
        null_dict = {k: v for k, v in null_dict.items() if len(v) > 1} # 0 will n t compute
        distRes = calcDists(phylo, null_dict, omes2dist = omes2dist, cpus = cpus)
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

        # sort to apply percentile threshold for each size OGx and for both the
        # border percentile and cluster percentile
        nullSizes.sort()
        ogxBordPercs[size] = nullSizes[round(ogx_perc * len(nullSizes) + .5)]
        ogxClusPercs[size] = nullSizes[round(clus_perc * len(nullSizes) + .5)]

    return ogxBordPercs, ogxClusPercs


def loadNulls(max_clus_size, wrk_dir, ogx_perc, clus_perc):

    ogxBordPercs, ogxClusPercs = {}, {}
    for size in range(3, max_clus_size + 1):
        nullSizes = []
        with open(wrk_dir + str(size) + '.null.txt', 'r') as raw:
            for line in raw:
                nullSizes.append(float(line.rstrip()))
        nullSizes.sort()
        ogxBordPercs[size] = nullSizes[round(ogx_perc * len(nullSizes) + .5)]
        ogxClusPercs[size] = nullSizes[round(clus_perc * len(nullSizes) + .5)]

    return ogxBordPercs, ogxClusPercs


def hash_ogx(gff_path, ome, ogX_genes, gene2og, clusplusminus):
    '''sliding window of size clusplusminus, identify sequences that may have
    ogs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus'''

    gff_list, gene_dict = gff2list(gff_path), {}
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))

    ogx_dict = {}
    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]):
            if seq0 in ogX_genes: # if the og is part of a significant seed
            # locus
                if i0 < clusplusminus: # is i0 - clusplusminus < 0 ?
                    locus = cds_dict[scaf][:i0+clusplusminus+1] # then gather
                    # all the beginning
                else:
                    locus = cds_dict[scaf][i0-clusplusminus:i0+clusplusminus+1]
                    # instead get the +/- and adjust for python
                og0 = gene2og[seq0]
                for ogX in ogX_genes[seq0]:
                    if ogX not in ogx_dict:
                        ogx_dict[ogX] = {og: [] for og in ogX}
                    ogx_dict[ogX][og0].append(seq0) # add the sequence to the ogX_dict
                    start, end = None, None
                    for i1, seq1 in enumerate(locus): # for each index and sequence
                    # in the locus
                        try:
                            og1 = gene2og[seq1]
                        except KeyError: # missing entry
                            continue
                        if og1 in set(ogX) and start is None: # if the og is
                        # in the sig clus and we haven't started
                            start = i1 # this is the start index
                            ogx_dict[ogX][og1].append(seq1)

                        elif og1 in set(ogX): # otherwise if it is in the
                        # sig clus label this the other border unless another is 
                        # found
                            end = i1 + 1
                            ogx_dict[ogX][og1].append(seq1)
                    for gene in locus[start:end]:
                        if gene not in gene_dict:
                            gene_dict[gene] = []
                        gene_dict[gene].append(ogX)

    gene_tup = tuple([
        tuple([gene, tuple(sorted(vals))]) \
        for gene, vals in gene_dict.items()
        ])
    # ((gene, {gene: [OGx...]}))

    ogx_tup = tuple([
        tuple([
            ogx, 
            tuple([tuple([og, tuple(sorted(set(seqs)))]) \
                for og, seqs in ogs.items()]) \
            ])
            for ogx, ogs in ogx_dict.items()
        ])

    return ome, gene_tup, ogx_tup

def findOGxPairSingle(gene2ogx, ome, ogx2i, pairsDict, minimum = 2):
    omePairs_dict = defaultdict(list)
    for gene, ogxs in gene2ogx.items():
        for pairTup in combinations(ogxs, 2):
            omePairs_dict[pairTup].append(gene) # gene centric 
#            try: #og approach
#                omePairs_dict[pairTup].append(gene2og[gene]) # og approach
 #           except KeyError: # gene not in an og/is singleton
   #             pass

    for pairTup, genes in omePairs_dict.items():
        if len(genes) >= minimum:
            ogXpair = tuple(sorted([ogx2i[pairTup[0]], ogx2i[pairTup[1]]]))
            pairsDict[ogXpair].append(ome)
    return pairsDict


def findOGxPairs(gene2ogx, ome, count = 2):
    omePairs_dict = defaultdict(list)
    for gene in gene2ogx:
        for combo in combinations(gene2ogx[gene], 2):
            try:
                omePairs_dict[combo].append(gene2og[gene])
            except KeyError:
                pass

    omePairs = set([x for x,v in omePairs_dict.items() if len(set(v)) >= count])
    # if there are two, independent overlapping OGs between the OGxs, consider
    # them a pair
    return ome, omePairs

def OGx2locSimilarity(ogx, ogxGenes, ome2i, ogx_dir):
    omes = []
    for og, genes in ogxGenes.items():
        for gene in genes:
            omes.append(ome2i[gene[:gene.find('_')]])
    omes = sorted(set(omes))
    ome2ome = {i: {og: [] for og in ogx} for i in combinations(omes, 2)}
    for og in ogx:
        with open(ogx_dir + str(og) + '.out', 'r') as raw:
            # unfortunately, I envisage a scenario where Diamond doesn't
            # compare two genes, even though they are part of related clusters
            # due to its heuristics; thus it will show up as a zero value later
            # on and vastly deplete the average overall. will just have to see
            # how this shakes out with known clusters
            genes = set(ogxGenes[og])
            for line in raw:
                q,s,val = line.split('\t')[:3]
                if {q,s}.issubset(genes):
                    ome0, ome1 = q[:q.find("_")], s[:s.find('_')]
                    if ome0 != ome1:
                        i = tuple(sorted([ome2i[ome0], ome2i[ome1]]))
                        ome2ome[i][og].append(float(val)/100)

    ome2ome_adj = {}
    for i, ogs in ome2ome.items():
        ome2ome_adj[i] = 0
        for og, vals in ogs.items():
            try:
                ome2ome_adj[i] += max(vals) # take the highest
            except ValueError: # no hits
                continue
        ome2ome_adj[i] /= len(ogx) # if an og is missing, its a zero
        # contribution to the average.
        #  does this need to be a distance for mcl?

    return tuple([tuple([i, v]) for i,v in ome2ome_adj.items()])


def OGxClan2cluster_families(
    adj_path, ome2ogx_path, clan, clanNum, ogxGenes, ome2i, i2ogx, ogx2i, wrk_dir,
    inflation = 1.5
    ):
    '''ogxGenes = {ogx: og: []}'''
    ogx_dir = wrk_dir + 'ogx/'


    node_adj, cmds, ome2ogx, sims = {}, [], {}, {}
    for ogxI in clan:
        ogx = i2ogx[ogxI]
        data = OGx2locSimilarity(
            ogx, ogxGenes[ogx], ome2i, ogx_dir
            )
        for i, av_sim in data:
            if i not in sims:
                sims[i] = []
            sims[i].append(av_sim)
        for og, genes in ogxGenes[ogx].items():
            for gene in genes:
                ome = gene[:gene.find('_')]
                iOme = ome2i[ome]
                if iOme not in ome2ogx:
                    ome2ogx[iOme] = []
                ome2ogx[iOme].append(ogx)

    i2i_adj = {i: max(sim) for i, sim in sims.items()}
    ome2ogx = {i: sorted(set(ogxs)) for i, ogxs in ome2ogx.items()}

#    i2i_adj = {i: max(vals) for i, vals in i2iSim.items()} 
    # take the highest similarity value as the representative for each
    # ome 2 ome comparison. ideally, this should be locus by locus
    # and a future iteration should implement a robust hash to make
    # that process more efficient

    # normalize - I dont think I should do this, perhaps MCL should
    # take care of this, and perhaps a large graph is needed to make
    # sure MCL normalizes with respect to everything else
#    minv_, maxv_ = min(i2i_adj.values()), max(i2i_adj.values())
 #   denom = maxv_ - minv_
  #  try: # should this be a distance for mcl?
   #     i2i_adj = {i: (v-minv_)/denom for i, v in i2i_adj.items()}
#    except ZeroDivisionError: # they all have the same
 #       i2i_adj = {i: 1 for i in i2i_adj.items()}
        # should just return as one cluster here

    with open(adj_path + '.tmp', 'w') as out:
        for omePair, val in i2i_adj.items():
            out.write(str(omePair[0]) + ' ' + str(omePair[1]) + ' ' + str(val) + '\n')
    os.rename(adj_path + '.tmp', adj_path)

    with open(ome2ogx_path + '.tmp', 'wb') as out:
        pickle.dump(ome2ogx, out)
    os.rename(ome2ogx_path + '.tmp', ome2ogx_path)

    return ome2ogx



def MCL(adj_path, clusFile, inflation = 1.5, cpus = 1):

    if cpus > 5:
        cpus = 5

    subprocess.call([
        'mcl', adj_path, '--abc', '-I', str(inflation),
        '-o', clusFile
        ], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
        )
#    mclFile = prep_mcl(adj_path + '.tmp', sym = True, output = None) 
    # symmetrical matrix, no verbose, mcxload must be in path
 #   run_mcl(mclFile, clusFile, inflation, cpus = cpus)
    # parallelization could be more efficient overall

    return clusFile


def parseMCL(clusFile, ome2ogx):
    clusters = []
    with open(clusFile, 'r') as raw:
        for line in raw:
            clusters.append([int(x) for x in line.rstrip().split()])
#            i, clus = [int(x) for x in line.rstrip().split()]
 #           if clus not in clus2i:
  #              clus2i[clus] = []
   #         clus2i[clus].append(i)

    moduleOGxs = []
    moduleOmes = []
    modules = []

    for iOmes in clusters:
        moduleOGxs.append([])
        moduleOmes.append(sorted(iOmes))
        modules.append({})
        for i in iOmes:
           for ogx in ome2ogx[i]:
               moduleOGxs[-1].extend(list(ogx))
               if ogx not in modules[-1]:
                   modules[-1][ogx] = []
               modules[-1][ogx].append(i)
        moduleOGxs[-1] = sorted(set(moduleOGxs[-1]))
    return modules, moduleOmes, moduleOGxs

def FamMNGR(
    args
#    clan, i, ogx2genes, ome2i, i2ogx, wrk_dir, inflation
    ):
    clan, i, ogx2genesT, ome2iT, i2ogxT, wrk_dir, inflation = args
    ogx2genes = {ogx: dict(vals) for ogx, vals in ogx2genesT}
    ome2i, i2ogx = dict(ome2iT), dict(i2ogxT)


    ogx2clus_dir = wrk_dir + 'ogx2clus/'
    ome2ogx = None # prime for loading if necessary
    adj_path = ogx2clus_dir + str(i) + '.mtx'
    clus_path = ogx2clus_dir + str(i) + '.clus'
    ome2ogx_path = ogx2clus_dir + str(i) + '.pickle'
    print('\t\t\tClan: ' + str(i), flush = True)
    ogx2i = {v: k for k, v in i2ogx.items()} # make here to decrease
    # serialization

    if not os.path.isfile(adj_path):
        if not os.path.isfile(clus_path):
            ome2ogx = OGxClan2cluster_families(
                adj_path, ome2ogx_path,
                clan, i, ogx2genes, ome2i, i2ogx, ogx2i, wrk_dir, inflation
                )
    if not os.path.isfile(clus_path):
        MCL(adj_path, clus_path, inflation = 1.5)
    if not ome2ogx:
        with open(ome2ogx_path, 'rb') as raw:
            ome2ogx = pickle.load(raw)
    tm, tmo, tmx = parseMCL(clus_path, ome2ogx)
    tuples = []
    for module in tm:
        tuples.append(tuple([tuple(module.items())]))
    tm = tuple(tuples)
    return tm, tmo, tmx


def WriteAdjMatrix(Q, out_file):
    with open(out_file, 'w') as out: # might be nice to compress this/write to MCL binary matrix
        x = True
        while x:
            x = Q.get()
            if x:
                out.write(x)
                out.flush() # shouldn't be a bottleneck
    
def BLASTclanOG(db, og, fileBase, genes, minid, diamond = 'diamond'):
    blast_ids = defaultdict(dict)
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
        os.rename(fileBase + '.out.tmp', fileBase + '.out')

    with open(fileBase + '.out', 'r') as raw:
        for line in raw:
            q, s, pident = line.split()
            if q not in blast_ids[og]:
                blast_ids[og][q] = {}
            blast_ids[og][q][s] = float(pident)/100 # adjust diamond to decimal

    return dict(blast_ids)

def AcquireClusFamSim(
    i0, i1, index, loc0, loc1, ogL0, ogL1, blast_ids, Q
    ):
    set0, set1 = set(ogL0), set(ogL1)
    if None in set0:
        set0.remove(None)
    if None in set1:
        set1.remove(None)
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
    # primary OGx into separate families

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



def mpClan2FamLoci(
    db, clanI, loci, ogLoci, ogxXloci, fam_dir, Q, index, 
    diamond = 'diamond', minid = 30, cpus = 1
    ):

    blast_hash = defaultdict(list)
    for i, locus in enumerate(loci):
        ogs = ogLoci[i]
        for i1, og in enumerate(ogs):
            if og is not None:
                blast_hash[og].append(locus[i1])
    
    blast_ids = defaultdict(dict)
    cmds = []
    for og, genes in blast_hash.items():
        if len(genes) > 1:
            fileBase = fam_dir + str(clanI) + '.' + str(og)
            cmds.append([
                db, og, fileBase, genes, minid, diamond
                ])
    with mp.Pool(processes = cpus) as pool:
        tblast_res = pool.starmap(BLASTclanOG, cmds)

    for tblast_ids in tblast_res:
        for og, qdict in tblast_ids.items():
            for q, sdict in qdict.items():
                if q in blast_ids[og]:
                    blast_ids[og][q] = {**blast_ids[og][q], **sdict}
                else:
                    blast_ids[og][q] = sdict

    blast_ids = dict(blast_ids)    
    if blast_ids:
        with mp.Pool(processes = cpus) as pool:
            pool.starmap(
                AcquireClusFamSim,
                ([i0, i1, index, loci[i0], loci[i1], ogLoci[i0], ogLoci[i1], blast_ids, Q] \
                for i0, i1 in combinations(range(len(loci)), 2))
                )

def Clan2FamLoci(
    db, clanI, loci, ogLoci, ogxXloci, fam_dir, Q, index, 
    diamond = 'diamond', minid = 30
    ):
    print(clanI, flush = True)
    blast_hash = defaultdict(list)
    for i, locus in enumerate(loci):
        ogs = ogLoci[i]
        for i1, og in enumerate(ogs):
            if og is not None:
                blast_hash[og].append(locus[i1])
    
    blast_ids = defaultdict(dict)
    for og, genes in blast_hash.items():
        if len(genes) > 1:
            fileBase = fam_dir + str(clanI) + '.' + str(og)
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
                os.rename(fileBase + '.out.tmp', fileBase + '.out')
        
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
            AcquireClusFamSim(i0, i1, index, loc0, loc1, ogL0, ogL1, blast_ids, Q)


def HashClanLoci(ome, gff_path, ome_sig_clus, gene2og, clusplusminus):
    '''sliding window of size clusplusminus, identify sequences that may have
    ogs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus'''

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
                            og = gene2og[seq1]
                        except KeyError:
                            continue
                        if og in sig_clus[0] and start is None: # if the og is
                        # in the sig clus and we haven't started
                            start = i1 # this is the start index
                        elif og in sig_clus[0]: # otherwise if it is in the
                        # sig clus label this the other border unless another is 
                        # found
                            end = i1 + 1
                    ogX = tuple(sorted(list(sig_clus[0])))
                    preclanLoci[clan].append([
                         locus[start:end], ogX
                         ])
                    # clus_out = [(og0, og1, ... ogn), [prot0, prot1, ... protn])]

    outclanLoci, outclanOGx = defaultdict(list), defaultdict(list)
    for clan, loci in preclanLoci.items():
        while loci: # exhaustively merge overlapping loci
            loc0, ogxs0, locIntersect = set(loci[0][0]), loci[0][1], None
            for i1, loc1d in enumerate(loci[1:]):
                loc1, ogxs1 = set(loc1d[0]), loc1d[1]
                locIntersect = loc0.intersection(loc1)
                if locIntersect: # check for overlap
                    newOGx = list(ogxs0)
                    newOGx.extend(list(ogxs1))
                    loci.append([list(loc1.union(loc0)), newOGx])
                    break
            if locIntersect: # if there was overlap, delete the overlappers
                del loci[0]
                del loci[i1 + 1]
            else: # no more overlap for this locus, add to the final output
                outclanLoci[clan].append(sorted(loc0))
                outclanOGx[clan].append(tuple(sorted(ogxs0)))
                del loci[0]

# blast each locus OG against itself
    outogLoci = {}
    for clan, loci in outclanLoci.items():
        outogLoci[clan] = []
        for locus in loci:
            outogLoci[clan].append([])
            for gene in locus:
                try:
                    outogLoci[clan][-1].append(gene2og[gene])
                except KeyError:
                    outogLoci[clan][-1].append(None)
    
    return ome, outclanLoci, outclanOGx, outogLoci




def groupOGx(
    ogx2loc, db, gene2og, i2ogx, ogx2i,
    phylo, clusScores, bordScores, ome2i, ogx2omes,
    wrk_dir, omes2dist = {}, clusplusminus = 3, 
    inflation = 1.5, minimum = 2, cpus = 1
    ):

    groupI = wrk_dir + 'group.I.pickle'
    groupII = wrk_dir + 'group.II.pickle'
    if not os.path.isfile(groupI):
        print('\tDiscovering OGxs with overlapping loci', flush = True)

        ogX_genes = {}
        for ogX in ogx2loc:
            if ogX in ogx2i: # if it passed the border threshold
                for gene in ogx2loc[ogX]:
                    ome = gene[:gene.find('_')]
                    if ome not in ogX_genes:
                        ogX_genes[ome] = {}
                    if gene not in ogX_genes[ome]:
                        ogX_genes[ome][gene] = []
                    ogX_genes[ome][gene].append(ogX)
    
        hashOgx_cmds = []
        start = datetime.now()
        for ome in ogX_genes:
            hashOgx_cmds.append([
                db[ome]['gff3'],
                ome, ogX_genes[ome], gene2og, clusplusminus
                ])
        with mp.Pool(processes = cpus) as pool:
            gene_tups = pool.starmap(hash_ogx, hashOgx_cmds)
    
        gene2ogx, ogx2genes = {}, {}
        for res in gene_tups:
            gene2ogx[ome2i[res[0]]] = {x[0]: set(x[1]) for x in res[1]}
            for ogx, ogs in res[2]:
                if ogx not in ogx2genes:
                    ogx2genes[ogx] = {og: [] for og in ogx}
                for og, seqs in ogs:
                    ogx2genes[ogx][og].extend(seqs)
            # ogx2genes = {ogx: og: [gene, gene]}
            # gene2ogx = {omeI: gene: {OGx}}
        with open(groupI, 'wb') as out:
            pickle.dump([gene2ogx, ogx2genes], out)
        print('\t\t\t' + str(datetime.now() - start))
    elif not os.path.isfile(groupII):
        print('\tLoading OGxs with overlapping loci', flush = True)
        with open(groupI, 'rb') as in_:
            gene2ogx, ogx2genes = pickle.load(in_)

    groupIII = wrk_dir + 'group.III.pickle'
    if not os.path.isfile(groupIII):
        print('\tIdentifying significantly overlapping OGxs', flush = True)
        tpairDict = defaultdict(list)
        for ome, gene2ogx_ome in gene2ogx.items():
            tpairDict = findOGxPairSingle(
                gene2ogx_ome, ome, ogx2i, tpairDict,
                minimum = minimum
                )

        pairDict = {}
        for id_, omes in tpairDict.items():
            omes_set = set(omes)
            if len(omes_set) > 1: # need at least one ome to calc a branch length
                pairDict[id_] = tuple(sorted(omes_set)) 


        omes2dist = updateScores(phylo, pairDict, cpus = cpus, omes2dist = omes2dist)

        print('\tClassifying OGx clans', flush = True)
        # populate a lil_matrix here, then use that network to identify modules
        print('\t\tBuilding binary OGx-OGx network', flush = True)
        matrix = lil_matrix((len(i2ogx), len(i2ogx)), dtype=bool)
        for idPair, omes in pairDict.items():
            i0, i1 = idPair[0], idPair[1]
            nullSize = max([len(i2ogx[x]) for x in idPair])
            if omes2dist[omes] >= bordScores[nullSize]:
                matrix[i0, i1] = True
                matrix[i1, i0] = True
    
        print('\t\tIsolating subgraphs', flush = True)
        network = nx.from_scipy_sparse_matrix(matrix)
        subgraphs = [network.subgraph(c) for c in nx.connected_components(network)]
        # use .copy() after iterator to copy
        clans, clanOmes, clanOGxs = [], [], []
        for sg in subgraphs:
            clans.append({})
            clanOmes.append([])
            clanOGxs.append([])
            for id_ in list(sg.nodes):
                ogX = i2ogx[id_]
                omes = ogx2omes[ogX]
                clans[-1][ogX] = omes # storing in a hash to allow future
                # manipulation
                clanOGxs[-1].extend(ogX)
                clanOmes[-1].extend(omes)
            clanOmes[-1] = tuple(sorted(set(clanOmes[-1])))
            clanOGxs[-1] = tuple(sorted(set(clanOGxs[-1])))

        omes2dist = updateScores(
            phylo, {clanOGxs[i]: omes for i, omes in enumerate(clanOmes)},
            cpus = cpus, omes2dist = omes2dist
            )
        with open(groupIII, 'wb') as out:
            pickle.dump([clans, clanOmes, clanOGxs], out)
    else:
        print('\tLoading OGx clans', flush = True)
        with open(groupIII, 'rb') as in_:
            clans, clanOmes, clanOGxs = pickle.load(in_)
   
    print('\t\t' + str(len(clans)) + ' clans', flush = True)

    print('\tClassifying OGx clans into cluster families', flush = True)
    clanFile = wrk_dir + 'clan2loci.'
    if not os.path.isfile(clanFile + 'og.json.gz'):
        print('\t\tPreparing loci extraction', flush = True)
        ome2clus2extract = defaultdict(dict)
        for i, clan in enumerate(clans):
            for ogx, omes in clan.items():
                for gene in ogx2loc[ogx]:
                    ome = gene[:gene.find('_')]
                    if gene not in ome2clus2extract[ome]:
                        ome2clus2extract[ome][gene] = [[set(ogx), i]]
                    else:
                        ome2clus2extract[ome][gene].append([set(ogx), i])
                    # {ome: {gene: [[set(ogx), clanI]]}}

        print('\t\tExtracting clan loci', flush = True)
        cmds = []
        for ome, clus2extract in ome2clus2extract.items():
            # prepare commands for locus extraction
            cmds.append([
                ome, db[ome]['gff3'],
                clus2extract, gene2og, clusplusminus 
                ])
        with mp.Pool(processes = cpus) as pool:
            hash_res = pool.starmap(HashClanLoci, cmds)
    
        clanLoci, clanOGx4fams, clanOGloci = defaultdict(list), defaultdict(list), defaultdict(list)
        for ome, outclanLoci, outclanOGx, outclanOGloci in hash_res:
            for clan in outclanLoci:
                clanLoci[clan].extend(outclanLoci[clan])
                clanOGx4fams[clan].extend(outclanOGx[clan])
                clanOGloci[clan].extend(outclanOGloci[clan])
        writeJson(clanFile + 'json.gz', clanLoci)
        writeJson(clanFile + 'ogx.json.gz', clanOGx4fams)
        writeJson(clanFile + 'og.json.gz', clanOGloci)
    else:
        print('\t\tLoading clan loci', flush = True)
        clanLoci = {
            int(k): tuple(v) for k,v in sorted(
                readJson(clanFile + 'json.gz').items(), key = lambda x: len(x[1]), reverse = True
            )}
        clanOGx4fams = {int(k): tuple([tuple(i) for i in v]) for k,v in readJson(clanFile + 'ogx.json.gz').items()}
        clanOGloci = {int(k): tuple(v) for k, v in readJson(clanFile + 'og.json.gz').items()}

    fam_dir = wrk_dir + 'fam/'
    if not os.path.isdir(fam_dir):
        os.mkdir(fam_dir)
   
    print('\t\tCalling cluster families', flush = True)
    print('\t\t\t' + str(sum([len(v) for v in list(clanLoci.values())])) + ' loci', flush = True)
    m, adj_mtr = mp.Manager(), fam_dir + 'loci.adj.tmp'

    Q = m.Queue()
    W = mp.Process(target=WriteAdjMatrix, args=(Q, adj_mtr))
    W.start()

    if cpus < 2:
        cpus = 2
    index, cmds = 0, []

    loci, ogxXloci = {}, {}
    for clanI in clanLoci:
        cmds.append([
            db, clanI, clanLoci[clanI], clanOGloci[clanI],
            clanOGx4fams[clanI], fam_dir, Q, index
            ])
        for i, locus in enumerate(clanLoci[clanI]):
            loci[index] = clanLoci[clanI][i]
            ogxXloci[index] = clanOGx4fams[clanI][i]
            index += 1

    del cmds[0] # remove the first one because it is huge enough to mp on its own
    with mp.Pool(processes = cpus - 1) as pool:
        pool.starmap(
            Clan2FamLoci, cmds
            )
    mpClan2FamLoci(
        db, 0, clanLoci[0], clanOGloci[0], 
        clanOGx4fams[0], fam_dir, Q, 0,
        diamond = 'diamond', minid = 30, cpus = cpus - 1
        )
    W.join()
    m.close()

    os.rename(adj_mtr, fam_dir + 'loci.adj')
    MCL(
        fam_dir + 'loci.adj', fam_dir + 'loci.clus', 
        inflation = 1.5, cpus = cpus*2
        )
    fams = []
    with open(fam_dir + 'loci.clus', 'r') as raw:
        for line in raw:
            d = line.rstrip().split('\t')
            if len(d) > 1:
                indices = [int(x) for x in d]
            else: # singletons won't pass thresholds, disregard them
                continue
            fams.append(indices)

    fams, famOGxs, famOmes = [], [], []
    for fam, locIs in enumerate(fams):
        fams.append(defaultdict(list))
        famOGx, famOme_list = [], []
        for locI in locIs:
            loc = loci[locI]
            ogxs = ogxXloci[locI]
            omeI = ome2i[loc[0][:loc[0].find('_')]]
            [fams[-1][ogx].append(omeI) for ogx in ogxs]
            famOme_list.append(omeI)
            famOGx.extend(ogxs)
        famOGxs.append(tuple(sorted(set(famOGx))))
        famOmes.append(sorted(set(famOme_list)))

    print('\t\t\t' + str(len(fams)) + ' fams', flush = True)

    omes2dist = updateScores(
        phylo, {famOGxs[i]: omes for i, omes in enumerate(famOmes)},
        cpus = cpus, omes2dist = omes2dist
        )

    return fams, famOGxs, famOmes, omes2dist


def exSigClus(
    clus_scores_dict, ogx2loc, top_ogxs, ome2i
    ):
    '''Create a hash to seed retrieving clusters based on
    their entry in ogx2loc.'''
 
    sig_clus = defaultdict(dict)
    for top in top_ogxs:
        ogX, clan, omeIs = top[0], top[1], top[2]
        for gene in ogx2loc[ogX]:
            ome = gene[:gene.find('_')]
            if ome2i[ome] in omeIs:
                if gene in sig_clus[ome]:            
                    sig_clus[ome][gene] = [[set(ogX), clan]]
                else:
                    sig_clus[ome][gene].append([set(ogX), clan])

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
    '''recursively add branches without the trait to the missing total'''

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

def calcPatchiness(phylo, omes):
    '''calculate the percent branch length that a trait is missing over the
    MRCA of where the trait is present'''
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


def makeBlastDB(og_file, out_file, makeblastdb):
#    og_int = int(re.sub(r'\.fa$', '', os.path.basename(og_file).replace('OG','')))
 #   out_file = blast_dir + str(og_int) + '.db'
    cmd = [
        makeblastdb, '-in', og_file, '-dbtype', 'prot',
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


def parseRes(res_file, ome_set, ogx, og):

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

    return results, ogx, og


def calcRBHcorr(ogx_res, ogx):

    omes, total = set(), 0
    for og in ogx_res:
        t_dict = {}
        for q in ogx_res[og]:
            ome = q[:q.find('_')]
            if ome not in t_dict:
                t_dict[ome] = 0
            if ogx_res[og][q] == 1: # if there are any hits for an ome w/multiple queries
                t_dict[ome] = 1
        try:
            total += sum(t_dict.values()) / len(t_dict)
        except ZeroDivisionError:
            print('\t\t' + str(og) + ' no hits', flush = True)
    try:
        total /= len(ogx_res)
    except ZeroDivisionError:
        print('\t\t' + str(ogx) + ' no hits', flush = True)

    return total, ogx


def runMakeDMNDdb(db, diamond, ogx_dir, og, genes):
    fa_dict = acc2fa(db, genes)
    makeDBcmd = subprocess.Popen([
        diamond, 'makedb', '--in', '-', '--db',
        ogx_dir + str(og) + '.dmnd'
        ], stdin = subprocess.PIPE, stdout = subprocess.DEVNULL,
        stderr = subprocess.DEVNULL
        )
    makeDBcmd.communicate(input=dict2fa(fa_dict).encode())[0]
    makeDBcmd.stdin.close()
    makeDBcmd.wait()

def ogBlast(db, ogs, og_dir, ogx_dir, diamond, og2gene, cpus = 1, printexit = False):


    cmds1, cmds2 = [], []
    for og in ogs:
        if not os.path.isfile(ogx_dir + str(og) + '.dmnd'):
            cmds1.append([
                db, diamond, ogx_dir, og, og2gene[og]
                ])
        if not os.path.isfile(ogx_dir + str(og) + '.out'):
            cmds2.append([
                diamond, 'blastp', '--query', og_file, 
                '--db', ogx_dir + str(og) + '.dmnd', 
                '-o', ogx_dir + str(og) + '.out'
                ])

    print('\tMaking ' + str(len(cmds1)) + ' OG diamond databases', flush = True)
    with mp.Pool(processes = cpus) as pool:
        pool.starmap(runMakeDMNDdb, cmds1)

    if printexit:
        with open(ogx_dir + '../../rbh.sh', 'w') as out:
            out.write('\n'.join([' '.join(x) for x in cmds2]))
        print('\nRBH commands outputted to `rbh.sh`', flush = True)
        sys.exit(0)

    print('\tRunning ' + str(len(cmds2)) + ' OG self v self diamonds', flush = True)
    multisub(cmds2, processes = cpus)


def retroOGxGeneGrab(
    gff_path, ome_loc, gene2og, clusplusminus
    ):

    gff_list = gff2list(gff_path)
    clus_ogs = defaultdict(list)
    cds_dict, cds_dict2 = compileCDS2(
        gff_list, os.path.basename(gff_path).replace('.gff3','')
        )

    # the goal here is to take each OGx and compile and OG by OG list of each
    # gene that is associated with it - I think I am going to need some sort of
    # while loop to accomodate searches outside of the plusminus range.
    # Essentially this will operate by first checking if the entirety of the
    # OGx is within the plusminus range. The next step is that, if not, go
    # ahead and use a while loop to simultaneously parse up and down beyond the
    # plusminus to grab the rest of the OGx and stop when acquired. This can be
    # accomplished by taking the set of OGx and one by one removing the hits
    # found. The while loop would then operate on what is contained within the
    # set.
    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]):
            if seq0 in ome_loc: # if the og is part of a significant seed
            # locus
                for ogX in ome_loc[seq0]:
                    ogs = set(ogX)
                    clus_ogs[ogX].append(defaultdict(list))
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
                            og = gene2og[seq1]
                        except KeyError:
                            continue
                        if og in ogs:
                            clus_ogs[ogX][-1][og].append(seq1)
                        # I'm envisioning a scenario where a P450 from outside
                        # the cluster joins in, even though a P450 OG was an
                        # original part of the OGx; the newly joining P450
                        # would not hit another organism with the cluster and
                        # thus impact the score. it is an assumption that
                        # this is minor.
    return clus_ogs
    # clus_ogs = {ogx: [{og: []}]} list is per locus

def rbhCalc(
    ogX, ogxDict, ogx_dir
    ):
    # ogxDict = {og: set(gene1, gene2)}

    res= {}
    for og in ogxDict:
        res[og] = {}
        with open(ogx_dir + str(og) + '.out', 'r') as raw:
            geneInfo = defaultdict(list)
            for line in raw:
                d = line.rstrip().split('\t')
                q = d[0]
                s = d[1]
                bit = float(d[-1])
                if q in ogxDict[og]:
                    geneInfo[q].append((s, bit))
        geneInfo = {
            k: sorted(v, key = lambda x: x[1], reverse = True) for k,v in
            geneInfo.items()
            }


        for gene in list(ogxDict[og]):
            ome = gene[:gene.find('_')]
            if ome not in res[og]:
                res[og][ome] = False
            if res[og][ome]:
                continue
            try:
                while geneInfo[gene][0][0][:geneInfo[gene][0][0].find('_')] == ome:
                    del geneInfo[gene][0]
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

    rbhScore = \
        sum([sum(omeScores[ome].values()) / len(omeScores[ome]) for ome in omeScores]) \
        / len(omeScores) # overall average of average binary positive per ome 

    return ogX, rbhScore


def rbhMngr2(
    ogs, omes, og_dir, ogx_dir, diamond, ogx2loc, 
    db, gene2og, clusplusminus, og2gene, cpus = 1
    ):

    ogBlast(db, ogs, og_dir, ogx_dir, diamond, og2gene, cpus = 1, printexit = False)
    ogxGene_cmds = []
    ome_locs = {ome: {} for ome in omes}
    for ogX in ogx2loc:
        for seq in ogx2loc[ogX]:
            ome = seq[:seq.find('_')]
            if seq not in ome_locs[ome]:
                ome_locs[ome][seq] = []
            ome_locs[ome][seq].append(ogX)

    for ome in ome_locs:
        gff = db['gff3'][ome]
        ogxGene_cmds.append([gff, ome_locs[ome], gene2og, clusplusminus])

    print('\tAssimilating loci with significant OGxs', flush = True)
    with mp.Pool(processes = cpus) as pool:
        clusOgs_prep = pool.starmap(retroOGxGeneGrab, ogxGene_cmds)
#    {ogx:[{og:[seq]}]}

    clusOgs = {}
    for res in clusOgs_prep:
        for ogX in res:
            if ogX not in clusOgs:
                clusOgs[ogX] = {x: [] for x in ogX}
            for locus in res[ogX]:
                for og in locus:
                    clusOgs[ogX][og].extend(locus[og])
    clusOgs = {
        ogX: {og: set(v) for og, v in ogs.items()} for ogX, ogs in clusOgs.items()
        } # make sets from it
    # {ogX: {og: set(seqs)}}

    print('\tCalculating RBH scores', flush = True)
    with mp.Pool(processes = cpus) as pool:
        rbhRes = pool.starmap(
            rbhCalc, [[ogX, clusOgs[ogX], ogx_dir] for ogX in clusOgs]
            )

    rbhScores = dict(rbhRes)
     
    return rbhScores        

def rbhMngr3(
    ogs, omes, og_dir, ogx_dir, diamond, ogx2loc, 
    db, gene2og, clusplusminus, og2gene, modules,
    moduleOGxs, ome2i, cpus = 1
    ):

    if not os.path.isfile(ogx_dir + 'clusOGs.pickle'):
        ogxGene_cmds = []
        ome_locs = {ome: defaultdict(list) for ome in omes}
    
        for i, ogx2omes in enumerate(modules):
            modOGx = moduleOGxs[i]
            for ogX, omes in ogx2omes.items():
                omes_set = set(omes)
                for seq in ogx2loc[ogX]:
                    ome = seq[:seq.find('_')]
                    if ome2i[ome] in omes_set:
    #                    if seq not in ome_locs[ome]:
     #                       ome_locs[ome][seq] = []
                        ome_locs[ome][seq].append(modOGx)
            
        for ome in ome_locs:
            gff = db[ome]['gff3']
            ogxGene_cmds.append([gff, ome_locs[ome], gene2og, clusplusminus])
    
        print('\tAssimilating loci with significant OGxs', flush = True)
        with mp.Pool(processes = cpus) as pool:
            clusOgs_prep = pool.starmap(retroOGxGeneGrab, ogxGene_cmds)

        with open(ogx_dir + 'clusOgs.pickle', 'wb') as out:
            pickle.dump(clusOgs_prep, out)
    else:
        with open(ogx_dir + 'clusOgs.pickle', 'rb') as raw:
            clusOGs_prep = pickle.load(raw)

    clusOgs = {modOGx: defaultdict(list) for modOGx in moduleOGxs}
    for res in clusOgs_prep:
        # clus_ogs = {ogx: [{og: []}]} list is per locus
        for ogX in res:
            for locus in res[ogX]:
                for og in locus:
                    clusOgs[ogX][og].extend(locus[og])
    clusOgs = {
        ogX: {og: set(v) for og, v in ogs.items()} for ogX, ogs in clusOgs.items()
        } # make sets from it
    # {ogX: {og: set(seqs)}}

    print('\tCalculating RBH scores', flush = True)
    with mp.Pool(processes = cpus) as pool:
        rbhRes = pool.starmap(
            rbhCalc, [[ogX, clusOgs[ogX], ogx_dir] for ogX in clusOgs]
            )

    rbhScores = dict(rbhRes)
     
    return rbhScores        





def dndsGeneGrab(
    gff_path, assembly_path, proteome_path,
    ome_sig_clus, gene2og, clusplusminus
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
                    og0 = gene2og[seq0]
                except KeyError:
                    continue
                for hit in ome_sig_clus[seq0]:
                    ogX = tuple(sorted(list(hit)))
                    if ogX not in clus_out:
                        clus_out[ogX] = [{og: [{}, {}] for og in ogX}]
                    clus_out[ogX][-1][og0][0][seq0] = prot_dict[seq0]
                    clus_out[ogX][-1][og0][1] = {
                        **clus_out[ogX][-1][og0][1], 
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
                            og = gene2og[seq1]
                        except KeyError:
                            continue
                        if og in hit: # if the og is
                        # in the sig clus
                            clus_out[ogX][-1][og][0][seq1] = prot_dict[seq1]
                            clus_out[ogX][-1][og][1] = {
                                **clus_out[ogX][-1][og][1],
                                **ntmain(cds_dict2[seq1], assem_dict)
                                }
    return clus_out


def dnds_preparation(db, omes4ogs, ogx2omes, gene2og, plusminus, ogx_dir, i2ome):

    gene_checks, ogX_check = {}, {}
    
    for ome in omes4ogs:
        gff_path = db[i2ome[ome]]['gff3'] 
        proteome_path = db[i2ome[ome]]['faa']
        assembly_path = db[i2ome[ome]['fna']
        res = dndsGeneGrab(
            gff_path, assembly_path, proteome_path,
            omes4ogs[ome], gene2og, plusminus
            )
        for ogX in res:
            if ogX not in ogX_check:
                ogX_check[ogX] = set()
            if ogX not in gene_checks:
                gene_checks[ogX] = {og: [{}, {}] for og in ogX}
            for loc in res[ogX]:
                for og in loc:
                    gene_checks[ogX][og][0] = {
                        **gene_checks[ogX][og][0], **loc[og][0]
                        }
                    gene_checks[ogX][og][1] = {
                        **gene_checks[ogX][og][1], **loc[og][1]
                        }
            ogX_check[ogX].add(ome)
            if ogX_check[ogX] == set(ogx2omes[ogX]):
                for og in ogX:
                    out_base = ogx_dir + '-'.join([str(x) for x in ogX]) + \
                        '.' + str(og)
                    with open(out_base + '.aa.fa', 'w') as out:
                        out.write(dict2fa(gene_checks[ogX][og][0]))
                    with open(out_base + '.nt.fa', 'w') as out:
                        out.write(dict2fa(gene_checks[ogX][og][1]))
                del gene_checks[ogX]

    if gene_checks:
        print('\t\tERROR: discrepancies with previous run', flush = True)
        for ogX in gene_checks:
            print('\t\t\t' + str(ogX), flush = True)
            for og in ogX:
                out_base = ogx_dir + '-'.join([str(x) for x in ogX]) + \
                    '.' + str(og)
                with open(out_base + '.aa.fa', 'w') as out:
                    out.write(gene_checks[ogX][og][0])
                with open(out_base + '.nt.fa', 'w') as out:
                    out.write(gene_checks[ogX][og][1])
        

def calc_dnds(mafft):

    og_catch = re.search(r'(^[^\.]+)\.(\d+)\.mafft$', os.path.basename(mafft))
    ogX = tuple([int(x) for x in og_catch[1].split('-')])
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
        return ogX, og, dnds
    except ZeroDivisionError:
        print('\t\t' + err + ' unknown error')
    except KeyError:
        print('\t\t' + err + ' ambiguous NTs')
    return None


def parse_dnds(mpRes):

    ogXdNdS_dict = {}
    for res in mpRes:
        if res:
            ogx, og, dnds = res[0], res[1], res[2]
            try:
                float(dnds)
            except ValueError:
                continue
            if ogx not in ogXdNdS_dict:
                ogXdNdS_dict[ogx] = {
                    'selection_total': 0,
                    'og_total': 0,
                    'ogs': {}
                    }

            if dnds < 1: # if purifying selection
                ogXdNdS_dict[ogx]['selection_total'] += dnds
            else: # if positive selection use reciprocal
                ogXdNdS_dict[ogx]['selection_total'] += 1/dnds
            ogXdNdS_dict[ogx]['og_total'] += 1
            ogXdNdS_dict[ogx]['ogs'][og] = dnds
    for ogx in ogXdNdS_dict:
        selection_coefficient = ogXdNdS_dict[ogx]['selection_total']/ogXdNdS_dict[ogx]['og_total']
        mean_dnds = sum(ogXdNdS_dict[ogx]['ogs'].values())/len(ogXdNdS_dict[ogx]['ogs'])
        ogXdNdS_dict[ogx] = [selection_coefficient, mean_dnds, ogXdNdS_dict[ogx]['ogs']]

    return ogXdNdS_dict


def hash_clusters(ome, gff_path, ome_sig_clus, gene2og, clusplusminus):
    '''sliding window of size clusplusminus, identify sequences that may have
    ogs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus'''

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
                            og = gene2og[seq1]
                        except KeyError:
                            continue
                        if og in sig_clus[0] and start is None: # if the og is
                        # in the sig clus and we haven't started
                            start = i1 # this is the start index
                        elif og in sig_clus[0]: # otherwise if it is in the
                        # sig clus label this the other border unless another is 
                        # found
                            end = i1 + 1
                    ogX = tuple(sorted(list(sig_clus[0])))
                    clus_out.append([
                         ogX, locus[start:end], {sig_clus[1]}, scaf
                         ])
                    # clus_out = [(og0, og1, ... ogn), [prot0, prot1, ... protn])]
    
    return ome, clus_out, cds_dict


def merge_clusters(clus_out):
    '''because we have significant clusters with borders, we can merge each
    cluster and assume that at a certain number of ogs in a combo will never
    happen via random sampling; therefore, their significance in the
    microsynteny adjustment is sufficient for justification for merging'''

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
                ogX = tuple(sorted(list(set(clus0[0]).union(set(clus1[0])))))
                # obtain the higher order og combination as the union of the
                # two cluster sets, then sort it and store as a tuple
                comp_clus.append([
                    ogX, list(loc0.union(loc1)), 
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


def write_clusters(ome_sig_clus, ome, out_file, gff_path, gene2og, clusplusminus = 10):
    # is multiprocessing this messing something and memory and that's why
    # output is wonky? (sometimes fewer proteins than should be in locus) 
    # or is output consistent fucking up elsewhere

    clus_out, cds_dict, ome = hash_clusters(ome, gff_path, ome_sig_clus, gene2og, clusplusminus)
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
                clusOGs.append(gene2og[gene])
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
        out.write('#name\togs\tgenes\tclans\n')
        for clus in clus_out:
            clus[2] = sorted(clus[2])
            ogs = ','.join([str(x) for x in clus[0]]) # og1,og2,...ogn
            genes = ','.join(clus[1]) # prot1,prot2,prot3,...protn
            clans= ';'.join([str(x) for x in list(clus[2])])
            out.write(clus[4] + '\t' + ogs + '\t' + genes + '\t' + clans + '\n') 
            out_genes.append(clus[1])

    return ome, out_genes


def write_cluster_scores(out_dir, clus_probs, i2ome, ogx2omes): #, rand_clus_probs):
    
    out_list, clus_len_dict = [], {}
    for ogX in clus_probs:
        out_list.append([
            ogX, clus_probs[ogX], 
            [i2ome[z] for z in ogx2omes[ogX]] #, rand_clus_probs[ogX]
            ])
#        clus_len_dict[ogX] = len(ogx2omes[ogX])
         
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
       os.rename(hmm_dir + 'pfam.tmp', hmm_dir + 'pfam.out')
    

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


def grabClus(genes_list, gff_path, prot_path, ome, ome_dir, gene2og, pfamRes = {}):

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
                ogStr = ';OG=' + str(gene2og[gene])
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
        out.write('#name\togs\tgenes\togx_clans\tpfams\n')
        with open(ome_dir + ome + '/info.out', 'r') as raw:
            for line in raw:
                if not line.startswith('#'):
                    clus = line.split('\t')[0]
                    out.write(line.rstrip() + '\t' + clus_anns[clus] + '\n')
    os.rename(ome_dir + ome + '/info.out.tmp', ome_dir + ome + '/info.out')

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
            'og_file\t' + log_dict['og_file'] + '\n' + \
            'plusminus\t' + str(log_dict['plusminus']) + '\n' + \
            'pair_percentile\t' + str(log_dict['pair_percentile']) + '\n' + \
            'ogx_percentile\t' + str(log_dict['ogx_percentile']) + '\n' + \
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
            log_res['ogx_percentile'] = False
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
        seed_arr = wrk_dir + 'ogpair.arr.npy'
        if os.path.isfile(seed_file):
            os.remove(seed_file)
        if os.path.isfile(seed_arr):
            os.remove(seed_arr)
        clus_pickle = wrk_dir + 'ogx2loc.pickle'
        ogx_pickle = wrk_dir + 'ogx_scores.pickle'
        ome_pickle = wrk_dir + 'ogx_omes.pickle'
        if os.path.isfile(clus_pickle):
            os.remove(clus_pickle)
        if os.path.isfile(ogx_pickle):
            os.remove(ogx_pickle)
        if os.path.isfile(ome_pickle):
            os.remove(ome_pickle)
    if not log_res['border_percentile']:
        row_file = wrk_dir + 'mtx/mcl.prep.rows'
        prep_file = wrk_dir + 'mtx/mcl.prep.gz'
        if os.path.isfile(row_file):
            os.remove(row_file)    
        if os.path.isfile(prep_file):
            os.remove(prep_file)
    if not log_res['ogx_percentile']:
        kern_file = out_dir + 'ogx_clans.tsv.gz'
        clus_file = out_dir + 'ogxs.tsv.gz'
        patch_pickle = wrk_dir + 'patchiness.scores.pickle'
        ome_dir = wrk_dir + 'ome/'
        ogx_dir = wrk_dir + 'ogx/'
        hmm_dir = wrk_dir + 'hmm/'
        if os.path.isfile(kern_file):
            os.remove(kern_file)
        if os.path.isfile(clus_file):
            os.remove(clus_file)
        if os.path.isdir(ome_dir):
            shutil.rmtree(ome_dir)
        if os.path.isdir(ogx_dir):
            shutil.rmtree(ogx_dir)
        if os.path.isdir(hmm_dir):
            shutil.rmtree(hmm_dir)
        if os.path.isfile(wrk_dir + 'ogx.tar.gz'):
            os.remove(wrk_dir + 'ogx.tar.gz')
        if os.path.isfile(patch_pickle):
            os.remove(patch_pickle)
    if not log_res['og_file']:
        shutil.rmtree(wrk_dir)
        os.mkdir(wrk_dir)
        for key in log_res:
            log_res[key] = False



def logCheck(log_dict, log_path, out_dir, wrk_dir):

    if not os.path.isfile(log_path):
        log_res = {x: False for x in log_dict}
        rmOldData(log_res, out_dir, wrk_dir)
        initLog(log_path, log_dict)
    log_res = readLog(log_path, log_dict)
    if any(not log_res[x] for x in log_res):
        print('\nInitializing new run', flush = True)
        rmOldData(log_res, out_dir, wrk_dir)
        initLog(log_path, log_dict)

    return log_res


def PatchMain(
    phylo, ogx2omes, ogxs, wrk_dir, 
    old_path = 'patchiness.scores.pickle', cpus = 1
    ):

    if not os.path.isfile(wrk_dir + old_path):
        if isinstance(ogx2omes, dict): # round 1 patchiness
            clusOmes = set([tuple([str(x) for x in ogx2omes[y]]) for y in ogxs])
        else: # round 2 is a list
            clusOmes = set([
                tuple([str(x) for x in y]) for i, y in enumerate(ogxs)
                ])
#        more = set([tuple(omes) for omes in moduleOmes])
#        allOGxs = list(clusOgxs.union(more))    
        with mp.Pool(processes = cpus) as pool:
            patch_res = pool.starmap(
                calcPatchiness, [(phylo, x) for x in clusOmes]
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

def RBHmain(
    ogx2loc, wrk_dir, ome2i, og_dir, ogx_dir,
    diamond, db, gene2og, plusminus, og2gene, 
    old_path = 'ogx2rbh.pickle',
    modules = None, moduleOGxs = None, cpus = 1
    ):

    ogs_list = []
    for ogX in ogx2loc:
        ogs_list.extend(list(ogX))
    ogs = set(ogs_list)
    if not os.path.isfile(wrk_dir + old_path):
        if not modules: # for ogxs
            ogx2rbh = rbhMngr2(
                list(ogs), list(ome2i.keys()), og_dir, ogx_dir, diamond, ogx2loc,
                db, gene2og, plusminus, og2gene, cpus = cpus
                )
        else: # run the kernel detection version
            ogx2rbh = rbhMngr3(
                list(ogs), list(ome2i.keys()), og_dir, ogx_dir, diamond, ogx2loc, 
                db, gene2og, plusminus, og2gene, modules,
                moduleOGxs, ome2i, cpus = cpus
                )
        with open(wrk_dir + old_path, 'wb') as pickout:
            pickle.dump(ogx2rbh, pickout)
        ogx_dirTar = mp.Process(target=tardir, args=(ogx_dir, True))
        ogx_dirTar.start() # when to join ...
    else:
        print('\tLoading previous coevolution results', flush = True)
        with open(wrk_dir + old_path, 'rb') as pickin:
            ogx2rbh = pickle.load(pickin)

    return ogx2rbh


def main(
    db, og_file, out_dir, plusminus = 1, og_dir = None,
    seed_perc = 0.2, clus_perc = 0.7, ogx_perc = 0.7,
    minimum_omes = 2, samples = 10000, pfam = None,
    tree_path = None, diamond = 'diamond',
    run_dnds = False, cpus = 1, n50thresh = None, 
    root = None, coevo_thresh = 0, patch_thresh = 0
    ):
    '''
    The general workflow:
    log management -> input data parsing -> orthogroup pair identification ->
    microsynteny distance thresholding -> OGx formation ->
    microsynteny distance border thresholding -> OGx grouping ->
    microsynteny distance cluster thresholding -> patchiness calculation ->
    coevolution calculation -> optional dN/dS calculations ->
    OGx data output -> cluster retrieving -> data output
    '''

    # initialize the log and working directory
    wrk_dir = out_dir + 'working/'
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)
    log_path = out_dir + 'log.txt'
    log_dict = {
        'og_file': og_file, 'plusminus': plusminus,
        'pair_percentile': seed_perc, 'ogx_percentile': clus_perc,
        'border_percentile': ogx_perc,
        'null_samples': samples, 'n50': n50thresh
        }
    log_res = logCheck(log_dict, log_path, out_dir, wrk_dir)

    # if not using an external tree
    if not tree_path:
        tree_path = out_dir + 'microsynt.newick'

    # obtain useable omes
    useableOmes, dbOmes = set(), set(db['ome'])
    print('\nI. Inputting data', flush = True)
    if n50thresh: # optional n50 threshold via mycotoolsDB
        assemblyPath = formatPath('$MYCODB/../data/assemblyStats.tsv')
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
    print('\tParsing orthogroups', flush = True)
    ome2i, gene2og, i2ome, og2gene = compileOGs(og_file, wrk_dir, useableOmes)

    cc_arr_path = wrk_dir + 'ogpair'
    if not os.path.isfile(out_dir + 'seed_scores.tsv.gz'):
        # compile cooccuring pairs of orthogroups in each genome
        print('\tCompiling all loci', flush = True)
        pairs = compileLoci(
            db, ome2i, gene2og, plusminus, 
            cpus = cpus
            )
    
        # assimilate cooccurrences across omes
        print('\tIdentifying cooccurences', flush = True) 
    
        seed_len = sum([len(pairs[x]) for x in pairs])
        print('\t\t' + str(seed_len) + ' initial OG-pairs', flush = True)
        cooccur_array, cooccur_dict, ogpairs, revogpairs = \
            formCooccurStructures(pairs, 2, cc_arr_path)
        max_ome = max([len(cooccur_dict[x]) for x in cooccur_dict])
        print('\t\t' + str(max_ome) + ' maximum organisms with OG-pair', flush = True)
        cooccur_array[cooccur_array > 0] = 1
        cooccur_array.astype(np.int8)
        print('\t\t' + str(sys.getsizeof(cooccur_array)/1000000) + ' MB', flush = True)
        np.save(cc_arr_path, cooccur_array)

    else:
        cooccur_array = np.load(cc_arr_path + '.npy')
        # do I even need this considering the above is with respect to seed score filtering

    db = db.set_index('ome')
 
    if not os.path.isfile(wrk_dir + '2.null.txt'):
        print('\tGenerating null distributions', flush = True)
        omes2dist, pair_null = gen_nulls(
            pairs,
            phylo, samples = samples, cpus = cpus
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


    microsynt_dict = {}
    if not os.path.isfile(tree_path):
        # create microsynteny distance matrix and make tree
        print('\tPreparing microsynteny matrix', flush = True)
        microsynt_dict = PrepMicrosyntDict(cooccur_array)
    print('\tCompiling microsynteny tree', flush = True)
    phylo = compileTree(
        microsynt_dict, i2ome, tree_path, root = root
        )

    # create null distribution for orthogroup pairs

    # seed clusters by calcuating total microsynteny distance for 
    # orthogroup pairs
    print('\nII. Seeding OG pairs', flush = True) 
    if not os.path.isfile(out_dir + 'seed_scores.tsv.gz'):
        print('\tCalculating seed OG-pair scores', flush = True)
        seed_score_start = datetime.now()

        results = calcDists(phylo, cooccur_dict, cpus)
        omes2dist, top_ogs = {x[1]: x[0] for x in results}, []
        with open(wrk_dir + 'ome_scores.pickle', 'wb') as out:
            pickle.dump(omes2dist, out)
        for ogpair in cooccur_dict:
            score = omes2dist[cooccur_dict[ogpair]]
            if score > min_score:
                i = revogpairs[ogpair]
                top_ogs.append([
                    ogpair[0], ogpair[1], len(cooccur_dict[ogpair]), score#, score/ome_dist
                    ])

        # write scores            
        print('\t\t' + str(datetime.now() - seed_score_start), flush = True)
        top_ogs = sorted(top_ogs, key = lambda x: x[3])
        print('\tWriting seed scores', flush = True)
        with gzip.open(out_dir + 'seed_scores.tsv.gz', 'wt') as out:
            out.write('#og0\tog1\tcooccurrences\tscore\tadj_score\n')
            for line in top_ogs:
                out.write('\t'.join([str(x) for x in line]) + '\n')
        print('\t\t' + str(len(top_ogs)) + ' significant seeds', flush = True)
    
    elif not os.path.isfile(wrk_dir + 'ogx_scores.pickle'): # load previous og pairs
        print('\tLoading previous seed OG pairs', flush = True)
        top_ogs = load_seedScores(out_dir + 'seed_scores.tsv.gz', min_score)
        print('\t\t' + str(len(top_ogs)) + ' significant seeds', flush = True)


    # begin sifting for OGxs using pairs as seeds for OGx detection
    print('\nIII. Sprouting high order OG combinations (OGx)', flush = True)
    if not os.path.isfile(wrk_dir + 'ogx_omes.pickle'):
        ogPair_dict = form_ogPairDict(top_ogs)
        print('\tForming OGxs', flush = True)
        form_clus_start = datetime.now()
        ogx2omes, ogx2loc = ogPair2ogX(
            db, ogPair_dict, gene2og, 
            ome2i, cpus, clusplusminus = plusminus
            )
        with open(wrk_dir + 'ogx2loc.pickle', 'wb') as pickout:
            pickle.dump(ogx2loc, pickout)
        with open(wrk_dir + 'ogx_omes.pickle', 'wb') as pickout:
            pickle.dump(ogx2omes, pickout)
        formOGxTime = datetime.now() - form_clus_start
        print('\t\t' + str(formOGxTime), flush = True)

    else: # or just load available structures
        with open(wrk_dir + 'ogx2loc.pickle', 'rb') as pickin:
            ogx2loc = pickle.load(pickin)
        with open(wrk_dir + 'ogx_omes.pickle', 'rb') as pickin:
            ogx2omes = pickle.load(pickin)

    ome_combos = set([tuple(sorted(list(x))) for x in list(ogx2omes.values())])
    if not os.path.isfile(wrk_dir + 'ogx_scores.pickle'):
        ogx2i, i2ogx, ogx_cooccur_dict = FormCooccurDict(
            ogx2omes
            ) # create ogx data structures

        print('\tCalculating OGx microsynteny distances', flush = True)
        clus_score_start = datetime.now()
        clus_obs = len(ogx2omes)
        print('\t\t' + str(clus_obs) + ' observed OGx', flush = True)

        # calculate OGx microsynteny distances
        results = calcDists(phylo, ogx_cooccur_dict, cpus, omes2dist = omes2dist)
        ogx2dist = {}
        omes2dist, top_ogs = {**omes2dist, **{x[1]: x[0] for x in results}}, []
        with open(wrk_dir + 'ome_scores.pickle', 'wb') as out:
            pickle.dump(omes2dist, out)
        for ogX in ogx_cooccur_dict:
            score = omes2dist[tuple(ogx_cooccur_dict[ogX])]
            ogx2dist[ogX] = score
    
        print('\t\t' + str(datetime.now() - clus_score_start), flush = True)
        with open(wrk_dir + 'ogx_scores.pickle', 'wb') as pickout:
            pickle.dump(ogx2dist, pickout)

    else: # just load previous ogx scores
        print('\tLoading previous OGxs', flush = True)
        with open(wrk_dir + 'ogx_scores.pickle', 'rb') as pickin:
            ogx2dist = pickle.load(pickin) 

    if not os.path.isfile(wrk_dir + 'ogx2pval.pickle'):
        print('\tEstimating OGx probability', flush = True)
        omes2og2genes, omes2genes = {}, defaultdict(list)
        for gene, og in gene2og.items():
            ome = gene[:gene.find('_')]
            omeI = ome2i[ome]
            if omeI not in omes2og2genes:
                omes2og2genes[omeI] = defaultdict(list)
            omes2og2genes[omeI][og].append(gene)
            omes2genes[omeI].append(gene)
        unadjOGx2pval = HypergeoMNGR(
            ogx2omes, omes2og2genes, omes2genes, (plusminus*2)-1, cpus = cpus
            )
   #     comparisons = len(ogx2omes)
  #      ogx2pval = {
 #           k: v * comparisons for k, v in unadjOGx2pval.items()
#            } # apply bonferroni correction
        with open(wrk_dir + 'ogx2pval.pickle', 'wb') as out:
            pickle.dump(ogx2pval, out)
    else:
        with open(wrk_dir + 'ogx2pval.pickle', 'rb') as raw:
            ogx2pval = pickle.load(raw)
    

    # prepare null distributions for each size (# of OGs) observed
    # in OGxs    
    max_ogx_size = max([len(x) for x in ogx2omes])
    if not os.path.isfile(wrk_dir + str(max_ogx_size) + '.null.txt'):
        print('\tPreparing OGx nulls', flush = True)
        bordScores, clusScores = genOGxNulls(
            list(db['gff3']), gene2og, max_ogx_size, plusminus, phylo,
            ogx_perc, clus_perc, wrk_dir,
            omes2dist, samples = samples, cpus = cpus
            )
    else: # or just load the old ones
        print('\tLoading OGx nulls', flush = True)
        bordScores, clusScores = loadNulls(max_ogx_size, wrk_dir, ogx_perc, clus_perc)

    # collect for normalizing relative to the whole dataset later
    # should I move into absolute space? it may just be better for comparing
    # datasets in the future and some algorithm is going to pick up that gene clusters
    # are in some absolute microsynteny distance, so long as its not normalized. 
#    max_dist = max(ogx2dist.values())
 #   min_dist = min(ogx2dist.values())

    i2ogx, ogx2i, togx2dist, count = {}, {}, {}, 0
    for ogx in list(ogx2dist.keys()):
        if ogx2dist[ogx] >= bordScores[len(ogx)]:
             togx2dist[ogx] = ogx2dist[ogx]
             i2ogx[count], ogx2i[ogx] = ogx, count
             count += 1
             # apply the border threshold for ALL pre-grouped ogxs
    print(
        '\t\t' + str(len(ogx2dist)) + ' OGx pass border threshold', 
        flush = True
        )
    ogx2dist = togx2dist

    print('\nIV. Quantifying OGx patchiness', flush = True)
    omes2patch = PatchMain(
        phylo, ogx2omes, ogx2dist, wrk_dir,
        old_path = 'patchiness.scores.pickle', cpus = cpus
        )

    ogx_dir = wrk_dir + 'ogx/'
    if not checkdir(ogx_dir, unzip = True, rm = True):
        os.mkdir(ogx_dir)

    print('\nV. Quantifying OGx coevolution', flush = True)
    ogx2rbh = RBHmain(
        ogx2loc, wrk_dir, ome2i, og_dir, ogx_dir,
        diamond, db, gene2og, plusminus, og2gene, 
        old_path = 'ogx2rbh.pickle', cpus = cpus
        )


    # Group ogxs
    print('\nVI. Inferring OGx clans', flush = True)
    if not os.path.isfile(wrk_dir + 'fams.pickle'): # need to add this to
    # log parsing
        fams, famOmes, famOGxs, omes2dist = groupOGx(
            ogx2loc, db, gene2og, i2ogx, ogx2i,
            phylo, clusScores, bordScores, ome2i,
            ogx2omes, wrk_dir, #ogx2dist,
            omes2dist = omes2dist, clusplusminus = plusminus, cpus = cpus
            )
        with open(wrk_dir + 'fams.pickle', 'wb') as pickout:
            pickle.dump([fams, famOmes, famOGxs], pickout)
        with open(wrk_dir + 'ome_scores.pickle', 'wb') as pickout:
            pickle.dump(omes2dist, pickout)
    else: # or just load old data
        with open(wrk_dir + 'fams.pickle', 'rb') as raw:
            fams, famOmes, famOGxs = pickle.load(raw)



    runOmes = [
        omes for omes in famOmes \
        if omes not in omes2patch
        ]
    runOGxs = [
        ogx for i, ogx in enumerate(famOGxs) \
        if famOmes[i] not in omes2patch
        ]
   
    print('\nVI. Quantifying cluster family patchiness', flush = True)
    omes2patch = {**PatchMain(
        phylo, runOGxs, runOmes, wrk_dir,
        old_path = 'patchiness.full.pickle', cpus = cpus
        ), **omes2patch} # could make more efficient by skipping redos

    ogx_dir = wrk_dir + 'ogx/'
    if not checkdir(ogx_dir, unzip = True, rm = True):
        os.mkdir(ogx_dir)

    print('\nVII. Quantifying cluster family gene evolution congruence', flush = True)
    ogx2rbh = {**ogx2rbh, **RBHmain(
        ogx2loc, wrk_dir, ome2i, og_dir, ogx_dir,
        diamond, db, gene2og, plusminus, og2gene, 
        old_path = 'ogx2rbh.full.pickle',
        fams = fams, famOGxs = famOGxs,
        cpus = cpus
        )}

    if coevo_thresh > 0 or patch_thresh > 0 or microsyn_thresh > 0:
        print('\tApplying thresholds', flush = True)
        print('\t\t' + str(len(fams)) + ' families before', flush = True)
        # edit fams, ogx2dist
 #       newOgx2dist, newFams, newFamOmes, newFamOGxs = {}, [], [], []
        newFams, newFamOmes, newFamOGxs = [], [], []
        for i0, fam in enumerate(fams):
            check = False # have we added a new list
            ogc = famOGxs[i0]
            omesc = famOmes[i0]
            for x, omes in fam.items():
                if ogx2rbh[ogc] >= coevo_thresh \
                    and omes2patch[omesc] >= patch_thresh:
                    if check:
                        newFams[-1][x] = omes
                    else:
                        newFams.append({x: omes})
                        check = True
            if check:
                newFamOmes.append(famOmes[i0])
                newFamOGxs.append(ogc)
        fams, famOmes, famOGxs = \
            newFams, newFamOmes, newFamOGxs
        print('\t\t' + str(len(fams)) + ' families after', flush = True)

    if run_dnds: # need to bring file naming to speed
        print('\nIIX. Quantifying cluster family dn/ds', flush = True)
        omes4ogs, ogx_files = {x: defaultdict(list) for x in range(len(ome2i))}, []
        for i, fam in enumerate(fams):
            for ogx in fam:
#                ogX = ogx2i[x]
                for gene in ogx2loc[ogx]:
                    omeI = ome2i[gene[:gene.find('_')]]
                    if omeI in set(fam[i][ogx]):
                        omes4ogs[omeI][gene].append(famOGxs[i])

        ogx_files = []
        for ogx in famOGxs:
            for og in ogx:
                out_base = ogx_dir + '-'.join([str(x) for x in ogx]) + '.' + str(og)
                ogx_files.append(out_base)

        maffts = [file + '.mafft' for file in ogx_files]
        if not all(os.path.isfile(x + '.aa.fa') for x in ogx_files):
            print("\tCompiling OGx kernels' genes", flush = True)
            dnds_preparation(db, omes4ogs, ogx2omes, gene2og, plusminus, ogx_dir, i2ome)
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
        with mp.Pool(processes = cpus) as pool:
            dnds_res = pool.map(calc_dnds, maffts)
        pool.join()
    
        ogx2dnds = parse_dnds(dnds_res)
    else:
        print('\nIIX. Skipping quantifying selection', flush = True)
        ogx2dnds = {}

    print('\nIX. Writing and annotating clusters', flush = True)
    print('\tWriting cluster scores', flush = True)
    fam_output, done, maxval = [], set(), max(ogx2dist.values())
    for fam, modOGx2omes in enumerate(fams):
        for ogX, omes in modOGx2omes.items():
            ogX_id = ogx2i[ogX]
            if ogX_id not in done:
                done.add(ogX_id)
            else:
                continue
            fam_output.append([
                ','.join([str(x) for x in ogX]), ogX_id,
                fam, ogx2dist[ogX]/maxval, ogx2rbh[ogX], omes2patch[tuple(ogx2omes[ogX])],
                ','.join([i2ome[x] for x in ogx2omes[ogX]])
                ]) # OGxs at this stage are not segregated into groups
    fam_output = sorted(fam_output, key = lambda x: x[3], reverse = True)
    with gzip.open(out_dir + 'ogxs.tsv.gz', 'wt') as out:
        out.write('#ogs\togX_id\tfam\tdistance\tcoevolution\tpatchiness\tomes')
        for entry in fam_output:
            out.write('\n' + '\t'.join([str(x) for x in entry]))

    kern_output, top_ogxs = [], []
    for i, ogX in enumerate(famOGxs):
        omesc = famOmes[i]
        try:
            if ogX not in dnds_dict:
                dnds_dict[ogX] = ['na' for x in range(3)]
        except NameError:
            dnds_dict = {ogX: ['na' for x in range(3)]}
        if ogX in ogx2rbh:
            kern_output.append([
                ','.join([str(x) for x in ogX]), i, 
                omes2dist[omesc]/maxval, ogx2rbh[ogX], omes2patch[omesc], 
                ','.join([str(i2ome[x]) for x in omesc])
#                dnds_dict[ogX][0], dnds_dict[ogX][1], str(dnds_dict[ogX][2]),
#                omes2dist[omesc]
                ])
        else:
            kern_output.append([
                ','.join([str(x) for x in ogX]), i,
                omes2dist[omesc]/maxval, ogx2rbh[ogX], omes2patch[omesc], 
                ','.join([str(i2ome[x]) for x in omesc]) 
#                dnds_dict[ogX][0], dnds_dict[ogX][1], str(dnds_dict[ogX][2]),
               # omes2dist[omesc]
                ])
        for sogx, omes in fams[i].items():
            top_ogxs.append([sogx, i, set(omes)])

    kern_output = sorted(kern_output, key = lambda x: x[4], reverse = True)


    with gzip.open(out_dir + 'ogx_fams.tsv.gz', 'wt') as out:
        out.write(
            '#ogs\tfam\tdistance\t' + \
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
    sig_clus = exSigClus(
        ogx2dist, ogx2loc, top_ogxs, ome2i
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
        write_clus_cmds.append([sig_clus[ome], ome, out_file, gff, gene2og, plusminus])
    with mp.Pool(processes = cpus) as pool:
        out_genes = pool.starmap(write_clusters, write_clus_cmds)
    pool.join()
    genes = {x[0]: x[1] for x in out_genes if x[1]}

#    ogx_dirTar.join()

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
            grabClus_cmds.append([genes[ome], gff_path, pro_path, ome, ome_dir, gene2og, pfamRes[ome]])
        except KeyError:
            grabClus_cmds.append([genes[ome], gff_path, pro_path, ome, ome_dir, gene2og])


    with mp.Pool(processes = cpus) as pool:
        clus_info = pool.starmap(grabClus, grabClus_cmds)
    pool.join()

#    print('\tCluster SVGs', flush = True)
 #   ome_svg_cmds = []
  #  svg_width = 10
   # for ome in clus_info:
    #    ome_svg_cmds.append([ome_dir + ome + '/'])
#    with mp.Pool(processes = cpus) as pool:
 #       pool.starmap(runGFF2SVG, ome_svg_cmds)


if __name__ == '__main__':
    description = \
    '''Pharmaceuticals are primarily derived from biologically-produced compounds, their
    derivatives, and synthetic compounds inspired by natural pharmacophores.
    Many natural product specialized metabolites are produced by gene clusters, or 
    regions of genomes that contain colocalized genes with concerted 
    biosynthetic function. 
    Genomes rearrange and recombine through macroevolutionary timescales, so the 
    persistence of colocalized gene families within clusters is indicative of 
    selection. 
    og2clus detects these regions by identifying unexpected conserved microsynteny 
    of orthogroups throughout a sample of organisms. 
    Some unexpectedly conserved regions are not directly involved in specialized 
    metabolite biosynthesis - og2clus filters results by implementing a
    multivariate model derived from 1) total distance between organisms with a 
    particular orthogroup combination, 2) phylogenetic topological correlation between the
    constituent genes of the orthogroup combination, and 3) phylogenetic patchiness of
    the orthogroup combination.'''
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-d', '--database', required = True, default = masterDB(), 
        help = 'MycotoolsDB. DEFAULT: masterdb')
    parser.add_argument('-og', '--orthogroups', required = True, 
        help = 'OrthoFinder output directory')
    parser.add_argument('-w', '--window', default = 5, type = int,
        help = 'Max genes +/- for each locus window. DEFAULT: 5 (11 gene window)')
    parser.add_argument('-r', '--root', 
        help = 'Ome code or 2 ome codes to root upon, e.g. psicya1, ' + \
        '"psicya1 psicub1"; DEFAULT: midpoint')
    parser.add_argument('-sp', '--seed_percentile', type = int, default = 75,
        help = 'Percentile of OG pair distances. DEFAULT: 75')
    parser.add_argument('-op', '--ogx_percentile', type = int,
        help = 'Percentile of OGx microsynteny distances. ' + \
        'Must be less than -cp. DEFAULT: -cp')
    parser.add_argument('-cp', '--clus_percentile', type = int, default = 80, 
        help = 'Percentile of OGx microsynteny distances . DEFAULT: 80')
    parser.add_argument('-pt', '--patch_threshold', default = 0, type = float,
        help = "Threshold [0 < value < 1] of high order OG combinations' patchiness scores.")
    parser.add_argument('-ct', '--coevo_threshold', default = 0, type = float, 
        help = "Threshold [0 < value < 1] of high order OG combinations' coevolution scores.")
    parser.add_argument('-p', '--pfam', required = True, help = 'Pfam-A.hmm path')
#    parser.add_argument('-p', '--phylogeny', 
 #       help = 'Species phylogeny (newick). DEFAULT: generate neighbor joining tree from ' + \
  #      'orthogroup pairs')
    parser.add_argument('-n', '--null_sample', type = int, default = 10000,
        help = 'Samples for null distributions. DEFAULT: 10,000')
    parser.add_argument('-s', '--dnds', action = 'store_true', help = 'Run dN/dS calculations')
    parser.add_argument('--n50', help = 'Minimum assembly N50')
    parser.add_argument('-o', '--output')
    parser.add_argument('--cpus', default = 1, type = int)
    args = parser.parse_args()

    if args.root:
        if '"' in args.root or "'" in args.root:
            args.root = args.root.replace('"','').replace("'",'')
        root = args.root.split()
        root_txt = ','.join(root)
    else:
        root = []
        root_txt = 'midpoint'

    if not args.ogx_percentile:
        args.ogx_percentile = args.clus_percentile # set the default

  
    args_dict = {
        'OrthoFinder results': args.orthogroups, 'MycotoolsDB': args.database, 
        'Root': root_txt, 'Pfam DB': args.pfam, 'Window': args.window*2+1,
        'Seed Percentile': args.seed_percentile, #'Precluster threshold': args.clus_threshold,
        'Cluster percentile': args.clus_percentile, 'OGx percentile': args.ogx_percentile,
        'Patchiness threshold': args.patch_threshold, 'Coevolution threshold': args.coevo_threshold,
        'Null samples': args.null_sample, 'Calculate dN/dS': args.dnds, 'Minimum N50': args.n50,
        'Processors': args.cpus, 'Output directory': args.output
        }

    pfam = formatPath(args.pfam)
    if not os.path.isfile(pfam):
        print('\nERROR: invalid Pfam-A.hmm path', flush = True)
        sys.exit(4)
    findExecs(
        ['mafft', 'hmmsearch', 'diamond', 'mcxload', 'mcxdump', 'mcl'], # mcl], 
        exit = {'mafft', 'hmmsearch', 'diamond', 'mcxload', 'mcxdump', 'mcl'} #, mcl}
        )
    if args.clus_percentile < args.ogx_percentile:
        print('\nERROR: ogx percentile is greater than cluster percentile',
            flush = True)
        sys.exit(3)

    start_time = intro('orthogroups to clusters - og2clus', args_dict, 'Zachary Konkel, Jason Slot')
    date = datetime.strftime(start_time, '%Y%m%d')

    if not args.output:
        out_dir = os.getcwd() + '/' + date + '_og2clus/'
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
    else:
        out_dir = formatPath(args.output)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
            out_dir += '/'

#    if args.phylogeny:
 #       tree_path = formatPath(args.phylogeny)
  #  else:
   #     tree_path = None
    if len(str(args.seed_percentile)) > 1:
        seed_perc = float('.' + str(args.seed_percentile))
    else:
        seed_perc = float('.0' + str(args.seed_percentile))
    if len(str(args.clus_percentile)) > 1:
        clus_perc = float('.' + str(args.clus_percentile))
    else:
        clus_perc = float('.0' + str(args.clus_percentile))   
    if len(str(args.ogx_percentile)) > 1:
        ogx_perc = float('.' + str(args.ogx_percentile))
    else:
        ogx_perc = float('.0' + str(args.ogx_percentile))   


    db = mtdb(formatPath(args.database))
    orthogroups = formatPath(args.orthogroups) + '/Orthogroups/Orthogroups.txt'
    og_dir = formatPath(args.orthogroups) + '/Orthogroup_Sequences/'
    main(
        db, orthogroups, out_dir, plusminus = args.window,
        cpus = args.cpus, og_dir = og_dir,
        seed_perc = seed_perc, #clus_thresh = args.clus_threshold,
        clus_perc = clus_perc, diamond = 'diamond',#seed_thresh = args.seed_threshold,
        ogx_perc = ogx_perc, pfam = pfam, samples = args.null_sample,
        run_dnds = args.dnds, n50thresh = args.n50,
        root = root, coevo_thresh = args.coevo_threshold, 
        patch_thresh = args.patch_threshold 
        )

    outro(start_time)
