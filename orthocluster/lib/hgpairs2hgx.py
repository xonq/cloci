import os
import pickle
import multiprocessing as mp
from datetime import datetime
from itertools import combinations
from collections import defaultdict
from mycotools.lib.biotools import gff2list
from orthocluster.orthocluster.lib import phylocalcs
from orthocluster.orthocluster.lib.input_parsing import compileCDS

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


def formulate_hgx_dicts(gen_clus_cmds, hgx2omes, hgx2loc, ome2i, cpus):

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


def id_hgx(db, hgpair_dict, gene2hg, ome2i, cpus, clusplusminus = 10):
    """form hgxs from hgpairs by extracting higher order combinations from
    overlapping hgx loci"""

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
                    
    # generate dictionaries of crude predicted hgxs
    hgx2omes, hgx2loc, gen_clus_cmds = {}, {}, []
    for hgpair in protohgx2omes:
        hgpair_set = set(hgpair)
        for loci in combinations(protohgx2omes[hgpair], 2):
            loc0, loc1 = loci[0], loci[1]
            gen_clus_cmds.append([loc0, loc1, hgpair_set])
        if len(gen_clus_cmds) > 5000000:
            count += len(gen_clus_cmds) 
            print('\t\t\t' + str(count), flush = True)
            hgx2omes, hgx2loc = formulate_hgx_dicts(
                gen_clus_cmds, hgx2omes, hgx2loc, ome2i, cpus
                )               
            gen_clus_cmds = []      
                            
    # run the last set of commands
    hgx2omes, hgx2loc = formulate_hgx_dicts(
        gen_clus_cmds, hgx2omes, hgx2loc, ome2i, cpus
        )

    print('\t\tRemoving subset HGxs', flush = True)
    hgx2omes, hgx2loc = rm_subsets(hgx2omes, hgx2loc, cpus)
    hgx2omes = {x: tuple(sorted(list(hgx2omes[x]))) for x in hgx2omes}

    return hgx2omes, hgx2loc


def hgpairs2hgx(db, wrk_dir, top_hgs, gene2hg, ome2i, 
                omes2dist, phylo, plusminus = 3, cpus = 1):
    if not os.path.isfile(wrk_dir + 'hgx_omes.pickle'):
        hgpair_dict = form_hgpairDict(top_hgs)
        print('\tForming HGxs', flush = True)
        form_clus_start = datetime.now()
        hgx2omes, hgx2loc = id_hgx(
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

    return hgx2omes, hgx2loc
