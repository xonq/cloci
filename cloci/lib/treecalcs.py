import os
import re
import random
import pickle
import multiprocessing as mp
from tqdm import tqdm
from cogent3 import PhyloNode
from itertools import chain
from collections import Counter, defaultdict
from mycotools.lib.kontools import eprint

def addPatch(phylo, omes_set):
    """recursively add branches with the trait to the presence total"""
    total = 0
    for sphylo in phylo:
        subOmes = [str(x)[:str(x).find(':')] for x in sphylo.tips(True)]
        if any(
            x in omes_set \
            for x in subOmes
            ): # if some descendents aren't in the set, search this branch
            if all(x in omes_set for x in subOmes):
                if len(subOmes) > 1:
                    total += sphylo.total_descending_branch_length()
                else:
                    total += sphylo.length
            else:
                total += addPatch(sphylo, omes_set)
    return total

def calc_pds(phylo, omes):
    """calculate the percent branch length that a trait is missing over the
    MRCA of where the trait is present"""
    omes_set = set(omes)
    try:
        mrca = phylo.lowest_common_ancestor(omes)
 #   except ValueError:
#        print(phylo, omes)
    except AttributeError: # 1 ome
        eprint('\t\t' + ','.join([str(x) for x in omes]) \
             + ' missing tip(s)', flush = True)
        print(phylo)

    mrca_omes = set(x.name for x in mrca.iter_tips())
    if mrca_omes == omes_set:
        return tuple([int(x) for x in omes]), 0
    else:
        totalDist = mrca.total_descending_branch_length()
        subDist = addPatch(mrca, omes_set)
        return tuple([int(x) for x in omes]), 1 - subDist/totalDist



def patch_main(
    phylo, omes, wrk_dir,
    old_path = 'pds.pickle', cpus = 1
    ):

    if os.path.isfile(wrk_dir + old_path):
        print('\tLoading previous PDS results', flush = True)
        with open(wrk_dir + old_path, 'rb') as in_pick:
            omes2patch = pickle.load(in_pick)
    else:
        omes2patch = {}

    clusOmes = set([
        tuple([str(x) for x in y]) for y in omes \
               if y not in omes2patch
        ])
    if clusOmes:
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            patch_res = pool.starmap(
                calc_pds, tqdm([(phylo, x) for x in clusOmes],
                                                 total = len(clusOmes))
                )
            pool.close()
            pool.join()
        omes2patch = {**omes2patch,
                      **{ome_tup: pds for ome_tup, pds in patch_res}}
        with open(wrk_dir + old_path, 'wb') as out:
            pickle.dump(omes2patch, out)
    
    return omes2patch

def calc_mmd(phylo, omes):
    mrca = phylo.lowest_common_ancestor([str(x) for x in omes])
    mmd = mrca.get_max_tip_tip_distance()[0]
    return mmd, tuple([int(x) for x in omes])

def calc_tmd(phylo, omes):
    """calculate descending branch length from cogent3 tree"""
    # need to verify wtf descending branch length v total of supplied nodes is
    omes = [str(x) for x in omes]
    omes_set = set(omes)

    print(phylo.get_tip_names(), flush = True)

    try: 
        mrca = phylo.lowest_common_ancestor(omes) # this is the failing step
        mrca_omes = set(x.name for x in mrca.iter_tips())
        tmd = addPatch(phylo.lowest_common_ancestor(omes),
                       omes_set)
        return tmd, tuple([int(i) for i in omes])
    except ValueError:
        eprint('\t\t' + ','.join([str(x) for x in omes]) \
             + ' missing/extraneous tip(s)', flush = True)
        return 0, tuple([int(i) for i in omes])
    except AttributeError:
        print(omes, '\n', phylo)

def calc_tmd_uniq_omes(phylo, omes, o_omes):
    """calculate descending branch length from cogent3 tree"""
    # need to verify wtf descending branch length v total of supplied nodes is
    omes = [str(x) for x in omes]
    omes_set = set(omes)
    try:
        mrca = phylo.lowest_common_ancestor(omes)
        mrca_omes = set(x.name for x in mrca.iter_tips())
        tmd = addPatch(phylo.lowest_common_ancestor(omes),
                       omes_set)
        return tmd, tuple([int(i) for i in o_omes])
    except ValueError:
        eprint('\t\t' + ','.join([str(x) for x in omes]) \
             + ' missing/extraneous tip(s)', flush = True)
        return 0, tuple([int(i) for i in o_omes])
    except AttributeError:
        print(omes, '\n', phylo)


def calc_branch_sim(phylo, omes0, omes1):
    tmd_union = calc_tmd(phylo, list(set(omes0).union(set(omes1))))[0]
    tmd_inter = calc_tmd(phylo, list(set(omes0).intersection(set(omes1))))[0]
    return tmd_inter/tmd_union


def get_uniq_spp(db, iomes, i2ome):
    omes = {i2ome[x]: x for x in iomes}
    ab_omes = [re.search(r'([^\d+])\d', x)[1] for x in omes]
    reps_set = set(k for k, v in Counter(ab_omes).items() if v > 1)
    # if there are duplicate omes, then check if the species is duplicate
    if reps_set:
        rep_spp = defaultdict(list)
        [rep_spp[db[x]['taxonomy']['species']].append(x) for x in omes \
         if re.search(r'([^\d+])\d', x)[1] in reps_set]
        # need the omes that weren't duplicated, the ones that were but aren't the 
        # same species, and a random choice of the duplicated ones
        # species and their omes that were duplicated
        dup_spp = {k: v for k, v in rep_spp.items() if len(v) > 1}
        dup_omes = chain(*list(dup_spp.values()))
        uniq_omes = list(set(omes.keys()).difference(set(dup_omes)))
        for sp, rep_omes in dup_spp.items():
            r_ome = rep_omes[random.randint(0, len(rep_omes) - 1)]
            uniq_omes.append(r_ome)
        return tuple(sorted([omes[i] for i in uniq_omes]))
    else:
        return iomes
        

def calc_dists(phylo, cooccur_dict, cpus = 1, omes2dist = {}, func = calc_tmd,
               uniq_sp = False, i2ome = None):
    # multiprocessing calculating only new distances for omes2dist
    if uniq_sp:
        with mp.get_context('forkserver').Pool(processes = cpus) as pool:
            results = pool.starmap(
                calc_tmd_uniq_omes,
                [(phylo, get_uniq_spp(uniq_sp, x, i2ome), x) \
                  for x in list(set(cooccur_dict.values())) \
                  if x not in omes2dist]
                )
            pool.close()
            pool.join()
    else:
        with mp.get_context('forkserver').Pool(processes = cpus) as pool:
            results = pool.starmap(
                func,
                [(phylo, x,) for x in list(set(cooccur_dict.values())) \
                if x not in omes2dist]
                )
            pool.close()
            pool.join()
        
    return results 


def update_dists(phylo, cooccur_dict, cpus = 1, omes2dist = {}, func = calc_tmd,
                 uniq_sp = [], i2ome = None):
    """update the omes2dist with a new set of distance calulations"""
    results = calc_dists(phylo, cooccur_dict, cpus, omes2dist = omes2dist, func = func,
                         uniq_sp = uniq_sp, i2ome = i2ome)
    omes2dist = {**omes2dist, **{x[1]: x[0] for x in results}}
    return omes2dist
