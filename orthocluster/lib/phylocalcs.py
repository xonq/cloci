import multiprocessing as mp
from cogent3 import PhyloNode
from mycotools.lib.kontools import eprint

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
 #   except ValueError:
#        print(phylo, omes)
    except AttributeError: # 1 ome
        eprint('\t\t' + ','.join([str(x) for x in omes]) + ' raised a tip not found error', flush = True)
        print(phylo)

    mrca_omes = set(x.name for x in mrca.iter_tips())
    if mrca_omes == omes_set:
        return tuple([int(x) for x in omes]), 1
    else:
        totalDist = mrca.total_descending_branch_length()
        subDist = addPatch(mrca, omes_set)
        return tuple([int(x) for x in omes]), subDist/totalDist


def calc_branch_len(phylo, omes):
    """calculate descending branch length from cogent3 tree"""
    # need to verify wtf descending branch length v total of supplied nodes is
    omes = [str(x) for x in omes]
    omes_set = set(omes)
    try:
        mrca = phylo.lowest_common_ancestor(omes)
        mrca_omes = set(x.name for x in mrca.iter_tips())
        if mrca_omes == omes_set:
            return mrca.total_descending_branch_length(), \
                   tuple([int(i) for i in omes])
        else:
            return addPatch(phylo.lowest_common_ancestor(
                omes
                ), set(omes)), tuple([int(i) for i in omes])
    except ValueError:
        eprint('\t\t' + ','.join([str(x) for x in omes]) + ' raised a tip not found error', flush = True)
        return 0, tuple([int(i) for i in omes])
    except AttributeError:
        print(omes, '\n', phylo)



def calc_dists(phylo, cooccur_dict, cpus = 1, omes2dist = {}):
    # multiprocessing calculating only new distances for omes2dist
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        results = pool.starmap(
            calc_branch_len,
            [(phylo, x,) for x in list(set(cooccur_dict.values())) \
            if x not in omes2dist]
            )
    
    return results 


def update_dists(phylo, cooccur_dict, cpus = 1, omes2dist = {}):
    """update the omes2dist with a new set of distance calulations"""
    results = calc_dists(phylo, cooccur_dict, cpus, omes2dist = omes2dist)
    omes2dist = {**omes2dist, **{x[1]: x[0] for x in results}}
    return omes2dist
