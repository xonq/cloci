import os
import sys
import shutil
import random
import pickle
import multiprocessing as mp
from tqdm import tqdm
from cogent3 import PhyloNode
from itertools import combinations
from collections import defaultdict, Counter
from mycotools.lib.biotools import gff2list
from mycotools.lib.kontools import eprint, collect_files
from orthocluster.orthocluster.lib import phylocalcs
from orthocluster.orthocluster.lib.input_parsing import compileCDS

def hash_4_nulls_bysize(
    gff_path, gene2hg, size, plusminus, ome
    ):

    gff_list = gff2list(gff_path) # open here for decreased serialization
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))
  
    nulls = []
    for scaf, seqs in cds_dict.items():
        for i, seq in enumerate(seqs):
            if i < plusminus: # compile the locus
                locus = seqs[:i+plusminus+1]
            else:
                locus = seqs[i-plusminus:i+plusminus+1]
            loc_og = set([gene2hg[x] for x in locus if x in gene2hg]) # translate to OGs
            nulls.extend((
                [tuple(sorted(x)) for x in combinations(loc_og, size)]
                ))

    return ome, set(nulls)


def hash_4_nulls(
    gff_path, gene2hg, max_size, plusminus, ome = None
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
    return ome, nulls


def form_cooccur_dict(cooccur_dict):

    # sort the values by the length of the combination
    cooccur_dict = {k: tuple(sorted(v)) for k, v in sorted(cooccur_dict.items(), key = lambda x: len(x[0]))}

    hgx2i, i2hgx = {}, {}
    for i, hgx in enumerate(list(cooccur_dict.keys())):
        i2hgx[i] = hgx
        hgx2i[hgx] = i

    return i2hgx, hgx2i, cooccur_dict


def gen_null_dict_omes(combo_dict, sample = 10000, omes = []):
    """combo_dict = {ome: [(OG0...OGn)]}"""
 
    nulls = []
    for ome in omes:
        combos = combo_dict[ome]
        nulls.extend(list(combos)) # extend all combos to the null
     
    try: # sample
        null_list = random.sample(nulls, sample)
    except ValueError:
        print('\t\t\t\tWARNING: requested null sample greater than HGx size', flush = True)
        null_list = nulls

    reps = {k: v for k, v in Counter(null_list).items() if v > 1}
    null_set = set(null_list)
                                            
    cooccur_dict = defaultdict(list)
    for ome, hgxs in combo_dict.items():
        intersect = list(hgxs.intersection(null_set))
        for hgx in intersect:
            cooccur_dict[hgx].append(ome)
    
    cooccur_dict = {k: tuple(v) for k, v in cooccur_dict.items() \
                    if len(v) > 1}
    # {(OG0,OGn): [omei0...]}
    reps = {k: v for k, v in reps.items() \
            if k in cooccur_dict}

    return cooccur_dict, reps


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


def load_partitions(partition_file, ome2i):
    omes = set(ome2i.keys())
    with open(partition_file, 'r') as raw:
        partition_omes = [[ome2i[y] for y in x.rstrip().split() if y in omes] \
                          for x in raw \
                          if x.rstrip() and not x.startswith('#')]
    return partition_omes


def partition_for_nulls(phylo, partition_file, ome2i):
    partition_omes_prep = load_partitions(partition_file, ome2i)
    complete_set = {int(x.name) for x in phylo.iter_tips()}
    check_set = set()
        
    partition_omes, ome2partition = [], {}
    for i, partitions in enumerate(partition_omes_prep):
        part_phy = phylo.lowest_common_ancestor([str(x) for x in partitions]) # MRCA
        part_phy_omes = [int(x.name) for x in part_phy.iter_tips()]
        if set(part_phy_omes).intersection(check_set):
            eprint('\tWARNINGS: MRCA of partitions overlap: ' \
                 + str(partitions), flush = True)
            sys.exit(13)
        check_set = check_set.union(set(part_phy_omes))
        partition_phylos.append(part_phy)
        partition_omes.append(part_phy_omes)
        ome2partition = {**ome2partition, **{x: i for x in partition_omes[-1]}}
                
    diff_list = list(complete_set.difference(check_set))
    if diff_list:   
        eprint('\t\tWARNING: omes missing from partitions: ' \
             + str(','.join([str(x) for x in diff_list])),
             flush = True)
        ome2partition = {**ome2partition, **{x: None for x in diff_list}}
    return partition_phylos, partition_omes, ome2partition


def gen_nulls(pairs, phylo, omes = [], samples = 10000, cpus = 1):
    """pairs = {ome: [(OG0,OG1)...]}"""
            
    null_dict, reps = gen_null_dict_omes(pairs, samples, omes = omes)
    oldLen = len(null_dict)
    null_dict = {k: v for k, v in null_dict.items() if len(v) > 1}
    results = phylocalcs.calc_dists(phylo, null_dict, cpus = cpus)
    omes2dist, pair_scores = {x[1]: x[0] for x in results}, []
    for k in null_dict:
        pair_scores.append(omes2dist[null_dict[k]])
    for k, iters in reps.items():
        for i in range(1, iters):
            pair_scores.append(omes2dist[null_dict[k]])
    for i in range(samples - len(pair_scores)):
        pair_scores.append(0)

    return omes2dist, sorted(pair_scores)


def gen_pair_nulls(phylo, ome2i, wrk_dir, nul_dir, seed_perc, ome2pairs,
                   i2ome, samples = 1000, partition_file = None, cpus = 1):
    if partition_file:
        print('\tPartitioning tree for null distributions', flush = True)
        partition_omes = load_partitions(partition_file, ome2i)
        ome2partition = {}
        for k, omes in enumerate(partition_omes):
            for ome in omes:
                ome2partition[ome] = k
        missing_omes = set(ome2i.values()).difference(set(ome2partition.keys()))
        ome2partition = {**ome2partition, **{k: None for k in missing_omes}}
#            partition_for_nulls(phylo, partition_file, ome2i)
    else: # spoof it
        partition_omes = [ome2i.values()]
        ome2partition = {i: 0 for i in ome2i.values()}

    if os.path.isfile(wrk_dir + 'omes2tmd.pickle'):
        with open(wrk_dir + 'omes2tmd.pickle', 'rb') as raw:
            omes2dist = pickle.load(raw)

    # create null distribution for orthogroup pairs
    final_partition = len(partition_omes) - 1 
    pair_nulls, omes2dist = [], {}
    if not os.path.isfile(f'{nul_dir}2.null.{final_partition}.txt'):
        print('\tGenerating null distributions', flush = True)
        for i, omes in enumerate(partition_omes):
            omes2dist, pair_null = gen_nulls(
                ome2pairs, phylo,
                omes, samples = samples, cpus = cpus
                )
            n_f = f'{nul_dir}2.null.{i}.txt'
            with open(f'{n_f}.tmp', 'w') as out:
                out.write('\n'.join([str(i0) for i0 in pair_null]))
            shutil.move(f'{n_f}.tmp', n_f)
            pair_nulls.append(pair_null)
        with open(wrk_dir + 'omes2tmd.pickle', 'wb') as out:
            pickle.dump(omes2dist, out)
    else: # or just load what is available
        if os.path.isfile(wrk_dir + 'omes2tmd.pickle'):
            with open(wrk_dir + 'omes2tmd.pickle', 'rb') as raw:
                omes2dist = pickle.load(raw)
        for i, omes in enumerate(partition_omes):
            with open(f'{nul_dir}2.null.{i}.txt', 'r') as raw:
               pair_nulls.append([float(x.rstrip()) for x in raw])
    min_pair_scores = []
    for i, pair_null in enumerate(pair_nulls):
        scores_i = round(seed_perc * len(pair_null) + .5)
        min_pair_scores.append(pair_null[scores_i])

    return partition_omes, ome2partition, \
           omes2dist, min_pair_scores


def load_null(nul_f):
    with open(nul_f, 'r') as raw:
        nullSizes = [float(x.rstrip()) for x in raw]
    nullSizes.sort()
    return nullSizes


def load_hgx_nulls(max_clus_size, nul_dir, hgx_perc, clus_perc, partition_num):

    hgxBordPercs, hgxClusPercs = {}, {}
    for size in range(3, max_clus_size + 1):
        nul_f = f'{nul_dir}{size}.null.{partition_num}.txt'
        if os.path.isfile(nul_f):
            with open(nul_f, 'r') as raw:
                nullSizes = [float(x.rstrip()) for x in raw]
            nullSizes.sort()
            hgxBordPercs[size] = nullSizes[round(hgx_perc * len(nullSizes) + .5)]
            hgxClusPercs[size] = nullSizes[round(clus_perc * len(nullSizes) + .5)]
        elif os.path.isfile(f'{nul_dir}skip.{size}-{partition_num}'):
            hgxBordPercs[size] = 0.0
            hgxClusPercs[size] = 0.0
        else:
            raise FileNotFoundError(f'missing null/skip file in {nul_dir}')
    
    return hgxBordPercs, hgxClusPercs


def parse_done_file(done_file):
    # NEED some accounting of changing bord percentile
    skip_i = {}
    done_data = done_file.split('.')
    if done_data:
        for d in done_data[1:]:
            k,v = d.split('-')
            skip_i[int(k)] = int(v)
    return skip_i


def parse_skip_files(skip_files):
    # NEED some accounting of changing bord percentile
    skip_i = defaultdict(int)
    for skip_f in skip_files:
        s, i = [int(x) for x in skip_f.replace('skip.','').split('-')]
        if s > skip_i[i]:
            skip_i[i] = s
    return {k: v for k,v in skip_i.items()}


def gen_hgx_nulls(
    gffs, gene2hg, max_clus_size, plusminus, phylo,
    hgx_perc, clus_perc, nul_dir, partition_omes,
    omes2dist = {}, skip_i = {}, samples = 10000, cpus = 1
    ):
    """Generates a null distribution of randomly sampled HGxs for each 
    # of homogroups observed in HGxs. Applies the percentile for each
    size as the minimum value for significance.
    Outputs a dictionary of the minimum values for each size HGx
    based on the inputted percentiles. {# of OGs: minimum value}"""

    # hashRes = [({size: [HGx]})...] by ome
    print('\tCalculating random sample TMD', flush = True)
    hgxBordPercs, hgxClusPercs = defaultdict(dict), defaultdict(dict)
    for size in range(3, max_clus_size + 1): # for all observed sizes
        if all(i in skip_i for i, v in enumerate(partition_omes)):
            if all(skip_i[i] <= size for i, v in enumerate(partition_omes)):
                for i, v in enumerate(partition_omes):
                    hgxBordPercs[i][size] = 0.0
                    hgxClusPercs[i][size] = 0.0
                continue
        print(f'\t\tHGx size {size}', flush = True)
        if not all(os.path.isfile(f'{nul_dir}{size}.null.{i}.txt') \
               or os.path.isfile(f'{nul_dir}skip.{size}-{i}') \
               for i in range(len(partition_omes))):
            print('\t\t\tRandomly sampling HGxs', flush = True)
            hash_null_cmds = [
                (x, gene2hg, size, plusminus, i) \
                for x, i in gffs.items()
                ]
            with mp.get_context('fork').Pool(processes = cpus) as pool:
                hashRes = pool.starmap(hash_4_nulls_bysize, 
                                       tqdm(hash_null_cmds, total = len(hash_null_cmds)))
            size_dict = {ome: nulls for ome, nulls in hashRes}
            size_dict = {k: size_dict[k] for k in sorted(size_dict.keys())}
        for i0, omes in enumerate(partition_omes):
            if i0 in skip_i:
                if skip_i[i0] <= size:
                    hgxBordPercs[i0][size] = 0.0
                    hgxClusPercs[i0][size] = 0.0
                    continue
            # multiprocessing
            print(f'\t\t\tLineage {i0}', flush = True)
            if os.path.isfile(f'{nul_dir}{size}.null.{i0}.txt'):
                nullSizes = load_null(f'{nul_dir}{size}.null.{i0}.txt')
            else:
                null_dict, reps = gen_null_dict_omes(size_dict, samples, omes)
                omes2dist = phylocalcs.update_dists(phylo, null_dict, 
                                                    omes2dist = omes2dist, cpus = cpus)
                nullSizes = []
                with open(f'{nul_dir}{size}.null.{i0}.txt.tmp', 'w') as out:
                    for omes in list(null_dict.values()):
                        nullSizes.append(omes2dist[omes])
                        out.write(str(omes2dist[omes]) + '\n')
                    count = 0
                    for k, v in reps.items():
                        for i1 in range(v - 1):
                            out.write(str(omes2dist[null_dict[k]]) + '\n')
                            count += 1
                    out.write('\n'.join([str(0) for x in range(samples-len(null_dict)-count)]))
                    # account for missing values
                shutil.move(f'{nul_dir}{size}.null.{i0}.txt.tmp',
                            f'{nul_dir}{size}.null.{i0}.txt')
        
            # sort to apply percentile threshold for each size HGx and for both the
            # border percentile and cluster percentile
                nullSizes.sort()
            hgxBordPercs[i0][size] = nullSizes[round(hgx_perc * len(nullSizes) + .5)]
            hgxClusPercs[i0][size] = nullSizes[round(clus_perc * len(nullSizes) + .5)]
            if not hgxBordPercs[i0][size]:
                skip_i[i0] = size
                for s in range(size, max_clus_size + 1):
                    with open(f'{nul_dir}skip.{s}-{i0}', 'w') as out:
                        pass
                print(f'\t\t\t\tWARNING: all HGx >= {size} HGs pass')

    skips = [f'{k}-{v}' for k,v in skip_i.items()]
    done_file = f'{nul_dir}done.{".".join(skips)}'
    with open(done_file, 'w') as out:
        pass

    return hgxBordPercs, hgxClusPercs


def partitions2hgx_nulls(db, partition_omes, ome2i, gene2hg, max_hgx_size,
                         plusminus, hgx_perc, clus_perc, nul_dir, omes2dist,
                         phylo, samples = 1000, cpus = 1):
    final_partition = len(partition_omes) - 1
    print('\tPreparing HGx nulls', flush = True)
    files = [os.path.basename(x) for x in collect_files(nul_dir, 'txt')]
    done_file = [x for x in files if x.startswith('done.')]
    if done_file:
        bordScores_list, clusScores_list = [], []
        skip_i = parse_done_file(done_file)
        for i, omes in enumerate(partition_omes):
            bordScores, clusScores = load_hgx_nulls(max_hgx_size, nul_dir, hgx_perc,
                                                    clus_perc, i)
            bordScores_list.append(bordScores)
            clusScores_list.append(clusScores)
    else:
        skip_files = [x for x in files if x.startswith('skip.')]
        skip_i = parse_skip_files(skip_files)
        bordScores_list, clusScores_list = gen_hgx_nulls(
            {db[k]['gff3']: ome2i[k] for k in list(db.keys()) \
             if k in ome2i},
            gene2hg, max_hgx_size, plusminus, phylo,
            hgx_perc, clus_perc, nul_dir, partition_omes,
            omes2dist, skip_i, samples = samples, cpus = cpus
            )

    return bordScores_list, clusScores_list
