import os
import sys
import shutil
import random
import pickle
import multiprocessing as mp
from tqdm import tqdm
from cogent3 import PhyloNode
from itertools import combinations, chain
from collections import defaultdict, Counter
from mycotools.lib.biotools import gff2list
from mycotools.lib.kontools import eprint, collect_files
from cloci.lib import treecalcs
from cloci.lib.input_parsing import compileCDS

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


def gen_null_dict_omes(db, i2ome, combo_dict, sample = 10000, omes = []):
    """combo_dict = {ome: [(OG0...OGn)]}"""

    # make a list of observed HGx combinations 
    nulls = []
    for ome in omes:
        combo_dict[ome] = list(combo_dict[ome])
        nulls.extend(combo_dict[ome]) # extend all combos to the null
    
    # randomly sample while accounting for oversampled spp
    # NO DO NOT DO THIS. This approach was designed to account for sampling bias
    # in null distributions that have oversampled taxa, e.g. Saccharomyces, Cryptococcus
    # However, lineages that are not monophyletic in the microsynteny tree will have
    # a chance of small distances being enough to be declared unexpected because
    # the sampling bias wasn't considered.
    spp2omes, null_list = defaultdict(list), []
    [spp2omes[db[i2ome[x]]['taxonomy']['species']].append(x) for x in omes]
    spp = list(spp2omes.keys())
    sp_omes = list(chain(*list(spp2omes.values())))
    for i in range(sample):
#        rand_sp = random.choice(spp)
#        rand_ome = random.choice(spp2omes[rand_sp])
        rand_ome = random.choice(sp_omes)
        null_list.append(random.choice(combo_dict[rand_ome]))
        
    reps = {k: v for k, v in Counter(null_list).items() if v > 1}
    null_set = set(null_list)
                                            
    cooccur_dict = defaultdict(list)
    for ome, hgxs in combo_dict.items():
        intersect = list(set(hgxs).intersection(null_set))
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


def rank_partition(db, wrk_dir, rank, ome2i, partition_file):
    lineages = defaultdict(list)
    for ome, row in db.items():
        lineage = row['taxonomy'][rank]
        if lineage:
            if ome in ome2i:
                lineages[row['taxonomy'][rank]].append(ome)
    with open(partition_file, 'w') as out:
        for lineage, omes in lineages.items():
            out.write(f'#{lineage}\n{" ".join(omes)}\n')


def load_partitions(partition_file, ome2i):
    omes = set(ome2i.keys())
    with open(partition_file, 'r') as raw:
        partition_omes = [[ome2i[y] for y in x.rstrip().split() if y in omes] \
                          for x in raw \
                          if x.rstrip() and not x.startswith('#')]
    return partition_omes


def gen_nulls(db, i2ome, pairs, phylo, omes = [], samples = 10000, 
              uniq_sp = False, dist_func = treecalcs.calc_tmd,
              cpus = 1):
    """pairs = {ome: [(OG0,OG1)...]}"""

    null_dict, reps = gen_null_dict_omes(db, i2ome, pairs, samples, omes = omes)
    oldLen = len(null_dict)
    null_dict = {k: v for k, v in null_dict.items() if len(v) > 1}
    if not uniq_sp:
        results = treecalcs.calc_dists(phylo, null_dict, cpus = cpus,
                                       func = dist_func)
    else:
        results = treecalcs.calc_dists(phylo, null_dict, cpus = cpus, i2ome = i2ome,
                                       func = treecalcs.calc_tmd, uniq_sp = uniq_sp)
    omes2dist, pair_scores = {x[1]: x[0] for x in results}, []
    for k in null_dict:
        pair_scores.append(omes2dist[null_dict[k]])
    for k, iters in reps.items():
        for i in range(1, iters):
            pair_scores.append(omes2dist[null_dict[k]])
    for i in range(samples - len(pair_scores)):
        pair_scores.append(0)

    return omes2dist, sorted(pair_scores)


def gen_pair_nulls(db, phylo, ome2i, wrk_dir, nul_dir, seed_perc, ome2pairs,
                   i2ome, samples = 1000, partition_file = None, partition_rank = None,
                   dist_func = treecalcs.calc_tmd, uniq_sp = False, cpus = 1):
    if not os.path.isdir(nul_dir):
        os.mkdir(nul_dir)
    if partition_file or partition_rank:
        if partition_rank:
            partition_file = f'{wrk_dir}{partition_rank}_partitions.txt'
            if not os.path.isfile(partition_file):
                rank_partition(db, wrk_dir, partition_rank, ome2i, partition_file)
        print('\tPartitioning null distributions', flush = True)
        partition_omes = load_partitions(partition_file, ome2i)
        ome2partition = {}
        for k, omes in enumerate(partition_omes):
            for ome in omes:
                ome2partition[ome] = k
        missing_omes = set(ome2i.values()).difference(set(ome2partition.keys()))
        ome2partition = {**ome2partition, **{k: None for k in missing_omes}}
    else: # spoof it
        partition_omes = [ome2i.values()]
        ome2partition = {i: 0 for i in ome2i.values()}

    if os.path.isfile(wrk_dir + 'omes2dist.pickle'):
        with open(wrk_dir + 'omes2dist.pickle', 'rb') as raw:
            omes2dist = pickle.load(raw)
    else:
        omes2dist = {}

    if not os.path.isfile(f'{wrk_dir}hgx2loc.pickle'):
        # create null distribution for orthogroup pairs
        final_partition = len(partition_omes) - 1 
        pair_nulls, omes2dist = [], {}
        if not os.path.isfile(f'{nul_dir}2.null.{final_partition}.txt'):
            print('\tGenerating null distributions', flush = True)
            for i, omes in enumerate(partition_omes):
                omes2dist, pair_null = gen_nulls(
                    db, i2ome, ome2pairs, phylo,
                    omes, samples = samples, 
                    dist_func = dist_func, uniq_sp = uniq_sp,
                    cpus = cpus
                    )
                n_f = f'{nul_dir}2.null.{i}.txt'
                with open(f'{n_f}.tmp', 'w') as out:
                    out.write('\n'.join([str(i0) for i0 in pair_null]))
                shutil.move(f'{n_f}.tmp', n_f)
                pair_nulls.append(pair_null)
            with open(wrk_dir + 'omes2dist.pickle', 'wb') as out:
                pickle.dump(omes2dist, out)
        else: # or just load what is available
            if os.path.isfile(wrk_dir + 'omes2dist.pickle'):
                with open(wrk_dir + 'omes2dist.pickle', 'rb') as raw:
                    omes2dist = pickle.load(raw)
            for i, omes in enumerate(partition_omes):
                with open(f'{nul_dir}2.null.{i}.txt', 'r') as raw:
                   pair_nulls.append([float(x.rstrip()) for x in raw])
        min_pair_scores = []
        for i, pair_null in enumerate(pair_nulls):
            scores_i = round(seed_perc * len(pair_null) + .5)
            min_pair_scores.append(pair_null[scores_i])
    else:
        min_pair_scores = None

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
        elif os.path.isfile(f'{nul_dir}skip.{size}-{partition_num}.txt'):
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
        s, i = [int(x) for x in skip_f.replace('skip.','').replace('.txt','').split('-')]
        if s > skip_i[i]:
            skip_i[i] = s
    return {k: v for k,v in skip_i.items()}


def gen_hgx_nulls(
    db, i2ome, gffs, gene2hg, max_clus_size, plusminus, phylo,
    hgx_perc, clus_perc, nul_dir, partition_omes,
    omes2dist = {}, skip_i = {}, samples = 10000, 
    dist_func = treecalcs.calc_tmd, uniq_sp = False, cpus = 1
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
                pool.close()
                pool.join() 
            size_dict = {ome: nulls for ome, nulls in hashRes}
            size_dict = {k: size_dict[k] for k in sorted(size_dict.keys())}
        for i0, omes in enumerate(partition_omes):
            if i0 in skip_i:
                if skip_i[i0] <= size:
                    hgxBordPercs[i0][size] = 0.0
                    hgxClusPercs[i0][size] = 0.0
                    continue
            # multiprocessing
            if os.path.isfile(f'{nul_dir}{size}.null.{i0}.txt'):
                nullSizes = load_null(f'{nul_dir}{size}.null.{i0}.txt')
            else:
                null_dict, reps = gen_null_dict_omes(db, i2ome, size_dict, samples, omes)
                omes2dist = treecalcs.update_dists(phylo, null_dict, func = dist_func,
                                                   uniq_sp = uniq_sp, i2ome = i2ome,
                                                   omes2dist = omes2dist, cpus = cpus)
                nullSizes = []
                with open(f'{nul_dir}{size}.null.{i0}.txt.tmp', 'w') as out:
                    for omes in list(null_dict.values()):
                        nullSizes.append(omes2dist[omes])
                        out.write(str(omes2dist[omes]) + '\n')
                    count = 0
                    for k, v in reps.items():
                        for i1 in range(v - 1):
                            t_dist = omes2dist[null_dict[k]]
                            out.write(str(t_dist) + '\n')
                            count += 1
                            nullSizes.append(t_dist)
                    z_vals = range(samples - len(null_dict) - count)
                    out.write('\n'.join([str(0) for x in z_vals]))
                    nullSizes.extend([0 for x in z_vals])
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
                    with open(f'{nul_dir}skip.{s}-{i0}.txt', 'w') as out:
                        pass
                print(f'\t\t\tWARNING: lineage {i0}: all HGx >= {size} HGs pass')

    #skips = [f'{k}-{v}' for k,v in skip_i.items()]
#    done_file = f'{nul_dir}done.{".".join(skips)}'
    done_file = f'{nul_dir}done.txt'
    with open(done_file, 'w') as out:
        pass

    return hgxBordPercs, hgxClusPercs


def partitions2hgx_nulls(db, partition_omes, ome2i, i2ome, gene2hg, max_hgx_size,
                         plusminus, hgx_perc, clus_perc, nul_dir, omes2dist,
                         phylo, samples = 1000, dist_func = treecalcs.calc_tmd,
                         uniq_sp = False, cpus = 1):
    final_partition = len(partition_omes) - 1
    print('\tPreparing HGx nulls', flush = True)
    files = [os.path.basename(x) for x in collect_files(nul_dir, 'txt')]
    done_file = [x for x in files if x.startswith('done.')]
    skip_files = [x for x in files if x.startswith('skip.')]
    skip_i = parse_skip_files(skip_files)
    if done_file:
        bordScores_list, clusScores_list = [], []
#        skip_i = parse_done_file(done_file)
        for i, omes in enumerate(partition_omes):
            bordScores, clusScores = load_hgx_nulls(max_hgx_size, nul_dir, hgx_perc,
                                                    clus_perc, i)
            bordScores_list.append(bordScores)
            clusScores_list.append(clusScores)
    else:
        bordScores_list, clusScores_list = gen_hgx_nulls(
            db, i2ome, {db[k]['gff3']: ome2i[k] for k in list(db.keys()) \
             if k in ome2i},
            gene2hg, max_hgx_size, plusminus, phylo,
            hgx_perc, clus_perc, nul_dir, partition_omes,
            omes2dist, skip_i, samples = samples, cpus = cpus,
            dist_func = dist_func, uniq_sp = uniq_sp
            )

    return bordScores_list, clusScores_list
