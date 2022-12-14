import random
import pickle
from cogent3 import PhyloNode
from mycotools.lib.biotools import gff2list
from mycotools.lib.kontools import eprint
from orthocluster.orthocluster.lib import phylocalcs
from orthocluster.orthocluster.lib.input_parsing import compileCDS

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



def form_cooccur_dict(cooccur_dict):

    # sort the values by the length of the combination
    cooccur_dict = {k: tuple(sorted(v)) for k, v in sorted(cooccur_dict.items(), key = lambda x: len(x[0]))}

    hgx2i, i2hgx = {}, {}
    for i, hgx in enumerate(list(cooccur_dict.keys())):
        i2hgx[i] = hgx
        hgx2i[hgx] = i

    return i2hgx, hgx2i, cooccur_dict


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


def partition_for_nulls(phylo, partition_file):
    with open(partition_file, 'r') as raw:
        partition_omes = [x.rstrip().split() \
                          for x in raw \
                          if x.rstrip() and not x.startswith('#')]
    complete_set = {x for x in phylo.iter_tips()}
    check_set = {}
        
    partition_phylos, partition_omes, ome2partition = [], [], {}
    for i, partitions in enumerate(partion_omes):
        part_phy = phylo.last_common_ancestor(partitions) # MRCA
        part_phy_omes = [x for x in part_phy.iter_tips()]
        if set(part_phy_omes).intersection(check_set):
            eprint('\tERROR: MRCA tree of partition overlaps previous partition: ' \
                 + str(partitions), flush = True)
            sys.exit(13)
        check_set = check_set.union(set(part_phy_omes))
        partition_phylos.append(part_phy)
        partition_omes.append(part_phy_omes)
        ome2partition = {**ome2partition, **{x: i for x in partition_omes[-1]}}
                
    diff_list = list(complete_set.difference(check_set))
    if diff_list:   
        eprint('\t\tWARNING: omes missing from partitions: ' \
             + str(','.join(diff_list)),
             flush = True)
        ome2partition = {**ome2partition, **{x: None for x in diff_list}}
    return partition_phylos, partition_omes, ome2partition


def gen_nulls(pairs, phylo, samples = 10000, cpus = 1):
    """pairs = {ome: [(OG0,OG1)...]}"""
            
    hgpair2i, i2hgpair, null_dict = gen_null_dict(pairs, samples)
    oldLen = len(null_dict)
    null_dict = {k: v for k, v in null_dict.items() if len(v) > 1}
    results = phylocalcs.calc_dists(phylo, null_dict, cpus = cpus)
    omes2dist, pair_scores = {x[1]: x[0] for x in results}, []
    for k in null_dict:
        pair_scores.append(omes2dist[null_dict[k]])
    for i in range(oldLen - len(null_dict)):
        pair_scores.append(0)

    return omes2dist, sorted(pair_scores)


def gen_pair_nulls(phylo, i2ome, wrk_dir, nul_dir, seed_perc,
                   samples = 1000, partition_file = None, cpus = 1):
    if partition_file:
        print('\tPartitioning tree for null distributions', flush = True)
        partition_phylos, partition_omes, ome2partition = \
            partition_for_nulls(phylo, partition_file)
    else: # spoof it
        partition_phylos = [phylo]
        partition_omes = [i2ome]
        ome2partition = {ome: 0 for ome in i2ome}


    # create null distribution for orthogroup pairs
    final_partition = len(partition_phylos) - 1 
    pair_nulls = []
    if not os.path.isfile(f'{nul_dir}2.null.{final_partition}.txt'):
        print('\tGenerating null distributions', flush = True)
        for i, pphylo in enumerate(partition_phylos):
                omes2dist, pair_null = gen_nulls(
                    {o: ome2pairs[o] for o in partition_omes[i]},
                    pphylo, samples = samples, cpus = cpus
                    )
                with open(f'{nul_dir}2.null.{i}.txt', 'w') as out:
                    out.write('\n'.join([str(i0) for i0 in pair_null]))
                pair_nulls.append(pair_null)
        with open(wrk_dir + 'ome_scores.pickle', 'wb') as out:
                pickle.dump(omes2dist, out)
    else: # or just load what is available
        with open(wrk_dir + 'ome_scores.pickle', 'rb') as raw:
            omes2dist = pickle.load(raw)
        for i, pphylo in enumerate(paritition_phylos):
            with open(f'{nul_dir}2.null.{i}.txt', 'r') as raw:
               pair_nulls.append([float(x.rstrip()) for x in raw])
    min_pair_scores = []
    for i, pair_null in enumerate(pair_nulls):
        scores_i = round(seed_perc * len(pair_null) + .5)
        min_pair_scores.append(pair_null[scores_i])

    return partition_phylos, partition_omes, ome2partition, \
           omes2dist, min_pair_scores


def load_hgx_nulls(max_clus_size, nul_dir, hgx_perc, clus_perc, partition_num):

    hgxBordPercs, hgxClusPercs = {}, {}
    for size in range(3, max_clus_size + 1):
        nullSizes = []
        with open(f'{nul_dir}{size}.{partition_num}.null.txt', 'r') as raw:
            for line in raw:
                nullSizes.append(float(line.rstrip()))
        nullSizes.sort()
        hgxBordPercs[size] = nullSizes[round(hgx_perc * len(nullSizes) + .5)]
        hgxClusPercs[size] = nullSizes[round(clus_perc * len(nullSizes) + .5)]

    return hgxBordPercs, hgxClusPercs


def gen_hgx_nulls(
    gffs, gene2hg, max_clus_size, plusminus, phylo,
    hgx_perc, clus_perc, nul_dir, partition_num,
    omes2dist = {}, samples = 10000, cpus = 1
    ): # NEED to adjust; this currently uses way too much memory
    """Generates a null distribution of randomly sampled HGxs for each 
    # of homogroups observed in HGxs. Applies the percentile for each
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
        distRes = phylocalcs.calc_dists(phylo, null_dict, omes2dist = omes2dist, cpus = cpus)
        omes2dist, scores = {
            **omes2dist, **{x[1]: x[0] for x in distRes}
            }, []
        
        nullSizes = []
        with open(f'{nul_dir}{size}.{partition_num}.null.txt', 'w') as out:
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


def partitions2hgx_nulls(db, partition_omes, ome2i, gene2hg, max_hgx_size,
                         plusminus, hgx_perc, clus_perc, nul_dir, omes2dist,
                         samples = 1000, cpus = 1):
    bordScores_list, clusScores_list = [], []
    print('\tPreparing HGx nulls', flush = True)
    for i, pphylo in enumerate(partition_phylos):
        if not os.path.isfile(f'{nul_dir}{i}.{final_partition}.null.txt'):
            bordScores, clusScores = gen_hgx_nulls(
                [db[k]['gff3'] for k in partition_omes[i] if k in ome2i],
                gene2hg, max_hgx_size, plusminus, pphylo,
                hgx_perc, clus_perc, nul_dir, i,
                omes2dist, samples = samples, cpus = cpus
                )
        else:
            bordScores, clusScores = load_hgx_nulls(max_hgx_size, nul_dir, hgx_perc,
                                                clus_perc, i)
        bordScores_list.append(bordScores)
        clusScores_list.append(clusScores)
    return bordScores_list, clusScores_list
