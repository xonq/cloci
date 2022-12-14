#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
import subprocess
import numpy as np
import multiprocessing as mp
from mycotools.db2files import soft_main as symlink_files
from mycotools.lib.kontools import format_path, mkOutput, \
                                   findExecs, intro, outro, \
                                   eprint
from mycotools.lib.dbtools import mtdb, masterDB
from orthocluster.orthocluster.lib.input_parsing import *

def run_mmseqs(db, wrk_dir, algorithm = 'mmseqs easy-cluster', 
                 min_id = 0.3, min_cov = 0.5, cpus = 1):
    symlink_files(['faa'], db, wrk_dir, verbose = False) # symlink proteomes
    cluster_res_file = wrk_dir + 'homolog_groups.tsv'
    if not os.path.isfile(cluster_res_file): # NEED to add to log removal
        # be extra cautious about shell injection because we need to glob
        int(cpus)
        float(min_id)
        float(min_cov)
        if not os.path.isdir(wrk_dir):
            raise OSError('invalid working directory')
        elif not algorithm in {'mmseqs easy-linclust', 'mmseqs easy-cluster'}:
            raise OSError('invalid mmseqs binary')
        mmseqs_cmd = subprocess.call(algorithm + ' ' + ' '.join([wrk_dir + 'faa/*faa',
                                      wrk_dir + 'cluster', wrk_dir + 'tmp/',
                                      '--min-seq-id', str(min_id), '--threads',
                                      str(cpus), '--compressed', '1',
                                      '--cov-mode', '0', '-c', str(min_cov)]),
                                      shell = True)
 #                                     stdout = subprocess.DEVNULL,
#                                      stderr = subprocess.DEVNULL)
        shutil.move(wrk_dir + 'cluster_cluster.tsv', cluster_res_file)
    elif os.path.getsize(cluster_res_file):
        mmseqs_cmd = 0
    else:
        mmseqs_cmd = 1
    if mmseqs_cmd:
        eprint('\tERROR: cluster failed')
        sys.exit(1)
    if os.path.isfile(wrk_dir + 'cluster_all_seqs.fasta'):
        os.remove(wrk_dir + 'cluster_all_seqs.fasta')
    if os.path.isfile(wrk_dir + 'cluster_rep_seq.fasta'):
        os.remove(wrk_dir + 'cluster_rep_seq.fasta')
    if os.path.isdir(wrk_dir + 'tmp/'):
        shutil.rmtree(wrk_dir + 'tmp/')
    return cluster_res_file

        
def parse_orthofinder(hg_file, useableOmes = set()):
    """
    imports orthofinder Orthogroups.txt "hg_file". outputs several data structures:
    ome_num = {ome: number}, gene2hg = {gene: og}, i2ome = [ome0, ome1, ome2]
    """

    gene2hg, ome_num, i2ome, hg2gene = \
        {}, {}, [], {}
    with open(hg_file, 'r') as raw:
        for line in raw:
            data = line.rstrip().split(' ') # OrthoFinder input
            og = int(data[0].replace(':','').replace('OG',''))
            hits = [x for x in data[1:] if x[:x.find('_')] in useableOmes]
            omes = [x[:x.find('_')] for x in hits]
            if len(set(omes)) < 2: # skip singleton organisms
                continue
            hg2gene[og] = hits
            for i, gene in enumerate(hits):
                ome = omes[i]
#                ome = gene[:gene.find('_')] # wtf, parsing omes raises and index error
                if ome not in ome_num:
                    i2ome.append(ome)
                    ome_num[ome] = len(i2ome) - 1
                gene2hg[gene] = og

    return ome_num, gene2hg, i2ome, hg2gene



def parse_1to1(hg_file, useableOmes = set()):
    derivations = defaultdict(list)
    with open(hg_file, 'r') as raw:
        for line in raw:
            k, v = line.rstrip().split() # all white space
            derivations[k].append(v)
  
    hgs_list = []
    for k, v in derivations.items():
        v.append(k) 
        hgs_list.append(sorted(set(v)))
        
    hg2gene = {i: v for i, v in enumerate(
                    sorted(hgs_list, key = lambda x: len(x),
                    reverse = True)
                    )}
    
    todel = []
    for og, genes in hg2gene.items():
        genes = [x for x in genes if x[:x.find('_')] in useableOmes]
        omes = set([x[:x.find('_')] for x in genes])
        if len(omes) < 2: # skip singleton omes
            todel.append(og)
    for og in todel:
        del hg2gene[og]    
    hg2gene = {k: v for k, v in sorted(hg2gene.items(), key = lambda x:
                                       len(x[1]), reverse = True)}
    
    gene2hg = {}
    for i, og_list in hg2gene.items():
        for gene in og_list:
            gene2hg[gene] = i
        
    i2ome = sorted(set([x[:x.find('_')] for x in list(gene2hg.keys())]))
    ome_num = {v: i for i, v in enumerate(i2ome)}
        
    return ome_num, gene2hg, i2ome, hg2gene

def compile_homolog_groups(hg_file, wrk_dir,
                           method = 'mmseqs easy-cluster',
                           useableOmes = set()):
    if method in {'mmseqs easy-cluster', 'mmseqs easy-linclust'}:
        ogInfo = parse_1to1(hg_file, useableOmes)
        hg2genes = ogInfo[-1]
        with open(wrk_dir + 'homolog_groups.txt', 'w') as out:
            for hg, genes in hg2genes.items():
                out.write(str(hg) + '\t' + ' '.join(genes) + '\n')
    elif method == 'orthofinder':
        ogInfo = parse_orthofinder(hg_file, useableOmes)

    with open(wrk_dir + 'ome2i.tsv', 'w') as out:
        out.write(
            '\n'.join([x + '\t' + str(ogInfo[0][x]) for x in ogInfo[0].keys()])
            )
    max_ome =  max([int(x) for x in ogInfo[1].values()])
    return ogInfo


def compile_loci(
    db, ome2i, gene2hg, plusminus, cpus = 1
    ):

    loci_hash_cmds = [
        [v['gff3'], ome2i[ome],
        gene2hg, plusminus]
        for ome, v in db.items() \
        if ome in set(ome2i)
        ]
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        loci_hashes = pool.starmap(parseLoci, loci_hash_cmds)
    pool.join()
    pairs = {x[0]: x[1] for x in loci_hashes}

    return pairs


def form_cooccur_array(cooccur_dict, ome_len):

    count, hgx2i, size_dict, cooccur_arrays, i2hgx = 0, {}, {}, {}, {}
    cooccur_dict = {k: tuple(sorted(v)) for k, v in sorted(cooccur_dict.items(), key = lambda x: len(x[0]))}
    cooccur_arrays = np.zeros([ome_len, len(cooccur_dict)], dtype = np.int32)

    old_len = len(list(cooccur_dict.keys())[0])
    for i, hgx in enumerate(list(cooccur_dict.keys())):
        i2hgx[i] = hgx
        hgx2i[hgx] = i
        for ome in cooccur_dict[hgx]:
#            cooccur_arrays[ome, i] = hgx_dict[ome][hgx]
            cooccur_arrays[ome, i] = 1

    return i2hgx, hgx2i, cooccur_dict, cooccur_arrays



def form_cooccur_structures(pairs, min_omes, ome_len, cc_arr_path = None):
    """
    Imports the out_dicts from parseLoci and creates index vectors that bear
    the ome_num's with a given  cooccurence. e.g. cooccur_dict[(og1, og2)] = [1, 2, 3]
    """

    cooccur_dict = defaultdict(list)
    for ome in pairs:
        for key in pairs[ome]:
            new_key = tuple(sorted(key))
            cooccur_dict[new_key].append(ome)

    cooccur_dict = {
        x: tuple(sorted(cooccur_dict[x])) for x in cooccur_dict \
        if len(cooccur_dict[x]) > min_omes
        }
            
    i2hgpair, hgpair2i, cooccur_dict, cooccur_array = form_cooccur_array(
        cooccur_dict, ome_len
        )   
            
    return cooccur_array, dict(cooccur_dict), hgpair2i, i2hgpair


def id_near_schgs(hg2gene, omes, max_hgs = 10000, max_median = 2, max_mean = 2):
    """Identify near single copy homology groups. Criteria:
       1) has all omes; 2) minimal overall size of homology group;
       3) lowest median"""
    near_schgs = []
    min_hg2gene = {k: v for k, v in sorted(hg2gene.items(),
                    key = lambda x: len(x[1])) \
                   if len(v) >= len(omes)}

    max_len = max_mean * len(omes) # dont need to compute means iteratively
    median_i = round((len(omes) / 2) - 0.5)

    for hg, genes in min_hg2gene.items():
        hg_omes = [x[:x.find('_')] for x in genes]
        if not omes.difference(hg_omes): # all omes are present
            count_omes = list(Counter(hg_omes).values()) # prepare for median calculation
            count_omes.sort()
            median = count_omes[median_i]
            print(median, max_median, len(genes), max_len)
            if median <= max_median and len(genes) <= max_len:
                print('\t', median, len(genes))
                near_schgs.append(hg)
            elif len(genes) >= max_len: # everything else will be higher
                break
            if len(near_schgs) == max_hgs:
                break

    return near_schgs


def extract_nschg_pairs(nschgs, hgpair2i, m_arr):
    nschg_set = set(nschgs)
    nschg_pairs = [i for hgpair, i in hgpair2i.items() \
                   if any(x in nschg_set for x in hgpair)]
    nschgs_arr = np.array(nschg_pairs)
    valid_arr = m_arr[:, nschgs_arr]
    return valid_arr


def align_microsynt_np(m_arr, i2ome, hg2gene, hgpair2i, wrk_dir, nschgs = None):
    if not nschgs:
        nschgs = id_near_schgs(hg2gene, set(i2ome), max_hgs = 100,
                           max_median = 4, max_mean = 3)
    trm_arr = extract_nschg_pairs(nschgs, hgpair2i, m_arr)
    with open(wrk_dir + 'microsynt.align.phy', 'w') as out:
        out.write(f'{trm_arr.shape[0]} {trm_arr.shape[1]}\n')
        [out.write(f'{i2ome[i]} {"".join([str(x) for x in trm_arr[i]])}\n') \
                   for i in range(trm_arr.shape[0])]
    return wrk_dir + 'microsynt.align.phy'


def remove_nulls(cc_arr):
    sum_arr = np.sum(cc_arr, axis = 1) # sum all ome rows
    null_i_list = list(np.where(sum_arr == 0)[0])
    del_list = sorted(null_i_list, reverse = True)
    for i in sorted(null_i_list, reverse = True):
        cc_arr = np.delete(cc_arr, i, axis = 0)
    return cc_arr, del_list


def run_tree(alignment, wrk_dir, constraint = False, iqtree = 'iqtree',
             model = 'GTR2+FO+ASC+R5', verbose = False, cpus = 1):
    
    tree_dir = wrk_dir + 'tree/'
    if not os.path.isdir(tree_dir):
        os.mkdir(tree_dir)
    
    prefix = tree_dir + 'microsynt'
    tree_cmd = [iqtree, '-s', alignment, '-m', model,
                '--prefix', prefix] 
    if constraint:
        tree_cmd.extend(['-te', constraint])
        comp_file = prefix + '.treefile'
    else: # use bootstrap
        tree_cmd.extend(['-B', '1000'])     
        comp_file = prefix + '.contree'
        
    if verbose:                             
        v = None                            
    else:
        v = subprocess.DEVNULL
    tree_res = subprocess.call(tree_cmd, stdout = v)
    shutil.copy(comp_file, wrk_dir + '../microsynt.newick')
        
    return tree_res


def main(db, hg_file, out_dir, wrk_dir, algorithm, 
         tree_path, plusminus = 3, min_cov = 0, min_id = 0.3, 
         n50thresh = None, near_single_copy_genes = [], 
         constraint = None, verbose = False, cpus = 1):

    # obtain useable omes 
    useableOmes, dbOmes = set(), set(db.keys())
    print('\nI. Inputting data', flush = True)
    if n50thresh: # optional n50 threshold via mycotoolsDB
        assemblyPath = format_path('$MYCODB/../data/assemblyStats.tsv')
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
    if not hg_file and not os.path.isfile(wrk_dir + 'homology_groups.tsv'):
        hg_file = run_mmseqs(db, wrk_dir, algorithm = algorithm,
                               min_id = min_id, cpus = cpus)
    print('\tParsing homology groups (HGs)', flush = True)
    ome2i, gene2hg, i2ome, hg2gene = compile_homolog_groups(hg_file, wrk_dir, 
                                                            algorithm, useableOmes)
    print('\t\tOmes:', len(ome2i), flush = True)
    print('\t\tHGs:', len(hg2gene), flush = True)
    print('\t\tGenes:', len(gene2hg), flush = True)
    
    cc_arr_path = wrk_dir + ''
    ome2pairs = compile_loci(
        db, ome2i, gene2hg, plusminus,
        cpus = cpus
        )

    if not os.path.isfile(out_dir + 'seed_scores.tsv.gz'):
        # compile cooccuring pairs of homogroups in each genome
        print('\tCompiling all loci', flush = True)
        
        # assimilate cooccurrences across omes
        print('\tIdentifying cooccurences', flush = True)
        
        seed_len = sum([len(ome2pairs[x]) for x in ome2pairs])
        print('\t\t' + str(seed_len) + ' initial HG-pairs', flush = True)
        cooccur_array, cooccur_dict, hgpair2i, i2hgpair = \
            form_cooccur_structures(ome2pairs, 2, len(ome2i), cc_arr_path)
        max_ome = max([len(cooccur_dict[x]) for x in cooccur_dict])
        print('\t\t' + str(max_ome) + ' maximum organisms with HG-pair', flush = True)
        cooccur_array[cooccur_array > 0] = 1
        cooccur_array.astype(np.int8)
        print('\t\t' + str(sys.getsizeof(cooccur_array)/1000000) + ' MB', flush = True)
        cooccur_array, del_omes = remove_nulls(cooccur_array)
        for i in del_omes:
            del i2ome[i]
        
        ome2i = {v: i for i, v in enumerate(i2ome)}
        with open(wrk_dir + 'ome2i.tsv', 'w') as out:
            out.write(
                '\n'.join([k + '\t' + str(v) for k, v in ome2i.items()])
                )
        np.save(cc_arr_path, cooccur_array)
#    elif not os.path.isfile(tree_path):
    else:
        cooccur_array = np.load(cc_arr_path + '.npy')
        cooccur_dict = None

    microsynt_dict = {}
    print('\nII. Building microsynteny tree', flush = True)
    if not os.path.isfile(tree_path):
        nschgs = []
        if near_single_copy_genes:
            try:
                for entry in near_single_copy_genes:
                    nschgs.append(int(entry))
                nschgs = sorted(set(nschgs))
            except ValueError:
                nschgs = sorted({gene2hg[gene] \
                             for gene in near_single_copy_genes \
                             if gene in gene2hg})
        # create microsynteny distance matrix and make tree
        print('\tPreparing microsynteny alignment', flush = True)
        align_file = align_microsynt_np(cooccur_array, i2ome, hg2gene,
                                        hgpair2i, wrk_dir, nschgs)
        print('\tBuilding microsynteny tree', flush = True)
        run_tree(align_file, wrk_dir, constraint = constraint, iqtree = 'iqtree',
                 model = 'GTR2+FO+ASC+R5', verbose = False, cpus = cpus)

    return ome2i, gene2hg, i2ome, hg2gene, ome2pairs, cooccur_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Generate a microsynteny ' \
      + 'tree that recapitulates the divergence of gene order within small ' \
      + 'loci.')
    parser.add_argument('-d', '--db', default = masterDB(), 
        help = 'MycotoolsDB. DEFAULT: masterdb')
    parser.add_argument('-f', '--focal_genes',
        help = 'File of genes for neighborhood extraction, e.g. SCOs')
    parser.add_argument('-of', '--orthofinder', 
        help = 'Precomputed OrthoFinder output directory')
    parser.add_argument('-i', '--input', 
        help = 'Precomputed cluster results file "gene\tcluster#"')
    parser.add_argument('-l', '--linclust', action = 'store_true',
        help = 'Run linclust instead of easy-cluster if lacking precomputed ' \
             + 'groups')
    parser.add_argument('-w', '--window', default = 5, type = int,
        help = 'Max genes +/- for locus window. DEFAULT: 5 (11 gene window)')
    parser.add_argument('-t', '--topology_constraint',
        help = 'Constrain microsynteny topology to species tree')
    parser.add_argument('-c', '--cpus', default = 1, type = int)
    parser.add_argument('-o', '--output')
    args = parser.parse_args()

    execs = ['iqtree']
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
    findExecs(execs, exit = set(execs))

    if not args.output:
        out_dir = mkOutput(os.getcwd() + '/', 'db2microsyntree')
    else:
        out_dir = mkOutput(format_path(args.output), 'db2microsyntree')

    args_dict = {'Database': args.db, 'Focal genes': args.focal_genes,
                 'Precomputed homogroups': homogroups,
                 'Window size': args.window,
                 'Constraint': args.topology_constraint,
                 'CPUS': args.cpus, 'output': out_dir}
    intro(args_dict, 'db2microsyntree', 'Zachary Konkel, Jason Slot')

    wrk_dir = out_dir + 'working/'
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)

    if args.focal_genes:
        with open(format_path(args.focal_genes), 'r') as raw:
            focal_genes = [x.rstrip() for x in raw.read().split() \
                           if x.rstrip()]

    tree_path = out_dir + 'microsynteny.newick'

    db = mtdb(format_path(args.database)).set_index()

    main(db, homogroups, out_dir, wrk_dir, method,
         tree_path, plusminus = args.window, min_cov = 0, min_id = 0.3,
         n50thresh = None, near_single_copy_genes = focal_genes,
         constraint = format_path(args.topology_constraint), verbose = False, 
         cpus = args.cpus)
