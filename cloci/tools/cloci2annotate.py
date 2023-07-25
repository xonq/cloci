#! /usr/bin/env python3

import os
import sys
import gzip
import shutil
import argparse
import multiprocessing as mp
from tqdm import tqdm
from itertools import chain
from cloci.lib.output_data import ann_mngr, annotate_clusters
from mycotools.db2search import compile_hmm_queries
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.lib.kontools import format_path, eprint, mkOutput, findExecs, \
    collect_files

# NEED variable hmmsearch args
# NEED to output annotations for hlg runs and gcfs runs and reference HLGs from gcf if present

def input_annotations(ann_file):
    # 0 = no change, 1 = and condition, 2 = or condition
    ann_change = False
    acc2ann = {}
    with open(ann_file, 'r') as raw:
        for line in raw:
           if not line.startswith('#'):
               d = line.rstrip().split('\t')
               if len(d) > 1:
                   ann_change = True
                   if ';' not in d[0] and '|' not in d[0]:
                       acc2ann[d[0]] = [None, d[1]]
                   else:
                       ands = d[0].split(';')
                       and_ors = [x.split('|') for x in ands]
                       for a in and_ors:
                           for o in a:
                               if ';' in o or '|' in o:
                                   eprint('\nERROR: multi-level nested annotations ' \
                                        + 'not supported', flush = True)
                                   sys.exit(5)
                               acc2ann[o] = [[set(x) for x in and_ors if x != a],
                                             d[1]]
                  
               elif len(d) == 1:
                   acc2ann[d[0]] = [None, d[0]]
    return acc2ann, ann_change


def apply_annotations(hlg_file, out_file, acc2ann, extract = False):
    data, count, passing, gcfs = [], 0, set(), []
    with open(hlg_file, 'r') as raw:
        for line in raw:
            clus_type = set()
            d = line.rstrip().split('\t')
            if line.startswith('#'):
                data.append(d)
                continue
            if len(d) == 6: # previous annotations
                prev_anns = [x.split('|') if x else [] \
                             for x in d[-2].split(';')]
            else:
                prev_anns = None
            anns = [x.split('|') if x else [] for x in d[-1].split(';')]
            a_set = set(chain(*anns))
            for i0, v in enumerate(anns):
                for i1, ann in enumerate(v):
                    if ann in acc2ann:
                        if acc2ann[ann][0] and acc2ann[ann][1] not in clus_type:
                            trues = 0
                            ands_len = len(acc2ann[ann][0])
                            for etc in acc2ann[ann][0]:
                                if any(x in a_set for x in etc):
                                    trues += 1
                            if trues == ands_len:
                                clus_type.add(acc2ann[ann][1])       
                        elif not acc2ann[ann][0]:
                            clus_type.add(acc2ann[ann][1]) 
                        anns[i0][i1] = acc2ann[ann][1]
            if any(x for x in clus_type):
                d[0] += '_' + '$'.join(sorted(clus_type))
                if '' in clus_type:
                    clus_type.remove('')
                passing.add(count)
                gcfs.append(int(d[3]))
            if prev_anns:
                for i, v in enumerate(anns):
                    if not v:
                        if prev_anns[i]:
                            anns[i] = prev_anns[i]
                d[-1] = ';'.join(['|'.join([x for x in y]) for y in anns])
                del d[-2]
            elif not prev_anns:
                d[-1] = ';'.join(['|'.join([x for x in y]) for y in anns])
            data.append(d)
            count += 1

    with open(out_file, 'w') as out:
        out.write('\t'.join(data[0]) + '\n')
        if extract:
            for i, v in enumerate(data[1:]):
                if i in passing:
                    out.write('\t'.join(v) + '\n')
        else:
            for d in data[1:]:
                out.write('\t'.join(v) + '\n')

    return sorted(set(gcfs))

def import_hlgs(ome, hlg_file):
    gene2hg, gene2ann, clus2gene, clus2gcf = {}, {}, {}, {}
    with open(hlg_file, 'r') as raw:
        for line in raw:
            if not line.startswith('#'):
                d = line.rstrip().split('\t')
                h, g = d[1].split(','), d[2].split(',')
                for i, gene in enumerate(g):
                    gene2hg[gene] = h[i]
    return ome, gene2hg


def main(hmm_paths, cloci_dir, in_dir, out_dir, db, prefix = 'hlg', 
         extract = False, ann_file = None, cpus = 1):

    db = db.set_index()

    print('\nCompiling input', flush = True)
    wrk_dir = out_dir + 'working/'
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)
    ome_dir = out_dir + 'ome/'
    if not os.path.isdir(ome_dir):
        os.mkdir(ome_dir)

    if ann_file:
        acc2ann, ann_change = input_annotations(ann_file)
    else:
        acc2ann, ann_change = {}, False


    q_hmm = wrk_dir + 'query.hmm'
    queries = compile_hmm_queries(hmm_paths, q_hmm)


    extract_cmds = []
    for hlg_file in collect_files(in_dir, 'tsv', recursive = True):
        if os.path.basename(hlg_file) == f'{prefix}.tsv':
            ome = os.path.basename(os.path.dirname(hlg_file))
            if ome not in db:
                continue
            if not os.path.isdir(ome_dir + ome):
                os.mkdir(ome_dir + ome)
            extract_cmds.append((ome, hlg_file))
            shutil.copy(hlg_file, ome_dir + ome + '/old.tsv')
            
    with mp.Pool(processes = cpus) as pool:
        out_res = pool.starmap(import_hlgs, tqdm(extract_cmds, total = len(extract_cmds)))

    ome2gene2hg = {}
    prot_paths = {}
    for ome, gene2hg in out_res:
        ome2gene2hg[ome] = gene2hg
        prot_paths[ome] = db[ome]['faa']

    genes_dict = {ome: [list(genes.keys())] for ome, genes in ome2gene2hg.items()}
    if not genes_dict:
        print('\nNothing to do', flush = True)
        sys.exit(0)

    print('\nAnnotating', flush = True)
    ann_res, failed_omes = ann_mngr(genes_dict, prot_paths, wrk_dir, q_hmm, False,
                                    evalue = 0.01,
                                    threshold = 0.35, cpus = cpus)
    eprint(f'\t\tWARNING: Failures: {",".join(sorted(failed_omes))}', flush = True)
    for ome in failed_omes:
        del ome2gene2hg[ome]


    print('\nApplying annotations', flush = True)
    annotation_cmds = []
    for ome, gene2hg in ome2gene2hg.items():
        annotation_cmds.append((list(gene2hg.keys()), ome, #f'{ome_dir}{ome}/old.tsv',
                                ome_dir, #f'{ome_dir}{ome}/{os.path.basename(hlg_file)}',
                                gene2hg, 'old', ann_res[ome]))

    with mp.Pool(processes = cpus) as pool:
        pool.starmap(annotate_clusters, tqdm(annotation_cmds, 
                                             total = len(annotation_cmds)))

    if extract or ann_change:
        print('\nConverting annotations', flush = True)
        with mp.Pool(processes = cpus) as pool:
            gcfs_res = pool.starmap(apply_annotations, 
                         tqdm(((f'{ome_dir}{ome}/old.tsv',
                               f'{ome_dir}{ome}/{prefix}.tsv',
                               acc2ann, extract) for ome in ome2gene2hg),
                              total = len(ome2gene2hg)))
        gcfs = set(chain(*sorted(gcfs_res)))
        out = ''
        with gzip.open(f'{cloci_dir}hlgs.tsv.gz', 'rt') as raw:
            for line in raw:
                if line.startswith('#'):
                    out += line
                else:
                    d = line.rstrip().split()
                    t_gcf = int(d[1])
                    if t_gcf in gcfs:
                        out += line
        with gzip.open(f'{out_dir}{prefix}s.tsv.gz', 'wt') as out_f:
            out_f.write(out)
    else:
        shutil.copy(f'{cloci_dir}hlgs.tsv.gz', 
                   f'{out_dir}hlgs.tsv.gz')
        

def cli():
    parser = argparse.ArgumentParser(
        description = 'Annotates HLGs/GCFs. Will overwrite gcfs if -e is specified and not -o')
    parser.add_argument('-i', '--input', help = 'HMM dir/file', required = True)
    parser.add_argument('-c', '--cloci', help = 'CLOCI|CLOCI/ome dir', required = True)
    parser.add_argument('-d', '--mtdb', default = primaryDB())
    parser.add_argument('-a', '--annotations', 
        help = 'File with annotation conversions: <ACC>\\t<ANN>; `;` = AND `|` = OR')
    parser.add_argument('-g', '--gcf', help = 'Extract from GCFs', action = 'store_true')
    parser.add_argument('-e', '--extract', help = 'Extract annotations', action = 'store_true')
    parser.add_argument('--cpus', type = int, default = 1)
    args = parser.parse_args()

    findExecs(['hmmsearch'], exit = ['hmmsearch'])

#    if os.path.basename(format_path(args.cloci)[:-1]) == 'ome':
 #       out_dir = mkOutput(format_path('./'), 'cloci2annotation') 
  #      ome_dir = format_path(args.cloci)
   # else:
    out_dir = mkOutput(format_path(args.cloci), 'cloci2annotation', suffix = False)
    ome_dir = format_path(args.cloci) + 'ome/'

    if os.path.isdir(format_path(args.input)):
        hmm_files = collect_files(format_path(args.input), 'hmm')
    elif os.path.isfile(format_path(args.input)):
        hmm_files = [format_path(args.input)]
    else:
        eprint('\nERROR: invalid -i')
        sys.exit(3)

    if args.gcf:
        prefix = 'gcf'
    else:
        prefix = 'hlg'    

    main(hmm_files, format_path(args.cloci), ome_dir, out_dir, 
         mtdb(format_path(args.mtdb)), prefix, 
         extract = args.extract, ann_file = format_path(args.annotations), 
         cpus = args.cpus)


if __name__ == '__main__':
    cli()
    sys.exit(0)
