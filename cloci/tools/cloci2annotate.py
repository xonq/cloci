#! /usr/bin/env python3

import os
import sys
import shutil
import argparse
from tqdm import tqdm
from cloci.lib.output_res import ann_mngr
from mycotools.db2search import compile_hmm_queries
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.lib.kontools import format_path, eprint, mkOutput, findExecs

# NEED output parsable gcfs.tsv.gz for hlg2biofile
# NEED variable hmmsearch args

def input_annotations(ann_file):
    with open(ann_file, 'r') as raw:
        for line in raw:
           if not line.startswith('#'):
               d = line.rstrip().split('\t')
               if len(d) > 1:
                   acc2ann[d[0]] = d[1]
               elif len(d) == 1:
                   acc2ann[d[0]] = d[1]
    return acc2ann

def import_hlgs(ome, hlg_file):
    gene2hg, gene2ann, clus2gene, clus2gcf = {}, {}, {}, {}
    with open(hlg_file, 'r') as raw:
        for line in raw:
            if not line.startswith('#'):
                d = line.rstrip().split('\t')
                c, h, g, f = d[0], d[1].split(','), d[2].split(','), int(d[3])
                if len(d) > 4:
                    a = d[4].split(';')
                    for i, gene in enumerate(g):
                        gene2hg[gene] = h[i]
                        gene2ann[gene] = a[i]
                else:
                    for i, gene in enumerate(g):
                        gene2hg[gene] = h[i]
                        gene2ann[gene] = ''
                clus2gene[c] = g
                clus2gcf[c] = f
    return ome, gene2hg, gene2ann, clus2gene, clus2gcf


def main(hmm_paths, ome_dir, out_dir, db, prefix = 'hlg', 
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
        acc2ann = import_annotations(ann_file)
    else:
        acc2ann = {}

    q_hmm = wrk_dir + 'query.hmm'
    queries = compile_hmm_queries(hmm_paths, q_hmm)


    extract_cmds = []
    for hlg_file in collect_files(ome_dir, 'tsv'):
        if os.path.basename(hlg_file) == f'{prefix}.tsv':
            ome = os.path.basename(os.path.dirname(hlg_file))
            if not os.path.isdir(ome_dir + ome):
                os.mkdir(ome_dir + ome)
            extract_cmds.append((ome, hlg_file))
            shutil.copy(hlg_file, ome_dir + ome + '/old.tsv')
            
    with mp.Pool(processes = cpus) as pool:
        out_res = pool.starmap(import_hlgs, tqdm(extract_cmds, total = len(extract_cmds))

    ome2gene2hg, ome2gene2ann, ome2clus2gene, ome2clus2gcf = \
        {}, {}, {}, {}
    prot_paths = []
    for ome, gene2hg, gene2ann, clus2gene, clus2gcf in out_res:
        ome2gene2hg[ome] = gene2hg
        ome2gene2ann[ome] = gene2ann
        ome2clus2gene[ome] = clus2gene
        ome2clus2gcf[ome] = clus2gcf
        prot_path.append(db[ome]['faa'])

    genes = list(chain(*[list(x.keys()) for x in ome2gene2hg.values()]))

    print('\nAnnotating', flush = True)
    ann_res, failed_omes = ann_mngr(genes, prot_paths, wrk_dir, evalue = 0.01,
                                    threshold = 0.35, cpus = cpus)


    print('\nApplying annotations', flush = True)
    annotation_cmds = []
    for ome, gene2hg in ome2gene2hgs.items():
        annotation_cmds.append((list(gene2hg.keys()), f'{ome_dir}{ome}/old.tsv',
                                f'{ome_dir}{ome}/{os.path.basename(hlg_file)}',
                                ome, gene2hg, ann_res[ome]))

    with mp.Pool(processes = cpus) as pool:
        pool.starmap(annotate_clusters, tqdm(annotation_cmds, 
                                             total = len(annotation_cmds)))
    

def cli():
    parser = argparse.ArgumentParser(
        description = 'Annotates HLGs/GCFs. Will overwrite gcfs if -e is specified and not -o')
    parser.add_argument('-i', '--input', help = 'HMM dir/file', reqired = True)
    parser.add_argument('-c', '--cloci', help = 'CLOCI|CLOCI/ome dir', required = True)
    parser.add_argument('-d', '--mtdb', default = primaryDB())
    parser.add_argument('-a', '--annotations', 
        help = 'File with annotation conversions: <ACC>\\t<ANN>')
    parser.add_argument('-g', '--gcf', 'Extract from GCFs', action = 'store_true')
    parser.add_argument('-e', '--extract', 'Extract annotations', action = 'store_true')
    parser.add_argument('--cpus', type = int, default = 1)

    findExecs(['hmmsearch'], exit = ['hmmsearch'])

    if os.path.basename(os.path.dirname(format_path(args.cloci)[:-1])) == 'ome':
        out_dir = mkOutput(format_path('./'), 'cloci2annotation') 
        ome_dir = format_path(args.cloci)
    else:
        out_dir = mkOutput(format_path(args.cloci), 'cloci2annotation', suffix = False)
        ome_dir = foarmat_path(args.cloci) + 'ome/'

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

    main(hmm_files, ome_dir, out_dir, mtdb(format_path(args.mtdb)), prefix, 
         extract = args.extract, ann_file = format_path(args.annotations), 
         cpus = args.cpus)


if __name__ == '__main__':
    cli()
    sys.exit(0)
