#! /usr/bin/env python3

import os
import sys
import gzip
import argparse
import multiprocessing as mp
from tqdm import tqdm
from itertools import chain
from collections import defaultdict
from mycotools.acc2gbk import ome_main as acc2gbk
from mycotools.acc2gff import gff_main as acc2gff
from mycotools.lib.dbtools import primaryDB, mtdb
from mycotools.lib.biotools import gff2list, fa2dict, list2gff, dict2fa
from mycotools.lib.kontools import stdin2str, format_path, eprint


def parse_hlgtsv(hlg_tsv):
    hlg2omes = {}
    if hlg_tsv.endswith('.gz'):
        with gzip.open(hlg_tsv, 'rt') as raw:
            for line in raw:
                if not line.startswith('#'):
                    d = line.rstrip().split()
                    hgs_s, hlg_s, hlc_s, tmd_s, p_s, g_s, mi_s, mp_s, ippr_s, o_s = d
                    hlg, omes = int(hlg_s), o_s.split(',')
                    hlg2omes[hlg] = omes

    omes2hlg = defaultdict(list)
    for hlg, omes in hlg2omes.items():
        for ome in omes:
            omes2hlg[ome].append(hlg)

    return hlg2omes, omes2hlg
                     

def parse_ome(ome, ome_file, hlgs = []):
    hlg_set = set(hlgs)

    hlg2genes = defaultdict(list)
    with open(ome_file, 'r') as raw:
        if not hlg_set:
            for line in raw:
                if not line.startswith('#'):
                    d = line.rstrip().split()
                    genes, hlg = d[2].split(','), int(d[3])
                    hlg2genes[hlg].append(genes)
        else:
            for line in raw:
                if not line.startswith('#'):
                    d = line.rstrip().split()
                    genes, hlg = d[2].split(','), int(d[3])
                    if hlg in hlg_set:
                        hlg2genes[hlg].append(genes)
    return ome, hlg2genes

        
def compile_ome_gffs(row, ome, hlg2genes, out_dir, by_ome = False, gbks = False, gffs = False):
    gff = gff2list(row['gff3'])
    hlg2gffs = defaultdict(list)
    for hlg, gene_ls in hlg2genes.items():
        for genes in gene_ls:
            acc2gff_dict = acc2gff(gff, genes)
            comp_gff = []
            for acc, acc_gff in acc2gff_dict.items():
                comp_gff.extend(acc_gff)
            hlg2gffs[hlg].append(comp_gff)

    if gffs:
        for hlg, gffs in hlg2gffs.items():
            for i, gff in enumerate(gffs):
                if by_ome:
                    with open(f'{out_dir}{ome}/{hlg}.{i}.gff3', 'w') as out:
                        out.write(list2gff(gff))
                else:
                    with open(f'{out_dir}{hlg}/{ome}.{i}.gff3', 'w') as out:
                        out.write(list2gff(gff))
            
    if gbks:
        try:
            hlg2gbks = acc2gbk(ome, hlg2gffs, row)
        except KeyError:
            eprint(f'\tERROR: {ome} discrepancy between GFF and FAA', flush = True)
            hlg2gbks = {}
        for hlg, gbks in hlg2gbks.items():
            for i, gbk in enumerate(gbks):
                gbk_data = gbk[list(gbk.keys())[0]]
                if len(gbk) > 1:
                    eprint('\nERROR: multiple contigs for one cluster; huh?', flush = True)
                    sys.exit(9)
                if by_ome:
                    with open(f'{out_dir}{ome}/{hlg}.{i}.gbk', 'w') as out:
                        out.write(gbk_data)
                else:
                    with open(f'{out_dir}{hlg}/{ome}.{i}.gbk', 'w') as out:
                        out.write(gbk_data)
    
def compile_ome_faas(faa_path, ome, hlg2genes, out_dir, by_ome = False):
    faa = fa2dict(faa_path)
    for hlg, gene_ls in hlg2genes.items():
        for i, genes in enumerate(gene_ls):
            hlg_faa = {k: faa[k] for k in genes}
            if by_ome:
                with open(f'{out_dir}{ome}/{hlg}.{i}.faa', 'w') as out:
                    out.write(dict2fa(hlg_faa))
            else:
                with open(f'{out_dir}{hlg}/{ome}.{i}.faa', 'w') as out:
                    out.write(dict2fa(hlg_faa))


def main(db, clo_dir, out_dir, omes = [], hlgs = [], gcf_output = False, 
         make_gbk = False, make_gff = False, make_faa = False, cpus = 1):
    """Generate a dictionary of the output file formats"""

    db = db.set_index('ome')

    ome_dir = clo_dir + 'ome/'
    if gcf_output:
        prefix = 'gcf'
    else:
        prefix = 'hlg'
    if hlgs:
        print(f'\nParsing {prefix.upper()}s', flush = True)
        hlg_tsv = f'{clo_dir}{prefix}s.tsv.gz'
        hlg2omes, omes2hlg = parse_hlgtsv(hlg_tsv)
        with mp.Pool(processes = cpus) as pool:
            ome_res = pool.starmap(parse_ome, tqdm(((ome, f'{ome_dir}{ome}/{prefix}.tsv', hlgs) \
                                               for ome in omes2hlg), total = len(omes2hlg)))
        for hlg in hlgs:
            if not os.path.isdir(out_dir + str(hlg)):
                os.mkdir(out_dir + str(hlg))
    if omes:
        ome_set = set(omes)
        print(f'\nParsing ome {prefix.upper()}s', flush = True)
        ome_dirs = [x for x in os.listdir(ome_dir) \
                    if os.path.isdir(ome_dir + x) and x in db and x in ome_set]
        with mp.Pool(processes = cpus) as pool:
            ome_res = pool.starmap(parse_ome, tqdm(((ome, f'{ome_dir}{ome}/{prefix}.tsv') \
                                               for ome in ome_dirs), total = len(ome_dirs)))
        for ome in omes:
            if not os.path.isdir(out_dir + ome):
                os.mkdir(out_dir + ome)

    if make_gbk or make_gff:
        print(f'\nGenerating GFFs/GBKs', flush = True)
        if cpus > 1:
            with mp.Pool(processes = cpus) as pool:
                pool.starmap(compile_ome_gffs, 
                             tqdm(((db[ome], ome, hlg2genes, out_dir, bool(omes), make_gbk, make_gff) \
                              for ome, hlg2genes in ome_res if hlg2genes), total = len(ome_res)))
        else:
            for ome, hlg2genes in ome_res:
                if hlg2genes:
                    compile_ome_gffs(db[ome], ome, hlg2genes, out_dir, bool(omes), make_gbk, make_gff)
    if make_faa:
        print(f'\nGenerating protein fastas', flush = True)
        with mp.Pool(processes = cpus) as pool:
            pool.starmap(compile_ome_faas,
                         tqdm(((db[ome]['faa'], ome, hlg2genes, out_dir, bool(omes)) \
                              for ome, hlg2genes in ome_res), total = len(ome_res)))
        


def cli():
    parser = argparse.ArgumentParser(description = 'Parse CLOCI output and ' \
           + 'generate biofiles. "-" for stdin; -g specifies GCFs, all HLGs by default')

    in_opt = parser.add_argument_group('Input parameters')
    in_opt.add_argument('-c', '--cloci', help = 'CLOCI output dir', required = True)
    in_opt.add_argument('-g', '--gcf', action = 'store_true', help = 'Generate files for GCFs')
    in_opt.add_argument('-o', '--ome', help = 'Generate files for all ome(s) HLG/GCFs')
    in_opt.add_argument('--hlg', help = 'Generate files for specific HLG/GCF(s)')
    in_opt.add_argument('-d', '--mtdb', help = 'Reference MTDB; DEFAULT: primaryDB',
                        default = primaryDB())

    out_opt = parser.add_argument_group('Output parameters')
    out_opt.add_argument('--gbk', action = 'store_true', 
                         help = 'Generate cluster GenBank files')
    out_opt.add_argument('--gff', action = 'store_true',
                         help = 'Generate cluster GFF files')
    out_opt.add_argument('--faa', action = 'store_true',
                         help = 'Generate cluster amino acid fastas')
    out_opt.add_argument('--out', help = 'Output directory; DEFAULT: -c/biofiles/')

    run_opt = parser.add_argument_group('Runtime parameters')
    run_opt.add_argument('--cpu', type = int, default = 1)
    args = parser.parse_args()
    
    if len([x for x in [args.ome, args.hlg] if x]) > 1:
        eprint('\nERROR: -o and --hlg are incompatible', flush = True)
        sys.exit(4)

    if not args.gcf and not args.ome and not args.hlg:
        print('\nOutputting all HLGs', flush = True)
    elif args.gcf and not args.ome and not args.hlg:
        print('\nOutputting all GCFs', flush = True)


    run_dir = format_path(args.cloci)
    if not run_dir.endswith('/'):
        eprint('\nERROR: -c does not exist', flush = True)
        sys.exit(5)
    elif not args.gbk and not args.gff and not args.faa:
        eprint('\nERROR: --gbk/--gff/--faa required', flush = True)
        sys.exit(6)

    db = mtdb(args.mtdb)

    in_data = None
    for arg in [args.ome, args.hlg]:
        if arg:
            if '-' == arg:
                in_data = stdin2str()
            break

    out_dir = format_path(args.out)
    if not out_dir:
        out_dir = run_dir + 'biofiles/'
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
    elif not out_dir.endswith('/'):
        os.mkdir(out_dir)
        out_dir += '/'

    omes, hlgs = [], []
    if args.ome:
        if in_data:
            omes = in_data.replace('"', '').replace("'", '').split().split(',')
        else:
            omes = [args.ome.replace('"', '').replace('"', '')]
        if args.gcf:
            print('\nOutputting GCFs for omes', flush = True)
        else:
            print('\nOutputting HLGs for omes', flush = True)
    elif args.hlg:
        if in_data:
            hlgs = [int(x) for x in in_data.replace('"', '').replace("'", '').split()]
        else:
            hlgs = [int(x) for x in args.hlg.replace('"', '').replace('"', '').split()]
        if args.gcf:
            print('\nOutputting specified GCFs', flush = True)
        else:
            print('\nOutputting specified HLGs', flush = True)

    main(db, run_dir, out_dir, omes, hlgs, args.gcf, 
         args.gbk, args.gff, args.faa, args.cpu)
    sys.exit(0)


if __name__ == '__main__':
    cli()
