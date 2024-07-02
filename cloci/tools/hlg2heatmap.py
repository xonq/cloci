#! /usr/bin/env python3

import sys
import gzip
import argparse
from itertools import chain
from collections import defaultdict
from mycotools.lib.kontools import format_path

def parse_hlgs(hlg_file):
    hlg2omes = {}
    with gzip.open(hlg_file, 'rt') as raw:
        for line in raw:
            if not line.startswith('#'):
                d = line.rstrip().split()
                hlg = int(d[1])
                omes = d[9].split(',')
                hlg2omes[hlg] = omes
    return hlg2omes

def gen_heatmap(hlg2omes, ref_omes):
    ome2hlgs = defaultdict(set)
    heatmap = defaultdict(list)
    if not ref_omes:
        hlgs = sorted(hlg2omes.keys())
        for hlg, omes in hlg2omes.items():
            for ome in omes:
                ome2hlgs[ome].add(hlg)
    else:
        set_omes = set(ref_omes)
        for hlg, omes in hlg2omes.items():
            for ome in omes:
                if ome in set_omes:
                    ome2hlgs[ome].add(hlg)
        hlgs = sorted(set(chain(*list(ome2hlgs.values()))))
        ome2hlgs = {o: ome2hlgs[o] for o in ref_omes}

    for ome, hlg_set in ome2hlgs.items():
        for hlg in hlgs:
            if hlg in hlg_set:
                heatmap[ome].append(1)
            else:
                heatmap[ome].append(0)
    return heatmap, hlgs
    
def write_output(heatmap, hlgs):
    print('#ome\t' + '\t'.join([str(x) for x in hlgs]) + '\n', flush = True)
    for ome, data in heatmap.items():
        print(f'{ome}\t' + '\t'.join([str(x) for x in data]) + '\n', 
              flush = True)

def main(cloci_dir, omes, gcfs = False):
    if gcfs:
        hlg_file = cloci_dir + 'gcfs.tsv.gz'
    else:
        hlg_file = cloci_dir + 'hlgs.tsv.gz'

    hlg2omes = parse_hlgs(hlg_file)
    return gen_heatmap(hlg2omes, omes)

def cli():
    # cloci output dir
    # gcf option
    # db/ome_list
    parser = argparse.ArgumentParser(description = 'Create table of HLG ' \
           + 'presence-absence')
    parser.add_argument('-c', '--cloci', help = 'CLOCI output dir', 
                        required = True)
    parser.add_argument('-g', '--gcf', help = 'Output GCFs instead of HLGs',
                        action = 'store_true')
    parser.add_argument('-o', '--omes', help = 'MTDB or list of omes to ' \
           + 'extract in order')
    args = parser.parse_args()

    if args.omes:
        with open(format_path(args.omes), 'r') as raw:
            omes = [x.split()[0] for x in raw \
                    if x.rstrip() and not x.startswith('#')]
    else:
        omes = []

    heatmap, hlgs = main(format_path(args.cloci), omes, args.gcf)
    write_output(heatmap, hlgs)

if __name__ == '__main__':
    cli()
    sys.exit(0)

