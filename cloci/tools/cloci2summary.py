#! /usr/bin/env python3

import os, sys, argparse, multiprocessing as mp
from mycotools.lib.kontools import format_path, collect_files
from mycotools.lib.dbtools import mtdb, primaryDB

def CalcMeanMedian(genesPerClus):

    genesPerClus.sort()
    meanGpC = sum(genesPerClus)/len(genesPerClus)

    if len(genesPerClus) % 2: # it's odd
        medianGpC = int(genesPerClus[round((len(genesPerClus)/2))])
    else:
        medianGpC = sum(
            genesPerClus[int((len(genesPerClus)/2)-1):int((len(genesPerClus)/2)+1)]
            )/2


    return meanGpC, medianGpC

def GrabGenes_cloci(info_path):
    genes, gene2cluster, genesPerClus = [], {}, []
    with open(info_path, 'r') as raw:
        for line in raw:
            d = line.strip().split('\t')
            clus, tGenes = d[0], d[2].split(',')
            genes.extend(tGenes)
            genesPerClus.append(len(tGenes))
            for gene in tGenes:
                gene2cluster[gene] = clus

    if not genesPerClus: # something failed
        return 

    meanGpC, medianGpC = CalcMeanMedian(genesPerClus)

    return set(genes), gene2cluster, meanGpC, medianGpC

def GrabGenes_antismash(gbks):
    genes, gene2cluster, genesPerClus = [], {}, []
    for gbk in gbks:
        clus = os.path.basename(os.path.abspath(gbk))[:-4]
        with open(gbk, 'r') as raw:
            tGenes = []
            for line in raw:
                line = line.lstrip().rstrip()
                if line.startswith('/Alias'):
                    tGene = line[line.find('"'):]
                    gene = tGene.replace('"','')
                    tGenes.append(gene)
                    gene2cluster[gene] = clus
            genesPerClus.append(len(tGenes))
            genes.extend(tGenes)

    if not genesPerClus:
        return
    meanGpC, medianGpC = CalcMeanMedian(genesPerClus)

    return set(genes), gene2cluster, meanGpC, medianGpC

def GrabHits(clociGenes, antismashGenes, gene2cluster):
    intersection = clociGenes.intersection(antismashGenes)

    hits = []
    for gene in list(intersection):
        hits.append(gene2cluster[gene])

    return list(set(hits))

def ome_clus_mngr(ome, cloci_dir, as_dir):
    print(ome, flush = True)
    try:
        o2cGenes, o2cGene2Cluster, o2cMean, o2cMedian = GrabGenes_cloci(
            cloci_dir + 'ome/' + ome + '/info.out'
            )
    except IndexError: # failed
        return

    as_check = False
    if as_dir:
        if os.path.isdir(as_dir + ome):
            as_dir = as_dir + ome
            as_check = True

    asHits, o2cHits, asGene2Cluster, asMean, asMedian = \
        [], [], {}, None, None
    if as_check:
        gbks = collect_files(as_dir, 'gbk')
        gbks = [
            x for x in gbks if not \
            os.path.basename(os.path.abspath(x)).startswith(ome)
            ] # remove full gbk
        if gbks:
            asGenes, asGene2Cluster, asMean, asMedian = GrabGenes_antismash(gbks)
            o2cHits = GrabHits(o2cGenes, asGenes, o2cGene2Cluster)
            asHits = GrabHits(o2cGenes, asGenes, asGene2Cluster)

    stats = tuple([
        tuple(['o2c', tuple([
            len(set(o2cGene2Cluster.values())),
            o2cMean, o2cMedian, len(o2cHits)
            ])]),
        tuple(['as', tuple([
            len(set(asGene2Cluster.values())),
            asMean, asMedian, len(asHits)
            ])])
        ])


    return ome, stats

def write_stat_results(results, out_file):
    data = {
        k: v for k,v in sorted({
            x[0]: x[1] for x in results
            }.items(), key = lambda y: y[0]
            )
        } # sort stats by ome
    header = [
        'ome', 'clus', 'mean_genes_per', 
        'median_genes_per' #,0 'o2cRecoveredByAS', 'asClus', 
#        'asMeanGenesPerClus', 'asMedianGenesPerClus', 'asRecoveredByOG2clus'
        ]
    header_str = '#' + '\t'.join(header) + '\n'
    with open(out_file, 'w') as out:
        out.write(header_str)
        for ome, stats in data.items():
            if not ome: # it failed
                continue
            d = [
                ome, stats[0][1][0], stats[0][1][1], stats[0][1][2]
                ]
#            if stats[1][1][0]: # if there are antismash clusters
 #               d.extend([
  #                  stats[0][1][3]/stats[0][1][0],
   #                 stats[1][1][0], #stats[1][1][1], stats[1][1][2],
#                    stats[1][1][3]/stats[1][1][0]
    #            ])
    #        else:
      #          d.extend(['', '', '', '', ''])
            out.write('\t'.join([
                str(x) for x in d
                ]) + '\n')

def calc_cluster_stats(db, cloci_dir, cpus = 1):
    calc_stats_cmds = []
    for ome in db.set_index():
        if os.path.isdir(cloci_dir + 'ome/' + ome):
            calc_stats_cmds.append(
                [ome, cloci_dir, None]
                )

    with mp.Pool(processes = cpus) as pool:
        res = pool.starmap(ome_clus_mngr, calc_stats_cmds)

    write_stat_results(res, cloci_dir + 'stats.tsv')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Summarize cloci output')
    parser.add_argument('-d', '--db', default = primaryDB())
    parser.add_argument('-c', '--cloci', required = True, help = 'cloci directory')
#    parser.add_argument('-a', '--antismash', help = 'antiSMASH directory')
    parser.add_argument('--cpu', default = 1, type = int)
    args = parser.parse_args()

    db = mtdb(format_path(args.db))
    o2c_dir = format_path(args.cloci)
 #   if args.antismash:
  #      as_dir = format_path(args.antismash)
   # else:
    #    as_dir = None

    calc_cluster_stats(db, o2c_dir, args.cpu)

    sys.exit(0)
