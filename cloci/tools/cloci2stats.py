import os
import re
import sys
import gzip
import argparse
import multiprocessing as mp
from math import log
from collections import defaultdict, Counter
from mycotools.lib.kontools import collect_files, eprint, format_path
from mycotools.lib.biotools import gff3Comps, gff2list
from mycotools.lib.dbtools import mtdb


def assimilate_proxies(prox_file):
    hlg2data = {}
    with gzip.open(prox_file, 'rt') as raw:
        for line in raw:
            if not line.startswith('#'):
                d = line.rstrip().split('\t')
                hlg, lntmd, pds = int(d[1]), float(d[3]), float(d[4])
                gcl, mmi, mmp, csb = float(d[5]), float(d[6]), float(d[7]), float(d[8])
                hlg2data[hlg] = (lntmd, gcl, mmi, mmp, csb, pds,)
    return hlg2data

def count_transcripts(gff_path):
    gff = gff2list(gff_path)
    error = False
    aliases, alia_comp = [], re.compile(gff3Comps()['Alias'])
    for line in gff:
        try:
            alias_init = line['attributes'].find('Alias=')
            alias_end = line['attributes'][alias_init:].find(';')
            if alias_end != -1:
                alias_end += alias_init
                alias = line['attributes'][alias_init + 7:alias_end]
            else:
                alias = line['attributes'][alias_init + 7:]
#            alias = alia_comp.search(line['attributes'])[1]
        except TypeError:
            if not error:
                eprint(f'WARNING: {gff_path} missing alia', flush = True)
            error = True
            continue
        if alias:
            aliases.append(alias)
    alia = set(aliases)
    return len(alia)

def calc_clus_shannon(top_anns):
    len_anns = len(top_anns)
    ann_counts = Counter(top_anns)
    return -sum([v/len_anns * log(v/len_anns) for v in ann_counts.values()])
            

def calc_stats(ome, f, gff_path, alia = None):
    genes_in_clus = []
    hlgs = []
    overall_anns = []
 #   ind_shannons = []
    with open(f, 'r') as raw:
        for line in raw:
            if not line.startswith('#'):
                try:
                    data = line.rstrip().split('\t')
                    clus_id, hgs, genes, t_hlgs = data[0],data[1],data[2],data[3]
                    if len(data) > 4: # there are annotations
                        raw_annotations = data[4]
                        top_anns = [x.split('|')[0] for x in raw_annotations.split(';') if x]
#                        clus_shannon = calc_clus_shannon(top_anns)
                    else:
#                        clus_shannon = None
                        top_anns = []
                except ValueError:
                    raise ValueError(line)
                    sys.exit(4)
                try:
                    hgs = [int(x) for x in hgs.split(',') if x]
                except ValueError:
                    eprint(f'WARNING: no HGs for: {line}')
                genes = genes.split(',')
                t_hlgs = [int(x) for x in t_hlgs.split(',')]
                genes_in_clus.append(len(genes))
                hlgs.extend(t_hlgs)
#                ind_shannons.append(clus_shannon)
                overall_anns.extend(top_anns)

    tot_gic = sum(genes_in_clus)
    mean_gic = tot_gic/len(genes_in_clus)
    genes_in_clus.sort()
    median_gic = genes_in_clus[round(len(genes_in_clus)/2) - 1]
    if overall_anns:
        shannon = calc_clus_shannon(overall_anns)
    else:
        shannon = None

    if not alia:
        alia = count_transcripts(gff_path)
    perc_clus = tot_gic/alia

    return ome, hlgs, mean_gic, median_gic, perc_clus, alia, shannon
                
def summarize_stats(ome_stats, rank, hlg2data, db, tmp_path):
    ome2alia = {}
    tax_dict = defaultdict(lambda: {'tmd': [], 'gcl': [],
                                    'mmi': [], 'mmp': [],
                                    'csb': [], 'pds': [], 'sd': [],
                                    'mean_g': [], 'med_g': [], 'pc': []})

    with open(tmp_path, 'w') as out:
        out.write('tip\ttmd\tgcl\tmmi\tmmp\tcsb\tpds\tpc\tdivers\ttaxon\n')
        if rank == 'ome':
            for ome, hlgs, mean_gic, median_gic, pc, alia, shan in ome_stats:
                prox = [hlg2data[hlg] for hlg in set(hlgs)]
                tax_dict[ome]['tmd'].extend([x[0] for x in prox])
                tax_dict[ome]['gcl'].extend([x[1] for x in prox])
                tax_dict[ome]['mmi'].extend([x[2] for x in prox])
                tax_dict[ome]['mmp'].extend([x[3] for x in prox])
                tax_dict[ome]['csb'].extend([x[4] for x in prox])
                tax_dict[ome]['pds'].extend([x[5] for x in prox])
                tax_dict[ome]['mean_g'].append(mean_gic)
                tax_dict[ome]['med_g'].append(median_gic)
                tax_dict[ome]['pc'].append(pc)
                tax_dict[ome]['sd'].append(shan)
                ome2alia[ome] = alia
                out.write(f'{ome}\t{sum(tax_dict[ome]["tmd"])/len(tax_dict[ome]["tmd"])}\t' \
                        + f'{sum(tax_dict[ome]["gcl"])/len(tax_dict[ome]["gcl"])}\t' \
                        + f'{sum(tax_dict[ome]["mmi"])/len(tax_dict[ome]["mmi"])}\t' \
                        + f'{sum(tax_dict[ome]["mmp"])/len(tax_dict[ome]["mmp"])}\t' \
                        + f'{sum(tax_dict[ome]["csb"])/len(tax_dict[ome]["csb"])}\t' \
                        + f'{sum(tax_dict[ome]["pds"])/len(tax_dict[ome]["pds"])}\t' \
                        + f'{pc}\t{sd}\t\n')
        else:
            for ome, hlgs, mean_gic, median_gic, pc, alia, shan in ome_stats:
                tax = db[ome]['taxonomy'][rank]
                prox = [hlg2data[hlg] for hlg in set(hlgs)]
                tmd = [x[0] for x in prox]
                gcl = [x[1] for x in prox]
                mmi = [x[2] for x in prox]
                mmp = [x[3] for x in prox]
                csb = [x[4] for x in prox]
                pds = [x[5] for x in prox]

                tax_dict[tax]['tmd'].extend(tmd)
                tax_dict[tax]['gcl'].extend(gcl)
                tax_dict[tax]['mmi'].extend(mmi)
                tax_dict[tax]['mmp'].extend(mmp)
                tax_dict[tax]['csb'].extend(csb)
                tax_dict[tax]['pds'].extend(pds)
                tax_dict[tax]['mean_g'].append(mean_gic)
                tax_dict[tax]['med_g'].append(median_gic)
                tax_dict[tax]['pc'].append(pc)
                tax_dict[tax]['sd'].append(shan)
                ome2alia[ome] = alia
                out.write(f'{ome}\t{sum(tmd)/len(tmd)}\t{sum(gcl)/len(gcl)}\t' \
                        + f'{sum(mmi)/len(mmi)}\t{sum(mmp)/len(mmp)}\t' \
                        + f'{sum(csb)/len(csb)}\t{sum(pds)/len(pds)}\t' \
                        + f'{pc}\t{shan}\t{tax}\n')
    return tax_dict, ome2alia

def output_stats(out_file, tax_dict):
    tax_dict = {k: v for k, v in sorted(tax_dict.items(), key = lambda x: x[0])}
    with open(out_file + '.mean.tsv', 'w') as out0, open(out_file + '.median.tsv', 'w') as out1:
        out0.write('#taxon\ttmd\tgcl\tmmi\tmmp\tcsb\tpds\tgenes\tperc_clustered\tdiversity\n')
        out1.write('#taxon\ttmd\tgcl\tmmi\tmmp\tcsb\tpds\tgenes\tperc_clustered\tdiversity\n')
        for tax, data in tax_dict.items():
            tmd = sum(data['tmd'])/len(data['tmd'])
            gcl = sum(data['gcl'])/len(data['gcl'])
            mmi = sum(data['mmi'])/len(data['mmi'])
            mmp = sum(data['mmp'])/len(data['mmp'])
            csb = sum(data['csb'])/len(data['csb'])
            pds = sum(data['pds'])/len(data['pds'])
            genes = sum(data['mean_g'])/len(data['mean_g'])
            pc = sum(data['pc'])/len(data['pc'])
            try:
                sd = sum(data['sd'])/len(data['sd'])
            except TypeError:
                sd = None
            line = f'{tax}\t{tmd}\t{gcl}\t{mmi}\t{mmp}\t{csb}\t{pds}\t{genes}\t{pc}\t{sd}\n'
            out0.write(line)
            
            data['tmd'].sort()
            data['gcl'].sort()
            data['mmi'].sort()
            data['mmp'].sort()
            data['csb'].sort()
            data['pds'].sort()
            data['med_g'].sort()
            data['pc'].sort()
            prox_len = len(data['tmd'])
            tmd = data['tmd'][round(prox_len/2) - 1]
            gcl = data['gcl'][round(prox_len/2) - 1]
            mmi = data['mmi'][round(prox_len/2) - 1]
            mmp = data['mmp'][round(prox_len/2) - 1]
            csb = data['csb'][round(prox_len/2) - 1]
            pds = data['pds'][round(prox_len/2) - 1]
            genes = data['med_g'][round(len(data['med_g'])/2 - 1)]
            pc = data['pc'][round(len(data['pc'])/2 - 1)]
            sd = data['sd'][round(len(data['sd'])/2 - 1)]
            line = f'{tax}\t{tmd}\t{gcl}\t{mmi}\t{mmp}\t{csb}\t{pds}\t{genes}\t{pc}\t{sd}\n'
            out1.write(line)
   

def main(cloci_dir, db, rank, cpus = 1):
    db = db.set_index()

    print('\nCollecting data for omes', flush = True)
    ome_dir = cloci_dir + 'ome/'
    tsvs = collect_files(ome_dir, 'tsv', recursive = True)
    hlgs = [x for x in tsvs if os.path.basename(x) == 'hlg.tsv']
    gcfs = [x for x in tsvs if os.path.basename(x) == 'gcf.tsv']
    omes = [os.path.basename(os.path.dirname(x)) for x in hlgs]
    omes = [x for x in omes if x in db]
    primary_hlg = cloci_dir + 'hlgs.tsv.gz'
    primary_gcf = cloci_dir + 'gcfs.tsv.gz'

    if not os.path.isfile(primary_hlg):
        eprint('\nERROR: `hlgs.tsv.gz` not detected', flush = True)
        sys.exit(1)

    print('\nAssimilating proxy values for HLGs', flush = True)
    hlg2data = assimilate_proxies(primary_hlg)
    
    print('\nCollecting stats by ome', flush = True)
    with mp.Pool(processes = cpus) as pool:
        ome_stats = pool.starmap(calc_stats, ((ome, ome_dir + ome + '/hlg.tsv', db[ome]['gff3']) \
                                              for ome in omes))

    print('\nSummarizing data by ' + rank, flush = True)
    tmp_file = cloci_dir + f'.cloci2summary.hlg.{rank}.tsv'
    hlg_stats, ome2alia = summarize_stats(ome_stats, rank, hlg2data, db, tmp_file)
#    if os.path.isfile(tmp_file):
 #       os.remove(tmp_file)

    out_f = cloci_dir + 'hlg_stats.' + rank

    print('\nWriting HLG stats', flush = True)
    output_stats(out_f, hlg_stats)

    if not os.path.isfile(primary_gcf):
        eprint('\nWARNING: `gcfs.tsv.gz` not detected; ignoring GCFs', flush = True)
 #       gcf2data = []
    else:
        print('\nAssimilating proxy values for GCFs', flush = True)
#        gcf2data = assimilate_proxies(primary_gcf)
        print('\nCollecting stats by ome', flush = True)
        with mp.Pool(processes = cpus) as pool:
            ome_stats = pool.starmap(calc_stats, ((ome, ome_dir + ome + '/gcf.tsv', db[ome]['gff3'], ome2alia[ome]) \
                                                  for ome in omes))
        print('\nSummarizing data', flush = True)
        tmp_file = cloci_dir + f'.cloci2summary.gcf.{rank}.tsv'
        gcf_stats, null = summarize_stats(ome_stats, rank, hlg2data, db, tmp_file)
  #      if os.path.isfile(tmp_file):
   #         os.remove(tmp_file)
        print('\nWriting GCF stats', flush = True)
        out_f = cloci_dir + 'gcf_stats.' + rank 
        output_stats(out_f, gcf_stats)


if __name__ == '__main__':
    ranks = ['ome', 'kingdom', 'phylum', 'subphylum', 'class', 'order', 'family', 'genus', 'species']
    parser = argparse.ArgumentParser(description = 'Summarize CLOCI output for taxonomic rank')
    parser.add_argument('-i', '--input', required = True, help = 'CLOCI input directory')
    parser.add_argument('-r', '--rank', help = f'{ranks}; DEFAULT: ome')
    parser.add_argument('-d', '--mtdb', help = 'MycotoolsDB')
    parser.add_argument('-c', '--cpu', type = int)
    args = parser.parse_args()

    rank = args.rank.lower().rstrip().lstrip()
    if rank not in set(ranks):
        eprint('\nERROR: invalid -r', flush = True)
        sys.exit(3)

    main(format_path(args.input), mtdb(format_path(args.mtdb)), rank, cpus = args.cpu)
    sys.exit(0)
