#! /usr/bin/env python3

import os
import re
import sys
import gzip
import argparse
import multiprocessing as mp
from math import log
from tqdm import tqdm
from collections import defaultdict, Counter
from cloci.lib.output_data import run_hmmsearch, parse_hmm_res
from mycotools.annotationStats import main as annotation_stats
from mycotools.assemblyStats import main as assembly_stats
from mycotools.lib.kontools import collect_files, eprint, format_path, findExecs, mkOutput
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

def calc_shannon(top_anns):
    len_anns = len(top_anns)
    ann_counts = Counter(top_anns)
    return -sum([v/len_anns * log(v/len_anns) for v in ann_counts.values()])
            

def calc_stats(ome, f, gff_path, alia = None):
    genes2hg = {}
    hlgs = []
    overall_anns = []
    if not os.path.isfile(f):
        return None, None, None, None, None, None, None, None, None
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
                    hgs = [int(x) if x else None for x in hgs.split(',')]
                except ValueError:
                    eprint(f'WARNING: no HGs for: {line}')
                    hgs = []
                genes = genes.split(',')
                t_hlgs = [int(x) for x in t_hlgs.split(',')]
                for i, g in enumerate(genes):
                    genes2hg[g] = hgs[i]
                hlgs.extend(t_hlgs)
#                ind_shannons.append(clus_shannon)
                overall_anns.extend(top_anns)

    genes_in_clus = sorted(set(genes2hg.keys()))
    tot_gic = len(genes_in_clus)
    mean_gic = tot_gic/len(genes_in_clus)
    genes_in_clus.sort()
    median_gic = genes_in_clus[round(len(genes_in_clus)/2) - 1]
    if overall_anns:
        alpha = calc_shannon(overall_anns)
    else:
        alpha = None

    if not alia:
        alia = count_transcripts(gff_path)
    perc_clus = tot_gic/alia

    return ome, hlgs, mean_gic, median_gic, \
           perc_clus, tot_gic, alia, alpha, genes2hg
    

def hg2alien(genes, hg_file, ome2tax):
    if not os.path.isfile(hg_file):
        return None
    genes = set(genes)
    g2h = defaultdict(list)
    with open(hg_file, 'r') as raw:
        for line in raw:
            d = line.rstrip().split()
            try:
                q, s, i = d[0], d[1], d[-1]
            except IndexError:
                continue
            if q in genes:
                g2h[q].append((s, float(i)))

    ome2alien = defaultdict(list)
    for q, hits in g2h.items():
        g2h[q] = sorted(hits, key = lambda x: [1], reverse = True)
        self_ome = q[:q.find('_')]
        try:
            self_tax = ome2tax[self_ome]
        except KeyError: # skip omes without taxonomy, remove later
            continue
        in_tax, out_tax = None, None
        while in_tax is None or out_tax is None:
            s, i = hits[0]
            s_tax = s[:s.find('_')]
            if s_tax == self_tax:
                if in_tax is None:
                    in_tax = i
            elif s_tax:
                if out_tax is None:
                    out_tax = i
            del hits[0]
            if not hits:
                break
        if in_tax is None and out_tax is None:
            continue
        elif in_tax is None:
            in_tax = 0
            ome2alien[self_ome].append(0 - log(out_tax))
        elif out_tax is None:
            out_tax = 0
            ome2alien[self_ome].append(log(in_tax) - 0)

    return ome2alien
           
def alien_mngr(genes_in_clus, algn_dir, ome2tax, cpus = 1):
    hg2genes = defaultdict(list)
    for ome, genes2hg in genes_in_clus.items():
        for gene, hg in genes2hg.items():
            if hg is not None:
                hg2genes[hg].append(gene)

#    alien_res = []
 #   for hg, genes in tqdm(hg2genes.items(), total = len(hg2genes)):
  #      alien_res.append(hg2alien(genes, f'{algn_dir}{hg}.out', ome2tax))
    with mp.Pool(processes = cpus) as pool:
        alien_res = pool.starmap(hg2alien, tqdm(((genes, f'{algn_dir}{hg}.out', ome2tax) \
                                     for hg, genes in hg2genes.items()), 
                                     total = len(hg2genes)))

    ome_res = defaultdict(list)
    for ome2alien in alien_res:
        if ome2alien:
            for ome, aliens in ome2alien.items():
                ome_res[ome].extend(aliens)
    
    return ome_res

def summarize_stats(ome_stats, rank, hlg2data, db, tmp_path, ass_stats):
    ome2alia = {}
    tax_dict = defaultdict(lambda: {'tmd': [], 'gcl': [], 'genes': [],
                                    'mmi': [], 'mmp': [], 'gic': [],
                                    'csb': [], 'pds': [], 'sd': [],
                                    'mean_g': [], 'med_g': [], 'pc': [],
                                    'n50': [], 'bp': [], 'mask': []})
    tax2genes2hg = {}
    with open(tmp_path, 'w') as out:
        out.write('tip\ttmd\tgcl\tmmi\tmmp\tcsb\tpds\tpc' \
                + '\tclus_genes\tgenes\tdivers\ttaxon\tn50\tbp\tmask\n')
        if rank == 'ome':
            for ome, hlgs, mean_gic, median_gic, pc, gic, alia, shan, genes in ome_stats:
                if not ome:
                    continue
                tax2genes2hg[ome] = genes
                prox = [hlg2data[hlg] for hlg in set(hlgs) if hlg in hlg2data]
                n50, alen, mask = ass_stats[ome]['n50'], ass_stats[ome]['len'], ass_stats[ome]['mask']
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
                tax_dict[ome]['genes'].append(alia)
                tax_dict[ome]['gic'].append(gic)
                tax_dict[ome]['n50'].append(ass_stats[ome]['n50'])
                tax_dict[ome]['bp'].append(ass_stats[ome]['len'])
                tax_dict[ome]['mask'].append(ass_stats[ome]['mask'])
                ome2alia[ome] = alia
                try:
                    out.write(f'{ome}\t{sum(tax_dict[ome]["tmd"])/len(tax_dict[ome]["tmd"])}\t' \
                            + f'{sum(tax_dict[ome]["gcl"])/len(tax_dict[ome]["gcl"])}\t' \
                            + f'{sum(tax_dict[ome]["mmi"])/len(tax_dict[ome]["mmi"])}\t' \
                            + f'{sum(tax_dict[ome]["mmp"])/len(tax_dict[ome]["mmp"])}\t' \
                            + f'{sum(tax_dict[ome]["csb"])/len(tax_dict[ome]["csb"])}\t' \
                            + f'{sum(tax_dict[ome]["pds"])/len(tax_dict[ome]["pds"])}\t' \
                            + f'{pc}\t{gic}\t{alia}\t{shan}\t{ome}\t{n50}\t{alen}\t{mask}\n')
                except ZeroDivisionError:
                    continue
        else:
            for ome, hlgs, mean_gic, median_gic, pc, gic, alia, shan, genes in ome_stats:
                if not ome:
                    continue
                tax2genes2hg[ome] = genes
                tax = db[ome]['taxonomy'][rank]
                prox = [hlg2data[hlg] for hlg in set(hlgs) if hlg in hlg2data]                        
                tmd = [x[0] for x in prox]
                gcl = [x[1] for x in prox]
                mmi = [x[2] for x in prox]
                mmp = [x[3] for x in prox]
                csb = [x[4] for x in prox]
                pds = [x[5] for x in prox]
                n50, alen, mask = ass_stats[ome]['n50'], ass_stats[ome]['len'], ass_stats[ome]['mask']

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
                tax_dict[tax]['genes'].append(alia)
                tax_dict[tax]['gic'].append(gic)
                tax_dict[tax]['n50'].append(n50)
                tax_dict[tax]['bp'].append(alen)
                tax_dict[tax]['mask'].append(mask)
                ome2alia[ome] = alia
                try:
                    out.write(f'{ome}\t{sum(tmd)/len(tmd)}\t{sum(gcl)/len(gcl)}\t' \
                        + f'{sum(mmi)/len(mmi)}\t{sum(mmp)/len(mmp)}\t' \
                        + f'{sum(csb)/len(csb)}\t{sum(pds)/len(pds)}\t' \
                        + f'{pc}\t{gic}\t{alia}\t{shan}\t{tax}\t{n50}\t{alen}\t{mask}\n')
                except ZeroDivisionError:
                    continue

    return tax_dict, ome2alia, tax2genes2hg

def output_stats(out_file, tax_dict, run_alien = False):
    tax_dict = {k: v for k, v in sorted(tax_dict.items(), key = lambda x: x[0])}
    with open(out_file + '.mean.tsv', 'w') as out0, open(out_file + '.median.tsv', 'w') as out1:
        if run_alien:
            out0.write('#taxon\ttmd\tgcl\tmmi\tmmp\tcsb\tpds' \
                     + '\tgenes_per_clus\tperc_clustered\t' \
                     + 'genes_in_clus\ttot_genes\tdiversity\talienness\t' \
                     + 'n50\tbase_pairs\tmask\n')
            out1.write('#taxon\ttmd\tgcl\tmmi\tmmp\tcsb\tpds' \
                     + '\tgenes_per_clus\tperc_clustered\t' \
                     + 'genes_in_clus\ttot_genes\tdiversity\talienness\t' \
                     + 'n50\tbase_pairs\tmask\n')
        else:
            out0.write('#taxon\ttmd\tgcl\tmmi\tmmp\tcsb\tpds' \
                     + '\tgenes_per_clus\tperc_clustered\t' \
                     + 'genes_in_clus\ttot_genes\tdiversity\t' \
                     + 'n50\tbase_pairs\tmask\n')
            out1.write('#taxon\ttmd\tgcl\tmmi\tmmp\tcsb\tpds\t' \
                     + 'genes_per_clus\tperc_clustered\t' \
                     + 'genes_in_clus\ttot_genes\tdiversity\t' \
                     + 'n50\tbase_pairs\tmask\n')
        for tax, data in tax_dict.items():
            if len(data['tmd']) == 0:
                eprint(f'\tWARNING: {tax} has no results', flush = True)
                continue
            tmd = sum(data['tmd'])/len(data['tmd'])
            gcl = sum(data['gcl'])/len(data['gcl'])
            mmi = sum(data['mmi'])/len(data['mmi'])
            mmp = sum(data['mmp'])/len(data['mmp'])
            csb = sum(data['csb'])/len(data['csb'])
            pds = sum(data['pds'])/len(data['pds'])
            genes = sum(data['mean_g'])/len(data['mean_g'])
            tot_genes = sum(data['genes'])/len(data['genes'])
            tot_cg = sum(data['gic'])/len(data['gic'])
            pc = sum(data['pc'])/len(data['pc'])
            n50 = sum(data['n50'])/len(data['n50'])
            bp = sum(data['bp'])/len(data['bp'])
            try:
                mask = sum(data['mask'])/len(data['mask'])
            except TypeError:
                mask = None
            try:
                sd = sum(data['sd'])/len(data['sd'])
            except TypeError:
                sd = None
            if run_alien:
                alien = sum(data['alien'])/len(data['alien'])
                line = f'{tax}\t{tmd}\t{gcl}\t{mmi}\t{mmp}\t{csb}\t' \
                     + f'{pds}\t{genes}\t{pc}\t{tot_cg}\t{tot_genes}\t' \
                     + f'{sd}\t{alien}\t{n50}\t{bp}\t{mask}\n'
            else:
                line = f'{tax}\t{tmd}\t{gcl}\t{mmi}\t{mmp}\t{csb}\t' \
                     + f'{pds}\t{genes}\t{pc}\t{tot_cg}\t{tot_genes}\t' \
                     + f'{sd}\t{n50}\t{bp}\t{mask}\n'
            out0.write(line)
            
            data['tmd'].sort()
            data['gcl'].sort()
            data['mmi'].sort()
            data['mmp'].sort()
            data['csb'].sort()
            data['pds'].sort()
            data['med_g'].sort()
            data['pc'].sort()
            data['genes'].sort()
            data['gic'].sort()
            data['n50'].sort()
            data['bp'].sort()
            data['mask'].sort()
            prox_len = len(data['tmd'])
            tmd = data['tmd'][round(prox_len/2) - 1]
            gcl = data['gcl'][round(prox_len/2) - 1]
            mmi = data['mmi'][round(prox_len/2) - 1]
            mmp = data['mmp'][round(prox_len/2) - 1]
            csb = data['csb'][round(prox_len/2) - 1]
            pds = data['pds'][round(prox_len/2) - 1]
            genes = data['med_g'][round(len(data['med_g'])/2 - 1)]
            tot_genes = data['genes'][round(len(data['genes'])/2 - 1)]
            tot_cg = data['gic'][round(len(data['gic'])/2 - 1)]
            pc = data['pc'][round(len(data['pc'])/2 - 1)]
            sd = data['sd'][round(len(data['sd'])/2 - 1)]
            n50 = data['n50'][round(len(data['n50'])/2 - 1)]
            bp = data['bp'][round(len(data['bp'])/2 - 1)]
            data['mask'] = [x for x in data['mask'] if x] # IGNORES TRUE 0% MASKED
            try:
                mask = data['mask'][round(len(data['mask'])/2) - 1]
            except IndexError:
                mask = None

            if run_alien:
                data['alien'].sort()
                alien = data['alien'][round(len(data['alien'])/2) - 1]
                line = f'{tax}\t{tmd}\t{gcl}\t{mmi}\t{mmp}\t{csb}\t' \
                     + f'{pds}\t{genes}\t{pc}\t{tot_cg}\t{tot_genes}\t' \
                     + f'{sd}\t{alien}\t{n50}\t{bp}\t{mask}\n'
            else:
                line = f'{tax}\t{tmd}\t{gcl}\t{mmi}\t{mmp}\t{csb}' \
                     + f'\t{pds}\t{genes}\t{pc}\t{tot_cg}\t{tot_genes}\t{sd}' \
                     + f'\t{n50}\t{bp}\t{mask}\n'
            out1.write(line)


def calc_gamma_diversity(ome, ann_dir):
    """Calculate genome-wide Shannon diversity referencing the top hit,
    excludes genes with missing annotations"""
    null, ann_res = parse_hmm_res(ome, ann_dir, 0.001, 0.5, 1)
    gamma_diversity = calc_shannon(list(ann_res.values()))
    return ome, gamma_diversity


def add_alien_data(db, stats_file, stats_dict, rank, ome2alien):
    if rank != 'ome':
        tax2alien = defaultdict(list)
        for ome, aliens in ome2alien.items():
            tax = db[ome]['taxonomy'][rank]
            tax2alien[tax].extend(aliens)
    else:
        tax2alien = ome2alien

    for tax, aliens in tax2alien.items():
        stats_dict[tax]['alien'] = aliens

    with open(stats_file + '.tmp', 'w') as out:
        with open(stats_file, 'r') as raw:
            for line in raw:
                if line.startswith('#'):
                    out.write(line.rstrip() + '\talienness\n')
                else:
                    taxon = line.split()[0]
                    if taxon in tax2alien:
                        out.write(line.rstrip() \
                        + f'\t{sum(tax2alien[taxon])/len(tax2alien[taxon])}\n')
                    else:
                        out.write(line.rstrip() + '\tna\n')
    os.rename(stats_file + '.tmp', stats_file)

    return stats_dict


def main(out_dir, ome_dir, db, rank, gamma = False, ann_dir = None, 
         pfam = False, alien = False, gcf_only = False, cpus = 1):

    db = db.set_index()

    if gamma:
        if not ann_dir:
            ann_dir = out_dir + 'genome_ann/'
        if not os.path.isdir(ann_dir):
            os.mkdir(ann_dir)
        annotations = collect_files(ann_dir, 'out')
        omes = [os.path.basename(x[:-4]) for x in annotations]
        to_run = sorted(set(db.keys()).difference(set(omes)))
        if to_run:
            print(f'\nAnnotating {len(to_run)} full genomes')
            failed_omes = run_hmmsearch(pfam, to_run, ann_dir, cpus, db = db)
            passing_omes = sorted(set(db.keys()).difference(set(failed_omes)))
            
            ann_res = compile_hmm_res(ann_dir, passing_omes, cpus = cpus, 
                                      max_hits = 1)
        

    print('\nPreparing assembly statistics', flush = True)
    assembly_stats('', log_path = out_dir + 'assembly_stats.tsv', cpus = cpus, db = db)

    with open(out_dir + 'assembly_stats.tsv', 'r') as raw:
        ome2ass_stats = {}
        for line in raw:
            if not line.startswith('#'):  
                d = line.rstrip().split()
                ome, n50, mask, ass_len = d[0], d[1], d[-1], d[-5]
                ome2ass_stats[ome] = {'n50': int(n50), 'mask': float(mask),
                                      'len': int(ass_len)}



    print('\nCollecting data for omes', flush = True)
    tsvs_p = collect_files(ome_dir, 'tsv', recursive = True)
    tsvs = [x for x in tsvs_p if os.path.basename(os.path.dirname(x)) in db]
    gcfs = [x for x in tsvs if os.path.basename(x) == 'gcf.tsv']
    omes = [os.path.basename(os.path.dirname(x)) for x in gcfs]
    if rank != 'ome':
        ome2tax = {o: db[o]['taxonomy'][rank] for o in omes}
    else:
        ome2tax = {o: o for o in omes}
    primary_gcf = ome_dir + '../gcfs.tsv.gz'

    if not gcf_only:
        hlgs = [x for x in tsvs if os.path.basename(x) == 'hlg.tsv']
        omes = [os.path.basename(os.path.dirname(x)) for x in hlgs]
        primary_hlg = out_dir + '../hlgs.tsv.gz'
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
        tmp_file = out_dir + f'.cloci2summary.hlg.{rank}.tsv'
        hlg_stats, ome2alia, genes_in_clus = summarize_stats(ome_stats, rank, hlg2data, 
                                                             db, tmp_file, ome2ass_stats)
    #    if os.path.isfile(tmp_file):
     #       os.remove(tmp_file)
        if alien:
            print('\nCalculating HLG alienness', flush = True)
            ome2alien = alien_mngr(genes_in_clus, alien, ome2tax, cpus = cpus)
            hlg_stats = add_alien_data(db, tmp_file, hlg_stats, rank, ome2alien)
        out_f = out_dir + 'hlg_stats.' + rank
    
        print('\nWriting HLG stats', flush = True)
        output_stats(out_f, hlg_stats, alien)
    else:
        hlg2data = assimilate_proxies(primary_gcf)
        ome2alia = None


    
    if not os.path.isfile(primary_gcf):
        eprint('\nWARNING: `gcfs.tsv.gz` not detected; ignoring GCFs', flush = True)
 #       gcf2data = []
    else:
        print('\nAssimilating proxy values for GCFs', flush = True)
#        gcf2data = assimilate_proxies(primary_gcf)
        print('\nCollecting stats by ome', flush = True)
        with mp.Pool(processes = cpus) as pool:
            if ome2alia:
                ome_stats = pool.starmap(calc_stats, tqdm(((ome, f'{ome_dir}{ome}/gcf.tsv', 
                                                       db[ome]['gff3'], ome2alia[ome]) \
                                                      for ome in omes), total = len(omes)))
            else:
                ome_stats = pool.starmap(calc_stats, tqdm(((ome, f'{ome_dir}{ome}/gcf.tsv',
                                                       db[ome]['gff3']) for ome in omes),
                                                       total = len(omes)))
        print('\nSummarizing data', flush = True)
        tmp_file = out_dir + f'.cloci2summary.gcf.{rank}.tsv'
        gcf_stats, null, genes_in_clus = summarize_stats(ome_stats, rank, hlg2data, 
                                                         db, tmp_file, ome2ass_stats)
        if not ome2alia:
            ome2alia = null
        if alien:
            ome2alien = alien_mngr(genes_in_clus, alien, ome2tax, cpus = cpus)
            gcf_stats = add_alien_data(db, tmp_file, gcf_stats, rank, ome2alien)
  #      if os.path.isfile(tmp_file):
   #         os.remove(tmp_file)
        print('\nWriting GCF stats', flush = True)
        out_f = out_dir + 'gcf_stats.' + rank 
        output_stats(out_f, gcf_stats, alien)


def cli():
    ranks = ['ome', 'kingdom', 'phylum', 'subphylum', 'class', 'order', 'family', 'genus', 'species']
    parser = argparse.ArgumentParser(description = 'Summarize CLOCI output for taxonomic rank')
    parser.add_argument('-c', '--cloci', required = True, help = 'CLOCI|CLOCI/ome output')
    parser.add_argument('-r', '--rank', help = f'{ranks}; DEFAULT: ome')
    parser.add_argument('-d', '--mtdb', help = 'MycotoolsDB')
    parser.add_argument('-a', '--alien',
        help = 'Directory of HG alignments to quantify Alien Index')
    parser.add_argument('-g', '--gcf', action = 'store_true',
        help = 'Only run calculations for GCFs')
#    parser.add_argument('-g', '--gamma', action = 'store_true',
   #     help = 'Calculate gamma and beta diversity')
 #   parser.add_argument('-a', '--annotations', 
   #     help = '[-g] Directory of tbl-formatted Pfam annotations for gamma diversity, labeled <ome>.out')
  #  parser.add_argument('-p', '--pfam', help = '[-g] Pfam.hmm path')
    parser.add_argument('--cpus', type = int, default = 1)
    args = parser.parse_args()

    rank = args.rank.lower().rstrip().lstrip()
    if rank not in set(ranks):
        eprint('\nERROR: invalid -r', flush = True)
        sys.exit(3)

#    if args.gamma:
 #       findExecs(['hmmsearch', exit = set('hmmsearch')])

    if os.path.basename(os.path.dirname(format_path(args.cloci))) == 'ome':
        out_dir = mkOutput(format_path('./'), 'cloci2stats')
        ome_dir = format_path(args.cloci)
    else:
        out_dir = mkOutput(format_path(args.cloci), 'cloci2stats', suffix = None)
        ome_dir = format_path(args.cloci) + 'ome/'

    main(out_dir, ome_dir, mtdb(format_path(args.mtdb)), rank, 
         alien = format_path(args.alien), gcf_only = args.gcf, cpus = args.cpus)
  #       gamma = False, ann_dir = format_path(args.annotations), 
   #      pfam = format_path(args.pfam),    
    sys.exit(0)


if __name__ == '__main__':
    cli()
