import os
import sys
import pickle
import datetime
import subprocess
import multiprocessing as mp
from tqdm import tqdm
from itertools import chain
from collections import defaultdict
from mycotools.acc2fa import dbmain as acc2fa
from mycotools.lib.biotools import gff2list
from mycotools.lib.kontools import multisub, tardir, collect_files, checkdir, eprint
from orthocluster.orthocluster.lib.input_parsing import compileCDS2, hg_fa_mngr
from orthocluster.orthocluster.lib.treecalcs import calc_tmd

# NEED an option to import OrthoFinder pairwise alignments
	# not recommended, need higher max target sequences than OrthoFinder 
def parse_algn_wpos(hgx_dir, hg, hgx_genes):
    qlen, hit_len = len(hgx_genes), 0
    with open(f'{hgx_dir}{hg}.out', 'r') as raw:
        # open the blast results
        gene2algn = defaultdict(list)
        for line in raw:
            d = line.rstrip().split('\t')
            # need to adjust to just use d when safe
            q, s, p, i = d[0], d[1], d[2], d[3]
            if q in hgx_genes:
                gene2algn[q].append((s, float(i), float(p)))
                hit_len = len(gene2algn)
            elif hit_len == qlen:
                break
    return gene2algn

def parse_algn_wopos(hgx_dir, hg, hgx_genes):
    qlen, hit_len = len(hgx_genes), 0
    with open(f'{hgx_dir}{hg}.out', 'r') as raw:
        # open the blast results
        gene2algn = defaultdict(list)
        for line in raw:
            d = line.rstrip().split('\t')
            # need to adjust to just use d when safe
            q, s, i = d[0], d[1], d[2]
            if q in hgx_genes:
                gene2algn[q].append((s, float(i)))
                hit_len = len(gene2algn)
            elif hit_len == qlen:
                break
    return gene2algn

def get_top_hits_wpos(gene, gcf, hgx_genes, hits, gene2algn, res, ome):
    ome_hits = set()
    while gene2algn[gene]:
        sbj_gene, sbj_id, sbj_pos = gene2algn[gene][0]
        sbj_ome = sbj_gene[:sbj_gene.find('_')]
        # if the subject is also in the checked queries and it is not
        # a recent paralog of the query species
        if sbj_gene in hgx_genes and sbj_ome != ome:
            # grab the hit similarity
            res[gcf][gene][1][sbj_gene] = (sbj_id, sbj_pos,)
            hits.add(sbj_gene)
            ome_hits.add(sbj_ome)
            if not hgx_genes.difference(hits):
                return hits, res, gene2algn
        # add one to the failed count
        else: 
            res[gcf][gene][0] += 1 
        # proceed to the next query
        try:
            gene2algn[gene].pop(0)
        except IndexError:
            break
    
    hgx_omes = set(x[:x.find('_')] for x in list(hgx_genes))
    if hgx_genes.difference(hits) and hgx_omes.difference(ome_hits):
        # missing gene and no paralog
        res[gcf][gene][2] = True 

    return hits, res, gene2algn
        
def get_top_hits_wopos(gene, gcf, hgx_genes, hits, gene2algn, res, ome):
    ome_hits = set()
    while gene2algn[gene]:
        sbj_gene, sbj_id = gene2algn[gene][0]
        sbj_ome = sbj_gene[:sbj_gene.find('_')]
        # if the subject is also in the checked queries and it is not
        # a recent paralog of the query species
        if sbj_gene in hgx_genes and sbj_ome != ome:
            # grab the hit similarity
            res[gcf][gene][1][sbj_gene] = sbj_id
            hits.add(sbj_gene)
            ome_hits.add(sbj_ome)

            if not hgx_genes.difference(hits):
                return hits, res, gene2algn   
        # add one to the failed count
        else: 
            res[gcf][gene][0] += 1 
        # proceed to the next query
        try:
            gene2algn[gene].pop(0)
        except IndexError:
            break
    
    hgx_omes = set(x[:x.find('_')] for x in list(hgx_genes))
    if hgx_genes.difference(hits) and hgx_omes.difference(ome_hits):
        # missing gene and no paralog
        # biased against convergent assimilation
        res[gcf][gene][2] = True 

    return hits, res, gene2algn   

def sim_calc_wpos(res):
    omeScores = defaultdict(dict)
    ids_y_pos = defaultdict(dict)
    for hg in res:
        for ome, data in res[hg].items():
            t_gcc, id_y_pos_data = data
            omeScores[ome][hg] = t_gcc
            ids = [x[0] for x in id_y_pos_data.values()]
            pos = [x[1] for x in id_y_pos_data.values()]
            if ids:
                ids_y_pos[hg][ome] = (min(ids), min(pos),)
            else:
                ids_y_pos[hg][ome] = (0, 0,)
    return omeScores, ids_y_pos

def sim_calc_wopos(res):
    omeScores = defaultdict(dict)
    ids_dict = defaultdict(dict)
    for hg in res:
        for ome, data in res[hg].items():
            t_gcc, id_data = data
            omeScores[ome][hg] = t_gcc
            ids = list(id_data.values())
            if ids:
                ids_dict[hg][ome] = min(ids)
            else:
                ids_dict[hg][ome] = 0
    return omeScores, ids_dict

def mms_calcs_wpos(ids_y_pos):
    # average minimum identity / minimum positive; for each gene homolog group
    gcf_min_id, gcf_min_pos, total = 0, 0, 0
    for ome_dict in ids_y_pos.values():
        gcf_min_id += sum([x[0] for x in ome_dict.values()])
        gcf_min_pos += sum([x[1] for x in ome_dict.values()])
        total += len(ome_dict)
    try:
        gcf_min_id /= total
        gcf_min_pos /= total
    except ZeroDivisionError:
        gcf_min_id, gcf_min_pos = 0, 0
    return gcf_min_id, gcf_min_pos

def mms_calcs_wopos(ids_dict):
    # average minimum identity; for each gene homolog group
    gcf_min_id, total = 0, 0
    for ome_dict in ids_dict.values():
        gcf_min_id += sum(ome_dict.values())
        total += len(ome_dict)
    try:
        gcf_min_id /= total
    except ZeroDivisionError:
        gcf_min_id = 0
    return gcf_min_id, None

def hg_sim_calc_wpos(res, min_id):
    output_data = {}
    for gcf, gene_info in res.items():
        output_data[gcf] = {}
        for gene, data in gene_info.items():
            ome = gene[:gene.find('_')]
            ome_gcc = data[0]
            if not data[2]:
                try:
                    ome_mi = min([v[0] for v in data[1].values()])
                    ome_mp = min([v[1] for v in data[1].values()])
                except ValueError:
                    ome_mi, ome_mp = 0, 0
            else:
                ome_mi, ome_mp = min_id, min_id
            if ome not in output_data[gcf]:
                output_data[gcf][ome] = (ome_gcc, ome_mi, ome_mp)
            # if there's a paralog, take the best gcc
            else:
                if ome_gcc > output_data[gcf][ome][0]:
                    output_data[gcf][ome] = (ome_gcc, ome_mi, ome_mp)
    return output_data

def hg_sim_calc_wopos(res, min_id):
    output_data = {}
    for gcf, gene_info in res.items():
        output_data[gcf] = {}
        for gene, data in gene_info.items():
            ome = gene[:gene.find('_')]
            ome_gcc = data[0]
            if not data[2]:
                try:
                    ome_mi = min(data[1].values())
                except ValueError:
                    ome_mi = 0
            else:
                ome_mi = min_id
            if ome not in output_data[gcf]:
                output_data[gcf][ome] = (ome_gcc, ome_mi)
            # if there's a paralog, take the best gcc
            else:
                if ome_gcc > output_data[gcf][ome][0]:
                    output_data[gcf][ome] = (ome_gcc, ome_mi)
    return output_data



def hg_parse_and_calc(hg, hg_dict, hgx_dir, phylo, ome2i, min_id = 30, pos_data = True):

    if pos_data:
       parse_func = parse_algn_wpos
       top_hits_func = get_top_hits_wpos
       hg_sim_func = hg_sim_calc_wpos
    else:
       parse_func = parse_algn_wopos
       top_hits_func = get_top_hits_wopos
       hg_sim_func = hg_sim_calc_wopos

    gene2gcf = {}
    for gcf, genes in hg_dict.items():
        # ignore singletons
        if len(genes) > 1:
            for gene in genes:
                gene2gcf[gene] = gcf

    hg_genes = set(gene2gcf.keys())

    try:
        gene2algn = parse_func(hgx_dir, hg, hg_genes)
    except IndexError:
        if not pos_data:
            eprint(f'\tERROR: malformatted alignment: {hg}', flush = True)
            sys.exit(123)
        else:
            parse_func = parse_algn_wopos
            top_hits_func = get_top_hits_wopos
            hg_sim_func = hg_sim_calc_wopos
        try:
            gene2algn = parse_func(hgx_dir, hg, hg_genes)
        except IndexError:
            eprint(f'\tERROR: malformatted alignment {hg}', flush = True)
        return None, None, None
    except FileNotFoundError:
        return None, None, None

    # sort each query preferably by percent positives
    try:
        gene2algn = {
            k: sorted(v, key = lambda x: x[2], reverse = True) \
            for k,v in gene2algn.items()
            } 
    # otherwise run it by positives
    except IndexError:
        gene2algn = {k: sorted(v, key = lambda x: x[1], reverse = True) \
                     for k, v in gene2algn.items()}
    res = defaultdict(dict)
    todel_genes = defaultdict(lambda: defaultdict(list))

    gcf2considered = {}
    for gcf, genes in hg_dict.items():
        failed_genes = True
        gene_set = set(genes)

        while failed_genes:
            failed_genes = []
            # disregard the self from the length
            hgx_gene_len = len(gene_set) - 1
            # for the weighted average, we can start at this step
            # only consider unique omes because we are only taking the max
            considered_omes = len(set(x[:x.find('_')] for x in genes))
            gcf2considered[gcf] = considered_omes
            for gene in genes:
                ome = gene[:gene.find('_')] # identify the ome
                if gene not in res[gcf]:
                    res[gcf][gene] = [0, {}, False]
        
                # to identify if the subject is in the family of omes
                hits = {gene}
                # while all cluster homologs aren't accounted for and there remain hits
                if gene not in gene2algn:
                    continue
                hits, res, gene2algn = top_hits_func(gene, gcf, gene_set, hits, 
                                                     gene2algn, res, ome)
        
                # if there are missing hits (which means we are disregarding
                # missing alignments)
                if hgx_gene_len - len(res[gcf][gene][1]) > 0:
                    # then get the percent of observed hits that are valid
                    hit_len = res[gcf][gene][0] + len(res[gcf][gene][1])
                    try:
                        res[gcf][gene][0] = hit_len * considered_omes \
                                          / (hit_len + res[gcf][gene][0])
                    except ZeroDivisionError:
                        # no valid hits, the gene should not be called in the cluster
                        # typically happens from extension module
                        todel_genes[gcf][ome].append(gene)
                        failed_genes.append(gene)
                else:
                    # the GCC for this HG and ome is the total of the genes in this
                    # HG (disregarding self) divided by the total + the failed genes
                    res[gcf][gene][0] = hgx_gene_len * considered_omes \
                                      / (hgx_gene_len + res[gcf][gene][0])
            # have to rerun the calculation if there are failed genes
            for failed_gene in failed_genes:
                gene_set.remove(failed_gene)
            genes = list(gene_set)
    
    return hg, hg_sim_func(res, min_id), \
           {gcf: dict(ome_dict) for gcf, ome_dict in todel_genes.items()}, \
           gcf2considered


def gcc_mngr(
    hgs, omes, hgx_dir, hgx2loc,
    db, gene2hg, clusplusminus, hg2gene, gcfs,
    gcf_omes, gcf_hgxs, ome2i, min_id, phylo,
    d2gcc, d2id_, d2pos, cpus = 1
    ):

    i2ome = {v: k for k, v in ome2i.items()}

    clus_hgs = {}
    for fam, loci in gcfs.items():
        gcf_hgx = gcf_hgxs[fam]
        omesc = gcf_omes[fam]
        if gcf_hgx not in d2id_:
            clus_hgs[fam] = [gcf_hgx, defaultdict(list)]
    # are some gcf_hgx HGs relics of the original hgx and thus only contain one gene?
    # would lower evo_conco scores
            gcf_hgx_set = set(gcf_hgx)
            for loc in loci:
                for gene in loc:
                    try:
                        hg = gene2hg[gene]
                    except KeyError:
                        continue
                    if hg in gcf_hgx_set:
                        clus_hgs[fam][1][hg].append(gene)
        elif omesc not in d2id_[gcf_hgx]:
            clus_hgs[fam] = [gcf_hgx, defaultdict(list)]
    # are some gcf_hgx HGs relics of the original hgx and thus only contain one gene?
    # would lower evo_conco scores
            gcf_hgx_set = set(gcf_hgx)
            for loc in loci:
                for gene in loc:
                    try:
                        hg = gene2hg[gene]
                    except KeyError:
                        continue
                    if hg in gcf_hgx_set:
                        clus_hgs[fam][1][hg].append(gene)

            # {fam: [hgx, {hg: set(seqs)}}

    clus_hgs = {
        fam: [d[0], {hg: set(v) for hg, v in d[1].items() if len(set(v)) > 1}] \
        for fam, d in clus_hgs.items()
        } # make sets from it
    hg2genes = defaultdict(dict)
    for gcf, hgs in clus_hgs.items():
        for hg, gene_set in hgs[1].items():
            hg2genes[hg][gcf] = tuple(sorted(gene_set))

    print('\tParsing alignments', flush = True)
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        hg_results = pool.starmap(hg_parse_and_calc,
                              tqdm(((hg, hg_dict, hgx_dir, phylo, ome2i, min_id) \
                               for hg, hg_dict in hg2genes.items()),
                               total = len(hg2genes)))
        pool.close()
        pool.join()

    print('\tQuantifying', flush = True)
    gcf_res = defaultdict(dict)
    deleted = 0
    gcf2tot_considered = defaultdict(int)
    gcfs = {gcf: list(locs) for gcf, locs in gcfs.items()}
    for hg, hg_res, todel_genes, gcf2considered in hg_results:
        todel = defaultdict(list)
        if hg_res:
            for gcf, considered in gcf2considered.items():
                gcf2tot_considered[gcf] += considered
            for gcf, ome_dict in hg_res.items():
                gcf_res[gcf][hg] = []
                for ome, res in ome_dict.items():
                    gcf_res[gcf][hg].append(res)
            for gcf, ome_fails in todel_genes.items():
                for ome, failed_genes in ome_fails.items():
                    fails = set(failed_genes)
                    for i, loc in enumerate(gcfs[gcf]):
                        if loc[0][:loc[0].find('_')] == ome:
                            # if all/all but 1 of the genes failed, then delete it
                            if len([x for x in loc if x in fails]) >= len(loc) - 1:
                                todel[gcf].append(i)
            for gcf, iz in todel.items():
                deleted += len(iz)
                for i in reversed(iz):
                    del gcfs[gcf][i]
    print(f'\t{deleted} clusters removed due to lacking homology', flush = True)
            
    hgx2omes2gcc = defaultdict(dict)
    hgx2omes2id = defaultdict(dict)
    hgx2omes2pos = defaultdict(dict)
    for gcf, hg_dict in tqdm(gcf_res.items(), total = len(gcf_res)):
        hgx, omes = gcf_hgxs[gcf], gcf_omes[gcf]
        gccs, ids, poss = [], [], []
        for hg, res in hg_dict.items():
            gccs.append(sum([x[0] for x in res])/len(res))
            ids.append(sum([x[1] for x in res])/len(res))
            try:
                poss.append(sum([x[2] for x in res])/len(res))
            except IndexError:
                poss.append(None)
            except TypeError:
                poss.append(None)
#        gcc = sum(gccs)/(len(gccs) * gcf2tot_considered[gcf])
        gcc = sum(gccs)/gcf2tot_considered[gcf]
        id_ = sum(ids)/len(ids)
        try:
            pos = sum(poss)/len(poss)
        except ValueError:
            pos = None
        except TypeError:
            pos = None

        hgx2omes2gcc[hgx][omes] = gcc
        hgx2omes2id[hgx][omes] = id_
        if pos is not None:
            hgx2omes2pos[hgx][omes] = pos

    hgx2omes2pos = dict(hgx2omes2pos)
    return {**d2gcc, **hgx2omes2gcc}, \
        {**d2id_, **hgx2omes2id}, \
        {**d2pos, **hgx2omes2pos}, \
        {g: tuple(l) for g, l in gcfs.items()}


def find_missing_algns(hgs, hgx_dir):
    alnd_hgs = set([int(os.path.basename(x[:-4])) \
                    for x in collect_files(hgx_dir, 'out')])
    missing_alns = set(hgs).difference(alnd_hgs)
    return missing_alns


def parse_failures(fail_file):
    if os.path.isfile(fail_file):
        with open(fail_file, 'r') as raw:
            failed_hgs = [int(x.rstrip()) for x in raw]
    else:
        failed_hgs = []
    return failed_hgs


def prep_blast_cmds(db, hgs, hg_dir, hgx_dir, 
                    minid = 30, algorithm = 'diamond',
                    sensitivity = '', hg2gene = None,
                    cpus = 1, rerun = False):

    missing_alns = find_missing_algns(hgs, hgx_dir)

    if not rerun:
        failed_hgs = parse_failures(hgx_dir + 'failed.txt')
        missing_alns = missing_alns.difference(set(failed_hgs))

    if os.path.basename(algorithm) == 'diamond':
        dmnd_dbs = set([int(os.path.basename(x[:-5])) \
                    for x in collect_files(hgx_dir, 'dmnd')])
        missing_dmnds = set(missing_alns).difference(dmnd_dbs)
        cmds1 = [(algorithm, 'makedb', '--db', f'{hgx_dir}{hg}.dmnd',
                  '--in', f'{hg_dir}{hg}.faa', '--threads', '2') \
                  for hg in missing_dmnds]
        if sensitivity:
            cmds2 = [((algorithm, 'blastp', '--query', f'{hg_dir}{hg}.faa',
                      '--db', f'{hgx_dir}{hg}.dmnd', '-o', f'{hgx_dir}{hg}.out.tmp',
                      '--outfmt', '6', 'qseqid', 'sseqid', 'ppos', 'pident',
                      '--threads', '2', '--id', str(minid), '--ultra-sensitive',
                      '--no-self-hits', '--max-target-seqs', str(len(hg2gene[hg])), '&&'), 
                      ('mv', f'{hgx_dir}{hg}.out.tmp', f'{hgx_dir}{hg}.out')) \
                     for hg in missing_alns]
        else:
            cmds2 = [((algorithm, 'blastp', '--query', f'{hg_dir}{hg}.faa',
                      '--db', f'{hgx_dir}{hg}.dmnd', '-o', f'{hgx_dir}{hg}.out.tmp',
                      '--outfmt', '6', 'qseqid', 'sseqid', 'ppos', 'pident',
                      '--threads', '2', '--id', str(minid), '--max-target-seqs', 
                      str(len(hg2gene[hg])), '--no-self-hits', '&&'), 
                      ('mv', f'{hgx_dir}{hg}.out.tmp', f'{hgx_dir}{hg}.out')) \
                     for hg in missing_alns]
    elif os.path.basename(algorithm) == 'blastp':
        cmds1 = []
        cmds2 = [((algorithm, '-query', f'{hg_dir}{hg}.faa', '-subject',
                   f'{hg_dir}{hg}.faa', '-out', f'{hgx_dir}{hg}.out.tmp',
                   '-outfmt', '6 "qseqid sseqid ppos pident"', '-num_threads',
                   str(cpus*2),
                   '-max_target_seqs', str(len(hg2gene[hg])), '&&'),
                   ('mv', f'{hgx_dir}{hg}.out.tmp', f'{hgx_dir}{hg}.out')) \
                  for hg in missing_alns]
    elif os.path.basename(algorithm) == 'mmseqs':
        mmseqs_dbs = set([int(os.path.basename(x[:-7])) \
                     for x in collect_files(hgx_dir, 'mmseqs')])
        missing_mmseqs = set(missing_alns).difference(mmseqs_dbs)
        cmds1 = [(algorithm , 'createdb', f'{hg_dir}{hg}.faa', f'{hgx_dir}{hg}.mmseqs',
                  '--createdb-mode', '1', '--shuffle', '0') for hg in missing_mmseqs]
        cmds2 = [((algorithm, 'search', f'{hgx_dir}{hg}.mmseqs', f'{hgx_dir}{hg}.mmseqs',
                  f'{hgx_dir}{hg}.raw', f'{hgx_dir}tmp{hg}', '--threads', str(cpus),
                  '--num-iterations', '3', '-s', '7.5', '-e', '0.001', 
                  '--max-seqs', str(len(hg2gene[hg])), '--max-rejected', '10', 
                  '--min-ungapped-score', '30', '&&'),
                  (algorithm, 'filterdb', f'{hgx_dir}{hg}.raw', f'{hgx_dir}{hg}.raw.filter', 
                   '--comparison-operator', 'ge', '--comparison-value', str(minid/100), 
                   '--filter-column', '2', '--threads', str(cpus), '&&'),
                  (algorithm, 'convertalis', f'{hgx_dir}{hg}.mmseqs', f'{hgx_dir}{hg}.mmseqs',
                   f'{hgx_dir}{hg}.raw.filter', f'{hgx_dir}{hg}.out.tmp', '--format-output', 
                   'query,target,pident', '--threads', str(cpus), '&&'),
                  ('mv', f'{hgx_dir}{hg}.out.tmp', f'{hgx_dir}{hg}.out', '&&'),
                  ('rm', '-rf', f'{hgx_dir}tmp{hg}', f'{hgx_dir}{hg}.raw', 
                   f'{hgx_dir}{hg}.raw.filter')) for hg in missing_alns]

    return cmds1, cmds2


def run_blast(hgs, db, hg_dir, hgx_dir, algorithm = 'diamond', 
              printexit = False, sensitivity = '', hg2gene = None,
              skipalgn = False, fallback = False, minid = 30, cpus = 1):
#    hgs = sorted(set(chain(*list(hgx2loc.keys()))))
    db_cmds, algn_cmds = prep_blast_cmds(db, hgs, hg_dir, hgx_dir, 
                                         minid = minid, rerun = not skipalgn,
                                         algorithm = algorithm, 
                                         sensitivity = sensitivity,
                                         hg2gene = hg2gene, cpus = cpus)
    if algn_cmds:
        if db_cmds:
            with open(hgx_dir + '../makedb.sh', 'w') as out:
                out.write('\n'.join([' '.join(x) for x in db_cmds]))
        with open(hgx_dir + '../srch.sh', 'w') as out:
            if len(algn_cmds[0]) == 2:
                for algn, tmp_mv in algn_cmds:
                    out.write(' '.join(algn) + '\n' + ' '.join(tmp_mv) + '\n')
            else:
                for algn, fltr, convert, tmp_mv, rm in algn_cmds:
                    out.write(' '.join(algn[:-1]) + ' && \n ' \
                            + ' '.join(fltr[:-1]) + ' && \n ' \
                            + ' '.join(convert[:-1]) + ' && \n ' \
                            + ' '.join(tmp_mv[:-1]) + ' && \n ' \
                            + ' '.join(rm) + ' \n\n')
        if printexit:
            print(f'\n{algorithm} db commands outputted to ' \
                 + f'{hgx_dir}../makedb.sh; run this first', flush = True)
            print(f'\n{algorithm} commands outputted to {hgx_dir}../srch.sh',
                  flush = True)
            sys.exit(0)

        elif not skipalgn:
            if db_cmds:
                print(f'\tBuilding {len(db_cmds)} aligner DBs', flush = True)
                multisub(db_cmds, verbose = 2, processes = cpus)
            print(f'\tAligning {len(algn_cmds)} HGs', flush = True)
            if algorithm == 'diamond':
                multisub(algn_cmds, verbose = 2, processes = cpus,
                         injectable = True)
            else: # leave the multithreading optimization to the program
    #            multisub(algn_cmds, verbose = 2, processes = round((cpus - 1)/2),
     #                    injectable = True)
                # launch 1 subprocess to minimize python overhead
                subprocess.run([os.environ['SHELL'], f'{hgx_dir}../srch.sh'],
                                stdout = subprocess.DEVNULL)
    
    missing_alns = find_missing_algns(hgs, hgx_dir)
    if not fallback:
        print(f'\t\t{len(missing_alns)} failed alignments will be skipped', 
              flush = True)
    elif missing_alns:
        print(f'\t\t{len(missing_alns)} failed alignments will fallback to diamond', 
              flush = True)
        db_cmds, algn_cmds = prep_blast_cmds(db, hgs, hg_dir, hgx_dir, 
                                         minid = minid, rerun = True,
                                         algorithm = 'diamond', 
                                         sensitivity = None,
                                         hg2gene = hg2gene, cpus = cpus)
        if db_cmds:
            print(f'\tBuilding {len(db_cmds)} diamond DBs', flush = True)
            multisub(db_cmds, verbose = 2, processes = cpus)
        print(f'\tAligning {len(algn_cmds)} HGs', flush = True)
        multisub(algn_cmds, verbose = 2, processes = cpus,
                 injectable = True)
        missing_alns = find_missing_algns(hgs, hgx_dir)

       
    with open(hgx_dir + 'failed.txt', 'w') as out:
        out.write('\n'.join([str(x) for x in missing_alns]))


def gcc_main(
    hgx2loc, wrk_dir, ome2i, hg_dir, hgx_dir,
    algorithm, db, gene2hg, plusminus, hg2gene, phylo,
    old_path = 'hgx2omes2gcc.pickle',
    gcfs = None, gcf_hgxs = None,
    gcf_omes = None, gcf2clan = {}, 
    cpus = 1, printexit = False,
    algn_sens = '', skipalgn = False, minid = 30,
    fallback = False
    ):


    if not checkdir(hgx_dir, unzip = True, rm = True):
        os.mkdir(hgx_dir)

    if os.path.isfile(wrk_dir + old_path):
        print('\tLoading previous coevolution results', flush = True)
        with open(wrk_dir + 'hgx2omes2gcc.full.pickle', 'rb') as pickin:
            d2gcc = pickle.load(pickin)
        with open(wrk_dir + 'hgx2omes2id.full.pickle', 'rb') as pickin:
            d2id_ = pickle.load(pickin)
        if os.path.isfile(wrk_dir + 'hgx2omes2pos.full.pickle'):
            with open(wrk_dir + 'hgx2omes2pos.full.pickle', 'rb') as pickin:
                d2pos = pickle.load(pickin)
        else:
            d2pos = {}
    else:
        d2gcc, d2id_, d2pos = {}, {}, {}
    
    hgs = list(chain(*list(gcf_hgxs.values())))
    hg_dir = hg_fa_mngr(wrk_dir, None, hgs, 
                        db, hg2gene, cpus = cpus,
                        low_mem = True)
    run_blast(hgs, db, hg_dir, hgx_dir, #minid = minid, let minid be at locus sim step
              algorithm = algorithm, printexit = printexit,
              sensitivity = algn_sens, hg2gene = hg2gene,
              skipalgn = skipalgn, fallback = fallback,
              cpus = cpus)

    d2gcc, d2id_, d2pos, gcfs = gcc_mngr(
        list(hgs), list(ome2i.keys()), hgx_dir, hgx2loc,
        db, gene2hg, plusminus, hg2gene, gcfs,
        gcf_omes, gcf_hgxs, ome2i, minid, phylo,
        d2gcc, d2id_, d2pos, cpus = cpus #d2pos, cpus = cpus
        )
#    hgx_dirTar = mp.Process(target=tardir, args=(hgx_dir, True))
 #   hgx_dirTar.start() # when to join ...
    with open(f'{wrk_dir}gcfs.pickle', 'wb') as out:
        pickle.dump([gcfs, gcf_omes, gcf_hgxs, gcf2clan], out)
    with open(wrk_dir + 'hgx2omes2gcc.full.pickle', 'wb') as pickout:
        pickle.dump(d2gcc, pickout)
    with open(wrk_dir + 'hgx2omes2id.full.pickle', 'wb') as pickout:
        pickle.dump(d2id_, pickout)
    if d2pos:
        with open(wrk_dir + 'hgx2omes2pos.full.pickle', 'wb') as pickout:
            pickle.dump(d2pos, pickout)

    return d2gcc, d2id_, d2pos, gcfs
