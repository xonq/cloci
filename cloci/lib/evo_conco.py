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
from cloci.lib.input_parsing import compileCDS2, hg_fa_mngr
from cloci.lib.treecalcs import calc_tmd

# NEED todel_hg to remove fully incompatible GCL HGs from HLG_HGXs
# NEED an option to import OrthoFinder pairwise alignments
	# not recommended, need higher max target sequences than OrthoFinder 
def parse_algn_wpos(hgx_dir, hg, hgx_genes):
    """Parse alignment data that contains percent positives"""
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
    """Parse alignment data that does not contain percent positives (MMseqs)"""
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


def get_top_hits_wpos(gene, hlg, hgx_genes, ome_set, hits, gene2algn, res, ome):
    """Identify the top hits of a query by sorting based on alignment
    identity"""
    ome_hits = set()
    ome_len = len(ome_set) - 1
    for sbj_gene, sbj_id, sbj_pos in gene2algn[gene]:
        sbj_ome = sbj_gene[:sbj_gene.find('_')]
    # if the subject is also in the checked queries and it is not
    # a recent paralog of the query species
        if sbj_gene in hgx_genes and sbj_ome != ome:
        # grab the hit similarity
            res[hlg][gene][1][sbj_gene] = (sbj_id, sbj_pos,)
            hits.add(sbj_gene)
            ome_hits.add(sbj_ome)
            if ome_len == len(ome_hits):
                return hits, res 
    # add one to the failed count
        elif sbj_ome not in ome_set: 
            res[hlg][gene][0] += 1 
    
    if len(ome_hits) < ome_len:
        # missing gene and no paralog
        res[hlg][gene][2] = True 

    return hits, res 
        
def get_top_hits_wopos(gene, hlg, hgx_genes, ome_set, hits, gene2algn, res, ome):
    """Identify the top hits of a query by sorting based on alignment
    identity (MMseqs)"""
    ome_len = len(ome_set) - 1
    ome_hits = set()
    for sbj_gene, sbj_id in gene2algn[gene]:
        sbj_ome = sbj_gene[:sbj_gene.find('_')]
        # if the subject is also in the checked queries and it is not
        # a recent paralog of the query species
        if sbj_gene in hgx_genes and sbj_ome != ome:
            # grab the hit similarity
            res[hlg][gene][1][sbj_gene] = sbj_id
            hits.add(sbj_gene)
            ome_hits.add(sbj_ome)
            if len(ome_hits) == ome_len:
                return hits, res
        # add one to the failed count
        elif sbj_ome not in ome_set: 
            res[hlg][gene][0] += 1 
        
    if len(ome_hits) < ome_len:
        # missing gene and no paralog
        # biased against convergent assimilation
        res[hlg][gene][2] = True 

    return hits, res


def sim_calc_wpos(res):
    """Calclate the minimum similarity for identity and positives across the genes
    within a shared HG"""
    omeScores = defaultdict(dict)
    ids_y_pos = defaultdict(dict)
    for hg in res:
        for ome, data in res[hg].items():
            t_gcl, id_y_pos_data = data
            omeScores[ome][hg] = t_gcl
            ids = [x[0] for x in id_y_pos_data.values()]
            pos = [x[1] for x in id_y_pos_data.values()]
            if ids:
                ids_y_pos[hg][ome] = (min(ids), min(pos),)
            else:
                ids_y_pos[hg][ome] = (0, 0,)
    return omeScores, ids_y_pos


def sim_calc_wopos(res):
    """Calclate the minimum similarity for identity across the genes
    within a shared HG"""
    omeScores = defaultdict(dict)
    ids_dict = defaultdict(dict)
    for hg in res:
        for ome, data in res[hg].items():
            t_gcl, id_data = data
            omeScores[ome][hg] = t_gcl
            ids = list(id_data.values())
            if ids:
                ids_dict[hg][ome] = min(ids)
            else:
                ids_dict[hg][ome] = 0
    return omeScores, ids_dict

# deprecated
def mms_calcs_wpos(ids_y_pos):
    # average minimum identity / minimum positive; for each gene homolog group
    hlg_min_id, hlg_min_pos, total = 0, 0, 0
    for ome_dict in ids_y_pos.values():
        hlg_min_id += sum([x[0] for x in ome_dict.values()])
        hlg_min_pos += sum([x[1] for x in ome_dict.values()])
        total += len(ome_dict)
    try:
        hlg_min_id /= total
        hlg_min_pos /= total
    except ZeroDivisionError:
        hlg_min_id, hlg_min_pos = 0, 0
    return hlg_min_id, hlg_min_pos

# deprecated
def mms_calcs_wopos(ids_dict):
    # average minimum identity; for each gene homolog group
    hlg_min_id, total = 0, 0
    for ome_dict in ids_dict.values():
        hlg_min_id += sum(ome_dict.values())
        total += len(ome_dict)
    try:
        hlg_min_id /= total
    except ZeroDivisionError:
        hlg_min_id = 0
    return hlg_min_id, None


def hg_sim_calc_wpos(res, min_id, hlg2tot_con):
    """Calculate and return the average GCL, MMI, and MMP for all HGs shared
    within and HLG for each ome"""
    output_data = {}
    for hlg, gene_info in res.items():
        output_data[hlg] = {}
        tot_con = hlg2tot_con[hlg]
        for gene, data in gene_info.items():
            ome = gene[:gene.find('_')]
            ome_gcl = data[0]
            if not data[2]:
                try:
                    mi_p = [v[0] for v in data[1].values()]
                    ome_mi = min(mi_p)
                    mi_i = mi_p.index(ome_mi)
                    mp_p = [v[1] for v in data[1].values()]
                    ome_mp = mp_p[mi_i]
                except ValueError:
                    ome_mi, ome_mp = 0, 0
            else:
                ome_mi, ome_mp = 0, 0
    #            ome_mi, ome_mp = min_id, min_id
            if ome not in output_data[hlg]:
                output_data[hlg][ome] = (ome_gcl, ome_mi * tot_con, 
                                         ome_mp * tot_con)
            # if there's a paralog, take the best gcl
            else:
                if ome_gcl > output_data[hlg][ome][0]:
                    output_data[hlg][ome] = (ome_gcl, ome_mi * tot_con, 
                                             ome_mp * tot_con)
    return output_data


# deprecated
def hg_sim_calc_wopos(res, min_id, hlg2tot_con):
    """Calculate and return the average GCl and MMI for all HGs shared
    within an HLG for each ome, without considering MMP"""
    output_data = {}
    for hlg, gene_info in res.items():
        output_data[hlg] = {}
        tot_con = hlg2tot_con[hlg]
        for gene, data in gene_info.items():
            ome = gene[:gene.find('_')]
            ome_gcl = data[0]
            if not data[2]:
                try:
                    ome_mi = min(data[1].values())
                except ValueError:
                    ome_mi = 0
            else:
                ome_mi = 0
#                ome_mi = min_id
            if ome not in output_data[hlg]:
                output_data[hlg][ome] = (ome_gcl, ome_mi * tot_con)
            # if there's a paralog, take the best gcl
            else:
                if ome_gcl > output_data[hlg][ome][0]:
                    output_data[hlg][ome] = (ome_gcl, ome_mi * tot_con)
    return output_data


def hg_parse_and_calc(hg, hg_dict, hgx_dir, ome2i, min_id = 30, pos_data = True):
    """Manage the parsing and calculation of commitment to the locus and
    minimum similarity (identity and positives) for a given HG and the HLGs it
    is a part of"""

    # if there is percent positive data then use the functions that correspond
    # to that
    if pos_data:
       parse_func = parse_algn_wpos
       top_hits_func = get_top_hits_wpos
       hg_sim_func = hg_sim_calc_wpos
    # deprecated: lack of percent positives data [MMseqs]
    else:
       parse_func = parse_algn_wopos
       top_hits_func = get_top_hits_wopos
       hg_sim_func = hg_sim_calc_wopos

    gene2hlg = {}
    for hlg, genes in hg_dict.items():
        # ignore singletons
        if len(genes) > 1:
            for gene in genes:
                gene2hlg[gene] = hlg

    # grab the genes in an hg
    hg_genes = set(gene2hlg.keys())

    # parse the self-alignment results
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
        return None, None, None, None
    except FileNotFoundError:
        return None, None, None, None

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
    res = defaultdict(lambda: defaultdict(lambda: [0, {}, False]))
    todel_genes = defaultdict(lambda: defaultdict(list))

    # create a data structure that maintains all of the considered genes of an
    # HG in a particular HLG
    hlg2considered = {}
    for hlg, genes in hg_dict.items():
        failed_genes = True
        gene_set = set(genes)
        omes_set = set(x[:x.find('_')] for x in genes)

        # if there is a failed gene (too low identity to call a homolog, then
        # continually rerun until those have been exhaustively removed from
        # consideration
        while failed_genes:
            failed_genes = []
            # disregard the self from the length
            hgx_ome_len = len(omes_set) - 1
            if hgx_ome_len == 0:
                ome = genes[0][:genes[0].find('_')]
                todel_genes[hlg][ome].extend(genes)
                failed_genes.extend(genes)
                break
            # for the weighted average, we can start at this step
            # only consider unique omes because we are only taking the max
            considered_genes = len(set(x[:x.find('_')] for x in genes))
            hlg2considered[hlg] = considered_genes
            for gene in genes:
                ome = gene[:gene.find('_')] # identify the ome
#                if gene not in res[hlg]:
 #                   res[hlg][gene] = [0, {}, False]
        
                # to identify if the subject is in the family of omes
                hits = {gene}
                # while all cluster homologs aren't accounted for and there remain hits
                if gene not in gene2algn:
                    continue
                hits, res = top_hits_func(gene, hlg, gene_set, omes_set, hits, 
                                                     gene2algn, res, ome)
        
                # if there are missing hits (which means we are disregarding
                # missing alignments)
                if hgx_ome_len - len(res[hlg][gene][1]) > 0:
                    # then get the percent of observed hits that are valid
                    hit_len = res[hlg][gene][0] + len(res[hlg][gene][1])
                    try:
#* considered_genes \

                        res[hlg][gene][0] = len(res[hlg][gene][1]) / hit_len
                    except ZeroDivisionError:
#                        eprint(hlg, gene, gene_set, omes_set, res[hlg][gene], gene2algn, flush = True)
                        # no valid hits, the gene should not be called in the cluster
                        # typically happens from extension module
                        todel_genes[hlg][ome].append(gene)
                        failed_genes.append(gene)
                else:
                    # the GCL for this HG and ome is the total of the genes in this
                    # HG (disregarding self) divided by the total + the failed genes
#* considered_genes \

                    res[hlg][gene][0] = hgx_ome_len \
                                      / (hgx_ome_len + res[hlg][gene][0])
            # have to rerun the calculation if there are failed genes
            for failed_gene in failed_genes:
                gene_set.remove(failed_gene)
            genes = list(gene_set)
            omes_set = set(x[:x.find('_')] for x in genes)
   
    return hg, hg_sim_func(res, min_id, hlg2considered), \
           {hlg: dict(ome_dict) for hlg, ome_dict in todel_genes.items()}, \
           hlg2considered


def gcl_mngr(
    hgs, omes, hgx_dir, hgx2loc,
    db, gene2hg, clusplusminus, hg2gene, hlgs,
    hlg_omes, hlg_hgxs, ome2i, min_id, phylo,
    d2gcl, d2id_, d2pos, cpus = 1
    ):
    """Manage the calculation of GCL and similarity measurements (MMI and MMP)
    for each HLG"""
    i2ome = {v: k for k, v in ome2i.items()}

    # identify all the HGs that need to be parsed and organize them by their
    # HLG
    clus_hgs = {}
    for hlg, loci in hlgs.items():
        hlg_hgx = hlg_hgxs[hlg]
        omesc = hlg_omes[hlg]
        if hlg_hgx not in d2id_:
            clus_hgs[hlg] = [hlg_hgx, defaultdict(list)]
            hlg_hgx_set = set(hlg_hgx)
            for loc in loci:
                for gene in loc:
                    try:
                        hg = gene2hg[gene]
                    except KeyError:
                        continue
                    if hg in hlg_hgx_set:
                        clus_hgs[hlg][1][hg].append(gene)
        elif omesc not in d2id_[hlg_hgx]:
            clus_hgs[hlg] = [hlg_hgx, defaultdict(list)]
            hlg_hgx_set = set(hlg_hgx)
            for loc in loci:
                for gene in loc:
                    try:
                        hg = gene2hg[gene]
                    except KeyError:
                        continue
                    if hg in hlg_hgx_set:
                        clus_hgs[hlg][1][hg].append(gene)

            # {fam: [hgx, {hg: set(seqs)}}

    clus_hgs = {
        hlg: [d[0], {hg: set(v) for hg, v in d[1].items() if len(set(v)) > 1}] \
        for hlg, d in clus_hgs.items()
        } # make sets from it
    hg2genes = defaultdict(dict)
    for hlg, hgs in clus_hgs.items():
        for hg, gene_set in hgs[1].items():
            hg2genes[hg][hlg] = tuple(sorted(gene_set))

    # parse the alignments and acquire the preliminary quantitations of
    # similarity and GCL measurements on an HG-by-HG basis
    print('\tParsing alignments', flush = True)
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        hg_results = pool.starmap(hg_parse_and_calc,
                              tqdm(((hg, hg_dict, hgx_dir, ome2i, min_id) \
                               for hg, hg_dict in hg2genes.items()),
                               total = len(hg2genes)))
        pool.close()
        pool.join()

    # prepare for quantifying the overall HLG GCL, MMI, MMP, and eventual CSB
    # by removing any loci and genomes that do not have sufficient
    # justification to warrant their retention following minimum similarity
    # thresholding
    print('\tQuantifying', flush = True)
    hlg_res = defaultdict(dict)
    deleted = 0
    hlg2tot_considered = defaultdict(int)
    hlgs = {hlg: list(locs) for hlg, locs in hlgs.items()}
    for hg, hg_res, todel_genes, hlg2considered in hg_results:
        todel = defaultdict(list)
        if hg_res:
            # acquire the total genes considered for each HLG
            for hlg, considered in hlg2considered.items():
                hlg2tot_considered[hlg] += considered
            for hlg, ome_dict in hg_res.items():
                hlg_res[hlg][hg] = []
                for ome, res in ome_dict.items():
                    hlg_res[hlg][hg].append(res)
            for hlg, ome_fails in todel_genes.items():
                for ome, failed_genes in ome_fails.items():
                    fails = set(failed_genes)
                    for i, loc in enumerate(hlgs[hlg]):
                        if loc[0][:loc[0].find('_')] == ome:
                            # if all/all but 1 of the genes failed, then delete it
                            if len([x for x in loc if x in fails]) >= len(loc) - 1:
                                todel[hlg].append(i)
            # remove any loci that now lack sufficient homology to justify them
            # being a part of an HLG, due to them have too low identity to be
            # considered
            for hlg, iz in todel.items():
                deleted += len(iz)
                for i in sorted(set(iz), reverse = True):
                    del hlgs[hlg][i]
    print(f'\t{deleted} loci removed due to lacking homology', flush = True)

    # delete any HLGs that have too few genomes to justify them
    todel = set()
    for hlg, locs in hlgs.items():
        omes = set(ome2i[loc[0][:loc[0].find('_')]] for loc in locs)
        if len(omes) > 2:        
            hlg_omes[hlg] = tuple(sorted(omes))
        else:
            todel.add(hlg)
    print(f'\t{len(todel)} HLGs removed for too few remaining omes', flush = True)
    for hlg in todel:
        del hlgs[hlg]
        del hlg_omes[hlg]
        del hlg_hgxs[hlg]
        
    # quantify the overall measurements for each HLG            
    hgx2omes2gcl = defaultdict(dict)
    hgx2omes2id = defaultdict(dict)
    hgx2omes2pos = defaultdict(dict)
    for hlg, hg_dict in tqdm(hlg_res.items(), total = len(hlg_res)):
        if hlg in todel:
            continue
        tot_considered = hlg2tot_considered[hlg]
        hgx, omes = hlg_hgxs[hlg], hlg_omes[hlg]
        gcls, ids, poss = [], [], []
        # calculate the average of each measurement for the overall HLG
        for hg, res in hg_dict.items():
            gcls.append(sum([x[0] for x in res]))
            ids.append(sum([x[1] for x in res])/len(res))
            try:
                poss.append(sum([x[2] for x in res])/len(res))
            except IndexError:
                poss.append(None)
            except TypeError:
                poss.append(None)
#        gcl = sum(gcls)/(len(gcls) * hlg2tot_considered[hlg])
        gcl = sum(gcls)/tot_considered
        id_ = sum(ids)/tot_considered
        try:
            pos = sum(poss)/tot_considered
        except ValueError:
            pos = None
        except TypeError:
            pos = None

        hgx2omes2gcl[hgx][omes] = gcl
        hgx2omes2id[hgx][omes] = id_
        if pos is not None:
            hgx2omes2pos[hgx][omes] = pos

    hgx2omes2pos = dict(hgx2omes2pos)
    return {**d2gcl, **hgx2omes2gcl}, \
        {**d2id_, **hgx2omes2id}, \
        {**d2pos, **hgx2omes2pos}, \
        {g: tuple(l) for g, l in hlgs.items()}, \
        {g: tuple(o) for g, o in hlg_omes.items()}


def find_missing_algns(hgs, hgx_dir):
    """Identify missing alignments that need to be run"""
    alnd_hgs = set([int(os.path.basename(x[:-4])) \
                    for x in collect_files(hgx_dir, 'out')])
    missing_alns = set(hgs).difference(alnd_hgs)
    return missing_alns


def parse_failures(fail_file):
    """Parse failed alignments file so that they are not continually rerun when
    it is unnecessary"""
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
    """Prepare commands for self-aligning HGs based on the inputted
    algorithm"""
    missing_alns = find_missing_algns(hgs, hgx_dir)

    # if failed runs are not to be rerun
    if not rerun:
        failed_hgs = parse_failures(hgx_dir + 'failed.txt')
        missing_alns = missing_alns.difference(set(failed_hgs))

    # prepare diamond commands
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
    # prepare BLASTp commands
    elif os.path.basename(algorithm) == 'blastp':
        cmds1 = []
        cmds2 = [((algorithm, '-query', f'{hg_dir}{hg}.faa', '-subject',
                   f'{hg_dir}{hg}.faa', '-out', f'{hgx_dir}{hg}.out.tmp',
                   '-outfmt', '6 "qseqid sseqid ppos pident"', '-num_threads',
                   str(cpus*2),
                   '-max_target_seqs', str(len(hg2gene[hg])), '&&'),
                   ('mv', f'{hgx_dir}{hg}.out.tmp', f'{hgx_dir}{hg}.out')) \
                  for hg in missing_alns]
    # deprecated - prepare MMseqs commands
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
    """Run HG self-alignment commands"""
#    hgs = sorted(set(chain(*list(hgx2loc.keys()))))
    db_cmds, algn_cmds = prep_blast_cmds(db, hgs, hg_dir, hgx_dir, 
                                         minid = minid, rerun = not skipalgn,
                                         algorithm = algorithm, 
                                         sensitivity = sensitivity,
                                         hg2gene = hg2gene, cpus = cpus)
    if algn_cmds:
        # make the databases if necessary (Diamond)
        if db_cmds:
            with open(hgx_dir + '../makedb.sh', 'w') as out:
                out.write('\n'.join([' '.join(x) for x in db_cmds]))
        # write and run the search commands
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
        # exit once the commands are outputted to allow for user-controlled
        # parallelization
        if printexit:
            print(f'\n{algorithm} db commands outputted to ' \
                 + f'{hgx_dir}../makedb.sh; run this first', flush = True)
            print(f'\n{algorithm} commands outputted to {hgx_dir}../srch.sh',
                  flush = True)
            sys.exit(0)

        # NEED to rewrite out what skipalgn does
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

    # identify missing alignments
    missing_alns = find_missing_algns(hgs, hgx_dir)
    # call these alignments failed
    if not fallback:
        print(f'\t\t{len(missing_alns)} failed alignments will be skipped', 
              flush = True)
    # if there are missing alignments and a fallback method, then run that
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

    # finally, output the failed alignments for future reruns to reference 
    with open(hgx_dir + 'failed.txt', 'w') as out:
        out.write('\n'.join([str(x) for x in missing_alns]))


def gcl_main(
    hgx2loc, wrk_dir, ome2i, hg_dir, hgx_dir,
    algorithm, db, gene2hg, plusminus, hg2gene, phylo,
    old_path = 'gcl.pickle',
    hlgs = None, hlg_hgxs = None,
    hlg_omes = None, hlg2clan = {}, 
    cpus = 1, printexit = False,
    algn_sens = '', skipalgn = False, minid = 30,
    fallback = False
    ):
    """The main function for controlling and executing the management of
    running GCL, MMI, MMP, and CSB quantitation"""

    # determine if the hgx directory exists, and unzip if necessary
    if not checkdir(hgx_dir, unzip = True, rm = True):
        os.mkdir(hgx_dir)

    # open existing data structures to rerun from what is already complete
    if os.path.isfile(wrk_dir + old_path):
        print('\tLoading previous GCL, MMI, and MMP results', flush = True)
        with open(wrk_dir + 'gcl.pickle', 'rb') as pickin:
            d2gcl = pickle.load(pickin)
        with open(wrk_dir + 'mmi.pickle', 'rb') as pickin:
            d2id_ = pickle.load(pickin)
        if os.path.isfile(wrk_dir + 'mmp.pickle'):
            with open(wrk_dir + 'mmp.pickle', 'rb') as pickin:
                d2pos = pickle.load(pickin)
        else:
            d2pos = {}
    else:
        d2gcl, d2id_, d2pos = {}, {}, {}

    # output fastas for any HGs that do not have them
    # NEED to figure out why this is necessary now when earlier steps should
    # have already completed this
    hgs = list(chain(*list(hlg_hgxs.values())))
    hg_dir = hg_fa_mngr(wrk_dir, None, hgs, 
                        db, hg2gene, cpus = cpus,
                        low_mem = True)

    # run HG self-alignment to set the stage for similarity calculations
    run_blast(hgs, db, hg_dir, hgx_dir, #minid = minid, let minid be at locus sim step
              algorithm = algorithm, printexit = printexit,
              sensitivity = algn_sens, hg2gene = hg2gene,
              skipalgn = skipalgn, fallback = fallback,
              cpus = cpus)

    # calculate the proxies of coordinated gene evolution
    d2gcl, d2id_, d2pos, hlgs, hlg_omes = gcl_mngr(
        list(hgs), list(ome2i.keys()), hgx_dir, hgx2loc,
        db, gene2hg, plusminus, hg2gene, hlgs,
        hlg_omes, hlg_hgxs, ome2i, minid, phylo,
        d2gcl, d2id_, d2pos, cpus = cpus #d2pos, cpus = cpus
        )
#    hgx_dirTar = mp.Process(target=tardir, args=(hgx_dir, True))
 #   hgx_dirTar.start() # when to join ...

    # output the data structures for future swift parsing and update the HLG
    # data based on any HLGs that were removed following determination of being
    # below minimum gene similarity thresholds
    with open(f'{wrk_dir}hlgs.pickle', 'wb') as out:
        pickle.dump([hlgs, hlg_omes, hlg_hgxs, hlg2clan], out)
    with open(wrk_dir + 'gcl.pickle', 'wb') as pickout:
        pickle.dump(d2gcl, pickout)
    with open(wrk_dir + 'mmi.pickle', 'wb') as pickout:
        pickle.dump(d2id_, pickout)
    if d2pos:
        with open(wrk_dir + 'mmp.pickle', 'wb') as pickout:
            pickle.dump(d2pos, pickout)

    return d2gcl, d2id_, d2pos, hlgs, hlg_omes
