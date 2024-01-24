import re
import os
import sys
import gzip
import shutil
import hashlib
import multiprocessing as mp
from tqdm import tqdm
from ast import literal_eval
from cogent3 import PhyloNode, load_tree
from collections import defaultdict
from itertools import combinations
from datetime import datetime
from mycotools.acc2fa import dbmain as acc2fa
from mycotools.lib.biotools import gff2list, dict2fa, fa2dict
from mycotools.lib.kontools import format_path, eprint, collect_files, collect_dirs

def init_log(
    log_file, log_dict
    ):
    with open(log_file, 'w') as out:
        out.write(f'mtdb\t{log_dict["mtdb"]}\n' \
                + f'focal_genes\t{log_dict["focal_genes"]}\n' \
                + f'microsynt_constraint\t{log_dict["microsynt_constraint"]}\n' \
                + f'microsynt_tree\t{log_dict["microsynt_tree"]}\n'
                + f'hg_file\t{log_dict["hg_file"]}\n' \
                + f'dist_type\t{log_dict["dist_type"]}\n' \
                + f'uniq_sp\t{log_dict["uniq_sp"]}\n' \
                + f'plusminus\t{log_dict["plusminus"]}\n' \
                + f'hgp_percentile\t{log_dict["hgp_percentile"]}\n' \
                + f'hgx_percentile\t{log_dict["hgx_percentile"]}\n' \
                + f'aligner\t{log_dict["aligner"]}\n' \
#                + f'hlg_percentile\t{log_dict["hlg_percentile"]}\n' \
                + f'gene_id\t{log_dict["gene_id"]}\n' \
                + f'hlg_id\t{log_dict["hlg_id"]}\n' \
                + f'hlg_sim\t{log_dict["hlg_sim"]}\n' \
                + f'orig_hlg_id\t{log_dict["hlg_id"]}\n' \
                + f'topology_merge\t{log_dict["topology_merge"]}\n' \
                + f'min_topology_sim\t{log_dict["min_topology_sim"]}\n' \
                + f'domain_inflation\t{log_dict["domain_inflation"]}\n' \
                + f'hlg_inflation\t{log_dict["hlg_inflation"]}\n' \
                + f'tuning\t{log_dict["tuning"]}\n' \
                + f'id_percent\t{log_dict["id_percent"]}\n' \
                + f'pos_percent\t{log_dict["pos_percent"]}\n' \
                + f'iipr_percent\t{log_dict["iipr_percent"]}\n' \
                + f'patch_threshold\t{log_dict["patch_threshold"]}\n' \
                + f'gcl_threshold\t{log_dict["gcl_threshold"]}\n' \
                + f'dist_threshold\t{log_dict["dist_threshold"]}\n' \
                + f'partition\t{log_dict["partition"]}\n' \
                + f'null_samples\t{log_dict["null_samples"]}\n' \
                + f'n50\t{log_dict["n50"]}\n' \
                + f'hg_dir\t{log_dict["hg_dir"]}\n' \
                + f'hgx_dir\t{log_dict["hgx_dir"]}\n' \
                + f'ipr_path\t{log_dict["ipr_path"]}\n' \
                + f'pfam_path\t{log_dict["pfam_path"]}')

def read_log(
    log_file, log_dict
    ):
    log_res = {}
    with open(log_file, 'r') as raw:
        for line in raw:
            key = line[:line.find('\t')]
            res = line[line.find('\t') + 1:].rstrip()
            if key not in log_dict: # old feature
                log_res[key] = False
                continue
            elif key == 'orig_hlg_id':
                if round(float(res) * 100) \
                   <= round(float(log_dict['hlg_id']) * 100):
                   log_res[key] = True
                else:
                   log_dict['orig_hlg_id'] = log_dict['hlg_id']
                   log_res[key] = False
            elif key == 'hgx_percentile':
                if round(float(res) * 100) > \
                    round(float(log_dict['hgx_percentile']) * 100):
                    log_res['hgx_percentile_I'] = False
                    log_res['hgx_percentile_II'] = False
                elif round(float(res) * 100) != round(float(log_dict['hgx_percentile']) * 100):
                    log_res['hgx_percentile_I'] = True
                    log_res['hgx_percentile_II'] = False
                else:
                    log_res['hgx_percentile_I'] = True
                    log_res['hgx_percentile_II'] = True
            elif key == 'hlg_id':
                continue
            elif key == 'hlg_sim':
                if res.lower() != str(log_dict[key]).lower():
                    log_res[key] = False
                else:
                    log_res[key] = True
            elif res != str(log_dict[key]):
                log_res[key] = False
                if key == 'tuning':
                    # if you're tuning and weren't before
                    if log_dict['tuning']:
                        # then the old inflation is irrelevant
                        log_res['domain_inflation'] = False
                        log_res['hlg_inflation'] = False
                    # elif you're not tuning, were before, and inflations are the same
                    elif log_res['domain_inflation'] and log_res['hlg_inflation']:
                        log_res['tuning'] = True
                        log_dict['tuning'] = literal_eval(res)
                    # elif the domain inflation is the same and the hlg is not
 #                   else:
#                        log_res['hlg_inflation'] = False
            else:
                log_res[key] = True
                if key == 'tuning' and log_dict['tuning']:
                    if inflation_rnd1:
                        log_res['domain_inflation'] = True
                    else:
                        inflation_rnd2 = None
                        log_res['domain_inflation'] = False
                    if inflation_rnd2:
                        log_res['hlg_inflation'] = True
                    else:
                        log_res['hlg_inflation'] = False
            # if there's an inflation, tuning is complete
            if key == 'domain_inflation':
                if log_dict['tuning']:
                    inflation_rnd1 = literal_eval(res)
                else: # else the inflation is the inflation in the log_dict
                    inflation_rnd1 = log_dict['domain_inflation']
            elif key == 'hlg_inflation':
                # if you're tuning 
                if log_dict['tuning'] and log_res['hlg_inflation']:
                    inflation_rnd2 = literal_eval(res)
                else:
                    inflation_rnd2 = log_dict['hlg_inflation']
             

    init_discrep = []
    try:
        if not log_res['mtdb']:
            init_discrep.append('-d')
            log_res['n50'] = False
    except KeyError:
        init_discrep.append('-d')
    try:
        if not log_res['partition']:
            init_discrep.append('-nr/-np')
            log_res['null_samples'] = False
    except KeyError:
        init_discrep.append('-nr/-np')
    try:
        if not log_res['focal_genes']:
            init_discrep.append('-f')
            log_res['microsynt_tree'] = False
    except KeyError:
        init_discrep.append('-f')
    try:
        if not log_res['microsynt_constraint']:
            init_discrep.append('-c')
            log_res['microsynt_tree'] = False
    except KeyError:
        init_discrep.append('-c')
    try:
        if not log_res['n50']:
            init_discrep.append('--n50')
            log_res['plusminus'] = False
    except KeyError:
        init_discrep.append('--n50')
    try:
        if not log_res['dist_type']:
            init_discrep.append('-mmd')
            log_res['uniq_sp'] = False
    except KeyError:
        init_discrep.append('-mmd')
    try:
        if not log_res['uniq_sp']:
            init_discrep.append('-u')
            log_res['null_samples'] = False
    except KeyError:
        init_discrep.append('-u')
    try:
        if not log_res['plusminus']:
            init_discrep.append('-w')
            log_res['null_samples'] = False
    except KeyError:
        init_discrep.append('-w')
    try:
        if not log_res['null_samples']:
            init_discrep.append('-ns')
            log_res['hgp_percentile'] = False
    except KeyError:
        init_discrep.append('-ns')
    try:
        if not log_res['hgp_percentile']:
            init_discrep.append('-hp')
            log_res['hgx_percentile_I'] = False
    except KeyError:
        init_discrep.append('-hp')
    try:
        if not log_res['hgx_percentile_I']:
            init_discrep.append('-xp')
            log_res['hgx_percentile_II'] = False
    except KeyError:
        init_discrep.append('-xp')
    try:
        if not log_res['hgx_percentile_II']:
            init_discrep.append('-xp')
            log_res['topology_merge'] = False
    except KeyError:
        init_discrep.append('-xp')
    try:
        if not log_res['topology_merge']:
            log_res['min_topology_sim'] = False
            init_discrep.append('-tm')
    except KeyError:
        init_discrep.append('-tm')
    try:
        if not log_res['min_topology_sim']:
            init_discrep.append('-ts')
            log_res['hlg_sim'] = False
    except KeyError:
        init_discrep.append('-ts')
    try:
        if not log_res['aligner']:
            init_discrep.append('-a')
            log_res['gene_id'] = False
    except KeyError:
        log_res['aligner'] = False
        init_discrep.append('-a')
    try:
        if not log_res['gene_id']:
            init_discrep.append('-mg')
            log_res['hlg_sim'] = False
    except KeyError:
        init_discrep.append('-mg')
    try:
        if not log_res['hlg_sim']:
            init_discrep.append('-s')
            log_res['orig_hlg_id'] = False
    except KeyError:
        init_discrep.append('-s')
    try:
        if not log_res['orig_hlg_id']:
            init_discrep.append('-ml')
            log_res['domain_inflation'] = False
    except KeyError:
        init_discrep.append('-ml')
    try:
        if not log_res['domain_inflation']:
            init_discrep.append('-I2/-T')
            log_res['hlg_inflation'] = False
    except KeyError:
        init_discrep.append('-I1/-T')
        print('\nERROR: corrupted log.txt.' + \
            '\nIf not rectified, future runs will completely overwrite the current\n')
        sys.exit(149)
    try:
        if not log_res['hlg_inflation']:
            init_discrep.append('-I2/-T')
#            log_res['hlg_percentile'] = False
    except KeyError:
        init_discrep.append('-I2/-T')
        print('\nERROR: corrupted log.txt\nFuture runs may overwrite the current\n')
        sys.exit(149)
#    try:
 #       if not log_res['hlg_percentile']:
  #          init_discrep.append('-fp')
   # except KeyError:
    #    init_discrep.append('-fp')

    # stop gap for legacy data
    if 'hg_dir' not in log_res:
        log_res['hg_dir'] = True
    if 'hgx_dir' not in log_res:
        log_res['hgx_dir'] = True

    return log_res, inflation_rnd1, inflation_rnd2, init_discrep


def rm_old_data(
    log_res, out_dir, wrk_dir
    ):

    todel, tosave, save_dir = [], [], None
    if not log_res['dist_type'] or not log_res['uniq_sp']:
        if os.path.isfile(out_dir + 'working/omes2dist.pickle'):
            os.remove(out_dir + 'working/omes2dist.pickle')
    if not log_res['null_samples']:
        if os.path.isdir(wrk_dir + 'null/'):
            shutil.rmtree(wrk_dir + 'null/')
    if not log_res['hgp_percentile']:
        seed_file = out_dir + 'hgps.tsv.gz'
        seed_arr = wrk_dir + 'microsynt.npy'
        if os.path.isfile(seed_file):
            tosave.append(seed_file)
        if os.path.isfile(seed_arr):
            os.remove(seed_arr)
        clus_pickle = wrk_dir + 'hgx2loc.pickle'
        if os.path.isfile(clus_pickle):
            os.remove(clus_pickle)
    if not log_res['hgx_percentile_I']:
        skip_files = collect_files(wrk_dir + 'null/', '*')
        skip_files = [x for x in skip_files \
                      if os.path.basename(x).startswith('skip.')]
        for skip_file in skip_files:
            os.remove(skip_file)
        done_file = [x for x in collect_files(f'{wrk_dir}null/', 'txt') \
                     if os.path.basename(x).startswith('done')]
        if done_file:
            os.remove(done_file[0])
    if not log_res['hgx_percentile_II']:
        hlg_dir = wrk_dir + 'hlg/'
        if os.path.isdir(hlg_dir):
            todel.append(hlg_dir)
        groups = ['hgcs', 'gene2hgx']
        for group in groups:
             if os.path.isfile(f'{wrk_dir}{group}.pickle'):
                 os.remove(f'{wrk_dir}{group}.pickle')
    if not log_res['min_topology_sim']:
        groups = ['group.III']
        clan_data = ['hg.json.gz', 'hgx.json.gz', 'json.gz']
        for file_ in clan_data:
            if os.path.isfile(f'{wrk_dir}clan2loci.{file_}'):
                os.remove(f'{wrk_dir}clan2loci.{file_}')
        for group in groups:
             if os.path.isfile(f'{wrk_dir}{group}.pickle'):
                 os.remove(f'{wrk_dir}{group}.pickle')       
        todel.append(f'{wrk_dir}hlg/')
    if not log_res['aligner']:
        hgx_dir = f'{wrk_dir}hgx/'
        eprint('\tWARNING: alignment software changed, old alignments retained:', flush = True)
        eprint(f'\t\t{hgx_dir}', flush = True)
 #       shutil.rmtree(hgx_dir)
    if not log_res['orig_hlg_id']:
        hlg_dir = f'{wrk_dir}hlg/'
        adj_rows = collect_files(hlg_dir, 'tmp.w')
        adj_rows.extend(collect_files(hlg_dir, 'tmp.r'))
        adj_rows.extend(collect_files(hlg_dir, 'tmp'))
        row_file = wrk_dir + 'hlg/loci.adj'
        mci_file = wrk_dir + 'hlg/loci.mci'
        mcl_rows = wrk_dir + 'hlg/mcl_rows.tsv'
        for f in [row_file, mci_file, mcl_rows]:
            if os.path.isfile(f):
                todel.append(f)

    if not log_res['hlg_inflation']:
        hlgs_file = f'{wrk_dir}hlgs.pickle'
        if os.path.isfile(hlgs_file):
            os.remove(hlgs_file)
        mcl_res = f'{wrk_dir}hlg/loci.clus'
        tosave.extend([mcl_res, mcl_res + '.tmp'])
        hgx_files = ['gcl.pickle',
                     'mmi.pickle',
                     'mmp.pickle']
        for file_ in hgx_files:
            if os.path.isfile(f'{wrk_dir}{file_}'):
                os.remove(f'{wrk_dir}{file_}')

#    if not log_res['hlg_percentile']:
        hlg_output = out_dir + 'hlgs.tsv.gz'
        gcf_output = out_dir + 'gcfs.tsv.gz'
        kern_file = out_dir + 'hgxs.tsv.gz'
        ome_dir = out_dir + 'ome/'
        ann_dir = wrk_dir + 'ann/'
        tosave.extend([hlg_output, gcf_output])
        if os.path.isdir(ome_dir):
            count = 0
            save_dir = f'{out_dir}run{count}/'
            while os.path.isdir(save_dir):
                count += 1
                save_dir = f'{out_dir}run{count}/'
            os.mkdir(save_dir)
            os.mkdir(save_dir + 'working/')
            os.mkdir(save_dir + 'working/hlg/')
            os.mkdir(save_dir + 'working/hlg/lg')
            shutil.move(ome_dir, save_dir)
            if os.path.isfile(kern_file):
                shutil.copy(kern_file, save_dir + os.path.basename(kern_file))
            if os.path.isfile(hlg_output):
                shutil.move(hlg_output, save_dir + os.path.basename(hlg_output))
            if os.path.isfile(out_dir + 'hlgs.html'):
                shutil.move(out_dir + 'hlgs.html', save_dir + 'hlgs.html')
            if os.path.isfile(out_dir + 'hlg_metrics.html'):
                shutil.move(out_dir + 'hlg_metrics.html', 
                            save_dir + 'hlg_metrics.html')
            if os.path.isfile(gcf_output):
                shutil.move(gcf_output, save_dir + os.path.basename(clus_file))
            if os.path.isfile(out_dir + 'gcfs.html'):
                shutil.move(out_dir + 'gcfs.html', save_dir + 'gcfs.html')
            if os.path.isfile(out_dir + 'gcf_metrics.html'):
                shutil.move(out_dir + 'gcf_metrics.html', 
                            save_dir + 'gcf_metrics.html')

            shutil.move(out_dir + 'log.txt', save_dir + 'log.txt')
            for f in collect_files(f'{wrk_dir}hlg/', '*'):
                shutil.copy(f, f'{save_dir}working/hlg/{os.path.basename(f)}')
            for f in collect_files(f'{wrk_dir}hlg/lg/', '*'):
                shutil.copy(f, f'{save_dir}working/hlg/lg/{os.path.basename(f)}')

#        else:
 #           if os.path.isfile(clus_file):
  #              os.remove(clus_file)
    
#        if os.path.isdir(hmm_dir):
 #           shutil.rmtree(hmm_dir)
  #      if os.path.isdir(f'{out_dir}net/'):
   #         shutil.move(f'{out_dir}net/', f'{save_dir}net/')

    if not log_res['hg_dir']:
        if os.path.isdir(f'{wrk_dir}hg/'):
            shutil.rmtree(f'{wrk_dir}hg/')
    if not log_res['hgx_dir']:
        if os.path.isdir(f'{wrk_dir}hgx/'):
            shutil.rmtree(f'{wrk_dir}hgx/')

    if save_dir:
        for f in tosave:
            try:
                shutil.move(f, save_dir + os.path.basename(f))
            except FileNotFoundError:
                continue
    else:
        for f in tosave:
            try:
                os.remove(f)
            except FileNotFoundError:
                continue
    for f in todel:
        if os.path.isdir(f):
            shutil.rmtree(f)
        elif os.path.isfile(f):
            os.remove(f)

    if not log_res['hg_file']:
        shutil.rmtree(wrk_dir)
        os.mkdir(wrk_dir)
        for key in log_res:
            log_res[key] = False


def log_check(log_dict, log_path, out_dir, wrk_dir, flag = True):

    gcf_thresh = False
    if not os.path.isfile(log_path):
        log_res = {x: False for x in log_dict}
        try:
            rm_old_data(log_res, out_dir, wrk_dir)
        except KeyError:
            eprint('\nWARNING: missing log, rerunning without deleting', flush = True)
        init_log(log_path, log_dict)
    log_res, inflation_1, inflation_2, init_discrep = read_log(log_path, log_dict)
    if any(not log_res[x] for x in log_res):
        failed = set([x for x, v in log_res.items() if not v])
        for skip in ['patch_threshold', 'gcl_threshold', 'iipr_percent',
                     'id_percent', 'pos_percent', 'dist_threshold', 'ipr_path',
                     'pfam_path']:
            if skip in failed:
                failed.remove(skip)
        if failed:
            try:
                if not flag:
                    print(f'\n{init_discrep[0]} changed; removing data to checkpoint', flush = True)
                    rm_old_data(log_res, out_dir, wrk_dir)
                    init_log(log_path, log_dict)
                else:
                    eprint('\nERROR: -n not called and incompatible parameters: \
                            \n\t' + init_discrep[0], flush = True)
                    sys.exit(15)
            # irrelevant missing discrepancy
            except IndexError:
                init_log(log_path, log_dict)
        else:
            if os.path.isdir(out_dir + 'ome/'):
                tsvs = set(collect_files(out_dir + 'ome/', 'tsv', recursive = True))
                # NEED an intelligent resume
                ome_dirs = [f'{out_dir}ome/{x}/' for x in os.listdir(out_dir + 'ome/') \
                            if os.path.isdir(f'{out_dir}ome/{x}/')]
                if all(f'{x}hlg.tsv' in tsvs for x in ome_dirs) or os.path.isfile(out_dir + 'hlgs.tsv.gz'):
                    print('\nBypassing to filtration', flush = True)
                    gcf_thresh = True
                else:
                    init_log(log_path, log_dict)
            else:
                init_log(log_path, log_dict)
    elif os.path.isdir(out_dir + 'ome/'):
        tsvs = set(collect_files(out_dir + 'ome/', 'tsv', recursive = True))
        ome_dirs = [f'{out_dir}ome/{x}/' for x in os.listdir(out_dir + 'ome/') \
                    if os.path.isdir(f'{out_dir}ome/{x}/')]
        if all(f'{x}hlg.tsv' in tsvs for x in ome_dirs):
            print('\nBypassing to filtration', flush = True)
            gcf_thresh = True


    return log_res, inflation_1, inflation_2, gcf_thresh


def compileCDS(gff_list, ome):
    """ 
    Inputs the gff_list and organism ome code. Compiles CDS entries for loci parsing and outputs:
    cds_dict = {contig: {protein: [[CDSstart, CDSstop]] } }
    Sorts the cds_dict for each protein from smallest to largest
    """     
                
    cds_dict, fail = defaultdict(dict), False
    for entry in gff_list:
        if entry['type'] == 'CDS':
            prot_prep_i0 = entry['attributes'].index(';Alias=') # grab the
            # mycotools accession index
            try: # grab the end of the Alias by grabbing the index for a ';',
            # that is after the accession tag's occurrence.
            # then add that to the start of the accession tag and add one to
            # remove that semicolon. If there is a value error, the accession 
            # is at the end of the attributes (gff col 8 (?))
                prot_prep_i1 = entry['attributes'][prot_prep_i0+1:].index(';') + prot_prep_i0+1
            except ValueError: 
                prot_prep_i1 = len(entry['attributes'])
            prot = entry['attributes'][prot_prep_i0:prot_prep_i1].replace(';Alias=','')
            # obtain the protein accession
            if prot not in cds_dict[entry['seqid']] and prot: # add the protein
            # to the cds_dict entry
                cds_dict[entry['seqid']][prot] = [] 
            elif not prot: # if there isn't a valid accession it may mean the
            # mycotools curation did not work or the user did not curate
            # correctly 
                print(entry['attributes'], prot_prep_i0, prot_prep_i1)
                if not fail:
                    print('\tWARNING: ' + ome + ' has proteins in gff with no Alias', flush = True)
                fail = True
                continue    
            cds_dict[entry['seqid']][prot].extend([int(entry['start']),
                int(entry['end'])]) # add the coordinates
                        
    for contig in cds_dict:
        for prot in cds_dict[contig]:
            cds_dict[contig][prot].sort() # sort the coordinates of the proteins
            # lowest to highest
        cds_dict[contig] = list(sorted(cds_dict[contig].keys(), key = lambda k:
            cds_dict[contig][k][0])) # sort the proteins in the contigs lowest
            # to highest coordinates

    return dict(cds_dict)

def compileCDS2(gff_list, ome):
    """
    Inputs the gff_list and organism ome code. Compiles CDS entries for loci parsing and outputs:
    cds_dict = {contig: {protein: [[CDSstart, CDSstop]] } }
    Sorts the cds_dict for each protein from smallest to largest
    """

    cds_dict, cds_dict2, fail = defaultdict(dict), {}, False
    for entry in gff_list:
        if entry['type'] == 'CDS':
            prot_prep_i0 = entry['attributes'].index(';Alias=') # grab the
            # mycotools accession index
            try: # grab the end of the Alias by grabbing the index for a ';',
            # that is after the accession tag's occurrence.
            # then add that to the start of the accession tag and add one to
            # remove that semicolon. If there is a value error, the accession
            # is at the end of the attributes (gff col 8 (?))
                prot_prep_i1 = entry['attributes'][prot_prep_i0+1:].index(';') + prot_prep_i0+1
            except ValueError:
                prot_prep_i1 = len(entry['attributes'])
            prot = entry['attributes'][prot_prep_i0:prot_prep_i1].replace(';Alias=','')
            # obtain the protein accession
            if prot not in cds_dict[entry['seqid']] and prot: # add the protein
            # to the cds_dict entry
                cds_dict[entry['seqid']][prot] = []
                cds_dict2[prot] = []
            elif not prot: # if there isn't a valid accession it may mean the
            # mycotools curation did not work or the user did not curate
            # correctly
                print(entry['attributes'], prot_prep_i0, prot_prep_i1)
                if not fail:
                    print('\tWARNING: ' + ome + ' has proteins in gff with no Alias', flush = True)
                fail = True
                continue
            cds_dict[entry['seqid']][prot].extend([int(entry['start']),
                int(entry['end'])]) # add the coordinates
            cds_dict2[prot].append(entry)

    for contig in cds_dict:
        for prot in cds_dict[contig]:
            cds_dict[contig][prot].sort() # sort the coordinates of the proteins
            # lowest to highest
        cds_dict[contig] = list(sorted(cds_dict[contig].keys(), key = lambda k:
            cds_dict[contig][k][0])) # sort the proteins in the contigs lowest
            # to highest coordinates

    return dict(cds_dict), cds_dict2


def parseLoci(
    gff_path, ome_num, gene2hg, plusminus = 6
    ):
    """obtain a set of tuples of HG pairs {(OG0, OG1)...}"""

    gff_list = gff2list(gff_path) # open here to improve pickling
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))
    pairs = []
    for scaf in cds_dict: # for each contig
        for i, seq in enumerate(cds_dict[scaf]): # for each gene
            if i < plusminus: # does the locus start at the contig border?
                locus = cds_dict[scaf][:i+plusminus+1]
            else: # or I extend it to the +/- or opposite border
                locus = cds_dict[scaf][i-plusminus:i+plusminus+1]
            loc_og = []
            for gene in locus:
                try: # failure should be relatively rare, so faster than `if`
                    loc_og.append(gene2hg[gene])
                except KeyError: # it's either a singleton or not in an OG
                    pass
            pairs.extend([tuple(sorted(x)) for x in combinations(set(loc_og), 2)]) # don't deal with
            # same OGs - tandem duplications will just naturally cooccur more often

    out_pairs = set(pairs) # unique pairs of OGs
    return ome_num, out_pairs


def load_seedScores(file_):#, seed_thresh):

    out_hgs = []
    with gzip.open(file_, 'rt') as raw:
        for line in raw:
            if not line.startswith('#'):
                d = line.rstrip().split()
                out_hgs.append((int(d[0]), int(d[1])))
    
    return out_hgs


def compile_tree(i2ome, tree_path, root = []):
    with open(tree_path, 'r') as raw:
        nwk = raw.read()

    phylo = load_tree(tree_path)
    omes = set(i2ome)
    tips = set(phylo.get_tip_names())
    missing = tips.difference(omes)
    if missing:
        eprint(f'\tWARNING: removing {len(missing)} ' \
              + 'missing tips from microsynteny tree',
               flush = True)
        todel = []
        for n in phylo.tips():
            if n.name in missing:
                todel.append(n)
        for n in todel:
            n.parent.remove(n)
            phylo.prune()
        
    if len(root) == 1:
        phylo = phylo.rooted_with_tip(root[0])
        phylo.write(tree_path, with_distances = True)
    elif len(root) > 1:
        nodes = {k: (v, len(v.get_tip_names())) \
                 for k, v in phylo.get_nodes_dict().items() \
                 if set(root).issubset(set(v.get_tip_names()))}
        mrca_tip_len = min([v[1] for v in list(nodes.values())])
        mrca_edge = [k for k, v in nodes.items() if v[1] == mrca_tip_len]
        phylo = phylo.rooted_at(mrca_edge)
        phylo.write(tree_path, with_distances = True)

    phylo.reassign_names({v: str(i) for i, v in enumerate(i2ome)})
    return phylo


def write_hg_fa(hg_file, big_dict, genes):
    fa_str = dict2fa({k: big_dict[k] for k in genes})
    with open(f'{hg_file}.tmp', 'w') as out:
        out.write(fa_str)
    os.rename(f'{hg_file}.tmp', hg_file)


def output_hg_fas(db, genes, hg_file):
    fa_str = dict2fa(acc2fa(
                db, genes, error = False, spacer = '\t\t',
                coord_check = False
                ))
    with open(hg_file + '.tmp', 'w') as out:
        out.write(fa_str)
    os.rename(hg_file + '.tmp', hg_file)


def cp_files(f0, f1):
    shutil.copy(f0, f1)

def symlink_files(f0, f1):
    os.symlink(f0, f1)


def load_ome2fa(db, cpus = 1):
#    big_fa = {}
#    with mp.Pool(processes = cpus) as pool:
 #       fa_dicts = pool.map(fa2dict, 
  #                          tqdm((row['faa'] for row in db.values()), 
   #                         total = len(db)))
    #    pool.close()
     #   pool.join()
#    for res in fa_dicts:
 #       big_fa = {**big_fa, **res}
    fas = []
    for row in tqdm(db.values(), total = len(db)):
        fas.append(fa2dict(row['faa']))
    big_fa = {}
    for fa in fas:
        for nam, seq in fa.items():
            big_fa[nam] = seq

    return big_fa

def load_ome2fa_mp(db, cpus = 1):
    big_fa = {}
    with mp.Pool(processes = cpus) as pool:
        fa_dicts = pool.map(fa2dict, 
                            tqdm((row['faa'] for row in db.values()), 
                            total = len(db)))
        pool.close()
        pool.join()
    for res in fa_dicts:
        big_fa = {**big_fa, **res}

    return big_fa

def big_acc2fa(db, hg_dir, hgs, hg2gene, big_fa, cpus = 1):
    for hg in tqdm(hgs, total = len(hgs)):
        write_hg_fa(f'{hg_dir}{hg}.faa', big_fa, hg2gene[hg])



def big_acc2fa_mp(db, hg_dir, hgs, hg2gene, big_fa, cpus = 1):
    big_dict = mp.Manager().dict(big_fa)
    with mp.Pool(processes = cpus) as pool:
        pool.starmap(write_hg_fa, tqdm(((f'{hg_dir}{hg}.faa',
                                         big_dict, hg2gene[hg]) \
                                         for hg in hgs),
                                       total = len(hgs)))
        pool.close()
        pool.join()

    
def hg_fa_mngr(wrk_dir, hg_dir, hgs,
               db, hg2gene, cpus = 1, low_mem = False):
    new_hg_dir = wrk_dir + 'hg/'
    if not os.path.isdir(new_hg_dir):
        os.mkdir(new_hg_dir)
    # extract only used HGs
    hg_files = set([int(os.path.basename(x[:-4])) \
                    for x in collect_files(new_hg_dir, 'faa')])
    missing_hgs = sorted(set(hgs).difference(hg_files))

    if not missing_hgs:
        return new_hg_dir

    if not hg_dir:
        if low_mem: # this process is way too slow
            with mp.Pool(processes = cpus) as pool:
                pool.starmap(output_hg_fas,
                             tqdm(((db, hg2gene[hg], f'{new_hg_dir}{hg}.faa') \
                             for hg in missing_hgs), total = len(missing_hgs)))
                pool.close()
                pool.join()
        else:
            big_fa = load_ome2fa(db, cpus = cpus)
            big_acc2fa(db, new_hg_dir, missing_hgs, hg2gene, big_fa, cpus)
    else: # predetermined hg input
        hg_fa_cmds = []
        orthofinder = False
        if not os.path.isfile(f'{hg_dir}{hgs[0]}.faa'):
             digits = len(str(hgs[0]))
             zeros = 7 - digits
             if os.path.isfile(f'{hg_dir}OG{"0" * zeros}{hgs[0]}.fa'):
                 orthofinder = True
                 ext = '.fa'
             elif os.path.isfile(f'{hg_dir}{hgs[0]}.fa'):
                 ext = '.fa'
             elif os.path.isfile(f'{hg_dir}{hgs[0]}.fasta'):
                 ext = '.faa'
        else:
             raise FileNotFoundError('HG fastas missing')

        if orthofinder:
             for hg in missing_hgs:
                 digits = len(str(hg))
                 zeros = 7 - digits
                 hg_file = (f'{hg_dir}OG{"0" * zeros}{hg}.fa')
                 hg_fa_cmds.append((hg_file, f'{new_hg_dir}{hg}.faa'))
        else:
            hg_fa_cmds = [(f'{hg_dir}{hg}{ext}', f'{new_hg_dir}{hg}.faa') \
                          for hg in missing_hgs]
        copy_hgs = False
        try:
            os.symlink(hg_fa_cmds[0][0], hg_fa_cmds[0][1])
            del hg_fa_cmds[0]
        except:
            copy_hgs = True

        if copy_hgs:
            with mp.Pool(processes = cpus) as pool:
                pool.starmap(cp_files, hg_fa_cmds)
                pool.close()
                pool.join()
        else:
            with mp.Pool(processes = cpus) as pool:
                pool.starmap(symlink_files, hg_fa_cmds)
                pool.close()
                pool.join()
    return new_hg_dir


def sha_tune_file(tune_file):
    with open(tune_file, 'r') as raw:
        tune_data = [x.rstrip() for x in raw if not x.startswith('#')]
    tune = {}
    for cluster in tune_data:
        try:
            name, rawgenes, rawomes, rawfalse = cluster.split('\t')
        except ValueError:
            name, rawgenes, rawomes = cluster.split('\t')
            rawfalse = None
        if name in tune:
            eprint(f'\nERROR: duplicate entry for {name} in tune file', 
                   flush = True)
        if ',' in rawgenes:
            genes = [x.rstrip().lstrip() for x in rawgenes.split(',')]
        else:
            genes = [x.rstrip().lstrip() for x in rawgenes.split()]
        if ',' in rawomes:
            omes = [x.rstrip().lstrip() for x in rawomes.split(',')]
        else:
            omes = [x.rstrip().lstrip() for x in rawomes.split()]
        if rawfalse:
            if ',' in rawfalse:
                false_omes = [x.rstrip().lstrip() for x in rawfalse.split(',')]
            else:
                false_omes = [x.rstrip().lstrip() for x in rawfalse.split()]
        else:
            false_omes = []
        tune[name] = [tuple(sorted(set(genes))), tuple(sorted(set(omes))),
                      tuple(sorted(set(false_omes)))]
        
    tune_sha = hashlib.sha256(str(sorted(tune)).encode('utf-8')).hexdigest()
    return tune_sha


def init_run(db, out_dir, near_single_copy_genes, constraint_path,
             tree_path, hg_file, plusminus, seed_perc,# clus_perc,
             hgx_perc, aligner, id_perc, pos_perc, iipr_perc,
             patch_thresh, gcl_thresh, dist_thresh,
             samples, n50thresh, flag, min_gene_id, min_hlg_id, inflation_rnd1, 
             inflation_rnd2, simfun,
             tune_file, dist_type, uniq_sp, partition, min_topology_sim,
             topology_merge, hg_dir, hgx_dir, ipr_path, pfam_path):

    wrk_dir = out_dir + 'working/'
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)

    nul_dir = wrk_dir + 'null/'
    if not os.path.isdir(nul_dir):
        os.mkdir(nul_dir)

    log_path = out_dir + 'log.txt'
    if not tree_path:
        tree_path = out_dir + 'microsynt.newick'

    if tune_file:
        tune_sha = sha_tune_file(tune_file)
    else:
        tune_sha = None

    log_dict = {'mtdb': hashlib.sha256(str(sorted(db.keys())).encode('utf-8')).hexdigest(),
        'focal_genes': hashlib.sha256(','.join(sorted(
                                     near_single_copy_genes
                                     )).encode('utf-8')).hexdigest(),
        'microsynt_constraint': format_path(constraint_path),
        'microsynt_tree': tree_path,
        'hg_file': hg_file, 'plusminus': plusminus, 'dist_type': dist_type,
        'uniq_sp': bool(uniq_sp),
        'hgp_percentile': seed_perc, 'hgx_percentile': hgx_perc, 
        'aligner': aligner, 'gene_id': min_gene_id,
        'orig_hlg_id': None, 'hlg_id': min_hlg_id, 'hlg_sim': str(simfun).lower(),
#        'hlg_percentile': clus_perc, 
        'id_percent': id_perc, 'pos_percent': pos_perc, 'iipr_percent': iipr_perc,
        'patch_threshold': patch_thresh, 'dist_threshold': dist_thresh,
        'gcl_threshold': gcl_thresh, 'domain_inflation': inflation_rnd1,
        'hlg_inflation': inflation_rnd2,
        'tuning': tune_sha, 'partition': partition, 
        'topology_merge': topology_merge, 'min_topology_sim': min_topology_sim,
        'null_samples': samples, 'n50': n50thresh, 'hg_dir': hg_dir,
        'hgx_dir': hgx_dir, 'ipr_path': ipr_path, 'pfam_path': pfam_path}

    log_res, inflation_1, inflation_2, gcf_bypass = log_check(log_dict, log_path, 
                                                              out_dir, wrk_dir, flag)

    # symlink files on an individual basis so the new dir is writeable without
    # propagating new changes back to the reference dir
    hg_dir, hgx_dir = wrk_dir + 'hg/', wrk_dir + 'hgx/'

    existing_hg = set(os.path.basename(x) for x in collect_files(hg_dir, 'faa'))
    existing_hgx = set(os.path.basename(x) for x in collect_files(hgx_dir, 'out'))
    if log_dict['hg_dir']:
        if not os.path.isdir(hg_dir):
            os.mkdir(wrk_dir + 'hg/')
        prev_hg = set(os.path.basename(x) for x in collect_files(log_dict['hg_dir'], 'faa'))
        missing_hg = prev_hg.difference(existing_hg)
        for i in list(missing_hg):
            os.symlink(log_dict['hg_dir'] + i, wrk_dir + 'hg/' + i)
    if log_dict['hgx_dir']:
        if not os.path.isdir(hgx_dir):
            os.mkdir(wrk_dir + 'hgx/')
        prev_hgx = set(os.path.basename(x) for x in collect_files(log_dict['hgx_dir'], 'out'))
        missing_hgx = prev_hgx.difference(existing_hgx)
        for i in list(missing_hgx):
            os.symlink(log_dict['hgx_dir'] + i, wrk_dir + 'hgx/' + i)

    return wrk_dir, nul_dir, inflation_1, inflation_2, gcf_bypass
