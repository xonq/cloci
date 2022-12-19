import os
import sys
import gzip
import hashlib
from cogent3 import PhyloNode, load_tree
from collections import defaultdict
from itertools import combinations
from mycotools.lib.biotools import gff2list
from mycotools.lib.kontools import format_path, eprint


def init_log(
    log_file, log_dict
    ):
    with open(log_file, 'w') as out:
        out.write(f'mtdb\t{log_dict["mtdb"]}\n' \
                + f'focal_genes\t{log_dict["focal_genes"]}' \
                + f'microsynt_constraint\t{log_dict["microsynt_constraint"]}' \
                + f'microsynt_tree\t{log_dict["microsynt_tree"]}\n'
                + f'hg_file\t{log_dict["hg_file"]}\n' \
                + f'plusminus\t{log_dict["plusminus"]}\n' \
                + f'pair_percentile\t{log_dict["pair_percentile"]}\n' \
                + f'hgx_percentile\t{log_dict["hgx_percentile"]}\n' \
                + f'border_percentile\t{log_dict["border_percentile"]}\n' \
                + f'id_percent\t{log_dict["id_percent"]}\n' \
                + f'pos_percent\t{log_dict["pos_percent"]}\n' \
                + f'patch_threshold\t{log_dict["patch_threshold"]}\n' \
                + f'coevo_threshold\t{log_dict["coevo_threshold"]}\n' \
                + f'null_samples\t{log_dict["null_samples"]}\n' \
                + f'n50\t{log_dict["n50"]}')

def read_log(
    log_file, log_dict
    ):
    log_res = {}
    with open(log_file, 'r') as raw:
        for line in raw:
            key = line[:line.find('\t')]
            res = line[line.find('\t') + 1:].rstrip()
            if res != str(log_dict[key]):
                log_res[key] = False
            else:
                log_res[key] = True
    try:
        if not log_res['mtdb']:
            log_res['n50'] = False
        if not log_res['focal_genes']:
            log_res['microsynt_tree'] = False
        if not log_res['microsynt_constraint']:
            log_res['microsynt_tree'] = False
        if not log_res['n50']:
            log_res['plusminus'] = False
        if not log_res['plusminus']:
            log_res['null_samples'] = False
        if not log_res['null_samples']:
            log_res['pair_percentile'] = False
        if not log_res['pair_percentile']:
            log_res['border_percentile'] = False
        if not log_res['border_percentile']:
            log_res['hgx_percentile'] = False
    except KeyError:
        print('\nERROR: corrupted log.txt.' + \
            '\nIf not rectified, future runs will completely overwrite the current\n')
        sys.exit(149)

    return log_res



def rm_old_data(
    log_res, out_dir, wrk_dir
    ):
    if not log_res['null_samples']:
        nulls = collect_files(wrk_dir + 'null/', 'null.txt')
        for null in nulls:
            os.remove(null)
    if not log_res['pair_percentile']:
        seed_file = out_dir + 'seed_scores.tsv.gz'
        seed_arr = wrk_dir + '.arr.npy'
        if os.path.isfile(seed_file):
            os.remove(seed_file)
        if os.path.isfile(seed_arr):
            os.remove(seed_arr)
        clus_pickle = wrk_dir + 'hgx2loc.pickle'
        hgx_pickle = wrk_dir + 'hgx_scores.pickle'
        ome_pickle = wrk_dir + 'hgx_omes.pickle'
        if os.path.isfile(clus_pickle):
            os.remove(clus_pickle)
        if os.path.isfile(hgx_pickle):
            os.remove(hgx_pickle)
        if os.path.isfile(ome_pickle):
            os.remove(ome_pickle)
    if not log_res['border_percentile']:
        row_file = wrk_dir + 'mtx/mcl.prep.rows'
        prep_file = wrk_dir + 'mtx/mcl.prep.gz'
        if os.path.isfile(row_file):
            os.remove(row_file)
        if os.path.isfile(prep_file):
            os.remove(prep_file)
    if not log_res['hgx_percentile']:
        kern_file = out_dir + 'hgx_clans.tsv.gz'
        clus_file = out_dir + 'hgxs.tsv.gz'
        patch_pickle = wrk_dir + 'patchiness.scores.pickle'
        ome_dir = wrk_dir + 'ome/'
        hgx_dir = wrk_dir + 'hgx/'
        hmm_dir = wrk_dir + 'hmm/'
        if os.path.isfile(kern_file):
            os.remove(kern_file)
        if os.path.isfile(clus_file):
            os.remove(clus_file)
        if os.path.isdir(ome_dir):
            shutil.rmtree(ome_dir)
        if os.path.isdir(hgx_dir):
            shutil.rmtree(hgx_dir)
        if os.path.isdir(hmm_dir):
            shutil.rmtree(hmm_dir)
        if os.path.isfile(wrk_dir + 'hgx.tar.gz'):
            os.remove(wrk_dir + 'hgx.tar.gz')
        if os.path.isfile(patch_pickle):
            os.remove(patch_pickle)
    if not log_res['hg_file']:
        shutil.rmtree(wrk_dir)
        os.mkdir(wrk_dir)
        for key in log_res:
            log_res[key] = False


def log_check(log_dict, log_path, out_dir, wrk_dir, flag = True):

    if not os.path.isfile(log_path):
        log_res = {x: False for x in log_dict}
        rm_old_data(log_res, out_dir, wrk_dir)
        init_log(log_path, log_dict)
    log_res = read_log(log_path, log_dict)
    if any(not log_res[x] for x in log_res):
        if not flag:
            print('\nInitializing new run', flush = True)
            rm_old_data(log_res, out_dir, wrk_dir)
            init_log(log_path, log_dict)
        else:
            failed = set([x for x, v in log_res.items() if not v])
            for skip in ['patch_threshold', 'coevo_threshold',
                         'pos_percent', 'id_percent']:
                if skip in failed:
                    failed.remove(skip)
            if failed:
                eprint('\nERROR: -n not called and incompatible parameters: \
                        \n\t' + ','.join([x for x,v in log_res.items() if not v]),
                   flush = True)
                sys.exit(15)

    return log_res


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
                data = [x.rstrip() for x in line.split('\t')]
#                if float(data[3]) > seed_thresh:
                out_hgs.append(line.split('\t'))
    
    return out_hgs


def compile_tree(i2ome, tree_path, root = []):
    phylo = load_tree(tree_path)
    if root:
      #  if len(root) > 1:
#            from ete3 import Tree
 #           t = Tree(tree_path)
  #          ancestor = t.get_common_ancestor(root[0], root[1])
   #         t.set_outgroup(ancestor)
    #        with open(tree_path, 'w') as out:
     #           out.write(t.write())
      #      phylo = load_tree(tree_path)
      #      phylo = phylo.rooted_at(
       #         phylo.get_edge_names(root[0], root[1])[0]
        #            )
      
       # else:
        phylo = phylo.rooted_with_tip(root[0])

        phylo.write(tree_path, with_distances = True)

    phylo.reassign_names({v: str(i) for i, v in enumerate(i2ome)})
    return phylo


def init_run(db, out_dir, near_single_copy_genes, constraint_path,
             tree_path, hg_file, plusminus, seed_perc, clus_perc,
             hgx_perc, id_perc, pos_perc, patch_thresh, coevo_thresh,
             samples, n50thresh, flag):
    wrk_dir = out_dir + 'working/'
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)

    nul_dir = wrk_dir + 'null/'
    if not os.path.isdir(nul_dir):
        os.mkdir(nul_dir)

    log_path = out_dir + 'log.txt'
    if not tree_path:
        tree_path = out_dir + 'microsynt.newick'
    log_dict = {'mtdb': hashlib.sha256(str(sorted(db.keys())).encode('utf-8')).hexdigest(),
        'focal_genes': hashlib.sha256(','.join(sorted(
                                     near_single_copy_genes
                                     )).encode('utf-8')).hexdigest(),
        'microsynt_constraint': format_path(constraint_path),
        'microsynt_tree': tree_path,
        'hg_file': hg_file, 'plusminus': plusminus,
        'pair_percentile': seed_perc, 'hgx_percentile': clus_perc,
        'border_percentile': hgx_perc, 'id_percent': id_perc,
        'pos_percent': pos_perc, 'patch_threshold': patch_thresh,
        'coevo_threshold': coevo_thresh,
        'null_samples': samples, 'n50': n50thresh}

    log_res = log_check(log_dict, log_path, out_dir, wrk_dir, flag)
    return wrk_dir, nul_dir
