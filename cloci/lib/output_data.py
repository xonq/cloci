import os
import re
import gzip
import shutil
import pickle
import subprocess
import multiprocessing as mp
import plotly.graph_objects as go
from math import log
from tqdm import tqdm
from itertools import chain
from collections import defaultdict
from plotly.subplots import make_subplots
from mycotools.acc2gff import grab_gff_acc
from mycotools.lib.kontools import eprint
from mycotools.lib.biotools import list2gff, gff2list, fa2dict, dict2fa
from cloci.lib import input_parsing

# NEED to add PCA output
# NEED locID output and clans
# NEED gbk output option
# NEED to add merge abutting clusters feature
# NEED annotation for gcfs resume and hlgs independent

def mk_subplots(labels, axes, axes_labels, alpha = 0.6):
    
    colors = [
        [240,163,255], [0,117,220], [153,63,0], [76,0,92], [25,25,25],
        [0,92,49], [43,206,72], [255,204,153], [128,128,128], [148,255,181],
        [143,124,0], [157,204,0], [194,0,136], [0,51,128], [255,164,5],
        [255,168,187], [66,102,0], [255,0,16], [94,241,242], [0,153,143],
        [224,255,102], [116,10,255], [153,0,0], [255,255,128], [255,255,0], [255,80,5]
        ]
#    random.shuffle(colors)
    pl_clr = [
        'rgba(' + ', '.join([str(y) for y in x]) + ', ' + str(alpha) + ')' \
        for x in colors
        ]
    fig = make_subplots(
        rows=len(axes) - 1, cols=len(axes) - 1,
        start_cell = "bottom-left"
        )

    clr_count = 0
    for i0, axis0 in enumerate(axes):
        for i1, axis1 in enumerate(axes[i0+1:]):
            fig.add_trace(
                go.Scatter(
                    mode = 'markers', x=axis0, y=axis1,
                    marker = dict(
                        color = pl_clr[clr_count],
                        size = 5,
                        opacity = alpha
                        ),
                    text = labels
                    ), row = i0+1, col = i0+i1+1,
                )
            fig.update_xaxes(
                title_text=axes_labels[i0], row = i0+1, col =
                i0+i1+1
                )
            fig.update_yaxes( 
                title_text=axes_labels[i0+i1+1], row = i0+1, col =
                i0+i1+1
                )
            clr_count += 1
    
    return fig 


def mk_3d(labels, axes, axes_labels, alpha = 0.6):

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        mode = 'markers',
        x = axes[0], y = axes[1], z = axes[2],
        marker = dict(
            color = 'rgba(163, 25, 169,' + str(alpha) + ')',
            size = 5, opacity = alpha
            ),
        text = labels, showlegend = False
        )
    )
    fig.update_layout(scene = dict(
                    xaxis_title=axes_labels[0],
                    yaxis_title=axes_labels[1],
                    zaxis_title=axes_labels[2]))#,
#                    width=700,
 #                   margin=dict(r=20, b=10, l=10, t=10))

    return fig


def extract_sig_clus(
    clus_scores_dict, hgx2loc, top_hgxs, ome2i
    ):
    """Create a hash to seed retrieving clusters based on
    their entry in hgx2loc."""
        
    sig_clus = defaultdict(dict)
    for top in top_hgxs:
        hgx, gcf, omeIs = top[0], top[1], top[2]
        for gene in hgx2loc[hgx]:
            ome = gene[:gene.find('_')]
            if ome2i[ome] in omeIs:
                if gene not in sig_clus[ome]:            
                    sig_clus[ome][gene] = [[set(hgx), gcf]]
                else:
                    sig_clus[ome][gene].append([set(hgx), gcf])
    
    return dict(sig_clus)


def hash_clusters(gff_path):
    """sliding window of size clusplusminus, identify sequences that may have
    hgs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus"""
    gff_list = gff2list(gff_path)
    cds_dict = input_parsing.compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))
    return cds_dict


def merge_clusters(clus_out, cds_dict):
    """because we have significant clusters with borders, we can merge each
    cluster and assume that at a certain number of hgs in a combo will never
    happen via random sampling; therefore, their significance in the
    microsynteny adjustment is sufficient for justification for merging"""
    
    comp_clus, change = defaultdict(list), False
    for scaf, clusters in clus_out.items():
        gene_list = cds_dict[scaf]
        while clusters: # while there are clusters
            any_intersect, toDel = False, [0] # we haven't found intersection
            clus0 = clusters[0] # grab the first tuple
            loc0 = clus0[0] # grab the genes in the cluster
            low_i0 = gene_list.index(loc0[0])
            hi_i0 = gene_list.index(loc0[-1])
            for i, clus1 in enumerate(clusters[1:]): # look at all other clusters
                loc1 = clus1[0] # grab their loci
                low_i1 = gene_list.index(loc1[0])
                hi_i1 = gene_list.index(loc1[-1])
                if low_i1 - hi_i0 == 1 and False:
                    any_intersect, change = True, True # we found overlap in
                    # clus_out and in this specific cluster comparison
 #                   hgx = tuple(sorted(list(set(clus0[0]).union(set(clus1[0])))))
                    # obtain the higher order og combination as the union of the
                    # two cluster sets, then sort it and store as a tuple
                    comp_clus[scaf].append([
#                        hgx, 
                        loc0 + loc1, 
                        clus0[1].union(clus1[1]) #, merge_x
                        ])
                    # append the og combo and the union of the two loci's genes
                    # protein order is random
                    toDel.append(i) # we can delete the intersected locus as well
                    # as the first
                elif low_i0 - hi_i1 == 1 and False:
                    any_intersect, change = True, True
                    comp_clus[scaf].append([loc1 + loc0, clus0[1].union(clus1[1])])
                    toDel.append(i)
                elif set(loc0).intersection(set(loc1)) and False:
                    any_intersect, change = True, True
                    comp_clus[scaf].append([set(loc0).union(set(loc1)), clus0[1].union(clus1[1])])
                    comp_clus[scaf][-1][0] = sorted(comp_clus[scaf][-1][0], 
                                                    key = cds_dict[scaf].index)   
                    toDel.append(i)
            if not any_intersect: # if there is not any intersection in the locus
                comp_clus[scaf].append(clus0) #, clus0[2]]) # simply append the original
                # result
            toDel.sort(reverse = True) # sort the todel from highest to lowest to
            # not mess up order when indexing
            for i in toDel:
                clusters.pop(i)
            
    return comp_clus, change


def clus2scaf(cds_dict, clusters):

    gene2scaf, index_map = {}, {}
    for scaf, genes in cds_dict.items():
        index_map[scaf] = {}
        for i, gene in enumerate(genes):
            index_map[scaf][gene] = i
            gene2scaf[gene] = scaf

    clus_out = defaultdict(list)
    for loc, gcf in clusters:
        gene0 = loc[0]
        scaf = gene2scaf[gene0]
        ord_loc = sorted(loc, key = lambda pair: index_map[scaf][pair])
#        hgs = []
 #       for gene in ord_loc:
  #          try:
   #             hgs.append(gene2hg[gene])
    #        except KeyError:
     #           hgs.append('')
      #  clus_info = [hgs, ord_loc, {gcf}]
        
        clus_out[scaf].append([ord_loc, gcf])

    sorted_clus_out = {k: sorted(v, key = lambda x: index_map[k][x[0][0]]) \
                       for k, v in clus_out.items()}

    return sorted_clus_out, index_map


def write_clusters(clusters, ome, out_file, gff_path, gene2hg, clusplusminus = 10):
    cds_dict = input_parsing.compileCDS(gff2list(gff_path), 
                                        os.path.basename(gff_path).replace('.gff3',''))

    clus_out, scaf2gene2i = clus2scaf(cds_dict, clusters)
    merge_out = defaultdict(list)
    for scaf, locs in clus_out.items():
        for loc, gcf in locs:
            merge_out[scaf].append([loc, {gcf}])
  #  change = True
#    while change: # so long as there is a single merge
 #       clus_out, change = merge_clusters(clus_out, cds_dict) # attempt to merge clusters
#    clus_out = sorted(clus_out, key = lambda x: len(x[0]), reverse = True)
    # sort the clusters by the size of the OG combination tuple highest to
    # lowest

    for scaf, clusters in merge_out.items():
        mapping = {v: i for i, v in enumerate(cds_dict[scaf])}
        clusters = sorted(clusters, key = lambda pair: mapping[pair[0][0]])
        for i, clus in enumerate(clusters):
            clus[0] = sorted(clus[0], key = cds_dict[scaf].index)
            clusHGs = []
            for gene in clus[0]:
                try:
                    clusHGs.append(gene2hg[gene])
                except KeyError: # gene not in an OG/is a singleton
                    clusHGs.append('')
            clus.append(clusHGs)
   #         scaf = clus[2]
    #        scaf = clus[3][:20]
     #       if clus[3] != scaf:
    #            scaf = clus[3][:20] + '..'
#            count = 0
 #           while scaf + '_' + str(count) in scafs:
  #              count += 1
            name = scaf + '_' + str(i)
#            scafs.add(name)
            clus.append(name)
        clus_out[scaf] = clusters
    

    out_genes = []
    with open(out_file, 'w') as out:
        out.write('#name\thgs\tgenes\tgcfs\n')
        for scaf, clusters in clus_out.items():
            for clus in clusters:
 #           clus[2] = sorted(clus[2])
                hgs = ','.join([str(x) for x in clus[2]]) # og1,og2,...ogn
                genes = ','.join(clus[0]) # prot1,prot2,prot3,...protn
                gcfs = ';'.join([str(x) for x in sorted(clus[1])])
                out.write(clus[3] + '\t' + hgs + '\t' + genes + '\t' + gcfs + '\n')
                out_genes.append(clus[0])

    return ome, out_genes


def read_ome_hlgs(ome, ome_hlg_file):
    out_genes = []
    try:
        with open(ome_hlg_file, 'r') as raw:
            for line in raw:
                if not line.startswith('#'):
#        out.write('#name\thgs\tgenes\tgcfs\n')
                    data = line.rstrip().split()
                    genes = data[2].split(',')
                    out_genes.append(genes)
        return ome, out_genes
    except FileNotFoundError:
        return ome, []




def output_ann_fas(ann_dir, ome, ome_genes, prot_path):
    faa = fa2dict(prot_path)
    genes = list(chain(*ome_genes))
    try:
        out_faa = {k: faa[k] for k in genes}
        with open(ann_dir + ome + '.faa', 'w') as out:
            out.write(dict2fa(out_faa))
        return ome, True
    except KeyError:
        eprint('\t\tERROR: ' + ome + ' accessions modified from input',
               flush = True)
        return ome, False


def output_res_faas(genes_in, prot_paths, ann_dir, cpus = 1):
    with mp.Pool(processes = cpus) as pool:
        res = pool.starmap(output_ann_fas, 
                           tqdm(((ann_dir, ome, ome_genes, prot_paths[ome]) \
                                 for ome, ome_genes in genes_in.items()),
                                total = len(genes_in)))

    return set(ome for ome, passing in res if not passing)

def call_hmmsearch(ome, pfam, ann_dir, faa_file, threads = 4):
    out_file = ann_dir + ome + '_pfam.out'
    if not os.path.isfile(out_file):
        hmmsearch = subprocess.call([
            'hmmsearch', '--cpu', str(threads), 
            '--domtblout', out_file + '.tmp', pfam, faa_file
            ], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
            )
        os.rename(out_file + '.tmp', out_file)
    else:
        hmmsearch = 0
    return ome, hmmsearch


def run_hmmsearch(pfam, passing_omes, ann_dir, cpus, 
                  db = None, threads = 4):
    """Manage running hmmsearch for Pfam annotations. Default runs just the
    Pfam hits, inputting a db runs for the full genome"""

    with mp.Pool(processes = round(cpus/2)) as pool:
        if db:
            hmm_res = pool.starmap(call_hmmsearch, 
                                   tqdm(((ome, pfam, ann_dir,
                                     db[ome]['faa'], threads) \
                                     for ome in list(passing_omes)),
                                    total = len(passing_omes)))
        else:
            hmm_res = pool.starmap(call_hmmsearch, 
                               tqdm(((ome, pfam, ann_dir,
                                     f'{ann_dir}{ome}.faa', threads) \
                                     for ome in list(passing_omes)),
                                    total = len(passing_omes)))

    failed_omes = []
    for ome, exit_code in hmm_res:
        if exit_code:
            eprint('\t\tERROR: ' + ome + ' failed hmmsearch', flush = True)
            failed_omes.append(ome)

    return failed_omes

def call_iprscan(ome, iprscan, ann_dir, threads = 4):
    out_file = ann_dir + ome + '_ipr.out'
    faa_file = ann_dir + ome + '.faa'
    if not os.path.isfile(out_file):
        ipr_code = subprocess.call([
            iprscan, '-i', faa_file, '-cpu', str(threads), 
            '-o', out_file + '.tmp', '-f', 'TSV'
            ], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
            )
    else:
        ipr_code = 0
    return ome, ipr_code 


def run_ipr_scan(iprscan, passing_omes, ann_dir, cpus):

    threads = 4
    with mp.Pool(processes = round(cpus/2)) as pool:
        ipr_res = pool.starmap(call_iprscan, 
                               tqdm(((ome, iprscan, ann_dir, threads) \
                                     for ome in list(passing_omes)),
                                    total = len(passing_omes)))

    failed_omes = []
    for ome, exit_code in ipr_res:
        if exit_code:
            eprint('\t\tERROR: ' + ome + ' failed InterProScan', flush = True)
            failed_omes.append(ome)

    return failed_omes




def parse_hmm_res(ome, ann_dir, evalue, threshold, max_hits = None):
    ome_out = ann_dir + ome + '_pfam.out'
#    lineComp = re.compile(r'([^ ]+)')
    res = defaultdict(list)
    with open(ome_out, 'r') as raw:
        for line in raw:
            if not line.startswith('#'):
#                data = lineComp.findall(line.rstrip())
                data = line.rstrip().split()
                hit_perc = (int(data[18]) - int(float(data[17])))/int(data[5])
                if hit_perc > threshold:
                    if float(data[6]) < evalue: # and float(data[7]) > score:
                        # pfam, name, cov_perc, evalue, bit, algn_start, algn_end
                        res[data[0]].append(( 
                            data[4], data[3], hit_perc,
                            float(data[21]), int(data[17]), int(data[18])
                            ))
    ome_res = {}
    for gene, gene_hits in res.items():
        gene_hits = sorted(gene_hits, key = lambda x: x[4], reverse = True)
        if max_hits:
            gene_hits = gene_hits[:max_hits]
        todel, hits = [], set()
        for i, hit in enumerate(gene_hits):
            if hit[0] in hits:
                todel.append(i)
            hits.add(hit)
        for i in reversed(todel):
            del gene_hits[i]
        ome_res[gene] = gene_hits

    return ome, ome_res


def compile_hmm_res(ann_dir, passing_omes, evalue = 0.001, 
                    threshold = 0.5, cpus = 1, max_hits = None):
        
    with mp.Pool(processes = cpus) as pool:
        omebyome_res = pool.starmap(parse_hmm_res, 
                               tqdm(((ome, ann_dir, evalue, threshold) \
                                     for ome in list(passing_omes)), 
                                    total = len(passing_omes)))
    ome_res = {}
    for ome, res in omebyome_res:
        ome_res[ome] = res

    return ome_res


def parse_ipr_res(ome, ann_dir, evalue, threshold):
    # NOT BUILT YET
    ome_out = ann_dir + ome + '_ipr.out'
    lineComp, res = re.compile(r'([^ ]+)'), defaultdict(list)
    with open(ome_out, 'r') as raw:
        for line in raw:
            if not line.startswith('#'):
                data = lineComp.findall(line.rstrip())
                hit_perc = (int(data[15]) - int(float(data[14])))/int(data[5])
                if hit_perc > threshold:
                    if float(data[6]) < evalue: # and float(data[7]) > score:
                        # pfam, name, cov_perc, evalue, bit, algn_start, algn_end
                        res[data[0]].append(( 
                            data[4], data[3], hit_perc, float(data[6]), 
                            float(data[7]), int(data[16]), int(data[17])
                            ))
    ome_res = {}
    for gene, gene_hits in res.items():
        gene_hits = sorted(gene_hits, key = lambda x: x[3], reverse = True)
        todel, hits = [], set()
        for i, hit in enumerate(gene_hits):
            if hit[0] in hits:
                todel.append(i)
            hits.add(hit)
        for i in reversed(todel):
            del gene_hits[i]
        ome_res[gene] = gene_hits

    return ome, ome_res




def compile_ipr_res(ann_dir, passing_omes, evalue = 0.01, threshold = 0.5, cpus = 1):
        
    with mp.Pool(processes = cpus) as pool:
        omebyome_res = pool.starmap(parse_ipr_res, 
                               tqdm(((ome, ann_dir, evalue, threshold) \
                                     for ome in list(passing_omes)), 
                                    total = len(passing_omes)))
    ome_res = {}
    for ome, res in omebyome_res:
        ome_res[ome] = res

    return ome_res




def ann_mngr(genes_list, prot_paths, wrk_dir, pfam, ipr_scan,
             evalue = 0.0001, threshold = 0.5, cpus = 1):

    ann_res = None
    ann_dir = wrk_dir + 'ann/'
    if not os.path.isdir(ann_dir):
        os.mkdir(ann_dir)

    failed_omes = output_res_faas(genes_list, prot_paths, ann_dir, cpus = cpus)
    passing_omes = set(genes_list.keys()).difference(failed_omes)

    if pfam:
        print("\tHmmsearch'ing Pfam database", flush = True)
        hmm_failed_omes = run_hmmsearch(pfam, passing_omes, ann_dir, cpus)
        passing_omes = passing_omes.difference(set(hmm_failed_omes))
        failed_omes = failed_omes.union(set(hmm_failed_omes))
        print('\tParsing hmmsearch output', flush = True)
        ann_res = compile_hmm_res(ann_dir, passing_omes, evalue, threshold, cpus)
    elif ipr_scan:
        print('\tRunning InterProScan', flush = True)
        ipr_failed_omes = run_ipr_scan(ipr_scan, passing_omes, ann_dir, cpus)
        passing_omes = passing_omes.difference(set(ipr_failed_omes))
        failed_omes = failed_omes.union(set(ipr_failed_omes))
        print('\tParsing InterProScan output', flush = True)
        pass

    return ann_res, failed_omes


def annotate_clusters(genes_list, #gff_path, prot_path, 
                      ome, ome_dir, gene2hg, prefix = 'hlg',
                      ann_res = {}):

    if not os.path.isfile(f'{ome_dir}{ome}/{prefix}.tsv'):
        return

    with open(ome_dir + ome + f'/{prefix}.tsv.tmp', 'w') as out:
        out.write('#name\thgs\tgenes\tgcfs\tpfams\n')
        with open(ome_dir + ome + f'/{prefix}.tsv', 'r') as raw:
            for line in raw:
                if not line.startswith('#'):
                    clus_d = line.split()
                    genes = clus_d[2].split(',')
                    anns_p = [ann_res[x] if x in ann_res else None for x in genes]
                    anns = ['|'.join([y[1] if y else '' for y in x]) \
                            if x else '' for x in anns_p]
                    ann_str = ';'.join(anns)
                    out.write('\t'.join(clus_d + [ann_str]) + '\n')
    shutil.move(f'{ome_dir}{ome}/{prefix}.tsv.tmp', 
                f'{ome_dir}{ome}/{prefix}.tsv')

#    return ome#, svg_dict
#    return ome, clus_gffs, clus_fas, svg_dict


def output_hgxs(hgx2dist, hgx2omes, hgx2i, i2ome, out_dir):
    maxhgxd = max(hgx2dist.values())
    minhgxd = min(hgx2dist.values())

    denom = maxhgxd - minhgxd
    hgx_output = []
    for hgx, omes in hgx2omes.items():
        hgx_id = hgx2i[hgx]
        hgx_output.append([
            ','.join([str(x) for x in hgx]), hgx_id,
            (log(hgx2dist[hgx]) - minhgxd)/denom,
            ','.join([i2ome[x] for x in omes])
            ]) # HGxs at this stage are not segregated into groups
    hgx_output = sorted(hgx_output, key = lambda x: x[3], reverse = True)
    with gzip.open(out_dir + 'hgxs.tsv.gz', 'wt') as out:
        out.write('#hgs\thgx_id\tnrm_log_tmd\tomes')#\tpdd\tomes')
        for entry in hgx_output:
            out.write('\n' + '\t'.join([str(x) for x in entry]))

def write_hlgs_txt_wpos(hlg_hgxs, hlg_omes, logg2d, 
                        hlg2clan, omes2patch, 
                        hgx2omes2gcl, hgx2omes2id, hgx2omes2pos, i2ome,
                        out_dir, group_name):
    hlg_output = []
    for i, omesc in hlg_omes.items():
        hlg_hgx = hlg_hgxs[i]
        try:
            hlg_output.append([
                ','.join([str(x) for x in hlg_hgx]), i, hlg2clan[i],
                logg2d[omesc], omes2patch[omesc], 
                hgx2omes2gcl[hlg_hgx][omesc],
                hgx2omes2id[hlg_hgx][omesc], hgx2omes2pos[hlg_hgx][omesc],
                100 * (1 - (hgx2omes2id[hlg_hgx][omesc]/hgx2omes2pos[hlg_hgx][omesc])), 
                ','.join([str(i2ome[x]) for x in omesc])#,
         #        dnds_dict[hgx][0], dnds_dict[hgx][1], str(dnds_dict[hgx][2]),
                ])
        except ZeroDivisionError:
            hlg_output.append([
                ','.join([str(x) for x in hlg_hgx]), i, hlg2clan[i],
                logg2d[omesc], omes2patch[omesc], 
                hgx2omes2gcl[hlg_hgx][omesc],
                hgx2omes2id[hlg_hgx][omesc], hgx2omes2pos[hlg_hgx][omesc],
                0, ','.join([str(i2ome[x]) for x in omesc])#,
         #        dnds_dict[hgx][0], dnds_dict[hgx][1], str(dnds_dict[hgx][2]),
                ])
        except KeyError:
            try:
                hlg_output.append([
                   ','.join([str(x) for x in hlg_hgx]), i, hlg2clan[i],
                   logg2d[omesc], omes2patch[omesc], 
                   hgx2omes2gcl[hlg_hgx][omesc],
                   hgx2omes2id[hlg_hgx][omesc], 'na', 'na',
                   ','.join([str(i2ome[x]) for x in omesc])#,
            #        dnds_dict[hgx][0], dnds_dict[hgx][1], str(dnds_dict[hgx][2]),
                   ])
            except KeyError:
                eprint(f'ERROR: missing HGx/omes from results: {hlg_hgx}, {omesc}', flush =  True)

    hlg_output = sorted(hlg_output, key = lambda x: x[3], reverse = True)
    with gzip.open(out_dir + f'{group_name}s.tsv.gz', 'wt') as out:
        out.write(f'#hgs\t{group_name}\thlc\tnrm_log_tmd\tpds' \
                + '\tgcl\tmmi\tmmp\tcsb\tomes') #\tmmp\tomes') #+ \
            #'selection_coef\tmean_dnds\tog_dnds\t' + \
         #   'total_dist'
#            )
        for entry in hlg_output:
            out.write('\n' + '\t'.join([str(x) for x in entry]))
    return hlg_output

def write_hlgs_txt_wopos(hlg_hgxs, hlg_omes, logg2d, 
                        hlg2clan, omes2patch, 
                        hgx2omes2gcl, hgx2omes2id, i2ome,
                        out_dir, group_name):
    hlg_output = []
    for i, hlg_hgx in hlg_hgxs.items():
        omesc = hlg_omes[i]
        hlg_output.append([
            ','.join([str(x) for x in hlg_hgx]), i, hlg2clan[i],
            logg2d[omesc], omes2patch[omesc], 
            hgx2omes2gcl[hlg_hgx][omesc],
            hgx2omes2id[hlg_hgx][omesc],
            ','.join([str(i2ome[x]) for x in omesc])#,
     #        dnds_dict[hgx][0], dnds_dict[hgx][1], str(dnds_dict[hgx][2]),
            ])

    hlg_output = sorted(hlg_output, key = lambda x: x[3], reverse = True)
    with gzip.open(out_dir + f'{group_name}s.tsv.gz', 'wt') as out:
        out.write(f'#hgs\t{group_name}\thlc\tnrm_log_tmd\tpds' \
                + '\tgcl\tmmi\tomes') #\tmmp\tomes') #+ \
            #'selection_coef\tmean_dnds\tog_dnds\t' + \
         #   'total_dist'
#            )
        for entry in hlg_output:
            out.write('\n' + '\t'.join([str(x) for x in entry]))
    return hlg_output


def calc_logg2d(hlg_omes, omes2dist):
    logg2d_prep = {}
    for omesc in hlg_omes.values():
        logg2d_prep[omesc] = log(omes2dist[omesc])
    maxhlgd = max(logg2d_prep.values()) # max observed, even of those truncated/removed
    minhlgd = min(logg2d_prep.values())
    denom = maxhlgd - minhlgd
    logg2d = {k: (v - minhlgd)/denom for k, v in logg2d_prep.items()}
    return logg2d


def threshold_hlg_by_dist(logg2d, dist_thresh, hlgs, hlg_hgxs, hlg_omes, omes2dist):
    if dist_thresh:
        gcfs, gcf_omes, gcf_hgxs = {}, {}, {}
        for hlg, locs in hlgs.items():
            hgxc = hlg_hgxs[hlg]
            omesc = hlg_omes[hlg]
            if logg2d[omesc] >= dist_thresh:
                gcfs[hlg] = locs
                gcf_omes[hlg] = omesc
                gcf_hgxs[hlg] = hgxc
        return gcfs, gcf_omes, gcf_hgxs
    else:
        return hlgs, hlg_omes, hlg_hgxs


def threshold_hlg_wpos(gcl_thresh, patch_thresh, id_perc, pos_perc, csb_thresh, hlgs,
                       hlg_hgxs, hlg_omes, omes2patch, hgx2omes2gcl,
                       hgx2omes2id, hgx2omes2pos):

    if any(x > 0 for x in [gcl_thresh, patch_thresh, id_perc, pos_perc, csb_thresh]):
        gcfs, gcf_omes, gcf_hgxs = {}, {}, {}
        for hlg, locs in hlgs.items():
            check = False # have we added a new list
            hgxc = hlg_hgxs[hlg]
            omesc = hlg_omes[hlg]
            try:
                patch = omes2patch[omesc]
                gcl = hgx2omes2gcl[hgxc][omesc]
                id_, pos = hgx2omes2id[hgxc][omesc], hgx2omes2pos[hgxc][omesc]
                if pos != 0 and id_ != 0:
                    csb = (1 - (id_/pos))
                else:
                    csb = 0
                if gcl >= gcl_thresh \
                    and patch >= patch_thresh \
                    and id_ >= id_perc \
                    and pos >= pos_perc \
                    and csb >= csb_thresh:
                    gcfs[hlg] = locs
                    gcf_omes[hlg] = omesc
                    gcf_hgxs[hlg] = hgxc
            except TypeError:
                 if hgx2omes2gcl[hgxc][omesc] >= gcl_thresh \
                    and omes2patch[omesc] >= patch_thresh \
                    and hgx2omes2id[hgxc][omesc] >= id_perc:
                    gcfs[hlg] = locs
                    gcf_omes[hlg] = omesc
                    gcf_hgxs[hlg] = hgxc
            except KeyError:
                eprint(f'\tWARNING: missing from proxies and removed: {omesc, hgxc}')
        return gcfs, gcf_omes, gcf_hgxs
    else:
        return hlgs, hlg_omes, hlg_hgxs

def threshold_hlg_wopos(gcl_thresh, patch_thresh, id_perc, hlgs,
                       hlg_hgxs, hlg_omes, omes2patch,
                       hgx2omes2gcl, hgx2omes2id):
    if any(x > 0 for x in [gcl_thresh, patch_thresh, id_perc]):
        print('\tApplying thresholds', flush = True)
        print('\t\t' + str(len(hlgs)) + ' HLGs before', flush = True)
        gcfs, gcf_omes, gcf_hgxs = {}, {}, {}
        for hlg, locs in hlgs.items():
            hgxc = hlg_hgxs[hlg]
            omesc = hlg_omes[hlg]
            if hgx2omes2gcl[hgxc][omesc] >= gcl_thresh \
                and omes2patch[omesc] >= patch_thresh \
                and hgx2omes2id[hgxc][omesc] >= id_perc:
                gcfs[hlg] = locs
                gcf_omes[hlg] = omesc
                gcf_hgxs[hlg] = hgxc

        return gcfs, gcf_omes, gcf_hgxs
    else:
        return hlgs, hlg_omes, hlg_hgxs


def output_thresholds(logg2d, dist_thresh, gcl_thresh, patch_thresh, id_perc, pos_perc,
                      csb_thresh, hlgs, hlg_hgxs, hlg_omes, omes2dist, omes2patch, hlg2clan,
                      hgx2omes2gcl, hgx2omes2id, hgx2omes2pos, out_dir, i2ome):
    print(f'\tFiltering HLGs for GCFs', flush = True)
    print('\tApplying thresholds', flush = True)
    print('\t\t' + str(len(hlgs)) + ' HLGs before', flush = True)
    hlgs, hlg_omes, hlg_hgxs = threshold_hlg_by_dist(logg2d, dist_thresh, hlgs, hlg_hgxs,
                                                     hlg_omes, omes2dist)
    if hgx2omes2pos:
        gcfs, gcf_omes, gcf_hgxs = threshold_hlg_wpos(gcl_thresh, patch_thresh, id_perc, pos_perc,
                                          csb_thresh, hlgs, hlg_hgxs, hlg_omes, omes2patch,
                                          hgx2omes2gcl, hgx2omes2id, hgx2omes2pos)
        gcf_output = write_hlgs_txt_wpos(gcf_hgxs, gcf_omes, logg2d, 
                            hlg2clan, omes2patch, 
                            hgx2omes2gcl, hgx2omes2id, hgx2omes2pos, i2ome,
                            out_dir, 'gcf')

    else:
        gcfs, gcf_omes, gcf_hgxs = threshold_gcf_wopos(gcl_thresh, patch_thres, 
                                           id_perc, gcfs,
                                           gcf_hgxs, gcf_omes, omes2patch,
                                           hgx2omes2gcl, hgx2omes2id)
        gcf_output = write_hlgs_txt_wopos(gcf_hgxs, gcf_omes, logg2d, 
                            hlg2clan, omes2patch, 
                            hgx2omes2gcl, hgx2omes2id, i2ome,
                            out_dir, 'gcf')


    print(f'\t\t{len(gcfs)} GCFs pass', flush = True)

    return gcf_output, gcfs, gcf_omes, gcf_hgxs


def output_figures(hlg_output, prefix, out_dir):
    # legacy, need to remove unless dn/ds is reimplemented
#    if dnds_dict:
 #       axes = [[],[],[],[],[], []]
  #  else:
    axes = [[],[],[], [], [], []]
    labels = []

    try:
        for entry in hlg_output:
            labels.append(
                ['HLG:', str(entry[1]) + ' | Omes: ' + entry[-1]]
                )
            axes[0].append(entry[5])
            axes[1].append(entry[4])
            axes[2].append(entry[3])
            axes[3].append(entry[6])
            axes[4].append(entry[7])
            axes[5].append(entry[8])
    except ValueError: # no positives
        labels = []
        axes = [[], [], [], []]
        for entry in hlg_output:
            labels.append(
                ['HLG:', str(entry[1]) + ' | Omes: ' + entry[-1]]
                )
            axes[0].append(entry[5])
            axes[1].append(entry[4])
            axes[2].append(entry[3])
            axes[3].append(entry[6])
            
    print(f'\tOutputting {prefix.upper()} plots', flush = True)
    if len(axes) == 4:
        axes_labels = ['PDS', 'GCL', 'nl_TMD',
                       'MMI']
    else:
        axes_labels = ['PDS', 'GCL', 
                       'nl_TMD', 'MMI', 
                       'MMP', 'CSB']

    fig = mk_subplots(labels, axes, axes_labels, alpha = 0.6)
    fig.write_html(out_dir + f'{prefix}_proxies.html')
    fig = mk_3d(labels, axes[:3], axes_labels[:3], alpha = 0.7)
    fig.write_html(out_dir + f'{prefix}s.html')


def quick_thresh(ome, ome_dir, hlgs):
    hlg_f = ome_dir + ome + '/hlg.tsv'
    gcf_f = ome_dir + ome + '/gcf.tsv'
    with open(hlg_f, 'r') as raw, open(gcf_f, 'w') as out:
        for line in raw:
            if line.startswith('#'):
                out.write(line)
            else:
                hlg = int(line.rstrip().split()[3])
                if hlg in hlgs:
                    out.write(line)


def threshold_gcf_quick(db, out_dir, ome_dir, 
                        dist_thresh, gcl_thresh, patch_thresh,
                        id_perc, pos_perc, csb_thresh, cpus = 1):

    ome2hlgs = defaultdict(set)
    print('\tLoading HLGs file', flush = True)
    with gzip.open(out_dir + 'hlgs.tsv.gz', 'rt') as raw, \
        gzip.open(out_dir + 'gcfs.tsv.gz', 'wt') as out:
        for line in raw:
            if line.startswith('#'):
                out.write(line)
            else:
                d = line.rstrip().split()
                hlg = int(d[1])
                tmd = float(d[3])
                pds = float(d[4])
                gcl = float(d[5])
                mmi = float(d[6])
                mmp = float(d[7])
                csb = float(d[8])/100
                omes = d[9].split(',')
                if tmd >= dist_thresh and gcl >= gcl_thresh \
                    and pds >= patch_thresh and mmi >= id_perc \
                    and mmp >= pos_perc and csb >= csb_thresh:
                    for ome in omes:
                        ome2hlgs[ome].add(hlg)
                    out.write(line)

    print('\tExtracting GCFs')    
    with mp.Pool(processes = cpus) as pool:
        pool.starmap(quick_thresh, ((ome, ome_dir, hlgs) \
                                    for ome, hlgs in ome2hlgs.items()))                        


def threshold_gcf_bypass(db, out_dir, wrk_dir, i2ome, gene2hg,
                         dist_thresh, gcl_thresh, patch_thresh,
                         id_perc, pos_perc, csb_thresh, ipr_path,
                         pfam_path, cpus = 1):

    print('\tLoading data structures', flush = True)
    with open(wrk_dir + 'hlgs.pickle', 'rb') as raw:
        hlgs, hlg_omes, hlg_hgxs, hlg2clan = pickle.load(raw)

    with open(wrk_dir + 'pds.full.pickle', 'rb') as in_pick:
        omes2patch = pickle.load(in_pick)

    with open(wrk_dir + 'gcl.pickle', 'rb') as pickin:
        hgx2omes2gcl = pickle.load(pickin)
    with open(wrk_dir + 'mmi.pickle', 'rb') as pickin:
        hgx2omes2id = pickle.load(pickin)
    if os.path.isfile(wrk_dir + 'mmp.pickle'):
        with open(wrk_dir + 'mmp.pickle', 'rb') as pickin:
            hgx2omes2pos = pickle.load(pickin)
    else:
        hgx2omes2pos = None

    with open(wrk_dir + 'omes2dist.pickle', 'rb') as raw:
        omes2dist = pickle.load(raw)

    logg2d = calc_logg2d(hlg_omes, omes2dist)
    gcf_output, gcfs, gcf_omes, gcf_hgxs = output_thresholds(logg2d, dist_thresh, 
                      gcl_thresh, patch_thresh, id_perc, pos_perc, csb_thresh,
                      hlgs, hlg_hgxs, hlg_omes, omes2dist, omes2patch, hlg2clan,
                      hgx2omes2gcl, hgx2omes2id, hgx2omes2pos, out_dir, i2ome)
    output_figures(gcf_output, 'gcf', out_dir)

    
    with mp.Pool(processes = cpus) as pool:
        hlg_gene_res = pool.starmap(read_ome_hlgs, ((ome, f'{out_dir}ome/hlg.tsv') \
                                                    for ome in i2ome))
    hlg_genes = {}
    for ome, ome_hlg_gene in hlg_gene_res:
        if ome_hlg_gene:
            hlg_genes[ome] = ome_hlg_gene

    gcf_genes = write_ome_output('gcf', gcfs, out_dir, db, cpus, gene2hg)
    annotation_mngr(hlg_genes, db, wrk_dir, pfam_path, ipr_path, 
                    gene2hg, out_dir, prefix = 'hlg', cpus = cpus)
    annotation_mngr(gcf_genes, db, wrk_dir, pfam_path, ipr_path, 
                    gene2hg, out_dir, prefix = 'gcf', cpus = cpus)

def write_ome_output(prefix, hlgs, out_dir, db, cpus, gene2hg):
    print(f'\tWriting {prefix.upper()} files', flush = True)
    write_clus_cmds = []
    ome2clusters = defaultdict(list)
    for hlg, loci in hlgs.items():
        for loc in loci:
            ome = loc[0][:loc[0].find('_')]
            ome2clusters[ome].append([loc, hlg])

    ome_dir = out_dir + 'ome/'
    if not os.path.isdir(ome_dir):
        os.mkdir(ome_dir)
    for ome, clusters in ome2clusters.items():
        if not os.path.isdir(ome_dir + ome):
            os.mkdir(ome_dir + ome)
        gff = db[ome]['gff3']
        out_file = ome_dir + ome + f'/{prefix}.tsv'
        write_clus_cmds.append([clusters, ome, out_file, gff, gene2hg])
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        out_genes = pool.starmap(write_clusters, 
                                 tqdm(write_clus_cmds, 
                                      total = len(write_clus_cmds)))
        pool.close()
        pool.join()
    genes = {x[0]: x[1] for x in out_genes if x[1]}
    return genes


def annotation_mngr(hlg_genes, db, wrk_dir, pfam_path, ipr_path, 
                    gene2hg, out_dir, prefix = 'hlg', cpus = 1):
    print('\tSearching', flush = True)
    prot_paths = {}
    for ome in hlg_genes:
        prot_paths[ome] = db[ome]['faa']

    ann_res, failedOmes = ann_mngr(
        hlg_genes, prot_paths, wrk_dir, pfam_path, ipr_path,
        evalue = 0.0001, threshold = 0.5, cpus = cpus
        )
    print('\tAnnotating', flush = True)
    annotate_clusters_cmds = []
    for ome, genes in hlg_genes.items():
        if ome in failedOmes:
            continue
        gff_path = db[ome]['gff3']
        pro_path = db[ome]['faa']
        try:
            annotate_clusters_cmds.append([genes, #gff_path, pro_path, 
                                  ome, out_dir + 'ome/', gene2hg, 
                                  prefix, ann_res[ome]])
        except KeyError:
            pass
#             annotate_clusters_cmds.append([genes, #gff_path, pro_path, 
 #                                  ome, out_dir + 'ome/', gene2hg,
   #                                prefix])


    with mp.get_context('fork').Pool(processes = cpus) as pool:
        pool.starmap(annotate_clusters, 
                             tqdm(annotate_clusters_cmds, 
                                  total = len(annotate_clusters_cmds)))
        pool.close()
        pool.join()


def output_hlgs(db, wrk_dir, hlgs, hlg_omes, i2ome, out_dir, hlg_hgxs,
         omes2dist, omes2patch, hgx2omes2gcl, hgx2omes2id, hgx2omes2pos, 
         gene2hg, plusminus, ome2i, hlg2clan, dist_thresh, gcl_thresh,
         patch_thresh, id_perc, pos_perc, csb_thresh, ipr_path = None,
         pfam_path = None, dnds_dict = {}, cpus = 1):

    print('\tWriting cluster scores', flush = True)

    logg2d = calc_logg2d(hlg_omes, omes2dist)
    if hgx2omes2pos:
        hlg_output = write_hlgs_txt_wpos(hlg_hgxs, hlg_omes, logg2d, 
                            hlg2clan, omes2patch, 
                            hgx2omes2gcl, hgx2omes2id, hgx2omes2pos, i2ome,
                            out_dir, 'hlg')
    else:
        hlg_output = write_hlgs_txt_wopos(hlg_hgxs, hlg_omes, logg2d, 
                            hlg2clan, omes2patch, 
                            hgx2omes2gcl, hgx2omes2id, i2ome,
                            out_dir, 'hlg')

    if dist_thresh or gcl_thresh or patch_thresh or id_perc or pos_perc:
        gcf_output, gcfs, gcf_omes, gcf_hgxs = output_thresholds(logg2d, 
                                   dist_thresh, gcl_thresh, patch_thresh, 
                                   id_perc, pos_perc, csb_thresh, hlgs, hlg_hgxs, 
                                   hlg_omes, omes2dist, omes2patch, hlg2clan,
                                   hgx2omes2gcl, hgx2omes2id, hgx2omes2pos, 
                                   out_dir, i2ome)
    else:
        gcf_output, gcfs, gcf_omes, gcf_hgxs = None, hlgs, hlg_omes, hlg_hgxs

    output_figures(hlg_output, 'hlg', out_dir)
    if gcf_output:
        output_figures(gcf_output, 'gcf', out_dir)

    prefix = 'hlg'
    hlg_genes = write_ome_output('hlg', hlgs, out_dir, db, cpus, gene2hg)
    if gcf_output:
#        prefix = 'gcf'
        gcf_genes = write_ome_output('gcf', gcfs, out_dir, db, cpus, gene2hg)
   
    if pfam_path or ipr_path:    
        print('\nX. Annotating clusters', flush = True)

    annotation_mngr(hlg_genes, db, wrk_dir, pfam_path, ipr_path, 
                    gene2hg, out_dir, prefix = 'hlg', cpus = cpus)
    if dist_thresh or gcl_thresh or patch_thresh or id_perc or pos_perc:
        annotation_mngr(hlg_genes, db, wrk_dir, pfam_path, ipr_path, 
                        gene2hg, out_dir, prefix = 'gcf', cpus = cpus)
