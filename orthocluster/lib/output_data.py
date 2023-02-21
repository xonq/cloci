import os
import re
import gzip
import shutil
import subprocess
import multiprocessing as mp
import plotly.graph_objects as go
from math import log
from tqdm import tqdm
from collections import defaultdict
from plotly.subplots import make_subplots
from mycotools.acc2gff import grabGffAcc
from mycotools.lib.kontools import eprint
from mycotools.lib.biotools import list2gff, gff2list, fa2dict, dict2fa
from orthocluster.orthocluster.lib import input_parsing

# NEED locID output and clans
# NEED to fix and optimize Pfam annotation
# NEED gbk output option
# NEED to add merge abutting clusters feature

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


def clus2scaf(cds_dict, clusters, gene2hg):

    gene2scaf, index_map = {}, {}
    for scaf, genes in cds_dict.items():
        index_map[scaf] = {v: i for i, v in enumerate(genes)}
        for gene in genes:
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
        
        clus_out[scaf].append([ord_loc, {gcf}])

    return clus_out 


def write_clusters(clusters, ome, out_file, gff_path, gene2hg, clusplusminus = 10):
    cds_dict = input_parsing.compileCDS(gff2list(gff_path), 
                                        os.path.basename(gff_path).replace('.gff3',''))

    clus_out = clus2scaf(cds_dict, clusters, gene2hg)
    change = True
    while change: # so long as there is a single merge
        clus_out, change = merge_clusters(clus_out, cds_dict) # attempt to merge clusters
#    clus_out = sorted(clus_out, key = lambda x: len(x[0]), reverse = True)
    # sort the clusters by the size of the OG combination tuple highest to
    # lowest

    for scaf, clusters in clus_out.items():
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


def setupHmmsearch(genes_in, prot_paths, hmm_dir):
    
    genes_out, todel = {}, []
    for ome in genes_in:
        prot_dict = fa2dict(prot_paths[ome])
        try:
            for sub_list in genes_in[ome]:
                for gene in sub_list: 
                    genes_out[gene] = prot_dict[gene]
        except KeyError:
            eprint('\t\tERROR: ' + ome + ' accessions modified from input', flush = True)
            todel.append(ome)
            continue
    
    genes = list(genes_out.keys())
    for ome in todel:
        omeKeys = [x for x in genes if x.startswith(ome + '_')]
        for key in omeKeys:
            del genes_out[key]
    
    with open(hmm_dir + 'complete.fa', 'w') as out:
        out.write(dict2fa(genes_out))
    
    return list(genes_out.keys()), hmm_dir + 'complete.fa', todel


def runHmmsearch(pfam, fa, hmm_dir, cpus):
    hmmsearch = subprocess.call([
        'hmmsearch', '--cpu', str(cpus), '--domtblout', hmm_dir + 'pfam.tmp', pfam, fa
        ], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
        )
    if hmmsearch:
        print('\t\tERROR: hmmsearch failed: ' + str(hmmsearch), flush = True)
    else:
       shutil.move(hmm_dir + 'pfam.tmp', hmm_dir + 'pfam.out')


def parseHmmRes(hmm_out, evalue = 0.01, threshold = 0.5):
    
    lineComp, res = re.compile(r'([^ ]+)'), defaultdict(list)
    with open(hmm_out, 'r') as raw:
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
    
    ome_res = defaultdict(dict)
    for gene in res:
        res[gene] = sorted(res[gene], key = lambda x: x[3], reverse = True)
        todel, hits = [], set()
        for i, hit in enumerate(res[gene]):
            if hit[0] in hits:
                todel.append(i)
            hits.add(hit)
        for i in reversed(todel):
            del res[gene][i]
        ome = gene[:gene.find('_')]
        ome_res[ome][gene] = res[gene]

    return ome_res


def pfamMngr(genes_list, prot_paths, wrk_dir, pfam, evalue = 0.01, threshold = 0.5, cpus = 1):

    hmm_dir = wrk_dir + 'hmm/'
    if not os.path.isdir(hmm_dir):
        os.mkdir(hmm_dir)

    genes, hmm_fa, failedOmes = setupHmmsearch(genes_list, prot_paths, hmm_dir)
    if not os.path.isfile(hmm_dir + 'pfam.out'):
        print("\tHmmsearch'ing Pfam database", flush = True)
        runHmmsearch(pfam, hmm_fa, hmm_dir, cpus)
    print('\tParsing hmmsearch output', flush = True)
    hmm_res = parseHmmRes(hmm_dir + 'pfam.out', evalue, threshold)

    return hmm_res, set(failedOmes)


def grabClus(genes_list, gff_path, prot_path, ome, ome_dir, gene2hg, pfamRes = {}):

    gff_dir, fa_dir = ome_dir + ome + '/gff/', ome_dir + ome + '/fa/'
    if not os.path.isdir(gff_dir):
        os.mkdir(gff_dir)
    if not os.path.isdir(fa_dir):
        os.mkdir(fa_dir)

    clus_gffs, clus_fas, clus_anns, scafs, svg_dict = [], [], {}, set(), {}
    gff_list, prot_dict = gff2list(gff_path), fa2dict(prot_path)
    print('\t' + ome, flush = True)
    for genes in genes_list:
        clus_gffs.append([[]])
        clus_fas.append({})
        clus_ann = ''
        for gene in genes:
            geneCoord = None
            geneGff = grabGffAcc(gff_list, gene)
            for entry in geneGff:
                if entry['type'].lower() == 'gene':
                    geneCoord = (int(entry['start']), int(entry['end']), entry['strand'])
            geneFa = prot_dict[gene]
            if gene in pfamRes:
                pfamStr = ';Pfam=' + '|'.join([
                    (hit[0] + '-' + hit[1]).replace('|','&') for hit in pfamRes[gene]
                    ])
                clus_ann += pfamStr[6:]
            else:
                pfamStr = ''
            clus_ann += ','
            try:
                ogStr = ';OG=' + str(gene2hg[gene])
            except KeyError: # gene not in an OG or is a singleton
                ogStr = ';OG='
            for entry in geneGff:
                if entry['type'] == 'gene' or entry['type'].lower() == 'mrna':
                    entry['attributes'] += pfamStr + ogStr
            geneFa['description'] = pfamStr[6:]
            clus_gffs[-1][0].extend(geneGff)
            clus_fas[-1][gene] = geneFa
        scaf_prep = clus_gffs[-1][0][0]['seqid']
        scaf = scaf_prep[:20]
        if scaf_prep != scaf:
            scaf = scaf_prep[:20] + '..'
        count = 0
        while scaf + '_' + str(count) in scafs:
            count += 1
        name = scaf + '_' + str(count)
        scafs.add(name)
        clus_anns[name] = clus_ann[:-1]
        clus_gffs[-1].append(name)
#        svg_dict[name] = copy.deepcopy(t_svg_dict)
        with open(gff_dir + name + '.gff3', 'w') as out:
            out.write(list2gff(clus_gffs[-1][0]))
        with open(fa_dir + name + '.fa', 'w') as out:
            out.write(dict2fa(clus_fas[-1]))

    with open(ome_dir + ome + '/info.out.tmp', 'w') as out:
        out.write('#name\thgs\tgenes\tgcfs\tpfams\n')
        with open(ome_dir + ome + '/info.out', 'r') as raw:
            for line in raw:
                if not line.startswith('#'):
                    clus = line.split('\t')[0]
                    out.write(line.rstrip() + '\t' + clus_anns[clus] + '\n')
    shutil.move(ome_dir + ome + '/info.out.tmp', ome_dir + ome + '/info.out')

    return ome#, svg_dict
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
        out.write('#hgs\thgx_id\tnrm_log_tmd\tomes')#\tpatchiness\tomes')
        for entry in hgx_output:
            out.write('\n' + '\t'.join([str(x) for x in entry]))


def output_gcfs(db, wrk_dir, gcfs, gcf_omes, i2ome, out_dir, gcf_hgxs,
         omes2dist, hgx2omes2gcc, omes2patch, hgx2omes2id,
         hgx2omes2pos, gene2hg, plusminus, ome2i,
         gcf2clan, pfam_path = None, dnds_dict = {}, dist_thresh = 0, cpus = 1):

    print('\tWriting cluster scores', flush = True)
    logg2d = {}
    for omesc in gcf_omes.values():
        logg2d[omesc] = log(omes2dist[omesc])
    maxgcfd = max(logg2d.values()) # max observed, even of those truncated/removed
    mingcfd = min(logg2d.values())
    denom = maxgcfd - mingcfd
    
    gcf_output = []
    for i, gcf_hgx in gcf_hgxs.items():
        omesc = gcf_omes[i]
        nml_log_tmd = (logg2d[omesc] - mingcfd)/denom
        if nml_log_tmd >= dist_thresh:
            gcf_output.append([
                ','.join([str(x) for x in gcf_hgx]), i, gcf2clan[i],
                nml_log_tmd, omes2patch[omesc], 
                hgx2omes2gcc[gcf_hgx][omesc],
                hgx2omes2id[gcf_hgx][omesc], hgx2omes2pos[gcf_hgx][omesc],
                ','.join([str(i2ome[x]) for x in omesc])#,
         #        dnds_dict[hgx][0], dnds_dict[hgx][1], str(dnds_dict[hgx][2]),
                ])

    gcf_output = sorted(gcf_output, key = lambda x: x[3], reverse = True)
    with gzip.open(out_dir + 'gcfs.tsv.gz', 'wt') as out:
        out.write('#hgs\tgcf\tclan\tnrm_log_tmd\tpatchiness' \
                + '\tgcc\tmmi\tmmp\tomes') #+ \
            #'selection_coef\tmean_dnds\tog_dnds\t' + \
         #   'total_dist'
#            )
        for entry in gcf_output:
            out.write('\n' + '\t'.join([str(x) for x in entry]))

    if dnds_dict:
        axes = [[],[],[],[],[]]
    else:
        axes = [[],[],[]]
    labels = []

    for entry in gcf_output:
        labels.append(
            ['GCF:', str(entry[1]) + ' | Omes: ' + entry[-1] + ' | HGs: ' + str(entry[0])]
            )
        axes[0].append(entry[5])
        axes[1].append(entry[4])
        axes[2].append(entry[3])
    print('\tOutputting scatter plots', flush = True)
    axes_labels = [
        'Distribution Patchiness', 'Gene Cluster Commitment', 'Log Microsynteny Distance' #,
        ]

    fig = mk_subplots(labels, axes, axes_labels, alpha = 0.6)
    fig.write_html(out_dir + 'metrics.html')
    fig = mk_3d(labels, axes[:3], axes_labels[:3], alpha = 0.7)
    fig.write_html(out_dir + 'gcfs.html')

    print('\tCompiling clusters from annotations', flush = True)
    write_clus_cmds = []
    ome2clusters = defaultdict(list)
    for gcf, loci in gcfs.items():
        for loc in loci:
            ome = loc[0][:loc[0].find('_')]
            ome2clusters[ome].append([loc, gcf])

    ome_dir = out_dir + 'ome/'
    if not os.path.isdir(ome_dir):
        os.mkdir(ome_dir)
    for ome, clusters in ome2clusters.items():
        if not os.path.isdir(ome_dir + ome):
            os.mkdir(ome_dir + ome)
        gff = db[ome]['gff3']
        out_file = ome_dir + ome + '/info.out'
        write_clus_cmds.append([clusters, ome, out_file, gff, gene2hg])
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        out_genes = pool.starmap(write_clusters, 
                                 tqdm(write_clus_cmds, 
                                      total = len(write_clus_cmds)))
        pool.close()
        pool.join()
    genes = {x[0]: x[1] for x in out_genes if x[1]}

#    hgx_dirTar.join()

    print('\tApplying Pfam annotations', flush = True)
    prot_paths = {}
    for ome in genes:
        prot_paths[ome] = db[ome]['faa']

    if pfam_path:
        pfamRes, failedOmes = pfamMngr(
            genes, prot_paths, wrk_dir, pfam_path, 
            evalue = 0.01, threshold = 0.5, cpus = cpus
            )
    else:
        pfamRes = {}
    

    print('\nX. Outputting clusters', flush = True)
    print('\tCluster biofiles', flush = True)
    grabClus_cmds = []
    for ome in genes:
        if ome in failedOmes:
            continue
        gff_path = db[ome]['gff3']
        pro_path = db[ome]['faa']
        try:
            grabClus_cmds.append([genes[ome], gff_path, pro_path, ome, ome_dir, gene2hg, pfamRes[ome]])
        except KeyError:
            grabClus_cmds.append([genes[ome], gff_path, pro_path, ome, ome_dir, gene2hg])


    with mp.get_context('fork').Pool(processes = cpus) as pool:
        clus_info = pool.starmap(grabClus, grabClus_cmds)
        pool.close()
        pool.join()
