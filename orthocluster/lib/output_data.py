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


def extrct_sig_clus(
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


def hash_clusters(ome, gff_path, ome_sig_clus, gene2hg, clusplusminus):
    """sliding window of size clusplusminus, identify sequences that may have
    hgs that may be significant, grab their window, and check for a significant
    higher order og combo in the locus"""
    
    gff_list, protoclus, clus_out = gff2list(gff_path), {}, []
    cds_dict = input_parsing.compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))
    
    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]):
            if seq0 in ome_sig_clus: # if the og is part of a significant seed
            # locus 
                for sig_clus in ome_sig_clus[seq0]:
                    start, end = None, None
                    if i0 < clusplusminus: # is i0 - clusplusminus < 0 ?
                        locus = cds_dict[scaf][:i0+clusplusminus+1] # then gather
                        # all the beginning
                    else:
                        locus = cds_dict[scaf][i0-clusplusminus:i0+clusplusminus+1]
                        # instead get the +/- and adjust for python
                    for i1, seq1 in enumerate(locus): # for each index and sequence
                    # in the locus
                        try:
                            og = gene2hg[seq1]
                        except KeyError:
                            continue
                        if og in sig_clus[0] and start is None: # if the og is
                        # in the sig clus and we haven't started
                            start = i1 # this is the start index
                        elif og in sig_clus[0]: # otherwise if it is in the
                        # sig clus label this the other border unless another is 
                        # found 
                            end = i1 + 1
                    hgx = tuple(sorted(list(sig_clus[0])))
                    clus_out.append([
                         hgx, locus[start:end], {sig_clus[1]}, scaf
                         ])
                    # clus_out = [(og0, og1, ... ogn), [prot0, prot1, ... protn])]
    
    return ome, clus_out, cds_dict


def merge_clusters(clus_out):
    """because we have significant clusters with borders, we can merge each
    cluster and assume that at a certain number of hgs in a combo will never
    happen via random sampling; therefore, their significance in the
    microsynteny adjustment is sufficient for justification for merging"""
    
    comp_clus, change = [], False
    while clus_out: # while there are clusters
        any_intersect, toDel = False, [0] # we haven't found intersection
        clus0 = clus_out[0] # grab the first tuple
        loc0 = set(clus0[1]) # grab the set of proteins in the cluster
        for i, clus1 in enumerate(clus_out[1:]): # look at all other clusters
            loc1 = set(clus1[1]) # grab their loci
            intersect = loc0.intersection(loc1) # check for overlap between
            # proteins
            if intersect:
                any_intersect, change = True, True # we found overlap in
                # clus_out and in this specific cluster comparison
                hgx = tuple(sorted(list(set(clus0[0]).union(set(clus1[0])))))
                # obtain the higher order og combination as the union of the
                # two cluster sets, then sort it and store as a tuple
                comp_clus.append([
                    hgx, list(loc0.union(loc1)), 
                    clus0[2].union(clus1[2]), clus0[3] #, merge_x
                    ])
                # append the og combo and the union of the two loci's proteins
                # protein order is random
                toDel.append(i) # we can delete the intersected locus as well
                # as the first
        
        if not any_intersect: # if there is not any intersection in the locus
            comp_clus.append(clus0) #, clus0[2]]) # simply append the original
            # result
        toDel.sort(reverse = True) # sort the todel from highest to lowest to
        # not mess up order when indexing
        for i in toDel:
            clus_out.pop(i)
    
    return comp_clus, change


def write_clusters(ome_sig_clus, ome, out_file, gff_path, gene2hg, clusplusminus = 10):
    # is multiprocessing this messing something and memory and that's why
    # output is wonky? (sometimes fewer proteins than should be in locus) 
    # or is output consistent fucking up elsewhere

    ome, clus_out, cds_dict = hash_clusters(ome, gff_path, ome_sig_clus, gene2hg, clusplusminus)
    change = True
    while change: # so long as there is a single merge
        clus_out, change = merge_clusters(clus_out) # attempt to merge clusters
    clus_out = sorted(clus_out, key = lambda x: len(x[0]), reverse = True)
    # sort the clusters by the size of the OG combination tuple highest to
    # lowest

    scafs = set()
    for clus in clus_out:
        clus[1] = sorted(clus[1], key = cds_dict[clus[3]].index)
        clusOGs = []
        for gene in clus[1]:
            try:
                clusOGs.append(gene2hg[gene])
            except KeyError: # gene not in an OG/is a singleton
                clusOGs.append('')
        clus[0] = clusOGs
        scaf = clus[3][:20]
        if clus[3] != scaf:
            scaf = clus[3][:20] + '..'
        count = 0
        while scaf + '_' + str(count) in scafs:
            count += 1
        name = scaf + '_' + str(count)
        scafs.add(name)
        clus.append(name)


    out_genes = []
    with open(out_file, 'w') as out:
        out.write('#name\thgs\tgenes\tgcfs\n')
        for clus in clus_out:
            clus[2] = sorted(clus[2])
            hgs = ','.join([str(x) for x in clus[0]]) # og1,og2,...ogn
            genes = ','.join(clus[1]) # prot1,prot2,prot3,...protn
            gcfs= ';'.join([str(x) for x in list(clus[2])])
            out.write(clus[4] + '\t' + hgs + '\t' + genes + '\t' + gcfs + '\n')
            out_genes.append(clus[1])

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


def output_res(db, wrk_dir, hgx2dist, gcfs, gcf_omes, i2ome, hgx2omes, out_dir, gcf_hgxs,
         omes2dist, hgx2omes2gcc, omes2patch, hgx2omes2id,
         hgx2omes2pos, hgx2loc, gene2hg, plusminus, ome2i,
         hgx2i, pfam_path = None, dnds_dict = {}, cpus = 1):
    print('\tWriting cluster scores', flush = True)
    omes2dist = {k: log(v) for k, v in omes2dist.items() if v}
    maxval = max(omes2dist.values()) # max observed, even of those truncated/removed
    minval = min(omes2dist.values())
    denom = maxval - minval
    gcf_output, done = [], set()
    for gcf, modHGx2omes in enumerate(gcfs):
        for hgx, omes in modHGx2omes.items():
            hgx_id = hgx2i[hgx]
            if hgx2dist[hgx] == 0:
                print(omes, '0 hgx2dist value', flush = True)
                continue
            if hgx_id not in done:
                done.add(hgx_id)
            else:
                continue
            gcf_output.append([
                ','.join([str(x) for x in hgx]), hgx_id,
                gcf, (log(hgx2dist[hgx]) - minval)/denom,
                ','.join([i2ome[x] for x in hgx2omes[hgx]])
                ]) # HGxs at this stage are not segregated into groups
    gcf_output = sorted(gcf_output, key = lambda x: x[3], reverse = True)
    with gzip.open(out_dir + 'hgxs.tsv.gz', 'wt') as out:
        out.write('#hgs\thgx_id\tgcf\tnrm_log_tmd\tomes')#\tpatchiness\tomes')
        for entry in gcf_output:
            out.write('\n' + '\t'.join([str(x) for x in entry]))

    
    kern_output, top_hgxs = [], []
    for i, ogc in enumerate(gcf_hgxs):
        omesc = gcf_omes[i]
        if omesc not in omes2patch:
            print(omesc, 'not in patch', flush = True)
            continue
        elif omes2dist[omesc] == 0:
            print(omesc, '0 value', flush = True)
            continue
        kern_output.append([
            ','.join([str(x) for x in ogc]), i,
            (omes2dist[omesc] - minval)/denom, omes2patch[omesc], 
            hgx2omes2gcc[ogc][omesc],
            hgx2omes2id[ogc][omesc], hgx2omes2pos[ogc][omesc],
            ','.join([str(i2ome[x]) for x in omesc])#,
     #        dnds_dict[hgx][0], dnds_dict[hgx][1], str(dnds_dict[hgx][2]),
           # omes2dist[omesc]
            ])
        for shgx, omes in gcfs[i].items():
            top_hgxs.append([shgx, i, set(omes)])

    kern_output = sorted(kern_output, key = lambda x: x[2], reverse = True)
    with gzip.open(out_dir + 'gcfs.tsv.gz', 'wt') as out:
        out.write('#hgs\tgcf\tnrm_log_tmd\tpatchiness' \
                + '\tgcc\tmmi\tmmp\tomes') #+ \
            #'selection_coef\tmean_dnds\tog_dnds\t' + \
         #   'total_dist'
#            )
        for entry in kern_output:
            out.write('\n' + '\t'.join([str(x) for x in entry]))

    if dnds_dict:
        axes = [[],[],[],[],[]]
    else:
        axes = [[],[],[]]
    labels = []

    for entry in kern_output:
        labels.append(
            ['GCF:', str(entry[1]) + ' | Omes: ' + entry[-1] + ' | HGs: ' + str(entry[0])]
            )
        axes[0].append(entry[4])
        axes[1].append(entry[3])
        axes[2].append(entry[2])
        if dnds_dict:
            axes[3].append(entry[6])
            axes[4].append(entry[7])
    print('\tOutputting scatter plots', flush = True)
    axes_labels = [
        'Distribution Patchiness', 'Gene BLAST Congruence', 'Log Microsynteny Distance' #,
#        'Selection Coefficient', 'Mean dN/dS'
        ]
    if dnds_dict:
        axes_labels.extend(['Selection Coefficient', 'Mean dN/dS'])
    fig = mk_subplots(labels, axes, axes_labels, alpha = 0.6)
    fig.write_html(out_dir + 'metrics.html')
    fig = mk_3d(labels, axes[:3], axes_labels[:3], alpha = 0.7)
    fig.write_html(out_dir + 'gcfs.html')

    print('\tCompiling clusters from annotations', flush = True)
    sig_clus = extrct_sig_clus(
        hgx2dist, hgx2loc, top_hgxs, ome2i
        )
    write_clus_cmds = []
    ome_dir = out_dir + 'ome/'
    if not os.path.isdir(ome_dir):
        os.mkdir(ome_dir)
    for ome in sig_clus:
        if not os.path.isdir(ome_dir + ome):
            os.mkdir(ome_dir + ome)
        gff = db[ome]['gff3']
        out_file = ome_dir + ome + '/info.out'
        write_clus_cmds.append([sig_clus[ome], ome, out_file, gff, gene2hg, plusminus])
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

