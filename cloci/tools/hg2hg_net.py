#! /usr/bin/env python3

import os
import re
import sys
import argparse
import numpy as np
import multiprocessing as mp
from math import log
from graph_tool import centrality, clustering
from graph_tool.all import *
from scipy.sparse import lil_matrix
from itertools import combinations, chain
from collections import Counter, defaultdict
from cloci.lib.hgx2hlgs import MCL, read_clus, read_mci_rows
from mycotools.lib.kontools import format_path, getColors, hex2rgb, \
    write_json, read_json, eprint, mkOutput
from mycotools.lib.dbtools import mtdb
from mycotools.lib.biotools import gff2list
from mycotools.utils.extractHmmAcc import grabAccs

# NEED to output json of the network modules
# NEED to organize output
def collectGffs(omes, top_dir):
    clusGFFs = [
        top_dir + ome + '/' + ome + '.gff3' for ome in omes \
        if os.path.isfile(top_dir + ome + '/' + ome + '.gff3')
        ]
    return clusGFFs

def mkPfam2name(pfam_file):
    # compile pfam 2 name
    with open(pfam_file, 'r') as raw:
        pfam_str = raw.read()
    pfams = grabAccs(pfam_str)
    pfam2name = {x[1].lstrip(): x[0].lstrip() for x in pfams}
    pfam2name = {k[:k.find('.')]: v for k,v in pfam2name.items()}
    return pfam2name


def extract_lineages(db, rank):
    # acquire all lineages that belong to rank within input database
    omesXlineage = {}
    for ome, row in db.set_index().items():
        lineage = row['taxonomy'][rank]
        if lineage not in omesXlineage:
            omesXlineage[lineage] = []
        omesXlineage[lineage].append(ome)
    return omesXlineage


def parse_info(info_path, ome):
    # acquire all homology_group pairs from each cluster, return ready for counter hash

    counts, hg_counts, hg2pfam = Counter(), Counter(), defaultdict(list)
    with open(info_path, 'r') as raw:
        for line in raw:
            todel = []
            if not line.startswith('#'):
                data = line.rstrip().split('\t')
                hgs, genes = data[1].split(','), data[2].split(',')
                if not len(data) == 5:
                    eprint(f'\tWARNING: skipping {ome} missing annotations', flush = True)
                    return ome, None, None, None
                pfams = [x.split('|') for x in data[-1].split(';')]
                gene2pfam = {k: pfams[i] for i, k in enumerate(genes)}
                for i, hg in enumerate(hgs):
                    try:
                        hgs[i] = int(hg)
                    except ValueError:
                        todel.append(i)
                if todel:
                    for i in reversed(todel):
                        del hgs[i]
                        del pfams[i]

                for i, hg in enumerate(hgs):
                    if pfams[i]:
                        hg2pfam[hg].extend(list(set(pfams[i])))
                sort_hgs = sorted(hgs)
                hg_counts += Counter(set(sort_hgs))
                pairs = combinations(sort_hgs, 2)
                counts += Counter(pairs)
                

    return ome, tuple([tuple([k, v]) for k,v in counts.items()]), \
        tuple([tuple([k,v]) for k,v in hg_counts.items()]), \
        tuple([tuple([k,tuple(v)]) for k,v in hg2pfam.items()])


def combine_counts(counts_res, db, rank = 'species'):
    # combine counts from multiprocessing
    db = db.set_index()

    counts = defaultdict(Counter)
    hg_counts = defaultdict(Counter)
    hg2pfam = defaultdict(dict)
    taxonLens = defaultdict(int)
    taxon_fails = defaultdict(int)
    for ome, tCounts, tHGcounts, tHG2Pfam in counts_res:
        taxon = db[ome]['taxonomy'][rank]
        if tCounts == None:
            taxon_fails[taxon] += 1
            continue
        taxonLens[taxon] += 1
        counts[taxon] += Counter({x[0]: x[1] for x in tCounts})
        hg_counts[taxon] += Counter({x[0]: x[1] for x in tHGcounts})
        for hg, pfams in tHG2Pfam:
            if hg not in hg2pfam[taxon]:
                hg2pfam[taxon][hg] = Counter()
            hg2pfam[taxon][hg] += Counter(pfams)

    fCounts, fHGcounts, fHG2Pfam = Counter(), Counter(), {}
    for taxon in counts:
        # average across the rank
        tCounts = {k: round(v/taxonLens[taxon] + 0.5) for k, v in counts[taxon].items()}
        tHGcounts = {k: round(v/taxonLens[taxon] + 0.5) for k, v in hg_counts[taxon].items()}
        tHG2Pfam = {}
        for hg, pfams in hg2pfam[taxon].items():
            tHG2Pfam[hg] = {k: round(v/taxonLens[taxon]) for k, v in pfams.items()}

        fCounts += Counter(tCounts)
        fHGcounts += Counter(tHGcounts)
        for hg, pfams in tHG2Pfam.items():
            if hg not in fHG2Pfam:
                fHG2Pfam[hg] = Counter()
            fHG2Pfam[hg] += pfams

    for hg in fHG2Pfam:
        fHG2Pfam[hg] = {k: v for k,v in sorted(fHG2Pfam[hg].items(), key = lambda x: x[1], reverse = True)}
            
    return fCounts, fHGcounts, fHG2Pfam, taxon_fails

def write_node_quants(hg_counts, out_path):

    hg_counts = [
        (k, v) for k,v in sorted(hg_counts.items(), key = lambda x: x[1],
        reverse = True)
        ]
    with open(out_path, 'w') as out:
        for k,v in hg_counts:
            out.write(str(k) + '\t' + str(v) + '\n')


def increase_float(val):
    # high fidelity float increase
    strVal = str(val)
    if strVal.endswith('9'):
        newVal = float(strVal + '1')
    else:
        factor = len(strVal[strVal.find('.')+1:])
        scale = 10**factor
        transVal = round(val*scale)
        newVal = float('0.' + str(transVal + 1))
    return newVal


def write_to_adj(oEdges, edges, out_path, oi2ni):
    minCount, maxCount = min(oEdges.values()), max(oEdges.values())
    denom = maxCount - minCount
    if denom == 0:
        eprint('\t\t\tWARNING: low connection weight throughout network', flush = True)
        denom = 1
    with open(out_path, 'w') as outF:
        for k, value in edges.items():
#            value = ((v-minCount)/denom) # normalize
            out_data = [oi2ni[k[0]], oi2ni[k[1]], value]
            outF.write('\t'.join([str(x) for x in out_data]) + '\n')



def write_adj_matr(oEdges, out_path, maximum = 500, duplicates = False):

    if not duplicates: # dont look at same HG-HG relationships
        oEdges = {k: v for k,v in oEdges.items() if k[0] != k[1]}

    data = sorted(list(oEdges.values()))
    edges = oEdges
    
#    write_to_adj(oEdges, edges, out_path + '.raw', oi2ni)

    percentile = 0.0
    while len(edges) > maximum: # while too many edges
        percentile = increase_float(percentile)
        minPerc = data[round((len(data)-1)*percentile)]
        edges = {k: v for k, v in oEdges.items() if v > minPerc}

    edge_ns = sorted(set(chain(*list(oEdges.keys()))))
    ni2oi = {int(k): int(v) for k, v in enumerate(edge_ns)}
    oi2ni = {v: k for k, v in ni2oi.items()}

    print('\t\tEdge #|percentile:', len(edges), percentile, flush = True)
    write_to_adj(oEdges, edges, out_path, oi2ni)
    
    return ni2oi, oi2ni    


def import_adj(adj_file, loclen, minimum = 0):
#    adj_arr = np.zeros([len(loci2grab), len(loci2grab)])
    adj_arr = lil_matrix((loclen, loclen))

    count = 0
    with open(adj_file, 'r') as raw:
        for line in raw:
            l0, l1, gid = line.rstrip().split()
            l0, l1 = int(l0), int(l1)
            val = float(gid)
            if val >= minimum:
                adj_arr[(l0, l1)] = val
                
    return adj_arr


def run_community_detection(wrk_dir, adj_f, inflation = 1.5, threads = 1):
    clus_f = adj_f + '.clus'
    rows_f = adj_f + '.rows'
    MCL(adj_f, clus_f, inflation, threads, rows_f)
    mci2pre = read_mci_rows(rows_f)
    clus = read_clus(clus_f, mci2pre)
    i2clus = {}
    for i,v in enumerate(clus):
        for i1 in v:
            i2clus[i1] = i
    return clus, i2clus


def make_network(
    wrk_dir, adj_file, quant_path, network_path, hg2pfam_path,
    ni2oi, oi2ni, calc_transitivity = False, calc_centrality = False,
    calc_local_transitivity = False,
    skip = False, overwrite = False, scale = 1.5, def_size = 1, annotate = False,
    inflation = 1.5, cpus = 1
    ):
    # open adjacency matrix and construct a network as an html


    clus, i2clus = run_community_detection(wrk_dir, adj_file, 
                                        inflation, cpus * 2 - 1)

    adj_arr = import_adj(adj_file, len(oi2ni))
    idx = adj_arr.nonzero()
    g = Graph(directed = False)
    g.add_edge_list(np.transpose(idx))
    weights = adj_arr[idx].toarray()
    reprop = g.new_edge_property('float')
    reprop.a = weights
    v = weights[:, 1]
    weights[:, 1] = log((v - v.min()) / (v.max() - v.min()) + def_size)
    eprop = g.new_edge_property('float')
    eprop.a = weights

    quants_p = {}
    with open(quant_path, 'r') as raw:
        for line in raw:
            k, v = line.rstrip().split('\t')
            if int(k) in oi2ni: # HGs with only duplicates may be excluded
                quants_p[oi2ni[int(k)]] = float(v)


    max_q, min_q = max(quants_p.values()), min(quants_p.values())
    denom = max_q - min_q
    quants = g.new_vertex_property('float')
    if denom:
        for k, v in quants_p.items():
            quants[k] = log((v - min_q)/denom + def_size) + def_size * 3
    else:
        for k, v in quants_p.items():
            quants[k] = 4

    hg2pfam_dict = {}
    hg2pfam = g.new_vertex_property('string')
    if os.path.isfile(hg2pfam_path):
        hg2pfamPrep = read_json(hg2pfam_path)
        for hg, pfams in hg2pfamPrep.items():
            hg = int(hg)
            if hg in oi2ni:
                if pfams:
                    v = list(pfams.keys())[0]
                    hg2pfam_dict[oi2ni[hg]] = f'{hg}|{v[v.find("-")+1:]}'
                    hg2pfam[oi2ni[hg]] = f'{hg}|{v[v.find("-")+1:]}'
                else:
                    hg2pfam[oi2ni[hg]] = str(hg)
                    hg2pfam_dict[oi2ni[hg]] = str(hg)
            
    else:
        for k, v in enumerate(quants):
            hg2pfam[k] = str(ni2oi[k])
            hg2pfam_dict[k] = str(ni2oi[k])


#    if len(quants) < 2:
 #       eprint('\t\tERROR: less than 2 nodes in network', flush = True)
  #      return
    count = 0
    v_fil_col = g.new_vertex_property('vector<float>')
    colors = getColors(len(clus), ignore = ['#ffffff'])
    for c, nodes in enumerate(clus):
        if len(nodes) > 1:
            try:
                color = colors[count]
                count += 1
            except IndexError:
                count = 0
                color = colors[count]
        else:
            color = '#ffffff'
        rgb = [x/255.0 for x in hex2rgb(color)] + [0.8]
        for hg in nodes:
#            if hg in oi2ni:
 #               ni = oi2ni[hg]
            v_fil_col[hg] = rgb
        

    g.vp['name'] = hg2pfam
    print('\t\tWriting', flush = True)
    if not skip:
        if overwrite or not os.path.isfile(network_path):
            if annotate:
                graph_draw(g, vertex_size = quants, vertex_text = g.vp['name'],
                   output = network_path, vertex_halo_color = v_fil_col, 
                   vertex_fill_color = v_fil_col,
                   edge_pen_width = eprop, edge_color = [0,0,0,1], vertex_font_family = 'arial',
                   vertex_font_size = 1)
            else:
                graph_draw(g, vertex_size = quants,
                   output = network_path, vertex_halo_color = v_fil_col, 
                   vertex_fill_color = v_fil_col,
                   edge_pen_width = eprop, edge_color = [0,0,0,1])
    
    print('\t\tNodes:', len(oi2ni), flush = True)
#    print('\t\tEdges:', len(edgeLabels), flush = True)


    print('\t\tQuantifying graph summary', flush = True)
    if calc_transitivity:
        print('\t\t\tQuantifying global graph transitivity', flush = True)
        print('\t\t\t\tNode:', flush = True)
        global_transitivity = clustering.global_clustering(g, reprop)
        stats = global_transitivity
        print(f'\t\t\t\t\t{global_transitivity[0]} stdev: {global_transitivity[1]}', 
              flush = True)
    else:
        stats = [None, None]

    density = len(weights[0])/((len(quants_p) - (len(quants_p) - 1))/2)
    stats = stats + (density,)
    print(f'\t\t\tDensity: {density}', flush = True)

    if calc_local_transitivity:
        local_transitivity = clustering.local_clustering(g, reprop)
        trans_path = adj_file[:-4] + '.transitivity'
        with open(trans_path, 'w') as out:
            out.write('#hg|pfam\ttransitivity\tnode_size\n')
            for k, v in sorted(enumerate(local_transitivity),
                               key = lambda x: x[1], reverse = True):
                if k in quants_p:
                    out.write(f'{hg2pfam_dict[k]}\t{v}\t{quants_p[k]}\n')
                elif k in hg2pfam_dict:
                    eprint(f'\t\t\tWARNING: {hg2pfam_dict[k]} missing from node sizes',
                           flush = True)
                else:
                    eprint(f'\t\t\tWARNING: rogue network node {k}', flush = True)
#        print('\t\t\t\tEdge:', flush = True)
 #       edge_transitivity = nx.edge_transitivity(G)
  #      print('\t\t\t\t\t' + str(edge_transitivity), flush = True)
    if calc_centrality:
        print('\t\t\tQuantifying betweenness centrality', flush = True)
        v_between, e_between = centrality.betweenness(g)
        central_path = adj_file[:-4] + '.centrality'
        with open(central_path, 'w') as out:
            out.write('#hg|pfam\tbetweenness\tnode_size\n')
            for k, v in sorted(enumerate(v_between),
                               key = lambda x: x[1], reverse = True):
                if v == 0:
                    continue
                if k in quants_p:
                    out.write(f'{hg2pfam_dict[k]}\t{v}\t{quants_p[k]}\n')
                elif k in hg2pfam_dict:
                    eprint(f'\t\t\tWARNING: {hg2pfam_dict[k]} missing from node sizes',
                           flush = True)
                else:
                    eprint(f'\t\t\tWARNING: rogue network node {k}', flush = True)

    return stats 

def main(
    db, rank, cloci_dir, transitivity = False, 
    calc_centrality = False, local_transitivity = False,
    max_edges = 500, hlg = 'hlg',
    overwrite = False, cpus = 1, skip = False
    ):

    if cloci_dir:
        ome_dir, net_dir = cloci_dir + 'ome/', cloci_dir + 'net/'
        if not os.path.isdir(ome_dir):
            eprint('\nERROR: invalid cloci output', flush = True)
            sys.exit(2)

    wrk_dir = net_dir + 'working/'
    if not os.path.isdir(net_dir):
        os.mkdir(net_dir)
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)
    
    omesXlineage = extract_lineages(db, rank)

    omesXlineage = {
        k: v for k,v in sorted(omesXlineage.items(), key = lambda x: len(x[1]))
        }
    lineage2stats = {}
    for lineage, omes in omesXlineage.items():
        print('\n' + lineage, flush = True)
        adj_path = wrk_dir + lineage + '.adj'
        quant_path = wrk_dir + lineage + '.tsv'
        hg2pfam_path = wrk_dir + lineage + '.hg2pfam.txt'
#        if not os.path.isfile(adj_path): # need to log the factor then
        print('\tAdjacency matrix', flush = True)
        info_cmds = [
            [f'{ome_dir}{ome}/{hlg}.tsv', ome] \
            for ome in omes \
            if os.path.isfile(f'{ome_dir}{ome}/{hlg}.tsv')
            ]
        if not info_cmds:
            eprint('\tWARNING: skipping; no output files found', flush = True)
            continue
        with mp.Pool(processes = cpus) as pool:
            counts_res = pool.starmap(parse_info, info_cmds)
        counts, hg_counts, hg2pfam, taxon_fails = combine_counts(
            counts_res, db, rank = 'species'
            )
        if not hg_counts:
            eprint('\tWARNING: skipping; no annotations', flush = True)
            continue
        write_json(hg2pfam, hg2pfam_path)
        write_node_quants(hg_counts, quant_path)
        ni2oi, oi2ni = write_adj_matr(counts, adj_path, max_edges)
    
        net_path = net_dir + lineage + '.svg'
        print('\tNetwork', flush = True)
        stats = make_network(wrk_dir, adj_path, quant_path, net_path, 
                         hg2pfam_path, ni2oi, oi2ni,
                         transitivity, calc_centrality, 
                         local_transitivity, skip,
                         overwrite)
        lineage2stats[lineage] = stats + (len(omes) - taxon_fails[lineage], taxon_fails[lineage],)

    with open(net_dir + rank + '.stats.tsv', 'w') as out:
        out.write('#lineage global_transitivity gl_tr_stdev density pass_genomes failed_genomes\n')
        for lineage, stats in lineage2stats.items():
            out.write(f'{lineage} {" ".join([str(x) for x in stats])}\n')
     
            

def cli():
    parser = argparse.ArgumentParser(
        description = 'Takes CLOCI output, creates networks of HG' + \
        ' relationships. Node sizes are normalized and edge sizes are ' + \
        'log-normalized during visualization'
        )
    parser.add_argument(
        '-i', '--input', required = True,
        help = 'CLOCI results directory'
        )
    parser.add_argument(
        '-g', '--gcfs', help = 'Make networks from GCFs; DEFAULT HLGs',
        action = 'store_true',
        )
    parser.add_argument(
        '-d', '--database', required = True, 
        help = 'Must be same as CLOCI input'
        )
    parser.add_argument(
        '-r', '--rank', help = 'Taxon rank; DEFAULT: class',
        default = 'class'
        )
    parser.add_argument('-s', '--skip', help = 'Skip writing network file', 
                        action = 'store_true')
    parser.add_argument(
        '-m', '--maximum', help = 'Find upper percentile of edges to fit ' + \
        'maximum edges',
        type = int
        )
    parser.add_argument(
        '--centrality', help = 'Calculate betweenness centrality',
        action = 'store_true'
        )
    parser.add_argument(
        '-gt', '--global_transitivity', help = 'Quantify global transitivity',
        action = 'store_true'
        )
    parser.add_argument('-lt', '--local_transitivity', help = 'Quantify local transitivity',
                        action = 'store_true')
    parser.add_argument(
        '-c', '--cpu', default = 1, type = int
        )
    parser.add_argument('--overwrite', help = 'Rerun if necessary',
        action = 'store_true')
    args = parser.parse_args()

    db = mtdb(format_path(args.database))
    cloci_dir = format_path(args.input)

    if args.gcfs:
        hlg = 'gcf'
    else:
        hlg = 'hlg'

    main(
        db, args.rank, cloci_dir, transitivity = args.global_transitivity, 
        calc_centrality = args.centrality, local_transitivity = args.local_transitivity,
        hlg = hlg,
        max_edges = args.maximum, cpus = args.cpu, overwrite = args.overwrite, 
        skip = args.skip
        )
    sys.exit(0)


if __name__ == '__main__':
    cli()
