#! /usr/bin/env python3

from collections import Counter, defaultdict
from itertools import combinations
import multiprocessing as mp, os, argparse, sys, re
import networkx as nx
from igraph import Graph
from pyvis.network import Network
from mycotools.lib.kontools import format_path, getColors, hex2rgb, \
    write_json, read_json, eprint, mkOutput
from mycotools.lib.dbtools import mtdb
from mycotools.lib.biotools import gff2list
from mycotools.utils.extractHmmAcc import grabAccs
import matplotlib.pyplot as plt

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

def parseClusGff(gff_list, ome):
    clus = {}
    for entry in gff_list:
        clusID = re.search(r'clusID=([^;]+)', entry['attributes'])[1]
        if clusID not in clus:
            clus[clusID] = []
        pfamSearch = re.search(r'PFAM=([^;]+)', entry['attributes'])
        if pfamSearch is not None:
            pfams = pfamSearch[1].replace('"','').split('|')
            clus[clusID].extend(pfams)

    clus = {
        k: sorted(list(set(v))) for k, v in sorted(clus.items(), key = lambda x: x[0])
        }

    counts, pfamCounts = Counter(), Counter()
    for clusID, pfams in clus.items():
        counts += Counter(combinations(pfams, 2))
        pfamCounts += Counter(pfams)

    return ome, tuple((k, v) for k, v in counts.items()), tuple((k, v) for k, v in pfamCounts.items())

#


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
                    break
                pfams = [x.split('|') for x in data[-1].split(';')]
                gene2pfam = {k: pfams[i] for i, k in enumerate(genes)}
                for i, hg in enumerate(hgs):
                    try:
                        int(hg)
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

def parse_pfam(ome, pfam_res):
    gene2pfam = defaultdict(dict)
    pfam2def = {}
    with open(pfam_res, 'r') as raw:
        for line in raw:
            if line.startswith('#'):
                continue
            data = line.rstrip().split()
            g, p, q, b = data[0], data[3], data[2], float(data[5])
            # this is where a post-hoc threshold would be applied
            if q not in gene2pfam[g]: # in case there is a duplicate
            # for the same gene, add in the better score 
                gene2pfam[g][q] = b
            else:
                if b > gene2pfam[g][q]:
                    gene2pfam[g][q] = b
            pfam2def[p[:p.find('.')]] = q

#    for g,ps in gene2pfam.items():
 #       ps = [k for k,v in sorted(ps.items(), key = lambda x: x[1], reverse = True)][0]
        # grab the best bit for the gene

    return ome, tuple([tuple([k, tuple(list(v.keys()))]) for k, v in gene2pfam.items()])



def combineGffs(results, db, pfam2name, rank = 'species'):
    db = db.set_index()

    counts, pfamCounts, taxonLens = {}, {}, {}
    for ome, mCounts, mPfams in results:
        taxon = db[ome]['taxonomy'][rank]
        if taxon not in counts: 
            counts[taxon], taxonLens[taxon], pfamCounts[taxon] = Counter(), 0, Counter()
        taxonLens[taxon] += 1
        tCounts = {}
        for codes, val in mCounts:
            names = []
            for acc in codes:
                if acc[:acc.find('.')] in pfam2name:
                    names.append(pfam2name[acc[:acc.find('.')]])
                else:
                    names.append(acc)
            tCounts[tuple(sorted(names))] = val


        tPfams = {}
        for pfam, val in mPfams:
            abbr = pfam[:pfam.find('.')]
            if abbr in pfam2name:
                name = pfam2name[abbr]
            else:
                name = pfam
            tPfams[name] = val
    
        pfamCounts[taxon] += Counter(tPfams)
        counts[taxon] += Counter(tCounts)

    fCounts, fPfams = Counter(), Counter()
    for taxon in counts:
        fCounts += Counter({
            k: round(v/taxonLens[taxon]) for k, v in counts[taxon].items()
            })
        fPfams += Counter({
            k: round(v/taxonLens[taxon]) for k, v in pfamCounts[taxon].items()
            })

    return fCounts, fPfams

def combine_counts(counts_res, db, rank = 'species'):
    # combine counts from multiprocessing
    db = db.set_index()

    counts, hg_counts, hg2pfam, taxonLens = {}, {}, {}, {}
    for ome, tCounts, tHGcounts, tHG2Pfam in counts_res:
        taxon = db[ome]['taxonomy'][rank]
        if taxon not in counts:
            counts[taxon], hg_counts[taxon] = Counter(), Counter()
            hg2pfam[taxon], taxonLens[taxon] = {}, 0
        taxonLens[taxon] += 1
        counts[taxon] += Counter({x[0]: x[1] for x in tCounts})
        hg_counts[taxon] += Counter({x[0]: x[1] for x in tHGcounts})
        for og, pfams in tHG2Pfam:
            if og not in hg2pfam[taxon]:
                hg2pfam[taxon][og] = Counter()
            hg2pfam[taxon][og] += Counter(pfams)

    fCounts, fHGcounts, fHG2Pfam = Counter(), Counter(), {}
    for taxon in counts:
        tCounts = {k: round(v/taxonLens[taxon]) for k, v in counts[taxon].items()}
        tHGcounts = {k: round(v/taxonLens[taxon]) for k, v in hg_counts[taxon].items()}
        tHG2Pfam = {}
        for og, pfams in hg2pfam[taxon].items():
            tHG2Pfam[og] = {k: round(v/taxonLens[taxon]) for k, v in pfams.items()}

        fCounts += Counter(tCounts)
        fHGcounts += Counter(tHGcounts)
        for og, pfams in tHG2Pfam.items():
            if og not in fHG2Pfam:
                fHG2Pfam[og] = Counter()
            fHG2Pfam[og] += pfams

    for og in fHG2Pfam:
        fHG2Pfam[og] = {k: v for k,v in sorted(fHG2Pfam[og].items(), key = lambda x: x[1], reverse = True)}
            
    return fCounts, fHGcounts, fHG2Pfam

def write_node_quants(hg_counts, out_path):

    minC, maxC = min(hg_counts.values()), max(hg_counts.values())
    denom = maxC - minC
    hg_counts = [
        (k, v) for k,v in sorted(hg_counts.items(), key = lambda x: x[1],
        reverse = True)
        ]
    with open(out_path, 'w') as out:
        for k,v in hg_counts:
            out.write(str(k) + '\t' + str(((v-minC)/denom) + 0.01) + '\n')


def increaseFloat(val):
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


def writeToAdj(oEdges, edges, out_path):
    minCount, maxCount = min(oEdges.values()), max(oEdges.values())
    denom = maxCount - minCount
    if denom == 0:
        eprint('\t\t\tWARNING: low connection weight throughout network', flush = True)
        denom = 1
    with open(out_path, 'w') as outF:
        for k, v in edges.items():
            value = ((v-minCount)/denom)
            out_data = [k[0], k[1], value]
            outF.write('\t'.join([str(x) for x in out_data]) + '\n')



def write_adj_matr(oEdges, out_path, maximum = 500, duplicates = False):
    # write adjancency matrix of counts normalized scaled with factor

    if not duplicates:
        oEdges = {k: v for k,v in oEdges.items() if k[0] != k[1]}

    data = sorted(list(oEdges.values()))
    edges = oEdges
    writeToAdj(oEdges, edges, out_path + '.raw')

    percentile = 0.0
    while len(edges) > maximum: # while too many edges
        percentile = increaseFloat(percentile)
        minPerc = data[round((len(data)-1)*percentile)]
        edges = {k: 4*v for k, v in oEdges.items() if v > minPerc}

    print('\t\tEdge #|percentile:', len(edges), percentile, flush = True)
    writeToAdj(oEdges, edges, out_path)


def write_net_file(
    adjMtx_path, quant_path, network_path, hg2pfam_path,
    edgeConnect = False, centrality = False
    ):
    # open adjacency matrix and construct a network as an html

    quants = {}
    with open(quant_path, 'r') as raw:
        for line in raw:
            k, v = line.rstrip().split('\t')
            quants[k] = round(float(v)*4750)


    if os.path.isfile(hg2pfam_path):
        hg2pfamPrep = read_json(hg2pfam_path)
        hg2pfam = {}
        for hg, pfams in hg2pfamPrep.items():
            if pfams:
                v = list(pfams.keys())[0]
                hg2pfam[hg] = hg + '|' + v[v.find('-')+1:]
            else:
                hg2pfam[hg] = hg 
    else:
        hg2pfam = {x: x for x in quants}

    if len(quants) < 2:
        eprint('\t\tERROR: less than 2 nodes in network', flush = True)
        return



    if os.path.isfile(adjMtx_path + '.raw'):
        print('\t\tQuantifying graph summary', flush = True)
        G = nx.read_weighted_edgelist(adjMtx_path + '.raw')

        if edgeConnect:
            print('\t\t\tQuantifying normalized graph connectivity', flush = True)
            print('\t\t\t\tNode:', flush = True)
            nodeConnectivity = nx.average_node_connectivity(G)
            print('\t\t\t\t\t' + str(nodeConnectivity), flush = True)
            print('\t\t\t\tEdge:', flush = True)
            edgeConnectivity = nx.edge_connectivity(G)
            print('\t\t\t\t\t' + str(edgeConnectivity), flush = True)
        elif centrality:
            print('\t\t\t\tQuantifying betweenness centrality', flush = True)
            cent_dict = nx.betweenness_centrality(G)
            central_path = adjMtx_path[:-4] + '.centrality'
            with open(central_path, 'w') as out:
                out.write('\n'.join([
                    str(hg2pfam[k]) + '\t' + str(v) + '\t' + str(v/(quants[k]/4750))\
                    for k,v in sorted(
                        cent_dict.items(), key = lambda x: x[1], 
                        reverse = True
                        ) if v > 0 and k in quants # why isnt k always in quants
                    ]))
    net = nx.read_weighted_edgelist(adjMtx_path)

#    bb = nx.betweenness_centrality(net)
#    nx.set_node_attributes(net, bb, "betweenness")
    nx.set_node_attributes(net, quants, "size") # verify node sizes are applied right
    pos = nx.nx_pydot.pydot_layout(net)
#    pos = nx.get_node_attributes(net, 'pos')
    nodeSize = nx.get_node_attributes(net, 'size')
    nodeQuant = len(nodeSize)
    edgeLabels = nx.get_edge_attributes(net, 'weight')
    edgeQuant = len(edgeLabels)
    edgeLabels = {k: 4*v for k, v in edgeLabels.items()}
    
    print('\t\tNodes:', nodeQuant, flush = True)
    print('\t\tEdges:', len(edgeLabels), flush = True)

    print('\t\tCommunity detection', flush = True)
    # modularity detection via multilevel algorithm
#    h = Graph.from_networkx(net)
    h = Graph.TupleList(net.edges(), directed = False)
    com = Graph.community_multilevel(h, list(edgeLabels.values())) # verify edgelabels are extracted right
    membership = com.membership
    modulesPrep = {k: v for k, v in zip(h.vs["name"], membership)}
    modules = {}
    for og, clus in modulesPrep.items():
        if clus not in modules:
            modules[clus] = []
        modules[clus].append(og)
    modules = {k: v for k,v in sorted(modules.items(), key = lambda x: len(x[1]), reverse = True)}

    print('\t\t\tModules w/> 1 HG:', len([v for v in modules.values() if len(v) > 1]), flush = True)

    dimScale = 1+round(edgeQuant**0.65)
    if dimScale > 650:
        print('\t\tWARNING: extremely large network may not format well', flush = True)
        deScale = 650/dimScale
        dimScale = 650
        quants = {k: v*deScale for k,v in quants.items()}
        nx.set_node_attributes(net, quants, "size")
        edgeLabels = {k: v*deScale for k,v in edgeLabels.items()}


    fig = plt.figure(figsize=(dimScale,dimScale))
    ax = plt.subplot(111)

# just draw it outright to get all edges in there
    nx.draw(
        net, pos, node_size=list(nodeSize.values()), node_color='white' #, font_size=8,
#        with_labels = True, font_weight='bold'
        )

    colors = getColors(len([k for k, v in modules.items() if len(v) > 1]), ignore = ['#ffffff'])
    count, edgeCount, edgeColor = 0, -1, 'black'
    for module, nodes in modules.items():
        if len(nodes) > 1:
            try:
                color = colors[count]
                count += 1
                if edgeCount > -1:
                    edgeColor = colors[edgeCount]
                else:
                    edgeColor = 'black'
            except IndexError: # exceeds color
                count = 0
                color = colors[count]
                edgeCount += 1
                edgeColor = colors[edgeCount]
        else:
            color = '#ffffff'
            edgeColor = 'black'
        nx.draw_networkx_nodes(
            nx, pos, nodelist = tuple(nodes), node_color = color,
            node_size = tuple([quants[k] for k in nodes])
            #, font_size = 8, font_weight = 'bold', with_labels = True
            ).set_edgecolor(edgeColor)
        rgb = hex2rgb(color)
        if (rgb[0]*0.299 + rgb[1]*0.587 + rgb[2]*0.114) > 186: 
           fontColor = '#000000' 
        else:
           fontColor = '#aaaaaa'
        nx.draw_networkx_labels(
            nx, pos, labels = {k: hg2pfam[k] for k in nodes}, font_size = 8, font_weight = 'bold', font_color = fontColor
            ) 

# edge labels
#    nx.draw_networkx_edge_labels(
 #       net, pos, edge_labels = edgeLabels, 
  #      font_size = 6
   #     )


    nx.draw_networkx_edges(
        net, pos, edgelist = edgeLabels.keys(), width = list(edgeLabels.values()),
        alpha = 0.4
        )
    try:
        plt.tight_layout()
    except ValueError: # too big
        ax.update_datalim(((0,0), (65535, 65535)))

    plt.show(block = True)
    plt.savefig(network_path, format="svg")
    plt.close()

# pyvis
#    nt = Network('500px', '500px')
 #   nt.from_nx(net)
  #  nt.show(network_path)

def main(
    db, rank, cloci_dir, connectivity = False, 
    centrality = False, max_edges = 500, hlg = 'hlg',
    cpus = 1
    ):

    if cloci_dir:
        ome_dir, net_dir = cloci_dir + 'ome/', cloci_dir + 'net/'
        if not os.path.isdir(ome_dir):
            eprint('\nERROR: invalid cloci output', flush = True)
            sys.exit(2)

    if not os.path.isdir(net_dir):
        os.mkdir(net_dir)
    omesXlineage = extract_lineages(db, rank)

    omesXlineage = {
        k: v for k,v in sorted(omesXlineage.items(), key = lambda x: len(x[1]))
        }
    for lineage, omes in omesXlineage.items():
        print('\n' + lineage, flush = True)
        adj_path = net_dir + lineage + '.adj'
        quant_path = net_dir + lineage + '.tsv'
        hg2pfam_path = net_dir + lineage + '.hg2pfam.txt'
        if not os.path.isfile(adj_path): # need to log the factor then
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
            counts, hg_counts, hg2pfam = combine_counts(
                counts_res, db, rank = 'species'
                )
            write_json(hg2pfam, hg2pfam_path)
            write_node_quants(hg_counts, quant_path)
            write_adj_matr(counts, adj_path, max_edges)
        
        net_path = adj_path[:-4] + '.svg'
        if not os.path.isfile(net_path):
            print('\tNetwork', flush = True)
            write_net_file(adj_path, quant_path, net_path, 
                           hg2pfam_path, connectivity, centrality)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Takes CLOCI output, creates networks of HG' + \
        ' relationships'
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
    parser.add_argument(
        '-m', '--maximum', help = 'Find upper percentile of edges to fit ' + \
        'maximum edges; DEFAULT: 500',
        default = 500, type = int
        )
    parser.add_argument(
        '--centrality', help = 'Calculate betweenness centrality',
        action = 'store_true'
        )
    parser.add_argument(
        '--connectivity', help = 'Quantify node and edge connectivity',
        action = 'store_true'
        )
    parser.add_argument(
        '-c', '--cpu', default = 1, type = int
        )
    args = parser.parse_args()

    db = mtdb(format_path(args.database))
    cloci_dir = format_path(args.input)

    if args.gcfs:
        hlg = 'gcf'
    else:
        hlg = 'hlg'

    main(
        db, args.rank, cloci_dir, connectivity = args.connectivity, 
        centrality = args.centrality, hlg = hlg,
        max_edges = args.maximum, cpus = args.cpu
        )
    sys.exit(0)
