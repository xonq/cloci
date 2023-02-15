import os
import sys
import gzip
import shutil
import pickle
import argparse
#import networkx as nx
import numpy as np
from graph_tool.all import *
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.sparse import lil_matrix
from mycotools.lib.kontools import eprint, format_path, stdin2str, getColors, hex2rgb
from mycotools.lib.dbtools import mtdb, masterDB

# NEED to deal with multiple familes
# NEED gcf2name conversion
# NEED node annotation option
     # has to space nodes effectively, probably graphvis
# NEED taxonomy coloring
# NEED to filter irrelevant omes
# NEED option to adjust vertex size by degree
# NEED option to filter for only one representative of same species per graph
    # ideally need a better way to address sampling bias, must be phylogeny-aware
# NEED centrality, connectedness calculations
# NEED to choose size of everything intelligently
    # likely nonlinearly proportional to edge number

def import_gcfs(gcf_in, out_dir):
    gcf_set, clans = set(gcf_in), set()
    with gzip.open(f'{out_dir}gcfs.tsv.gz', 'rt') as raw:
        for line in raw:
            if not line.startswith('#'):
                hgs, gcf, clan, tmd, dp, gcc, mmi, mmp, omes = \
                    line.rstrip().split()
                gcf = int(gcf)
                clan = int(clan)
                if gcf in gcf_set:
                    clans.add(clan)
    return clans
       
def import_loci(loci_file, clans, gcfs):
    locI2gcf = {}
    with open(loci_file, 'r') as raw:
        if gcfs and clans:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, gcf, loc = line.rstrip().split()
                    clan, gcf = int(clan), int(gcf)
                    if clan in clans and gcf in gcfs:
                        locI = int(locI)
                        gcf = int(gcf)
                        ome = loc[:loc.find('_')]
                        locI2gcf[locI] = (gcf, ome)
        elif clans:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, gcf, loc = line.rstrip().split()
                    clan = int(clan)
                    if clan in clans:
                        locI = int(locI)
                        gcf = int(gcf)
                        ome = loc[:loc.find('_')]
                        locI2gcf[locI] = (gcf, ome)
        elif gcfs:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, gcf, loc = line.rstrip().split()
                    gcf = int(gcf)
                    if gcf in gcfs:
                        locI = int(locI)
                        ome = loc[:loc.find('_')]
                        locI2gcf[locI] = (gcf, ome)
        else:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, gcf, loc = line.rstrip().split()
                    clan = int(clan)
                    locI = int(locI)
                    gcf = int(gcf)
                    ome = loc[:loc.find('_')]
                    locI2gcf[locI] = (gcf, ome)
    return {k: v for k, v in sorted(locI2gcf.items(), key = lambda x: x[0])}

def import_adj(adj_file, loci2grab, minimum, ol2nl):
#    adj_arr = np.zeros([len(loci2grab), len(loci2grab)])
    loclen = len(loci2grab)
    adj_arr = lil_matrix((loclen, loclen))

    count = 0
    with open(adj_file, 'r') as raw:
        for line in raw:
            l0, l1, gid = line.rstrip().split()
            l0, l1 = int(l0), int(l1)
            if l0 in loci2grab and l1 in loci2grab:
#                val = 50**float(gid)
                val = float(gid)
                if val >= minimum:
                    adj_arr[(ol2nl[l0], ol2nl[l1])] = val

    return adj_arr

def make_network(net_file, adj_arr, locIdata, nl2ol, img, annotate = False):
    print('\tLoading adjacency matrix', flush = True)
#    G = nx.read_weighted_edgelist(net_file)
 #   pos = nx.nx_pydot.pydot_layout(G)
    g = Graph(directed = False)
    idx = adj_arr.nonzero()
    weights = adj_arr[idx].toarray()
    g.add_edge_list(np.transpose(adj_arr.nonzero()))

 #   dim = len(idx[0]) 
#    weights *= 10

    # normalize 
#    max_weight = np.max(weights)
 #   opacity = weights/max_weight

  #  n_opacity = np.select([True], [np.array([0, 0, 0, opacity])], opacity)
   # print(n_opacity)


    # old_gcf2new_locI
    modules = defaultdict(list)
    for locI, data in locIdata.items():
        modules[data[0]].append(locI)

    # old gcf 2 new gcf
    op2np = {v: i for i, v in enumerate(sorted(modules.keys()))}
    np2op = {v: k for k, v in op2np.items()}

    # new_gcf2new_locI
    modules = {op2np[k]: v for k, v in sorted(modules.items(), key = lambda x: x[0])}
    # new_locI2new_gcf
    nl2np = {}
    for p, nls in modules.items():
        for nl in nls:
            nl2np[nl] = p

    colors = getColors(len(modules), ignore = ['#ffffff'])

    print(f'\tDrawing {len(locIdata)} nodes', flush = True)
    count, edgeCount, edgeColor, modcol, hexcol = 0, -1, 'black', [], []
    for module, nodes in modules.items():
        nodes = [str(x) for x in nodes]
        if len(nodes) > 1:
            try:
                color = colors[count]
                count += 1
#                if edgeCount > -1:
 #                   edgeColor = colors[edgeCount]
  #              else:
   #                 edgeColor = 'black'
            except IndexError: # exceeds color
                count = 0
                color = colors[count]
    #            edgeCount += 1
     #           edgeColor = colors[edgeCount]
        else:
            color = '#ffffff'
            edgeColor = 'black'
#        nx.draw_networkx_nodes(
 #           G, pos, nodelist = tuple(nodes), node_color = color,
            #, font_size = 8, font_weight = 'bold', with_labels = True
  #          ).set_edgecolor(edgeColor)
        rgb = [x/255.0 for x in hex2rgb(color)]
        if (rgb[0]*0.299 + rgb[1]*0.587 + rgb[2]*0.114) > 186:
           fontColor = '#000000'
        else:
           fontColor = '#aaaaaa'
        modcol.append(rgb)
        hexcol.append(color)
   #     nx.draw_networkx_labels(
    #        G, pos, labels = {k: f'{gcfs[int(k)][1]}_{nl2ol[int(k)]}' for k in nodes}, 
     #       font_size = 8, font_weight = 'bold', font_color = fontColor
      #      )
    vname = g.new_vertex_property('string')
    vprop = g.new_vertex_property('vector<float>')
    vprop2 = g.new_vertex_property('string')

 #   for locI, data in locIdata.items():
#        v = g.add_vertex()
    if not annotate:
        for v in g.vertices():
            ngcf = nl2np[int(v)]
            vprop[v] = modcol[ngcf]
            if modules[ngcf][-1] == int(v):
    #            ol = nl2ol[int(v)]
    #            name= locIdata[v][1] + '_' + str(ol)
                vname[v] = str(np2op[ngcf])
        eprop = g.new_edge_property('float')
 #   eop = g.new_edge_property('float')
        eprop.a = weights * 2.5
  #  eop.a = opacity
        g.ep['edge_weight'] = eprop
    else:
        eprop = g.new_edge_property('float')
 #   eop = g.new_edge_property('float')
        eprop.a = weights * 2.5
  #  eop.a = opacity
        g.ep['edge_weight'] = eprop

        for v in g.vertices():
            ngcf = nl2np[int(v)]
            vprop[v] = modcol[ngcf]
            vprop2[v] = hexcol[ngcf]
            ol = nl2ol[int(v)]
#            name = locIdata[int(v)][1] + '_' + str(ol)
            vname[v] = str(ol)
    g.vp['color'] = vprop
    g.vp['name'] = vname
    g.vp['hex'] = vprop2


  #  g.vertex_properties['name'] = vname
#    vprops_dict = {'fill_color': vprop, 'font_family': 'DejaVu Sans',
 #                  'text': vname, 'font_size': 8}
    if img:
         graph_draw(g, vertex_fill_color = [1,1,1,1], vertex_color = g.vp['color'],
 #                   vertex_text_position = 1, vertex_aspect = 1,
                    vertex_text = g.vp['name'], #vertex_size = g.degree_property_map('total'),
     #                  vertex_font_size=12, vertex_text = g.vp['name'], 
#                               vertex_font_size=18, vertex_size = 5, output = (dim, dim),
                   edge_color = [0, 0, 0, 1], output = net_file, edge_pen_width = eprop)
#         graphviz_draw(g, vcolor = g.vp['hex'], ecolor = '#000000', overlap = False,
  #                     output = net_file, penwidth = eprop, vprops = {'vertex_text': g.vp['name']}) 
    else:
        g.save(net_file)



def parse_gcfs(gcf_file):
    clans, gcfs = set(), set()
    with gzip.open(gcf_file, 'rt') as raw:
        for line in raw:
            if not line.startswith('#'):
                hgs, gcf, clan, nlt, p, gcc, mmi, mmp, omes = line.rstrip().split()
                gcf, clan = int(gcf), int(clan)
                clans.add(clan)
                gcfs.add(gcf)
    return gcfs, clans


def main(clans, gcfs, res_dir, minimum, passing = False, db = None, rank = None,
         ext = 'pdf', img = True, annotate = False):
    net_dir = res_dir + 'net/'
    if not os.path.isdir(net_dir):
        os.mkdir(net_dir)

    run = 0
    while os.path.isfile(f'{net_dir}{run}.{ext}'):
        run += 1
    out_file = f'{net_dir}{run}.{ext}'

    if passing:
        print('\nImporting passing GCFs and clans', flush = True)
        pass_gcfs, pass_clans = parse_gcfs(f'{res_dir}gcfs.tsv.gz')
    else:
        pass_clans, pass_gcfs = set(), set()
    if pass_clans and clans:
        diff = set(clans).difference(pass_clans)
        if set(clans).difference(pass_clans):
            eprint(f'\tWARNING: nonpassing clans: {",".join([str(x) for x in diff])}',
                   flush = True)
            clans = list(pass_clans.intersection(set(clans)))
    elif pass_clans:
        clans = list(pass_clans)
    else:
        pass_clans, pass_gcfs = list(clans), set(gcfs)

    print('\nImporting loci', flush = True)
    if clans:
        locI2gcf = import_loci(f'{res_dir}working/gcf/loci.txt', set(clans), pass_gcfs)
    else:
        locI2gcf = import_loci(f'{res_dir}working/gcf/loci.txt', set(), pass_gcfs)
    ol2nl = {locI: i for i, locI in enumerate(list(locI2gcf.keys()))}
    nl2ol = {v: k for k, v in ol2nl.items()}

    if not locI2gcf:
        eprint(f'\tERROR: clans {clans} were not recovered', flush = True)
        sys.exit(5)

    print('\nBuilding adjacency matrix', flush = True)
    adj_arr = import_adj(f'{res_dir}working/gcf/loci.adj', 
                         set(locI2gcf.keys()), minimum, ol2nl)

    print('\nMaking network', flush = True)
    modules = {i: v for i, v in enumerate(list(locI2gcf.values()))}
    make_network(out_file, adj_arr, modules, nl2ol, img, annotate)
    print(f'\nNetwork outputted to {out_file}', flush = True)


if __name__ == '__main__':
    ranks = ['kingdom', 'phylum', 'subphylum', 'class', 'order', 'family', 'genus']
    img_exts = ['pdf', 'png', 'ps', 'svg']
    net_exts = ['dot', 'gml', 'graphml', 'gt', 'xml']
    parser = argparse.ArgumentParser(
        description = 'hg2clus locus-locus similarity network. Will output all clans if no -c/-f'
        )
    parser.add_argument('-i', '--input', help = 'Completed hg2clus directory', required = True)
    parser.add_argument('-c', '--clans', help = 'Input clans; "-" for stdin')
    parser.add_argument('-f', '--gcfs', help = 'Input GCFs; "-" for stdin')
    parser.add_argument('-p', '--passing', help = 'Passing GCFs only',
                        action = 'store_true')
    parser.add_argument('-m', '--min', type = float, default = 0,
        help = '0 < -m < 1: Minimum GCF similarity for network edges')
    parser.add_argument('-a', '--annotate', action = 'store_true',
        help = 'Annotate nodes with their locusID')
    parser.add_argument('-d', '--database', help = 'MTDB for taxonomy',
        default = masterDB())
    parser.add_argument('-r', '--rank', 
        help = f'Rank for taxonomy node coloring: {ranks}')
    parser.add_argument('-e', '--extension', help = f'Network file format: {img_exts}, {net_exts}',
        default = 'pdf')
    args = parser.parse_args()

    if args.rank:
        rank = args.rank.lower()
        if rank not in set(ranks):
            eprint(f'\nERROR: unacceptable rank "{rank}"', flush = True)
            sys.exit(3)
    else:
        rank = None

    if not args.extension.lower() in set(img_exts).union(set(net_exts)):
        eprint(f'\nERROR: invalid file extension "{args.extension}"', flush = True)
    else:
        ext = args.extension.lower()
        if ext in set(img_exts):
            img = True
            if ext == 'svg':
                eprint('\nWARNING: svg for large networks can fail due to Cairo memory management', flush = True)
        else:
            img = False

    gcf_str, clan_str = '', ''
    if args.clans and args.gcfs:
        eprint('\nERROR: -c or -f only', flush = True)
        sys.exit(2)
    elif args.clans == '-':
        clan_str = stdin2str(args.clans)
    elif args.gcfs == '-':
        gcf_str = stdin2str(args.gcfs)
    elif args.clans:
        clan_str = args.clans
    elif args.gcfs:
        gcf_str = args.gcfs

    in_dir = format_path(args.input)

    if gcf_str:
        gcf_p = gcf_str.replace('"','').replace("'",'')
        if ',' in gcf_str:
            gcfs = [int(x) for x in gcf_p.split(',')]
        else:
            gcfs = [int(x) for x in gcf_p.split()]
        if args.passing:
            clans = import_gcfs(gcfs, in_dir)
        else:
            clans = []
    elif clan_str:
        clan_p = clan_str.replace('"','').replace("'",'')
        if ',' in clan_str:
            clans = [int(x) for x in clan_p.split(',')]
        else:
            clans = [int(x) for x in clan_p.split()]
    else:
        clans, gcfs = None, None

    db = mtdb(format_path(args.database))
    
    main(clans, gcfs, in_dir, args.min, args.passing, db, rank, ext = ext, img = img,
         annotate = args.annotate)
