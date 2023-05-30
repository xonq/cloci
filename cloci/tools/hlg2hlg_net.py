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
from mycotools.lib.dbtools import mtdb, primaryDB

# NEED to deal with multiple familes
# NEED hlg2name conversion
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

def import_hlgs(hlg_in, out_dir):
    hlg_set, clans = set(hlg_in), set()
    with gzip.open(f'{out_dir}hlgs.tsv.gz', 'rt') as raw:
        for line in raw:
            if not line.startswith('#'):
                hgs, hlg, clan, tmd, dp, gcc, mmi, mmp, omes = \
                    line.rstrip().split()
                hlg = int(hlg)
                clan = int(clan)
                if hlg in hlg_set:
                    clans.add(clan)
    return clans
       
def import_loci(loci_file, clans, hlgs):
    locI2hlg = {}
    with open(loci_file, 'r') as raw:
        if hlgs and clans:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, hlg, loc = line.rstrip().split()
                    clan, hlg = int(clan), int(hlg)
                    if clan in clans and hlg in hlgs:
                        locI = int(locI)
                        hlg = int(hlg)
                        ome = loc[:loc.find('_')]
                        locI2hlg[locI] = (hlg, ome)
        elif clans:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, hlg, loc = line.rstrip().split()
                    clan = int(clan)
                    if clan in clans:
                        locI = int(locI)
                        hlg = int(hlg)
                        ome = loc[:loc.find('_')]
                        locI2hlg[locI] = (hlg, ome)
        elif hlgs:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, hlg, loc = line.rstrip().split()
                    hlg = int(hlg)
                    if hlg in hlgs:
                        locI = int(locI)
                        ome = loc[:loc.find('_')]
                        locI2hlg[locI] = (hlg, ome)
        else:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, hlg, loc = line.rstrip().split()
                    clan = int(clan)
                    locI = int(locI)
                    hlg = int(hlg)
                    ome = loc[:loc.find('_')]
                    locI2hlg[locI] = (hlg, ome)
    return {k: v for k, v in sorted(locI2hlg.items(), key = lambda x: x[0])}

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

def make_network(net_file, adj_arr, locIdata, nl2ol, img, 
                 annotate = False, locI2color = {}, locI2name = {}, scale = 1):
    print('\tLoading adjacency matrix', flush = True)
    g = Graph(directed = False)
    idx = adj_arr.nonzero()
    weights = np.log(adj_arr[idx].toarray())
    v = weights[:, 1]
    weights[:, 1] = (v - v.min()) / (v.max() - v.min()) * scale

    g.add_edge_list(np.transpose(adj_arr.nonzero()))

    # old_hlg2new_locI
    modules = defaultdict(list)
    for locI, data in locIdata.items():
        modules[data[0]].append(locI)

    # old hlg 2 new hlg
    op2np = {v: i for i, v in enumerate(sorted(modules.keys()))}
    np2op = {v: k for k, v in op2np.items()}

    # new_hlg2new_locI
    modules = {op2np[k]: v for k, v in sorted(modules.items(), key = lambda x: x[0])}
    # new_locI2new_hlg
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
        rgb = [x/255.0 for x in hex2rgb(color)]
        if (rgb[0]*0.299 + rgb[1]*0.587 + rgb[2]*0.114) > 186:
           fontColor = '#000000'
        else:
           fontColor = '#aaaaaa'
        modcol.append(rgb)
        hexcol.append(color)
    vname = g.new_vertex_property('string')
    vprop = g.new_vertex_property('vector<float>')
    vprop2 = g.new_vertex_property('string')
    v_fil_col = g.new_vertex_property('vector<float>')

 #   for locI, data in locIdata.items():
#        v = g.add_vertex()
    if not annotate:
        for v in g.vertices():
            nlocI = int(v)
            locI = nl2ol[nlocI]
            nhlg = nl2np[nlocI]
   #         vprop[v] = modcol[nhlg]
            if modules[nhlg][-1] == nlocI and not locI2name:
    #            ol = nl2ol[int(v)]
    #            name= locIdata[v][1] + '_' + str(ol)
                vname[v] = str(np2op[nhlg])
            elif locI in locI2name:
                vname[v] = locI2name[locI]
            if locI in locI2color:
                v_fil_col[v] = locI2color[locI]
            else:
                v_fil_col[v] = [1.0, 1.0, 1.0, 0.3]  
        eprop = g.new_edge_property('float')
 #   eop = g.new_edge_property('float')
        eprop.a = weights #* scale
  #  eop.a = opacity
        g.ep['edge_weight'] = eprop
    else:
        eprop = g.new_edge_property('float')
 #   eop = g.new_edge_property('float')
        eprop.a = weights #* scale
  #  eop.a = opacity
        g.ep['edge_weight'] = eprop

        for v in g.vertices():
            nlocI = int(v)
            locI = nl2ol[nlocI]
            nhlg = nl2np[nlocI]
            vprop[v] = modcol[nhlg]
            vprop2[v] = hexcol[nhlg]
            if locI in locI2color:
                v_fil_col[v] = [x for x in locI2color[locI]]
  #              vprop[v] = locI2color[locI]
            else:
                v_fil_col[v] = (1,1,1,0.3)  
 #               vprop[v] = (255,255,255,1)  

            ol = nl2ol[int(v)]
#            name = locIdata[int(v)][1] + '_' + str(ol)
            if locI in locI2name:
                vname[v] = locI2name[locI]
            else:
                vname[v] = str(ol)
#    g.vp['color'] = vprop
    g.vp['name'] = vname
    g.vp['hex'] = vprop2
    g.vp['fill'] = v_fil_col


    if img:
        graph_draw(g, vertex_halo_color = g.vp['fill'], vertex_color = [0,0,0,1],
                   vertex_fill_color = g.vp['fill'], 
 #                   vertex_text_position = 1, vertex_aspect = 1,
    #                vertex_text = g.vp['name'], #vertex_size = g.degree_property_map('total'),
     #                  vertex_font_size=12, vertex_text = g.vp['name'], 
#                               vertex_font_size=18, vertex_size = 5, output = (dim, dim),
                   edge_color = [0, 0, 0, 1], output = net_file, edge_pen_width = eprop)
#         graphviz_draw(g, vcolor = g.vp['hex'], ecolor = '#000000', overlap = False,
  #                     output = net_file, penwidth = eprop, vprops = {'vertex_text': g.vp['name']}) 
    else:
        g.save(net_file)



def parse_hlgs(hlg_file):
    clans, hlgs = set(), set()
    with gzip.open(hlg_file, 'rt') as raw:
        for line in raw:
            if not line.startswith('#'):
                hgs, hlg, clan, nlt, p, gcc, mmi, mmp, omes = line.rstrip().split()
                hlg, clan = int(hlg), int(clan)
                clans.add(clan)
                hlgs.add(hlg)
    return hlgs, clans


def main(clans, hlgs, res_dir, minimum, passing = False, db = None, rank = None,
         ext = 'pdf', img = True, annotate = False, highlight = set(), scale = 0.25):
    net_dir = res_dir + 'net/'
    if not os.path.isdir(net_dir):
        os.mkdir(net_dir)

    run = 0
    while os.path.isfile(f'{net_dir}{run}.{ext}'):
        run += 1
    out_file = f'{net_dir}{run}.{ext}'

    if passing:
        print('\nImporting passing GCFs and clans', flush = True)
        pass_hlgs, pass_clans = parse_hlgs(f'{res_dir}hlgs.tsv.gz')
    else:
        pass_clans, pass_hlgs = set(), set()
    if pass_clans and clans:
        diff = set(clans).difference(pass_clans)
        if set(clans).difference(pass_clans):
            eprint(f'\tWARNING: nonpassing clans: {",".join([str(x) for x in diff])}',
                   flush = True)
            clans = list(pass_clans.intersection(set(clans)))
    elif pass_clans:
        clans = list(pass_clans)
    else:
        pass_clans, pass_hlgs = list(clans), set(hlgs)

    print('\nImporting loci', flush = True)
    if clans:
        locI2hlg = import_loci(f'{res_dir}working/hlg/loci.txt', set(clans), pass_hlgs)
    else:
        locI2hlg = import_loci(f'{res_dir}working/hlg/loci.txt', set(), pass_hlgs)
    ol2nl = {locI: i for i, locI in enumerate(list(locI2hlg.keys()))}
    nl2ol = {v: k for k, v in ol2nl.items()}

    if not locI2hlg:
        eprint(f'\tERROR: clans {clans} were not recovered', flush = True)
        sys.exit(5)

    print('\nBuilding adjacency matrix', flush = True)
    adj_arr = import_adj(f'{res_dir}working/hlg/loci.adj', 
                         set(locI2hlg.keys()), minimum, ol2nl)

    print('\nMaking network', flush = True)
    modules = {i: v for i, v in enumerate(list(locI2hlg.values()))}
    if highlight:
        locI2color = {x: (250/255, 128/255, 114/255, 0.8) for x in list(highlight)}
        locI2name = {x: locI2hlg[x][1] for x in list(highlight)}
    else:
        locI2color = {}
        locI2name = {}
    make_network(out_file, adj_arr, modules, nl2ol, img, 
                 annotate, locI2color, locI2name, scale)
    print(f'\nNetwork outputted to {out_file}', flush = True)


if __name__ == '__main__':
    ranks = ['kingdom', 'phylum', 'subphylum', 'class', 'order', 'family', 'genus']
    img_exts = ['pdf', 'png', 'ps', 'svg']
    net_exts = ['dot', 'gml', 'graphml', 'gt', 'xml']
    parser = argparse.ArgumentParser(
        description = 'Cloci locus-locus similarity network. Will output all clans if no -c/-f'
        )
    parser.add_argument('-i', '--input', help = 'Completed Cloci directory', required = True)
    parser.add_argument('-c', '--clans', help = 'Input clans; "-" for stdin')
    parser.add_argument('-f', '--hlgs', help = 'Input GCFs; "-" for stdin')
    parser.add_argument('-l', '--locids', help = 'LocIDs for highlighting')
    parser.add_argument('-p', '--passing', help = 'Passing GCFs only',
                        action = 'store_true')
    parser.add_argument('-m', '--min', type = float, default = 0,
        help = '0 < -m < 1: Minimum GCF similarity for network edges')
    parser.add_argument('-a', '--annotate', action = 'store_true',
        help = 'Annotate nodes with their locusID')
    parser.add_argument('-d', '--database', help = 'MTDB for taxonomy',
        default = primaryDB())
    parser.add_argument('-r', '--rank', 
        help = f'Rank for taxonomy node coloring: {ranks}')
    parser.add_argument('-e', '--extension', help = f'Network file format: {img_exts}, {net_exts}',
        default = 'pdf')
    parser.add_argument('-s', '--scale', type = float, default = 0.1,
        help = 'Edge scaling coefficient; DEFAULT: 0.1')
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

    hlg_str, clan_str = '', ''
    if args.clans and args.hlgs:
        eprint('\nERROR: -c or -f only', flush = True)
        sys.exit(2)
    elif args.clans == '-':
        clan_str = stdin2str(args.clans)
    elif args.hlgs == '-':
        hlg_str = stdin2str(args.hlgs)
    elif args.clans:
        clan_str = args.clans
    elif args.hlgs:
        hlg_str = args.hlgs

    in_dir = format_path(args.input)

    if hlg_str:
        hlg_p = hlg_str.replace('"','').replace("'",'')
        if ',' in hlg_str:
            hlgs = [int(x) for x in hlg_p.split(',')]
        else:
            hlgs = [int(x) for x in hlg_p.split()]
        if args.passing:
            clans = import_hlgs(hlgs, in_dir)
        else:
            clans = []
    elif clan_str:
        clan_p = clan_str.replace('"','').replace("'",'')
        if ',' in clan_str:
            clans = [int(x) for x in clan_p.split(',')]
        else:
            clans = [int(x) for x in clan_p.split()]
    else:
        clans, hlgs = None, None

    if args.locids:
        locids = set(int(x) for x in \
                     args.locids.replace('"','').replace("'",'').replace(',',' ').split())
    else:
        locids = set()

    db = mtdb(format_path(args.database))
    
    main(clans, hlgs, format_path(in_dir, force_dir = True), 
         args.min, args.passing, db, rank, ext = ext, img = img,
         annotate = args.annotate, highlight = locids, scale = args.scale)
