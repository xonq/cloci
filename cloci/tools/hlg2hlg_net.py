#! /usr/bin/env python3

import os
import re
import sys
import gzip
import shutil
import pickle
import random
import argparse
#import networkx as nx
import numpy as np
from graph_tool import clustering, centrality
from graph_tool.all import *
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.sparse import lil_matrix
from mycotools.lib.kontools import eprint, format_path, stdin2str, getColors, \
    hex2rgb, split_input, collect_files
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.extract_mtdb import main as extract_mtdb

# NEED to deal with multiple familes
# NEED hlg2name conversion
# NEED node annotation option
     # has to space nodes effectively, probably graphvis
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


def import_loci_cat(loci_file, clans, hlgs, ome2gene2cat, ex_omes):
    locI2hlg, locI2cat = {}, {}
    with open(loci_file, 'r') as raw:
        if hlgs and clans:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, hlg, loc = line.rstrip().split()
                    a_g = loc.split(',')[0]
                    clan, hlg = int(clan), int(hlg)
                    if clan in clans and hlg in hlgs:
                        locI = int(locI)
                        hlg = int(hlg)
                        ome = loc[:loc.find('_')]
                        if ex_omes:
                            if ome not in ex_omes:
                                continue
                        locI2hlg[locI] = (hlg, ome)
                        if ome in ome2gene2cat:
                            if a_g in ome2gene2cat[ome]:
                                locI2cat[locI] = ome2gene2cat[ome][a_g]
        elif clans:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, hlg, loc = line.rstrip().split()
                    a_g = loc.split(',')[0]
                    clan = int(clan)
                    if clan in clans:
                        locI = int(locI)
                        hlg = int(hlg)
                        ome = loc[:loc.find('_')]
                        if ex_omes:
                            if ome not in ex_omes:
                                continue
                        locI2hlg[locI] = (hlg, ome)
                        if ome in ome2gene2cat:
                            if a_g in ome2gene2cat[ome]:
                                locI2cat[locI] = ome2gene2cat[ome][a_g]

        elif hlgs:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, hlg, loc = line.rstrip().split()
                    a_g = loc.split(',')[0]
                    hlg = int(hlg)
                    if hlg in hlgs:
                        locI = int(locI)
                        ome = loc[:loc.find('_')]
                        if ex_omes:
                            if ome not in ex_omes:
                                continue
                        locI2hlg[locI] = (hlg, ome)
                        if ome in ome2gene2cat:
                            if a_g in ome2gene2cat[ome]:
                                locI2cat[locI] = ome2gene2cat[ome][a_g]

        else:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, hlg, loc = line.rstrip().split()
                    a_g = loc.split(',')[0]
                    clan = int(clan)
                    locI = int(locI)
                    hlg = int(hlg)
                    ome = loc[:loc.find('_')]
                    if ex_omes:
                        if ome not in ex_omes:
                            continue
                    locI2hlg[locI] = (hlg, ome)
                    if ome in ome2gene2cat:
                        if a_g in ome2gene2cat[ome]:
                            locI2cat[locI] = ome2gene2cat[ome][a_g]

    return {k: v for k, v in sorted(locI2hlg.items(), key = lambda x: x[0])}, \
           locI2cat


       
def import_loci(loci_file, clans, hlgs, ex_omes = None):
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
                        if ex_omes:
                           if ome not in ex_omes:
                               continue
   
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
                        if ex_omes:
                           if ome not in ex_omes:
                               continue

                        locI2hlg[locI] = (hlg, ome)
        elif hlgs:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, hlg, loc = line.rstrip().split()
                    hlg = int(hlg)
                    if hlg in hlgs:
                        locI = int(locI)
                        ome = loc[:loc.find('_')]
                        if ex_omes:
                           if ome not in ex_omes:
                               continue

                        locI2hlg[locI] = (hlg, ome)
        else:
            for line in raw:
                if not line.startswith('#'):
                    i, locI, clan, hlg, loc = line.rstrip().split()
                    clan = int(clan)
                    locI = int(locI)
                    hlg = int(hlg)
                    ome = loc[:loc.find('_')]
                    if ex_omes:
                       if ome not in ex_omes:
                           continue

                    locI2hlg[locI] = (hlg, ome)
    return {k: v for k, v in sorted(locI2hlg.items(), key = lambda x: x[0])}

def import_adj_hlg(adj_file, minimum, locI2hlg):
    # collapse all loci in each hlg to a binary network
    hlg2grab = set([x[0] for x in locI2hlg.values()])
    hlgs = len(hlg2grab)

    oh2nh = {v: i for i, v in enumerate(sorted(hlg2grab))}
    nh2oh = {v: k for k, v in oh2nh.items()}

    adj_arr = lil_matrix((hlgs, hlgs))

    count = 0
    with open(adj_file, 'r') as raw:
        for line in raw:
            l0, l1, gid = line.rstrip().split()
            l0, l1 = int(l0), int(l1)
            if l0 in locI2hlg and l1 in locI2hlg:
                h0, h1 = locI2hlg[l0][0], locI2hlg[l1][0]
#                val = 50**float(gid)
                if h0 != h1:
#                val = float(gid)
 #               if val >= minimum:
                    adj_arr[(oh2nh[h0], oh2nh[h1])] = 1

    hlg2omes = defaultdict(list)
    for hlg, ome in locI2hlg.values():
        hlg2omes[hlg].append(ome)

    return adj_arr, oh2nh, nh2oh, hlg2omes


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

def make_network(net_dir, ext, adj_arr, locIdata, nl2ol, img, 
                 annotate = False, font_size = 20,
                 locI2color = {}, locI2name = {}, scale = 1,
                 prefix = None, bcentrality = False, binary = False):
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
            if not binary:
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


    if len(locIdata) < 131060:
        out_size = round(len(locIdata)/4)
    else:
        out_size = 32760
    out_size = 5000

    if not prefix:
        run = 0
        while os.path.isfile(f'{net_dir}{run}.{ext}'):
            run += 1
        net_file = f'{net_dir}{run}.{ext}'
    else:
        net_file = f'{net_dir}{prefix}.{ext}'
    

    if bcentrality:
        v_between, e_between = centrality.betweenness(g)
        central_path = net_file[:-4] + '.centrality'
        with open(central_path, 'w') as out:
            out.write('#hlg\tbetweenness\n')
            for k, v in sorted(enumerate(v_between),
                               key = lambda x: x[1], reverse = True):
                out.write(f'{nl2ol[k]}\t{v}\n')


    if not annotate:
        if img:
#
            graph_draw(g, vertex_halo_color = g.vp['fill'], vertex_color = g.vp['fill'],
                       vertex_fill_color = g.vp['fill'], output_size = [out_size, out_size],
     #                   vertex_text_position = 1, vertex_aspect = 1,
        #                vertex_text = g.vp['name'], #vertex_size = g.degree_property_map('total'),
         #                  vertex_font_size=12, vertex_text = g.vp['name'], 
    #                               vertex_font_size=18, vertex_size = 5, output = (dim, dim),
                       edge_color = [0, 0, 0, 1], output = net_file, edge_pen_width = 1) #eprop
    #        graphviz_draw(g, vcolor = g.vp['hex'], ecolor = '#000000', overlap = False,
     #                     output = net_file + '.1', penwidth = 1, size = [100, 100]) #vprops = {'vertex_text': g.vp['name']}) 
        else:
            g.save(net_file)
    else:
        if img:
            graph_draw(g, vertex_halo_color = g.vp['fill'], vertex_color = [0,0,0,1],
                       vertex_fill_color = g.vp['fill'], vertex_font_size = font_size,
                       vertex_font_family = 'arial', output_size = [out_size, out_size],
                        vertex_text = g.vp['name'], #vertex_size = g.degree_property_map('total'),
                       edge_color = [0, 0, 0, 1], output = net_file, edge_pen_width = eprop)
        else:
            g.save(net_file)

    global_transitivity = clustering.global_clustering(g, eprop)
    print(f'\tTransitivity: {global_transitivity[0]}, stdev: {global_transitivity[1]}')
    return net_file

def parse_hlgs(hlg_file):
    clans, hlgs = set(), set()
    with gzip.open(hlg_file, 'rt') as raw:
        for line in raw:
            if not line.startswith('#'):
                hgs, hlg, clan, nlt, p, gcc, mmi, mmp, csb, omes = line.rstrip().split()
                hlg, clan = int(hlg), int(clan)
                clans.add(clan)
                hlgs.add(hlg)
    return hlgs, clans


def import_annotations(ome, hlg_file):
    gene2cat = {}
    with open(hlg_file, 'r') as raw:
        for line in raw:
            if not line.startswith('#'):
                d = line.rstrip().split()
                n, h, g, G, a = d
                cat_prep = re.search(r'_([^_]+)$', n)[1]
                cats = cat_prep.split('$')
                if len(cats) > 1:
                    cat = 'hybrid'
                else:
                    cat = cats[0]
                a_g = g.split(',')[0]
                G = int(G)
                gene2cat[a_g] = cat
    return ome, gene2cat
                  


def main(clans, hlgs, res_dir, minimum, passing = False, db = None, rank = None,
         ext = 'pdf', img = True, annotate_loc = False, annotate_hlg = False,
         highlight = set(), scale = 0.25, chrono = False, font_size = 20,
         colors = [], func_dir = None, ex_omes = None, out_prefix = None,
         binary = False, bcentrality = False):
    net_dir = res_dir + 'net/'
    if not os.path.isdir(net_dir):
        os.mkdir(net_dir)

    if passing:
        prefix = 'gcf'
    else:
        prefix = 'hlg'


    if passing:
        print('\nImporting passing GCFs and clans', flush = True)
        pass_hlgs, pass_clans = parse_hlgs(f'{res_dir}{prefix}s.tsv.gz')
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

    if func_dir:
        print('\nImporting annotations', flush = True)
        hlg_files_p = collect_files(func_dir, 'tsv', recursive = True)
        hlg_fs = [x for x in hlg_files_p if os.path.basename(x).startswith(prefix)]
        ome2gene2cat = {}
        for f in hlg_fs:
            ome = os.path.basename(os.path.dirname(f))
            if ex_omes:
                if ome not in ex_omes:
                    continue
            ome, gene2cat = import_annotations(ome, f)
            ome2gene2cat[ome] = gene2cat
        print('\nImporting loci', flush = True)
        if clans:
            locI2hlg, locI2cat = import_loci_cat(f'{res_dir}working/hlg/loci.txt', set(clans), 
                                       pass_hlgs, ome2gene2cat, ex_omes)
        else:
            locI2hlg, locI2cat = import_loci_cat(f'{res_dir}working/hlg/loci.txt', set(), 
                                       pass_hlgs, ome2gene2cat, ex_omes)
    else:
        locI2cat = {}
        print('\nImporting loci', flush = True)
        if clans:
            locI2hlg = import_loci(f'{res_dir}working/hlg/loci.txt', set(clans), 
                                   pass_hlgs, ex_omes)
        else:
            locI2hlg = import_loci(f'{res_dir}working/hlg/loci.txt', set(), 
                                   pass_hlgs, ex_omes)
    ol2nl = {locI: i for i, locI in enumerate(list(locI2hlg.keys()))}
    nl2ol = {v: k for k, v in ol2nl.items()}

    if not locI2hlg:
        eprint(f'\tERROR: clans {clans} were not recovered', flush = True)
        sys.exit(5)

    print('\nBuilding adjacency matrix', flush = True)
    if not binary:
        adj_arr = import_adj(f'{res_dir}working/hlg/loci.adj', 
                             set(locI2hlg.keys()), minimum, ol2nl)
    else:
        adj_arr, ol2nl, nl2ol, hlg2omes = import_adj_hlg(f'{res_dir}working/hlg/loci.adj', 
                                                         minimum, locI2hlg)

    print('\nMaking network', flush = True)
    locI2color = {}
    locI2name = {}

    if rank:
        locI2rank = {}
        if not binary:
            for locI, d in locI2hlg.items():
                o = d[1]
                tax = db[o]['taxonomy'][rank]
                locI2rank[locI] = tax
            tax2color = {k: None for k in locI2rank.values()}
            if not colors:
                colors = getColors(len(tax2color), ignore = ['#ffffff', '#000000'])
                random.shuffle(colors)
            for i, k in enumerate(list(tax2color.keys())):
                tax2color[k] = [x/255.0 for x in hex2rgb(colors[i])]
                tax2color[k].append(0.8)
    
    
        else: # hardcode for study
            for hlg, omes in hlg2omes.items():
                locI2rank[hlg] = set()
                for ome in omes:
                    tax = db[ome]['taxonomy'][rank]
                    locI2rank[hlg].add(tax)
                tax2color = {'Mucoromycota': '#004D40', 'Basidiomycota': '#1E88E5',
                             'Ascomycota': '#D81B60', 'Hybrid': '#000000'}
                for hlg, tax in locI2rank.items():
                    if len(tax) > 1:
                        locI2rank[hlg] = 'Hybrid'
                    else:
                        locI2rank[hlg] = list(tax)[0]
            for i, k in enumerate(list(tax2color.keys())):
                tax2color[k] = [x/255.0 for x in hex2rgb(tax2color[k])]
                tax2color[k].append(0.8)
        

        locI2color = {k: tax2color[v] for k, v in locI2rank.items()}
    if locI2cat:
        cats = ['hybrid', 'NI-siderophore', 'NRPS', 'PKS', 'RiPP', 'terpene']
        # hardcode to conserve colors for study
#        cats = sorted(set(locI2cat.values()))
        if not colors:
            colors = getColors(len(cats), ignore = ['#ffffff', '#000000'])
        cat2col = {}
        for i, k in enumerate(cats):
            cat2col[k] = [x/255.0 for x in hex2rgb(colors[i])]
            cat2col[k].append(0.8)
        print(f'\tColors: {cat2col}', flush = True)
        if not binary:
            locI2color = {k: cat2col[v] for k, v in locI2cat.items()}
        else:
            hlg2cat = defaultdict(set)
            for locI, cat in locI2cat.items():
                hlg = locI2hlg[locI][0]
                hlg2cat[hlg].add(cat)
            for hlg, cat in hlg2cat.items():
                if len(cat) > 1:
                    hlg2cat[hlg] = 'hybrid'
                else:
                    hlg2cat[hlg] = list(cat)[0]
            locI2color = {k: cat2col[v] for k, v in hlg2cat.items()}

    if not binary:
        modules = {i: v for i, v in enumerate(list(locI2hlg.values()))}
    else:
        modules = {k: [k] for k in ol2nl}
        

    if highlight:
        locI2color = {x: (250/255, 128/255, 114/255, 0.8) for x in list(highlight)}
        locI2name = {x: locI2hlg[x][1] for x in list(highlight)}
    if annotate_hlg:
        locI2name = {k: str(v[0]) for k, v in locI2hlg.items()}
        annotate = True
    elif annotate_loc:
        annotate = True
    else:
        annotate = False

    if chrono:
        locI2name = {k: int(v) for k, v in locI2name.items()}
        order = sorted(set(locI2name.values()))
        i02i1 = {}
        for i, v in enumerate(order):
            i02i1[v] = i
        locI2name = {k: i02i1[v] for k, v in locI2name.items()}

    out_file = make_network(net_dir, ext, adj_arr, modules, nl2ol, img, 
                 annotate, font_size, locI2color, locI2name, scale, prefix = out_prefix,
                 bcentrality = bcentrality, binary = binary)
    print(f'\nNetwork outputted to {out_file}', flush = True)


def cli():
    ranks = ['kingdom', 'phylum', 'subphylum', 'class', 'order', 'family', 'genus']
    img_exts = ['pdf', 'png', 'ps', 'svg']
    net_exts = ['dot', 'gml', 'graphml', 'gt', 'xml']
    parser = argparse.ArgumentParser(
        description = 'CLOCI locus-locus similarity network. Will output all clans if no -c/-hlg'
        )
    parser.add_argument('-i', '--input', help = 'CLOCI dir', required = True)
    parser.add_argument('-c', '--clans', help = 'Input clans; "-" for stdin')
    parser.add_argument('-hlg', '--hlgs', help = 'Input HLGs; "-" for stdin')
    parser.add_argument('-l', '--locids', help = 'LocIDs for highlighting')
    parser.add_argument('-p', '--passing', help = 'GCFs only',
                        action = 'store_true')
    parser.add_argument('-m', '--min', type = float, default = 0,
        help = '0 < -m < 1: Minimum GCF similarity for network edges')
    parser.add_argument('-d', '--database', help = 'MTDB for taxonomy',
        default = primaryDB())
    parser.add_argument('-li', '--lineages_in',
        help = f'Only consider inputted lineages')
    parser.add_argument('-r', '--rank', 
        help = f'Rank for taxonomy node coloring: {ranks}')
    parser.add_argument('-f', '--function',
        help = f'Ome directory for annotated cluster names from cloci2annotate')
    parser.add_argument('-al', '--annotate_loc', action = 'store_true',
        help = 'Annotate nodes with their locusID')
    parser.add_argument('-ah', '--annotate_hlg', action = 'store_true',
        help = 'Annotate nodes with the HLG')
    parser.add_argument('--gnet', action = 'store_true',
        help = 'Collapse HLGs into single nodes for a binary network')
    parser.add_argument('--convert', action = 'store_true', 
        help = 'Convert annotation IDs into chronological numbers')
    parser.add_argument('--font_size', type = int, default = 20)
    parser.add_argument('--colors', help = 'List of HEX colors for node annotation')
    parser.add_argument('-e', '--extension', help = f'Network file format: {img_exts}, {net_exts}',
        default = 'pdf')
    parser.add_argument('-s', '--scale', type = float, default = 0.1,
        help = 'Edge scaling coefficient; DEFAULT: 0.1')
    parser.add_argument('-o', '--out_name', help = 'Network output prefix')
    args = parser.parse_args()

    db = mtdb(format_path(args.database)).set_index()


    if args.rank:
        rank = args.rank.lower()
        if rank not in set(ranks):
            eprint(f'\nERROR: invalid rank "{rank}"', flush = True)
            sys.exit(3)
    else:
        rank = None

    if args.lineages_in:
        lineages = split_input(args.lineages_in)
        ex_db = extract_mtdb(db, lineage_list = lineages)
        ex_omes = set(ex_db.set_index().keys())
        if not ex_omes:
            eprint(f'\nERROR: lineages did not return omes "{lineages}"', 
                   flush = True)
            sys.exit(53)
    else:
        ex_omes = set()

    if args.function:
        func_dir = format_path(args.function)
    else:
        func_dir = None

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

    if args.colors:
        colors_p = args.colors.replace('"','').replace("'",'').replace(',', ' ').split()
        colors = ['#' + color if not color.startswith('#') else color for color in colors_p]
    else:
        colors = []

    if args.out_name:
        prefix = args.out_name
    else:
        prefix = None    

    main(clans, hlgs, format_path(in_dir, force_dir = True), 
         args.min, args.passing, db, rank, ext = ext, img = img,
         annotate_loc = args.annotate_loc, annotate_hlg = args.annotate_hlg,
         highlight = locids, scale = args.scale, chrono = args.convert,
         font_size = args.font_size, colors = colors, func_dir = func_dir,
         ex_omes = ex_omes, out_prefix = prefix, binary = args.gnet,
         bcentrality = args.gnet)


if __name__ == '__main__':
    cli()
