import os
import re
import sys
import shutil
import pickle
import numpy as np
import multiprocessing as mp
from ast import literal_eval
from tqdm import tqdm
from datetime import datetime
from collections import defaultdict
from scipy import sparse
try:
    from numba import njit
    from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
    @njit
    def gen_clusters_arr(loci, i2hg, clan_arr, check):
        """input is an HGp's 2D array of loci represented as a array of HGs
        presence absence.
        identify rows with > 2 overlapping HGs by multiplying the inputted
        2D np.array against its transposed self.
        if the overlapping loci are from different omes, report their
        overlapping HGs by pairwise loci summation and extracting coordinates
        with HG == 2 - indicating overlap.
        output these overlapping HGs (HGx), omes, and loci anchor genes""" 
        output = []
      #  adj_arr = clan_arr @ clan_arr.transpose() # sparse
        # e.g. hgs x loc     loc x hgs      loc x loc
        # [[1, 0, 1, 1],   [[1, 1, 1],    [[3, 3, 2],
        #  [1, 1, 1, 1], @  [0, 1, 1],  =  [3, 4, 3],
        #  [1, 1, 0, 1]]    [1, 1, 0],     [2, 3, 3]]
        #                   [1, 1, 1]]
        
        #  clan_arr         clan_arr.T      adj_arr
    
        # so in this example, loc 0 has an hgx (overlap > 2) with loc 0 and 1
        # loc 1 hgx with loc 0, 1, 2; and loc 2 an hgx with loc 1, 2
    
        # removing same locus and redundant overlap, there is thus a shared hgx
    
    
        adj_arr = clan_arr @ clan_arr.T # dense identify overlap
        for i0 in range(adj_arr.shape[0]): # for each row
            row = adj_arr[i0, :] # overlap row
            row0 = clan_arr[i0, :] # hg representation of locus
            loc0 = loci[i0] # the actual locus 
            ome0 = loc0[:loc0.find('_')] # the ome corresponding with the locus
    
            # these are loci coordinates non-self and non-redundant w/> 2 overlap
            overlap_loci = [x for x in np.where(row > 2)[0] if x > i0] # dense
    
            # for each of these overlaps
            for i1 in overlap_loci:
                loc1 = loci[i1]
                ome1 = loc1[:loc1.find('_')]
                if ome0 != ome1:
                    row1 = clan_arr[i1, :]
                    row_sum = row0 + row1
                    # if the sum of the row is 2, it indicates there is overlap
                    # at that particular HG
                    hgx = [i2hg[x] for x in np.where(row_sum >= 2)[0]] # dense
                    output.append((hgx, loc0, loc1))
        return output
except ImportError: # no numba 
    def gen_clusters_arr(loci, i2hg, clan_arr, check):
        """input is an HGp's 2D array of loci represented as a array of HGs
        presence absence.
        identify rows with > 2 overlapping HGs by multiplying the inputted
        2D np.array against its transposed self.
        if the overlapping loci are from different omes, report their
        overlapping HGs by pairwise loci summation and extracting coordinates
        with HG == 2 - indicating overlap.
        output these overlapping HGs (HGx), omes, and loci anchor genes""" 
        output = []
      #  adj_arr = clan_arr @ clan_arr.transpose() # sparse
        # e.g. hgs x loc     loc x hgs      loc x loc
        # [[1, 0, 1, 1],   [[1, 1, 1],    [[3, 3, 2],
        #  [1, 1, 1, 1], @  [0, 1, 1],  =  [3, 4, 3],
        #  [1, 1, 0, 1]]    [1, 1, 0],     [2, 3, 3]]
        #                   [1, 1, 1]]
        
        #  clan_arr         clan_arr.T      adj_arr
    
        # so in this example, loc 0 has an hgx (overlap > 2) with loc 0 and 1
        # loc 1 hgx with loc 0, 1, 2; and loc 2 an hgx with loc 1, 2
    
        # removing same locus and redundant overlap, there is thus a shared hgx
    
    
        adj_arr = clan_arr @ clan_arr.T # dense identify overlap
        for i0 in range(adj_arr.shape[0]): # for each row
            row = adj_arr[i0, :] # overlap row
            row0 = clan_arr[i0, :] # hg representation of locus
            loc0 = loci[i0] # the actual locus 
            ome0 = loc0[:loc0.find('_')] # the ome corresponding with the locus
    
            # these are loci coordinates non-self and non-redundant w/> 2 overlap
            overlap_loci = [x for x in np.where(row > 2)[0] if x > i0] # dense
    
            # for each of these overlaps
            for i1 in overlap_loci:
                loc1 = loci[i1]
                ome1 = loc1[:loc1.find('_')]
                if ome0 != ome1:
                    row1 = clan_arr[i1, :]
                    row_sum = row0 + row1
                    # if the sum of the row is 2, it indicates there is overlap
                    # at that particular HG
                    hgx = [i2hg[x] for x in np.where(row_sum >= 2)[0]] # dense
                    output.append((hgx, loc0, loc1))
        return output


import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
from itertools import combinations, chain
from mycotools.lib.biotools import gff2list
from mycotools.lib.kontools import write_json, collect_files, read_json
from cloci.lib import treecalcs
from cloci.lib.input_parsing import compileCDS

def form_hgpDict(out_hgs):

    hgs = sorted(set(chain(*out_hgs)))
    hgp_dict = defaultdict(set)
    for hg0, hg1 in out_hgs:
        hgp_dict[hg0].add(hg1) # {og0: set(og1, ..., ogn), } or 
        # {og0: set(og0, og1, ..., ogn)}    
        hgp_dict[hg1].add(hg0)
                                            
    return hgp_dict


def hash_protohgxs(gff_path, hgp_dict, gene2hg, clusplusminus = 10):
    """parse a gff and compile its organized CDS_dict. identify where og-pairs
    co-occur and retrieve the set of hgs for the locus with the locus seed
    protein"""

    gff_list, protohgx = gff2list(gff_path), defaultdict(list)
    cds_dict = compileCDS(gff_list, os.path.basename(gff_path).replace('.gff3',''))

    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]): # for each index and protein
            try:
                og0 = gene2hg[seq0]
            except KeyError: # gene not in OG / is a singleton
                continue
            if og0 in hgp_dict: # if this og is part of a significant seed
                if i0 < clusplusminus: # check if index - clusplusminus < 0
                    locus = cds_dict[scaf][:i0+clusplusminus+1] # obtain locus
                    # as all the sequences from the beginning of the contig to
                    # i0 plus the locus size tolerance (add 1 for python)
                else:
                    locus = cds_dict[scaf][i0-clusplusminus:i0+clusplusminus+1]
                    # otherwise obtain the locus from the +/-
                for seq1 in locus:
                    if seq1 != seq0: # if the protein isn't the first queried
                        try:
                            og1 = gene2hg[seq1]
                        except KeyError:
                            continue
                        if og1 in hgp_dict[og0]: # if it is in the set of
                        # hgs for the first queried og
                            og_loc_list = []
                            for gene in locus:
                                try:
                                    og_loc_list.append(gene2hg[gene])
                                except KeyError: # singleton or missing
                                    pass
                            og_loc = set(og_loc_list)
                            # the set of hgs in this locus
                            hgp = tuple(sorted([og0, og1])) # sort the hgs
                            protohgx[hgp].append([og_loc, seq0, seq1])
                            # {(og0, og1}: [[set(hgx, .. ogy), seq0, seq1], ...]}

    return dict(protohgx)


def merge_protos(protohgx_res):
    """combine protohgx results from each organism"""

    hgps2hgall = defaultdict(list)
    protohgx2omes = defaultdict(list)
    for ome_protohgx in protohgx_res: # for each protohgx from the mp
    # results
        for hgp, loc_d in ome_protohgx.items(): 
            protohgx2omes[hgp].extend(loc_d)
            for loc in loc_d:
                hgps2hgall[hgp].extend(list(loc[0]))
                    
            # {(hg0, hg1)}: {hg: i}}
    hgp2hgs = {hgp: {v: i \
                     for i, v in enumerate(sorted(set(hgs)))} \
               for hgp, hgs in hgps2hgall.items()}

    return {k: tuple(v) for k,v in protohgx2omes.items()}, hgp2hgs



def load_prehgx(wrk_dir, ome2i):
    h2l = defaultdict(list)
    pre_files = [os.path.basename(x) for x in collect_files(wrk_dir, 'txt')]
    preh_fs = [x for x in pre_files \
               if x.startswith('prehgx2data')]
    for preh_f in tqdm(preh_fs, total = len(preh_fs)):
        with open(f'{wrk_dir}{preh_f}', 'r') as raw:
            for line in raw:
                d = line.rstrip().split()
                hgx = tuple((int(x) for x in d[0].split(',')))
                gs = d[1].split(',')
                h2l[hgx].extend(gs)

    h2l = {k: tuple(sorted(set(v))) for k, v in h2l.items()}
    h2o = {k: tuple(sorted(set(ome2i[y[:y.find('_')]] for y in v))) \
           for k, v in h2l.items()}

    return h2o, h2l


def hgx2mngr(Q, wrk_dir, suffix):
    """Receive from gen_clusters and formulate_clus_dict shared Manager Queue
    Output the data to drive storage for parsing and collection later"""
    x = True
    hgp_file = f'{wrk_dir}prehgx.{suffix}.txt'
    dat_file = f'{wrk_dir}prehgx2data.{suffix}.txt'
    with open(hgp_file, 'a') as hgp_out, open(dat_file, 'a') as dat_out:
        while x:
            x = Q.get()
            if x:
                for mhgx, g0, g1 in x[1]:
                    out_str = ','.join([str(x) for x in mhgx])
                    out_str += f' {g0},{g1}\n'
                    dat_out.write(out_str)
    
                hgp = x[0]
                hgp_out.write(f'{hgp[0]} {hgp[1]}' + '\n')
                hgp_out.flush()
                dat_out.flush()
    Q.put(None)


def formulate_hgx_dicts(hgp, loci, hgs2i, Q):
    hgp_set = set(hgp)
    # make a matrix to identify overlap hgs v loci
    # loci must a be sorted by hgs
    clan_arr = np.zeros([len(loci), len(hgs2i)],
                        dtype = np.int8)
    for i, loc in enumerate(loci):
        hg_loc = np.array(list([hgs2i[x] for x in loc[0]]))
        clan_arr[i, hg_loc] = 1
    loci = [x[1] for x in loci]

#    if len(set(hgs2i.keys()).intersection({10, 1025, 6380, 10267, 12567, 63630})) > 2:
 #       print(loci, flush = True)
    check = False
    if hgp in {(10, 10267), (10, 12567), (10267, 12567), (1025, 63630), (10, 63630)}:
        check = True
    i2hg = np.array([k for k, v in hgs2i.items()])
    output = gen_clusters_arr(loci, i2hg, clan_arr.astype(float), check)
#                output.append([hgx, (ome0, ome1), (loc0, loc1)])
    Q.put((hgp, tuple((tuple(x[0]), x[1], x[2],) \
           for x in output),))


def par_rm(hgx, hgx2omes):

#    check_set, out_p = {10, 1025, 6380, 10267, 12567, 63630}, False
    t_comb_sets, t_combs, todel = [set(hgx)], [hgx], []
    # for each length of combinations possible
    for comb_len in reversed(range(3, len(hgx))):
        # for each combination at that length
        for comb in combinations(hgx, comb_len):
            # if this combination already exists 
            if comb in hgx2omes:
                # get the set of this subhgx and check if its a subset
                set_comb = set(comb)
                for i, set_hgx in enumerate(t_comb_sets):
                    # if it is a subset
                    if set_comb < set_hgx:
                        # and if it has the same omes as the superset
                        if hgx2omes[comb] == hgx2omes[t_combs[i]]:
                            todel.append(comb)
                            break
                    t_combs.append(comb)
                    t_combs_sets.append(set_comb)

    return todel


# I think most subsets will be accounted from the initial comparison,
# and the ones that aren't weren't part of some significant HGp.... soo
# they weren't significant - this would be the most cautious approach,
# but a lot of overhead for something

# in other words, this would only discover few HGxs that are part of
# loci with an HGp that is not part of the HGx; so it recoups some
# of the significant HGxs overlooked from the HGp seeding algorithm,
# but at the end of the day HGxs that are part of nonsignificant seed
# loci are still overlooked. 

# the subset HGxs that do exist will be populated because the loci will
# be overlapped - the performance hit doesnt justify this. if you want
# these loci then decrease the HGp percentile
def add_subsets(hgx2omes, hgx2loc):
    for mhgx in list(hgx2omes.keys()):
        omes = hgx2omes[mhgx]
        for c_len in range(3, len(mhgx)):
            for hgx in combinations(mhgx, c_len):
                hgx2omes[hgx].update(omes)
                hgx2loc[hgx].update(hgx2loc[hgx])

    return hgx2omes, hgx2loc

def rm_subsets(hgx2omes, hgx2loc, cpus = 1):

#    hgx2omes, hgx2loc = add_subsets(hgx2omes, hgx2loc)
    max_len = max((len(x) for x in list(hgx2omes.keys())))
    for hgx_len in tqdm(reversed(range(4, max_len + 1)),
                        total = max_len + 1 - 4):
        i2hgx = list(hgx2omes.keys())
        rm_cmds = [[hgx, hgx2omes] for hgx in i2hgx if len(hgx) == hgx_len]
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            todels = pool.starmap(par_rm, rm_cmds)
            pool.close()
            pool.join()
        
        for todel in todels:
            for hgx in todel:
                try:
                    del hgx2omes[hgx]
                    del hgx2loc[hgx]
                except KeyError:
                    continue

    return hgx2omes, hgx2loc


def id_hgx(db, hgp_dict, gene2hg, ome2i, wrk_dir, cpus, clusplusminus = 10):
    """form hgxs from hgps by extracting higher order combinations from
    overlapping hgx loci"""

    hash_protohgx_cmds = []
    gffs = [v['gff3'] for k, v in db.items() if k in ome2i]
    print('\t\tForming proto-HGx hashes', flush = True)
    for gff in gffs: # prepare for protohgx hashing by organism
        hash_protohgx_cmds.append((gff, hgp_dict, gene2hg, clusplusminus,))
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        protohgx_res = pool.starmap(hash_protohgxs, 
                                    tqdm(hash_protohgx_cmds,
                                         total = len(hash_protohgx_cmds)))
        pool.close()
        pool.join()
                
    print('\t\tForming proto-HGxs', flush = True)
    protohgx2omes, hgp2hgs = merge_protos(protohgx_res) # merge mp results
            
    print('\t\tIdentifying HGxs', flush = True)
    count, ite, proc = 0, 1, None
    
    pre_files = [os.path.basename(x) for x in collect_files(wrk_dir, 'txt')]
    if any(x.startswith('prehgx') for x in pre_files):
        print('\t\t\tLoading checkpoint', flush = True)
        preh_fs = [x for x in pre_files if x.startswith('prehgx.')]
        sxs = [int(re.search(r'prehgx\.(\d+)\.txt', x)[1]) for x in preh_fs]
        sx = max(sxs)
        hgps = []
        for preh_f in tqdm(preh_fs, total = len(preh_fs)):
            with open(f'{wrk_dir}{preh_f}', 'r') as raw:
                for line in raw:
                    d = line.rstrip().split()
                    if len(d) > 1:
                        hgps.append((int(d[0]), int(d[1]),))
        hgps = list(set(hgps))
        for hgp in hgps:
            if hgp in protohgx2omes:
                del protohgx2omes[hgp]
    else:
        sx = -1

    if protohgx2omes:
#        prog = mp.Process(target=progress_report, args=(wrk_dir, 0, start_len))
 #       prog.start()
        # could make this intelligent via a test
        if cpus == 1:
            eprint('ERROR: need more than one CPU', flush = True)
        if cpus == 2:
            c_cpus = 1
            w_cpus = 1
        elif cpus < 3:
            w_cpus = 1
            c_cpus = cpus - 1
        else:
            w_cpus = round(cpus/3)
            c_cpus = cpus - w_cpus
    
        Q = mp.Manager().Queue()
        procs = []
        for i in range(w_cpus):
            sx += 1
            procs.append(mp.Process(target=hgx2mngr, args=(Q, wrk_dir, sx)))
            procs[-1].start()

        params = ((hgp, loci, hgp2hgs[hgp], Q) for hgp, loci in protohgx2omes.items()) 
        with mp.Pool(processes = c_cpus) as pool:
            pool.starmap(formulate_hgx_dicts, tqdm(params, total = len(protohgx2omes)))
            pool.close()
            pool.join()
                
        Q.put(None)
        for proc in procs:
            proc.join()
        Q.get()

   # print('joining', flush = True)
  #  prog.join()
    print('\t\tReading results', flush = True)
    hgx2omes, hgx2loc = load_prehgx(wrk_dir, ome2i)

 #   print('\t\tRemoving subset HGxs', flush = True)
#    hgx2omes, hgx2loc = rm_subsets(hgx2omes, hgx2loc, cpus)
 #   hgx2omes = {x: tuple(sorted(hgx2omes[x])) for x in hgx2omes}

    return hgx2omes, hgx2loc


def hgp2hgx(db, wrk_dir, top_hgs, gene2hg, ome2i, 
                phylo, plusminus = 3, cpus = 1):
    if not os.path.isfile(wrk_dir + 'hgx2loc.pickle'):
        hgp_dict = form_hgpDict(top_hgs)
        print('\tForming HGxs', flush = True)
        form_clus_start = datetime.now()
        hgx2omes, hgx2loc = id_hgx(
            db, hgp_dict, gene2hg,  
            ome2i, wrk_dir, cpus, clusplusminus = plusminus
            )
#        write_json({str(k): list(v) for k, v in hgx2loc.items()},
 #                  f'{wrk_dir}hgx2loc.json.tmp.gz')
        with open(wrk_dir + 'hgx2loc.pickle', 'wb') as pickout:
            pickle.dump(hgx2loc, pickout)
 #       write_json({str(k): list(v) for k, v in hgx2omes.items()},
  #                 f'{wrk_dir}hgx2omes.json.tmp.gz')
   #     shutil.move(f'{wrk_dir}hgx2loc.json.tmp.gz', 
    #                f'{wrk_dir}hgx2loc.json.gz')
     #   shutil.move(f'{wrk_dir}hgx2omes.json.tmp.gz',
      #              f'{wrk_dir}hgx2omes.json.gz')
        with open(wrk_dir + 'hgx2omes.pickle', 'wb') as pickout:
            pickle.dump(hgx2omes, pickout)
        formHGxTime = datetime.now() - form_clus_start
        pre_files = [os.path.basename(x) for x in collect_files(wrk_dir, 'txt')]
        todel = [x for x in pre_files if x.startswith('prehgx')]
        for f in todel:
            os.remove(wrk_dir + f)
        print('\t\t' + str(formHGxTime), flush = True)
    
    else: # or just load available structures 
#        hgx2loc = literal_load(f'{wrk_dir}hgx2loc.json.gz')
 #       hgx2omes = literal_load(f'{wrk_dir}hgx2omes.json.gz')
        with open(wrk_dir + 'hgx2loc.pickle', 'rb') as pickin:
            hgx2loc = pickle.load(pickin)
        hgx2omes = {k: tuple(sorted(set(ome2i[y[:y.find('_')]] for y in v))) \
                    for k, v in hgx2loc.items()}
#        with open(wrk_dir + 'hgx2omes.pickle', 'rb') as pickin:
 #           hgx2omes = pickle.load(pickin)

    return hgx2omes, hgx2loc
