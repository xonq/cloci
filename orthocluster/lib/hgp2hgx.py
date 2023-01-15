import os
import re
import shutil
import pickle
import numpy as np
import multiprocessing as mp
from ast import literal_eval
from datetime import datetime
from collections import defaultdict
from scipy import sparse
from itertools import combinations, chain
from mycotools.lib.biotools import gff2list
from mycotools.lib.kontools import write_json, collect_files, read_json
from orthocluster.orthocluster.lib import phylocalcs
from orthocluster.orthocluster.lib.input_parsing import compileCDS

def form_hgpDict(out_hgs):
    
    hgs = set([int(x[0]) for x in out_hgs]) # make a set of the first og in an
    # og-pair
    hgs = hgs.union(set([int(x[1]) for x in out_hgs])) # add in the second og
    hgp_dict = {x: set() for x in hgs} # prepare the keys to bypass `if`'s
    for i in out_hgs:
        hgp_dict[int(i[0])].add(int(i[1])) # {og0: set(og1, ..., ogn), } or 
        # {og0: set(og0, og1, ..., ogn)}    
        hgp_dict[int(i[1])].add(int(i[0]))
                                            
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
            # {(hg0, hg1)}: [[{hgy, hgz}, seq], ...]

    hgp2hgs = {hgp: {v: i \
                     for i, v in enumerate(sorted(set(hgs)))} \
               for hgp, hgs in hgps2hgall.items()}

    return dict(protohgx2omes), hgp2hgs


def gen_clusters(loc0, loc1, hgp_set):

    ome0 = loc0[1][:loc0[1].find('_')] # extract ome from gene acc
    ome1 = loc1[1][:loc1[1].find('_')]
    if ome0 != ome1:
        loc0vloc1 = loc0[0].intersection(loc1[0])
        if loc0vloc1 != hgp_set: #if there are more overlapping hgs
            hgx = tuple(sorted(list(loc0vloc1))) # sort the loci (not by
            return hgx, ome0, ome1, loc0[1], loc1[1]

def gen_clusters_arr(loci, hgs2i, clan_arr, ome2i):

    output = []
    clan_arr = clan_arr.toarray() # to dense
    i2hg = [k for k, v in hgs2i.items()]
  #  adj_arr = clan_arr @ clan_arr.transpose() # sparse
    adj_arr = np.matmul(clan_arr, clan_arr.T) # dense
    for i0 in range(adj_arr.shape[0]):
        row = adj_arr[i0, :]
        row0 = clan_arr[i0, :]
        loc0 = loci[i0][1]
        ome0 = loc0[:loc0.find('_')]
        # identify loci with > 2 overlapping hgs
        prep = row > 2
#        overlap_loci = [x for x in sparse.find(prep)[0]] # sparse
        overlap_loci = [x for x in np.where(row > 2)[0] if x > i0] # dense
        for i1 in overlap_loci:
            loc1 = loci[i1][1]
            ome1 = loc1[:loc1.find('_')]
            if ome0 != ome1:
                row1 = clan_arr[i1, :]
                row_sum = row0 + row1
    #            where = row_sum == 2
   #             hgx = tuple([i2hg[x] for x in sparse.find(where)[0]])
                hgx = tuple([i2hg[x] \
                             for x in np.where(row_sum == 2)[0]]) # dense
                output.append([hgx, (ome2i[ome0], ome2i[ome1]), (loc0, loc1)])

    return output


def formulate_hgx_dict_solo(cmd, hgx2omes, hgx2loc, ome2i):

    for i, res in enumerate(gen_clusters_arr(*cmd)):
        if res:
            for comb_len in range(3, len(res[0]) + 1):
                hgx = res[0]
                hgx2omes[hgx] = hgx2omes[hgx].union({
                     ome2i[res[1]], ome2i[res[2]]
                     })
                hgx2loc[hgx] = hgx2loc[hgx].union({res[3], res[4]})

    return hgx2omes, hgx2loc


def add2hgxdicts(Q, rQ):
    hgx2omes, hgx2loc = defaultdict(set), defaultdict(set)
    x = True
    while x:
        x = Q.get()
        if x:
            mhgx, i0, i1, g1, g2 = x
            hgx2omes[mhgx].update((
                i0, i1
                ))
            hgx2loc[mhgx].update((
                g1, g2
                ))
    rQ.put([hgx2omes, hgx2loc])


def literal_load(f):
    d = read_json(f)
    return {literal_eval(k): literal_eval(v) for k, v in d.items()}


def load_prehgx(wrk_dir, suffix = ''):
    h2o = literal_load(f'{wrk_dir}prehgx2omes{suffix}.json.gz')
    h2l = literal_load(f'{wrk_dir}prehgx2loc{suffix}.json.gz')
    with open(f'{wrk_dir}prehgx{suffix}.txt', 'r') as raw:
        phgps = [x.rstrip().replace('(','').replace(')','').replace(',','').split() \
                for x in raw if x.rstrip()]
    hgps = [(int(x[0]), int(x[1])) for x in phgps]
    return defaultdict(set, h2o), defaultdict(set, h2l), hgps


def write_prehgx(hgx2omes, hgx2loc, hgps, wrk_dir, suffix):
    write_json({str(k): str(v) for k,v in hgx2omes.items()}, 
               f'{wrk_dir}prehgx2omes{suffix}.json.tmp.gz')
    write_json({str(k): str(v) for k,v in hgx2loc.items()}, 
               f'{wrk_dir}prehgx2loc{suffix}.json.tmp.gz')
    shutil.move(f'{wrk_dir}prehgx2omes{suffix}.json.tmp.gz', 
                f'{wrk_dir}prehgx2omes{suffix}.json.gz')
    shutil.move(f'{wrk_dir}prehgx2loc{suffix}.json.tmp.gz', 
                f'{wrk_dir}prehgx2loc{suffix}.json.gz')
    with open(f'{wrk_dir}prehgx{suffix}.txt', 'a') as out:
        out.write('\n'.join([str(x) for x in hgps]) + '\n')

def formulate_hgx_dicts(gen_clus_cmds, hgx2omes, hgx2loc, prehgx, wrk_dir, suffix = '', cpus = 1):

    with mp.get_context('fork').Pool(processes = cpus) as pool:
        clus_res = pool.starmap(gen_clusters_arr, gen_clus_cmds)
    clus_res = chain(*clus_res)
    for res in clus_res:
        hgx, iS, gs = res
        hgx2omes[hgx].update(iS)
        hgx2loc[hgx].update(gs)

    proc = mp.Process(target=write_prehgx, args=[hgx2omes, hgx2loc, prehgx, wrk_dir, suffix])
    proc.start()

    return hgx2omes, hgx2loc, proc

def par_rm(hgx, hgx2omes):

    t_comb_sets, t_combs, todel = [set(hgx)], [hgx], []
    for comb_len in reversed(range(3, len(hgx))):
        for comb in combinations(hgx, comb_len):
            if comb in hgx2omes:
                set_comb = set(comb)
                for i, set_hgx in enumerate(t_comb_sets):
                    if set_comb < set_hgx:
                        if hgx2omes[comb] == hgx2omes[t_combs[i]]:
                            todel.append(comb)
                            break
                    t_combs.append(set_comb)
    
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
    max_len = len(list(hgx2omes.keys())[0])
    for hgx_len in reversed(range(4, max_len + 1)):
        i2hgx = list(hgx2omes.keys())
        rm_cmds = [[hgx, hgx2omes] for hgx in i2hgx if len(hgx) == hgx_len]
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            todels = pool.starmap(par_rm, rm_cmds)
        
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
        protohgx_res = pool.starmap(hash_protohgxs, hash_protohgx_cmds)
    pool.join()
                
    print('\t\tForming proto-HGxs', flush = True)
    protohgx2omes, hgp2hgs = merge_protos(protohgx_res) # merge mp results
            
    print('\t\tIdentifying HGxs', flush = True)
    print(f'\t\t\t{len(protohgx2omes)} to run', flush = True)
    count, ite, proc = 0, 1, None
                    
    # generate dictionaries of crude predicted hgxs
    gen_clus_cmds, prehgxs = [], []
    hgx2omes, hgx2loc = defaultdict(set), defaultdict(set)
    # need at least two cpus at this point
#    if cpus == 2:
     #   d_cpu = 1
    #    m_cpu = 1
   # elif cpus == 3:
  #      d_cpu = 1
 #       m_cpu = 2
#    else:
    #    m_cpu = round(0.7*cpus)
   #     d_cpu = cpus - m_cpu
  #  Q = mp.Manager().Queue()
 #   rQ = mp.Manager().Queue()
#    procs = [mp.Process(target=add2hgxdicts, args=[Q, rQ]) for i in range(d_cpu)]
#    for proc in procs:
 #       proc.start()
    # load checkpoint if available
    if os.path.isfile(wrk_dir + 'prehgx.txt'):
        print('\t\t\tLoading checkpoint', flush = True)
        pre_files = [os.path.basename(x) for x in collect_files(wrk_dir, 'txt')]
        preh_fs = [re.search(r'prehgx(\d*)\.txt', x)[1]
                   for x in pre_files \
                   if re.search(r'prehgx\d*\.txt', x) is not None]
        try:
            sx = max([int(x) for x in preh_fs \
                      if x and os.path.isfile(f'{wrk_dir}prehgx2loc{x}.json.gz')])
        except ValueError:
            sx = ''

        hgx2omes, hgx2loc, ran_hgp = load_prehgx(wrk_dir, sx)
        for hgp in ran_hgp:
            if hgp in protohgx2omes:
                del protohgx2omes[hgp]

        print(f'\t\t\t{len(protohgx2omes)} to run', flush = True)

        if sx:
            sx += 1
        else:
            sx = 1

        with open(f'{wrk_dir}prehgx{sx}.txt', 'w') as out:
            out.write('\n'.join([str(x) for x in ran_hgp]))
    else:
        sx = ''


    for hgp, loci in protohgx2omes.items():
        hgp_set = set(hgp)
        hgs2i = hgp2hgs[hgp]
        # make a matrix to identify overlap hgs v loci
        # loci must a be sorted by hgs
        clan_arr = np.zeros([len(loci), len(hgs2i)],
                            dtype = np.int8)
        for i, loc in enumerate(loci):
            hg_loc = np.array(list([hgs2i[x] for x in loc[0]]))
            clan_arr[i, hg_loc] = 1
        count += 1
        gen_clus_cmds.append([loci, hgs2i, sparse.csr_matrix(clan_arr), ome2i]) # sparse
        prehgxs.append(hgp)
#        gen_clus_cmds.append([loci, hgs2i, clan_arr, ome2i]) # dense

        if count > 5000 * ite:
            ite += 1
            print('\t\t\t' + str(count), flush = True)
            if proc:
                proc.join()
            hgx2omes, hgx2loc, proc = formulate_hgx_dicts(
                gen_clus_cmds, hgx2omes, hgx2loc, prehgxs, wrk_dir, cpus = cpus, suffix = sx
                )               
            gen_clus_cmds, prehgxs = [], []
                            
    # run the last set of commands
    if proc:
        proc.join()
    hgx2omes, hgx2loc, proc = formulate_hgx_dicts(
        gen_clus_cmds, hgx2omes, hgx2loc, prehgxs, wrk_dir, cpus = cpus, suffix = sx
        )
    proc.join()

    print('\t\tRemoving subset HGxs', flush = True)
    hgx2omes, hgx2loc = rm_subsets(hgx2omes, hgx2loc, cpus)
    hgx2omes = {x: tuple(sorted(list(hgx2omes[x]))) for x in hgx2omes}

    return hgx2omes, hgx2loc


def hgp2hgx(db, wrk_dir, top_hgs, gene2hg, ome2i, 
                phylo, plusminus = 3, cpus = 1):
    if not os.path.isfile(wrk_dir + 'hgx2omes.json.gz'):
        hgp_dict = form_hgpDict(top_hgs)
        print('\tForming HGxs', flush = True)
        form_clus_start = datetime.now()
        hgx2omes, hgx2loc = id_hgx(
            db, hgp_dict, gene2hg,  
            ome2i, wrk_dir, cpus, clusplusminus = plusminus
            )
        write_json({str(k): str(v) for k, v in hgx2loc.items()},
                   f'{wrk_dir}hgx2loc.json.tmp.gz')
#        with open(wrk_dir + 'hgx2loc.json.tmp.gz', 'wb') as pickout:
 #           write(hgx2loc, pickout)
        write_json({str(k): str(v) for k, v in hgx2omes.items()},
                   f'{wrk_dir}hgx2omes.json.tmp.gz')
        shutil.move(f'{wrk_dir}hgx2loc.json.tmp.gz', 
                    f'{wrk_dir}hgx2loc.json.gz')
        shutil.move(f'{wrk_dir}hgx2omes.json.tmp.gz',
                    f'{wrk_dir}hgx2omes.json.gz')
#        with open(wrk_dir + 'hgx2omes.pickle', 'wb') as pickout:
 #           pickle.dump(hgx2omes, pickout)
        formHGxTime = datetime.now() - form_clus_start
        pre_files = [os.path.basename(x) for x in collect_files(wrk_dir, 'txt')]
        todel = [x for x in pre_files if x.startswith('prehgx')]
        for f in todel:
            os.remove(wrk_dir + f)
        print('\t\t' + str(formHGxTime), flush = True)
    
    else: # or just load available structures 
        hgx2loc = literal_load(f'{wrk_dir}hgx2loc.json.gz')
        hgx2omes = literal_load(f'{wrk_dir}hgx2omes.json.gz')
#        with open(wrk_dir + 'hgx2loc.pickle', 'rb') as pickin:
 #           hgx2loc = pickle.load(pickin)
  #      with open(wrk_dir + 'hgx2omes.pickle', 'rb') as pickin:
   #         hgx2omes = pickle.load(pickin)

    return hgx2omes, hgx2loc
