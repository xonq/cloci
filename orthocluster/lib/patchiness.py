def patch_main( 
    phylo, hgx2omes, hgxs, wrk_dir, 
    old_path = 'patchiness.scores.pickle', cpus = 1
    ):
    
    if not os.path.isfile(wrk_dir + old_path):
        if isinstance(hgx2omes, dict): # round 1 patchiness
            clusOmes = set([tuple([str(x) for x in hgx2omes[y]]) for y in hgxs])
        else: # round 2 is a list
            clusOmes = set([
                tuple([str(x) for x in y]) for y in hgxs
                ])
#        more = set([tuple(omes) for omes in moduleOmes])
#        allHGxs = list(clusOgxs.union(more))    
        with mp.get_context('fork').Pool(processes = cpus) as pool:
            patch_res = pool.starmap(
                calc_patchiness, [(phylo, x) for x in clusOmes]
                )           
        pool.join()         
        omes2patch = {ome_tup: patchiness for ome_tup, patchiness in patch_res}
                            
        with open(wrk_dir + old_path, 'wb') as out:
            pickle.dump(omes2patch, out)
    
    else:
        print('\tLoading previous patchiness results', flush = True)
        with open(wrk_dir + old_path, 'rb') as in_pick:
            omes2patch = pickle.load(in_pick)
        
    return omes2patch
