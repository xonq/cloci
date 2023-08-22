#! /usr/bin/env python3

# hit percentage > 50%, evalue < 0.01
# index 15, 14, 5, and 6
# total genes for taxon of interest, total each pfam for taxon of interest
    # acquire

# NEED to make enrichment analysis consider all Pfam hits for a gene if it isn't already
# NEED top hit argument
# NEED specific GO term argument

from scipy.stats import hypergeom
from mycotools.lib.biotools import fa2dict, gff3Comps, gff2list
from mycotools.lib.kontools import format_path, eprint, mkOutput
from mycotools.lib.dbtools import mtdb
from collections import Counter, defaultdict
from tqdm import tqdm
import multiprocessing as mp, os, argparse, sys, re

def collectGffs(omes, top_dir):
    clusGFFs = {
        ome: top_dir + ome + '/' + ome + '.gff3' for ome in omes \
        if os.path.isfile(top_dir + ome + '/' + ome + '.gff3')
        }
    return tuple(clusGFFs.items())

def parseGff(ome, gff_path, aliasComp = re.compile(gff3Comps()['Alias'])):
    genes = {}
    gff = gff2list(gff_path)
    for entry in gff:
        clusID = re.search(r'clusID=([^;]+)', entry['attributes'])[1]
        alias = aliasComp.search(entry['attributes'])[1]
        genes[alias] = clusID
    return ome, tuple(genes.items())

def extractGO( file_, goTerm = 'GO:0008152' ):
    '''From the GO file, compile all GO terms. For each GO, add to `godict` as the key and the set of GO terms it is as value.
    For each GO term of each key in `godict`, check all of their sub GO terms to see if `goTerm` is within one. If so, append
    it to the output.'''

    goSet = set()
    goSet.add(goTerm)
    count = 1

    godict = {}
    goterm = re.compile(r'GO\:\d+')
    with open( file_, 'r' ) as raw:
        for line in raw:
            if line.startswith('[Term]'):
                goid = None
            elif line.startswith('id:'):
                goid = goterm.search( line )
                if goid:
                    goid = goid[0]
                    godict[goid] = set()
            elif line.startswith('is_a:'):
                goisa = goterm.search( line )
                if goisa:
                    goisa = goisa[0]
                    godict[goid].add( goisa )

    todel = []
    for goid in godict:
        for sub_goid in godict[goid]:
            temp_goSet = set(sub_goid)
            temp_goList = [sub_goid]
            count = 1
            while count == 1:
                del_list = []
                for index in range(len(temp_goList)):
                    del_list.append(index)
                    id1 = str(temp_goList[index])
                    temp_goList.extend(godict[id1])
                    temp_goSet = temp_goSet | set(godict[id1])
                    if goTerm in temp_goSet:
                        goSet = goSet | {goid, sub_goid}
                        count = 0
                        break
                del_list.sort(reverse = True)
                for del_term in del_list:
                    del temp_goList[del_term]
                if len(temp_goList) == 0:
                    count = 0

    return goSet


def parsePfam2GO(pfam2go_path, pfam2def, go_sets):
    pfam2go, go2onto = defaultdict(list), {}
    go2pfam = defaultdict(list)
    with open(pfam2go_path, 'r') as raw:
        for line in raw:
            try:
                pfam = re.search(r'PF\d+', line)[0]
                go = re.search(r'GO\:\d+', line)[0]
                onto = re.search(r'GO\:([^;]+)', line)[1]
                try:
                     definition = pfam2def[pfam]
                except KeyError: # not used, who cares
                    continue
                for go_term, go_set in go_sets.items():
                    if go in go_set:
                        go2pfam[go_term].append(definition)
                pfam2go[definition].append(go)
                go2onto[go] = onto
            except TypeError:
                continue
    return pfam2go, go2onto, dict(go2pfam)
            

def countGenes(ome, proteome_path):
    return ome, len(fa2dict(proteome_path))

def clusterGenes(ome, info_path):
    genes = {}
    try:
        with open(info_path, 'r') as raw:
            for line in raw:
                if line.startswith('#'):
                    continue
                d = line.split('\t')
                clus, genesT = d[0], d[2].split(',')
                genes = {**genes, **{x: clus for x in genesT}}
    except FileNotFoundError:
        pass
    return ome, tuple(genes.items())

def parsePfam(ome, pfam_res, clusGenes, hlg_genes, top_hit = True):
    gene2pfam = defaultdict(dict)
    clusGenes = dict(clusGenes)
    hlgGenes = dict(hlg_genes)
    pfam2def = {}
    try:
        with open(pfam_res, 'r') as raw:
            for line in raw:
                if line.startswith('#'):
                    continue
                data = line.rstrip().split()
                g, p, q, b = data[0], data[4], data[3], float(data[7])
                
                # this is where a post-hoc threshold would be applied
                if q not in gene2pfam[g]: # in case there is a duplicate
                # for the same gene, add in the better score 
                    gene2pfam[g][q] = b
                else:
                    if b > gene2pfam[g][q]:
                        gene2pfam[g][q] = b
                pfam2def[p[:p.find('.')]] = q
    except FileNotFoundError:
        return ome, None, None, None, None
    
    if top_hit:
        gene2pfam = {k: dict(sorted(v.items(), key = lambda x: x[1], reverse = True)) \
                     for k, v in gene2pfam.items()}
        gene2pfam = {k: {list(v.keys())[0]: v[list(v.keys())[0]]} for k, v in gene2pfam.items()}
    gene2pfam, clusGenes2pfam = dict(gene2pfam), {}
    for gene, clus in clusGenes.items():
        try:
            clusGenes2pfam[gene] = tuple([gene2pfam[gene], clus])
        except KeyError: # no pfam acc found, should be rarer
            pass
    hlg_genes2pfam = {}
    for gene, clus in hlgGenes.items():
        try:
            hlg_genes2pfam[gene] = tuple([gene2pfam[gene], clus])
        except KeyError:
            pass
    
    return ome, tuple(gene2pfam.items()), \
        tuple(clusGenes2pfam.items()), tuple(hlg_genes2pfam.items()), tuple(pfam2def.items())
#        tuple([tuple([gene, pfams.items()]) for gene, pfams in clusGenes2pfam.items()])


def calcP(pfam, totCount, clusCount, totalGenes, totalClusGenes, pfamInClus):
    pval = hypergeom.sf(
        clusCount-1, totalGenes, totCount, totalClusGenes
        )
    return pfam, pval, totCount/totalGenes, clusCount/totalClusGenes


def GetPopResultsGO(
    omes2genes2pfams, omes2hlgs2pfam, omes2clus2pfam, totalGenes, total_hlg_genes, totalClusGenes, 
    pfam2go, go2pfam, pool
    ):
    totalPfams, clusPfams, pfam2clusters = Counter(), Counter(), defaultdict(list)
    hlg_pfams, pfam2hlgs = Counter(), defaultdict(list)
    totalGOs, clusGOs, go2clusters = Counter(), Counter(), defaultdict(list)
    hlg_gos, go2hlgs = Counter(), defaultdict(list)
    total_ref_gos, ref_go2clusters = Counter(), defaultdict(list)
    ref_go2hlgs = defaultdict(list)
    sp_gos = Counter()
    for ome, genes2pfam in tqdm(omes2genes2pfams.items(), total = len(omes2genes2pfams)):
        clusGenes2pfam = omes2clus2pfam[ome]
        hlgGenes2pfam = omes2hlgs2pfam[ome]
        tPfams, tClusPfams, t_hlg_pfams = [], [], []
        tGOs, tClusGOs, t_hlg_gos = [], [], []
        t_s_gos, t_s_clus_gos, t_s_hlg_gos = defaultdict(int), defaultdict(int), defaultdict(int)
        # count all instances
        for pfams in list(genes2pfam.values()):
            tPfams.extend(pfams)
            for pfam in pfams:
                go_terms = pfam2go[pfam]
                if not go2pfam:
                    tGOs.extend(go_terms)
                else:
                    for ref_go, pfam_set in go2pfam.items():
                        if pfam in pfam_set:
                            t_s_gos[ref_go] += 1
             
        for pfams, hlg in list(hlgGenes2pfam.values()):
            for pfam in pfams:
                pfam2hlgs[pfam].append(f'{ome}_{hlg}')
                if not go2pfam:
                    for go in pfam2go[pfam]:
                        go2hlgs[go].append(f'{ome}_{hlg}')
                        t_hlg_gos.append(go)
                else:
                    for ref_go, pfam_set in go2pfam.items():
                        if pfam in pfam_set:
                            t_s_hlg_gos[ref_go] += 1
                            ref_go2hlgs[ref_go].append(f'{ome}_{hlg}')
                t_hlg_pfams.append(pfam)   

        for pfams, clus in list(clusGenes2pfam.values()):
#            if any('F-box' in x for x in pfams) and any('WD40' in x for x in pfams):
 #               print(clus, flush = True)
            for pfam in pfams:
                pfam2clusters[pfam].append(ome + '_' + clus)
                if not go2pfam:
                    for go in pfam2go[pfam]:
                        go2clusters[go].append(ome + '_' + clus)
                        tClusGOs.append(go)
                else:
                    for ref_go, pfam_set in go2pfam.items():
                        if pfam in pfam_set:
                            t_s_clus_gos[ref_go] += 1
                            ref_go2clusters[ref_go].append(ome + '_' + clus)
                tClusPfams.append(pfam)
        totalPfams += Counter(tPfams)
        clusPfams += Counter(tClusPfams)
        hlg_pfams += Counter(t_hlg_pfams)
        totalGOs += Counter(tGOs)
        clusGOs += Counter(tClusGOs)
        hlg_gos += Counter(t_hlg_gos)
        total_ref_gos += Counter(t_s_gos)

    cvt_c, cvh_c, hvt_c = [], [], []
#    p_cvt_pval, p_cvh_pval = [], []
    for pfam, clusCount in clusPfams.items():
        totCount = totalPfams[pfam]
        hlgCount = hlg_pfams[pfam]
        pfamInClus = len(set(pfam2clusters[pfam]))
 #       p_cvt_pval.append(calcP(pfam, totCount, clusCount, totalGenes, totalClusGenes, pfamInClus))
  #      p_cvh_pval.append(calcP(pfam, hlgCount, clusCount, total_hlg_genes, totalClusGenes, pfamInClus))
        cvt_c.append([pfam, totCount, clusCount, totalGenes, totalClusGenes, pfamInClus])
        cvh_c.append([pfam, hlgCount, clusCount, total_hlg_genes, totalClusGenes, pfamInClus])
    p_cvt_pval = pool.starmap(calcP, cvt_c)
    p_cvh_pval = pool.starmap(calcP, cvh_c)

#    p_hvt_pval = []    
    for pfam, hlgCount in hlg_pfams.items():
        totCount = totalPfams[pfam]
        hlgCount = hlg_pfams[pfam]
        pfam_in_hlg = len(set(pfam2hlgs[pfam]))
#        p_hvt_pval.append(calcP(pfam, totCount, hlgCount, totalGenes, total_hlg_genes, pfam_in_hlg))
        hvt_c.append([pfam, totCount, hlgCount, totalGenes, total_hlg_genes, pfam_in_hlg])
    p_hvt_pval = pool.starmap(calcP, hvt_c)

    cvt_c, cvh_c, hvt_c = [], [], []
    if not go2pfam:
        for go, clusCount in clusGOs.items():
            totCount = totalGOs[go]
            hlgCount = hlg_gos[go]
            goInClus = len(set(go2clusters[go]))
            cvt_c.append([go, totCount, clusCount, totalGenes, totalClusGenes, goInClus])
            cvh_c.append([go, hlgCount, clusCount, total_hlg_genes, totalClusGenes, goInClus])
        g_cvt_pval = pool.starmap(calcP, cvt_c)
        g_cvh_pval = pool.starmap(calcP, cvh_c)
        for go, hlgCount in hlg_gos.items():
            totCount = totalGOs[go]
            hlgCount = hlg_gos[go]
            go_in_hlg = len(set(go2hlgs[go]))
            hvt_c.append([go, totCount, hlgCount, totalGenes, total_hlg_genes, go_in_hlg])
        g_hvt_pval = pool.starmap(calcP, hvt_c)
        r_cvt_pval = g_cvt_pval
        r_cvh_pval = g_cvh_pval
        r_hvt_pval = g_hvt_pval
    else: 
 #       r_cvt_pval, r_cvh_pval, r_hvt_pval = [], [], []
        for ref_go, gos in ref_go2clusters.items():
            totCount = total_ref_gos[ref_go]
            ref_in_clus = len(set(gos))
            ref_in_hlg = len(set(ref_go2hlgs[ref_go]))
 #           r_cvt_pval.append(calcP(ref_go, totCount, ref_in_clus, totalGenes, totalClusGenes, ref_in_clus))
#            r_cvh_pval.append(calcP(ref_go, ref_in_hlg, ref_in_clus, total_hlg_genes, totalClusGenes, ref_in_clus))
            cvt_c.append([ref_go, totCount, ref_in_clus, totalGenes, totalClusGenes, ref_in_clus])
            cvh_c.append([ref_go, ref_in_hlg, ref_in_clus, total_hlg_genes, totalClusGenes, ref_in_clus])
        r_cvt_pval = pool.starmap(calcP, cvt_c)
        r_cvh_pval = pool.starmap(calcP, cvh_c)
    
        for ref_go, gos in ref_go2hlgs.items():
            totCount = total_ref_gos[ref_go]
            ref_in_hlg = len(set(gos))
#            r_hvt_pval.append(calcP(ref_go, totCount, ref_in_hlg, totalGenes, total_hlg_genes, ref_in_hlg))
            hvt_c.append([ref_go, totCount, ref_in_hlg, totalGenes, total_hlg_genes, ref_in_hlg])
        r_hvt_pval = pool.starmap(calcP, hvt_c)
    

    return [p_cvt_pval, p_cvh_pval, p_hvt_pval], \
           [None, None, None], \
           [r_cvt_pval, r_cvh_pval, r_hvt_pval]
#           [g_cvt_pval, g_cvh_pval, g_hvt_pval], \



# need total genes, per pfam counts, pfams in clusters counts, genes 
def GetPopResults(omes2genes2pfams, omes2clus2pfam, totalGenes, totalClusGenes, pool):
    totalPfams, clusPfams, pfam2clusters = Counter(), Counter(), defaultdict(list)
    for ome, genes2pfam in omes2genes2pfams.items():
        clusGenes2pfam = omes2clus2pfam[ome]
        tPfams, tClusPfams = [], []
        for pfams in list(genes2pfam.values()):
            tPfams.extend(pfams)
        for pfams, clus in list(clusGenes2pfam.values()):
            for pfam in pfams:
#                print(pfam, flush = True)
                if pfam == 'DUF3435' or pfam.lower() == 'plavaka':
                    print(ome, clus, pfam, flush = True)
                pfam2clusters[pfam].append(ome + '_' + clus)
                tClusPfams.append(pfam)
        totalPfams += Counter(tPfams)
        clusPfams += Counter(tClusPfams)

    cmds = []
    for pfam, clusCount in clusPfams.items():
        totCount = totalPfams[pfam]
        pfamInClus = len(set(pfam2clusters[pfam]))
        cmds.append([pfam, totCount, clusCount, totalGenes, totalClusGenes, pfamInClus])

    pvalRes = pool.starmap(calcP, cmds)

    return pvalRes

def WriteOutput(scoresDict, thresh, out_file, go2onto = None, go2pfam = False, sort = False):
# out_dir + taxon.lower() + '.scores.tsv'
    try:
        alpha = thresh/len(scoresDict) # bonferroni correction
    except ZeroDivisionError:
        eprint('\t\tWARNING: no hits', flush = True)
        return None

    sig, insig = [], []
    with open(out_file + '.tmp', 'w') as out:
        if go2onto:
            out.write(
                '#go\ttotal\tclusters\tp\tsignificance (alpha =' + \
                str(alpha) + ')\tontology\n'
                )
        elif go2pfam:
            out.write(f'#go\tcompared\tcomparee\tp\tsignificance (alpha ={alpha})\n')
        else:
            out.write('#pfam\tcompared\tcomparee\tp\tsignificance (alpha=' + str(alpha) + ')\n')
        for q, val, tot, clusTot in scoresDict:
            if val <= alpha: # significant
                out.write(f'{q}\t{tot}\t{clusTot}\t{val}\t1')
#                sig.append([q, tot, clusTot, val, 1])
            else: # not significant
                out.write(f'{q}\t{tot}\t{clusTot}\t{val}\t0')
            if go2onto:
                out.write(f'\t{go2onto[q]}\n')
            else:
                out.write('\n')

    os.rename(out_file + '.tmp', out_file)

    if sort:
        sig = []
        insig = []
        with open(out_file, 'r') as raw:
            for line in raw:
                d = line.rstrip().split()
                if d[4] == '1':
                    sig.append(d)
                else:
                    insig.append(d)
        with open(out_file + '.tmp', 'w') as out: 
            sig = sorted(sig, key = lambda x: float(x[2]), reverse = True)
            out.write('\n'.join(['\t'.join([str(x) for x in y]) for y in sig]) + '\n')
            out.write('\n'.join(['\t'.join([str(x) for x in y]) for y in insig]))
        os.rename(out_file + '.tmp', out_file)


def main(
    db, thresh, rank, ome_dir, pfam_dir, out_dir,
    pfam2go_path, go_file, go_terms = [], top_hit = True, 
    pool = mp.Pool(processes = 1), minimum = 0
    ):

    db = db.set_index()
    ranks = defaultdict(list)
    for ome, row in db.items():
        taxon = row['taxonomy'][rank]
        ranks[taxon].append(ome)

    ranks = {k: v for k, v in ranks.items() if len(v) > minimum}

    print('\nChecking input', flush = True)
    for taxon, omes in ranks.items():
        todel = []
        for i, ome in enumerate(omes):
 #           if ome_dir:
            if not os.path.isfile(ome_dir + ome + '/hlg.tsv'):
                todel.append(i)
            elif not os.path.isfile(ome_dir + ome + '/gcf.tsv'):
                eprint(f'\tWARNING: {ome} gcf.tsv missing, assuming no GCFs', flush = True)
#                sys.exit(14)
#            elif not os.path.isfile(as_dir + ome + '/' + ome + '.gff3'):
 #               todel.append(i)
        for i in reversed(todel):
            del omes[i]

    go_sets = {}
    for go_term in go_terms:
        go_sets[go_term] = extractGO(go_file, go_term)

    ranks = {k: v for k, v in sorted(ranks.items(), key = lambda x: len(x[1]))}
    for taxon, omes in ranks.items():
        if not omes:
            continue
        print('\n' + taxon + ': ' + str(len(omes)) + ' omes', flush = True)
        print('\tGenes', flush = True)
        taxonGenes, taxon_hlg_genes, taxonClusGenes = 0, 0, 0
        countres = pool.starmap(
            countGenes, 
            [[ome, db[ome]['faa']] for ome in omes]
            )
        for ome, tot in countres:
            taxonGenes += tot

      #  if not as_dir:
        print('\tHLGs', flush = True)
        hlgres = pool.starmap(
            clusterGenes,
            [[ome, ome_dir + ome + '/hlg.tsv'] for ome in omes]
            )

        print('\tGCFs', flush = True)
        gcfres = pool.starmap(
            clusterGenes,
            [[ome, f'{ome_dir}{ome}/gcf.tsv'] for ome in omes]
            )
#        elif as_dir:
 #           gffs = collectGffs(omes, as_dir)
  #          clusres = pool.starmap(
   #             parseGff,
    #            [[ome, gff_path] for ome, gff_path in gffs]
     #           )

        print('\tPfams', flush = True)
        ome2clus_genes, ome2hlg_genes = {}, {}
        for ome, genes2clus in hlgres:
            ome2hlg_genes[ome] = genes2clus
            taxon_hlg_genes += len(genes2clus)
        for ome, genes2clus in gcfres:
            ome2clus_genes[ome] = genes2clus
            taxonClusGenes += len(genes2clus)

        pfamres = pool.starmap(
            parsePfam,
            [[ome, f'{pfam_dir}{ome}_pfam.out', ome2clus_genes[ome], ome2hlg_genes[ome], top_hit] \
            for ome in omes]
            )

        ome2genes2pfam, ome2clus_genes2pfam, ome2hlgGenes2pfam, pfam2def = {}, {}, {}, {}
        for ome, tgene2pfam, tclusGenes2pfam, thlgGenes2pfam, tpfam2def in pfamres:
            if not tgene2pfam:
                eprint(f'\t\tWARNING: {ome} missing annotation', flush = True)
                continue
            ome2genes2pfam[ome] = {x[0]: dict(x[1]) for x in tgene2pfam}
            ome2clus_genes2pfam[ome] = {x[0]: x[1] for x in tclusGenes2pfam}
            ome2hlgGenes2pfam[ome] = {x[0]: x[1] for x in thlgGenes2pfam}
            pfam2def = {**pfam2def, **dict(tpfam2def)}

#        if as_dir:
 #           out_base = out_dir + taxon.lower() + '.as'
  #      else:
        out_base = out_dir + taxon.lower()

        print(f'\tEnrichment {out_base}*', flush = True)
        if pfam2go_path:
            pfam2go, go2onto, go2pfam = parsePfam2GO(pfam2go_path, pfam2def, go_sets)
            pfamScores, goScores, ref_scores = GetPopResultsGO(
                ome2genes2pfam, ome2hlgGenes2pfam, ome2clus_genes2pfam, taxonGenes, 
                taxon_hlg_genes, taxonClusGenes, pfam2go, go2pfam, pool
                )
            print(f'\tWriting', flush = True)
            WriteOutput(pfamScores[0], thresh, out_base + '_cvt.pfam.tsv')
            WriteOutput(pfamScores[1], thresh, out_base + '_cvh.pfam.tsv')
            WriteOutput(pfamScores[2], thresh, out_base + '_hvt.pfam.tsv')
#            WriteOutput(goScores[0], thresh, out_base + '_cvt.go.tsv', go2onto)
 #           WriteOutput(goScores[1], thresh, out_base + '_cvh.go.tsv', go2onto)
  #          WriteOutput(goScores[2], thresh, out_base + '_hvt.go.tsv', go2onto)
            if go2pfam:
                WriteOutput(ref_scores[0], thresh, out_base + '_cvt.ref.go.tsv', go2pfam = True)
                WriteOutput(ref_scores[1], thresh, out_base + '_cvh.ref.go.tsv', go2pfam = True)
                WriteOutput(ref_scores[2], thresh, out_base + '_hvt.ref.go.tsv', go2pfam = True)
        else:
            scoresDict = GetPopResults(ome2genes2pfam, ome2clus_genes2pfam, 
                                       taxonGenes, taxonClusGenes, pool)
            print(f'\tWriting', flush = True)

            WriteOutput(scoresDict, thresh, out_base + '.gcf.pfam.tsv', go2onto = None, sort = True)
            scoresDict = GetPopResults(ome2genes2pfam, ome2hlgGenes2pfam, 
                                       taxonGenes, taxonClusGenes, pool)

            WriteOutput(scoresDict, thresh, out_base + '.hlg.pfam.tsv', go2onto = None, sort = True)

def cli():

    parser = argparse.ArgumentParser(description = 'cloci2pfam enrichment')
    parser.add_argument('-a', '--alpha', help = 'Overall Bonferonni corrected alpha; DEFAULT: 0.05',
        type = float, default = 0.05)
    parser.add_argument('-d', '--mtdb', help = 'CLOCI mycotools db input', required = True)
    parser.add_argument('-r', '--rank', help = 'taxon rank', required = True)
    parser.add_argument('-c', '--cloci', help = 'CLOCI output dir', required = True)
#    parser.add_argument('-a', '--antismash_dir', help = 'antiSMASH dir with ome subfolders')
    parser.add_argument('-p', '--pfam_dir', help = 'Pfam tbl output dir', required = True)
    parser.add_argument('-m', '--minimum', help = 'Minimum omes to consider rank', 
                        type = int, default = 0)
    parser.add_argument('-p2g', '--pfam2go', help = 'pfam2go file for GO')
    parser.add_argument('-g', '--go_terms', help = '[-gf] Specific GO terms to mine')
    parser.add_argument('-gf', '--go_obo', help = 'GO ontology file')
    parser.add_argument('-t', '--top_hit', action = 'store_true', 
                        help = 'Only consider top Pfam hit')
    parser.add_argument('-o', '--out_dir')
    parser.add_argument('--cpu', type = int, default = 1)
    args = parser.parse_args()

    taxon_set = {'kingdom', 'phylum', 'subphylum', 'class', 'order', 'family', 'genus', 'species'}
    if args.rank.lower() not in taxon_set:
        raise ValueError('Taxon must be in ' + str(taxon_set))
#    elif (args.cloci and args.antismash_dir) or (not args.cloci and not args.antismash_dir):
 #       raise ValueError('-o2c OR -as must be inputted')

    db = mtdb(format_path(args.mtdb))
    #if args.cloci:
    ome_dir = format_path(args.cloci + 'ome/', force_dir = True)
   #     as_dir = None
  #  else:
 #       ome_dir = None
#        as_dir = format_path(args.antismash_dir)
    pfam_dir = format_path(args.pfam_dir)

    if not args.out_dir:
#        if os.path.basename(os.path.dirname(format_path(args.cloci)[:-1])) == 'ome':
 #           out_dir = mkOutput(format_path('./'), 'cloci2enrichment')
  #      else:
        out_dir = mkOutput(format_path(args.cloci), 'cloci2enrichment', 
                             suffix = None)
    else:
        out_dir = mkOutput(format_path(out_dir), 'cloci2enrichment')

    if args.pfam2go:
        pfam2go_path = format_path(args.pfam2go)
    else:
        pfam2go_path = None

    if args.go_terms:
        go_terms = args.go_terms.replace('"','').replace("'",'').replace(',',' ').split()
        go_file = format_path(args.go_obo)
    else:
        go_terms, go_file = [], None

    pool = mp.Pool(processes = args.cpu)
    main(db, args.alpha, args.rank.lower(), ome_dir, pfam_dir, out_dir, pfam2go_path, 
         go_file, go_terms, pool = pool, top_hit = args.top_hit, minimum = args.minimum)
    pool.close()
    pool.join()

    sys.exit(0)


if __name__ == '__main__':
    cli()
