#! /usr/bin/env python3

# hit percentage > 50%, evalue < 0.01
# index 15, 14, 5, and 6
# total genes for taxon of interest, total each pfam for taxon of interest
    # acquire

# NEED to make enrichment analysis consider all Pfam hits for a gene if it isn't already

from scipy.stats import hypergeom
from mycotools.lib.biotools import fa2dict, gff3Comps, gff2list
from mycotools.lib.kontools import formatPath
from mycotools.lib.dbtools import mtdb
from collections import Counter, defaultdict
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

def parsePfam2GO(pfam2go_path, pfam2def):
    pfam2go, go2onto = defaultdict(list), {}
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
                pfam2go[definition].append(go)
                go2onto[go] = onto
            except TypeError:
                continue
    return pfam2go, go2onto
            

def countGenes(ome, proteome_path):
    return ome, len(fa2dict(proteome_path))

def clusterGenes(ome, info_path):
    genes = {}
    with open(info_path, 'r') as raw:
        for line in raw:
            d = line.split('\t')
            clus, genesT = d[0], d[2].split(',')
            genes = {**genes, **{x: clus for x in genesT}}
    return ome, tuple(genes.items())

def parsePfam(ome, pfam_res, clusGenes):
    gene2pfam = defaultdict(dict)
    clusGenes = dict(clusGenes)
    pfam2def = {}
    with open(pfam_res, 'r') as raw:
        for line in raw:
            if line.startswith('#'):
                continue
            data = line.rstrip().split()
            g, p, q, b = data[0], data[3], data[2], float(data[5])
            
            if 'F-box' in q:
                if g in clusGenes:
                    print(g, flush=True)

            # this is where a post-hoc threshold would be applied
            if q not in gene2pfam[g]: # in case there is a duplicate
            # for the same gene, add in the better score 
                gene2pfam[g][q] = b
            else:
                if b > gene2pfam[g][q]:
                    gene2pfam[g][q] = b
            pfam2def[p[:p.find('.')]] = q
    gene2pfam, clusGenes2pfam = dict(gene2pfam), {}
    for gene, clus in clusGenes.items():
        try:
            clusGenes2pfam[gene] = tuple([gene2pfam[gene], clus])
        except KeyError: # no pfam acc found, should be rarer
            pass

    return ome, tuple([tuple([gene, tuple(pfams.items())]) for gene, pfams in gene2pfam.items()]), \
        tuple(clusGenes2pfam.items()), tuple(pfam2def.items())
#        tuple([tuple([gene, pfams.items()]) for gene, pfams in clusGenes2pfam.items()])


def calcP(pfam, totCount, clusCount, totalGenes, totalClusGenes, pfamInClus):
    pval = hypergeom.sf(
        clusCount-1, totalGenes, totCount, totalClusGenes
        )
    return pfam, pval, clusCount, pfamInClus


def GetPopResultsGO(
    omes2genes2pfams, omes2clus2pfam, totalGenes, totalClusGenes, 
    pfam2go, pool
    ):
    totalPfams, clusPfams, pfam2clusters = Counter(), Counter(), defaultdict(list)
    totalGOs, clusGOs, go2clusters = Counter(), Counter(), defaultdict(list)
    for ome, genes2pfam in omes2genes2pfams.items():
        clusGenes2pfam = omes2clus2pfam[ome]
        tPfams, tClusPfams = [], []
        tGOs, tClusGOs = [], []
        for pfams in list(genes2pfam.values()):
            tPfams.extend(pfams)
            for pfam in pfams:
                tGOs.extend(pfam2go[pfam])
        for pfams, clus in list(clusGenes2pfam.values()):
#            if any('F-box' in x for x in pfams) and any('WD40' in x for x in pfams):
 #               print(clus, flush = True)
            for pfam in pfams:
                pfam2clusters[pfam].append(ome + '_' + clus)
                for go in pfam2go[pfam]:
                    go2clusters[go].append(ome + '_' + clus)
                    tClusGOs.append(go)
                tClusPfams.append(pfam)
        totalPfams += Counter(tPfams)
        clusPfams += Counter(tClusPfams)
        totalGOs += Counter(tGOs)
        clusGOs += Counter(tClusGOs)

    cmds = []
    for pfam, clusCount in clusPfams.items():
        totCount = totalPfams[pfam]
        pfamInClus = len(set(pfam2clusters[pfam]))
        cmds.append([pfam, totCount, clusCount, totalGenes, totalClusGenes, pfamInClus])
    p_pval = pool.starmap(calcP, cmds)

    cmds =  []
    for go, clusCount in clusGOs.items():
        totCount = totalGOs[go]
        goInClus = len(set(go2clusters[go]))
        cmds.append([go, totCount, clusCount, totalGenes, totalClusGenes, goInClus])
    g_pval = pool.starmap(calcP, cmds)

    return p_pval, g_pval


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
                print(pfam, flush = True)
                if pfam == 'DUF3435':
                    print(ome, clus, flush = True)
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

def WriteOutput(scoresDict, thresh, out_file, go2onto = None):            
# out_dir + taxon.lower() + '.scores.tsv'
    alpha = thresh/len(scoresDict) # bonferroni correction

    sig, insig = [], []
    with open(out_file, 'w') as out:
        if go2onto:
            out.write(
                '#go\ttotal\tclusters\tp\tsignificance (alpha =' + \
                str(alpha) + ')\tontology\n'
                )
        else:
            out.write('#pfam\ttotal\tclusters\tp\tsignificance (alpha=' + str(alpha) + ')\n')
        for q, val, tot, clusTot in scoresDict:
            if val <= alpha: # significant
                sig.append([q, tot, clusTot, val, 1])
                if go2onto:
                    sig[-1].append(go2onto[q])
            else: # not significant
                insig.append([q, tot, clusTot, val, 0])
                if go2onto:
                    insig[-1].append(go2onto[q])
        sig = sorted(sig, key = lambda x: (x[-1], x[-3]), reverse = True)
        out.write('\n'.join(['\t'.join([str(x) for x in y]) for y in sig]) + '\n')
        out.write('\n'.join(['\t'.join([str(x) for x in y]) for y in insig]))



def main(
    db, thresh, rank, og2clus_dir, pfam_dir, out_dir,
    pfam2go_path, as_dir = None, pool = mp.Pool(processes = 1)
    ):

    db = db.set_index()
    ranks = defaultdict(list)
    for ome, row in db.items():
        taxon = {r.lower(): t for r, t in row['taxonomy'].items()}[rank]
        ranks[taxon].append(ome)

    for taxon, omes in ranks.items():
        todel = []
        for i, ome in enumerate(omes):
            if og2clus_dir:
                if not os.path.isfile(og2clus_dir + 'ome/' + ome + '/info.out'):
                    todel.append(i)
            elif not os.path.isfile(as_dir + ome + '/' + ome + '.gff3'):
                todel.append(i)
        for i in reversed(todel):
            del omes[i]

    ranks = {k: v for k, v in sorted(ranks.items(), key = lambda x: len(x[1]))}
    for taxon, omes in ranks.items():
        print('\n' + taxon + ': ' + str(len(omes)) + ' omes', flush = True)
        taxonGenes, taxonClusGenes = 0, 0
        countres = pool.starmap(
            countGenes, 
            [[ome, formatPath('$MYCOFAA/' + ome + '.aa.fa')] for ome in omes]
            )
        for ome, tot in countres:
            taxonGenes += tot

        if not as_dir:
            clusres = pool.starmap(
                clusterGenes,
                [[ome, og2clus_dir + 'ome/' + ome + '/info.out'] for ome in omes]
                )
        elif as_dir:
            gffs = collectGffs(omes, as_dir)
            clusres = pool.starmap(
                parseGff,
                [[ome, gff_path] for ome, gff_path in gffs]
                )

        ome2clusGenes = {}
        for ome, genes2clus in clusres:
            ome2clusGenes[ome] = genes2clus
            taxonClusGenes += len(genes2clus)

        pfamres = pool.starmap(
            parsePfam,
            [[ome, pfam_dir + ome + '.out', ome2clusGenes[ome]] for ome in omes]
            )

        ome2genes2pfam, ome2clusGenes2pfam, pfam2def = {}, {}, {}
        for ome, tgene2pfam, tclusGenes2pfam, tpfam2def in pfamres:
            ome2genes2pfam[ome] = {x[0]: dict(x[1]) for x in tgene2pfam}
            ome2clusGenes2pfam[ome] = {x[0]: x[1] for x in tclusGenes2pfam}
            pfam2def = {**pfam2def, **dict(tpfam2def)}

        if as_dir:
            out_base = out_dir + taxon.lower() + '.as'
        else:
            out_base = out_dir + taxon.lower() + '.o2c'

        if pfam2go_path:
            pfam2go, go2onto = parsePfam2GO(pfam2go_path, pfam2def)
            pfamScores, goScores = GetPopResultsGO(
                ome2genes2pfam, ome2clusGenes2pfam, taxonGenes, taxonClusGenes, 
                pfam2go, pool
                )
            WriteOutput(pfamScores, thresh, out_base + '.pfam.tsv')
            WriteOutput(goScores, thresh, out_base + '.go.tsv', go2onto)
        else:
            scoresDict = GetPopResults(ome2genes2pfam, ome2clusGenes2pfam, taxonGenes, taxonClusGenes, pool)
            WriteOutput(scoresDict, thresh, out_base + '.pfam.tsv', go2onto = None)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'og2clus2pfam enrichment')
    parser.add_argument('-a', '--alpha', help = 'Pre-Bonferonni corrected alpha. DEFAULT: 0.05',
        type = float, default = 0.05)
    parser.add_argument('-d', '--db', help = 'og2clus mycotools db input', required = True)
    parser.add_argument('-r', '--rank', help = 'taxon rank', required = True)
    parser.add_argument('-o2c', '--og2clus_dir', help = 'og2clus output dir')
    parser.add_argument('-as', '--antismash_dir', help = 'antiSMASH dir with ome subfolders')
    parser.add_argument('-p', '--pfam_dir', help = 'pfam tbl output dir', required = True)
    parser.add_argument('-p2g', '--pfam2go', help = 'pfam2go file for GO')
    parser.add_argument('-o', '--out_dir')
    parser.add_argument('-c', '--cpu', type = int, default = 1)
    args = parser.parse_args()

    taxon_set = {'phylum', 'subphylum', 'class', 'order', 'family', 'genus', 'species'}
    if args.rank.lower() not in taxon_set:
        raise ValueError('Taxon must be in ' + str(taxon_set))
    elif (args.og2clus_dir and args.antismash_dir) or (not args.og2clus_dir and not args.antismash_dir):
        raise ValueError('-o2c OR -as must be inputted')

    db = mtdb(formatPath(args.db))
    if args.og2clus_dir:
        og2clus_dir = formatPath(args.og2clus_dir)
        as_dir = None
    else:
        og2clus_dir = None
        as_dir = formatPath(args.antismash_dir)
    pfam_dir = formatPath(args.pfam_dir)

    if not args.out_dir:
        out_dir = os.getcwd() +'/enrichment/'
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
    else:
        out_dir = os.getcwd() + '/'

    if args.pfam2go:
        pfam2go_path = formatPath(args.pfam2go)

    pool = mp.Pool(processes = args.cpu)
    main(db, args.alpha, args.rank.lower(), og2clus_dir, pfam_dir, out_dir, pfam2go_path, as_dir, pool)
    pool.close()
    pool.join()

    sys.exit(0)
