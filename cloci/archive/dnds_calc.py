import os
import warnings
import numpy as np
import multiprocessing as mp
from Bio import AlignIO, SeqIO, BiopythonWarning, codonalign
from collections import defaultdict
from mycotools.gff2seq import ntmain
from mycotools.lib.kontools import multisub
from mycotools.lib.biotools import gff2list, fa2dict, dict2fa
from orthocluster.orthocluster.lib import input_parsing
warnings.simplefilter('ignore', BiopythonWarning)

def dndsGeneGrab(
    gff_path, assembly_path, proteome_path,
    ome_sig_clus, gene2hg, clusplusminus
    ):
    
    gff_list, prot_dict = gff2list(gff_path), fa2dict(proteome_path)
    assem_dict, clus_out = fa2dict(assembly_path), {}
    cds_dict, cds_dict2 = input_parsing.compileCDS2(
        gff_list, os.path.basename(gff_path).replace('.gff3','')
        )
    
    for scaf in cds_dict:
        for i0, seq0 in enumerate(cds_dict[scaf]):
            if seq0 in ome_sig_clus: # if the og is part of a significant seed
            # locus
                try:
                    og0 = gene2hg[seq0]
                except KeyError:
                    continue
                for hit in ome_sig_clus[seq0]:
                    hgx = tuple(sorted(list(hit)))
                    if hgx not in clus_out:
                        clus_out[hgx] = [{og: [{}, {}] for og in hgx}]
                    clus_out[hgx][-1][og0][0][seq0] = prot_dict[seq0]
                    clus_out[hgx][-1][og0][1] = {
                        **clus_out[hgx][-1][og0][1], 
                        **ntmain(cds_dict2[seq0], assem_dict)
                        }
                    
                    start, end = None, None
                    if i0 < clusplusminus: # is i0 - clusplusminus < 0 ?
                        locus = cds_dict[scaf][:i0+clusplusminus+1] # then gather
                        # all the beginning
                    else:
                        locus = cds_dict[scaf][i0-clusplusminus:i0+clusplusminus+1]
                        # instead get the +/- and adjust for python
                    for i1, seq1 in enumerate(locus): # for each index and sequence
                    # in the locus
                        try:
                            og = gene2hg[seq1]
                        except KeyError:
                            continue
                        if og in hit: # if the og is
                        # in the sig clus
                            clus_out[hgx][-1][og][0][seq1] = prot_dict[seq1]
                            clus_out[hgx][-1][og][1] = {
                                **clus_out[hgx][-1][og][1],
                                **ntmain(cds_dict2[seq1], assem_dict)
                                }
    return clus_out



def dnds_preparation(db, omes4hgs, hgx2omes, gene2hg, plusminus, hgx_dir, i2ome):
    
    gene_checks, hgx_check = {}, {}
    
    for ome in omes4hgs:
        gff_path = db[i2ome[ome]]['gff3'] 
        proteome_path = db[i2ome[ome]]['faa']
        assembly_path = db[i2ome[ome]]['fna']
        res = dndsGeneGrab(
            gff_path, assembly_path, proteome_path,
            omes4hgs[ome], gene2hg, plusminus
            )
        for hgx in res:
            if hgx not in hgx_check:
                hgx_check[hgx] = set()
            if hgx not in gene_checks:
                gene_checks[hgx] = {og: [{}, {}] for og in hgx}
            for loc in res[hgx]:
                for og in loc:
                    gene_checks[hgx][og][0] = {
                        **gene_checks[hgx][og][0], **loc[og][0]
                        }
                    gene_checks[hgx][og][1] = {
                        **gene_checks[hgx][og][1], **loc[og][1]
                        }
            hgx_check[hgx].add(ome)
            if hgx_check[hgx] == set(hgx2omes[hgx]):
                for og in hgx: 
                    out_base = hgx_dir + '-'.join([str(x) for x in hgx]) + \
                        '.' + str(og)
                    with open(out_base + '.aa.fa', 'w') as out:
                        out.write(dict2fa(gene_checks[hgx][og][0]))
                    with open(out_base + '.nt.fa', 'w') as out:
                        out.write(dict2fa(gene_checks[hgx][og][1]))
                del gene_checks[hgx]
    
    if gene_checks:
        print('\t\tERROR: discrepancies with previous run', flush = True)
        for hgx in gene_checks:
            print('\t\t\t' + str(hgx), flush = True)
            for og in hgx: 
                out_base = hgx_dir + '-'.join([str(x) for x in hgx]) + \
                    '.' + str(og)
                with open(out_base + '.aa.fa', 'w') as out:
                    out.write(gene_checks[hgx][og][0]) 
                with open(out_base + '.nt.fa', 'w') as out:
                    out.write(gene_checks[hgx][og][1])



def calc_dnds(mafft):
    
    og_catch = re.search(r'(^[^\.]+)\.(\d+)\.mafft$', os.path.basename(mafft))
    hgx = tuple([int(x) for x in og_catch[1].split('-')])
    og = int(og_catch[2])
    
    err = os.path.basename(mafft).replace('.mafft','')
    nt_path = mafft.replace('.mafft', '.nt.fa')
    try:
        prot_align = AlignIO.read(mafft, 'fasta')
    except ValueError: 
        print('\t\t' + err + ' empty alignment?', flush = True)
        return None
    
    nt_fa = SeqIO.parse(nt_path, 'fasta')
    
    try:
        aln = codonalign.build(prot_align, nt_fa) #, ids)
    except KeyError: 
        print('\t\t' + err + ' ambiguous AAs', flush = True)
        return
    except RuntimeError:
        print('\t\t' + err + ' BioPython error', flush = True)
        return
    except IndexError: 
        print('\t\t' + err + ' unknown error', flush = True)
        return
    try:
        dn_matr, ds_matr = aln.get_dn_ds_matrix()
        np_dn, np_ds = np.array(dn_matr), np.array(ds_matr)
        dnds_matrNulls = np.divide(np_dn, np_ds)
        dnds_matr = np.ma.array(dnds_matrNulls, mask=np.isnan(dnds_matrNulls))
        dnds_matr = np.abs(dnds_matr)
        dnds = np.mean(dnds_matr)
        return hgx, og, dnds
    except ZeroDivisionError:
        print('\t\t' + err + ' unknown error')
    except KeyError: 
        print('\t\t' + err + ' ambiguous NTs')
    return None


def parse_dnds(mpRes):

    hgxdNdS_dict = {}
    for res in mpRes:
        if res:
            hgx, og, dnds = res[0], res[1], res[2]
            try:
                float(dnds)
            except ValueError:
                continue
            if hgx not in hgxdNdS_dict:
                hgxdNdS_dict[hgx] = {
                    'selection_total': 0,
                    'og_total': 0,
                    'hgs': {}
                    }

            if dnds < 1: # if purifying selection
                hgxdNdS_dict[hgx]['selection_total'] += dnds
            else: # if positive selection use reciprocal
                hgxdNdS_dict[hgx]['selection_total'] += 1/dnds
            hgxdNdS_dict[hgx]['og_total'] += 1
            hgxdNdS_dict[hgx]['hgs'][og] = dnds
    for hgx in hgxdNdS_dict:
        selection_coefficient = hgxdNdS_dict[hgx]['selection_total']/hgxdNdS_dict[hgx]['og_total']
        mean_dnds = sum(hgxdNdS_dict[hgx]['hgs'].values())/len(hgxdNdS_dict[hgx]['hgs'])
        hgxdNdS_dict[hgx] = [selection_coefficient, mean_dnds, hgxdNdS_dict[hgx]['hgs']]

    return hgxdNdS_dict

def dnds_main(db, ome2i, gcfs, hgx_dir, i2ome, hx2omes, gene2hg, plusminus):
    print('\nIX. Quantifying GCF dn/ds', flush = True)
    omes4hgs, hgx_files = {x: defaultdict(list) for x in range(len(ome2i))}, []
    for i, gcf in enumerate(gcfs):
        for hgx in gcf:
            hgx = hgx2i[x]
            for gene in hgx2loc[hgx]:
                omeI = ome2i[gene[:gene.find('_')]]
                if omeI in set(gcf[i][hgx]):
                    omes4hgs[omeI][gene].append(gcf_hgxs[i])

    hgx_files = []
    for hgx in gcf_hgxs:
        for og in hgx:
            out_base = hgx_dir + '-'.join([str(x) for x in hgx]) + '.' + str(og)
            hgx_files.append(out_base)

    maffts = [file + '.mafft' for file in hgx_files]
    if not all(os.path.isfile(x + '.aa.fa') for x in hgx_files):
        print("\tCompiling HGx kernels' genes", flush = True)
        dnds_preparation(db, omes4hgs, hgx2omes, gene2hg, plusminus, hgx_dir, i2ome)
    if not all(os.path.isfile(x) for x in maffts):
        print('\tAligning proteins', flush = True)
        maffts = [x for x in maffts if not os.path.isfile(x)]
        mafft_cmds = [
            [['mafft', '--auto', '--thread', '2', x.replace('.mafft', '.aa.fa')], x,
            x.replace('.mafft', '.aa.fa')] \
            for x in maffts
            ]
        mafft_res = multisub(mafft_cmds, processes = cpus, stdout = True, rm = True)
        maffts = [x['output'] for x in mafft_res]

    print('\tCalculating dn/ds', flush = True)
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        dnds_res = pool.map(calc_dnds, maffts)
        pool.close()
        pool.join()

    hgx2dnds = parse_dnds(dnds_res)
