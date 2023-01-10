import pickle
from mycotools.lib.kontools import read_json


# NEED clan2gcf in hgx2gcfs
# NEED to filter irrelevant omes
# NEED to open and obtain cluster id from ome/info.out
# NEED to build network
	# NEED to color code

def import_data(gcf_file, clan2loci_file, adj_file, clans):
    clanLoci = {int(k): tuple(v) for k,v in sorted(
        read_json(clan2loci_file).items(),
        key = lambda x: len(x[1]), reverse = True)}

    index, loci2grab = 0, []
    for i0, clanI in enumerate(sorted(clans)):
        if i0 != 0:
            to_acct = range(clans[i0 - 1] + 1, clanI)
            for i1 in to_acct:
                index += len(clanLoci[i1])
        else:
            to_acct = range(0, clanI)
            for i1 in to_acct:
                index += len(clanLoci[i1])

        for loc in clanLoci[clanI]:
            loci2grab.append(index)
            index += 1

    loci2grab, adj = set(loci2grab), {}
    with open(adj_file, 'r') as raw:
        for line in raw:
            l0, l1, gid = line.rstrip().split()
            if int(l0) in loci2grab and int(l1) in loci2grab:
                adj[(l0, l1)] = float(gid)
        


    with open(gcf_file, 'rb') as raw:
         gcf_data = pickle.load(raw)
