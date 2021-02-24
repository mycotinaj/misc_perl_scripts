##########################
# Written by Kevin Amses
# amsesk@umich.edu
##########################

import sys
import os
from ete3 import NCBITaxa
from collections import OrderedDict

def best_blast_hit (tabular, bitcol=11):
    best = OrderedDict()
    for line in open(tabular).readlines():
        spl = line.split("\t")
        spl = map(str.strip,spl)
        label = spl[0]
        bit = float(spl[bitcol])
        if label in best.keys():
            if float(best[label][bitcol]) < bit:
                best[label] = spl
        else:
            best[label] = spl
    return best

if sys.argv[1] in ['-h','--help']:
    print '''
ncbi_taxrpt.py
USAGE: python ncbi_taxrpt.py [blast_output]\n
    REQUIREMENT IN ORIGINAL BLAST CALL:
    \t-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"
    '''
    sys.exit(0)
ncbi = NCBITaxa()
ids = OrderedDict()
rank_of_int = ('superkingdom','kingdom','phylum','class','order','family','genus','species')
#print '\t'.join(['query','hit','evalue'])+'\t'+'\t'.join(rank_of_int)

for key,best in best_blast_hit(sys.argv[1]).iteritems():
    name = best[0]
    hit = best[1]
    tid = best[12]
    evalue = best[10]
    if name not in ids.keys():
        ids[name] = {'hit': '',
                    'lineage': None
                    }
    ids[name]['hit'] = hit
    ids[name]['evalue'] = evalue

    if tid == "N/A":
        ids[name]['lineage'] = {"superkingdom":"No_associated_taxid"}
    else:
        try:
            if len(tid.split(';')) > 1:
                for i in tid.split(';'):
                    i = int(i)
                    if name in ids.keys():
                        lineage = ncbi.get_taxid_translator(ncbi.get_lineage(i))
                        ids[name]['lineage'].update(lineage)
                    else:
                        lineage = ncbi.get_taxid_translator(ncbi.get_lineage(i))
                        ids[name]['lineage'] = lineage
            else:
                tid = int(tid)
                lineage = ncbi.get_taxid_translator(ncbi.get_lineage(tid))
                ids[name]['lineage'] = lineage
        except:
            ids[name]['lineage'] = {"superkingdom":"Not_in_taxdb"}

print "{}\t{}".format('\t'.join(['query','hit','evalue']), '\t'.join(rank_of_int))
for key,p in ids.iteritems():
    ranks = ncbi.get_rank(p['lineage'].keys())
    ids_at_rank = [[i,r] for [i,r] in ranks.iteritems() if r in rank_of_int]
    for lvl in ids_at_rank:
        lvl[0] = ncbi.get_taxid_translator([lvl[0]])
        lvl[0] = lvl[0][lvl[0].keys()[0]]
    ordered = []
    for r in rank_of_int:
        found = False
        for lvl in ids_at_rank:
            if lvl[1] == r:
                ordered.append(lvl[0])
                found = True
        if not found:
            ordered.append("NA")
    print '\t'.join([key,p['hit'],p['evalue']])+'\t'+'\t'.join(ordered)
