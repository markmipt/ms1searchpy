import csv
import os
from sys import argv
import operator
from collections import defaultdict
from pyteomics import auxiliary as aux, fasta
import utils


settings = utils.settings(argv[-1])
prefix = settings.get('input', 'decoy prefix')
isdecoy = lambda x: x[0].startswith(prefix)
isdecoy_key = lambda x: x.startswith(prefix)
escore = lambda x: -x[1]
protsV = set()
path_to_valid_fasta = settings.get('input', 'valid proteins')
if path_to_valid_fasta:
    for prot in fasta.read(path_to_valid_fasta):
        protsV.add(prot[0].split(' ')[0])

prots_spc = defaultdict(float)
for filename in argv[1:-1]:
    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        csvreader.next()
        for row in csvreader:
            dbname = row[0]
            protein_score = float(row[1])
            prots_spc[dbname] += protein_score

sortedlist_spc = sorted(prots_spc.iteritems(), key=operator.itemgetter(1))[::-1]

checked = set()
for k, v in prots_spc.items():
    if k not in checked:
        if isdecoy_key(k):
            if prots_spc.get(k.replace(prefix, ''), -1e6) > v:
                del prots_spc[k]
                checked.add(k.replace(prefix, ''))
        else:
            if prots_spc.get(prefix + k, -1e6) > v:
                del prots_spc[k]
                checked.add(prefix + k)

filtered_prots = aux.filter(prots_spc.items(), fdr=0.01, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1,
                            full_output=True)

identified_proteins = 0
identified_proteins_valid = 0

for x in filtered_prots:
    if x[0] in protsV:
        identified_proteins_valid += 1
    identified_proteins += 1

for x in filtered_prots[:5]:
    print x[0], x[1]
print 'results:%s;number of identified proteins = %d;number of valid proteins = %d' % (
'union', identified_proteins, identified_proteins_valid)
with open(os.path.join(os.path.dirname(argv[1]), 'union_proteins.csv'), 'w') as output:
    output.write('dbname\tscore\n')
    for x in filtered_prots:
        output.write('\t'.join((x[0], str(x[1]))) + '\n')
