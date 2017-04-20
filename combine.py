import csv
import os
import operator
from collections import defaultdict
from pyteomics import auxiliary as aux
import argparse

parser = argparse.ArgumentParser(
    description='Search proteins using LC-MS spectra',
    epilog='''

Example usage
-------------
$ search.py results_1_proteins_full.csv ... results_n_proteins_full.csv
-------------
''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('files', help='input .csv proteins_full files', nargs='+')
parser.add_argument('-prefix', help='decoy prefix', default='DECOY_')
args = vars(parser.parse_args())
prefix = args['prefix']

isdecoy = lambda x: x[0].startswith(prefix)
isdecoy_key = lambda x: x.startswith(prefix)
escore = lambda x: -x[1]

prots_spc = defaultdict(float)
for filename in args['files']:
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
    identified_proteins += 1

for x in filtered_prots[:5]:
    print x[0], x[1]
print 'results:%s;number of identified proteins = %d' % ('union', identified_proteins)
with open(os.path.join(os.path.dirname(args['files'][0]), 'union_proteins.csv'), 'w') as output:
    output.write('dbname\tscore\n')
    for x in filtered_prots:
        output.write('\t'.join((x[0], str(x[1]))) + '\n')
