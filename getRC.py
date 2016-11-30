from os import path
import argparse
import cPickle
from pyteomics import achrom, auxiliary as aux

parser = argparse.ArgumentParser(
    description='Create file with Retention Coefficients trained on input peptides.',
    epilog='''

Example usage
-------------
  $ getRC.py input.txt -o output.pickle -d tab
-------------

input file must contains lines with peptide sequence and Retention times:
PEPTIDE 10.5
PEPTIAS 30.1
KRPEEQ  59.1
...
''',
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('file', help='input txt file with peptides and RT experimental')
parser.add_argument('-o', '--out', help='path to output file')
parser.add_argument('-d', '--delimiter', default='tab', help='delimiter, can be: tab, comma, space')
args = parser.parse_args()

infile = args.file
outfile = args.out
if not outfile:
    outfile = path.splitext(infile)[0] + '.pickle'
delmap = {
    'tab': '\t',
    'comma': ',',
    'space': ' '
}
delimiter = delmap[args.delimiter]

seqs = []
RTexp = []
with open(infile) as inp:
    for line in inp:
        if line.strip():
            seq, RT = line.strip().split(delimiter)
            seqs.append(seq)
            RTexp.append(float(RT))
RC_def = achrom.RCs_gilar_rp
aa_labels = set(RC_def['aa'].keys())
xdict = {}
for key, val in RC_def['aa'].items():
    xdict[key] = [val, None]
RC_dict = achrom.get_RCs_vary_lcp(seqs, RTexp, labels=aa_labels)

RTpred = [achrom.calculate_RT(seq, RC_dict) for seq in seqs]

_, _, Rcoef, sigma = aux.linear_regression(RTpred, RTexp)

for key, val in RC_dict['aa'].items():
    try:
        xdict[key][1] = val
    except:
        xdict[key] = [None, val]
a, b, _, _ = aux.linear_regression([x[0] for x in xdict.values() if all(v is not None for v in x)],
                                   [x[1] for x in xdict.values() if all(v is not None for v in x)])
for key, x in xdict.items():
    if x[1] is None:
        x[1] = x[0] * a + b
    RC_dict['aa'][key] = x[1]
if 'C' not in RC_dict['aa']:
    RC_dict['aa']['C'] = RC_dict['aa']['C*']
cPickle.dump(RC_dict, open(outfile, 'w'))
print 'Model results: R^2 = %.2f, stdev = %.2f' % (Rcoef ** 2, sigma)
print 'Retention coefficients were saved.'
