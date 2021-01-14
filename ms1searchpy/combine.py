from .main import final_iteration
import pandas as pd
from collections import defaultdict
import argparse


def run():
    parser = argparse.ArgumentParser(
        description='Combine DirectMS1 search results',
        epilog='''

    Example usage
    -------------
    $ search.py file1_PFMs_ML.csv ... filen_PFMs_ML.csv
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', nargs='+', help='input csv PFMs_ML files')
    parser.add_argument('-out', help='prefix output file names', default='combined')
    parser.add_argument('-prots_full', help='path to any of *_proteins_full.csv file. By default this file will be searched in the folder with PFMs_ML files', default='')
    parser.add_argument('-fdr', help='protein fdr filter in %%', default=1.0, type=float)
    parser.add_argument('-prefix', help='decoy prefix', default='DECOY_')
    parser.add_argument('-nproc',   help='number of processes', default=1, type=int)
    args = vars(parser.parse_args())

    df1 = False
    for idx, filen in enumerate(args['file']):
        df3 = pd.read_csv(filen, sep='\t')
        df3['ids'] = df3['ids'].apply(lambda x: '%d:%s' % (idx, str(x)))
        if df1 is False:
            df1 = df3
            if args['prots_full']:
                df2 = pd.read_csv(args['prots_full'], sep='\t')
            else:
                try:
                    df2 = pd.read_csv(filen.replace('_PFMs_ML.csv', '_proteins_full.csv'), sep='\t')
                except:
                    print('Proteins_full file is missing!')
                    break

        else:
            df1 = df1.append(df3, ignore_index=True)

    pept_prot = defaultdict(set)
    for seq, prots in df1[['seqs', 'proteins']].values:
        for dbname in prots.split(';'):
            pept_prot[seq].add(dbname)

    protsN = dict()
    for dbname, theorpept in df2[['dbname', 'theoretical peptides']].values:
        protsN[dbname] = theorpept


    prefix = args['prefix']
    isdecoy = lambda x: x[0].startswith(prefix)
    isdecoy_key = lambda x: x.startswith(prefix)
    escore = lambda x: -x[1]
    fdr = float(args['fdr']) / 100

    resdict = dict()

    resdict['qpreds'] = df1['qpreds'].values
    resdict['preds'] = df1['preds'].values
    resdict['seqs'] = df1['peptide'].values
    resdict['ids'] = df1['ids'].values
    resdict['iorig'] = df1['iorig'].values

    mass_diff = resdict['qpreds']
    rt_diff = resdict['qpreds']

    base_out_name = args['out']
    final_iteration(resdict, mass_diff, rt_diff, pept_prot, protsN, base_out_name, prefix, isdecoy, isdecoy_key, escore, fdr, args['nproc'])


if __name__ == '__main__':
    run()
