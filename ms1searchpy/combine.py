from .main import final_iteration, filter_results
import pandas as pd
from collections import defaultdict
import argparse
import logging

logger = logging.getLogger(__name__)

def run():
    parser = argparse.ArgumentParser(
        description='Combine DirectMS1 search results',
        epilog='''

    Example usage
    -------------
    $ ms1combine.py file1_PFMs_ML.tsv ... filen_PFMs_ML.tsv
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', nargs='+', help='input tsv PFMs_ML files')
    parser.add_argument('-out', help='prefix output file names', default='combined')
    parser.add_argument('-prots_full', help='path to any of *_proteins_full.tsv file. By default this file will be searched in the folder with PFMs_ML files', default='')
    parser.add_argument('-fdr', help='protein fdr filter in %%', default=1.0, type=float)
    parser.add_argument('-prefix', help='decoy prefix', default='DECOY_')
    parser.add_argument('-nproc', help='number of processes', default=1, type=int)
    args = vars(parser.parse_args())
    logging.basicConfig(format='%(levelname)9s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=logging.INFO)


    d_tmp = dict()

    df1 = None
    for idx, filen in enumerate(args['file']):
        df3 = pd.read_csv(filen, sep='\t')
        df3['ids'] = df3['ids'].apply(lambda x: '%d:%s' % (idx, str(x)))
        df3['fidx'] = idx

        df3 = df3[df3['qpreds'] <= 10]


        qval_ok = 0
        for qval_cur in range(10):
            if qval_cur != 10:
                df1ut = df3[df3['qpreds'] == qval_cur]
                decoy_ratio = df1ut['decoy'].sum() / len(df1ut)
                d_tmp[(idx, qval_cur)] = decoy_ratio
                print(filen, qval_cur, decoy_ratio)

        if df1 is None:
            df1 = df3
            if args['prots_full']:
                df2 = pd.read_csv(args['prots_full'], sep='\t')
            else:
                try:
                    df2 = pd.read_csv(filen.replace('_PFMs_ML.tsv', '_proteins_full.tsv'), sep='\t')
                except:
                    logging.critical('Proteins_full file is missing!')
                    break

        else:
            df1 = pd.concat([df1, df3], ignore_index=True)

    d_tmp = [z[0] for z in sorted(d_tmp.items(), key=lambda x: x[1])]
    qdict = {}
    for idx, val in enumerate(d_tmp):
        qdict[val] = int(idx / len(args['file']))
    df1['qpreds'] = df1.apply(lambda x: qdict[(x['fidx'], x['qpreds'])], axis=1)


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

    # mass_diff = resdict['qpreds']
    # rt_diff = resdict['qpreds']

    base_out_name = args['out']

    e_ind = resdict['qpreds'] <= 9
    resdict = filter_results(resdict, e_ind)

    mass_diff = resdict['qpreds']
    rt_diff = resdict['qpreds']
    
    final_iteration(resdict, mass_diff, rt_diff, pept_prot, protsN, base_out_name, prefix, isdecoy, isdecoy_key, escore, fdr, args['nproc'])


if __name__ == '__main__':
    run()
