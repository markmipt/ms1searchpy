from __future__ import division
import argparse
import pkg_resources
import pandas as pd
import ast
import subprocess
import numpy as np
from scipy.stats import binom, ttest_ind
from collections import defaultdict
import itertools
from pyteomics import auxiliary as aux

def calc_sf_all(v, n, p):
    sf_values = -np.log10(binom.sf(v-1, n, p))
    sf_values[np.isinf(sf_values)] = max(sf_values[~np.isinf(sf_values)]) * 2
    return sf_values

def run():
    parser = argparse.ArgumentParser(
        description='run DirectMS1quant for ms1searchpy results',
        epilog='''

    Example usage
    -------------
    $ directms1quant -S1 sample1_1_proteins_full.tsv sample1_n_proteins_full.tsv -S2 sample2_1_proteins_full.tsv sample2_n_proteins_full.tsv
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-S1', nargs='+', help='input files for S1 sample', required=True)
    parser.add_argument('-S2', nargs='+', help='input files for S2 sample', required=True)
    # parser.add_argument('-allowed_prots', help='path to allowed prots', default='')
    parser.add_argument('-out', help='name of DirectMS1quant output file', default='directms1quant_out')
    parser.add_argument('-use_filt', help='use_filt', action='store_true')

    # parser.add_argument('-norm', help='normalization method. Can be average, median, GMM or None', default='None')
    # parser.add_argument('-impute_threshold', help='impute_threshold for missing values fraction', default='0.75')
    # parser.add_argument('-min_samples', help='minimum number of samples for peptide usage', default='3')
    # parser.add_argument('-new', help='new algo', action='store_true')
    args = vars(parser.parse_args())

    replace_label = '_proteins_full.tsv'

    df_final = False
    protsN = False

    all_s_lbls = {}

    allowed_prots = set()
    # allowed_peptides = set()

    for i in range(1, 3, 1):
        sample_num = 'S%d' % (i, )
        if args[sample_num]:
            for z in args[sample_num]:
                if args['use_filt']:
                    df0 = pd.read_table(z.replace('_proteins_full.tsv', '_proteins.tsv'))
                    allowed_prots.update(df0['dbname'])
                    allowed_prots.update(['DECOY_' + z for z in df0['dbname'].values])
                else:
                    df0 = pd.read_table(z)
                    allowed_prots.update(df0['dbname'])
                if not protsN:
                    df0 = pd.read_table(z)
                    protsN = df0.set_index('dbname')['theoretical peptides'].to_dict()


    # else:
    #     for prot in open(args['allowed_prots'], 'r'):
    #         allowed_prots.add(prot.strip())


    for i in range(1, 3, 1):
        sample_num = 'S%d' % (i, )
        if args.get(sample_num, 0):
            all_s_lbls[sample_num] = []
            for z in args[sample_num]:
                label = z.replace(replace_label, '')
                all_s_lbls[sample_num].append(label)
                df1 = pd.read_table(z)
                df3 = pd.read_table(z.replace(replace_label, '_PFMs.tsv'))
                print(z)
                df3 = df3[df3['proteins'].apply(lambda x: any(z in allowed_prots for z in x.split(';')))]
                df3['proteins'] = df3['proteins'].apply(lambda x: ';'.join([z for z in x.split(';') if z in allowed_prots]))

                df3['sequence'] = df3['sequence'] + df3['charge'].astype(int).astype(str)

                df3 = df3.sort_values(by='Intensity', ascending=False)

                df3 = df3.drop_duplicates(subset='sequence')
                # df3 = df3.explode('proteins')



                df3[label] = df3['Intensity']
                df3['protein'] = df3['proteins']
                df3['peptide'] = df3['sequence']
                df3 = df3[['peptide', 'protein', label]]
                    
                if df_final is False:
                    label_basic = label
                    df_final = df3.reset_index(drop=True)
                else:
                    df_final = df_final.reset_index(drop=True).merge(df3.reset_index(drop=True), on='peptide', how='outer')
                    df_final.protein_x.fillna(value=df_final.protein_y, inplace=True)
                    df_final['protein'] = df_final['protein_x']

                    df_final = df_final.drop(columns=['protein_x', 'protein_y'])

            
    df_final = df_final.assign(protein=df_final['protein'].str.split(';')).explode('protein').reset_index(drop=True)

    df_final = df_final.set_index('peptide')
    df_final['proteins'] = df_final['protein']
    df_final = df_final.drop(columns=['protein'])
    # cols = df_final.columns.tolist()
    cols = [z for z in df_final.columns.tolist() if not z.startswith('mz_')]
    cols.remove('proteins')
    cols.insert(0, 'proteins')
    df_final = df_final[cols]

    all_lbls = all_s_lbls['S1'] + all_s_lbls['S2']

    for cc in all_lbls:
        df_final[cc] = df_final[cc].fillna(df_final[cc].min() / 2)

    df_final['p-value'] = list(ttest_ind(df_final[all_s_lbls['S1']].values.astype(float), df_final[all_s_lbls['S2']].values.astype(float), axis=1, nan_policy='omit', equal_var=True)[1])
    df_final['p-value'] = df_final['p-value'].astype(float)
    df_final['p-value'] = df_final['p-value'].fillna(2.0)

    df_final['S2_mean'] = df_final[all_s_lbls['S2']].mean(axis=1)
    df_final['S1_mean'] = df_final[all_s_lbls['S1']].mean(axis=1)

    df_final['FC'] = np.log2(df_final['S2_mean']/df_final['S1_mean'])

    FC_max = df_final['FC'].max()
    FC_min = df_final['FC'].min()

    df_final.loc[(pd.isna(df_final['S2_mean'])) & (~pd.isna(df_final['S1_mean'])), 'FC'] = FC_min * 2
    df_final.loc[(~pd.isna(df_final['S2_mean'])) & (pd.isna(df_final['S1_mean'])), 'FC'] = FC_max * 2

    df_final['decoy'] = df_final['proteins'].apply(lambda x: all(z.startswith('DECOY_') for z in x.split(';')))

    t_p = [0.01, 0.05, 0.25, 0.5]
    t_fc = [0.3, 1.0, 2.0]

    # t_p = [0.05, ]
    # t_fc = [1.0, ]

    total_up = defaultdict(float)
    total_down = defaultdict(float)

    for p_val_threshold, FC_threshold in itertools.product(t_p, t_fc):
        print(p_val_threshold, FC_threshold)
        
        df_final['sign'] = df_final['p-value'] <= p_val_threshold
        
        df_final['up'] = df_final['sign'] * df_final['FC'] >= FC_threshold
        df_final['down'] = df_final['sign'] * df_final['FC'] <= -FC_threshold
        
        protsN = df_final.groupby('proteins')['up'].count().to_dict()
        prots_up = df_final.groupby('proteins')['up'].sum()
        prots_down = df_final.groupby('proteins')['down'].sum()

        N_decoy_total = df_final['decoy'].sum()
        changed_decoy_total = df_final[(df_final['p-value'] <= p_val_threshold) & (df_final['decoy'])].shape[0]
        upreg_decoy_total = df_final[(df_final['p-value'] <= p_val_threshold) & (df_final['decoy']) & (df_final['FC'] >= FC_threshold)].shape[0]
        downreg_decoy_total = df_final[(df_final['p-value'] <= p_val_threshold) & (df_final['decoy']) & (df_final['FC'] <= -FC_threshold)].shape[0]


        p_up = upreg_decoy_total / N_decoy_total
        p_down = downreg_decoy_total / N_decoy_total
        print(N_decoy_total, changed_decoy_total, upreg_decoy_total, downreg_decoy_total, p_up, p_down)
        
        
        names_arr = np.array(list(protsN.keys()))
        v_arr = np.array(list(prots_up[k] for k in names_arr))
        n_arr = np.array(list(protsN[k] for k in names_arr))

        all_pvals = calc_sf_all(v_arr, n_arr, p_up)
        
        v_arr = np.array(list(prots_down[k] for k in names_arr))
        all_pvals_down = calc_sf_all(v_arr, n_arr, p_up)
        
        for z, dbname in zip(all_pvals, names_arr):
            total_up[dbname] += z
        for z, dbname in zip(all_pvals_down, names_arr):
            total_down[dbname] += z

    all_pvals = [total_up[dbname] for dbname in names_arr]
    all_pvals_down = [total_down[dbname] for dbname in names_arr]

    df_out = pd.DataFrame()
    df_out['score'] = all_pvals
    df_out['dbname'] = names_arr
    df_out['decoy'] = df_out['dbname'].str.startswith('DECOY_')
    df_out_f = aux.filter(df_out, fdr=0.05, key='score', is_decoy='decoy', reverse=True)
    df_out_f.to_csv(path_or_buf=args['out']+'_upreg.tsv', sep='\t', index=False)

    df_out = pd.DataFrame()
    df_out['score'] = all_pvals_down
    df_out['dbname'] = names_arr
    df_out['decoy'] = df_out['dbname'].str.startswith('DECOY_')
    df_out_f = aux.filter(df_out, fdr=0.05, key='score', is_decoy='decoy', reverse=True)
    df_out_f.to_csv(path_or_buf=args['out']+'_downreg.tsv', sep='\t', index=False)

if __name__ == '__main__':
    run()
