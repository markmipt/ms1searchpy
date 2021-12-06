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
    sf_values[v <= 1] = 0
    sf_values[np.isinf(sf_values)] = 10#max(sf_values[~np.isinf(sf_values)]) * 2
    sf_values[n == 0] = 0
    # sf_values = binom.sf(v-1, n, p)
    # sf_values[np.isneginf(sf_values)] = min(sf_values[~np.isinf(sf_values)])
    # sf_values[n <= 3] = 1.0
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
    parser.add_argument('-out', help='name of DirectMS1quant output file', default='directms1quant_out')
    parser.add_argument('-min_samples', help='minimum number of samples for peptide usage. 0 means 50%% of input files', default=0)
    parser.add_argument('-fold_change', help='FC threshold in log2 scale', default=0.5)
    parser.add_argument('-qval', help='qvalue threshold', default=0.05)
    parser.add_argument('-intensity_norm', help='use intensity normalization', action='store_true')
    parser.add_argument('-all_proteins', help='use all proteins instead of FDR controlled', action='store_true')
    parser.add_argument('-all_pfms', help='use all PFMs instead of ML controlled', action='store_true')
    args = vars(parser.parse_args())

    replace_label = '_proteins_full.tsv'

    fold_change = float(args['fold_change'])

    df_final = False

    all_s_lbls = {}

    allowed_prots = set()
    all_peptides = set()
    allowed_peptides = set()

    for i in range(1, 3, 1):
        sample_num = 'S%d' % (i, )
        if args[sample_num]:
            for z in args[sample_num]:
                if not args['all_proteins']:
                    df0 = pd.read_table(z.replace('_proteins_full.tsv', '_proteins.tsv'))
                    allowed_prots.update(df0['dbname'])
                    allowed_prots.update(['DECOY_' + z for z in df0['dbname'].values])
                else:
                    df0 = pd.read_table(z)
                    allowed_prots.update(df0['dbname'])

                df0 = pd.read_table(z.replace('_proteins_full.tsv', '_PFMs_ML.tsv'))

                all_peptides.update(df0['seqs'])
                if not args['all_pfms']:
                    df0 = df0[df0['qpreds'] <= 10]
                    allowed_peptides.update(df0['seqs'])
                else:
                    allowed_peptides = all_peptides

    print('Total number of target proteins %d' % (len(allowed_prots)/2, ))
    print('Total number of peptides %d' % (len(allowed_peptides), ))


    for i in range(1, 3, 1):
        sample_num = 'S%d' % (i, )
        if args.get(sample_num, 0):
            all_s_lbls[sample_num] = []
            for z in args[sample_num]:
                label = z.replace(replace_label, '')
                all_s_lbls[sample_num].append(label)
                df1 = pd.read_table(z)
                df3 = pd.read_table(z.replace(replace_label, '_PFMs.tsv'))

                df3 = df3[df3['sequence'].apply(lambda x: x in allowed_peptides)]

                df3['plen'] = df3['sequence'].apply(lambda x: len(x))

                df3 = df3[df3['proteins'].apply(lambda x: any(z in allowed_prots for z in x.split(';')))]
                df3['proteins'] = df3['proteins'].apply(lambda x: ';'.join([z for z in x.split(';') if z in allowed_prots]))

                df3['sequence'] = df3['sequence']# + df3['charge'].astype(int).astype(str)
                # df3['sequence'] = df3['sequence'] + df3['charge'].astype(int).astype(str) + df3['ion_mobility'].astype(str)

                # df3['Intensity'] = df3.groupby('sequence')['Intensity'].transform('sum')

                df3 = df3.sort_values(by='Intensity', ascending=False)
                # df3 = df3.sort_values(by=('qpreds', 'Intensity'), ascending=(True, False))
                # df3 = df3.sort_values(by='nScans', ascending=False)

                df3 = df3.drop_duplicates(subset='sequence')
                # df3 = df3.explode('proteins')

                # df3['Intensity'] = np.log2(df3['Intensity'])
                # df3['Intensity'] = df3['Intensity'] - df3['Intensity'].median()

                # df3['Intensity'] = df3['Intensity'] / df3['Intensity'].median()

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

    df_final_copy = df_final.copy()

    df_res_list = []

    custom_min_samples = int(args['min_samples'])
    if custom_min_samples == 0:
        custom_min_samples = int(len(all_lbls)/2)

    df_final = df_final_copy.copy()

    max_missing = len(all_lbls) - custom_min_samples

    print('Allowed max number of missing values: %d' % (max_missing, ))

    df_final['nonmissing'] = df_final.isna().sum(axis=1) <= max_missing

    df_final = df_final[df_final['nonmissing']]
    print('Total number of peptides passed missing values threshold %d' % (len(df_final), ))

    if args['intensity_norm']:
        for cc in all_lbls:
            df_final[cc] = df_final[cc] / df_final[cc].median()

    for cc in all_lbls:
        df_final[cc] = df_final[cc].fillna(df_final[cc].min())

    df_final['p-value'] = list(ttest_ind(df_final[all_s_lbls['S1']].values.astype(float), df_final[all_s_lbls['S2']].values.astype(float), axis=1, nan_policy='omit', equal_var=True)[1])
    df_final['p-value'] = df_final['p-value'].astype(float)
    df_final['p-value'] = df_final['p-value'].fillna(1.0)

    df_final['S2_mean'] = df_final[all_s_lbls['S2']].mean(axis=1)
    df_final['S1_mean'] = df_final[all_s_lbls['S1']].mean(axis=1)

    df_final['FC'] = np.log2(df_final['S2_mean']/df_final['S1_mean'])

    FC_max = df_final['FC'].max()
    FC_min = df_final['FC'].min()

    df_final.loc[(pd.isna(df_final['S2_mean'])) & (~pd.isna(df_final['S1_mean'])), 'FC'] = FC_min * 2
    df_final.loc[(~pd.isna(df_final['S2_mean'])) & (pd.isna(df_final['S1_mean'])), 'FC'] = FC_max * 2

    df_final['decoy'] = df_final['proteins'].apply(lambda x: all(z.startswith('DECOY_') for z in x.split(';')))

    from scipy.stats import scoreatpercentile
    from scipy.optimize import curve_fit
    from scipy import exp
    def noisygaus(x, a, x0, sigma, b):
        return a * exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b

    def calibrate_mass(bwidth, mass_left, mass_right, true_md):

        bbins = np.arange(-mass_left, mass_right, bwidth)
        H1, b1 = np.histogram(true_md, bins=bbins)
        b1 = b1 + bwidth
        b1 = b1[:-1]


        popt, pcov = curve_fit(noisygaus, b1, H1, p0=[1, np.median(true_md), 1, 1])
        mass_shift, mass_sigma = popt[1], abs(popt[2])
        return mass_shift, mass_sigma, pcov[0][0]

    FC_l = -fold_change
    FC_r = fold_change

    # df_final.to_csv(path_or_buf=args['out']+'_final_unfilt.tsv', sep='\t', index=False)

    total_up = defaultdict(float)
    total_down = defaultdict(float)

    p_val_threshold = 0.05
        
    df_final['sign'] = df_final['p-value'] <= p_val_threshold
    
    df_final['up'] = df_final['sign'] * (df_final['FC'] >= FC_r)
    df_final['down'] = df_final['sign'] * (df_final['FC'] <= FC_l)

    up_dict = df_final.groupby('proteins')['up'].sum().to_dict()
    down_dict = df_final.groupby('proteins')['down'].sum().to_dict()

    ####### !!!!!!! #######
    df_final['up'] = df_final.apply(lambda x: x['up'] if up_dict.get(x['proteins'], 0) >= down_dict.get(x['proteins'], 0) else x['down'], axis=1)

    protsN = df_final.groupby('proteins')['up'].count().to_dict()

    prots_up = df_final.groupby('proteins')['up'].sum()
    prots_down = df_final.groupby('proteins')['down'].sum()

    N_decoy_total = df_final['decoy'].sum()
    changed_decoy_total = df_final[(df_final['p-value'] <= p_val_threshold) & (df_final['decoy'])].shape[0]

    upreg_decoy_total = df_final[df_final['decoy']]['up'].sum()
    downreg_decoy_total = df_final[df_final['decoy']]['down'].sum()


    p_up = upreg_decoy_total / N_decoy_total
    p_down = downreg_decoy_total / N_decoy_total
    print(N_decoy_total, changed_decoy_total, upreg_decoy_total, downreg_decoy_total, p_up, p_down)
    
    
    names_arr = np.array(list(protsN.keys()))
    v_arr = np.array(list(prots_up.get(k, 0) for k in names_arr))
    n_arr = np.array(list(protsN.get(k, 0) for k in names_arr))

    all_pvals = calc_sf_all(v_arr, n_arr, p_up)
    
    for z, dbname in zip(all_pvals, names_arr):
        total_up[dbname] += z

    all_pvals = [total_up[dbname] for dbname in names_arr]

    total_set = set()

    FC_up_dict = df_final[df_final['up']>0].groupby('proteins')['FC'].median().to_dict()

    df_out = pd.DataFrame()
    df_out['score'] = all_pvals
    df_out['dbname'] = names_arr

    df_out['FC'] = df_out['dbname'].apply(lambda x: FC_up_dict.get(x))

    df_res_list.append(df_out)

    df_out = False

    FC_dict_final = defaultdict(list)
    score_dict_final = defaultdict(float)

    for df_tmp in df_res_list:
        for dbname, score, FC in df_tmp[['dbname', 'score', 'FC']].values:
            FC_dict_final[dbname].append(FC)
            score_dict_final[dbname] += score

    df_out = pd.DataFrame()
    names_arr = [k for k in score_dict_final.keys()]
    all_pvals = [score_dict_final[k] for k in names_arr]
    all_FCs = [np.mean(FC_dict_final[k]) for k in names_arr]
    df_out['score'] = all_pvals
    df_out['v_arr'] = v_arr
    df_out['n_arr'] = n_arr
    df_out['dbname'] = names_arr
    df_out['FC'] = all_FCs

    df_out['decoy'] = df_out['dbname'].str.startswith('DECOY_')
    print((df_out['decoy']).sum(), (~df_out['decoy']).sum())

    df_out = df_out[~df_out['decoy']]

    try:
        FC_mean, FC_std, covvalue_cor = calibrate_mass(0.1, -df_out['FC'].min(), df_out['FC'].max(), df_out['FC'])
    except:
        FC_mean, FC_std, covvalue_cor = calibrate_mass(0.3, -df_out['FC'].min(), df_out['FC'].max(), df_out['FC'])
    print('df_out_FC', FC_mean, FC_std)

    df_out = df_out[df_out['FC'].abs() >= fold_change]

    qval_threshold = args['qval']

    df_out = df_out.sort_values(by='score', ascending=False)
    df_out['BH_threshold'] = -np.log10(df_out['score'].rank(ascending=False, method='max') * qval_threshold / len(df_out))
    df_out['pass'] = df_out['score'] > df_out['BH_threshold']
    df_out['p-value'] = 10**(-df_out['score'])
    score_threshold = df_out[~df_out['pass']]['score'].max()

    df_out.to_csv(path_or_buf=args['out']+'_proteins_full.tsv', sep='\t', index=False)

    df_out_f = df_out[df_out['score'] >= score_threshold]

    df_out_f.to_csv(path_or_buf=args['out']+'.tsv', sep='\t', index=False)

    total_set.update([z.split('|')[1] for z in set(df_out_f['dbname'])])

    f1 = open(args['out'] + '_proteins_for_stringdb.txt', 'w')
    for z in total_set:
        f1.write(z + '\n')
    f1.close()

if __name__ == '__main__':
    run()
