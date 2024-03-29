from __future__ import division
import argparse
import pandas as pd
import numpy as np
from scipy.stats import binom, ttest_ind
from scipy.optimize import curve_fit
import logging
from pyteomics import fasta
from collections import Counter, defaultdict
import random
random.seed(42)

logger = logging.getLogger(__name__)


def get_df_final(args, replace_label, allowed_peptides, allowed_prots_all, pep_RT=False, RT_threshold=False):

    df_final = False
    for i in range(1, 3, 1):
        sample_num = 'S%d' % (i, )
        if args.get(sample_num, 0):
            for z in args[sample_num]:
                label = sample_num + '_' + z.replace(replace_label, '')
                df3 = pd.read_table(z.replace(replace_label, '_PFMs.tsv'), usecols=['sequence', 'proteins', 'charge', 'ion_mobility', 'Intensity', 'RT'])

                if not args['allowed_peptides']:
                    df3['tmpseq'] = df3['sequence']
                    df3 = df3[df3['tmpseq'].apply(lambda x: x in allowed_peptides)]
                else:
                    df3 = df3[df3['sequence'].apply(lambda x: x in allowed_peptides)]


                df3 = df3[df3['proteins'].apply(lambda x: any(z in allowed_prots_all for z in x.split(';')))]
                df3['proteins'] = df3['proteins'].apply(lambda x: ';'.join([z for z in x.split(';') if z in allowed_prots_all]))

                df3['origseq'] = df3['sequence']
                df3['sequence'] = df3['sequence'] + df3['charge'].astype(int).astype(str) + df3['ion_mobility'].astype(str)

                if pep_RT is not False:
                    df3 = df3[df3.apply(lambda x: abs(pep_RT[x['sequence']] - x['RT']) <= 3 * RT_threshold, axis=1)]

                df3 = df3.sort_values(by='Intensity', ascending=False)

                df3 = df3.drop_duplicates(subset='sequence')

                df3[label] = df3['Intensity']
                df3['protein'] = df3['proteins']
                df3['peptide'] = df3['sequence']
                if pep_RT is False:
                    df3['RT_'+label] = df3['RT']
                    df3 = df3[['origseq', 'peptide', 'protein', label, 'RT_'+label]]
                else:
                    df3 = df3[['origseq', 'peptide', 'protein', label]]


                if df_final is False:
                    df_final = df3.reset_index(drop=True)
                else:
                    df_final = df_final.reset_index(drop=True).merge(df3.reset_index(drop=True), on='peptide', how='outer')
                    df_final.protein_x = df_final.protein_x.fillna(value=df_final.protein_y)
                    df_final.origseq_x = df_final.origseq_x.fillna(value=df_final.origseq_y)
                    df_final['protein'] = df_final['protein_x']
                    df_final['origseq'] = df_final['origseq_x']

                    df_final = df_final.drop(columns=['protein_x', 'protein_y'])
                    df_final = df_final.drop(columns=['origseq_x', 'origseq_y'])

    return df_final


def calc_sf_all(v, n, p):
    sf_values = -np.log10(binom.sf(v-1, n, p))
    sf_values[v <= 2] = 0
    sf_values[np.isinf(sf_values)] = 20
    sf_values[n == 0] = 0
    return sf_values

def noisygaus(x, a, x0, sigma, b):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b

def calibrate_mass(bwidth, mass_left, mass_right, true_md):

    bbins = np.arange(-mass_left, mass_right, bwidth)
    H1, b1 = np.histogram(true_md, bins=bbins)
    b1 = b1 + bwidth
    b1 = b1[:-1]


    popt, pcov = curve_fit(noisygaus, b1, H1, p0=[1, np.median(true_md), 1, 1])
    mass_shift, mass_sigma = popt[1], abs(popt[2])
    return mass_shift, mass_sigma, pcov[0][0]

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
    parser.add_argument('-fold_change', help='FC threshold standard deviations', default=2.0, type=float)
    parser.add_argument('-fold_change_abs', help='Use absolute log2 scale FC threshold instead of standard deviations', action='store_true')
    parser.add_argument('-bp', help='Experimental. Better percentage', default=80, type=int)
    parser.add_argument('-minl', help='Min peptide length for quantitation', default=7, type=int)
    parser.add_argument('-qval', help='qvalue threshold', default=0.05, type=float)
    parser.add_argument('-intensity_norm', help='Intensity normalization: 0-none, 1-median, 2-sum 1000 most intense peptides (default)', default=2, type=int)
    parser.add_argument('-all_proteins', help='use all proteins instead of FDR controlled', action='store_true')
    parser.add_argument('-all_pfms', help='use all PFMs instead of ML controlled', action='store_true')
    parser.add_argument('-allowed_peptides', help='path to allowed peptides')
    parser.add_argument('-allowed_proteins', help='path to allowed proteins')
    parser.add_argument('-d', '-db', help='path to uniprot fasta file for gene annotation')
    parser.add_argument('-prefix', help='Decoy prefix. Default DECOY_', default='DECOY_', type=str)
    args = vars(parser.parse_args())
    logging.basicConfig(format='%(levelname)9s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=logging.INFO)

    process_files(args)


def process_files(args):
    replace_label = '_proteins_full.tsv'
    decoy_prefix = args['prefix']

    fold_change = float(args['fold_change'])

    all_s_lbls = {}

    logger.info('Starting analysis...')

    allowed_prots = set()
    allowed_prots_all = set()
    allowed_peptides = set()

    # cnt0 = Counter()

    cnt_file = 0

    for i in range(1, 3, 1):
        sample_num = 'S%d' % (i, )
        if args[sample_num]:


            all_s_lbls[sample_num] = []

            for z in args[sample_num]:

                cnt_file += 1
                logger.debug('Processing file %d', cnt_file)

                label = sample_num + '_' + z.replace(replace_label, '')
                all_s_lbls[sample_num].append(label)

                if not args['allowed_proteins']:
                    if not args['all_proteins']:
                        df0 = pd.read_table(z.replace('_proteins_full.tsv', '_proteins.tsv'), usecols=['dbname', ])
                        allowed_prots.update(df0['dbname'])
                        allowed_prots.update([decoy_prefix + z for z in df0['dbname'].values])
                    else:
                        df0 = pd.read_table(z, usecols=['dbname', ])
                        allowed_prots.update(df0['dbname'])

                if not args['allowed_peptides']:
                    df0 = pd.read_table(z.replace('_proteins_full.tsv', '_PFMs_ML.tsv'), usecols=['seqs', 'qpreds', 'plen'])


                    if not args['all_pfms']:
                        df0 = df0[df0['qpreds'] <= 10]

                    df0 = df0[df0['plen'] >= args['minl']]
                    # df0['seqs'] = df0['seqs']
                    allowed_peptides.update(df0['seqs'])
                    # cnt0.update(df0['seqs'])

    if args['allowed_proteins']:
        allowed_prots = set(z.strip() for z in open(args['allowed_proteins'], 'r').readlines())
        allowed_prots.update([decoy_prefix + z for z in allowed_prots])

    if args['allowed_peptides']:
        allowed_peptides = set(z.strip() for z in open(args['allowed_peptides'], 'r').readlines())
    else:
        allowed_peptides = allowed_peptides

    logger.info('Total number of TARGET protein GROUPS: %d', len(allowed_prots) / 2)

    for i in range(1, 3, 1):
        sample_num = 'S%d' % (i, )
        if args.get(sample_num, 0):
            for z in args[sample_num]:
                df3 = pd.read_table(z.replace(replace_label, '_PFMs.tsv'), usecols=['sequence', 'proteins', ])
                df3 = df3[df3['sequence'].apply(lambda x: x in allowed_peptides)]

                df3_tmp = df3[df3['proteins'].apply(lambda x: any(z in allowed_prots for z in x.split(';')))]
                for dbnames in set(df3_tmp['proteins'].values):
                    for dbname in dbnames.split(';'):
                        allowed_prots_all.add(dbname)

    df_final = get_df_final(args, replace_label, allowed_peptides, allowed_prots_all, pep_RT=False, RT_threshold=False)

    logger.info('Total number of peptide sequences used in quantitation: %d', len(set(df_final['origseq'])))

    cols = [z for z in df_final.columns.tolist() if not z.startswith('mz_') and not z.startswith('RT_')]
    df_final = df_final[cols]

    df_final = df_final.set_index('peptide')


    all_lbls = all_s_lbls['S1'] + all_s_lbls['S2']

    df_final_copy = df_final.copy()

    custom_min_samples = int(args['min_samples'])
    if custom_min_samples == 0:
        custom_min_samples = int(len(all_lbls)/2)

    df_final = df_final_copy.copy()

    max_missing = len(all_lbls) - custom_min_samples

    logger.info('Allowed max number of missing values: %d', max_missing)

    df_final['nummissing'] = df_final.isna().sum(axis=1)
    df_final['nummissing_S1'] = df_final[all_s_lbls['S1']].isna().sum(axis=1)
    df_final['nummissing_S2'] = df_final[all_s_lbls['S2']].isna().sum(axis=1)
    df_final['nonmissing_S1'] = len(all_s_lbls['S1']) - df_final['nummissing_S1']
    df_final['nonmissing_S2'] = len(all_s_lbls['S2']) - df_final['nummissing_S2']
    df_final['nonmissing'] = df_final['nummissing'] <= max_missing

    df_final = df_final[df_final['nonmissing']]
    logger.info('Total number of PFMs passed missing values threshold: %d', len(df_final))


    df_final['S2_mean'] = df_final[all_s_lbls['S2']].mean(axis=1)
    df_final['S1_mean'] = df_final[all_s_lbls['S1']].mean(axis=1)
    df_final['FC_raw'] = np.log2(df_final['S2_mean']/df_final['S1_mean'])

    FC_max = df_final['FC_raw'].max()
    FC_min = df_final['FC_raw'].min()

    df_final.loc[(pd.isna(df_final['S2_mean'])) & (~pd.isna(df_final['S1_mean'])), 'FC_raw'] = FC_min
    df_final.loc[(~pd.isna(df_final['S2_mean'])) & (pd.isna(df_final['S1_mean'])), 'FC_raw'] = FC_max

    if args['intensity_norm'] == 2:
        for cc in all_lbls:
            df_final[cc] = df_final[cc] / df_final[cc].nlargest(1000).sum()

    elif args['intensity_norm'] == 1:
        for cc in all_lbls:
            df_final[cc] = df_final[cc] / df_final[cc].median()

    for slbl in ['1', '2']:
        S_len_current = len(all_s_lbls['S%s' % (slbl, )])
        df_final['S%s_mean' % (slbl, )] = df_final[all_s_lbls['S%s' % (slbl, )]].mean(axis=1)
        df_final['S%s_std' % (slbl, )] = np.log2(df_final[all_s_lbls['S%s' % (slbl, )]]).std(axis=1)

    df_final['S1_std'] = df_final['S1_std'].fillna(df_final['S2_std'])
    df_final['S2_std'] = df_final['S2_std'].fillna(df_final['S1_std'])

    idx_to_calc_initial_pval = df_final[['nonmissing_S1', 'nonmissing_S2']].min(axis=1) >= 2
    df_final.loc[idx_to_calc_initial_pval, 'p-value'] = list(ttest_ind(np.log10(df_final.loc[idx_to_calc_initial_pval, all_s_lbls['S1']].values.astype(float)), np.log10(df_final.loc[idx_to_calc_initial_pval, all_s_lbls['S2']].values.astype(float)), axis=1, nan_policy='omit', equal_var=True)[1])
    df_final['p-value'] = df_final['p-value'].astype(float)

    for cc in all_lbls:
        df_final[cc] = df_final[cc].fillna(df_final[cc].min())

    idx_missing_pval = pd.isna(df_final['p-value'])

    df_final.loc[idx_missing_pval, 'p-value'] = list(ttest_ind(np.log10(df_final.loc[idx_missing_pval, all_s_lbls['S1']].values.astype(float)), np.log10(df_final.loc[idx_missing_pval, all_s_lbls['S2']].values.astype(float)), axis=1, nan_policy='omit', equal_var=True)[1])

    df_final['p-value'] = df_final['p-value'].fillna(1.0)
    p_val_threshold = 0.1

    for cc in all_lbls:
        df_final[cc] = df_final[cc].fillna(df_final[cc].min())

    df_final['intensity_median'] = df_final[['S1_mean', 'S2_mean']].max(axis=1)
    df_final['iq'] = pd.qcut(df_final['intensity_median'], 5, labels=range(5)).fillna(0).astype(int)

    df_final['FC'] = np.log2(df_final['S2_mean']/df_final['S1_mean'])


    FC_max = df_final['FC'].max()
    FC_min = df_final['FC'].min()

    df_final_for_calib = df_final.copy()
    df_final_for_calib = df_final_for_calib[~pd.isna(df_final_for_calib['S1_mean'])]
    df_final_for_calib = df_final_for_calib[df_final_for_calib['FC'] <= FC_max/2]
    df_final_for_calib = df_final_for_calib[df_final_for_calib['FC'] >= FC_min/2]
    df_final_for_calib = df_final_for_calib[~pd.isna(df_final_for_calib['S2_mean'])]

    df_final.loc[(pd.isna(df_final['S2_mean'])) & (~pd.isna(df_final['S1_mean'])), 'FC'] = FC_min
    df_final.loc[(~pd.isna(df_final['S2_mean'])) & (pd.isna(df_final['S1_mean'])), 'FC'] = FC_max

    tmp = df_final_for_calib['FC']

    try:
        FC_mean, FC_std, covvalue_cor = calibrate_mass(0.05, -tmp.min(), tmp.max(), tmp)
        FC_mean2, FC_std2, covvalue_cor2 = calibrate_mass(0.1, -tmp.min(), tmp.max(), tmp)
        if not np.isinf(covvalue_cor2) and abs(FC_mean2) <= abs(FC_mean) / 10:
            FC_mean = FC_mean2
            FC_std = FC_std2
    except:
        FC_mean, FC_std, covvalue_cor = calibrate_mass(0.3, -tmp.min(), tmp.max(), tmp)

    if not args['fold_change_abs']:
        fold_change = FC_std * fold_change
    logger.info('Absolute FC threshold = %.2f +- %.2f', FC_mean, fold_change)

    df_final['decoy'] = df_final['protein'].apply(lambda x: all(z.startswith(decoy_prefix) for z in x.split(';')))

    df_final = df_final.assign(protein=df_final['protein'].str.split(';')).explode('protein').reset_index(drop=False)
    df_final['proteins'] = df_final['protein']
    df_final = df_final.drop(columns=['protein'])

    df_final = df_final.sort_values(by=['nummissing', 'intensity_median'], ascending=(True, False))
    df_final = df_final.drop_duplicates(subset=('origseq', 'proteins'))

    df_final['FC_corrected'] = df_final['FC'] - FC_mean
    df_final['FC_abs'] = df_final['FC_corrected'].abs()
    df_final = df_final.sort_values(by='FC_abs').reset_index(drop=True)
    df_final['FC_abs'] = df_final['FC_corrected']

    idx_stat = df_final['p-value'] <= p_val_threshold
    df_final['FC_gr_mean'] = df_final.groupby('proteins', group_keys=False)['FC_abs'].apply(lambda x: x.expanding().mean())

    neg_idx = (df_final['FC_corrected'] < 0)
    pos_idx = (df_final['FC_corrected'] >= 0)

    pos_idx_real = df_final[pos_idx].index
    neg_idx_real = df_final[neg_idx].index

    pos_idx_real_set = set(pos_idx_real)
    neg_idx_real_set = set(neg_idx_real)

    df_final_decoy = df_final[df_final['decoy']]

    FC_pools = dict()
    FC_pools['common'] = dict()

    df1_decoy_grouped_common = df_final.groupby('iq')
    for group_name, df_group in df1_decoy_grouped_common:
        FC_pools['common'][group_name] = list(df_group['FC_abs'].values)

    df_final['sign'] = False

    df1_grouped = df_final.groupby('proteins')

    for group_name, df_group in df1_grouped:

        prot_idx = df_group.index

        idx = sorted(list(prot_idx))
        idx_len = len(idx)

        loc_pos_values = df_final.loc[idx, 'FC_gr_mean'].values

        loc_pos_values = np.abs(loc_pos_values)

        pos_missing_list = list(df_final.loc[idx, 'iq'].values)
        better_res = np.array([0] * idx_len)
        for _ in range(100):
            random_list = [random.choice(FC_pools['common'][nm]) for nm in pos_missing_list]
            list_to_compare_current = np.cumsum(random_list) / np.arange(1, idx_len+1, 1)

            list_to_compare_current = np.abs(list_to_compare_current)

            better_res += loc_pos_values >= list_to_compare_current
        df_final.loc[idx, 'sign'] = better_res > args['bp']
        df_final.loc[idx, 'sign'] = df_final.loc[idx, 'sign'][::-1].cummin()[::-1]

    df_final.loc[df_final['p-value'] > p_val_threshold, 'sign'] = False

    df_final['up'] = df_final['sign'] * (df_final['FC_corrected'] > 0)
    df_final['down'] = df_final['sign'] * (df_final['FC_corrected'] < 0)

    cols = [z for z in df_final.columns.tolist() if not z.startswith('mz_') and not z.startswith('RT_')]
    cols.remove('proteins')
    cols.insert(0, 'proteins')
    df_final = df_final[cols]

    df_final.to_csv(path_or_buf=args['out']+'_quant_peptides.tsv', sep='\t', index=False)

    df_final = df_final.sort_values(by=['nummissing', 'intensity_median'], ascending=(True, False))
    df_final = df_final.drop_duplicates(subset=('origseq', 'proteins'))

    prot_to_peps = defaultdict(str)
    for dbname, pepseq in df_final.sort_values(by='origseq')[['proteins', 'origseq']].values:
        prot_to_peps[dbname] += pepseq

    all_peps_cnt = Counter(list(prot_to_peps.values()))
    peps_more_than_2 = set([k for k, v in all_peps_cnt.items() if v >= 2])

    pep_groups = {}
    protein_groups = {}
    cur_group = 1
    for dbname, pepseq in prot_to_peps.items():
        if pepseq not in peps_more_than_2:
            protein_groups[dbname] = cur_group
            cur_group += 1
        else:
            if pepseq not in pep_groups:
                pep_groups[pepseq] = cur_group
                protein_groups[dbname] = cur_group
                cur_group += 1
            else:
                protein_groups[dbname] = pep_groups[pepseq]

    del pep_groups
    del prot_to_peps
    del peps_more_than_2
    del all_peps_cnt

    up_dict = df_final.groupby('proteins')['up'].sum().to_dict()
    down_dict = df_final.groupby('proteins')['down'].sum().to_dict()

    ####### !!!!!!! #######
    df_final['up'] = df_final.apply(lambda x: x['up'] if up_dict.get(x['proteins'], 0) >= down_dict.get(x['proteins'], 0) else x['down'], axis=1)
    protsN = df_final.groupby('proteins')['up'].count().to_dict()

    prots_up = df_final.groupby('proteins')['up'].sum()
    decoy_df = df_final[df_final['decoy']].drop_duplicates(subset='origseq')

    N_decoy_total = len(decoy_df)
    upreg_decoy_total = decoy_df['up'].sum()

    N_nondecoy_total = (~df_final['decoy']).sum()

    p_up = upreg_decoy_total / N_decoy_total

    names_arr = np.array(list(protsN.keys()))

    logger.info('Total number of proteins used in quantitation: %d', sum(not z.startswith(decoy_prefix) for z in names_arr))
    logger.info('Total number of peptides: %d', len(df_final))
    logger.info('Total number of decoy peptides: %d', N_decoy_total)
    logger.info('Probability of random peptide to be differentially expressed: %.3f', p_up)

    v_arr = np.array(list(prots_up.get(k, 0) for k in names_arr))
    n_arr = np.array(list(protsN.get(k, 0) for k in names_arr))

    all_pvals = calc_sf_all(v_arr, n_arr, p_up)

    total_set = set()
    total_set_genes = set()

    FC_up_dict_basic = df_final.groupby('proteins')['FC_corrected'].median().to_dict()
    FC_up_dict_raw_basic = df_final.groupby('proteins')['FC_raw'].median().to_dict()

    df_final_up_idx = (df_final['up']>0)

    df_final.loc[df_final_up_idx, 'bestmissing'] = df_final.loc[df_final_up_idx, :].groupby('proteins')['nummissing'].transform('min')

    FC_up_dict2 = df_final.loc[df_final_up_idx, :].groupby('proteins')['FC_corrected'].median().to_dict()
    FC_up_dict_raw2 = df_final.loc[df_final_up_idx, :].groupby('proteins')['FC_raw'].median().to_dict()

    df_out = pd.DataFrame()
    df_out['score'] = all_pvals
    df_out['dbname'] = names_arr

    df_out['log2FoldChange(S2/S1)'] = df_out['dbname'].apply(lambda x: FC_up_dict2.get(x))
    df_out['log2FoldChange(S2/S1) no normalization'] = df_out['dbname'].apply(lambda x: FC_up_dict_raw2.get(x))

    df_out.loc[pd.isna(df_out['log2FoldChange(S2/S1)']), 'log2FoldChange(S2/S1)'] = df_out.loc[pd.isna(df_out['log2FoldChange(S2/S1)']), 'dbname'].apply(lambda x: FC_up_dict_basic.get(x))
    df_out.loc[pd.isna(df_out['log2FoldChange(S2/S1) no normalization']), 'log2FoldChange(S2/S1) no normalization'] = df_out.loc[pd.isna(df_out['log2FoldChange(S2/S1) no normalization']), 'dbname'].apply(lambda x: FC_up_dict_raw_basic.get(x))


    df_out.loc[:, 'log2FoldChange(S2/S1) using all peptides'] = df_out.loc[:, 'dbname'].apply(lambda x: FC_up_dict_basic.get(x))
    df_out.loc[:, 'log2FoldChange(S2/S1) using all peptides and no normalization'] = df_out.loc[:, 'dbname'].apply(lambda x: FC_up_dict_raw_basic.get(x))

    df_out['differentially expressed peptides'] = v_arr
    df_out['identified peptides'] = n_arr

    df_out['decoy'] = df_out['dbname'].str.startswith(decoy_prefix)

    df_out = df_out[~df_out['decoy']]

    df_out['protname'] = df_out['dbname'].apply(lambda x: x.split('|')[1] if '|' in x else x)
    df_out['protein_quant_group'] = df_out['dbname'].apply(lambda x: protein_groups[x])

    genes_map = {}
    if args['d']:
        for prot, protseq in fasta.read(args['d']):
            if decoy_prefix not in prot:
                try:
                    prot_name = prot.split('|')[1]
                except:
                    prot_name = prot
                try:
                    gene_name = prot.split('GN=')[1].split(' ')[0]
                except:
                    gene_name = prot
                genes_map[prot_name] = gene_name

        df_out['gene'] = df_out['protname'].apply(lambda x: genes_map[x])

    else:
        df_out['gene'] = df_out['protname']


    qval_threshold = args['qval']

    min_matched_peptides = 3

    df_out = df_out.sort_values(by='score', ascending=False).reset_index(drop=True)

    df_out['FC_pass'] = False
    df_out['FC_pass'] = df_out['log2FoldChange(S2/S1)'].abs() >= fold_change


    BH_idx = (df_out['identified peptides'] >= min_matched_peptides) & (df_out['FC_pass'])

    BH_idx_pos = (df_out['log2FoldChange(S2/S1)'] >= 0) & (df_out['identified peptides'] >= min_matched_peptides)
    BH_idx_neg = (df_out['log2FoldChange(S2/S1)'] < 0) & (df_out['identified peptides'] >= min_matched_peptides)

    # df_out['p-value'] = 1.0
    df_out['p-value'] = 10**(-df_out['score'])
    df_out['BH_pass'] = False

    df_out_BH_multiplier = len(set(df_out[BH_idx]['protein_quant_group']))
    lbl_to_use = 'protein_quant_group'

    current_rank = 0
    BH_threshold_array = []
    added_groups = set()
    for z in df_out[BH_idx][lbl_to_use].values:
        if z not in added_groups:
            added_groups.add(z)
            current_rank += 1
        BH_threshold_array.append(-np.log10(current_rank * qval_threshold / df_out_BH_multiplier))
    df_out.loc[BH_idx, 'BH_threshold'] = BH_threshold_array

    df_out.loc[BH_idx, 'BH_pass'] = df_out.loc[BH_idx, 'score'] >= df_out.loc[BH_idx, 'BH_threshold']
    df_out.loc[BH_idx, 'FDR_pass'] = df_out.loc[BH_idx, 'score'] >= -np.log10(args['qval'])
    # df_out.loc[BH_idx, 'BH_pass'] = df_out.loc[BH_idx, 'score'] >= -np.log10(args['qval'])

    # score_threshold = df_out.loc[(df_out['BH_pass']) & (BH_idx)]['score'].min()
    # df_out.loc[BH_idx, 'BH_pass'] = df_out.loc[BH_idx, 'score'] >= score_threshold

    df_out = df_out.drop(columns = {'decoy'})

    df_out = df_out[['score', 'p-value', 'dbname', 'log2FoldChange(S2/S1)', 'differentially expressed peptides',
                    'identified peptides', 'log2FoldChange(S2/S1) no normalization', 'log2FoldChange(S2/S1) using all peptides',
                    'log2FoldChange(S2/S1) using all peptides and no normalization', 'protname', 'protein_quant_group', 'gene', 'FC_pass', 'FDR_pass', 'BH_pass',]]# 'BH_threshold']]

    df_out.to_csv(path_or_buf=args['out']+'_quant_full.tsv', sep='\t', index=False)

    # df_out_f = df_out[(df_out['FDR_pass']) & (df_out['FC_pass'])]
    df_out_f = df_out[(df_out['BH_pass']) & (df_out['FC_pass'])]

    df_out_f.to_csv(path_or_buf=args['out']+'.tsv', sep='\t', index=False)

    for z in set(df_out_f['dbname']):
        try:
            prot_name = z.split('|')[1]
        except:
            prot_name = z

        gene_name = genes_map.get(prot_name, prot_name)

        total_set.add(prot_name)
        total_set_genes.add(gene_name)

    logger.info('Total number of significantly changed proteins: %d', len(total_set))
    logger.info('Total number of significantly changed genes: %d', len(total_set_genes))

    f1 = open(args['out'] + '_proteins_for_stringdb.txt', 'w')
    for z in total_set:
        f1.write(z + '\n')
    f1.close()

    f1 = open(args['out'] + '_genes_for_stringdb.txt', 'w')
    for z in total_set_genes:
        f1.write(z + '\n')
    f1.close()

if __name__ == '__main__':
    run()
