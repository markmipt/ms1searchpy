from __future__ import division
import argparse
import pandas as pd
import numpy as np
from scipy.stats import binom, ttest_ind, scoreatpercentile, percentileofscore, norm
from scipy.optimize import curve_fit
import logging
from pyteomics import fasta
from collections import Counter, defaultdict
import random
import subprocess
import os
random.seed(42)

logger = logging.getLogger(__name__)


def weighted_quantiles_interpolate(values, weights, quantiles=0.5):
    i = np.argsort(values)
    c = np.cumsum(weights[i])
    q = np.searchsorted(c, quantiles * c[-1])
    q_plus1 = np.clip(q + 1, a_min=None, a_max=values.shape[0] - 1)
    return np.where(c[q]/c[-1] == quantiles, 0.5 * (values[i[q]] + values[i[q_plus1]]), values[i[q]])


def get_df_final(args, replace_label, allowed_peptides, allowed_prots_all, pep_RT=False, RT_threshold=False, prot_spc=False):
    df_final = False
    for i in range(1, 3, 1):
        sample_num = 'S%d' % (i, )
        if args.get(sample_num, 0):
            for z in args[sample_num]:
                label = sample_num + '_' + z.replace(replace_label, '')
                df3 = pd.read_table(z.replace(replace_label, '_PFMs_ML.tsv'), usecols=['seqs', 'proteins', 'ch', 'im', 'Is', 'rt', 'qpreds', 'preds'])
                df3 = df3.rename(columns={'seqs': 'sequence', 'ch': 'charge', 'im': 'ion_mobility', 'Is': 'Intensity', 'rt': 'RT'})

                if pep_RT is False:
                    df3 = df3[df3['qpreds'] <= 10]

                if not args['allowed_peptides']:
                    df3['tmpseq'] = df3['sequence']
                    df3 = df3[df3['tmpseq'].apply(lambda x: x in allowed_peptides)]
                else:
                    df3 = df3[df3['sequence'].apply(lambda x: x in allowed_peptides)]


                df3 = df3[df3['proteins'].apply(lambda x: any(z in allowed_prots_all for z in x.split(';')))]
                df3['proteins'] = df3['proteins'].apply(lambda x: ';'.join([z for z in x.split(';') if z in allowed_prots_all]))

                df3['origseq'] = df3['sequence']
                df3['sequence'] = df3['sequence'] + df3['charge'].astype(int).astype(str) + df3['ion_mobility'].astype(str)


                if pep_RT is False:
                    df3 = df3.sort_values(by='preds')
                    df3 = df3.drop_duplicates(subset='sequence')

                if pep_RT is not False:
                    df3 = df3[df3['sequence'].apply(lambda x: x in pep_RT)]
                    df3['RT diff'] = df3.apply(lambda x: pep_RT[x['sequence']] - x['RT'], axis=1)
                    RT_shift, RT_threshold_l, RT_threshold_r = RT_threshold['RT_'+label]
                    df3['RT diff'] = df3['RT diff'] - RT_shift
                    df3 = df3[df3.apply(lambda x: RT_threshold_l <= (pep_RT[x['sequence']] - x['RT'] - RT_shift) <= RT_threshold_r, axis=1)]

                df3 = df3.sort_values(by='Intensity', ascending=False)
                df3 = df3.drop_duplicates(subset='sequence')

                df3[label] = df3['Intensity']
                df3['protein'] = df3['proteins']
                df3['peptide'] = df3['sequence']
                if pep_RT is False:
                    df3['RT_'+label] = df3['RT']
                    df3 = df3[['origseq', 'peptide', 'protein', label, 'RT_'+label]]
                else:

                    df3['RT_'+label] = df3['RT']
                    df3 = df3[['origseq', 'peptide', 'protein', label, 'RT_'+label]]


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


def calc_sf_all(v, n, p, min_peptides=3):
    sf_values = -np.log10(binom.sf(v-1, n, p))
    sf_values[v <= min_peptides-1] = 0
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
    parser.add_argument('-min_samples_group', help='minimum number of samples "best" group for peptide usage', default=0)
    parser.add_argument('-fold_change', help='FC threshold standard deviations', default=2.5, type=float)
    parser.add_argument('-fold_change_abs', help='Use absolute log2 scale FC threshold instead of standard deviations', action='store_true')
    parser.add_argument('-minl', help='Min peptide length for quantitation', default=7, type=int)
    parser.add_argument('-qval', help='qvalue threshold', default=0.05, type=float)
    parser.add_argument('-intensity_norm', help='Intensity normalization: 0-none, 1-median, 2-sum 1000 most intense peptides, 3-median by RT interval (default)', default=3, type=int)
    parser.add_argument('-all_proteins', help='use all proteins instead of FDR controlled', action='store_true')
    parser.add_argument('-all_pfms', help='use all PFMs instead of ML controlled', action='store_true')
    parser.add_argument('-do_not_remove_missing', help='do_not_remove_missing', action='store_true')
    parser.add_argument('-do_not_use_RT_alignment', help='do_not_use_RT_alignment', action='store_true')
    parser.add_argument('-allowed_peptides', help='path to allowed peptides')
    parser.add_argument('-allowed_proteins', help='path to allowed proteins')
    parser.add_argument('-protein_shifts', help='path to protein shifts (Used to normalize TPP experiments at initial temperature)')
    parser.add_argument('-d', '-db', help='path to uniprot fasta file for gene annotation')
    parser.add_argument('-prefix', help='Decoy prefix. Default DECOY_', default='DECOY_', type=str)
    parser.add_argument('-min_matched_peptides', help='Min significant peptides for reported DEPs', default=3, type=int)
    parser.add_argument('-legacy', help='Use legacy Directms1Quant workflow', action='store_true')
    args = vars(parser.parse_args())
    logging.basicConfig(format='%(levelname)9s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=logging.INFO)

    process_files(args)


def process_files(args):
    replace_label = '_proteins_full.tsv'
    decoy_prefix = args['prefix']

    fold_change = float(args['fold_change'])

    min_matched_peptides = args['min_matched_peptides']

    all_s_lbls = {}

    logger.info('Starting analysis...')

    if not args['legacy']:
        logger.info('Starting Directms1Quant2 preparations...')
        from . import combine_proteins
        all_files_list = []
        for i in range(1, 3, 1):
            sample_num = 'S%d' % (i, )
            if args[sample_num]:
                for z in args[sample_num]:
                    all_files_list.append(z)

        new_args = dict()
        outname_for_combined = args['out'] + '_union_IDs'
        new_args['out'] = outname_for_combined
        new_args['file'] = all_files_list
        new_args['fdr'] = 5.0
        new_args['prefix'] = args['prefix']

        combine_proteins.base_func(new_args, logger)

        args['allowed_proteins'] = outname_for_combined + '.features_proteins.tsv'

        tmp = pd.DataFrame(columns=['File Name', 'group', 'condition', 'BatchMS', 'vs'])

        for i in range(1, 3, 1):
            sample_num = 'S%d' % (i, )
            if args[sample_num]:
                for z in args[sample_num]:
                    tmp.loc[len(tmp), :] = [os.path.basename(z).replace('.features_proteins_full.tsv', ''), sample_num, '1', '1', '1']

        outname_for_multi = args['out'] + '_directms1quant_multi'
        outname_for_multi_samples = args['out'] + '_directms1quant_multi_samples.tsv'
        tmp[['File Name', 'group', 'condition', 'BatchMS', 'vs']].to_csv(outname_for_multi_samples, index=False, sep='\t')

        new_args = dict()
        pdir_tmp = os.path.dirname(z)
        if not pdir_tmp:
            pdir_tmp = os.getcwd()
        new_args['pdir'] = pdir_tmp
        print(new_args['pdir'])
        new_args['samples'] = outname_for_multi_samples
        new_args['out'] = outname_for_multi
        new_args['norm'] = 1
        new_args['proteins_for_figure'] = ''
        new_args['figdir'] = ''
        new_args['max_missing'] = 0.5
        new_args['prefix'] = args['prefix']
        new_args['start_stage'] = 1
        new_args['plot_figures'] = 0

        from . import directms1quantmulti
        directms1quantmulti.process_files(new_args, logger)

        dfx = pd.read_table(os.path.join(pdir_tmp, outname_for_multi) + '_proteins_LFQ.tsv')
        tmp_S1 = dfx[dfx['group'] == 'S1']
        tmp_S2 = dfx[dfx['group'] == 'S2']
        from scipy.stats import ttest_ind

        prots_pval2 = dict()

        banned_set = set(['File Name', 'group', 'condition', 'BatchMS', 'vs', 'sample', 'replicate', 'sample+condition'])
        for cc in dfx.columns:
            if cc not in banned_set:
                ar1 = tmp_S1[cc].values
                ar2 = tmp_S2[cc].values
                prots_pval2[cc] = ttest_ind(np.power(2, ar1), np.power(2, ar2))[1]
        logger.info('Directms1Quant2 preparations were finished...')

    allowed_prots = set()
    allowed_prots_all = set()
    allowed_peptides = set()
    cnt_file = 0
    prot_scores = defaultdict(list)
    s_koeff = 0

    for i in range(1, 3, 1):
        sample_num = 'S%d' % (i, )
        if args[sample_num]:
            all_s_lbls[sample_num] = []
            for z in args[sample_num]:
                cnt_file += 1
                logger.debug('Processing file %d', cnt_file)

                label = sample_num + '_' + z.replace(replace_label, '')
                all_s_lbls[sample_num].append(label)

                s_koeff += 1
                df0 = pd.read_table(z, usecols=['dbname', 'score'])
                for dbname, score in df0.values:
                    prot_scores[dbname].append(score)
                    
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
                    allowed_peptides.update(df0['seqs'])

    if args['allowed_proteins']:
        try:
            ap = pd.read_table(args['allowed_proteins'], usecols=['dbname', ])
            allowed_prots = set(ap['dbname'])
        except:
            allowed_prots = set(z.strip() for z in open(args['allowed_proteins'], 'r').readlines())
        allowed_prots.update([decoy_prefix + z for z in allowed_prots])

    if args['allowed_peptides']:
        allowed_peptides = set(z.strip() for z in open(args['allowed_peptides'], 'r').readlines())
    else:
        allowed_peptides = allowed_peptides

    logger.info('Total number of TARGET protein GROUPS: %d', len(allowed_prots) / 2)

    prot_spc = dict()
    for k, v in prot_scores.items():
        prot_spc[k] = sum(v) / s_koeff

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

    df_final = get_df_final(args, replace_label, allowed_peptides, allowed_prots_all, pep_RT=False, RT_threshold=False, prot_spc=prot_spc)

    if not args['do_not_use_RT_alignment']:
        rt_ccols = [z for z in df_final.columns.tolist() if z.startswith('RT_')]
        pep_RT = df_final.set_index('peptide')[rt_ccols].median(axis=1).to_dict()
        RT_threshold = dict()

        for cc in rt_ccols:
            dfqqq = df_final[~pd.isna(df_final[cc])].copy()
            dfqqq['RT diff'] = dfqqq.apply(lambda x: pep_RT[x['peptide']] - x[cc], axis=1)
            RT_shift = dfqqq['RT diff'].median()
            dfqqq['RT diff'] = dfqqq['RT diff'] - RT_shift
            RT_threshold_l = scoreatpercentile(dfqqq['RT diff'], 10) * 2
            RT_threshold_r = scoreatpercentile(dfqqq['RT diff'], 90) * 2
            RT_threshold[cc] = (RT_shift, RT_threshold_l, RT_threshold_r)

        df_final = get_df_final(args, replace_label, allowed_peptides, allowed_prots_all, pep_RT=pep_RT, RT_threshold=RT_threshold, prot_spc=prot_spc)

    logger.info('Total number of peptide sequences used in quantitation: %d', len(set(df_final['origseq'])))



    all_lbls = all_s_lbls['S1'] + all_s_lbls['S2']


    if args['intensity_norm'] == 3:
        df_final['RT_median'] = df_final[rt_ccols].median(axis=1)
        num_rt_groups = int(len(df_final) / 250)
        num_rt_groups = min(50, num_rt_groups)
        num_rt_groups = max(1, num_rt_groups)
        df_final['q_RT'] = pd.qcut(df_final['RT_median'], num_rt_groups, labels=range(num_rt_groups)).astype(int)

        norm_dict = dict()

        df_to_use = df_final
        non_missing_peptides_best = 0
        for cc in all_lbls:
            
            non_missing_peptides = (~pd.isna(df_final[cc])).sum()
            if non_missing_peptides >= non_missing_peptides_best:
                cc1 = cc
                non_missing_peptides_best = non_missing_peptides

        for cc in all_lbls:
            tmp_df = df_to_use
            ar1 = tmp_df[cc1].values
            ar2 = tmp_df[cc].values
            ar_ratio = np.array(tmp_df[cc] / tmp_df[cc1])
            idx_non_missing = (~np.isnan(ar_ratio))
            ar_ratio = ar_ratio[idx_non_missing]
            ar2 = ar2[idx_non_missing]
            koef2_base = weighted_quantiles_interpolate(ar_ratio, np.sqrt(ar2), 0.5)
            print('median intensity for sample %s: %.3f' % (cc, koef2_base))

            for RT_int in set(df_to_use['q_RT']):
                tmp_df = df_to_use[df_to_use['q_RT'] == RT_int]
                ar1 = tmp_df[cc1].values
                ar2 = tmp_df[cc].values
                ar_ratio = np.array(tmp_df[cc] / tmp_df[cc1])
                idx_non_missing = (~np.isnan(ar_ratio))
                ar_ratio = ar_ratio[idx_non_missing]
                ar2 = ar2[idx_non_missing]
                if len(ar_ratio) == 0:
                    norm_dict[(cc, RT_int)] = koef2_base
                else:
                    koef2 = weighted_quantiles_interpolate(ar_ratio, np.sqrt(ar2), 0.5)
                    norm_dict[(cc, RT_int)] = koef2



    df_final.to_csv(path_or_buf=args['out']+'_quant_peptides_raw.tsv', sep='\t', index=False, float_format="%.4g")

    cols = [z for z in df_final.columns.tolist() if not z.startswith('mz_') and not z.startswith('RT_')]
    df_final = df_final[cols]

    df_final = df_final.set_index('peptide')

    df_final_copy = df_final.copy()

    custom_min_samples = int(args['min_samples'])
    if custom_min_samples == 0:
        custom_min_samples = int(len(all_lbls)/2)

    custom_min_samples_group = int(args['min_samples_group'])

    df_final = df_final_copy.copy()

    max_missing = len(all_lbls) - custom_min_samples

    logger.info('Allowed max number of missing values: %d', max_missing)

    df_final['nummissing'] = df_final.isna().sum(axis=1)
    df_final['nummissing_S1'] = df_final[all_s_lbls['S1']].isna().sum(axis=1)
    df_final['nummissing_S2'] = df_final[all_s_lbls['S2']].isna().sum(axis=1)
    df_final['nonmissing_S1'] = len(all_s_lbls['S1']) - df_final['nummissing_S1']
    df_final['nonmissing_S2'] = len(all_s_lbls['S2']) - df_final['nummissing_S2']
    df_final['nonmissing'] = df_final['nummissing'] <= max_missing

    if not args['do_not_remove_missing']:
        df_final = df_final[df_final['nonmissing']]
    logger.info('Total number of PFMs: %d', len(df_final))
    logger.info('Total number of PFMs passed missing values threshold: %d', len(df_final[df_final['nonmissing']]))


    df_final['S2_mean'] = df_final[all_s_lbls['S2']].mean(axis=1)
    df_final['S1_mean'] = df_final[all_s_lbls['S1']].mean(axis=1)
    df_final['FC_raw'] = np.log2(df_final['S2_mean']/df_final['S1_mean'])

    FC_max = df_final[df_final['nonmissing']]['FC_raw'].max()
    FC_min = df_final[df_final['nonmissing']]['FC_raw'].min()

    df_final.loc[(pd.isna(df_final['S2_mean'])) & (~pd.isna(df_final['S1_mean'])), 'FC_raw'] = FC_min
    df_final.loc[(~pd.isna(df_final['S2_mean'])) & (pd.isna(df_final['S1_mean'])), 'FC_raw'] = FC_max

    if args['intensity_norm'] == 3:
        for cc in all_lbls:
            df_final[cc] = df_final.apply(lambda x: x[cc] / norm_dict[(cc, x['q_RT'])], axis=1)
    elif args['intensity_norm'] == 2:
        for cc in all_lbls:
            df_final[cc] = df_final[cc] / df_final[df_final['nonmissing']][cc].nlargest(1000).sum()
    elif args['intensity_norm'] == 1:
        for cc in all_lbls:
            df_final[cc] = df_final[cc] / df_final[df_final['nonmissing']][cc].median()

    for slbl in ['1', '2']:
        S_len_current = len(all_s_lbls['S%s' % (slbl, )])
        df_final['S%s_mean' % (slbl, )] = df_final[all_s_lbls['S%s' % (slbl, )]].mean(axis=1)
        df_final['S%s_std' % (slbl, )] = np.log2(df_final[all_s_lbls['S%s' % (slbl, )]]).std(axis=1)

    df_final['S1_std'] = df_final['S1_std'].fillna(df_final['S2_std'])
    df_final['S2_std'] = df_final['S2_std'].fillna(df_final['S1_std'])

    df_final['intensity_median'] = df_final[['S1_mean', 'S2_mean']].max(axis=1)
    df_final['iq'] = df_final['nummissing'].astype(str) + pd.qcut(df_final['intensity_median'], 5, labels=range(5)).fillna(0).astype(str)
    df_final['FC'] = np.log2(df_final['S2_mean']/df_final['S1_mean'])

    from scipy.stats import ttest_ind
    idx_to_calc_initial_pval = (df_final['nonmissing'])
    df_final.loc[idx_to_calc_initial_pval, 't-value'] = list(ttest_ind(np.log10(df_final.loc[idx_to_calc_initial_pval, all_s_lbls['S1']].values.astype(float)), np.log10(df_final.loc[idx_to_calc_initial_pval, all_s_lbls['S2']].values.astype(float)), axis=1, nan_policy='omit', equal_var=True)[0])
    df_final['t-value'] = df_final['t-value'].astype(float)
    df_final.loc[(pd.isna(df_final['S2_mean'])) & (~pd.isna(df_final['S1_mean'])), 't-value'] = 10
    df_final.loc[(~pd.isna(df_final['S2_mean'])) & (pd.isna(df_final['S1_mean'])), 't-value'] = -10
    df_final['t-value'] = df_final['t-value'].fillna(0.0)
    df_final['t-value'] = df_final['t-value'].clip(-10, 10)

    FC_max = df_final[df_final['nonmissing']]['FC'].max()
    FC_min = df_final[df_final['nonmissing']]['FC'].min()

    df_final_for_calib = df_final.copy()

    df_final_for_calib = df_final_for_calib[df_final_for_calib['nonmissing']]
    
    df_final_for_calib = df_final_for_calib[~pd.isna(df_final_for_calib['S1_mean'])]
    df_final_for_calib = df_final_for_calib[df_final_for_calib['FC'] <= FC_max/2]
    df_final_for_calib = df_final_for_calib[df_final_for_calib['FC'] >= FC_min/2]
    df_final_for_calib = df_final_for_calib[~pd.isna(df_final_for_calib['S2_mean'])]

    df_final.loc[(pd.isna(df_final['S2_mean'])) & (~pd.isna(df_final['S1_mean'])), 'FC'] = FC_min
    df_final.loc[(~pd.isna(df_final['S2_mean'])) & (pd.isna(df_final['S1_mean'])), 'FC'] = FC_max


    tmp1 = df_final.groupby('protein')['origseq'].count()
    proteins_best50percent_by_num_quantified_peptides = set(tmp1[tmp1 >= np.median(tmp1)].index)

    tmp2 = df_final[df_final['protein'].apply(lambda x: x in proteins_best50percent_by_num_quantified_peptides)].groupby('protein')['FC'].median().abs()
    proteins_best50percent_by_stable_FC = set(tmp2[tmp2 <= np.median(tmp2)].index)
    df_final_for_calib = df_final_for_calib[df_final_for_calib['protein'].apply(lambda x: x in proteins_best50percent_by_stable_FC)]

    tmp = df_final_for_calib['FC']

    try:
        FC_mean, FC_std, covvalue_cor = calibrate_mass(0.05, -tmp.min(), tmp.max(), tmp)
        FC_mean2, FC_std2, covvalue_cor2 = calibrate_mass(0.1, -tmp.min(), tmp.max(), tmp)
        if not np.isinf(covvalue_cor2) and abs(FC_mean2) <= abs(FC_mean) / 10:
            FC_mean = FC_mean2
            FC_std = FC_std2
    except:
        FC_mean, FC_std, covvalue_cor = calibrate_mass(0.3, -tmp.min(), tmp.max(), tmp)



    fc_dict_by_missing = dict()
    for num_missing in range(df_final_for_calib['nummissing'].max()+1):
        xtmp = df_final_for_calib[df_final_for_calib['nummissing'] == num_missing]['FC']
        try:
            try:
                xFC_mean, xFC_std, xcovvalue_cor = calibrate_mass(0.05, -xtmp.min(), xtmp.max(), xtmp)
                xFC_mean2, xFC_std2, xcovvalue_cor2 = calibrate_mass(0.1, -xtmp.min(), xtmp.max(), xtmp)
                if not np.isinf(xcovvalue_cor2) and abs(xFC_mean2) <= abs(xFC_mean) / 10:
                    xFC_mean = xFC_mean2
                    xFC_std = xFC_std2
            except:
                xFC_mean, xFC_std, xcovvalue_cor = calibrate_mass(0.3, -xtmp.min(), xtmp.max(), xtmp)

            fc_dict_by_missing[num_missing] = xFC_std
        except:
            try:
                fc_dict_by_missing[num_missing] = fc_dict_by_missing[num_missing-1]
            except:
                fc_dict_by_missing[num_missing] = False
    for num_missing in list(range(df_final_for_calib['nummissing'].max()+1))[::-1]:
        if fc_dict_by_missing[num_missing] is False:
            try:
                fc_dict_by_missing[num_missing] = fc_dict_by_missing[num_missing+1]
            except:
                fc_dict_by_missing[num_missing] = FC_std

    min_val = 0
    for num_missing in range(df_final_for_calib['nummissing'].max()+1):
        fc_dict_by_missing[num_missing] = max(fc_dict_by_missing[num_missing], min_val)
        min_val = fc_dict_by_missing[num_missing]

    fc_dict_by_missing_base = dict()
    if not args['fold_change_abs']:
        for num_missing in range(df_final_for_calib['nummissing'].max()+1):
            fc_dict_by_missing_base[num_missing] = float(fc_dict_by_missing[num_missing])
            fc_dict_by_missing[num_missing] = fc_dict_by_missing[num_missing] * 2# * fold_change

    if not args['fold_change_abs']:
        fold_change = FC_std * fold_change
    logger.info('Absolute FC threshold for peptides = %.2f +- %.2f', FC_mean, fold_change)

    df_final['decoy'] = df_final['protein'].apply(lambda x: all(z.startswith(decoy_prefix) for z in x.split(';')))

    df_final = df_final.assign(protein=df_final['protein'].str.split(';')).explode('protein').reset_index(drop=False)
    df_final['proteins'] = df_final['protein']
    df_final = df_final.drop(columns=['protein'])

    df_final = df_final.sort_values(by=['nummissing', 'intensity_median'], ascending=(True, False))
    df_final = df_final.drop_duplicates(subset=('origseq', 'proteins'))


    df_final['FC_corrected'] = df_final['FC'] - FC_mean

    if args['protein_shifts']:
        df_shifts = pd.read_table(args['protein_shifts'])
        if 'FC shift' in df_shifts.columns:
            shifts_map = df_shifts.set_index('dbname')['FC shift'].to_dict()
        else:
            df_shifts = df_shifts[df_shifts['identified peptides'] >= 3]
            shifts_map = df_shifts.set_index('dbname')['log2FoldChange(S2/S1) using all peptides'].to_dict()
        df_final['FC_corrected'] = df_final.apply(lambda x: x['FC_corrected'] - shifts_map.get(x['proteins'], 0), axis=1)


        for cc in all_s_lbls['S2']:
            df_final[cc] = df_final[cc] / df_final['proteins'].apply(lambda x: 2**shifts_map.get(x, 0))

    df_final['FC_abs'] = df_final['FC_corrected'].abs()
    df_final = df_final.sort_values(by='FC_abs').reset_index(drop=True)
    df_final['FC_abs'] = df_final['FC_corrected']

    idx_to_calc_initial_pval = (df_final[['nonmissing_S1', 'nonmissing_S2']].min(axis=1) >= 2) & (df_final['nonmissing'])

    df_final.loc[idx_to_calc_initial_pval, 'p-value'] = list(ttest_ind(np.log10(df_final.loc[idx_to_calc_initial_pval, all_s_lbls['S1']].values.astype(float)), np.log10(df_final.loc[idx_to_calc_initial_pval, all_s_lbls['S2']].values.astype(float)), axis=1, nan_policy='omit', equal_var=True)[1])
    df_final['p-value'] = df_final['p-value'].astype(float)


    for cc in all_lbls:
        df_final[cc] = df_final[cc].fillna(df_final[cc].min())

    idx_missing_pval = pd.isna(df_final['p-value'])

    df_final.loc[idx_missing_pval, 'p-value'] = list(ttest_ind(np.log10(df_final.loc[idx_missing_pval, all_s_lbls['S1']].values.astype(float)), np.log10(df_final.loc[idx_missing_pval, all_s_lbls['S2']].values.astype(float)), axis=1, nan_policy='omit', equal_var=True)[1])

    df_final['p-value'] = df_final['p-value'].fillna(1.0)

    df_final['sign'] = df_final.apply(lambda x: np.abs(x['FC_corrected']) >= fc_dict_by_missing[x['nummissing']], axis=1)
    df_final['up'] = df_final['sign'] * (df_final['FC_corrected'] > 0)
    df_final['down'] = df_final['sign'] * (df_final['FC_corrected'] < 0)

    cols = [z for z in df_final.columns.tolist() if not z.startswith('mz_') and not z.startswith('RT_')]
    cols.remove('proteins')
    cols.insert(0, 'proteins')
    df_final = df_final[cols]

    df_final.to_csv(path_or_buf=args['out']+'_quant_peptides.tsv', sep='\t', index=False, float_format="%.4g")

    df_final = df_final.sort_values(by=['nummissing', 'intensity_median'], ascending=(True, False))
    df_final = df_final.drop_duplicates(subset=('origseq', 'proteins'))




    prot_to_peps = defaultdict(str)
    prot_pep_map = defaultdict(set)
    for dbname, pepseq in df_final.sort_values(by='origseq')[['proteins', 'origseq']].values:
        prot_to_peps[dbname] += pepseq
        prot_pep_map[dbname].add(pepseq)


    genes_map = {}
    pep_pos_map = {}
    prot_len_map = {}
    if args['d']:
        for prot, protseq in fasta.read(args['d']):
            if decoy_prefix not in prot:
                try:
                    prot_name = prot.split('|')[1].split(' ')[0]
                except:
                    prot_name = prot.split(' ')[0]
                try:
                    gene_name = prot.split('GN=')[1].split(' ')[0]
                except:
                    gene_name = prot.split(' ')[0]
                genes_map[prot_name] = gene_name

            dbname = prot.split(' ')[0]
            if dbname in prot_pep_map:
                pep_pos_map[dbname] = dict()
                for pepseq in prot_pep_map[dbname]:
                    pep_pos_map[dbname][pepseq] = protseq.find(pepseq) + 1
                prot_len_map[dbname] = len(protseq)


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

    def get_difregmap(x):
        out = ''
        flag = 0
        last_pos = 0
        out += '%d/' % (x['prot_pos'].values[0], )
        for enum_idx, z in enumerate(x['up'].values):
            if not z:
                out += '_'
            else:
                if not flag:
                    out += '(%d)' % (x['prot_pos'].values[enum_idx], )
                    flag = 1
                last_pos = x['prot_pos'].values[enum_idx] + len(x['origseq'].values[enum_idx])-1
                out += '*'
        if flag:
            out_tmp = out[::-1].split('*', 1)
            out = out_tmp[1][::-1] + '*(%d)' % (last_pos, ) + out_tmp[0][::-1]
        out += '/%d/%d' % (x['prot_pos'].values[-1]+len(x['origseq'].values[-1])-1, x['prot_len'].values[-1])
        return out

    if args['d']:
        df_final['prot_pos'] = df_final.apply(lambda x: pep_pos_map[x['proteins']][x['origseq']], axis=1)
        df_final['prot_len'] = df_final.apply(lambda x: prot_len_map[x['proteins']], axis=1)
        df_final = df_final.sort_values(by='prot_pos')
        difregmap = df_final.reset_index(drop=True).groupby('proteins').apply(get_difregmap).to_dict()
    else:
        difregmap = {}


    prots_up = df_final.groupby('proteins')['up'].sum()
    prots_missing = df_final.groupby('proteins')['nummissing'].sum()
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
    all_pvals = calc_sf_all(v_arr, n_arr, p_up, min_matched_peptides)

    df_final = df_final[df_final['nonmissing']]

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

    lbl_FC_to_use = 'log2FoldChange(S2/S1) using all peptides'

    df_out = df_out[~df_out['decoy']]

    df_out['protname'] = df_out['dbname'].apply(lambda x: x.split('|')[1] if '|' in x else x)
    df_out['protein_quant_group'] = df_out['dbname'].apply(lambda x: protein_groups[x])

    if args['d']:
        df_out['gene'] = df_out['protname'].apply(lambda x: genes_map[x])
    else:
        df_out['gene'] = df_out['protname']


    qval_threshold = args['qval']

    df_out = df_out.sort_values(by='score', ascending=False).reset_index(drop=True)
    df_out['FC_pass'] = False
    df_out['FC_pass'] = df_out[lbl_FC_to_use].abs() >= fold_change

    df_out['p-value'] = 10**(-df_out['score'])

    if not args['legacy']:
        from scipy.stats import combine_pvalues
        df_out['p-value1'] = df_out['p-value']
        df_out['p-value2'] = df_out['dbname'].apply(lambda x: prots_pval2.get(x, 1))
        df_out['p-value2'] = df_out['p-value2'].fillna(1.0)
        df_out['p-value'] = df_out[['p-value1', 'p-value2']].apply(lambda x: combine_pvalues([x['p-value1'], x['p-value2']])[1], axis=1)

        df_out['p-value'] = df_out['p-value'].clip(1e-60, 1.0)
        df_out['score'] = -np.log10(df_out['p-value'])
        df_out = df_out.sort_values(by='score', ascending=False).reset_index(drop=True)

    BH_idx = (df_out['identified peptides'] >= min_matched_peptides) & (df_out['FC_pass'])
    df_out['BH_pass'] = False

    from scipy.stats import false_discovery_control

    df_out['p-adj'] = 1.0
    df_out.loc[BH_idx, 'p-adj'] = false_discovery_control(df_out.loc[BH_idx, 'p-value'])
    df_out['BH_pass'] = df_out['p-adj'] <= args['qval']

    df_out.loc[BH_idx, 'FDR_pass'] = df_out.loc[BH_idx, 'score'] >= -np.log10(args['qval'])

    df_out = df_out.drop(columns = {'decoy'})

    df_out['difregmap'] = df_out['dbname'].apply(lambda x: difregmap.get(x, ''))

    if not args['legacy']:
        df_out = df_out[['score', 'p-value', 'p-value1', 'p-value2', 'dbname', 'log2FoldChange(S2/S1)', 'differentially expressed peptides',
                        'identified peptides', 'log2FoldChange(S2/S1) no normalization', 'log2FoldChange(S2/S1) using all peptides',
                        'log2FoldChange(S2/S1) using all peptides and no normalization', 'protname', 'protein_quant_group', 'gene', 'FC_pass', 'FDR_pass', 'BH_pass', 'difregmap', 'p-adj']]# 'BH_threshold']]

    else:
        df_out = df_out[['score', 'p-value', 'dbname', 'log2FoldChange(S2/S1)', 'differentially expressed peptides',
                        'identified peptides', 'log2FoldChange(S2/S1) no normalization', 'log2FoldChange(S2/S1) using all peptides',
                        'log2FoldChange(S2/S1) using all peptides and no normalization', 'protname', 'protein_quant_group', 'gene', 'FC_pass', 'FDR_pass', 'BH_pass', 'difregmap', 'p-adj']]# 'BH_threshold']]

    df_out.to_csv(path_or_buf=args['out']+'_quant_full.tsv', sep='\t', index=False, float_format="%.4g")

    # df_out_f = df_out[(df_out['FDR_pass']) & (df_out['FC_pass'])]
    df_out_f = df_out[(df_out['BH_pass']) & (df_out['FC_pass'])]

    df_out_f.to_csv(path_or_buf=args['out']+'.tsv', sep='\t', index=False, float_format="%.4g")

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
