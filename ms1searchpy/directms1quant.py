from __future__ import division
import argparse
import pandas as pd
import numpy as np
from scipy.stats import binom, ttest_ind
import logging
from pyteomics import fasta

def calc_sf_all(v, n, p):
    sf_values = -np.log10(binom.sf(v-1, n, p))
    sf_values[v <= 1] = 0
    sf_values[np.isinf(sf_values)] = 20
    sf_values[n == 0] = 0
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
    parser.add_argument('-fold_change', help='FC threshold standard deviations', default=3.0, type=float)
    parser.add_argument('-fold_change_abs', help='Use absolute log2 scale FC threshold instead of standard deviations', action='store_true')
    parser.add_argument('-qval', help='qvalue threshold', default=0.05, type=float)
    parser.add_argument('-intensity_norm', help='Intensity normalization: 0-none, 1-median', default=1, type=int)
    parser.add_argument('-all_proteins', help='use all proteins instead of FDR controlled', action='store_true')
    parser.add_argument('-all_pfms', help='use all PFMs instead of ML controlled', action='store_true')
    parser.add_argument('-output_peptides', help='Add output table with peptides', action='store_true')
    parser.add_argument('-d', '-db', help='path to uniprot fasta file for gene annotation')
    args = vars(parser.parse_args())
    logging.basicConfig(format='%(levelname)9s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=logging.INFO)
    logger = logging.getLogger(__name__)


    replace_label = '_proteins_full.tsv'

    fold_change = float(args['fold_change'])

    df_final = False

    all_s_lbls = {}


    allowed_prots = set()
    allowed_prots_all = set()
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

                if not args['all_pfms']:

                    df0 = df0[df0['qpreds'] <= 10]
                allowed_peptides.update(df0['seqs'])

    logger.info('Total number of TARGET protein GROUPS: %d', len(allowed_prots) / 2)

    for i in range(1, 3, 1):
        sample_num = 'S%d' % (i, )
        if args.get(sample_num, 0):
            for z in args[sample_num]:
                df3 = pd.read_table(z.replace(replace_label, '_PFMs.tsv'))
                df3 = df3[df3['sequence'].apply(lambda x: x in allowed_peptides)]

                df3_tmp = df3[df3['proteins'].apply(lambda x: any(z in allowed_prots for z in x.split(';')))]
                for dbnames in set(df3_tmp['proteins'].values):
                    for dbname in dbnames.split(';'):
                        allowed_prots_all.add(dbname)


    for i in range(1, 3, 1):
        sample_num = 'S%d' % (i, )
        if args.get(sample_num, 0):
            all_s_lbls[sample_num] = []
            for z in args[sample_num]:
                label = z.replace(replace_label, '')
                all_s_lbls[sample_num].append(label)
                df3 = pd.read_table(z.replace(replace_label, '_PFMs.tsv'))

                df3 = df3[df3['sequence'].apply(lambda x: x in allowed_peptides)]

                # allowed_prots2 = set()
                # df3_tmp = df3[df3['proteins'].apply(lambda x: any(z in allowed_prots for z in x.split(';')))]
                # for dbnames in set(df3_tmp['proteins'].values):
                #     for dbname in dbnames.split(';'):
                #         allowed_prots2.add(dbname)
                # print('!', len(allowed_prots), len(allowed_prots2))

                # df3 = df3[df3['proteins'].apply(lambda x: any(z in allowed_prots for z in x.split(';')))]

                df3 = df3[df3['proteins'].apply(lambda x: any(z in allowed_prots_all for z in x.split(';')))]
                df3['proteins'] = df3['proteins'].apply(lambda x: ';'.join([z for z in x.split(';') if z in allowed_prots_all]))
                ### df3['proteins'] = df3['proteins'].apply(lambda x: ';'.join([z for z in x.split(';') if z in allowed_prots]))

                df3['origseq'] = df3['sequence']
                df3['sequence'] = df3['sequence'] + df3['charge'].astype(int).astype(str) + df3['ion_mobility'].astype(str)

                df3 = df3.sort_values(by='Intensity', ascending=False)

                df3 = df3.drop_duplicates(subset='sequence')

                df3[label] = df3['Intensity']
                df3['protein'] = df3['proteins']
                df3['peptide'] = df3['sequence']
                df3 = df3[['origseq', 'peptide', 'protein', label]]

                if df_final is False:
                    df_final = df3.reset_index(drop=True)
                else:
                    df_final = df_final.reset_index(drop=True).merge(df3.reset_index(drop=True), on='peptide', how='outer')
                    df_final.protein_x.fillna(value=df_final.protein_y, inplace=True)
                    df_final.origseq_x.fillna(value=df_final.origseq_y, inplace=True)
                    df_final['protein'] = df_final['protein_x']
                    df_final['origseq'] = df_final['origseq_x']

                    df_final = df_final.drop(columns=['protein_x', 'protein_y'])
                    df_final = df_final.drop(columns=['origseq_x', 'origseq_y'])


    logger.info('Total number of peptide sequences used in quantitation: %d', len(set(df_final['origseq'])))
    # print('Total number of proteins used in quantitation %d' % (len(allowed_prots_all), ))


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

    custom_min_samples = int(args['min_samples'])
    if custom_min_samples == 0:
        custom_min_samples = int(len(all_lbls)/2)

    df_final = df_final_copy.copy()

    max_missing = len(all_lbls) - custom_min_samples

    logger.info('Allowed max number of missing values: %d', max_missing)

    df_final['nummissing'] = df_final.isna().sum(axis=1)
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


    if args['intensity_norm'] == 1:
        for cc in all_lbls:
            # print(cc, df_final[cc].median())
            df_final[cc] = df_final[cc] / df_final[cc].median()


    df_final['S2_mean'] = df_final[all_s_lbls['S2']].mean(axis=1)
    df_final['S1_mean'] = df_final[all_s_lbls['S1']].mean(axis=1)

    for cc in all_lbls:
        df_final[cc] = df_final[cc].fillna(df_final[cc].min())

    df_final['p-value'] = list(ttest_ind(df_final[all_s_lbls['S1']].values.astype(float), df_final[all_s_lbls['S2']].values.astype(float), axis=1, nan_policy='omit', equal_var=True)[1])
    df_final['p-value'] = df_final['p-value'].astype(float)
    df_final['p-value'] = df_final['p-value'].fillna(1.0)

    p_val_threshold = 0.05

    df_final['sign'] = df_final['p-value'] <= p_val_threshold

    df_final['intensity_median'] = df_final[all_s_lbls['S1'] + all_s_lbls['S2']].median(axis=1)

    df_final['FC'] = np.log2(df_final['S2_mean']/df_final['S1_mean'])

    df_final_for_calib = df_final.copy()
    df_final_for_calib = df_final_for_calib[~pd.isna(df_final_for_calib['S1_mean'])]
    df_final_for_calib = df_final_for_calib[~pd.isna(df_final_for_calib['S2_mean'])]
    df_final_for_calib = df_final_for_calib[~df_final_for_calib['sign']]

    FC_max = df_final['FC'].max()
    FC_min = df_final['FC'].min()

    df_final.loc[(pd.isna(df_final['S2_mean'])) & (~pd.isna(df_final['S1_mean'])), 'FC'] = FC_min
    df_final.loc[(~pd.isna(df_final['S2_mean'])) & (pd.isna(df_final['S1_mean'])), 'FC'] = FC_max

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

    try:
        FC_mean, FC_std, covvalue_cor = calibrate_mass(0.05, -df_final_for_calib['FC'].min(), df_final_for_calib['FC'].max(), df_final_for_calib['FC'])
    except:
        FC_mean, FC_std, covvalue_cor = calibrate_mass(0.1, -df_final_for_calib['FC'].min(), df_final_for_calib['FC'].max(), df_final_for_calib['FC'])
    # print('df_final_FC', FC_mean, FC_std)

    # FC_l = FC_mean-fold_change
    # FC_r = FC_mean+fold_change

    if not args['fold_change_abs']:
        fold_change = FC_std * fold_change
    logger.info('Absolute FC threshold = %.2f', fold_change)
    FC_l = -fold_change
    FC_r = fold_change


    df_final['up'] = df_final['sign'] * (df_final['FC'] >= FC_r)
    df_final['down'] = df_final['sign'] * (df_final['FC'] <= FC_l)

    df_final = df_final.sort_values(by=['nummissing', 'intensity_median'], ascending=(True, False))
    df_final = df_final.drop_duplicates(subset=('origseq', 'proteins'))

    up_dict = df_final.groupby('proteins')['up'].sum().to_dict()
    down_dict = df_final.groupby('proteins')['down'].sum().to_dict()

    ####### !!!!!!! #######
    df_final['up'] = df_final.apply(lambda x: x['up'] if up_dict.get(x['proteins'], 0) >= down_dict.get(x['proteins'], 0) else x['down'], axis=1)
    protsN = df_final.groupby('proteins')['up'].count().to_dict()

    prots_up = df_final.groupby('proteins')['up'].sum()

    N_decoy_total = df_final['decoy'].sum()

    upreg_decoy_total = df_final[df_final['decoy']]['up'].sum()

    p_up = upreg_decoy_total / N_decoy_total

    names_arr = np.array(list(protsN.keys()))

    logger.info('Total number of proteins used in quantitation: %d', sum(not z.startswith('DECOY_') for z in names_arr))
    logger.info('Total number of peptides: %d', len(df_final))
    logger.info('Total number of decoy peptides: %d', N_decoy_total)
    logger.info('Total number of significantly changed decoy peptides: %d', upreg_decoy_total)
    logger.info('Probability of random peptide to be significantly changed: %.3f', p_up)
    # print(N_decoy_total, upreg_decoy_total, p_up)

    if args['output_peptides']:
        df_final.to_csv(path_or_buf=args['out']+'_quant_peptides.tsv', sep='\t', index=False)

    v_arr = np.array(list(prots_up.get(k, 0) for k in names_arr))
    n_arr = np.array(list(protsN.get(k, 0) for k in names_arr))

    all_pvals = calc_sf_all(v_arr, n_arr, p_up)

    total_set = set()
    total_set_genes = set()

    FC_up_dict_basic = df_final.groupby('proteins')['FC'].median().to_dict()
    FC_up_dict_raw_basic = df_final.groupby('proteins')['FC_raw'].median().to_dict()

    df_final = df_final[df_final['up']>0]

    df_final['bestmissing'] = df_final.groupby('proteins')['nummissing'].transform('min')

    FC_up_dict = df_final[df_final['bestmissing']==df_final['nummissing']].groupby('proteins')['FC'].median().to_dict()
    FC_up_dict_raw = df_final[df_final['bestmissing']==df_final['nummissing']].groupby('proteins')['FC_raw'].median().to_dict()

    # FC_up_dict = df_final.groupby('proteins')['FC'].median().to_dict()


    df_out = pd.DataFrame()
    df_out['score'] = all_pvals
    df_out['dbname'] = names_arr

    df_out['FC'] = df_out['dbname'].apply(lambda x: FC_up_dict.get(x))
    df_out['FC_raw'] = df_out['dbname'].apply(lambda x: FC_up_dict_raw.get(x))

    df_out.loc[pd.isna(df_out['FC']), 'FC'] = df_out.loc[pd.isna(df_out['FC']), 'dbname'].apply(lambda x: FC_up_dict_basic.get(x))
    df_out.loc[pd.isna(df_out['FC_raw']), 'FC_raw'] = df_out.loc[pd.isna(df_out['FC_raw']), 'dbname'].apply(lambda x: FC_up_dict_raw_basic.get(x))

    df_out['v_arr'] = v_arr
    df_out['n_arr'] = n_arr

    df_out['decoy'] = df_out['dbname'].str.startswith('DECOY_')

    df_out = df_out[~df_out['decoy']]

    df_out['FC_pass'] = (df_out['FC'].abs() >= fold_change) & (df_out['v_arr'] > 0)

    df_out_BH_multiplier = df_out['FC_pass'].sum()

    qval_threshold = args['qval']

    df_out['p-value'] = 1.0
    df_out['BH_pass'] = False

    df_out = df_out.sort_values(by='score', ascending=False)
    df_out.loc[df_out['FC_pass'], 'BH_threshold'] = -np.log10(df_out.loc[df_out['FC_pass'], 'score'].rank(ascending=False, method='max') * qval_threshold / df_out_BH_multiplier)
    df_out.loc[df_out['FC_pass'], 'BH_pass'] = df_out.loc[df_out['FC_pass'], 'score'] > df_out.loc[df_out['FC_pass'], 'BH_threshold']
    df_out.loc[df_out['FC_pass'], 'p-value'] = 10**(-df_out.loc[df_out['FC_pass'], 'score'])
    score_threshold = df_out[df_out['BH_pass']]['score'].min()
    df_out.loc[df_out['FC_pass'], 'BH_pass'] = df_out.loc[df_out['FC_pass'], 'score'] >= score_threshold

    df_out.to_csv(path_or_buf=args['out']+'_quant_full.tsv', sep='\t', index=False)

    df_out_f = df_out[(df_out['BH_pass']) & (df_out['FC_pass'])]

    df_out_f.to_csv(path_or_buf=args['out']+'.tsv', sep='\t', index=False)

    genes_map = {}
    if args['d']:
        for prot, protseq in fasta.read(args['d']):
            try:
                prot_name = prot.split('|')[1]
            except:
                prot_name = prot
            try:
                gene_name = prot.split('GN=')[1].split(' ')[0]
            except:
                gene_name = prot
            genes_map[prot_name] = gene_name


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
