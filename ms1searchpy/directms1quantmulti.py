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

from . import directms1quant
from os import path, listdir, makedirs
from copy import copy
import matplotlib.pyplot as plt
import seaborn as sb

logger = logging.getLogger(__name__)

def run():
    parser = argparse.ArgumentParser(
        description='Prepare LFQ protein table for project with multiple conditions (time-series,\
        concentration, thermal profiling, etc.',
        epilog='''

    Example usage
    -------------
    $ directms1quantmulti -S1 sample1_1_proteins_full.tsv sample1_n_proteins_full.tsv -S2 sample2_1_proteins_full.tsv sample2_n_proteins_full.tsv
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-pdir', help='path to project folder with search results', required=True)
    # parser.add_argument('-d', '-db', help='path to uniprot fasta file used for ms1searchpy', required=True)
    parser.add_argument('-samples', help='tsv table with sample details', required=True)
    parser.add_argument('-out', help='name of DirectMS1quant output files', default='DQmulti')
    parser.add_argument('-norm', help='LFQ normalization: (1) using median of CONTROL group, default; (2) using median of all groups', default=1, type=int)
    parser.add_argument('-proteins_for_figure', help='path to proteins for figure plotting', default='', type=str)
    parser.add_argument('-figdir', help='path to output folder for figures', default='')
    # parser.add_argument('-min_samples', help='minimum number of non-missing peptide intensities for DirectMS1Quant pairwise comparisons. 0 means 50%% of input files', default=1)
    # parser.add_argument('-pep_min_non_missing_samples', help='minimum fraction of files with non missing values for peptide', default=0.5, type=float)
    parser.add_argument('-max_missing', help='maximum fraction of missing values. default = 0.5', default=0.5, type=float)
    # parser.add_argument('-min_signif_for_pept', help='minimum number of pairwise DE results where peptide should be significant', default=999, type=int)
    parser.add_argument('-prefix', help='Decoy prefix. Default DECOY_', default='DECOY_', type=str)
    parser.add_argument('-start_stage', help='Can be 1, 2 or 3 to skip any stage which were already done', default=1, type=int)
    args = vars(parser.parse_args())
    logging.basicConfig(format='%(levelname)9s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=logging.INFO)

    process_files(args, logger)


def process_files(args, logger):

    infolder = args['pdir']
    ms1folder = args['pdir']
    # path_to_fasta = args['d']

    df1 = pd.read_table(args['samples'], dtype={'group': str, 'condition': str, 'vs': str})
    df1['sample'] = df1['group']
    if 'condition' not in df1.columns:
        df1['condition'] = ''
        
    if 'replicate' not in df1.columns:
        df1['replicate'] = 1
        
    if 'BatchMS' not in df1.columns:
        df1['BatchMS'] = 1

    if 'vs' not in df1.columns:
        df1['vs'] = df1['condition']
    df1['vs'] = df1['vs'].fillna(df1['condition'])

    df1_filenames = set(df1['File Name'])

    f_dict_map = {}

    for fname, sname, ttime in df1[['File Name', 'sample', 'condition']].values:
        f_dict_map[fname] = (sname, ttime)
        
        
    control_label = df1['group'].values[0]


    all_conditions = df1[df1['group']!=control_label].set_index(['group', 'condition'])['vs'].to_dict()

    outlabel = args['out']

    replace_label = '_proteins_full.tsv'
    decoy_prefix = args['prefix']

    names_map = {}

    args_local = {'S1': []}
    args_local['allowed_peptides'] = ''
    args_local['allowed_proteins'] = ''
    if infolder:
        for z in listdir(infolder):
            if z.endswith(replace_label):
                zname = z.split('.')[0]
                if zname in df1_filenames:
                    args_local['S1'].append(path.join(infolder, z))
                    names_map[args_local['S1'][-1]] = zname
    else:
        for z in df1_filenames:
            zname = z
            args_local['S1'].append(z+'.features_proteins_full.tsv')
            names_map[args_local['S1'][-1]] = zname

    all_s_lbls = {}

    allowed_prots = set()
    allowed_prots_all = set()
    allowed_peptides = set()

    cnt0 = Counter()

    if args['start_stage'] > 2:

        sample_num = 'S1'

        all_s_lbls[sample_num] = []

        for z in args_local[sample_num]:
            label = sample_num + '_' + z.replace(replace_label, '')
            all_s_lbls[sample_num].append(label)

        all_lbls = all_s_lbls['S1']
        out_name = path.join(ms1folder, '%s_peptide_LFQ.tsv' % (outlabel, ))
        df_final = pd.read_table(out_name)

        logger.info('Skipping Stage 2: Prepare full LFQ peptide table...')

    else:
        logger.info('Starting Stage 2: Prepare full LFQ peptide table...')
                
        sample_num = 'S1'

        all_s_lbls[sample_num] = []

        for z in args_local[sample_num]:
            label = sample_num + '_' + z.replace(replace_label, '')
            all_s_lbls[sample_num].append(label)

            df0 = pd.read_table(z.replace('_proteins_full.tsv', '_proteins.tsv'), usecols=['dbname', ])
            allowed_prots.update(df0['dbname'])
            allowed_prots.update([decoy_prefix + z for z in df0['dbname'].values])

            df0 = pd.read_table(z.replace('_proteins_full.tsv', '_PFMs_ML.tsv'), usecols=['seqs', 'qpreds'])#, 'ch', 'im'])
            df0 = df0[df0['qpreds'] <= 10]
            df0['seqs'] = df0['seqs']# + df0['ch'].astype(str) + df0['im'].astype(str)
            allowed_peptides.update(df0['seqs'])
            cnt0.update(df0['seqs'])

        custom_min_samples = 1

        allowed_peptides = set()
        for k, v in cnt0.items():
            if v >= custom_min_samples:
                allowed_peptides.add(k)


        logger.info('Total number of TARGET protein GROUPS: %d', len(allowed_prots) / 2)

        sample_num = 'S1'

        if args_local.get(sample_num, 0):
            for z in args_local[sample_num]:

                df3 = pd.read_table(z.replace(replace_label, '_PFMs.tsv'), usecols=['sequence', 'proteins', 'charge', 'ion_mobility'])

                df3['tmpseq'] = df3['sequence']# + df3['charge'].astype(str) + df3['ion_mobility'].astype(str)
                df3 = df3[df3['tmpseq'].apply(lambda x: x in allowed_peptides)]

                df3_tmp = df3[df3['proteins'].apply(lambda x: any(z in allowed_prots for z in x.split(';')))]
                for dbnames in set(df3_tmp['proteins'].values):
                    for dbname in dbnames.split(';'):
                        allowed_prots_all.add(dbname)

        df_final = directms1quant.get_df_final(args_local, replace_label, allowed_peptides, allowed_prots_all, pep_RT=False, RT_threshold=False)

        logger.info('Total number of peptide sequences used in quantitation: %d', len(set(df_final['origseq'])))

        all_lbls = all_s_lbls['S1']



        rt_ccols = [z for z in df_final.columns.tolist() if z.startswith('RT_')]
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
            print('non-missing peptides for %s: %d' % (cc, non_missing_peptides))
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
            koef2_base = directms1quant.weighted_quantiles_interpolate(ar_ratio, np.sqrt(ar2), 0.5)
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
                    koef2 = directms1quant.weighted_quantiles_interpolate(ar_ratio, np.sqrt(ar2), 0.5)
                    norm_dict[(cc, RT_int)] = koef2



        cols = [z for z in df_final.columns.tolist() if not z.startswith('mz_') and not z.startswith('RT_')]
        df_final = df_final[cols]

        df_final = df_final.set_index('peptide')

        df_final_copy = df_final.copy()

        df_final = df_final_copy.copy()

        max_missing = len(all_lbls) - custom_min_samples

        df_final['nummissing'] = df_final.isna().sum(axis=1)
        df_final['nonmissing'] = df_final['nummissing'] <= max_missing

        df_final = df_final[df_final['nonmissing']]
        logger.info('Total number of PFMs: %d', len(df_final))

        out_name = path.join(ms1folder, '%s_peptide_LFQ.tsv' % (outlabel, ))
        df_final.to_csv(out_name, sep='\t', index=False)






    if args['start_stage'] <= 3:

        logger.info('Starting Stage 3: Prepare full LFQ protein table...')


        all_lbls_by_batch = defaultdict(list)
        all_filenames_control = set(df1[df1['group'] == control_label]['File Name'])
        all_lbls_control = set()

        bdict = df1.set_index('File Name')['BatchMS'].to_dict()

        for cc in all_lbls:
            try:
                origfn = names_map[path.join(infolder, cc.split('/')[-1] + replace_label)]
            except:
                try:
                    origfn = names_map[cc + replace_label]
                except:
                    origfn = cc.split('_', 1)[-1] + replace_label
            try:
                all_lbls_by_batch[bdict[origfn]].append(cc)
            except:
                all_lbls_by_batch[bdict[origfn.replace('.features' + replace_label, '')]].append(cc)

            if origfn in all_filenames_control:
                all_lbls_control.add(cc)
            

        for lbl_name, small_lbls in all_lbls_by_batch.items():
            lbl_len = len(small_lbls)
            idx_to_keep = df_final[small_lbls].isna().sum(axis=1) <= args['max_missing'] * lbl_len
            df_final = df_final[idx_to_keep]

        # for cc in all_lbls:
        #     df_final[cc] = df_final[cc] / df_final[cc].nlargest(1000).sum()


        for cc in all_lbls:
            df_final[cc] = df_final.apply(lambda x: x[cc] / norm_dict[(cc, x['q_RT'])], axis=1)



        # idx_to_keep = df_final['origseq'].apply(lambda x: x in allowed_peptides_base_only_sequences)
        # df_final = df_final[idx_to_keep]

        
        for small_lbls in all_lbls_by_batch.values():
            m = df_final[small_lbls].min(axis=1)
            for col in df_final[small_lbls]:
                df_final.loc[:, col] = df_final.loc[:, col].fillna(m)
                
        for lbl_key, small_lbls in all_lbls_by_batch.items():

            if args['norm'] == 2:
                m = df_final[small_lbls].median(axis=1)
            elif args['norm'] == 1:
                small_lbls_control = [z for z in small_lbls if z in all_lbls_control]
                m = df_final[small_lbls_control].median(axis=1)
            else:
                logger.error('norm option must be 1 or 2!')
                return -1
            
            
            for col in df_final[small_lbls]:
                df_final.loc[:, col] = np.log2(df_final.loc[:, col] / m)
                
        df_final = df_final.fillna(-10)


        df_final = df_final.assign(protein=df_final['protein'].str.split(';')).explode('protein').reset_index(drop=True)
        df_final['proteins'] = df_final['protein']
        df_final = df_final.drop(columns=['protein'])
        cols = [z for z in df_final.columns.tolist() if not z.startswith('mz_') and not z.startswith('RT_')]
        cols.remove('proteins')
        cols.insert(0, 'proteins')
        df_final = df_final[cols]

        # idx_to_keep = df_final.apply(lambda x: (x['origseq'], x['proteins']) in allowed_peptides_base, axis=1)
        # df_final = df_final[idx_to_keep].copy()



        # idx_to_keep = df_final['origseq'].apply(lambda x: x in allowed_peptides_up)
        # idx_to_keep = df_final.apply(lambda x: (x['origseq'], x['proteins']) in allowed_peptides_up, axis=1)
        # df_final_accurate = df_final[idx_to_keep].copy()
        # df_final_accurate = df_final_accurate[df_final_accurate.groupby('proteins')['origseq'].transform('count') > 1]
        # accurate_proteins = set(df_final_accurate['proteins'])
        # df_final = pd.concat([df_final[df_final['proteins'].apply(lambda x: x not in accurate_proteins)], df_final_accurate])


        df_final['S1_mean'] = df_final[all_lbls].median(axis=1)
        df_final['intensity_median'] = df_final['S1_mean']

        df_final = df_final.sort_values(by=['nummissing', 'intensity_median'], ascending=(True, False))
        df_final = df_final.drop_duplicates(subset=('origseq', 'proteins'))

        df_final = df_final[df_final['proteins'].apply(lambda x: not x.startswith('DECOY_'))]

        df_final = df_final[df_final.groupby('proteins')['S1_mean'].transform('count') > 1]

        origfnames = []
        for cc in all_lbls:
            if infolder:
                origfn = cc.split('/')[-1].split('.')[0]
            else:
                origfn = cc.replace('S1_', '').split('.')[0]
            origfnames.append(origfn)
            
        dft = pd.DataFrame.from_dict([df_final.groupby('proteins')[cc].median().to_dict() for cc in all_lbls])
        dft2 = pd.DataFrame.from_dict([df_final.groupby('origseq')[cc].median().to_dict() for cc in all_lbls])

        dft['File Name'] = origfnames
        dft2['File Name'] = origfnames

        # print(df1['File Name'].values)
        # print(dft['File Name'].values)

        df1x = pd.merge(df1, dft2, on='File Name', how='left')
        df1 = pd.merge(df1, dft, on='File Name', how='left')

        df1['sample+condition'] = df1['sample'].apply(lambda x: x + ' ') + df1['condition']
        df1x['sample+condition'] = df1x['sample'].apply(lambda x: x + ' ') + df1x['condition']


        out_name = path.join(ms1folder, '%s_proteins_LFQ.tsv' % (outlabel, ))
        out_namex = path.join(ms1folder, '%s_peptide_LFQ_processed.tsv' % (outlabel, ))
        df1.to_csv(out_name, sep='\t', index=False)
        df1x.to_csv(out_namex, sep='\t', index=False)

    else:

        logger.info('Skipping Stage 3: Prepare full LFQ protein table...')
        out_name = path.join(ms1folder, '%s_proteins_LFQ.tsv' % (outlabel, ))
        df1 = pd.read_table(out_name)

    if args['start_stage'] <= 4:

        logger.info('Starting Stage 4: Plot figures for selected proteins...')

        warning_msg_1 = 'Provide file with proteins selected for figures.\
            It should be a tsv table with dbname column. For example, it could be the standard output table of directms1quant\
            Note, that protein names should be provided in uniprot format. For example, sp|P28838|AMPL_HUMAN'

        if not args['proteins_for_figure']:
            logger.warning(warning_msg_1)
        else:

            df_prots = pd.read_table(args['proteins_for_figure'])
            if 'dbname' not in df_prots.columns:
                logger.warning('dbname column is missing in proteins file')
            else:

                allowed_proteins_for_figures = set(df_prots['dbname'])

                for k in allowed_proteins_for_figures:
                    if k not in df1.columns:
                        logger.info('Protein %s was not quantified', k)
                    else:


                        try:
                            gname = k.split('|')[1]
                        except:
                            gname = k


                        plt.figure()
                        prot_name = k

                        all_one_char_colors = ['m', 'r', 'c', 'g', 'b', 'y', ]

                        color_custom = {}
                        for s_idx, s_val in enumerate(set(df1['sample'])):
                            s_idx_sm = s_idx
                            while s_idx_sm >= 6:
                                s_idx_sm -= 6
                            color_custom[s_val] = all_one_char_colors[s_idx_sm]

                        my_pal = dict()
                        for x_val, s_val in df1[['sample+condition', 'sample']].values:
                            my_pal[x_val] = color_custom[s_val]

                        ax = sb.boxplot(data=df1, x='sample+condition', hue = 'sample+condition', y = prot_name, palette=my_pal, legend=False)
                        xlabels_custom = [z for z in ax.get_xticklabels()]
                        ax.set_xticks(ax.get_xticks())
                        ax.set_xticklabels(xlabels_custom, rotation=90, size=12)
                        plt.title(ax.get_ylabel())
                        ax.set_ylabel('log2 LFQ', size=14)
                        plt.subplots_adjust()
                        plt.tight_layout()
                        if args['figdir']:
                            out_figdir = args['figdir']
                            if not path.isdir(out_figdir):
                                makedirs(out_figdir)
                        else:
                            out_figdir = infolder
                        plt.savefig(path.join(out_figdir, '%s_%s.png' % (args['out'], gname, )))


if __name__ == '__main__':
    run()
