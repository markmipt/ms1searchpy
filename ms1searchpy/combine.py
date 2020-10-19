from .utils import calc_sf_all, keywithmaxval
from .main import filter_results, final_iteration
import pandas as pd
import numpy as np
from pyteomics import auxiliary as aux
from copy import copy, deepcopy
from collections import defaultdict
from scipy.stats import binom
from scipy.stats import zscore, spearmanr


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
    fdr = args['fdr']

    resdict = dict()

    resdict['qpreds'] = df1['qpreds'].values
    resdict['preds'] = df1['preds'].values
    resdict['seqs'] = df1['peptide'].values
    resdict['ids'] = df1['ids'].values

    mass_diff = resdict['qpreds']
    rt_diff = resdict['qpreds']

    base_out_name = args['out']
    final_iteration(resdict, mass_diff, rt_diff, pept_prot, protsN, base_out_name, prefix, isdecoy, isdecoy_key, escore, fdr, args['nproc'])



    # p1 = set(resdict['seqs'])

    # prots_spc2 = defaultdict(set)
    # for pep, proteins in pept_prot.items():
    #     if pep in p1:
    #         for protein in proteins:
    #             prots_spc2[protein].add(pep)

    # for k in protsN:
    #     if k not in prots_spc2:
    #         prots_spc2[k] = set([])
    # prots_spc = dict((k, len(v)) for k, v in prots_spc2.items())

    # names_arr = np.array(list(prots_spc.keys()))
    # v_arr = np.array(list(prots_spc.values()))
    # n_arr = np.array([protsN[k] for k in prots_spc])

    # top100decoy_score = [prots_spc.get(dprot, 0) for dprot in protsN if isdecoy_key(dprot)]
    # top100decoy_N = [val for key, val in protsN.items() if isdecoy_key(key)]
    # p = np.mean(top100decoy_score) / np.mean(top100decoy_N)
    # print('p=%s' % (np.mean(top100decoy_score) / np.mean(top100decoy_N)))

    # prots_spc = dict()
    # all_pvals = calc_sf_all(v_arr, n_arr, p)
    # for idx, k in enumerate(names_arr):
    #     prots_spc[k] = all_pvals[idx]


    # p_out = defaultdict(float)
    # for i in range(0, 10, 1):
        
        
    #     prots_spc_basic2 = copy(prots_spc)
    #     prots_spc_final = dict()
    #     prots_spc_final2 = dict()


    #     mass_diff = resdict['qpreds']
    #     rt_diff = resdict['qpreds']

    #     mass_koef = i
    #     rtt_koef = mass_koef
    #     e_ind = mass_diff <= mass_koef
    #     resdict2 = filter_results(resdict, e_ind)


    #     prots_spc_final2 = {}
    #     prots_spc_IntensityNorm = defaultdict(list)
    #     prots_spc_IntensityNorm_predicted = defaultdict(list)

    #     for pep, Inorm, Inorm_predicted in zip(resdict2['seqs'], resdict2['IntensityNorm'], resdict2['IntensityNorm_predicted']):
    #         for bprot in pept_prot[pep]:
    #             prots_spc_IntensityNorm[bprot].append(Inorm)
    #             prots_spc_IntensityNorm_predicted[bprot].append(Inorm_predicted)
    #     for bprot in prots_spc_IntensityNorm:
    #         cor_basic = spearmanr(prots_spc_IntensityNorm[bprot], prots_spc_IntensityNorm_predicted[bprot])[0]
    #         prots_spc_final2[bprot] = cor_basic

    #     features_dict = dict()
    #     for pep in set(resdict2['seqs']):
    #         for bprot in pept_prot[pep]:
    #             prot_score = prots_spc_basic2[bprot]
    #             if prot_score > features_dict.get(pep, [-1, ])[-1]:
    #                 features_dict[pep] = (bprot, prot_score)

    #     prots_spc_basic = dict()

    #     p1 = set(resdict2['seqs'])

    #     pep_pid = defaultdict(set)
    #     pid_pep = defaultdict(set)
    #     banned_dict = dict()
    #     for pep, pid in zip(resdict2['seqs'], resdict2['ids']):
    #         pep_pid[pep].add(pid)
    #         pid_pep[pid].add(pep)
    #         if pep in banned_dict:
    #             banned_dict[pep] += 1
    #         else:
    #             banned_dict[pep] = 1

    #     if len(p1):
    #         prots_spc_final = dict()
    #         prots_spc_copy = False
    #         prots_spc2 = False
    #         unstable_prots = set()
    #         p0 = False

    #         names_arr = False
    #         tmp_spc_new = False
    #         decoy_set = False

    #         while 1:
    #             if not prots_spc2:

    #                 best_match_dict = dict()
    #                 n_map_dict = defaultdict(list)
    #                 for k, v in protsN.items():
    #                     n_map_dict[v].append(k)

    #                 decoy_set = set()
    #                 for k in protsN:
    #                     if isdecoy_key(k):
    #                         decoy_set.add(k)
    #                 decoy_set = list(decoy_set)


    #                 prots_spc2 = defaultdict(set)
    #                 for pep, proteins in pept_prot.items():
    #                     if pep in p1:
    #                         for protein in proteins:
    #                             if protein == features_dict[pep][0]:
    #                                 prots_spc2[protein].add(pep)

    #                 for k in protsN:
    #                     if k not in prots_spc2:
    #                         prots_spc2[k] = set([])
    #                 prots_spc2 = dict(prots_spc2)
    #                 unstable_prots = set(prots_spc2.keys())

    #                 top100decoy_N = sum([val for key, val in protsN.items() if isdecoy_key(key)])

    #                 names_arr = np.array(list(prots_spc2.keys()))
    #                 n_arr = np.array([protsN[k] for k in names_arr])

    #                 tmp_spc_new = dict((k, len(v)) for k, v in prots_spc2.items())


    #                 top100decoy_score_tmp = [tmp_spc_new.get(dprot, 0) for dprot in decoy_set]
    #                 top100decoy_score_tmp_sum = float(sum(top100decoy_score_tmp))

    #             tmp_spc = tmp_spc_new
    #             prots_spc = tmp_spc_new
    #             if not prots_spc_copy:
    #                 prots_spc_copy = deepcopy(prots_spc)

    #             for idx, v in enumerate(decoy_set):
    #                 if v in unstable_prots:
    #                     top100decoy_score_tmp_sum -= top100decoy_score_tmp[idx]
    #                     top100decoy_score_tmp[idx] = prots_spc.get(v, 0)
    #                     top100decoy_score_tmp_sum += top100decoy_score_tmp[idx]
    #             p = float(sum(top100decoy_score_tmp)) / top100decoy_N
    #             p = top100decoy_score_tmp_sum / top100decoy_N
    #             if not p0:
    #                 p0 = float(p)
    #                 print(p0)

    #             n_change = set(protsN[k] for k in unstable_prots)
    #             for n_val in n_change:
    #                 for k in n_map_dict[n_val]:
    #                     v = prots_spc[k]
    #                     if n_val not in best_match_dict or v > prots_spc[best_match_dict[n_val]]:
    #                         best_match_dict[n_val] = k
    #             n_arr_small = []
    #             names_arr_small = []
    #             v_arr_small = []
    #             for k, v in best_match_dict.items():
    #                 n_arr_small.append(k)
    #                 names_arr_small.append(v)
    #                 v_arr_small.append(prots_spc[v])

    #             prots_spc_basic = dict()
    #             all_pvals = calc_sf_all(v_arr_small, n_arr_small, p)
    #             for idx, k in enumerate(names_arr_small):
    #                 prots_spc_basic[k] = all_pvals[idx]

    #             best_prot = keywithmaxval(prots_spc_basic)

    #             best_score = prots_spc_basic[best_prot]
    #             unstable_prots = set()
    #             if best_prot not in prots_spc_final:
    #                 prots_spc_final[best_prot] = best_score
    #                 banned_pids = set()
    #                 for pep in prots_spc2[best_prot]:
    #                     for pid in pep_pid[pep]:
    #                         banned_pids.add(pid)
    #                 for pid in banned_pids:
    #                     for pep in pid_pep[pid]:
    #                         banned_dict[pep] -= 1
    #                         if banned_dict[pep] == 0:
    #                             best_prot_val = features_dict[pep][0]
    #                             for bprot in pept_prot[pep]:
    #                                 if bprot == best_prot_val:
    #                                     tmp_spc_new[bprot] -= 1
    #                                     unstable_prots.add(bprot)
    #             else:

    #                 v_arr = np.array([prots_spc[k] for k in names_arr])
    #                 all_pvals = calc_sf_all(v_arr, n_arr, p)
    #                 for idx, k in enumerate(names_arr):
    #                     prots_spc_basic[k] = all_pvals[idx]

    #                 for k, v in prots_spc_basic.items():
    #                     if k not in prots_spc_final:
    #                         prots_spc_final[k] = v

    #                 break

    #             try:
    #                 prot_fdr = aux.fdr(prots_spc_final.items(), is_decoy=isdecoy)
    #             except ZeroDivisionError:
    #                 prot_fdr = 100.0
    #             if prot_fdr >= 12.5 * fdr:

    #                 v_arr = np.array([prots_spc[k] for k in names_arr])
    #                 all_pvals = calc_sf_all(v_arr, n_arr, p)
    #                 for idx, k in enumerate(names_arr):
    #                     prots_spc_basic[k] = all_pvals[idx]

    #                 for k, v in prots_spc_basic.items():
    #                     if k not in prots_spc_final:
    #                         prots_spc_final[k] = v
    #                 break

    #     keys_array = prots_spc_final.keys()
    #     values_array = np.array(list(prots_spc_final.values()))
    #     # values_array = zscore(values_array)
        
    #     val_median = np.median(values_array)
        
    #     for k, val, zval in zip(keys_array, values_array, zscore(values_array)):
    #         prots_spc_final[k] = val
    #         p_out[k] += val

    #     filtered_prots = aux.filter(prots_spc_final.items(), fdr=fdr, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1, full_output=True, correction=0)
    #     print(mass_koef, len(filtered_prots))
        
    #     filtered_prots = aux.filter(p_out.items(), fdr=fdr, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1, full_output=True, correction=0)
    #     print('z', len(filtered_prots))

if __name__ == '__main__':
    run()
