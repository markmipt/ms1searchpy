import os
import pickle
from . import utils
import numpy as np
from scipy.stats import scoreatpercentile
from scipy.optimize import curve_fit
from scipy import exp
import operator
from copy import copy, deepcopy
from collections import defaultdict, Counter
import re
from pyteomics import parser, mass, fasta, auxiliary as aux, achrom
try:
    from pyteomics import cmass
except ImportError:
    cmass = mass
import subprocess
from sklearn import linear_model
import tempfile
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Queue, Process, cpu_count
from itertools import chain
try:
    import seaborn
    seaborn.set(rc={'axes.facecolor':'#ffffff'})
    seaborn.set_style('whitegrid')
except:
    pass

from .utils import calc_sf_all, recalc_spc
# from .utils import recalc_spc

import lightgbm as lgb
import pandas as pd
from sklearn.model_selection import train_test_split
from scipy.stats import zscore, spearmanr






# from __future__ import division
import pandas as pd
from pyteomics import pepxml, achrom, auxiliary as aux, mass, fasta, mzid, parser
from pyteomics import electrochem
import numpy as np
import random
SEED = 42
# from catboost import CatBoostClassifier, CatBoostRegressor
from sklearn.model_selection import train_test_split
from os import path, mkdir
from collections import Counter, defaultdict
import warnings
import pylab as plt
warnings.formatwarning = lambda msg, *args, **kw: str(msg) + '\n'

import pandas as pd
# from catboost import CatBoostClassifier, CatBoostRegressor
from sklearn.model_selection import train_test_split, KFold
import os
from collections import Counter, defaultdict
from scipy.stats import scoreatpercentile
from sklearn.isotonic import IsotonicRegression
import warnings
import numpy as np

import matplotlib
import numpy
import pandas
import random
import sklearn
# import xgboost
import matplotlib.pyplot as plt

from sklearn import (
    feature_extraction, feature_selection, decomposition, linear_model,
    model_selection, metrics, svm
)

import scipy
# import eli5
from scipy.stats import rankdata
import pickle
from copy import deepcopy
import csv

from scipy.stats import rankdata
import lightgbm as lgb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import chain
import time as timemodule
import ast
from sklearn import metrics




from scipy.stats import binom
def calc_sf_all_2(v, n, p):
    sf_values = -np.log10(binom.sf(v, n, p))
    sf_values[np.isinf(sf_values)] = max(sf_values[~np.isinf(sf_values)])
    return sf_values

def agg_func(x):
    return np.average(x.im, weights=x.Is)

        
SEED = 50

def get_cat_model(df, hyperparameters, feature_columns, train, test):
    
#     train = df[df['era'].isin(train_eras)]
#     test = df[df['era'].isin(test_eras)]
    feature_columns = list(feature_columns)
    dtrain = lgb.Dataset(get_X_array(train, feature_columns), get_Y_array(train), feature_name=feature_columns, free_raw_data=False)
    dvalid = lgb.Dataset(get_X_array(test, feature_columns), get_Y_array(test), feature_name=feature_columns, free_raw_data=False)
    np.random.seed(SEED)
    evals_result = {}
    model = lgb.train(hyperparameters, dtrain, num_boost_round=20000, valid_sets=(dvalid,), valid_names=('valid',), verbose_eval=False,
                early_stopping_rounds=100, evals_result=evals_result)
    return model

def get_X_array(df, feature_columns):
    return df.loc[:, feature_columns].values

def get_Y_array(df):
    return df.loc[:, 'mass diff'].values

def get_features(dataframe):
    feature_columns = dataframe.columns
    columns_to_remove = []
    allowed_features = {
        'mz',
        'RT',
        'Intensity',
    }
    for feature in feature_columns:
        if feature not in allowed_features:
            columns_to_remove.append(feature)
    feature_columns = feature_columns.drop(columns_to_remove)
    return feature_columns

def worker_RT(qin, qout, shift, step, RC=False, elude_path=False, ns=False, nr=False, win_sys=False):



    pepdict = dict()    
    if elude_path:
        outtrain = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
        outres = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
        outres_name = outres.name
        outres.close()
        for seq, RT in zip(ns, nr):
            outtrain.write(seq + '\t' + str(RT) + '\n')
        outtrain.flush()

        outtest = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')

        maxval = len(qin)
        start = 0
        while start + shift < maxval:
            item = qin[start+shift]
            outtest.write(item + '\n')
            start += step
        outtest.flush()

        subprocess.call([elude_path, '-t', outtrain.name, '-e', outtest.name, '-a', '-o', outres_name])
        for x in open(outres_name).readlines()[3:]:
            seq, RT = x.strip().split('\t')
            pepdict[seq] = float(RT)
        outtest.close()
        outtrain.close()
    else:
        maxval = len(qin)
        start = 0
        while start + shift < maxval:
            item = qin[start+shift]
            pepdict[item] = achrom.calculate_RT(item, RC)
            start += step
    if win_sys:
        return pepdict
    else:
        qout.put(pepdict)
        qout.put(None)

def final_iteration(resdict, mass_diff, rt_diff, pept_prot, protsN, base_out_name, prefix, isdecoy, isdecoy_key, escore, fdr, nproc, fname=False):
    # n = args['nproc']
    n = nproc


    prots_spc_basic = dict()

    p1 = set(resdict['seqs'])

    pep_pid = defaultdict(set)
    pid_pep = defaultdict(set)
    banned_dict = dict()
    for pep, pid in zip(resdict['seqs'], resdict['ids']):
        pep_pid[pep].add(pid)
        pid_pep[pid].add(pep)
        if pep in banned_dict:
            banned_dict[pep] += 1
        else:
            banned_dict[pep] = 1

    if len(p1):
        prots_spc_final = dict()
        prots_spc_copy = False
        prots_spc2 = False
        unstable_prots = set()
        p0 = False

        names_arr = False
        tmp_spc_new = False
        decoy_set = False

        while 1:
            if not prots_spc2:

                best_match_dict = dict()
                n_map_dict = defaultdict(list)
                for k, v in protsN.items():
                    n_map_dict[v].append(k)

                decoy_set = set()
                for k in protsN:
                    if isdecoy_key(k):
                        decoy_set.add(k)
                decoy_set = list(decoy_set)
                

                prots_spc2 = defaultdict(set)
                for pep, proteins in pept_prot.items():
                    if pep in p1:
                        for protein in proteins:
                            prots_spc2[protein].add(pep)

                for k in protsN:
                    if k not in prots_spc2:
                        prots_spc2[k] = set([])
                prots_spc2 = dict(prots_spc2)
                unstable_prots = set(prots_spc2.keys())

                top100decoy_N = sum([val for key, val in protsN.items() if isdecoy_key(key)])

                names_arr = np.array(list(prots_spc2.keys()))
                n_arr = np.array([protsN[k] for k in names_arr])

                tmp_spc_new = dict((k, len(v)) for k, v in prots_spc2.items())


                top100decoy_score_tmp = [tmp_spc_new.get(dprot, 0) for dprot in decoy_set]
                top100decoy_score_tmp_sum = float(sum(top100decoy_score_tmp))

            tmp_spc = tmp_spc_new
            prots_spc = tmp_spc_new
            if not prots_spc_copy:
                prots_spc_copy = deepcopy(prots_spc)

            for idx, v in enumerate(decoy_set):
                if v in unstable_prots:
                    top100decoy_score_tmp_sum -= top100decoy_score_tmp[idx]
                    top100decoy_score_tmp[idx] = prots_spc.get(v, 0)
                    top100decoy_score_tmp_sum += top100decoy_score_tmp[idx]
            p = float(sum(top100decoy_score_tmp)) / top100decoy_N
            p = top100decoy_score_tmp_sum / top100decoy_N
            if not p0:
                p0 = float(p)

            n_change = set(protsN[k] for k in unstable_prots)
            for n_val in n_change:
                for k in n_map_dict[n_val]:
                    v = prots_spc[k]
                    if n_val not in best_match_dict or v > prots_spc[best_match_dict[n_val]]:
                        best_match_dict[n_val] = k
            n_arr_small = []
            names_arr_small = []
            v_arr_small = []
            for k, v in best_match_dict.items():
                n_arr_small.append(k)
                names_arr_small.append(v)
                v_arr_small.append(prots_spc[v])

            prots_spc_basic = dict()
            all_pvals = calc_sf_all(v_arr_small, n_arr_small, p)
            for idx, k in enumerate(names_arr_small):
                prots_spc_basic[k] = all_pvals[idx]

            best_prot = utils.keywithmaxval(prots_spc_basic)

            best_score = prots_spc_basic[best_prot]
            unstable_prots = set()
            if best_prot not in prots_spc_final:
                prots_spc_final[best_prot] = best_score
                banned_pids = set()
                for pep in prots_spc2[best_prot]:
                    for pid in pep_pid[pep]:
                        banned_pids.add(pid)
                for pid in banned_pids:
                    for pep in pid_pep[pid]:
                        banned_dict[pep] -= 1
                        if banned_dict[pep] == 0:
                            for bprot in pept_prot[pep]:
                                tmp_spc_new[bprot] -= 1
                                unstable_prots.add(bprot)
            else:

                v_arr = np.array([prots_spc[k] for k in names_arr])
                all_pvals = calc_sf_all(v_arr, n_arr, p)
                for idx, k in enumerate(names_arr):
                    prots_spc_basic[k] = all_pvals[idx]

                for k, v in prots_spc_basic.items():
                    if k not in prots_spc_final:
                        prots_spc_final[k] = v

                break

            prot_fdr = aux.fdr(prots_spc_final.items(), is_decoy=isdecoy)
            if prot_fdr >= 12.5 * fdr:

                v_arr = np.array([prots_spc[k] for k in names_arr])
                all_pvals = calc_sf_all(v_arr, n_arr, p)
                for idx, k in enumerate(names_arr):
                    prots_spc_basic[k] = all_pvals[idx]

                for k, v in prots_spc_basic.items():
                    if k not in prots_spc_final:
                        prots_spc_final[k] = v
                break

    prots_spc_basic2 = copy(prots_spc_final)
    prots_spc_final = dict()
    prots_spc_final2 = dict()

    if n == 0:
        try:
            n = cpu_count()
        except NotImplementedError:
            n = 1

    if n == 1 or os.name == 'nt':
        qin = []
        qout = []
        for mass_koef in range(10):
            rtt_koef = mass_koef
            qin.append((mass_koef, rtt_koef))
        # for mass_koef in np.arange(1.0, 0.2, -0.33):
        #     for rtt_koef in np.arange(1.0, 0.2, -0.33):
                # qin.append((mass_koef, rtt_koef))
        qout = worker(qin, qout, mass_diff, rt_diff, resdict, protsN, pept_prot, isdecoy_key, isdecoy, fdr, prots_spc_basic2, True)
        
        for item, item2 in qout:
            if item2:
                prots_spc_copy = item2
            for k in protsN:
                # if k not in prots_spc_final:
                #     prots_spc_final[k] = [item.get(k, 0.0), ]
                # else:
                #     prots_spc_final[k].append(item.get(k, 0.0))
                if k not in prots_spc_final:
                    prots_spc_final[k] = [item.get(k, [0, 0])[0], ]
                    prots_spc_final2[k] = [item.get(k, [0, 0])[1], ]
                else:
                    prots_spc_final[k].append(item.get(k, [0, 0])[0])
                    prots_spc_final2[k].append(item.get(k, [0, 0])[1])

    else:
        qin = Queue()
        qout = Queue()

        # for mass_koef in np.arange(1.0, 0.2, -0.33):
        #     for rtt_koef in np.arange(1.0, 0.2, -0.33):
        #         qin.put((mass_koef, rtt_koef))
        for mass_koef in range(10):
            rtt_koef = mass_koef
            qin.put((mass_koef, rtt_koef))

        for _ in range(n):
            qin.put(None)

        procs = []
        for proc_num in range(n):
            p = Process(target=worker, args=(qin, qout, mass_diff, rt_diff, resdict, protsN, pept_prot, isdecoy_key, isdecoy, fdr, prots_spc_basic2))
            p.start()
            procs.append(p)

        for _ in range(n):
            for item, item2 in iter(qout.get, None):
                if item2:
                    prots_spc_copy = item2
                for k in protsN:
                    # if k not in prots_spc_final:
                    #     prots_spc_final[k] = [item.get(k, 0.0), ]
                    # else:
                    #     prots_spc_final[k].append(item.get(k, 0.0))
                    if k not in prots_spc_final:
                        prots_spc_final[k] = [item.get(k, [0, 0])[0], ]
                        prots_spc_final2[k] = [item.get(k, [0, 0])[1], ]
                    else:
                        prots_spc_final[k].append(item.get(k, [0, 0])[0])
                        prots_spc_final2[k].append(item.get(k, [0, 0])[1])

        for p in procs:
            p.join()
        # worker(qin, qout, mass_diff, rt_diff, resdict, protsN, pept_prot, isdecoy_key, isdecoy, fdr, prots_spc_basic2)
        # prots_spc_final, prots_spc_copy = qout.get()


    pickle.dump(protsN, open('/home/mark/protsN.pickle', 'wb'))
    pickle.dump(prots_spc_final, open('/home/mark/prots_spc_final.pickle', 'wb'))
    pickle.dump(prots_spc_final2, open('/home/mark/prots_spc_final2.pickle', 'wb'))


    all_dbnames = []
    all_theor_peptides = []
    for item in protsN.items():
        all_dbnames.append(item[0])
        all_theor_peptides.append(item[1])
        
        
    df2 = pd.DataFrame()
    df2['dbname'] = all_dbnames
    df2['theor peptides'] = all_theor_peptides

    df2['decoy'] = df2['dbname'].apply(lambda x: x.startswith('DECOY_'))

    zero_ar = [-10 for _ in range(n)]
    for qr in range(n):
        df2['qr_%d' % (qr, )] = df2['dbname'].apply(lambda x: prots_spc_final.get(x, zero_ar)[qr])


    for qr in range(n):
        df2['qr_%d_corr' % (qr, )] = df2['dbname'].apply(lambda x: prots_spc_final2.get(x, zero_ar)[qr])

    # Hyperparameter grid
    param_grid_prots = {
    #     'boosting_type': ['gbdt', 'goss', 'dart'],
    #     'boosting_type': ['gbdt', 'goss'],
        'boosting_type': ['gbdt', ],
    #     'boosting_type': ['dart', ],
        'num_leaves': list(range(10, 1000)),
        'learning_rate': list(np.logspace(np.log10(0.01), np.log10(0.3), base = 10, num = 1000)),
        'subsample_for_bin': list(range(1, 500, 5)),
        'min_child_samples': list(range(1, 150, 1)),
        'reg_alpha': list(np.linspace(0, 1)),
        'reg_lambda': list(np.linspace(0, 1)),
        'colsample_bytree': list(np.linspace(0.01, 1, 100)),
        'subsample': list(np.linspace(0.01, 1, 100)),
        'is_unbalance': [True, False],
        'metric': ['rmse', ],
        'verbose': [-1, ],
        'num_threads': [5, ],
    }


    # print('Start Machine Learning on proteins...')
    # MAX_EVALS = 25
    # out_file = 'test_randomCV_proteins2.csv'
    # of_connection = open(out_file, 'w')
    # writer = csv.writer(of_connection)

    # # Write column names
    # headers = ['auc', 'params', 'iteration']
    # writer.writerow(headers)
    # of_connection.close()

    # # df = df.reset_index(drop=True)

    # random_results = random_search_prots(df2, param_grid_prots, out_file, MAX_EVALS)
    # random_results = pd.read_csv('test_randomCV_proteins2.csv')
    # random_results = random_results[random_results['auc'] != 'auc']
    # random_results['params'] = random_results['params'].apply(lambda x: ast.literal_eval(x))
    # random_results['boosts'] = random_results['params'].apply(lambda x: int(x['n_estimators']))
    # convert_dict = {'auc': float, 
    #             } 
    # random_results = random_results.astype(convert_dict) 

    # random_results['auc_per_boost'] = random_results['auc'].values / random_results['boosts'].values

    # bestparams = random_results.sort_values(by='auc',ascending=False)['params'].values[0]
    # best_auc = random_results.sort_values(by='auc',ascending=False)['auc'].values[0]
    # print(best_auc)

    # # if best_auc >= 1:

    # kf = KFold(n_splits=3, shuffle=True, random_state=SEED)

    # test_res = []
    # for sp1 in kf.split(df2):
    #     train_df = df2.iloc[sp1[0]]
    #     test_df = df2.iloc[sp1[1]]

    #     feature_columns = list(get_features_prots(train_df))
    #     model = get_cat_model_final_prots(train_df, bestparams, feature_columns)
        
    #     df2.loc[sp1[1], 'preds'] = model.predict(get_X_array(test_df, feature_columns))

        
    #     df2['preds'] = df2['preds'].max() - df2['preds']

    # else:
    #     print('Skipping ML for proteins, use simple binomial scores')
    #     cols_for_sum = []
    #     for cc in df2.columns:
    #         if cc.startswith('sf_peptide'):
    #             cols_for_sum.append(cc)
    #     df2['preds'] = df2[cols_for_sum].sum(axis=1)


    # df2 = df2.rename({'preds': 'score', 'peptides_total': 'matched peptides', 'theor peptides': 'theoretical peptides'}, axis='columns')
    # df2.to_csv(base_out_name + '_proteins_full.csv', sep='\t', index=False, columns=['dbname', 'score', 'matched peptides', 'theoretical peptides'])

    # df2['basedbname'] = df2['dbname'].apply(lambda x: x.replace('DECOY_', ''))
    # df2 = df2.sort_values(by='score', ascending=False)
    # df2 = df2.drop_duplicates(subset=['basedbname'])

    # df2f = aux.filter(df2, fdr=fdr, key='score', is_decoy='decoy', reverse=True)
    # df2f.to_csv(base_out_name + '_proteins.csv', sep='\t', index=False, columns=['dbname', 'score', 'matched peptides', 'theoretical peptides'])

    # prots_spc_final = {}
    # for z in df2[['dbname', 'preds']].values:
    #     prots_spc_final[z[0]] = z[1]


    for k in prots_spc_final.keys():
        prots_spc_final[k] = np.mean(prots_spc_final[k])

    prots_spc = prots_spc_final
    sortedlist_spc = sorted(prots_spc.items(), key=operator.itemgetter(1))[::-1]
    with open(base_out_name + '_proteins_full.csv', 'w') as output:
        output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\n')
        for x in sortedlist_spc:
            # output.write('\t'.join((x[0], str(x[1]), str(protsN[x[0]]))) + '\n')
            output.write('\t'.join((x[0], str(x[1]), str(prots_spc_copy[x[0]]), str(protsN[x[0]]))) + '\n')

    checked = set()
    for k, v in list(prots_spc.items()):
        if k not in checked:
            if isdecoy_key(k):
                if prots_spc.get(k.replace(prefix, ''), -1e6) > v:
                    del prots_spc[k]
                    checked.add(k.replace(prefix, ''))
            else:
                if prots_spc.get(prefix + k, -1e6) > v:
                    del prots_spc[k]
                    checked.add(prefix + k)

    filtered_prots = aux.filter(prots_spc.items(), fdr=fdr, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1, full_output=True, correction=1)
    if len(filtered_prots) < 1:
        filtered_prots = aux.filter(prots_spc.items(), fdr=fdr, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1, full_output=True, correction=0)
    identified_proteins = 0

    for x in filtered_prots:
        identified_proteins += 1

    print('TOP 5 identified proteins:')
    print('dbname\tscore\tnum matched peptides\tnum theoretical peptides')
    for x in filtered_prots[:5]:
        # print '\t'.join((str(x[0]), str(x[1]), str(protsN[x[0]])))
        print('\t'.join((str(x[0]), str(x[1]), str(int(prots_spc_copy[x[0]])), str(protsN[x[0]]))))
    print('results:%s;number of identified proteins = %d' % (base_out_name, identified_proteins, ))
    # print('R=', r)
    with open(base_out_name + '_proteins.csv', 'w') as output:
        output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\n')
        for x in filtered_prots:
            # output.write('\t'.join((x[0], str(x[1]), str(protsN[x[0]]))) + '\n')
            output.write('\t'.join((x[0], str(x[1]), str(prots_spc_copy[x[0]]), str(protsN[x[0]]))) + '\n')


    if fname:
        fig = plt.figure(figsize=(16, 12))
        DPI = fig.get_dpi()
        fig.set_size_inches(2000.0/float(DPI), 2000.0/float(DPI))

        df0 = pd.read_table(os.path.splitext(fname)[0].replace('.features', '') + '.features' + '.tsv')

        # Features RT distribution
        # TODO add matched features and matched to 1% FDR proteins features
        ax = fig.add_subplot(3, 1, 1)
        bns = np.arange(0, df0['rtApex'].max() + 1, 1)
        ax.hist(df0['rtApex'], bins = bns)
        ax.set_xlabel('RT, min', size=16)
        ax.set_ylabel('# features', size=16)

        # Features mass distribution

        # TODO add matched features and matched to 1% FDR proteins features
        ax = fig.add_subplot(3, 1, 2)
        bns = np.arange(0, df0['massCalib'].max() + 6, 5)
        ax.hist(df0['massCalib'], bins = bns)
        ax.set_xlabel('neutral mass, Da', size=16)
        ax.set_ylabel('# features', size=16)

        # Features intensity distribution

        # TODO add matched features and matched to 1% FDR proteins features
        ax = fig.add_subplot(3, 1, 3)
        bns = np.arange(np.log10(df0['intensityApex'].min()) - 0.5, np.log10(df0['intensityApex'].max()) + 0.5, 0.5)
        ax.hist(np.log10(df0['intensityApex']), bins = bns)
        ax.set_xlabel('log10(Intensity)', size=16)
        ax.set_ylabel('# features', size=16)

        plt.savefig(base_out_name + '.png')

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

def calibrate_RT_gaus(bwidth, mass_left, mass_right, true_md):

    bbins = np.arange(-mass_left, mass_right, bwidth)
    H1, b1 = np.histogram(true_md, bins=bbins)
    b1 = b1 + bwidth
    b1 = b1[:-1]


    popt, pcov = curve_fit(noisygaus, b1, H1, p0=[1, np.median(true_md), bwidth * 5, 1])
    mass_shift, mass_sigma = popt[1], abs(popt[2])
    return mass_shift, mass_sigma, pcov[0][0]


def get_RCs2(sequences, RTs, lcp = -0.21,
            term_aa = False, **kwargs):
    labels = kwargs.get('labels')

    # Make a list of all amino acids present in the sample.
    peptide_dicts = [
            parser.amino_acid_composition(peptide, False, term_aa,
                               allow_unknown_modifications=True,
                               labels=labels)
            if not isinstance(peptide, dict) else peptide
        for peptide in sequences]

    detected_amino_acids = {aa for peptide_dict in peptide_dicts
                                for aa in peptide_dict}

    composition_array = []
    for pdict in peptide_dicts:
        loglen = np.log(parser.length(pdict))
        composition_array.append([pdict.get(aa, 0.)
             * (1. + lcp * loglen)
               for aa in detected_amino_acids] + [1.])

    if term_aa:
        for term_label in ['nterm', 'cterm']:
            normalizing_peptide = []
            for aa in detected_amino_acids:
                if aa.startswith(term_label):
                    normalizing_peptide.append(1.0)
                elif (term_label+aa) in detected_amino_acids:
                    normalizing_peptide.append(-1.0)
                else:
                    normalizing_peptide.append(0.0)
            normalizing_peptide.append(0.0)
            composition_array.append(normalizing_peptide)
            RTs.append(0.0)
    
    model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression(n_jobs=12), min_samples=0.5, max_trials=5000, random_state=42)
    model_ransac.fit(np.array(composition_array), np.array(RTs))
    RCs = model_ransac.estimator_.coef_
    if term_aa:
        for term_label in ['nterm', 'cterm']:
            RTs.pop()

    # Form output.
    RC_dict = {}
    RC_dict['aa'] = dict(
        zip(list(detected_amino_acids),
            RCs[:len(detected_amino_acids)]))
    RC_dict['aa'][parser.std_nterm] = 0.0
    RC_dict['aa'][parser.std_cterm] = 0.0
    RC_dict['const'] = RCs[len(detected_amino_acids)]
    RC_dict['lcp'] = lcp

    # Find remaining terminal RCs.
    if term_aa:
        for term_label in ['nterm', 'cterm']:
            # Check if there are terminal RCs remaining undefined.
            undefined_term_RCs = [aa for aa in RC_dict['aa']
                                if aa[1:5] != 'term'
                                and term_label + aa not in RC_dict['aa']]
            if not undefined_term_RCs:
                continue

            # Find a linear relationship between internal and terminal RCs.
            defined_term_RCs = [aa for aa in RC_dict['aa']
                              if aa[1:5] != 'term'
                              and term_label + aa in RC_dict['aa']]

            a, b, r, stderr = linear_regression(
                [RC_dict['aa'][aa] for aa in defined_term_RCs],
                [RC_dict['aa'][term_label+aa] for aa in defined_term_RCs])

            # Define missing terminal RCs using this linear equation.
            for aa in undefined_term_RCs:
                RC_dict['aa'][term_label + aa] = a * RC_dict['aa'][aa] + b
    inlier_mask = model_ransac.inlier_mask_
    outlier_mask = np.logical_not(inlier_mask)
    return RC_dict, outlier_mask


def process_file(args):
    # fname = args['file']
    # ftype = fname.rsplit('.', 1)[-1].lower()
    utils.seen_target.clear()
    utils.seen_decoy.clear()
    args = utils.prepare_decoy_db(args)
    return process_peptides(args)


def peptide_processor(peptide, **kwargs):
    seqm = peptide
    results = []
    m = cmass.fast_mass(seqm, aa_mass=kwargs['aa_mass']) + kwargs['aa_mass'].get('Nterm', 0) + kwargs['aa_mass'].get('Cterm', 0)
    acc_l = kwargs['acc_l']
    acc_r = kwargs['acc_r']
    dm_l = acc_l * m / 1.0e6
    if acc_r == acc_l:
        dm_r = dm_l
    else:
        dm_r = acc_r * m / 1.0e6
    start = nmasses.searchsorted(m - dm_l)
    end = nmasses.searchsorted(m + dm_r)
    for i in range(start, end):
        peak_id = ids[i]
        I = Is[i]
        massdiff = (m - nmasses[i]) / m * 1e6
        massdiff3 = (m - imraw[i]) / m * 1e6
        results.append((seqm, massdiff, rts[i], peak_id, I, Scans[i], Isotopes[i], mzraw[i], avraw[i], charges[i], imraw[i], massdiff3))
    return results


def prepare_peptide_processor(fname, args):
    global nmasses
    global rts
    global charges
    global ids
    global Is
    global Scans
    global Isotopes
    global mzraw
    global avraw
    global imraw
    nmasses = []
    rts = []
    charges = []
    ids = []
    Is = []
    Scans = []
    Isotopes = []
    mzraw = []
    avraw = []
    imraw = []

    min_ch = args['cmin']
    max_ch = args['cmax']

    min_isotopes = args['i']
    min_scans = args['sc']

    print('Reading spectra ...')
    for m, RT, c, peak_id, I, nScans, nIsotopes, mzr, avr, im in utils.iterate_spectra(fname, min_ch, max_ch, min_isotopes, min_scans):
        nmasses.append(m)
        rts.append(RT)
        charges.append(c)
        ids.append(peak_id)
        Is.append(I)
        Scans.append(nScans)
        Isotopes.append(nIsotopes)
        mzraw.append(mzr)
        avraw.append(avr)
        imraw.append(im)

    print(len(nmasses))

    # print(imraw)

    i = np.argsort(nmasses)
    nmasses = np.array(nmasses)[i]
    rts = np.array(rts)[i]
    charges = np.array(charges)[i]
    ids = np.array(ids)[i]
    Is = np.array(Is)[i]
    Scans = np.array(Scans)[i]
    Isotopes = np.array(Isotopes)[i]
    mzraw = np.array(mzraw)[i]
    avraw = np.array(avraw)[i]
    imraw = np.array(imraw)[i]

    fmods = args['fmods']
    aa_mass = mass.std_aa_mass
    if fmods:
        for mod in fmods.split(','):
            m, aa = mod.split('@')
            if aa == '[':
                aa_mass['Nterm'] = float(m)
            elif aa == ']':
                aa_mass['Cterm'] = float(m)
            else:
                aa_mass[aa] += float(m)

    acc_l = args['ptol']
    acc_r = args['ptol']

    return {'aa_mass': aa_mass, 'acc_l': acc_l, 'acc_r': acc_r, 'args': args}


def peptide_processor_iter_isoforms(peptide, **kwargs):
    out = []
    out.append(peptide_processor(peptide, **kwargs))
    return out

def get_results(ms1results):
    resdict = dict()
    labels = [
        'seqs',
        'md',
        'rt',
        'ids',
        'Is',
        'Scans',
        'Isotopes',
        'mzraw',
        'av',
        'ch',
        'im',
        'md3',
    ]
    for label, val in zip(labels, zip(*ms1results)):
        resdict[label] = np.array(val)
    return resdict

def filter_results(resultdict, idx):
    tmp = dict()
    for label in resultdict:
        tmp[label] = resultdict[label][idx]
    return tmp

def process_peptides(args):
    fname = args['file']
    fdr = args['fdr'] / 100
    min_isotopes_calibration = args['ci']
    try:
        outpath = args['outpath']
    except:
        outpath = False


    if outpath:
        base_out_name = os.path.splitext(os.path.join(outpath, os.path.basename(fname)))[0]
    else:
        base_out_name = os.path.splitext(fname)[0]

    out_log = open(base_out_name + '_log.txt', 'w')
    print(base_out_name + '_log.txt')
    out_log.close()
    out_log = open(base_out_name + '_log.txt', 'w')

    elude_path = args['elude']
    elude_path = elude_path.strip()

    deeplc_path = args['deeplc']
    deeplc_path = deeplc_path.strip()

    calib_path = args['pl']
    calib_path = calib_path.strip()

    if calib_path and args['ts']:
        args['ts'] = 0
        print('Two-stage RT prediction does not work with list of MS/MS identified peptides...')

    args['enzyme'] = utils.get_enzyme(args['e'])

    ms1results = []
    peps = utils.peptide_gen(args)
    kwargs = prepare_peptide_processor(fname, args)
    func = peptide_processor_iter_isoforms
    print('Running the search ...')
    for y in utils.multimap(1, func, peps, **kwargs):
        for result in y:
            if len(result):
                ms1results.extend(result)

    prefix = args['prefix']
    protsN, pept_prot = utils.get_prot_pept_map(args)

    resdict = get_results(ms1results)
    del ms1results

    resdict['mc'] = np.array([parser.num_sites(z, args['enzyme']) for z in resdict['seqs']])
    # resdict['mc'] = resdict['im']

    isdecoy = lambda x: x[0].startswith(prefix)
    isdecoy_key = lambda x: x.startswith(prefix)
    escore = lambda x: -x[1]

    # e_ind = resdict['mc'] == 0
    # resdict2 = filter_results(resdict, e_ind)

    e_ind = resdict['Isotopes'] >= min_isotopes_calibration
    # e_ind = resdict['Isotopes'] >= 1
    resdict2 = filter_results(resdict, e_ind)

    print(len(resdict2['seqs']))

    e_ind = resdict2['mc'] == 0
    resdict2 = filter_results(resdict2, e_ind)

    print(len(resdict2['seqs']))

    p1 = set(resdict2['seqs'])

    if len(p1):
        prots_spc2 = defaultdict(set)
        for pep, proteins in pept_prot.items():
            if pep in p1:
                for protein in proteins:
                    prots_spc2[protein].add(pep)

        for k in protsN:
            if k not in prots_spc2:
                prots_spc2[k] = set([])
        prots_spc = dict((k, len(v)) for k, v in prots_spc2.items())

        names_arr = np.array(list(prots_spc.keys()))
        v_arr = np.array(list(prots_spc.values()))
        n_arr = np.array([protsN[k] for k in prots_spc])

        top100decoy_score = [prots_spc.get(dprot, 0) for dprot in protsN if isdecoy_key(dprot)]
        top100decoy_N = [val for key, val in protsN.items() if isdecoy_key(key)]
        p = np.mean(top100decoy_score) / np.mean(top100decoy_N)
        print('p=%s' % (np.mean(top100decoy_score) / np.mean(top100decoy_N)))

        prots_spc = dict()
        all_pvals = calc_sf_all(v_arr, n_arr, p)
        for idx, k in enumerate(names_arr):
            prots_spc[k] = all_pvals[idx]

        checked = set()
        for k, v in list(prots_spc.items()):
            if k not in checked:
                if isdecoy_key(k):
                    if prots_spc.get(k.replace(prefix, ''), -1e6) > v:
                        del prots_spc[k]
                        checked.add(k.replace(prefix, ''))
                else:
                    if prots_spc.get(prefix + k, -1e6) > v:
                        del prots_spc[k]
                        checked.add(prefix + k)

        filtered_prots = aux.filter(prots_spc.items(), fdr=fdr, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1,
                                    full_output=True)

        print(sorted(prots_spc.items(), key=lambda x: -x[1])[:15])

        print(len(filtered_prots))

        identified_proteins = 0

        for x in filtered_prots:
            identified_proteins += 1
        print('results for default search: number of identified proteins = %d' % (identified_proteins, ))

        print('Running mass recalibration...')


        e_ind = resdict['mc'] == 0
        resdict2 = filter_results(resdict, e_ind)

        print(len(resdict2['seqs']))


        true_md = []
        true_isotopes = []
        true_seqs = []
        true_prots = set(x[0] for x in filtered_prots)
        for pep, proteins in pept_prot.items():
            if any(protein in true_prots for protein in proteins):
                true_seqs.append(pep)

        # true_seqs = list(resdict2['seqs'])

        e_ind = np.in1d(resdict2['seqs'], true_seqs)

        true_seqs = resdict2['seqs'][e_ind]
        true_md.extend(resdict2['md'][e_ind])
        true_md = np.array(true_md)
        true_isotopes.extend(resdict2['Isotopes'][e_ind])
        true_isotopes = np.array(true_isotopes)
        true_intensities = np.array(resdict2['Is'][e_ind])
        true_rt = np.array(resdict2['rt'][e_ind])
        true_mz = np.array(resdict2['mzraw'][e_ind])

        import pickle
        pickle.dump(true_md, open('/home/mark/true_md.pickle', 'wb'))
        pickle.dump(true_intensities, open('/home/mark/true_intensities.pickle', 'wb'))
        pickle.dump(true_rt, open('/home/mark/true_rt.pickle', 'wb'))
        pickle.dump(true_mz, open('/home/mark/true_mz.pickle', 'wb'))

        df1 = pd.DataFrame()
        df1['mass diff'] = true_md# * true_mz / 1e6
        df1['mz'] = true_mz
        df1['RT'] = true_rt
        df1['Intensity'] = true_intensities
        df1['seqs'] = true_seqs
        df1['orig_md'] = true_md

        mass_left = args['ptol']
        mass_right = args['ptol']

        try:
            mass_shift, mass_sigma, covvalue = calibrate_mass(0.001, mass_left, mass_right, true_md)
        except:
            mass_shift, mass_sigma, covvalue = calibrate_mass(0.01, mass_left, mass_right, true_md)

        print('Calibrated mass shift: ', mass_shift)
        print('Calibrated mass sigma in ppm: ', mass_sigma)

        out_log.write('Calibrated mass shift: %s\n' % (mass_shift, ))
        out_log.write('Calibrated mass sigma in ppm: %s\n' % (mass_sigma, ))

        e_all = abs(resdict['md'] - mass_shift) / (mass_sigma)
        r = 3.0
        e_ind = e_all <= r
        resdict = filter_results(resdict, e_ind)

        zs_all = e_all[e_ind] ** 2



        # e_all = abs(resdict['md3'] - mass_shift3) / (mass_sigma3)
        # r = 3.0
        # e_ind = e_all <= r
        # resdict = filter_results(resdict, e_ind)

        # zs_all = zs_all[e_ind]


        e_ind = resdict['mc'] == 0
        resdict2 = filter_results(resdict, e_ind)

        e_ind = resdict2['Isotopes'] >= min_isotopes_calibration
        resdict2 = filter_results(resdict2, e_ind)

        p1 = set(resdict2['seqs'])

        prots_spc2 = defaultdict(set)
        for pep, proteins in pept_prot.items():
            if pep in p1:
                for protein in proteins:
                    prots_spc2[protein].add(pep)

        for k in protsN:
            if k not in prots_spc2:
                prots_spc2[k] = set([])
        prots_spc = dict((k, len(v)) for k, v in prots_spc2.items())

        names_arr = np.array(list(prots_spc.keys()))
        v_arr = np.array(list(prots_spc.values()))
        n_arr = np.array([protsN[k] for k in prots_spc])

        top100decoy_score = [prots_spc.get(dprot, 0) for dprot in protsN if isdecoy_key(dprot)]
        top100decoy_N = [val for key, val in protsN.items() if isdecoy_key(key)]
        p = np.mean(top100decoy_score) / np.mean(top100decoy_N)
        print('p=%s' % (np.mean(top100decoy_score) / np.mean(top100decoy_N)))

        prots_spc = dict()
        all_pvals = calc_sf_all(v_arr, n_arr, p)
        for idx, k in enumerate(names_arr):
            prots_spc[k] = all_pvals[idx]

        checked = set()
        for k, v in list(prots_spc.items()):
            if k not in checked:
                if isdecoy_key(k):
                    if prots_spc.get(k.replace(prefix, ''), -1e6) > v:
                        del prots_spc[k]
                        checked.add(k.replace(prefix, ''))
                else:
                    if prots_spc.get(prefix + k, -1e6) > v:
                        del prots_spc[k]
                        checked.add(prefix + k)

        filtered_prots = aux.filter(prots_spc.items(), fdr=fdr, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1,
                                    full_output=True)

        identified_proteins = 0

        for x in filtered_prots:
            identified_proteins += 1
        print('results for default search after mass calibration: number of identified proteins = %d' % (identified_proteins, ))



        print('Running RT prediction...')


        e_ind = resdict['Isotopes'] >= min_isotopes_calibration
        # e_ind = resdict['Isotopes'] >= 1
        resdict2 = filter_results(resdict, e_ind)

        e_ind = resdict2['mc'] == 0
        resdict2 = filter_results(resdict2, e_ind)


        true_seqs = []
        true_rt = []
        true_isotopes = []
        true_prots = set(x[0] for x in filtered_prots)#[:5])
        for pep, proteins in pept_prot.items():
            if any(protein in true_prots for protein in proteins):
                true_seqs.append(pep)
        e_ind = np.in1d(resdict2['seqs'], true_seqs)


        true_seqs = resdict2['seqs'][e_ind]
        true_rt.extend(resdict2['rt'][e_ind])
        true_rt = np.array(true_rt)
        true_isotopes.extend(resdict2['Isotopes'][e_ind])
        true_isotopes = np.array(true_isotopes)

        e_all = abs(resdict2['md'][e_ind] - mass_shift) / (mass_sigma)
        # e_all3 = abs(resdict2['md3'][e_ind] - mass_shift3) / (mass_sigma3)
        zs_all_tmp = e_all ** 2# + e_all3 ** 2

        # e_ind = true_isotopes >= 1
        e_ind = true_isotopes >= min_isotopes_calibration
        true_seqs = true_seqs[e_ind]
        true_rt = true_rt[e_ind]
        true_isotopes = true_isotopes[e_ind]
        zs_all_tmp = zs_all_tmp[e_ind]

        e_ind = np.argsort(zs_all_tmp)
        true_seqs = true_seqs[e_ind]
        true_rt = true_rt[e_ind]
        true_isotopes = true_isotopes[e_ind]

        true_seqs = true_seqs[:2500]
        true_rt = true_rt[:2500]
        true_isotopes = true_isotopes[:2500]

        best_seq = defaultdict(list)
        newseqs = []
        newRTs = []
        for seq, RT in zip(true_seqs, true_rt):
            best_seq[seq].append(RT)
        for k, v in best_seq.items():
            newseqs.append(k)
            newRTs.append(np.median(v))
        true_seqs = np.array(newseqs)
        true_rt = np.array(newRTs)

        if calib_path:
            df1 = pd.read_csv(calib_path, sep='\t')
            true_seqs2 = df1['peptide'].values
            true_rt2 = df1['RT exp'].values

        else:
            true_seqs2 = true_seqs
            true_rt2 = true_rt
        # RC, outmask = get_RCs2(true_seqs, true_rt)

# ~/virtualenv_deeplc/bin/deeplc --file_pred ms1_all.txt --file_cal ms1_top.txt --file_pred_out test_out.txt
        if args['ts'] != 2 and deeplc_path:


            outtrain = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
            outcalib = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
            outres = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
            outres_name = outres.name
            outres.close()
            # ns = true_seqs[~outmask]
            # nr = true_rt[~outmask]
            ns = true_seqs
            nr = true_rt
            print(len(ns))
            ns2 = true_seqs2
            nr2 = true_rt2
            print(len(ns2))
            # ll = len(ns)
            # ns = ns[:ll]
            # nr = nr[:ll]

            # train_dict = {}

            outtrain.write('seq,modifications,tr\n')
            for seq, RT in zip(ns2, nr2):
                # train_dict[seq] = RT
                mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
                # df1['modifications'] = df1['seq'].apply(lambda x: '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(x) if aa == 'C']))
                outtrain.write(seq + ',' + str(mods_tmp) + ',' + str(RT) + '\n')
            outtrain.flush()

            outcalib.write('seq,modifications,tr\n')
            for seq, RT in zip(ns, nr):
                # train_dict[seq] = RT
                mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
                # df1['modifications'] = df1['seq'].apply(lambda x: '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(x) if aa == 'C']))
                outcalib.write(seq + ',' + str(mods_tmp) + ',' + str(RT) + '\n')
            outcalib.flush()

            subprocess.call([deeplc_path, '--file_pred', outcalib.name, '--file_cal', outtrain.name, '--file_pred_out', outres_name])
            # subprocess.call([elude_path, '-t', outtrain.name, '-e', outtrain.name, '-a', '-g', '-o', outres_name])
            pepdict = dict()
            train_RT = []
            train_seq = []
            for x in open(outres_name).readlines()[1:]:
                # print(x.strip())
                _, seq, _, RTexp, RT = x.strip().split(',')
                pepdict[seq] = float(RT)
                train_seq.append(seq)
                train_RT.append(float(RTexp))
                # train_RT.append(float(train_dict[seq]))


            train_RT = np.array(train_RT)
            RT_pred = np.array([pepdict[s] for s in train_seq])

            rt_diff_tmp = RT_pred - train_RT
            RT_left = -min(rt_diff_tmp)
            RT_right = max(rt_diff_tmp)

            # import random
            # random.shuffle(rt_diff_tmp)

            try:
                start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 100
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
            except:
                start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 50
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
            if np.isinf(covvalue):
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(0.1, RT_left, RT_right, rt_diff_tmp)
            if np.isinf(covvalue):
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(1.0, RT_left, RT_right, rt_diff_tmp)
            print('Calibrated RT shift: ', XRT_shift)
            print('Calibrated RT sigma: ', XRT_sigma)

            # ns = np.array(ns)
            # nr = np.array(nr)
            # idx = np.abs((rt_diff_tmp) - XRT_shift) <= 3 * XRT_sigma
            # ns = ns[idx]
            # nr = nr[idx]
            # RT_pred = RT_pred[idx]
            # train_RT = train_RT[idx]

            aa, bb, RR, ss = aux.linear_regression(RT_pred, train_RT)

        else:

            if args['ts'] != 2 and elude_path:
                outtrain = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
                outcalib = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
                outres = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
                outres_name = outres.name
                outres.close()
                # ns = true_seqs[~outmask]
                # nr = true_rt[~outmask]
                ns = true_seqs
                nr = true_rt
                print(len(ns))
                ns2 = true_seqs2
                nr2 = true_rt2
                print(len(ns2))
                # ll = len(ns)
                # ns = ns[:ll]
                # nr = nr[:ll]
                for seq, RT in zip(ns, nr):
                    outtrain.write(seq + '\t' + str(RT) + '\n')
                outtrain.flush()
                for seq, RT in zip(ns, nr):
                    outcalib.write(seq + '\t' + str(RT) + '\n')
                outcalib.flush()

                subprocess.call([elude_path, '-t', outtrain.name, '-e', outcalib.name, '-a', '-g', '-o', outres_name])
                pepdict = dict()
                train_RT = []
                train_seq = []
                for x in open(outres_name).readlines()[3:]:
                    seq, RT, RTexp = x.strip().split('\t')
                    pepdict[seq] = float(RT)
                    train_seq.append(seq)
                    train_RT.append(float(RTexp))
                train_RT = np.array(train_RT)
                RT_pred = np.array([pepdict[s] for s in train_seq])

                rt_diff_tmp = RT_pred - train_RT
                RT_left = -min(rt_diff_tmp)
                RT_right = max(rt_diff_tmp)

                # import random
                # random.shuffle(rt_diff_tmp)

                start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 50
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
                if np.isinf(covvalue):
                    XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(0.1, RT_left, RT_right, rt_diff_tmp)
                if np.isinf(covvalue):
                    XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(1.0, RT_left, RT_right, rt_diff_tmp)
                print('Calibrated RT shift: ', XRT_shift)
                print('Calibrated RT sigma: ', XRT_sigma)

                aa, bb, RR, ss = aux.linear_regression(RT_pred, train_RT)
            else:
                ns = true_seqs
                nr = true_rt
                ns2 = true_seqs2
                nr2 = true_rt2
                # RC = achrom.get_RCs_vary_lcp(true_seqs[~outmask], true_rt[~outmask])
                RC = achrom.get_RCs_vary_lcp(ns2, nr2)
                RT_pred = np.array([achrom.calculate_RT(s, RC) for s in ns])
                train_RT = nr
                aa, bb, RR, ss = aux.linear_regression(RT_pred, nr)

                rt_diff_tmp = RT_pred - nr
                RT_left = -min(rt_diff_tmp)
                RT_right = max(rt_diff_tmp)

                start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 50
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
                if np.isinf(covvalue):
                    XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(0.1, RT_left, RT_right, rt_diff_tmp)
                if np.isinf(covvalue):
                    XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(1.0, RT_left, RT_right, rt_diff_tmp)
                print('Calibrated RT shift: ', XRT_shift)
                print('Calibrated RT sigma: ', XRT_sigma)

        print(aa, bb, RR, ss)



        best_sigma = XRT_sigma#ss
        RT_sigma = XRT_sigma#best_sigma

    else:
        print('No matches found')

    if args['ts']:

        print('Running second stage RT prediction...')


        ns = np.array(ns)
        nr = np.array(nr)
        idx = np.abs((rt_diff_tmp) - XRT_shift) <= 3 * XRT_sigma
        ns = ns[idx]
        nr = nr[idx]

        if deeplc_path:


            outtrain = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
            outres = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
            outres_name = outres.name
            outres.close()
            print(len(ns))
            ll = len(ns)
            ns = ns[:ll]
            nr = nr[:ll]

            # train_dict = {}

            outtrain.write('seq,modifications,tr\n')
            for seq, RT in zip(ns, nr):
                # train_dict[seq] = RT
                mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
                # df1['modifications'] = df1['seq'].apply(lambda x: '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(x) if aa == 'C']))
                outtrain.write(seq + ',' + str(mods_tmp) + ',' + str(RT) + '\n')
            outtrain.flush()

            subprocess.call([deeplc_path, '--file_pred', outtrain.name, '--file_cal', outtrain.name, '--file_pred_out', outres_name])
            # subprocess.call([elude_path, '-t', outtrain.name, '-e', outtrain.name, '-a', '-g', '-o', outres_name])
            pepdict = dict()
            train_RT = []
            train_seq = []
            for x in open(outres_name).readlines()[1:]:
                # print(x.strip())
                _, seq, _, RTexp, RT = x.strip().split(',')
                pepdict[seq] = float(RT)
                train_seq.append(seq)
                train_RT.append(float(RTexp))
                # train_RT.append(float(train_dict[seq]))


            train_RT = np.array(train_RT)
            RT_pred = np.array([pepdict[s] for s in train_seq])

            rt_diff_tmp = RT_pred - train_RT
            RT_left = -min(rt_diff_tmp)
            RT_right = max(rt_diff_tmp)

            # import random
            # random.shuffle(rt_diff_tmp)

            try:
                start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 100
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
            except:
                start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 50
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
            if np.isinf(covvalue):
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(0.1, RT_left, RT_right, rt_diff_tmp)
            if np.isinf(covvalue):
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(1.0, RT_left, RT_right, rt_diff_tmp)
            print('Calibrated RT shift: ', XRT_shift)
            print('Calibrated RT sigma: ', XRT_sigma)

            # ns = np.array(ns)
            # nr = np.array(nr)
            # idx = np.abs((rt_diff_tmp) - XRT_shift) <= 3 * XRT_sigma
            # ns = ns[idx]
            # nr = nr[idx]
            # RT_pred = RT_pred[idx]
            # train_RT = train_RT[idx]

            aa, bb, RR, ss = aux.linear_regression(RT_pred, train_RT)

        else:

            if elude_path:
                outtrain = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
                outres = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
                outres_name = outres.name
                outres.close()
                # ns = true_seqs[~outmask]
                # nr = true_rt[~outmask]
                # ns = true_seqs
                # nr = true_rt
                print(len(ns))
                ll = len(ns)
                ns = ns[:ll]
                nr = nr[:ll]
                for seq, RT in zip(ns, nr):
                    outtrain.write(seq + '\t' + str(RT) + '\n')
                outtrain.flush()

                subprocess.call([elude_path, '-t', outtrain.name, '-e', outtrain.name, '-a', '-g', '-o', outres_name])
                pepdict = dict()
                train_RT = []
                train_seq = []
                for x in open(outres_name).readlines()[3:]:
                    seq, RT, RTexp = x.strip().split('\t')
                    pepdict[seq] = float(RT)
                    train_seq.append(seq)
                    train_RT.append(float(RTexp))
                train_RT = np.array(train_RT)
                RT_pred = np.array([pepdict[s] for s in train_seq])

                rt_diff_tmp = RT_pred - train_RT
                RT_left = -min(rt_diff_tmp)
                RT_right = max(rt_diff_tmp)

                # import random
                # random.shuffle(rt_diff_tmp)

                start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 50
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
                if np.isinf(covvalue):
                    XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(0.1, RT_left, RT_right, rt_diff_tmp)
                if np.isinf(covvalue):
                    XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(1.0, RT_left, RT_right, rt_diff_tmp)
                print('Calibrated RT shift: ', XRT_shift)
                print('Calibrated RT sigma: ', XRT_sigma)

                aa, bb, RR, ss = aux.linear_regression(RT_pred, train_RT)
            else:
                # RC = achrom.get_RCs_vary_lcp(true_seqs[~outmask], true_rt[~outmask])
                RC = achrom.get_RCs_vary_lcp(ns, nr)
                RT_pred = np.array([achrom.calculate_RT(s, RC) for s in ns])
                aa, bb, RR, ss = aux.linear_regression(RT_pred, nr)

                rt_diff_tmp = RT_pred - nr
                RT_left = -min(rt_diff_tmp)
                RT_right = max(rt_diff_tmp)

                start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 50
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
                if np.isinf(covvalue):
                    XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(0.1, RT_left, RT_right, rt_diff_tmp)
                if np.isinf(covvalue):
                    XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(1.0, RT_left, RT_right, rt_diff_tmp)
                print('Calibrated RT shift: ', XRT_shift)
                print('Calibrated RT sigma: ', XRT_sigma)

        print(aa, bb, RR, ss)



        best_sigma = XRT_sigma#ss
        RT_sigma = XRT_sigma#best_sigma


    out_log.write('Calibrated RT shift: %s\n' % (XRT_shift, ))
    out_log.write('Calibrated RT sigma: %s\n' % (XRT_sigma, ))
    out_log.close()

    p1 = set(resdict['seqs'])

    n = args['nproc']

    if deeplc_path:


        pepdict = dict()

        outtrain = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
        outtrain2 = open('/home/mark/deeplc_train.txt', 'w')
        outres = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
        outres_name = outres.name
        outres.close()


        outtrain.write('seq,modifications,tr\n')
        outtrain2.write('seq,modifications,tr\n')
        for seq, RT in zip(ns, nr):
            mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
            outtrain.write(seq + ',' + str(mods_tmp) + ',' + str(RT) + '\n')
            outtrain2.write(seq + ',' + str(mods_tmp) + ',' + str(RT) + '\n')
        outtrain.flush()

        outtest = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')


        outtest.write('seq,modifications\n')
        for seq in p1:
            mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
            outtest.write(seq + ',' + str(mods_tmp) + '\n')
        outtest.flush()

        subprocess.call([deeplc_path, '--file_pred', outtest.name, '--file_cal', outtrain.name, '--file_pred_out', outres_name])
        # subprocess.call([elude_path, '-t', outtrain.name, '-e', outtrain.name, '-a', '-g', '-o', outres_name])
        for x in open(outres_name).readlines()[1:]:
            _, seq, _, RT = x.strip().split(',')
            pepdict[seq] = float(RT)

        outtest.close()
        outtrain.close()
        outtrain2.close()
    else:

        if n == 1 or os.name == 'nt':
            qin = list(p1)
            qout = []
            if elude_path:
                pepdict = worker_RT(qin, qout, 0, 1, False, elude_path, ns, nr, True)
            else:
                pepdict = worker_RT(qin, qout, 0, 1, RC, False, False, False, True)
        else:
            qin = list(p1)
            qout = Queue()
            procs = []
            for i in range(n):
                if elude_path:
                    p = Process(target=worker_RT, args=(qin, qout, i, n, False, elude_path, ns, nr))
                else:
                    p = Process(target=worker_RT, args=(qin, qout, i, n, RC, False, False, False))
                p.start()
                procs.append(p)

            pepdict = dict()
            for _ in range(n):
                for item in iter(qout.get, None):
                    for k, v in item.items():
                        pepdict[k] = v

            for p in procs:
                p.join()


    rt_pred = np.array([pepdict[s] for s in resdict['seqs']])
    rt_diff = resdict['rt'] - rt_pred
    # random.shuffle(rt_diff)
    e_all = (rt_diff) ** 2 / (RT_sigma ** 2)
    zs_all = zs_all + e_all
    r = 9.0#16.0#9.0
    e_ind = e_all <= r
    resdict = filter_results(resdict, e_ind)
    rt_diff = rt_diff[e_ind]
    zs_all = zs_all[e_ind]
    rt_pred = rt_pred[e_ind]


    with open(base_out_name + '_protsN.csv', 'w') as output:
        output.write('dbname\ttheor peptides\n')
        for k, v in protsN.items():
            output.write('\t'.join((k, str(v))) + '\n')
   
    with open(base_out_name + '_PFMs.csv', 'w') as output:
        output.write('sequence\tmass diff\tRT diff\tpeak_id\tIntensity\tnScans\tnIsotopes\tproteins\tm/z\tRT\taveragineCorr\tcharge\tion_mobility\n')
        for seq, md, rtd, peak_id, I, nScans, nIsotopes, mzr, rtr, av, ch, im in zip(resdict['seqs'], resdict['md'], rt_diff, resdict['ids'], resdict['Is'], resdict['Scans'], resdict['Isotopes'], resdict['mzraw'], resdict['rt'], resdict['av'], resdict['ch'], resdict['im']):
            output.write('\t'.join((seq, str(md), str(rtd), str(peak_id), str(I), str(nScans), str(nIsotopes), ';'.join(pept_prot[seq]), str(mzr), str(rtr), str(av), str(ch), str(im))) + '\n')
            
    # mass_diff = np.abs(resdict['md'] - mass_shift) / (mass_sigma)
    # rt_diff = np.abs(resdict['rt'] - rt_pred) / RT_sigma
    mass_diff = (resdict['md'] - mass_shift) / (mass_sigma)
    rt_diff = (resdict['rt'] - rt_pred) / RT_sigma

    import pickle
    pickle.dump(resdict, open('/home/mark/resdict.pickle', 'wb'))
    pickle.dump(mass_diff, open('/home/mark/mass_diff.pickle', 'wb'))
    pickle.dump(rt_diff, open('/home/mark/rt_diff.pickle', 'wb'))
    pickle.dump(pept_prot, open('/home/mark/pept_prot.pickle', 'wb'))
    pickle.dump(protsN, open('/home/mark/protsN.pickle', 'wb'))



    prefix = 'DECOY_'
    isdecoy = lambda x: x[0].startswith(prefix)
    isdecoy_key = lambda x: x.startswith(prefix)
    escore = lambda x: -x[1]
    fdr = 0.01


    SEED = 42

    # Hyperparameter grid
    param_grid = {
    #     'boosting_type': ['gbdt', 'goss', 'dart'],
    #     'boosting_type': ['gbdt', 'goss'],
        'boosting_type': ['gbdt', ],
    #     'boosting_type': ['dart', ],
        'num_leaves': list(range(10, 1000)),
        'learning_rate': list(np.logspace(np.log10(0.001), np.log10(0.05), base = 10, num = 1000)),
        # 'subsample_for_bin': list(range(10, 10000, 100)),
        # 'min_child_samples': list(range(1, 2500, 5)),

        # # 'subsample_for_bin': list(range(1, 1000, 10)),
        # # 'min_child_samples': list(range(1, 250, 5)),

        # 'reg_alpha': list(np.linspace(0, 1)),
        # 'reg_lambda': list(np.linspace(0, 1)),
        # 'colsample_bytree': list(np.linspace(0.1, 1, 10)),
        # 'subsample': list(np.linspace(0.1, 1, 10)),
        # 'is_unbalance': [True, False],
        'metric': ['rmse', ],
        'verbose': [-1, ],
        'num_threads': [5, ],
        'n_estimators': [100, ],
    }

    def get_X_array(df, feature_columns):
        return df.loc[:, feature_columns].values

    def get_Y_array_pfms(df):
        return df.loc[:, 'decoy'].values

    def get_features_Is(dataframe):
        feature_columns = dataframe.columns
        columns_to_remove = []
        for feature in feature_columns:
            if feature not in ['plen', 'mass', 'pI', 'charge_theor', 'ch']:
                if not feature.startswith('c_'):
                    columns_to_remove.append(feature)
        feature_columns = feature_columns.drop(columns_to_remove)
        return feature_columns

    def get_features_FAIMS(dataframe):
        feature_columns = dataframe.columns
        columns_to_remove = []
        for feature in feature_columns:
            if feature not in ['plen', 'mass', 'pI', 'charge_theor', 'ch', 'mzraw']:
                if not feature.startswith('c_'):
                    columns_to_remove.append(feature)
        feature_columns = feature_columns.drop(columns_to_remove)
        return feature_columns

    def get_Y_array_Is(df):
        return df.loc[:, 'IntensityNorm'].values

    def get_Y_array_FAIMS(df):
        return df.loc[:, 'im'].values

    def get_features_pfms(dataframe):
        feature_columns = dataframe.columns
        columns_to_remove = []
        banned_features = {
            'ids',
            'seqs',
            'decoy',
            'preds',
            'av',
            'best_prot_rank',
            # 'Is',
            'Scans',
            # 'Isotopes',
            'IntensityNorm',
            'IntensityNorm_predicted',
            'proteins',
            'peptide',
            'bestprotein',
            'IntensityNorm_diff',
            # 'charge_theor',
            # 'pI',
            # 'mass',
            # 'mzraw',
            'im_predicted',
            # 'im',
            # 'im_diff',

            
        }
        for feature in feature_columns:
            if feature in banned_features or feature.startswith('c_'):
                columns_to_remove.append(feature)
        feature_columns = feature_columns.drop(columns_to_remove)
        return feature_columns

    # def get_features_pfms(dataframe):
    #     feature_columns = dataframe.columns
    #     columns_to_remove = []
    #     banned_features = {
    #         'rt_diff',
    #         'mass_diff',
            
    #     }
    #     for feature in feature_columns:
    #         if feature not in banned_features:
    #             columns_to_remove.append(feature)
    #     feature_columns = feature_columns.drop(columns_to_remove)
    #     return feature_columns

    def objective_pfms(df, hyperparameters, iteration):
        """Objective function for grid and random search. Returns
        the cross validation score from a set of hyperparameters."""
        
        all_res = []
        # all_iters = []
        
        # Number of estimators will be found using early stopping
        # if 'n_estimators' in hyperparameters.keys():
        #     del hyperparameters['n_estimators']

        groups = df['peptide']
        # ix = np.arange(len(groups))
        ix = df.index.values
        unique = np.unique(groups)
        np.random.RandomState(SEED).shuffle(unique)
        result = []
        for split in np.array_split(unique, 3):
            mask = groups.isin(split)
            train, test = ix[~mask], ix[mask]
            train_df = df.iloc[train]
            test_df = df.iloc[test]
            # result.append((train, test))

            feature_columns = get_features_pfms(df)
            # print(feature_columns)
        #     kf = KFold(n_splits=3, shuffle=True, random_state=SEED+1)
        
        # for sp1 in kf.split(df):
        #     train_df = df.iloc[sp1[0]]
        #     test_df = df.iloc[sp1[1]]
            # model = get_cat_model_pfms(df, hyperparameters, feature_columns, train_df, test_df)
            model = get_cat_model_final_pfms(train_df, hyperparameters, feature_columns)
            
            # all_iters.append(model.best_iteration)
            
            test_df['preds'] = model.predict(get_X_array(test_df, feature_columns))
            all_res.append(test_df)
            # all_res.append(len(aux.filter(test_df, fdr=0.1, key='preds', is_decoy='decoy')))
        test_df6 = pd.concat(all_res)
        shr_v = len(aux.filter(test_df6, fdr=0.25, key='preds', is_decoy='decoy'))#np.min(all_res)
        # hyperparameters['n_estimators'] = int(np.mean(all_iters))
        
        return [shr_v, hyperparameters, iteration]

    def random_search_pfms(df, param_grid, out_file, max_evals):
        """Random search for hyperparameter optimization. 
        Writes result of search to csv file every search iteration."""
        
        
        # Dataframe for results
        results = pd.DataFrame(columns = ['sharpe', 'params', 'iteration'],
                                    index = list(range(max_evals)))
        for i in range(max_evals):

            print('%d/%d' % (i+1, max_evals))
            
            # Choose random hyperparameters
            random_params = {k: random.sample(v, 1)[0] for k, v in param_grid.items()}
            # random_params['subsample'] = 1.0 if random_params['boosting_type'] == 'goss' else random_params['subsample']

            # Evaluate randomly selected hyperparameters
            eval_results = objective_pfms(df, random_params, i)
            results.loc[i, :] = eval_results

            # open connection (append option) and write results
            of_connection = open(out_file, 'a')
            writer = csv.writer(of_connection)
            writer.writerow(eval_results)
            
            # make sure to close connection
            of_connection.close()
            
        # Sort with best score on top
        results.sort_values('sharpe', ascending = False, inplace = True)
        results.reset_index(inplace = True)

        return results 

    def random_search_prots(df, param_grid, out_file, max_evals):
        """Random search for hyperparameter optimization. 
        Writes result of search to csv file every search iteration."""
        
        
        # Dataframe for results
        results = pd.DataFrame(columns = ['sharpe', 'params', 'iteration'],
                                    index = list(range(max_evals)))
        for i in range(max_evals):

            print('%d/%d' % (i+1, max_evals))
            
            # Choose random hyperparameters
            random_params = {k: random.sample(v, 1)[0] for k, v in param_grid.items()}
            # random_params['subsample'] = 1.0 if random_params['boosting_type'] == 'goss' else random_params['subsample']

            # Evaluate randomly selected hyperparameters
            eval_results = objective_prots(df, random_params, i)
            results.loc[i, :] = eval_results

            # open connection (append option) and write results
            of_connection = open(out_file, 'a')
            writer = csv.writer(of_connection)
            writer.writerow(eval_results)
            
            # make sure to close connection
            of_connection.close()
            
        # Sort with best score on top
        results.sort_values('sharpe', ascending = False, inplace = True)
        results.reset_index(inplace = True)

        return results 

    def get_cat_model_pfms(df, hyperparameters, feature_columns, train, test):
        
    #     train = df[df['era'].isin(train_eras)]
    #     test = df[df['era'].isin(test_eras)]
        feature_columns = list(feature_columns)
        dtrain = lgb.Dataset(get_X_array(train, feature_columns), get_Y_array_pfms(train), feature_name=feature_columns, free_raw_data=False)
        dvalid = lgb.Dataset(get_X_array(test, feature_columns), get_Y_array_pfms(test), feature_name=feature_columns, free_raw_data=False)
        np.random.seed(SEED)
        evals_result = {}
        model = lgb.train(hyperparameters, dtrain, num_boost_round=20000, valid_sets=(dvalid,), valid_names=('valid',), verbose_eval=False,
                    early_stopping_rounds=100, evals_result=evals_result)
        return model

    def get_cat_model_prots(df, hyperparameters, feature_columns, train, test):
        
    #     train = df[df['era'].isin(train_eras)]
    #     test = df[df['era'].isin(test_eras)]
        feature_columns = list(feature_columns)
        dtrain = lgb.Dataset(get_X_array(train, feature_columns), get_Y_array_prots(train), feature_name=feature_columns, free_raw_data=False)
        dvalid = lgb.Dataset(get_X_array(test, feature_columns), get_Y_array_prots(test), feature_name=feature_columns, free_raw_data=False)
        np.random.seed(SEED)
        evals_result = {}
        model = lgb.train(hyperparameters, dtrain, num_boost_round=20000, valid_sets=(dvalid,), valid_names=('valid',), verbose_eval=False,
                    early_stopping_rounds=100, evals_result=evals_result)
        return model

    def get_cat_model_final_pfms(df, hyperparameters, feature_columns):
        feature_columns = list(feature_columns)
        train = df
        dtrain = lgb.Dataset(get_X_array(train, feature_columns), get_Y_array_pfms(train), feature_name=feature_columns, free_raw_data=False)
        np.random.seed(SEED)
        model = lgb.train(hyperparameters, dtrain)
        return model

    def get_cat_model_final_prots(df, hyperparameters, feature_columns):
        feature_columns = list(feature_columns)
        train = df
        dtrain = lgb.Dataset(get_X_array(train, feature_columns), get_Y_array_prots(train), feature_name=feature_columns, free_raw_data=False)
        np.random.seed(SEED)
        model = lgb.train(hyperparameters, dtrain)
        return model

    def get_cat_model_final_Is(df, hyperparameters, feature_columns):
        feature_columns = list(feature_columns)
        train = df
        dtrain = lgb.Dataset(get_X_array(train, feature_columns), get_Y_array_Is(train), feature_name=feature_columns, free_raw_data=False)
        np.random.seed(SEED)
        model = lgb.train(hyperparameters, dtrain)
        return model

    def get_cat_model_final_FAIMS(df, hyperparameters, feature_columns):
        feature_columns = list(feature_columns)
        train = df
        dtrain = lgb.Dataset(get_X_array(train, feature_columns), get_Y_array_FAIMS(train), feature_name=feature_columns, free_raw_data=False)
        np.random.seed(SEED)
        model = lgb.train(hyperparameters, dtrain)
        return model

    def get_Y_array_prots(df):
        return df.loc[:, 'decoy'].values

    def get_features_prots(dataframe):
        feature_columns = dataframe.columns
        columns_to_remove = []
        banned_features = {
            'dbname',
            'decoy',
            'preds',
        }
        for feature in feature_columns:
            # if feature in banned_features or not feature.startswith('sf_peptide'):
            if not feature.startswith('qr_'):
                columns_to_remove.append(feature)
        feature_columns = feature_columns.drop(columns_to_remove)
        return feature_columns

    TARGET_NAME = 'decoy'
    PREDICTION_NAME = 'preds'

    def objective_prots(df, hyperparameters, iteration):
        """Objective function for grid and random search. Returns
        the cross validation score from a set of hyperparameters."""
        
        all_res = []
        all_iters = []
        
        # Number of estimators will be found using early stopping
        if 'n_estimators' in hyperparameters.keys():
            del hyperparameters['n_estimators']
        
        feature_columns = get_features_prots(df)
        print(feature_columns)
        kf = KFold(n_splits=3, shuffle=True, random_state=SEED+1)
        
        for sp1 in kf.split(df):
            train_df = df.iloc[sp1[0]]
            test_df = df.iloc[sp1[1]]
            model = get_cat_model_prots(df, hyperparameters, feature_columns, train_df, test_df)
            all_iters.append(model.best_iteration)
            
            test_df[PREDICTION_NAME] = model.predict(get_X_array(test_df, feature_columns))
            # all_res.append(len(aux.filter(test_df, fdr=0.1, key='preds', is_decoy='decoy')))



            all_res.append(test_df)
        test_df6 = pd.concat(all_res)
        shr_v = len(aux.filter(test_df6, fdr=0.25, key='preds', is_decoy='decoy'))#np.min(all_res)


    #         all_res.append(test_df[[TARGET_NAME, PREDICTION_NAME]])
    #     test_df6 = pd.concat(all_res)
        
    #     fpr, tpr, thresholds = metrics.roc_curve(test_df6[TARGET_NAME], test_df6[PREDICTION_NAME])
    #     shr_v = metrics.auc(fpr, tpr)
        # shr_v = np.min(all_res)
        hyperparameters['n_estimators'] = int(np.mean(all_iters))
        
        return [shr_v, hyperparameters, iteration]


    df1 = pd.DataFrame()
    for k in resdict.keys():
        # if k == 'Is':
        #     df1[k] = np.log10(resdict[k])
        # else:
        df1[k] = resdict[k]
    df1['mass_diff'] = mass_diff
    df1['rt_diff'] = rt_diff
    df1['decoy'] = df1['seqs'].apply(lambda x: all(z.startswith(prefix) for z in pept_prot[x]))

    print(sum(df1['decoy']), sum(~df1['decoy']))

    # df1['ids_count'] = df1.groupby('ids')['seqs'].transform('count')
    df1['peptide'] = df1['seqs']

    # df1 = df1.set_index('peptide')

    # tmp = df1.groupby('peptide').apply(agg_func)
    # df1['im'] = tmp
    
    # df1 = df1.reset_index()


    all_dbnames = []
    all_theor_peptides = []
    for item in protsN.items():
        all_dbnames.append(item[0])
        all_theor_peptides.append(item[1])
        
        
    df2 = pd.DataFrame()
    df2['dbname'] = all_dbnames
    df2['theor peptides'] = all_theor_peptides
    df2 = df2[df2['theor peptides'] >= 5]

    df3 = df1.copy()

    df3['proteins'] = df3['seqs'].apply(lambda x: pept_prot[x])

    df3 = df3.assign(proteins=df3.proteins).explode('proteins').reset_index(drop=True)

    tmp = df3.groupby('proteins')['seqs'].agg(['count', ])
    tmp2 = df3.groupby('proteins')['seqs'].nunique()

    df2['decoy'] = df2['dbname'].apply(lambda x: x.startswith('DECOY_'))

    df2 = df2.set_index(keys='dbname')

    df2['PSMs_total'] = tmp['count']
    df2['peptides_total'] = tmp2

    p_psms = df2[df2['decoy']]['PSMs_total'].sum() / df2[df2['decoy']]['theor peptides'].sum()
    p_peptides = df2[df2['decoy']]['peptides_total'].sum() / df2[df2['decoy']]['theor peptides'].sum()
    df2['sf_psm_total'] = calc_sf_all_2(df2['PSMs_total'].values, df2['theor peptides'].values, p_psms)
    df2['sf_peptide_total'] = calc_sf_all_2(df2['peptides_total'].values, df2['theor peptides'].values, p_peptides)

    df2['dbname'] = df2.index

    protein_ranks = dict()
    cur_r = 1
    for z in df2.sort_values(by='sf_peptide_total', ascending=False)[['dbname']].values:
        protein_ranks[z[0]] = cur_r
        cur_r += 1
    worst_rank = cur_r + 1

    def get_protein_with_max_score(proteins):
        maxscore = 0
        best_prot = 'Unknown'
        for prot in proteins:
            curscore = protein_ranks.get(prot, 0)
            if curscore > maxscore:
                maxscore = curscore
                best_prot = prot
        return best_prot



    def recalc_rank(raw):
        dbn = raw['bestprotein']
        theor_peps = protsN.get(dbn, 1e6)
        return float(theor_peps + 1 - raw['IntensityNorm'])/theor_peps

    df1['proteins'] = df1['seqs'].apply(lambda x: pept_prot[x])
    df1['bestprotein'] = df1['proteins'].apply(get_protein_with_max_score)
    # df1["IntensityNorm"] = df1['Is'].values / df1.groupby("bestprotein")['Is'].transform(lambda x: x.median())
    df1['IntensityNorm'] = df1.groupby('bestprotein').Is.rank(method='average', ascending=False)
    df1['IntensityNorm'] = df1.apply(recalc_rank, axis=1)
    mass_dict = {}
    pI_dict = {}
    charge_dict = {}
    for pep in set(df1['peptide']):
        try:
            mass_dict[pep] = mass.fast_mass2(pep)
            pI_dict[pep] = electrochem.pI(pep)
            charge_dict[pep] = electrochem.charge(pep, pH=7.0)
        except:
            mass_dict[pep] = 0
            pI_dict[pep] = 0
            charge_dict[pep] = 0


    for aa in mass.std_aa_mass:
        df1['c_%s' % (aa, )] = df1['peptide'].apply(lambda x: x.count(aa))
    df1['c_DP'] = df1['peptide'].apply(lambda x: x.count('DP'))
    df1['c_KP'] = df1['peptide'].apply(lambda x: x.count('KP'))
    df1['c_RP'] = df1['peptide'].apply(lambda x: x.count('RP'))

    df1['plen'] = df1['peptide'].apply(lambda z: len(z))
    df1['mass'] = df1['peptide'].apply(lambda x: mass_dict[x])
    df1['pI'] = df1['peptide'].apply(lambda x: pI_dict[x])
    df1['charge_theor'] = df1['peptide'].apply(lambda x: charge_dict[x])

    df1['rt_diff_abs'] = df1['rt_diff'].abs()
    df1['rt_diff_abs_pdiff'] = df1['rt_diff_abs'] - df1.groupby('ids')['rt_diff_abs'].transform('median')
    df1['rt_diff_abs_pnorm'] = df1['rt_diff_abs'] / (df1.groupby('ids')['rt_diff_abs'].transform('sum') + 1e-2)
    df1['id_count'] = df1.groupby('ids')['mass_diff'].transform('count')
    df1['seq_count'] = df1.groupby('peptide')['mass_diff'].transform('count')

    df1t5 = df1.sort_values(by='Is', ascending=False).copy()
    df1t5 = df1t5.drop_duplicates(subset='peptide', keep='first')    
    # df1t5 = df1[df1.groupby('bestprotein').peptide.transform('count')>5]

    # df1t5 = df1t5.assign(proteins=df1t5.proteins).explode('proteins').reset_index(drop=True)
    # df1t5['IntensityNorm'] = df1t5.groupby('proteins').Intensity.rank(method='average', ascending=False)
    # df1t5['IntensityNorm'] = df1t5.apply(recalc_rank, axis=1)

    random_results = pd.read_csv('/home/mark/notebooks/test_randomCV_MS1Intensity2.csv')
    random_results = random_results[random_results['auc'] != 'auc']
    random_results['params'] = random_results['params'].apply(lambda x: ast.literal_eval(x))
    random_results['boosts'] = random_results['params'].apply(lambda x: int(x['n_estimators']))
    convert_dict = {'auc': float, 
                } 
    random_results = random_results.astype(convert_dict) 

    random_results['auc_per_boost'] = random_results['auc'].values / random_results['boosts'].values

    bestparams = random_results.sort_values(by='auc',ascending=True)['params'].values[0]
    # bestparams['n_estimators'] = int(bestparams['n_estimators'] * 1.5)
    bestparams['num_threads'] = 5




    # print('Start Machine Learning on Is...')

    # feature_columns_Is = get_features_Is(df1t5)

    # print(feature_columns_Is)

    # groups = df1['peptide']
    # ix = df1.index.values
    # unique = np.unique(groups)
    # np.random.RandomState(SEED).shuffle(unique)
    # result = []
    # for split in np.array_split(unique, 3):
    #     mask = groups.isin(split)
    #     train, test = ix[~mask], ix[mask]
    #     train_df = df1.iloc[train]
    #     test_df = df1.iloc[test]
    #     model = get_cat_model_final_Is(train_df, bestparams, feature_columns_Is)
        
    #     df1.loc[test, 'IntensityNorm_predicted'] = model.predict(get_X_array(test_df, feature_columns_Is))
    # df1['IntensityNorm_diff'] = df1['IntensityNorm_predicted'] - df1['IntensityNorm']



    print('Start Machine Learning on Is...')


    feature_columns_Is = get_features_Is(df1t5)

    kf = KFold(n_splits=3, shuffle=True, random_state=SEED)

    test_res = []
    for sp1 in kf.split(df1):
        train_df = df1.iloc[sp1[0]]
        test_df = df1.iloc[sp1[1]]

        model = get_cat_model_final_Is(train_df, bestparams, feature_columns_Is)
        
        df1.loc[sp1[1], 'IntensityNorm_predicted'] = model.predict(get_X_array(test_df, feature_columns_Is))

    # model = get_cat_model_final_Is(df1t5, bestparams, feature_columns_Is)
    # df1['IntensityNorm_predicted'] = model.predict(get_X_array(df1, feature_columns_Is))
    df1['IntensityNorm_diff'] = df1['IntensityNorm_predicted'] - df1['IntensityNorm']


    # print('Start calculation of spearman rank corr')

    # t_counter = 0

    # for dbname in set(df1['bestprotein']):
    #     t_counter += 1
    #     if t_counter % 100 == 0:
    #         print(t_counter)
    #     idx_mask = df1['bestprotein'] == dbname
    #     tmp = df1[idx_mask]
    #     ar1 = tmp['IntensityNorm'].values
    #     ar2 = rankdata(tmp['IntensityNorm_predicted'].values)

    #     ar3 = []

    #     mask = np.ones(ar1.shape,dtype=bool)
    #     cor_basic = spearmanr(ar1[mask], ar2[mask])[0]

    #     l_m = len(mask)
    #     for idx in range(l_m):
    #         if idx > 0:
    #             mask[idx-1] = 1
    #         mask[idx] = 0
    #         if l_m > 2:
    #             cor_v = spearmanr(ar1[mask], ar2[mask])[0]
    #         else:
    #             cor_v = cor_basic
    #         ar3.append(cor_v)

    #     # print(l_m, len(ar3), sum(idx_mask))

    #     df1.loc[idx_mask, 'cor_diff'] = np.array(ar3) - cor_basic
    #     df1.loc[idx_mask, 'cor_basic'] = cor_basic
    #     # prots_cor[dbname] = spearmanr(tmp['IntensityNorm'], tmp['IntensityNorm_predicted'])[0]

    # # idx = df1['proteins'] == 'DECOY_sp|O75146|HIP1R_HUMAN'

    # # # df1['plen'] = df1['seqs'].apply(lambda x: len(x))

    # print('Start Machine Learning on FAIMS...')


    # random_results = pd.read_csv('/home/mark/notebooks/test_randomCV_FAIMS.csv')
    # random_results = random_results[random_results['auc'] != 'auc']
    # random_results['params'] = random_results['params'].apply(lambda x: ast.literal_eval(x))
    # random_results['boosts'] = random_results['params'].apply(lambda x: int(x['n_estimators']))
    # convert_dict = {'auc': float, 
    #             } 
    # random_results = random_results.astype(convert_dict) 

    # random_results['auc_per_boost'] = random_results['auc'].values / random_results['boosts'].values

    # bestparams = random_results.sort_values(by='auc',ascending=True)['params'].values[0]
    # # bestparams['n_estimators'] = int(bestparams['n_estimators'] * 1.5)
    # bestparams['num_threads'] = 6

    # feature_columns_FAIMS = get_features_FAIMS(df1)

    # kf = KFold(n_splits=3, shuffle=False, random_state=SEED)

    # test_res = []
    # for sp1 in kf.split(df1):
    #     train_df = df1.iloc[sp1[0]]

    #     model = get_cat_model_final_FAIMS(train_df, bestparams, feature_columns_FAIMS)
        
    #     df1.loc[sp1[1], 'im_predicted'] = model.predict(get_X_array(df1.loc[sp1[1], :], feature_columns_FAIMS))

    # # df1_top = df1[df1['bestprotein'].apply(lambda x: protein_ranks.get(x, 1e6) <= 100)]

    # # model = get_cat_model_final_FAIMS(df1_top, bestparams, feature_columns_FAIMS)
    # # df1['im_predicted'] = model.predict(get_X_array(df1, feature_columns_FAIMS))
    # df1['im_diff'] = df1['im_predicted'] - df1['im']

    # # df1['plen'] = df1['seqs'].apply(lambda x: len(x))

    if args['ml']:

        print('Start Machine Learning on PFMs...')

        print(get_features_pfms(df1))

        MAX_EVALS = 25
        out_file = 'test_randomCV_PFMs_2.csv'
        of_connection = open(out_file, 'w')
        writer = csv.writer(of_connection)

        # Write column names
        headers = ['auc', 'params', 'iteration']
        writer.writerow(headers)
        of_connection.close()

        # df = df.reset_index(drop=True)

        random_results = random_search_pfms(df1, param_grid, out_file, MAX_EVALS)

        random_results = pd.read_csv(out_file)
        random_results = random_results[random_results['auc'] != 'auc']
        random_results['params'] = random_results['params'].apply(lambda x: ast.literal_eval(x))
        random_results['boosts'] = random_results['params'].apply(lambda x: int(x['n_estimators']))
        convert_dict = {'auc': float, 
                    } 
        random_results = random_results.astype(convert_dict) 

        random_results['auc_per_boost'] = random_results['auc'].values / random_results['boosts'].values

        bestparams = random_results.sort_values(by='auc',ascending=False)['params'].values[0]
        bestparams['num_threads'] = 5
        print(random_results.sort_values(by='auc',ascending=False)['auc'].values[0])



        groups = df1['peptide']
        # ix = np.arange(len(groups))
        ix = df1.index.values
        unique = np.unique(groups)
        np.random.RandomState(SEED).shuffle(unique)
        result = []
        for split in np.array_split(unique, 3):
            mask = groups.isin(split)
            train, test = ix[~mask], ix[mask]
            train_df = df1.iloc[train]
            test_df = df1.iloc[test]

        # kf = KFold(n_splits=3, shuffle=True, random_state=SEED)

        # test_res = []
        # for sp1 in kf.split(df1):
        #     train_df = df1.iloc[sp1[0]]
        #     test_df = df1.iloc[sp1[1]]

            feature_columns = list(get_features_pfms(train_df))
            model = get_cat_model_final_pfms(train_df, bestparams, feature_columns)
            
            df1.loc[test, 'preds'] = model.predict(get_X_array(test_df, feature_columns))

    else:
        df1['preds'] = np.power(df1['mass_diff'], 2) + np.power(df1['rt_diff'], 2)

    # df1['pc'] = df1.groupby('ids')['preds'].transform('count')
    # df1.loc[df1['pc'] > 1, 'preds'] = df1.loc[df1['pc'] > 1, 'preds'] - df1.loc[df1['pc'] > 1, :].groupby('ids')['preds'].transform('median')
        
    df1['qpreds'] = pd.qcut(df1['preds'], 10, labels=range(10)) 
    df1['proteins'] = df1['seqs'].apply(lambda x: ';'.join(pept_prot[x]))

    df1.to_csv(base_out_name + '_PFMs_ML.csv', sep='\t', index=False)
     
        
    # all_dbnames = []
    # all_theor_peptides = []
    # for item in protsN.items():
    #     all_dbnames.append(item[0])
    #     all_theor_peptides.append(item[1])
        
        
    # df2 = pd.DataFrame()
    # df2['dbname'] = all_dbnames
    # df2['theor peptides'] = all_theor_peptides
    # df2 = df2[df2['theor peptides'] >= 5]

    # df3 = df1.copy()

    # df3['proteins'] = df3['seqs'].apply(lambda x: pept_prot[x])

    # df3 = df3.assign(proteins=df3.proteins).explode('proteins').reset_index(drop=True)

    # tmp = df3.groupby('proteins')['preds'].agg(['mean', 'count'])
    # tmp2 = df3.groupby('proteins')['seqs'].nunique()

    # df2['decoy'] = df2['dbname'].apply(lambda x: x.startswith('DECOY_'))

    # df2 = df2.set_index(keys='dbname')



    # df2['av_score_total'] = tmp['mean']
    # df2['PSMs_total'] = tmp['count']
    # df2['peptides_total'] = tmp2

    # p_psms = df2[df2['decoy']]['PSMs_total'].sum() / df2[df2['decoy']]['theor peptides'].sum()
    # p_peptides = df2[df2['decoy']]['peptides_total'].sum() / df2[df2['decoy']]['theor peptides'].sum()
    # df2['sf_psm_total'] = calc_sf_all_2(df2['PSMs_total'].values, df2['theor peptides'].values, p_psms)
    # df2['sf_peptide_total'] = calc_sf_all_2(df2['peptides_total'].values, df2['theor peptides'].values, p_peptides)


    # df2['dbname'] = df2.index

    # protein_ranks = dict()
    # cur_r = 1
    # for z in df2.sort_values(by='sf_peptide_total', ascending=False)[['dbname']].values:
    #     protein_ranks[z[0]] = cur_r
    #     cur_r += 1
    # worst_rank = cur_r + 1

        
    # df2 = pd.DataFrame()
    # df2['dbname'] = all_dbnames
    # df2['theor peptides'] = all_theor_peptides
    # df2 = df2[df2['theor peptides'] >= 5]

    # df3 = df1.copy()

    # df3['best_prot_rank'] = df3['seqs'].apply(lambda x: min(protein_ranks.get(z, worst_rank) for z in pept_prot[x]))

    # df3 = df3.sort_values(by='best_prot_rank', ascending=True)
    # df3 = df3.drop_duplicates(subset=['ids'])
    # df3 = df3.reset_index()
    # df3 = df3.drop(columns='index')

    # df3['proteins'] = df3['seqs'].apply(lambda x: pept_prot[x])

    # df3 = df3.assign(proteins=df3.proteins).explode('proteins').reset_index(drop=True)

    # tmp = df3.groupby('proteins')['preds'].agg(['mean', 'count'])
    # tmp2 = df3.groupby('proteins')['seqs'].nunique()

    # df2['decoy'] = df2['dbname'].apply(lambda x: x.startswith('DECOY_'))

    # df2 = df2.set_index(keys='dbname')

    # # from scipy.stats import binom
    # # def calc_sf_all(v, n, p):
    # #     sf_values = -np.log10(binom.sf(v, n, p))
    # #     sf_values[np.isinf(sf_values)] = max(sf_values[~np.isinf(sf_values)])
    # #     return sf_values

    # df2['av_score_total'] = tmp['mean']
    # df2['PSMs_total'] = tmp['count']
    # df2['peptides_total'] = tmp2

    # p_psms = df2[df2['decoy']]['PSMs_total'].sum() / df2[df2['decoy']]['theor peptides'].sum()
    # p_peptides = df2[df2['decoy']]['peptides_total'].sum() / df2[df2['decoy']]['theor peptides'].sum()
    # df2['sf_psm_total'] = calc_sf_all_2(df2['PSMs_total'].values, df2['theor peptides'].values, p_psms)
    # df2['sf_peptide_total'] = calc_sf_all_2(df2['peptides_total'].values, df2['theor peptides'].values, p_peptides)


    # for cc in df2.columns:
    #     if cc.startswith('PSMs_'):
    #         df2[cc] = df2[cc].fillna(value=0)
    #     elif cc.startswith('peptides_'):
    #         df2[cc] = df2[cc].fillna(value=0)
    #     elif cc.startswith('av_score'):
    #         df2[cc] = df2[cc].fillna(value=df2[cc].max())
    #     elif cc.startswith('sf_'):
    #         df2[cc] = df2[cc].fillna(value=0)

    # df2['sf_peptide_total'] = zscore(df2['sf_peptide_total'])

    # for qval in range(10):
    #     qqq = df3[df3['qpreds'] <= qval]
    #     tmp = qqq.groupby('proteins')['preds'].agg(['mean', 'count'])
    #     tmp2 = qqq.groupby('proteins')['seqs'].nunique()
        
    #     df2['PSMs_%s' % (qval, )] = tmp['count']
    #     df2['peptides_%s' % (qval, )] = tmp2
    #     df2['av_score_%s' % (qval, )] = tmp['mean']
        
    #     for cc in df2.columns:
    #         if cc.startswith('PSMs_'):
    #             df2[cc] = df2[cc].fillna(value=0)
    #         elif cc.startswith('peptides_'):
    #             df2[cc] = df2[cc].fillna(value=0)
    #         elif cc.startswith('av_score'):
    #             df2[cc] = df2[cc].fillna(value=df2[cc].max())
    #         elif cc.startswith('sf_'):
    #             df2[cc] = df2[cc].fillna(value=0)
        
    #     p_psms = df2[df2['decoy']]['PSMs_%s' % (qval, )].sum() / df2[df2['decoy']]['theor peptides'].sum()
    #     p_peptides = df2[df2['decoy']]['peptides_%s' % (qval, )].sum() / df2[df2['decoy']]['theor peptides'].sum()
    #     df2['sf_psm_%s' % (qval, )] = calc_sf_all_2(df2['PSMs_%s' % (qval, )].values, df2['theor peptides'].values, p_psms)
    #     df2['sf_peptide_%s' % (qval, )] = calc_sf_all_2(df2['peptides_%s' % (qval, )].values, df2['theor peptides'].values, p_peptides)
    #     df2['sf_peptide_%s' % (qval, )] = zscore(df2['sf_peptide_%s' % (qval, )])
        
    # for cc in df2.columns:
    #     if cc.startswith('PSMs_'):
    #         df2[cc] = df2[cc].fillna(value=0)
    #     elif cc.startswith('peptides_'):
    #         df2[cc] = df2[cc].fillna(value=0)
    #     elif cc.startswith('av_score'):
    #         df2[cc] = df2[cc].fillna(value=df2[cc].max())
    #     elif cc.startswith('sf_'):
    #         df2[cc] = df2[cc].fillna(value=0)
            
    # df2['dbname'] = df2.index

    # df2 = df2.reset_index(drop=True)

    # from sklearn.preprocessing import quantile_transform
    # feature_columns = list(get_features_prots(df2))
    # # for fc in feature_columns:
    # #     df2[fc] = quantile_transform(df2[[fc]], n_quantiles=300, random_state=0)



    # # Hyperparameter grid
    # param_grid_prots = {
    # #     'boosting_type': ['gbdt', 'goss', 'dart'],
    # #     'boosting_type': ['gbdt', 'goss'],
    #     'boosting_type': ['gbdt', ],
    # #     'boosting_type': ['dart', ],
    #     'num_leaves': list(range(10, 1000)),
    #     'learning_rate': list(np.logspace(np.log10(0.01), np.log10(0.3), base = 10, num = 1000)),
    #     'subsample_for_bin': list(range(1, 500, 5)),
    #     'min_child_samples': list(range(1, 150, 1)),
    #     'reg_alpha': list(np.linspace(0, 1)),
    #     'reg_lambda': list(np.linspace(0, 1)),
    #     'colsample_bytree': list(np.linspace(0.01, 1, 100)),
    #     'subsample': list(np.linspace(0.01, 1, 100)),
    #     'is_unbalance': [True, False],
    #     'metric': ['rmse', ],
    #     'verbose': [-1, ],
    #     'num_threads': [5, ],
    # }

    # print('Start Machine Learning on proteins...')
    # MAX_EVALS = 1
    # out_file = 'test_randomCV_proteins2.csv'
    # of_connection = open(out_file, 'w')
    # writer = csv.writer(of_connection)

    # # Write column names
    # headers = ['auc', 'params', 'iteration']
    # writer.writerow(headers)
    # of_connection.close()

    # # df = df.reset_index(drop=True)

    # random_results = random_search_prots(df2, param_grid_prots, out_file, MAX_EVALS)
    # random_results = pd.read_csv('test_randomCV_proteins2.csv')
    # random_results = random_results[random_results['auc'] != 'auc']
    # random_results['params'] = random_results['params'].apply(lambda x: ast.literal_eval(x))
    # random_results['boosts'] = random_results['params'].apply(lambda x: int(x['n_estimators']))
    # convert_dict = {'auc': float, 
    #             } 
    # random_results = random_results.astype(convert_dict) 

    # random_results['auc_per_boost'] = random_results['auc'].values / random_results['boosts'].values

    # bestparams = random_results.sort_values(by='auc',ascending=False)['params'].values[0]
    # bestparams['num_threads'] = 6
    # best_auc = random_results.sort_values(by='auc',ascending=False)['auc'].values[0]
    # print(best_auc)

    # if best_auc >= 100000:

    #     kf = KFold(n_splits=3, shuffle=True, random_state=SEED)

    #     test_res = []
    #     for sp1 in kf.split(df2):
    #         train_df = df2.iloc[sp1[0]]
    #         test_df = df2.iloc[sp1[1]]

    #         feature_columns = list(get_features_prots(train_df))
    #         model = get_cat_model_final_prots(train_df, bestparams, feature_columns)
            
    #         df2.loc[sp1[1], 'preds'] = model.predict(get_X_array(test_df, feature_columns))

        
    #     df2['preds'] = df2['preds'].max() - df2['preds']

    # else:
    #     print('Skipping ML for proteins, use simple binomial scores')
    #     cols_for_sum = []
    #     for cc in df2.columns:
    #         if cc.startswith('sf_peptide'):
    #             cols_for_sum.append(cc)
    #     df2['preds'] = df2[cols_for_sum].sum(axis=1)


    # df2 = df2.rename({'preds': 'score', 'peptides_total': 'matched peptides', 'theor peptides': 'theoretical peptides'}, axis='columns')
    # df2.to_csv(base_out_name + '_proteins_full.csv', sep='\t', index=False, columns=['dbname', 'score', 'matched peptides', 'theoretical peptides'])

    # df2['basedbname'] = df2['dbname'].apply(lambda x: x.replace('DECOY_', ''))
    # df2 = df2.sort_values(by='score', ascending=False)
    # df2 = df2.drop_duplicates(subset=['basedbname'])

    # df2f = aux.filter(df2, fdr=fdr, key='score', is_decoy='decoy', reverse=True)
    # df2f.to_csv(base_out_name + '_proteins.csv', sep='\t', index=False, columns=['dbname', 'score', 'matched peptides', 'theoretical peptides'])


    resdict['qpreds'] = df1['qpreds'].values
    resdict['IntensityNorm'] = df1['IntensityNorm'].values
    resdict['IntensityNorm_predicted'] = df1['IntensityNorm_predicted'].values
    mass_diff = resdict['qpreds']
    rt_diff = resdict['qpreds']

    pickle.dump(resdict, open('/home/mark/resdict.pickle', 'wb'))

    p1 = set(resdict['seqs'])

    prots_spc2 = defaultdict(set)
    for pep, proteins in pept_prot.items():
        if pep in p1:
            for protein in proteins:
                prots_spc2[protein].add(pep)

    for k in protsN:
        if k not in prots_spc2:
            prots_spc2[k] = set([])
    prots_spc = dict((k, len(v)) for k, v in prots_spc2.items())

    names_arr = np.array(list(prots_spc.keys()))
    v_arr = np.array(list(prots_spc.values()))
    n_arr = np.array([protsN[k] for k in prots_spc])

    top100decoy_score = [prots_spc.get(dprot, 0) for dprot in protsN if isdecoy_key(dprot)]
    top100decoy_N = [val for key, val in protsN.items() if isdecoy_key(key)]
    p = np.mean(top100decoy_score) / np.mean(top100decoy_N)
    print('p=%s' % (np.mean(top100decoy_score) / np.mean(top100decoy_N)))

    prots_spc = dict()
    all_pvals = calc_sf_all(v_arr, n_arr, p)
    for idx, k in enumerate(names_arr):
        prots_spc[k] = all_pvals[idx]

    # checked = set()
    # for k, v in list(prots_spc.items()):
    #     if k not in checked:
    #         if isdecoy_key(k):
    #             if prots_spc.get(k.replace(prefix, ''), -1e6) > v:
    #                 del prots_spc[k]
    #                 checked.add(k.replace(prefix, ''))
    #         else:
    #             if prots_spc.get(prefix + k, -1e6) > v:
    #                 del prots_spc[k]
    #                 checked.add(prefix + k)

    # filtered_prots = aux.qvalues(prots_spc.items(), key=escore, is_decoy=isdecoy, remove_decoy=False, formula=1,
    #                             full_output=True)

    # pqmap = dict()
    # for z in filtered_prots:
    #     qval = z[2]
    #     pdbname = z[3][0]
    #     pqmap[pdbname] = qval
    # #     print(qval, pdbname)
    # #     break

    # df1 = pd.DataFrame()
    # for k in resdict.keys():
    #     if k == 'Is':
    #         df1[k] = np.log10(resdict[k])
    #     else:
    #         df1[k] = resdict[k]
    # df1['mass_diff'] = mass_diff
    # df1['rt_diff'] = rt_diff
    # df1['decoy'] = df1['seqs'].apply(lambda x: all(z.startswith(prefix) for z in pept_prot[x]))

    # df1['peak_rank_md'] = df1.groupby("ids")["mass_diff"].rank("dense", ascending=True)
    # df1['peak_rank_rd'] = df1.groupby("ids")["rt_diff"].rank("dense", ascending=True)

    # def get_X_array(df, feature_columns):
    #     return df.loc[:, feature_columns].values

    # def get_Y_array(df):
    #     return df.loc[:, 'decoy'].values

    # def get_features(dataframe):
    #     feature_columns = dataframe.columns
    #     columns_to_remove = []
    #     banned_features = {
    #         'ids',
    #         'seqs',
    #         'decoy',
    #         'preds',
    #     }
    #     for feature in feature_columns:
    #         if feature in banned_features:
    #             columns_to_remove.append(feature)
    #     feature_columns = feature_columns.drop(columns_to_remove)
    #     return feature_columns

    # bestparams = {'boosting_type': 'gbdt',
    # 'num_leaves': 850,
    # 'learning_rate': 0.01534511205602959,
    # #  'subsample_for_bin': 6910,
    # #  'min_child_samples': 306,
    # 'reg_alpha': 0.26530612244897955,
    # 'reg_lambda': 0.24489795918367346,
    # 'colsample_bytree': 0.54,
    # 'subsample': 0.3,
    # 'is_unbalance': False,
    # 'metric': 'rmse',
    # 'n_estimators': 289}

    # SEED = 50

    # def get_cat_model_final(df, hyperparameters, feature_columns):
    #     feature_columns = list(feature_columns)
    #     train = df
    #     dtrain = lgb.Dataset(get_X_array(train, feature_columns), get_Y_array(train), feature_name=feature_columns, free_raw_data=False)
    #     np.random.seed(SEED)
    #     model = lgb.train(hyperparameters, dtrain)
    #     return model

    # # kf = KFold(n_splits=3, shuffle=True, random_state=SEED)

    # # test_res = []
    # # for sp1 in kf.split(df1):
    # #     train_df = df1.iloc[sp1[0]]
    # #     test_df = df1.iloc[sp1[1]]

    # #     feature_columns = list(get_features(train_df))
    # #     model = get_cat_model_final(train_df, bestparams, feature_columns)
        
    # #     df1.loc[sp1[1], 'preds'] = model.predict(get_X_array(test_df, feature_columns))
    # #     test_df['searchScore'] = model.predict(get_X_array(test_df, feature_columns))

    # mass_diff = df1['preds'].values


    # i_diff = resdict['Is']


    final_iteration(resdict, mass_diff, rt_diff, pept_prot, protsN, base_out_name, prefix, isdecoy, isdecoy_key, escore, fdr, args['nproc'], fname)


def worker(qin, qout, mass_diff, rt_diff, resdict, protsN, pept_prot, isdecoy_key, isdecoy, fdr, prots_spc_basic2, win_sys=False):

    for item in (iter(qin.get, None) if not win_sys else qin):
        print(item)
        mass_koef, rtt_koef = item

        # # rtt_koef = 1.0

        # m_k = scoreatpercentile(mass_diff, mass_koef * 100)
        # e_ind = mass_diff <= m_k
        # resdict2 = filter_results(resdict, e_ind)
        # rt_diff2 = rt_diff[e_ind]
        # # i_diff2 = i_diff[e_ind]

        # # r_k = scoreatpercentile(rt_diff2, rtt_koef * 100)
        # # e_ind = rt_diff2 <= r_k
        # # resdict2 = filter_results(resdict2, e_ind)
        # # i_diff3 = i_diff2[e_ind]

        # # r_k = scoreatpercentile(-i_diff3, i_koef * 100)
        # # e_ind = -i_diff3 <= r_k
        # # resdict2 = filter_results(resdict2, e_ind)



        # # m_k = scoreatpercentile(mass_diff, mass_koef * 100)
        # # r_k = scoreatpercentile(rt_diff, rtt_koef * 100)
        # # i_k = scoreatpercentile(-i_diff, i_koef * 100)
        # # e_ind1 = mass_diff <= m_k
        # # e_ind2 = rt_diff <= r_k
        # # e_ind3 = -i_diff <= i_k
        # # e_ind = e_ind1 & e_ind2 & e_ind3

        # # print(len(e_ind), sum(e_ind))

        # # resdict2 = filter_results(resdict, e_ind)



        # # resdict2 = filter_results(resdict, e_ind)
        # # rt_diff2 = rt_diff[e_ind]
        # # i_diff2 = i_diff[e_ind]

        # r_k = scoreatpercentile(rt_diff2, rtt_koef * 100)
        # e_ind = rt_diff2 <= r_k
        # resdict2 = filter_results(resdict2, e_ind)
        # # i_diff3 = i_diff2[e_ind]

        # # r_k = scoreatpercentile(-i_diff3, i_koef * 100)
        # # e_ind = -i_diff3 <= r_k
        # # resdict2 = filter_results(resdict2, e_ind)


        # m_k = scoreatpercentile(mass_diff, mass_koef * 100)
        e_ind = mass_diff <= mass_koef
        resdict2 = filter_results(resdict, e_ind)


        prots_spc_final2 = {}
        prots_spc_IntensityNorm = defaultdict(list)
        prots_spc_IntensityNorm_predicted = defaultdict(list)

        for pep, Inorm, Inorm_predicted in zip(resdict2['seqs'], resdict2['IntensityNorm'], resdict2['IntensityNorm_predicted']):
            for bprot in pept_prot[pep]:
                prots_spc_IntensityNorm[bprot].append(Inorm)
                prots_spc_IntensityNorm_predicted[bprot].append(Inorm_predicted)
        for bprot in prots_spc_IntensityNorm:
            cor_basic = spearmanr(prots_spc_IntensityNorm[bprot], prots_spc_IntensityNorm_predicted[bprot])[0]
            prots_spc_final2[bprot] = cor_basic

        features_dict = dict()
        for pep in set(resdict2['seqs']):
            for bprot in pept_prot[pep]:
                prot_score = prots_spc_basic2[bprot]
                if prot_score > features_dict.get(pep, [-1, ])[-1]:
                    features_dict[pep] = (bprot, prot_score)

        prots_spc_basic = dict()

        p1 = set(resdict2['seqs'])

        # seqs_ar = resdict2['seqs']
        # ids_ar = resdict2['ids']
        # RT_ar = resdict2['RT']
        # idx = np.argsort(RT_ar)
        # seqs_ar = seqs_ar[idx]

        pep_pid = defaultdict(set)
        pid_pep = defaultdict(set)
        banned_dict = dict()
        for pep, pid in zip(resdict2['seqs'], resdict2['ids']):
            pep_pid[pep].add(pid)
            pid_pep[pid].add(pep)
            if pep in banned_dict:
                banned_dict[pep] += 1
            else:
                banned_dict[pep] = 1

        if len(p1):
            prots_spc_final = dict()
            prots_spc_copy = False
            prots_spc2 = False
            unstable_prots = set()
            p0 = False

            names_arr = False
            tmp_spc_new = False
            decoy_set = False

            while 1:
                if not prots_spc2:

                    best_match_dict = dict()
                    n_map_dict = defaultdict(list)
                    for k, v in protsN.items():
                        n_map_dict[v].append(k)

                    decoy_set = set()
                    for k in protsN:
                        if isdecoy_key(k):
                            decoy_set.add(k)
                    decoy_set = list(decoy_set)
                    

                    prots_spc2 = defaultdict(set)
                    for pep, proteins in pept_prot.items():
                        if pep in p1:
                            for protein in proteins:
                                if protein == features_dict[pep][0]:
                                    prots_spc2[protein].add(pep)

                    for k in protsN:
                        if k not in prots_spc2:
                            prots_spc2[k] = set([])
                    prots_spc2 = dict(prots_spc2)
                    unstable_prots = set(prots_spc2.keys())

                    top100decoy_N = sum([val for key, val in protsN.items() if isdecoy_key(key)])

                    names_arr = np.array(list(prots_spc2.keys()))
                    n_arr = np.array([protsN[k] for k in names_arr])

                    tmp_spc_new = dict((k, len(v)) for k, v in prots_spc2.items())


                    top100decoy_score_tmp = [tmp_spc_new.get(dprot, 0) for dprot in decoy_set]
                    top100decoy_score_tmp_sum = float(sum(top100decoy_score_tmp))

                tmp_spc = tmp_spc_new
                prots_spc = tmp_spc_new
                if not prots_spc_copy:
                    prots_spc_copy = deepcopy(prots_spc)

                for idx, v in enumerate(decoy_set):
                    if v in unstable_prots:
                        top100decoy_score_tmp_sum -= top100decoy_score_tmp[idx]
                        top100decoy_score_tmp[idx] = prots_spc.get(v, 0)
                        top100decoy_score_tmp_sum += top100decoy_score_tmp[idx]
                p = float(sum(top100decoy_score_tmp)) / top100decoy_N
                p = top100decoy_score_tmp_sum / top100decoy_N
                if not p0:
                    p0 = float(p)

                n_change = set(protsN[k] for k in unstable_prots)
                for n_val in n_change:
                    for k in n_map_dict[n_val]:
                        v = prots_spc[k]
                        if n_val not in best_match_dict or v > prots_spc[best_match_dict[n_val]]:
                            best_match_dict[n_val] = k
                n_arr_small = []
                names_arr_small = []
                v_arr_small = []
                for k, v in best_match_dict.items():
                    n_arr_small.append(k)
                    names_arr_small.append(v)
                    v_arr_small.append(prots_spc[v])

                prots_spc_basic = dict()
                all_pvals = calc_sf_all(v_arr_small, n_arr_small, p)
                for idx, k in enumerate(names_arr_small):
                    prots_spc_basic[k] = all_pvals[idx]

                best_prot = utils.keywithmaxval(prots_spc_basic)

                best_score = prots_spc_basic[best_prot]
                unstable_prots = set()
                if best_prot not in prots_spc_final:
                    prots_spc_final[best_prot] = best_score
                    banned_pids = set()
                    for pep in prots_spc2[best_prot]:
                        for pid in pep_pid[pep]:
                            banned_pids.add(pid)
                    for pid in banned_pids:
                        for pep in pid_pep[pid]:
                            banned_dict[pep] -= 1
                            if banned_dict[pep] == 0:
                                best_prot_val = features_dict[pep][0]
                                for bprot in pept_prot[pep]:
                                    if bprot == best_prot_val:
                                        tmp_spc_new[bprot] -= 1
                                        unstable_prots.add(bprot)
                else:

                    v_arr = np.array([prots_spc[k] for k in names_arr])
                    all_pvals = calc_sf_all(v_arr, n_arr, p)
                    for idx, k in enumerate(names_arr):
                        prots_spc_basic[k] = all_pvals[idx]

                    for k, v in prots_spc_basic.items():
                        if k not in prots_spc_final:
                            prots_spc_final[k] = v

                    break

                try:
                    prot_fdr = aux.fdr(prots_spc_final.items(), is_decoy=isdecoy)
                except ZeroDivisionError:
                    prot_fdr = 100.0
                if prot_fdr >= 12.5 * fdr:

                    v_arr = np.array([prots_spc[k] for k in names_arr])
                    all_pvals = calc_sf_all(v_arr, n_arr, p)
                    for idx, k in enumerate(names_arr):
                        prots_spc_basic[k] = all_pvals[idx]

                    for k, v in prots_spc_basic.items():
                        if k not in prots_spc_final:
                            prots_spc_final[k] = v
                    break

        keys_array = prots_spc_final.keys()
        values_array = np.array(list(prots_spc_final.values()))
        # values_array = zscore(values_array)
        for k, val in zip(keys_array, values_array):
            prots_spc_final[k] = [val, prots_spc_final2.get(k, -5)]

        # if mass_koef + rtt_koef >= 1.99:
        if mass_koef == 9:
            item2 = prots_spc_copy
        else:
            item2 = False
        if not win_sys:
            qout.put((prots_spc_final, item2))
        else:
            qout.append((prots_spc_final, item2))
    if not win_sys:
        qout.put(None)
    else:
        return qout
