import os
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
import lightgbm as lgb
import pandas as pd
from sklearn.model_selection import train_test_split
from scipy.stats import zscore, spearmanr
import pandas as pd
from pyteomics import pepxml, achrom, auxiliary as aux, mass, fasta, mzid, parser
from pyteomics import electrochem
import numpy as np
import random
SEED = 42
from sklearn.model_selection import train_test_split
from os import path, mkdir
from collections import Counter, defaultdict
import warnings
import pylab as plt
warnings.formatwarning = lambda msg, *args, **kw: str(msg) + '\n'

import pandas as pd
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
import matplotlib.pyplot as plt

from sklearn import (
    feature_extraction, feature_selection, decomposition, linear_model,
    model_selection, metrics, svm
)

import scipy
from scipy.stats import rankdata
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

SEED = 50

def get_cat_model(df, hyperparameters, feature_columns, train, test):
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
        outtrain_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
        outtrain = open(outtrain_name, 'w')
        outres_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
        for seq, RT in zip(ns, nr):
            outtrain.write(seq + '\t' + str(RT) + '\n')
        outtrain.close()

        outtest_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
        outtest = open(outtest_name, 'w')

        maxval = len(qin)
        start = 0
        while start + shift < maxval:
            item = qin[start+shift]
            outtest.write(item + '\n')
            start += step
        outtest.close()

        subprocess.call([elude_path, '-t', outtrain_name, '-e', outtest_name, '-a', '-o', outres_name])
        for x in open(outres_name).readlines()[3:]:
            seq, RT = x.strip().split('\t')
            pepdict[seq] = float(RT)
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

            if not p0:
                p0 = float(p)

                prots_spc_tmp = dict()
                v_arr = np.array([prots_spc[k] for k in names_arr])
                all_pvals = calc_sf_all(v_arr, n_arr, p)
                for idx, k in enumerate(names_arr):
                    prots_spc_tmp[k] = all_pvals[idx]

                sortedlist_spc = sorted(prots_spc_tmp.items(), key=operator.itemgetter(1))[::-1]
                with open(base_out_name + '_proteins_full_noexclusion.csv', 'w') as output:
                    output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\n')
                    for x in sortedlist_spc:
                        output.write('\t'.join((x[0], str(x[1]), str(prots_spc_copy[x[0]]), str(protsN[x[0]]))) + '\n')

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
        qout = worker(qin, qout, mass_diff, rt_diff, resdict, protsN, pept_prot, isdecoy_key, isdecoy, fdr, prots_spc_basic2, True)
        
        for item, item2 in qout:
            if item2:
                prots_spc_copy = item2
            for k in protsN:
                if k not in prots_spc_final:
                    prots_spc_final[k] = [item.get(k, 0.0), ]
                else:
                    prots_spc_final[k].append(item.get(k, 0.0))

    else:
        qin = Queue()
        qout = Queue()

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
                    if k not in prots_spc_final:
                        prots_spc_final[k] = [item.get(k, 0.0), ]
                    else:
                        prots_spc_final[k].append(item.get(k, 0.0))

        for p in procs:
            p.join()

    for k in prots_spc_final.keys():
        prots_spc_final[k] = np.mean(prots_spc_final[k])

    prots_spc = prots_spc_final
    sortedlist_spc = sorted(prots_spc.items(), key=operator.itemgetter(1))[::-1]
    with open(base_out_name + '_proteins_full.csv', 'w') as output:
        output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\n')
        for x in sortedlist_spc:
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
        print('\t'.join((str(x[0]), str(x[1]), str(int(prots_spc_copy[x[0]])), str(protsN[x[0]]))))
    print('results:%s;number of identified proteins = %d' % (base_out_name, identified_proteins, ))
    # print('R=', r)
    with open(base_out_name + '_proteins.csv', 'w') as output:
        output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\n')
        for x in filtered_prots:
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

def process_file(args):
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
        # massdiff3 = (m - imraw[i]) / m * 1e6
        # results.append((seqm, massdiff, rts[i], peak_id, I, Scans[i], Isotopes[i], mzraw[i], avraw[i], charges[i], imraw[i], massdiff3))
        results.append((seqm, massdiff, i))#rts[i], peak_id, I, Scans[i], Isotopes[i], mzraw[i], avraw[i], charges[i], imraw[i], massdiff3))
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

    print('Number of peptide isotopic clusters: %d' % (len(nmasses), ))

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
        'iorig',
        # 'rt',
        # 'ids',
        # 'Is',
        # 'Scans',
        # 'Isotopes',
        # 'mzraw',
        # 'av',
        # 'ch',
        # 'im',
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

    # resdict['mc'] = np.array([parser.num_sites(z, args['enzyme']) for z in resdict['seqs']])
    # resdict['mc'] = resdict['im']

    isdecoy = lambda x: x[0].startswith(prefix)
    isdecoy_key = lambda x: x.startswith(prefix)
    escore = lambda x: -x[1]

    # e_ind = resdict['mc'] == 0
    # resdict2 = filter_results(resdict, e_ind)

    e_ind = np.array([Isotopes[iorig] for iorig in resdict['iorig']]) >= min_isotopes_calibration
    # e_ind = resdict['Isotopes'] >= min_isotopes_calibration
    # e_ind = resdict['Isotopes'] >= 1
    resdict2 = filter_results(resdict, e_ind)

    # e_ind = resdict2['mc'] == 0
    # resdict2 = filter_results(resdict2, e_ind)

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

        identified_proteins = 0

        for x in filtered_prots:
            identified_proteins += 1
        print('results for default search: number of identified proteins = %d' % (identified_proteins, ))

        print('Running mass recalibration...')


        # e_ind = resdict['mc'] == 0
        # resdict2 = filter_results(resdict, e_ind)
        resdict2 = resdict

        true_md = []
        true_isotopes = []
        true_seqs = []
        true_prots = set(x[0] for x in filtered_prots)
        for pep, proteins in pept_prot.items():
            if any(protein in true_prots for protein in proteins):
                true_seqs.append(pep)

        e_ind = np.in1d(resdict2['seqs'], true_seqs)

        true_seqs = resdict2['seqs'][e_ind]
        true_md.extend(resdict2['md'][e_ind])
        true_md = np.array(true_md)
        # true_isotopes.extend(resdict2['Isotopes'][e_ind])
        true_isotopes.extend(np.array([Isotopes[iorig] for iorig in resdict2['iorig']])[e_ind])
        true_isotopes = np.array(true_isotopes)
        true_intensities = np.array([Is[iorig] for iorig in resdict2['iorig']])[e_ind]
        # true_intensities = np.array(resdict2['Is'][e_ind])
        # true_rt = np.array(resdict2['rt'][e_ind])
        # true_mz = np.array(resdict2['mzraw'][e_ind])
        true_rt = np.array([rts[iorig] for iorig in resdict2['iorig']])[e_ind]
        true_mz = np.array([mzraw[iorig] for iorig in resdict2['iorig']])[e_ind]

        df1 = pd.DataFrame()
        df1['mass diff'] = true_md
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

        # e_ind = resdict['mc'] == 0
        # resdict2 = filter_results(resdict, e_ind)

        # e_ind = resdict2['Isotopes'] >= min_isotopes_calibration
        e_ind = np.array([Isotopes[iorig] for iorig in resdict['iorig']]) >= min_isotopes_calibration
        resdict2 = filter_results(resdict, e_ind)

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


        e_ind = np.array([Isotopes[iorig] for iorig in resdict['iorig']]) >= min_isotopes_calibration
        # e_ind = resdict['Isotopes'] >= min_isotopes_calibration
        # e_ind = resdict['Isotopes'] >= 1
        resdict2 = filter_results(resdict, e_ind)

        # e_ind = resdict2['mc'] == 0
        # resdict2 = filter_results(resdict2, e_ind)


        true_seqs = []
        true_rt = []
        true_isotopes = []
        true_prots = set(x[0] for x in filtered_prots)#[:5])
        for pep, proteins in pept_prot.items():
            if any(protein in true_prots for protein in proteins):
                true_seqs.append(pep)
        e_ind = np.in1d(resdict2['seqs'], true_seqs)


        true_seqs = resdict2['seqs'][e_ind]

        true_rt.extend(np.array([rts[iorig] for iorig in resdict2['iorig']])[e_ind])
        # true_rt.extend(resdict2['rt'][e_ind])
        true_rt = np.array(true_rt)
        true_isotopes.extend(np.array([Isotopes[iorig] for iorig in resdict2['iorig']])[e_ind])
        # true_isotopes.extend(resdict2['Isotopes'][e_ind])
        true_isotopes = np.array(true_isotopes)

        e_all = abs(resdict2['md'][e_ind] - mass_shift) / (mass_sigma)
        zs_all_tmp = e_all ** 2

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

        if args['ts'] != 2 and deeplc_path:

            
            outtrain_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
            outtrain = open(outtrain_name, 'w')
            outcalib_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
            outcalib = open(outcalib_name, 'w')
            outres_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
            ns = true_seqs
            nr = true_rt
            print('Peptides used for RT prediction: %d' % (len(ns), ))
            ns2 = true_seqs2
            nr2 = true_rt2

            outtrain.write('seq,modifications,tr\n')
            for seq, RT in zip(ns2, nr2):
                mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
                outtrain.write(seq + ',' + str(mods_tmp) + ',' + str(RT) + '\n')
            outtrain.close()

            outcalib.write('seq,modifications,tr\n')
            for seq, RT in zip(ns, nr):
                mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
                outcalib.write(seq + ',' + str(mods_tmp) + ',' + str(RT) + '\n')
            outcalib.close()


            subprocess.call([deeplc_path, '--file_pred', outcalib_name, '--file_cal', outtrain_name, '--file_pred_out', outres_name])
            pepdict = dict()
            train_RT = []
            train_seq = []
            for x in open(outres_name).readlines()[1:]:
                _, seq, _, RTexp, RT = x.strip().split(',')
                pepdict[seq] = float(RT)
                train_seq.append(seq)
                train_RT.append(float(RTexp))


            train_RT = np.array(train_RT)
            RT_pred = np.array([pepdict[s] for s in train_seq])

            rt_diff_tmp = RT_pred - train_RT
            RT_left = -min(rt_diff_tmp)
            RT_right = max(rt_diff_tmp)

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

            aa, bb, RR, ss = aux.linear_regression(RT_pred, train_RT)

        else:

            if args['ts'] != 2 and elude_path:


                outtrain_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
                outtrain = open(outtrain_name, 'w')
                outcalib_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
                outcalib = open(outcalib_name, 'w')
                outres_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())

                ns = true_seqs
                nr = true_rt
                print('Peptides used for RT prediction: %d' % (len(ns), ))
                ns2 = true_seqs2
                nr2 = true_rt2
                for seq, RT in zip(ns, nr):
                    outtrain.write(seq + '\t' + str(RT) + '\n')
                outtrain.close()
                for seq, RT in zip(ns, nr):
                    outcalib.write(seq + '\t' + str(RT) + '\n')
                outcalib.close()

                subprocess.call([elude_path, '-t', outtrain_name, '-e', outcalib_name, '-a', '-g', '-o', outres_name])
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



        best_sigma = XRT_sigma
        RT_sigma = XRT_sigma

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



            outtrain_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
            outtrain = open(outtrain_name, 'w')
            outres_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())

            print('Peptides used for RT prediction: %d' % (len(ns), ))
            ll = len(ns)
            ns = ns[:ll]
            nr = nr[:ll]

            outtrain.write('seq,modifications,tr\n')
            for seq, RT in zip(ns, nr):
                mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
                outtrain.write(seq + ',' + str(mods_tmp) + ',' + str(RT) + '\n')
            outtrain.close()

            subprocess.call([deeplc_path, '--file_pred', outtrain_name, '--file_cal', outtrain_name, '--file_pred_out', outres_name])
            pepdict = dict()
            train_RT = []
            train_seq = []
            for x in open(outres_name).readlines()[1:]:
                _, seq, _, RTexp, RT = x.strip().split(',')
                pepdict[seq] = float(RT)
                train_seq.append(seq)
                train_RT.append(float(RTexp))


            train_RT = np.array(train_RT)
            RT_pred = np.array([pepdict[s] for s in train_seq])

            rt_diff_tmp = RT_pred - train_RT
            RT_left = -min(rt_diff_tmp)
            RT_right = max(rt_diff_tmp)

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

            aa, bb, RR, ss = aux.linear_regression(RT_pred, train_RT)

        else:

            if elude_path:


                outtrain_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
                outtrain = open(outtrain_name, 'w')
                outres_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())

                print(len(ns))
                ll = len(ns)
                ns = ns[:ll]
                nr = nr[:ll]
                for seq, RT in zip(ns, nr):
                    outtrain.write(seq + '\t' + str(RT) + '\n')
                outtrain.close()

                subprocess.call([elude_path, '-t', outtrain_name, '-e', outtrain_name, '-a', '-g', '-o', outres_name])
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



        best_sigma = XRT_sigma
        RT_sigma = XRT_sigma


    out_log.write('Calibrated RT shift: %s\n' % (XRT_shift, ))
    out_log.write('Calibrated RT sigma: %s\n' % (XRT_sigma, ))
    out_log.close()

    p1 = set(resdict['seqs'])

    n = args['nproc']

    if deeplc_path:


        pepdict = dict()

        outtrain_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
        outtrain = open(outtrain_name, 'w')
        outres_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())

        outtrain.write('seq,modifications,tr\n')
        for seq, RT in zip(ns, nr):
            mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
            outtrain.write(seq + ',' + str(mods_tmp) + ',' + str(RT) + '\n')
        outtrain.close()


        outtest_name = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
        outtest = open(outtest_name, 'w')


        outtest.write('seq,modifications\n')
        for seq in p1:
            mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
            outtest.write(seq + ',' + str(mods_tmp) + '\n')
        outtest.close()

        subprocess.call([deeplc_path, '--file_pred', outtest_name, '--file_cal', outtrain_name, '--file_pred_out', outres_name])
        for x in open(outres_name).readlines()[1:]:
            _, seq, _, RT = x.strip().split(',')
            pepdict[seq] = float(RT)

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
    rt_diff = np.array([rts[iorig] for iorig in resdict['iorig']]) - rt_pred
    # rt_diff = resdict['rt'] - rt_pred
    e_all = (rt_diff) ** 2 / (RT_sigma ** 2)
    zs_all = zs_all + e_all
    r = 9.0
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
        # for seq, md, rtd, peak_id, I, nScans, nIsotopes, mzr, rtr, av, ch, im in zip(resdict['seqs'], resdict['md'], rt_diff, resdict['ids'], resdict['Is'], resdict['Scans'], resdict['Isotopes'], resdict['mzraw'], resdict['rt'], resdict['av'], resdict['ch'], resdict['im']):
        for seq, md, rtd, iorig in zip(resdict['seqs'], resdict['md'], rt_diff, resdict['iorig']):
            peak_id = ids[iorig]
            I = Is[iorig]
            nScans = Scans[iorig]
            nIsotopes = Isotopes[iorig]
            mzr = mzraw[iorig]
            rtr = rts[iorig]
            av = avraw[iorig]
            ch = charges[iorig]
            im = imraw[iorig]
            output.write('\t'.join((seq, str(md), str(rtd), str(peak_id), str(I), str(nScans), str(nIsotopes), ';'.join(pept_prot[seq]), str(mzr), str(rtr), str(av), str(ch), str(im))) + '\n')
            
    mass_diff = (resdict['md'] - mass_shift) / (mass_sigma)

    rt_diff = (np.array([rts[iorig] for iorig in resdict['iorig']]) - rt_pred) / RT_sigma
    # rt_diff = (resdict['rt'] - rt_pred) / RT_sigma

    prefix = 'DECOY_'
    isdecoy = lambda x: x[0].startswith(prefix)
    isdecoy_key = lambda x: x.startswith(prefix)
    escore = lambda x: -x[1]

    SEED = 42

    # Hyperparameter grid
    param_grid = {
        'boosting_type': ['gbdt', ],
        'num_leaves': list(range(10, 1000)),
        'learning_rate': list(np.logspace(np.log10(0.001), np.log10(0.05), base = 10, num = 1000)),
        'metric': ['rmse', ],
        'verbose': [-1, ],
        'num_threads': [5, ],
        'n_estimators': [100, ],
    }

    def get_X_array(df, feature_columns):
        return df.loc[:, feature_columns].values

    def get_Y_array_pfms(df):
        return df.loc[:, 'decoy'].values

    def get_features_pfms(dataframe):
        feature_columns = dataframe.columns
        columns_to_remove = []
        banned_features = {
            'iorig',
            'ids',
            'seqs',
            'decoy',
            'preds',
            'av',
            'Scans',
            'proteins',
            'peptide',
        }
        for feature in feature_columns:
            if feature in banned_features or feature.startswith('c_'):
                columns_to_remove.append(feature)
        feature_columns = feature_columns.drop(columns_to_remove)
        return feature_columns

    def objective_pfms(df, hyperparameters, iteration):
        """Objective function for grid and random search. Returns
        the cross validation score from a set of hyperparameters."""
        
        all_res = []

        groups = df['peptide']
        ix = df.index.values
        unique = np.unique(groups)
        np.random.RandomState(SEED).shuffle(unique)
        result = []
        for split in np.array_split(unique, 3):
            mask = groups.isin(split)
            train, test = ix[~mask], ix[mask]
            train_df = df.iloc[train]
            test_df = df.iloc[test]

            feature_columns = get_features_pfms(df)
            model = get_cat_model_final_pfms(train_df, hyperparameters, feature_columns)
            test_df['preds'] = model.predict(get_X_array(test_df, feature_columns))
            all_res.append(test_df)
        test_df6 = pd.concat(all_res)
        shr_v = len(aux.filter(test_df6, fdr=0.25, key='preds', is_decoy='decoy'))
        
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

    def get_cat_model_pfms(df, hyperparameters, feature_columns, train, test):
        feature_columns = list(feature_columns)
        dtrain = lgb.Dataset(get_X_array(train, feature_columns), get_Y_array_pfms(train), feature_name=feature_columns, free_raw_data=False)
        dvalid = lgb.Dataset(get_X_array(test, feature_columns), get_Y_array_pfms(test), feature_name=feature_columns, free_raw_data=False)
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

    df1 = pd.DataFrame()
    for k in resdict.keys():
        df1[k] = resdict[k]

    df1['ids'] = df1['iorig'].apply(lambda x: ids[x])
    df1['Is'] = df1['iorig'].apply(lambda x: Is[x])
    # df1['Scans'] = df1['iorig'].apply(lambda x: Scans[x])
    df1['Isotopes'] = df1['iorig'].apply(lambda x: Isotopes[x])
    df1['mzraw'] = df1['iorig'].apply(lambda x: mzraw[x])
    df1['rt'] = df1['iorig'].apply(lambda x: rts[x])
    # df1['av'] = df1['iorig'].apply(lambda x: avraw[x])
    df1['ch'] = df1['iorig'].apply(lambda x: charges[x])
    df1['im'] = df1['iorig'].apply(lambda x: imraw[x])

    df1['mass_diff'] = mass_diff
    df1['rt_diff'] = rt_diff
    df1['decoy'] = df1['seqs'].apply(lambda x: all(z.startswith(prefix) for z in pept_prot[x]))

    df1['peptide'] = df1['seqs']
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

    if args['ml']:

        print('Start Machine Learning on PFMs...')

        print('Features used for MachineLearning: ', get_features_pfms(df1))

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
        ix = df1.index.values
        unique = np.unique(groups)
        np.random.RandomState(SEED).shuffle(unique)
        result = []
        for split in np.array_split(unique, 3):
            mask = groups.isin(split)
            train, test = ix[~mask], ix[mask]
            train_df = df1.iloc[train]
            test_df = df1.iloc[test]

            feature_columns = list(get_features_pfms(train_df))
            model = get_cat_model_final_pfms(train_df, bestparams, feature_columns)
            
            df1.loc[test, 'preds'] = model.predict(get_X_array(test_df, feature_columns))

    else:
        df1['preds'] = np.power(df1['mass_diff'], 2) + np.power(df1['rt_diff'], 2)

    df1['qpreds'] = pd.qcut(df1['preds'], 10, labels=range(10)) 
    df1['proteins'] = df1['seqs'].apply(lambda x: ';'.join(pept_prot[x]))

    df1.to_csv(base_out_name + '_PFMs_ML.csv', sep='\t', index=False)

    resdict['qpreds'] = df1['qpreds'].values
    resdict['ids'] = df1['ids'].values
    mass_diff = resdict['qpreds']
    rt_diff = resdict['qpreds']

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

    final_iteration(resdict, mass_diff, rt_diff, pept_prot, protsN, base_out_name, prefix, isdecoy, isdecoy_key, escore, fdr, args['nproc'], fname)


def worker(qin, qout, mass_diff, rt_diff, resdict, protsN, pept_prot, isdecoy_key, isdecoy, fdr, prots_spc_basic2, win_sys=False):

    for item in (iter(qin.get, None) if not win_sys else qin):
        mass_koef, rtt_koef = item
        e_ind = mass_diff <= mass_koef
        resdict2 = filter_results(resdict, e_ind)

        features_dict = dict()
        for pep in set(resdict2['seqs']):
            for bprot in pept_prot[pep]:
                prot_score = prots_spc_basic2[bprot]
                if prot_score > features_dict.get(pep, [-1, ])[-1]:
                    features_dict[pep] = (bprot, prot_score)

        prots_spc_basic = dict()

        p1 = set(resdict2['seqs'])

        pep_pid = defaultdict(set)
        pid_pep = defaultdict(set)
        banned_dict = dict()
        for pep, pid in zip(resdict2['seqs'], resdict2['ids']):
        # for pep, pid in zip(resdict2['seqs'], [ids[iorig] for iorig in resdict2['iorig']]):
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
