import os
import numpy as np
from scipy.stats import scoreatpercentile, rankdata
from scipy.optimize import curve_fit
import operator
from copy import copy, deepcopy
from collections import defaultdict
from pyteomics import parser, mass, auxiliary as aux, achrom
try:
    from pyteomics import cmass
except ImportError:
    cmass = mass
import subprocess
import tempfile
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from multiprocessing import Queue, Process, cpu_count
try:
    import seaborn
    seaborn.set(rc={'axes.facecolor':'#ffffff'})
    seaborn.set_style('whitegrid')
except:
    pass

from . import utils
from .utils_figures import plot_outfigures
import lightgbm as lgb
SEED = 50
import warnings
warnings.formatwarning = lambda msg, *args, **kw: str(msg) + '\n'

import logging
import numpy
import pandas
from sklearn import metrics
import csv
import ast


logger = logging.getLogger(__name__)

def worker_RT(qin, qout, shift, step, RC=False, ns=False, nr=False, win_sys=False):
    pepdict = dict()
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




def calc_protein_scores(p1, pept_prot, protsN, isdecoy_key, prefix, best_base_results=False, p=False):

    if best_base_results is not False:

        sorted_proteins_names = sorted(best_base_results.items(), key=lambda x: -x[1])
        sorted_proteins_names_dict = {}
        for idx, k in enumerate(sorted_proteins_names):
            sorted_proteins_names_dict[k[0]] = idx
    else:
        sorted_proteins_names_dict = False


    prots_spc2 = defaultdict(set)
    for pep, proteins in pept_prot.items():
        if pep in p1:
            if sorted_proteins_names_dict is not False:
                best_pos = min(sorted_proteins_names_dict[k] for k in proteins)
                for protein in proteins:
                    if sorted_proteins_names_dict[protein] == best_pos:
                        prots_spc2[protein].add(pep)
            else:
                for protein in proteins:
                    prots_spc2[protein].add(pep)


    for k in protsN:
        if k not in prots_spc2:
            prots_spc2[k] = set([])
    prots_spc = dict((k, len(v)) for k, v in prots_spc2.items())

    names_arr = np.array(list(prots_spc.keys()))
    v_arr = np.array(list(prots_spc.values()))
    n_arr = np.array([protsN[k] for k in prots_spc])

    if p is False:
        top100decoy_score = [prots_spc.get(dprot, 0) for dprot in protsN if isdecoy_key(dprot)]
        top100decoy_N = [val for key, val in protsN.items() if isdecoy_key(key)]
        p = np.mean(top100decoy_score) / np.mean(top100decoy_N)
        logger.info('Probability of random match for theoretical peptide = %.3f', (np.mean(top100decoy_score) / np.mean(top100decoy_N)))


    prots_spc = dict()
    all_pvals = utils.calc_sf_all(v_arr, n_arr, p)
    for idx, k in enumerate(names_arr):
        prots_spc[k] = all_pvals[idx]

    if best_base_results is not False:
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

    return prots_spc, p




def calibrate_RT_gaus_full(rt_diff_tmp):
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
    return XRT_shift, XRT_sigma, covvalue

def final_iteration(resdict, mass_diff, rt_diff, pept_prot, protsN, base_out_name, prefix, isdecoy, isdecoy_key, escore, fdr, nproc, out_log=False, fname=False, prots_spc_basic2=False, output_all=True, separate_figures=False):
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

        cnt_target_final = 0
        cnt_decoy_final = 0

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
                # decoy_set = list(decoy_set)


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


                # top100decoy_score_tmp = [tmp_spc_new.get(dprot, 0) for dprot in decoy_set]
                # top100decoy_score_tmp_sum = float(sum(top100decoy_score_tmp))
                top100decoy_score_tmp = {dprot: tmp_spc_new.get(dprot, 0) for dprot in decoy_set}
                top100decoy_score_tmp_sum = float(sum(top100decoy_score_tmp.values()))

            prots_spc = tmp_spc_new
            if not prots_spc_copy:
                prots_spc_copy = deepcopy(prots_spc)

            for v in decoy_set.intersection(unstable_prots):
                top100decoy_score_tmp_sum -= top100decoy_score_tmp[v]
                top100decoy_score_tmp[v] = prots_spc.get(v, 0)
                top100decoy_score_tmp_sum += top100decoy_score_tmp[v]
            p = top100decoy_score_tmp_sum / top100decoy_N

            n_change = set(protsN[k] for k in unstable_prots)
            if len(best_match_dict) == 0:
                n_change = sorted(n_change)
            for n_val in n_change:
                for k in n_map_dict[n_val]:
                    v = prots_spc[k]
                    if n_val not in best_match_dict or v > prots_spc[best_match_dict[n_val]]:
                        best_match_dict[n_val] = k
            n_arr_small = []
            names_arr_small = []
            v_arr_small = []
            # for k, v in best_match_dict.items():
            #     n_arr_small.append(k)
            #     names_arr_small.append(v)
            #     v_arr_small.append(prots_spc[v])



            max_k = 0

            for k, v in best_match_dict.items():
                num_matched = prots_spc[v]
                if num_matched >= max_k:

                    max_k = num_matched

                    n_arr_small.append(k)
                    names_arr_small.append(v)
                    v_arr_small.append(prots_spc[v])


            prots_spc_basic = dict()
            all_pvals = utils.calc_sf_all(np.array(v_arr_small), n_arr_small, p)
            for idx, k in enumerate(names_arr_small):
                prots_spc_basic[k] = all_pvals[idx]

            if not p0:
                p0 = float(p)

                prots_spc_tmp = dict()
                v_arr = np.array([prots_spc[k] for k in names_arr])
                all_pvals = utils.calc_sf_all(v_arr, n_arr, p)
                for idx, k in enumerate(names_arr):
                    prots_spc_tmp[k] = all_pvals[idx]

                sortedlist_spc = sorted(prots_spc_tmp.items(), key=operator.itemgetter(1))[::-1]
                if output_all:
                    with open(base_out_name + '_proteins_full_noexclusion.tsv', 'w') as output:
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


                if best_prot in decoy_set:
                    cnt_decoy_final += 1
                else:
                    cnt_target_final += 1


            else:

                v_arr = np.array([prots_spc[k] for k in names_arr])
                all_pvals = utils.calc_sf_all(v_arr, n_arr, p)
                for idx, k in enumerate(names_arr):
                    prots_spc_basic[k] = all_pvals[idx]

                for k, v in prots_spc_basic.items():
                    if k not in prots_spc_final:
                        prots_spc_final[k] = v

                break
            try:
                prot_fdr = cnt_decoy_final / cnt_target_final
                # prot_fdr = aux.fdr(prots_spc_final.items(), is_decoy=isdecoy)
            except:
                prot_fdr = 100.0
            if prot_fdr >= 12.5 * fdr:

                v_arr = np.array([prots_spc[k] for k in names_arr])
                all_pvals = utils.calc_sf_all(v_arr, n_arr, p)
                for idx, k in enumerate(names_arr):
                    prots_spc_basic[k] = all_pvals[idx]

                for k, v in prots_spc_basic.items():
                    if k not in prots_spc_final:
                        prots_spc_final[k] = v
                break

    if prots_spc_basic2 is False:
        prots_spc_basic2 = copy(prots_spc_final)
    else:
        prots_spc_basic2 = prots_spc_basic2
        for k in prots_spc_final:
            if k not in prots_spc_basic2:
                prots_spc_basic2[k] = 0
    prots_spc_final = dict()

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

    prots_spc = deepcopy(prots_spc_final)
    sortedlist_spc = sorted(prots_spc.items(), key=operator.itemgetter(1))[::-1]
    if output_all:
        with open(base_out_name + '_proteins_full.tsv', 'w') as output:
            output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\tdecoy\n')
            for x in sortedlist_spc:
                output.write('\t'.join((x[0], str(x[1]), str(prots_spc_copy[x[0]]), str(protsN[x[0]]), str(isdecoy(x)))) + '\n')

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

    logger.info('TOP 5 identified proteins:')
    logger.info('dbname\tscore\tmatched peptides\ttheoretical peptides')
    for x in filtered_prots[:5]:
        logger.info('\t'.join((str(x[0]), str(x[1]), str(int(prots_spc_copy[x[0]])), str(protsN[x[0]]))))
    logger.info('Final stage search: identified proteins = %d', identified_proteins)

    if output_all:
        with open(base_out_name + '_proteins.tsv', 'w') as output:
            output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\tdecoy\n')
            for x in filtered_prots:
                output.write('\t'.join((x[0], str(x[1]), str(prots_spc_copy[x[0]]), str(protsN[x[0]]), str(isdecoy(x)))) + '\n')

        if out_log or (fname and identified_proteins > 10):
            df1_peptides = pd.read_table(base_out_name + '_PFMs.tsv')
            df1_peptides['decoy'] = df1_peptides['proteins'].apply(lambda x: any(isdecoy_key(z) for z in x.split(';')))

            df1_proteins = pd.read_table(base_out_name + '_proteins_full.tsv')
            df1_proteins_f = pd.read_table(base_out_name + '_proteins.tsv')
            top_proteins = set(df1_proteins_f['dbname'])
            df1_peptides_f = df1_peptides[df1_peptides['proteins'].apply(lambda x: any(z in top_proteins for z in x.split(';')))]

            if fname and identified_proteins > 10:

                df0 = pd.read_table(os.path.splitext(fname)[0] + '.tsv')

                plot_outfigures(df0, df1_peptides, df1_peptides_f,
                    base_out_name, df_proteins=df1_proteins,
                    df_proteins_f=df1_proteins_f, prefix=prefix, separate_figures=separate_figures)

            if out_log is not False:
                df1_peptides_f = df1_peptides_f.sort_values(by='Intensity')
                df1_peptides_f = df1_peptides_f.drop_duplicates(subset='sequence', keep='last')
                dynamic_range_estimation = np.log10(df1_peptides_f['Intensity'].quantile(0.99)) - np.log10(df1_peptides_f['Intensity'].quantile(0.01))
                out_log.write('Estimated dynamic range in Log10 scale: %.1f\n' % (dynamic_range_estimation, ))
                out_log.write('Matched peptides for top-scored proteins: %d\n' % (len(df1_peptides_f), ))
                out_log.write('Identified proteins: %d\n' % (identified_proteins, ))
                out_log.close()

            return top_proteins

    else:
        top_proteins_fullinfo = []
        for x in filtered_prots:
            top_proteins_fullinfo.append((x[0], str(x[1]), str(prots_spc_copy[x[0]]), str(protsN[x[0]]), str(isdecoy(x))))


        return top_proteins_fullinfo

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
    for filename in args['files']:

        # Temporary for pyteomics <= Version 4.5.5 bug
        from pyteomics import mass
        if 'H-' in mass.std_aa_mass:
            del mass.std_aa_mass['H-']
        if '-OH' in mass.std_aa_mass:
            del mass.std_aa_mass['-OH']

        try:
            args['file'] = filename
            process_peptides(deepcopy(args))
        except Exception as e:
            logger.error(e)
            logger.error('Search is failed for file: %s', filename)
    return 1



def prepare_peptide_processor(fname, args):
    global nmasses
    global rts
    global charges
    global ids
    global Is
    global Isums
    global Scans
    global Isotopes
    global mzraw
    global avraw
    global imraw

    min_ch = args['cmin']
    max_ch = args['cmax']

    min_isotopes = args['i']
    min_scans = args['sc']

    logger.info('Reading file %s', fname)

    df_features = utils.iterate_spectra(fname, min_ch, max_ch, min_isotopes, min_scans, args['nproc'], args['check_unique'])

    # Sort by neutral mass
    df_features = df_features.sort_values(by='massCalib')

    nmasses = df_features['massCalib'].values
    rts = df_features['rtApex'].values
    charges = df_features['charge'].values
    ids = df_features['id'].values
    Is = df_features['intensityApex'].values
    if 'intensitySum' in df_features.columns:
        Isums = df_features['intensitySum'].values
    else:
        Isums = df_features['intensityApex'].values
        logger.info('intensitySum column is missing in peptide features. Using intensityApex instead')

    Scans = df_features['nScans'].values
    Isotopes = df_features['nIsotopes'].values
    mzraw = df_features['mz'].values
    avraw = np.zeros(len(df_features))
    if len(set(df_features['FAIMS'])) > 1:
        imraw = df_features['FAIMS'].values
    else:
        imraw = df_features['im'].values

    logger.info('Number of peptide isotopic clusters passed filters: %d', len(nmasses))

    aa_mass, aa_to_psi = utils.get_aa_mass_with_fixed_mods(args['fmods'], args['fmods_legend'])

    acc_l = args['ptol']
    acc_r = args['ptol']

    return {'aa_mass': aa_mass, 'aa_to_psi': aa_to_psi, 'acc_l': acc_l, 'acc_r': acc_r, 'args': args}, df_features

def get_resdict(it, **kwargs):

    resdict = {
        'seqs': [],
        'md': [],
        'mods': [],
        'iorig': [],
    }

    for seqm in it:
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
        end = nmasses.searchsorted(m + dm_r, side='right')
        for i in range(start, end):
            massdiff = (m - nmasses[i]) / m * 1e6
            mods = 0

            resdict['seqs'].append(seqm)
            resdict['md'].append(massdiff)
            resdict['mods'].append(mods)
            resdict['iorig'].append(i)


    for k in list(resdict.keys()):
        resdict[k] = np.array(resdict[k])

    return resdict



def filter_results(resultdict, idx):
    tmp = dict()
    for label in resultdict:
        tmp[label] = resultdict[label][idx]
    return tmp

def process_peptides(args):
    logger.info('Starting search...')

    fname_orig = args['file']
    if fname_orig.lower().endswith('mzml'):
        fname = os.path.splitext(fname_orig)[0] + '.features.tsv'
    else:
        fname = fname_orig

    fdr = args['fdr'] / 100
    min_isotopes_calibration = args['ci']
    min_scans_calibration = args['csc']
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
    out_log.write('Filename: %s\n' % (fname, ))
    out_log.write('Used parameters: %s\n' % (args.items(), ))

    deeplc_path = args['deeplc']
    if deeplc_path:
        from deeplc import DeepLC
        logging.getLogger('deeplc').setLevel(logging.ERROR)

    deeplc_model_path = args['deeplc_model_path']
    deeplc_model_path = deeplc_model_path.strip()

    if len(deeplc_model_path) > 0:
        if deeplc_model_path.endswith('.hdf5'):
            path_model = deeplc_model_path
        else:
            path_model = [os.path.join(deeplc_model_path,f) for f in os.listdir(deeplc_model_path) if f.endswith(".hdf5")]

    else:
        path_model = None

    if args['deeplc_library']:

        path_to_lib = args['deeplc_library']
        if not os.path.isfile(path_to_lib):
            lib_file = open(path_to_lib, 'w')
            lib_file.close()
        write_library = True

    else:
        path_to_lib = None
        write_library = False



    calib_path = args['pl']
    calib_path = calib_path.strip()

    if calib_path and args['ts']:
        args['ts'] = 0
        logger.info('Two-stage RT prediction does not work with list of MS/MS identified peptides...')

    args['enzyme'] = utils.get_enzyme(args['e'])


    prefix = args['prefix']
    protsN, pept_prot, ml_correction = utils.get_prot_pept_map(args)

    kwargs, df_features = prepare_peptide_processor(fname_orig, args)
    logger.info('Running the search ...')
    out_log.write('Number of features: %d\n' % (len(df_features), ))

    resdict = get_resdict(pept_prot, **kwargs)

    aa_to_psi = kwargs['aa_to_psi']

    if args['mc'] > 0:
        resdict['mc'] = np.array([parser.num_sites(z, args['enzyme']) for z in resdict['seqs']])

    isdecoy = lambda x: x[0].startswith(prefix)
    isdecoy_key = lambda x: x.startswith(prefix)
    escore = lambda x: -x[1]



    e_ind = np.array([Isotopes[iorig] for iorig in resdict['iorig']]) >= min_isotopes_calibration
    resdict2 = filter_results(resdict, e_ind)

    e_ind = np.array([Scans[iorig] for iorig in resdict2['iorig']]) >= min_scans_calibration
    resdict2 = filter_results(resdict2, e_ind)

    e_ind = resdict2['mods'] == 0
    resdict2 = filter_results(resdict2, e_ind)

    if args['mc'] > 0:
        e_ind = resdict2['mc'] == 0
        resdict2 = filter_results(resdict2, e_ind)

    p1 = set(resdict2['seqs'])

    if len(p1):

        # Calculate basic protein scores including homologues
        prots_spc, p = calc_protein_scores(p1, pept_prot, protsN, isdecoy_key, prefix, best_base_results=False, p=False)
        # Calculate basic protein scores excluding homologues
        prots_spc, p = calc_protein_scores(p1, pept_prot, protsN, isdecoy_key, prefix, best_base_results=prots_spc, p=p)




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
        # logger.info('Stage 0 search: probability of random match for theoretical peptide = %.3f', (np.mean(top100decoy_score) / np.mean(top100decoy_N)))

        # prots_spc = dict()
        # all_pvals = utils.calc_sf_all(v_arr, n_arr, p)
        # for idx, k in enumerate(names_arr):
        #     prots_spc[k] = all_pvals[idx]

        # # checked = set()
        # # for k, v in list(prots_spc.items()):
        # #     if k not in checked:
        # #         if isdecoy_key(k):
        # #             if prots_spc.get(k.replace(prefix, ''), -1e6) > v:
        # #                 del prots_spc[k]
        # #                 checked.add(k.replace(prefix, ''))
        # #         else:
        # #             if prots_spc.get(prefix + k, -1e6) > v:
        # #                 del prots_spc[k]
        # #                 checked.add(prefix + k)





        # for k in protsN:
        #     if k not in prots_spc:
        #         prots_spc[k] = 0

        # sorted_proteins_names = sorted(prots_spc.items(), key=lambda x: -x[1])
        # sorted_proteins_names_dict = {}
        # for idx, k in enumerate(sorted_proteins_names):
        #     sorted_proteins_names_dict[k[0]] = idx

        # prots_spc2 = defaultdict(set)
        # for pep, proteins in pept_prot.items():
        #     if pep in p1:
        #         best_pos = min(sorted_proteins_names_dict[k] for k in proteins)
        #         for protein in proteins:
        #             if sorted_proteins_names_dict[protein] == best_pos:
        #                 prots_spc2[protein].add(pep)

        # for k in protsN:
        #     if k not in prots_spc2:
        #         prots_spc2[k] = set([])
        # prots_spc = dict((k, len(v)) for k, v in prots_spc2.items())

        # names_arr = np.array(list(prots_spc.keys()))
        # v_arr = np.array(list(prots_spc.values()))
        # n_arr = np.array([protsN[k] for k in prots_spc])

        # prots_spc = dict()
        # all_pvals = utils.calc_sf_all(v_arr, n_arr, p)
        # for idx, k in enumerate(names_arr):
        #     prots_spc[k] = all_pvals[idx]

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



















        filtered_prots = aux.filter(prots_spc.items(), fdr=0.05, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1,
                                    full_output=True)

        identified_proteins = 0

        for x in filtered_prots:
            identified_proteins += 1
        logger.info('Stage 0 search: identified proteins = %d', identified_proteins)
        if identified_proteins <= 25:
            logger.info('Low number of identified proteins, using first 25 top scored proteins for calibration...')
            filtered_prots = sorted(prots_spc.items(), key=lambda x: -x[1])[:25]

        logger.info('Running mass recalibration...')

        df1 = pd.DataFrame()
        df1['mass diff'] = resdict['md']
        df1['mc'] = (resdict['mc'] if args['mc'] > 0 else 0)
        df1['iorig'] = resdict['iorig']
        df1['seqs'] = resdict['seqs']

        df1['mods'] = resdict['mods']

        df1['nIsotopes'] = [Isotopes[iorig] for iorig in df1['iorig'].values]
        df1['nScans'] = [Scans[iorig] for iorig in df1['iorig'].values]

        # df1['orig_md'] = true_md


        true_seqs = set()
        true_prots = set(x[0] for x in filtered_prots)
        for pep, proteins in pept_prot.items():
            if any(protein in true_prots for protein in proteins):
                true_seqs.add(pep)

        df1['top_peps'] = (df1['mc'] == 0) & (df1['seqs'].apply(lambda x: x in true_seqs) & (df1['nIsotopes'] >= min_isotopes_calibration) & (df1['nScans'] >= min_scans_calibration))

        mass_calib_arg = args['mcalib']

        assert mass_calib_arg in [0, 1, 2]

        if mass_calib_arg:
            df1['RT'] = rts[df1['iorig'].values]#df1['iorig'].apply(lambda x: rts[x])

            if mass_calib_arg == 2:
                df1['im'] = imraw[df1['iorig'].values]#df1['iorig'].apply(lambda x: imraw[x])
            elif mass_calib_arg == 1:
                df1['im'] = 0

            im_set = set(df1['im'].unique())
            if len(im_set) <= 5:
                df1['im_qcut'] = df1['im']
                for im_value in im_set:
                    idx1 = df1['im'] == im_value
                    df1.loc[idx1, 'qpreds'] = str(im_value) + pd.qcut(df1.loc[idx1, 'RT'], 5, labels=range(5)).astype(str)
            else:
                df1['im_qcut'] = pd.qcut(df1['im'], 5, labels=range(5)).astype(str)
                im_set = set(df1['im_qcut'].unique())
                for im_value in set(df1['im_qcut'].unique()):
                    idx1 = df1['im_qcut'] == im_value
                    df1.loc[idx1, 'qpreds'] = str(im_value) + pd.qcut(df1.loc[idx1, 'RT'], 5, labels=range(5)).astype(str)

            # df1['qpreds'] = pd.qcut(df1['RT'], 10, labels=range(10))#.astype(int)

            cor_dict = df1[df1['top_peps']].groupby('qpreds')['mass diff'].median().to_dict()

            rt_q_list = list(range(5))
            for im_value in im_set:
                for rt_q in rt_q_list:
                    lbl_cur = str(im_value) + str(rt_q)
                    if lbl_cur not in cor_dict:

                        best_diff = 1e6
                        best_val = 0
                        for rt_q2 in rt_q_list:
                            cur_diff = abs(rt_q - rt_q2)
                            if cur_diff != 0:
                                lbl_cur2 = str(im_value) + str(rt_q2)
                                if lbl_cur2 in cor_dict:
                                    if cur_diff < best_diff:
                                        best_diff = cur_diff
                                        best_val = cor_dict[lbl_cur2]

                        cor_dict[lbl_cur] = best_val

            df1['mass diff q median'] = df1['qpreds'].apply(lambda x: cor_dict[x])
            df1['mass diff corrected'] = df1['mass diff'] - df1['mass diff q median']

        else:
            df1['qpreds'] = 0
            df1['mass diff q median'] = 0
            df1['mass diff corrected'] = df1['mass diff']




        mass_left = args['ptol']
        mass_right = args['ptol']

        try:
            mass_shift_cor, mass_sigma_cor, covvalue_cor = calibrate_mass(0.001, mass_left, mass_right, df1[df1['top_peps']]['mass diff corrected'])
        except:
            mass_shift_cor, mass_sigma_cor, covvalue_cor = calibrate_mass(0.01, mass_left, mass_right, df1[df1['top_peps']]['mass diff corrected'])

        if mass_calib_arg:

            try:
                mass_shift, mass_sigma, covvalue = calibrate_mass(0.001, mass_left, mass_right, df1[df1['top_peps']]['mass diff'])
            except:
                mass_shift, mass_sigma, covvalue = calibrate_mass(0.01, mass_left, mass_right, df1[df1['top_peps']]['mass diff'])

            logger.info('Uncalibrated mass shift: %.3f ppm', mass_shift)
            logger.info('Uncalibrated mass sigma: %.3f ppm', mass_sigma)

        logger.info('Estimated mass shift: %.3f ppm', mass_shift_cor)
        logger.info('Estimated mass sigma: %.3f ppm', mass_sigma_cor)

        out_log.write('Estimated mass shift: %.3f ppm\n' % (mass_shift_cor, ))
        out_log.write('Estimated mass sigma: %.3f ppm\n' % (mass_sigma_cor, ))

        resdict['md'] = df1['mass diff corrected'].values

        mass_shift = mass_shift_cor
        mass_sigma = mass_sigma_cor

        e_all = abs(resdict['md'] - mass_shift) / (mass_sigma)
        r = 3.0
        e_ind = e_all <= r
        resdict = filter_results(resdict, e_ind)


        e_ind = np.array([Isotopes[iorig] for iorig in resdict['iorig']]) >= min_isotopes_calibration
        resdict2 = filter_results(resdict, e_ind)


        e_ind = np.array([Scans[iorig] for iorig in resdict2['iorig']]) >= min_scans_calibration
        resdict2 = filter_results(resdict2, e_ind)

        e_ind = resdict2['mods'] == 0
        resdict2 = filter_results(resdict2, e_ind)


        if args['mc'] > 0:
            e_ind = resdict2['mc'] == 0
            resdict2 = filter_results(resdict2, e_ind)

        p1 = set(resdict2['seqs'])


        # Calculate basic protein scores including homologues
        prots_spc, p = calc_protein_scores(p1, pept_prot, protsN, isdecoy_key, prefix, best_base_results=False, p=False)
        # Calculate basic protein scores excluding homologues
        prots_spc, p = calc_protein_scores(p1, pept_prot, protsN, isdecoy_key, prefix, best_base_results=prots_spc, p=p)





        filtered_prots = aux.filter(prots_spc.items(), fdr=0.05, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1,
                                    full_output=True)

        identified_proteins = 0

        for x in filtered_prots:
            identified_proteins += 1
        logger.info('Stage 1 search: identified proteins = %d', identified_proteins)
        if identified_proteins <= 25:
            logger.info('Low number of identified proteins, using first 25 top scored proteins for calibration...')
            filtered_prots = sorted(prots_spc.items(), key=lambda x: -x[1])[:25]



        logger.info('Running RT prediction...')


        e_ind = np.array([Isotopes[iorig] for iorig in resdict['iorig']]) >= 1
        resdict2 = filter_results(resdict, e_ind)

        e_ind = resdict2['mods'] == 0
        resdict2 = filter_results(resdict2, e_ind)

        if args['mc'] > 0:
            e_ind = resdict2['mc'] == 0
            resdict2 = filter_results(resdict2, e_ind)


        true_seqs = []
        true_rt = []
        true_isotopes = []
        true_prots = set(x[0] for x in filtered_prots)
        for pep, proteins in pept_prot.items():
            if any(protein in true_prots for protein in proteins):
                true_seqs.append(pep)
        e_ind = np.in1d(resdict2['seqs'], true_seqs)


        true_seqs = resdict2['seqs'][e_ind]

        true_rt.extend(np.array([rts[iorig] for iorig in resdict2['iorig']])[e_ind])
        true_rt = np.array(true_rt)
        true_isotopes.extend(np.array([Isotopes[iorig] for iorig in resdict2['iorig']])[e_ind])
        true_isotopes = np.array(true_isotopes)

        e_all = abs(resdict2['md'][e_ind] - mass_shift) / (mass_sigma)
        zs_all_tmp = e_all ** 2

        zs_all_tmp += (true_isotopes.max() - true_isotopes) * 100

        e_ind = np.argsort(zs_all_tmp)
        true_seqs = true_seqs[e_ind]
        true_rt = true_rt[e_ind]

        true_seqs = true_seqs[:2500]
        true_rt = true_rt[:2500]

        per_ind = np.random.RandomState(seed=SEED).permutation(len(true_seqs))
        true_seqs = true_seqs[per_ind]
        true_rt = true_rt[per_ind]

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
            true_seqs = df1['peptide'].values
            true_rt = df1['RT exp'].values

            ll = len(true_seqs)
            true_seqs2 = true_seqs[int(ll/2):]
            true_rt2 = true_rt[int(ll/2):]
            true_seqs = true_seqs[:int(ll/2)]
            true_rt = true_rt[:int(ll/2)]

        else:

            ll = len(true_seqs)

            true_seqs2 = true_seqs[int(ll/2):]
            true_rt2 = true_rt[int(ll/2):]
            true_seqs = true_seqs[:int(ll/2)]
            true_rt = true_rt[:int(ll/2)]

        ns = true_seqs
        nr = true_rt
        ns2 = true_seqs2
        nr2 = true_rt2


        logger.info('First-stage peptides used for RT prediction: %d', len(true_seqs))

        # if args['ts'] != 2 and deeplc_path:

        #     dlc = DeepLC(verbose=False, batch_num=args['deeplc_batch_num'], path_model=path_model, write_library=write_library, use_library=path_to_lib, pygam_calibration=False)


        #     df_for_calib = pd.DataFrame({
        #         'seq': ns2,
        #         'modifications': [utils.mods_for_deepLC(seq, aa_to_psi) for seq in ns2],
        #         'tr': nr2,
        #     })

        #     df_for_check = pd.DataFrame({
        #         'seq': ns,
        #         'modifications': [utils.mods_for_deepLC(seq, aa_to_psi) for seq in ns],
        #         'tr': nr,
        #     })

        #     try:
        #         dlc.calibrate_preds(seq_df=df_for_calib, check_df=df_for_check)
        #     except:
        #         dlc.calibrate_preds(seq_df=df_for_calib)

        #     df_for_check['pr'] =  dlc.make_preds(seq_df=df_for_check)
        #     df_for_calib['pr'] =  dlc.make_preds(seq_df=df_for_calib)
        #     nr2_pred = df_for_calib['pr']

        #     rt_diff_tmp = df_for_check['pr'] - df_for_check['tr']

        #     XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus_full(rt_diff_tmp)

        # else:
        RC = achrom.get_RCs_vary_lcp(ns2, nr2, metric='mae')
        nr2_pred = np.array([achrom.calculate_RT(s, RC) for s in ns2])
        nr_pred = np.array([achrom.calculate_RT(s, RC) for s in ns])

        rt_diff_tmp = nr_pred - nr

        XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus_full(rt_diff_tmp)


        logger.info('First-stage calibrated RT shift: %.3f min', XRT_shift)
        logger.info('First-stage calibrated RT sigma: %.3f min', XRT_sigma)

        RT_sigma = XRT_sigma

    else:
        logger.info('No matches found')





    if args['ts']:


        if args['es']:
            logger.info('Use extra Stage 2 search...')

            qin = list(set(resdict['seqs']))
            qout = []
            pepdict = worker_RT(qin, qout, 0, 1, RC, False, False, True)

            rt_pred = np.array([pepdict[s] for s in resdict['seqs']])
            rt_diff = np.array([rts[iorig] for iorig in resdict['iorig']]) - rt_pred - XRT_shift
            e_all = (rt_diff) ** 2 / (RT_sigma ** 2)
            r = 9.0
            e_ind = e_all <= r
            resdict = filter_results(resdict, e_ind)




            e_ind = np.array([Isotopes[iorig] for iorig in resdict['iorig']]) >= min_isotopes_calibration
            resdict2 = filter_results(resdict, e_ind)


            e_ind = np.array([Scans[iorig] for iorig in resdict2['iorig']]) >= min_scans_calibration
            resdict2 = filter_results(resdict2, e_ind)

            e_ind = resdict2['mods'] == 0
            resdict2 = filter_results(resdict2, e_ind)


            if args['mc'] > 0:
                e_ind = resdict2['mc'] == 0
                resdict2 = filter_results(resdict2, e_ind)

            p1 = set(resdict2['seqs'])


            # Calculate basic protein scores including homologues
            prots_spc, p = calc_protein_scores(p1, pept_prot, protsN, isdecoy_key, prefix, best_base_results=False, p=False)
            # Calculate basic protein scores excluding homologues
            prots_spc, p = calc_protein_scores(p1, pept_prot, protsN, isdecoy_key, prefix, best_base_results=prots_spc, p=p)





            filtered_prots = aux.filter(prots_spc.items(), fdr=0.05, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1,
                                        full_output=True)

            identified_proteins = 0

            for x in filtered_prots:
                identified_proteins += 1
            logger.info('Stage 2 search: identified proteins = %d', identified_proteins)
            if identified_proteins <= 25:
                logger.info('Low number of identified proteins, using first 25 top scored proteins for calibration...')
                filtered_prots = sorted(prots_spc.items(), key=lambda x: -x[1])[:25]







            e_ind = np.array([Isotopes[iorig] for iorig in resdict['iorig']]) >= 1
            resdict2 = filter_results(resdict, e_ind)

            e_ind = resdict2['mods'] == 0
            resdict2 = filter_results(resdict2, e_ind)

            if args['mc'] > 0:
                e_ind = resdict2['mc'] == 0
                resdict2 = filter_results(resdict2, e_ind)


            true_seqs = []
            true_rt = []
            true_isotopes = []
            true_prots = set(x[0] for x in filtered_prots)
            for pep, proteins in pept_prot.items():
                if any(protein in true_prots for protein in proteins):
                    true_seqs.append(pep)
            e_ind = np.in1d(resdict2['seqs'], true_seqs)


            true_seqs = resdict2['seqs'][e_ind]

            true_rt.extend(np.array([rts[iorig] for iorig in resdict2['iorig']])[e_ind])
            true_rt = np.array(true_rt)
            true_isotopes.extend(np.array([Isotopes[iorig] for iorig in resdict2['iorig']])[e_ind])
            true_isotopes = np.array(true_isotopes)

            e_all = abs(resdict2['md'][e_ind] - mass_shift) / (mass_sigma)
            zs_all_tmp = e_all ** 2

            zs_all_tmp += (true_isotopes.max() - true_isotopes) * 100

            e_ind = np.argsort(zs_all_tmp)
            true_seqs = true_seqs[e_ind]
            true_rt = true_rt[e_ind]

            idx_limit = len(true_seqs)
            cnt_pep = len(set(true_seqs))
            while cnt_pep > 2500:
                idx_limit -= 1
                cnt_pep = len(set(true_seqs[:idx_limit]))
            true_seqs = true_seqs[:idx_limit]
            true_rt = true_rt[:idx_limit]

            per_ind = np.random.RandomState(seed=SEED).permutation(len(true_seqs))
            true_seqs = true_seqs[per_ind]
            true_rt = true_rt[per_ind]

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
                true_seqs = df1['peptide'].values
                true_rt = df1['RT exp'].values

                ll = len(true_seqs)
                true_seqs2 = true_seqs[int(ll/2):]
                true_rt2 = true_rt[int(ll/2):]
                true_seqs = true_seqs[:int(ll/2)]
                true_rt = true_rt[:int(ll/2)]

            else:

                ll = len(true_seqs)

                true_seqs2 = true_seqs[int(ll/2):]
                true_rt2 = true_rt[int(ll/2):]
                true_seqs = true_seqs[:int(ll/2)]
                true_rt = true_rt[:int(ll/2)]

            ns = true_seqs
            nr = true_rt
            ns2 = true_seqs2
            nr2 = true_rt2




        else:
            ns = np.array(ns)
            nr = np.array(nr)
            idx = np.abs((rt_diff_tmp) - XRT_shift) <= 3 * XRT_sigma
            ns = ns[idx]
            nr = nr[idx]

            rt_diff_tmp2 = nr2_pred - nr2
            ns2 = np.array(ns2)
            nr2 = np.array(nr2)
            idx = np.abs((rt_diff_tmp2) - XRT_shift) <= 3 * XRT_sigma
            ns2 = ns2[idx]
            nr2 = nr2[idx]

        logger.info('Second-stage peptides used for RT prediction: %d', len(ns))

        if deeplc_path:

            dlc = DeepLC(verbose=False, batch_num=args['deeplc_batch_num'], path_model=path_model, write_library=write_library, use_library=path_to_lib, pygam_calibration=False)


            df_for_calib = pd.DataFrame({
                'seq': ns2,
                'modifications': [utils.mods_for_deepLC(seq, aa_to_psi) for seq in ns2],
                'tr': nr2,
            })

            df_for_check = pd.DataFrame({
                'seq': ns,
                'modifications': [utils.mods_for_deepLC(seq, aa_to_psi) for seq in ns],
                'tr': nr,
            })


            try:
                dlc.calibrate_preds(seq_df=df_for_calib, check_df=df_for_check)
            except:
                dlc.calibrate_preds(seq_df=df_for_calib)

            df_for_check['pr'] =  dlc.make_preds(seq_df=df_for_check)

            rt_diff_tmp = df_for_check['pr'] - df_for_check['tr']

            XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus_full(rt_diff_tmp)

        else:

            RC = achrom.get_RCs_vary_lcp(ns, nr, metric='mae')
            RT_pred = np.array([achrom.calculate_RT(s, RC) for s in ns])

            rt_diff_tmp = RT_pred - nr

            XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus_full(rt_diff_tmp)

        RT_sigma = XRT_sigma

    logger.info('Second-stage calibrated RT shift: %.3f min', XRT_shift)
    logger.info('Second-stage calibrated RT sigma: %.3f min', XRT_sigma)

    out_log.write('Calibrated RT shift: %.3f min\n' % (XRT_shift, ))
    out_log.write('Calibrated RT sigma: %.3f min\n' % (XRT_sigma, ))

    p1 = set(resdict['seqs'])

    n = args['nproc']


    def divide_chunks(l, n):
        for i in range(0, len(l), n):
            yield l[i:i + n]

    if deeplc_path:

        pepdict = dict()

        if args['save_calib']:
            with open(base_out_name + '_calib.tsv', 'w') as output:
                output.write('peptide\tRT exp\n')
                for seq, RT in zip(ns, nr):
                    output.write('%s\t%s\n' % (seq, str(RT)))
                for seq, RT in zip(ns2, nr2):
                    output.write('%s\t%s\n' % (seq, str(RT)))

        seqs_batch = list(p1)

        df_for_check = pd.DataFrame({
            'seq': seqs_batch,
            'modifications': [utils.mods_for_deepLC(seq, aa_to_psi) for seq in seqs_batch],
        })

        df_for_check['pr'] =  dlc.make_preds(seq_df=df_for_check)

        pepdict_batch = df_for_check.set_index('seq')['pr'].to_dict()

        pepdict.update(pepdict_batch)


    else:

        qin = list(p1)
        qout = []
        pepdict = worker_RT(qin, qout, 0, 1, RC, False, False, True)

    rt_pred = np.array([pepdict[s] for s in resdict['seqs']])
    # rt_diff = np.array([rts[iorig] for iorig in resdict['iorig']]) - rt_pred
    rt_diff = np.array([rts[iorig] for iorig in resdict['iorig']]) - rt_pred - XRT_shift
    # rt_diff = resdict['rt'] - rt_pred
    e_all = (rt_diff) ** 2 / (RT_sigma ** 2)
    r = 9.0
    e_ind = e_all <= r

    resdict = filter_results(resdict, e_ind)
    rt_diff = rt_diff[e_ind]
    rt_pred = rt_pred[e_ind]


    logger.info('RT prediction was finished')


    with open(base_out_name + '_protsN.tsv', 'w') as output:
        output.write('dbname\ttheor peptides\n')
        for k, v in protsN.items():
            output.write('\t'.join((k, str(v))) + '\n')

    with open(base_out_name + '_PFMs.tsv', 'w') as output:
        output.write('sequence\tmass diff\tRT diff\tpeak_id\tIntensity\tIntensitySum\tnScans\tnIsotopes\tproteins\tm/z\tRT\taveragineCorr\tcharge\tion_mobility\n')
        for seq, md, rtd, iorig in zip(resdict['seqs'], resdict['md'], rt_diff, resdict['iorig']):
            peak_id = ids[iorig]
            I = Is[iorig]
            Isum = Isums[iorig]
            nScans = Scans[iorig]
            nIsotopes = Isotopes[iorig]
            mzr = mzraw[iorig]
            rtr = rts[iorig]
            av = avraw[iorig]
            ch = charges[iorig]
            im = imraw[iorig]
            output.write('\t'.join((seq, str(md), str(rtd), str(peak_id), str(I), str(Isum), str(nScans), str(nIsotopes), ';'.join(pept_prot[seq]), str(mzr), str(rtr), str(av), str(ch), str(im))) + '\n')

    mass_diff = (resdict['md'] - mass_shift) / (mass_sigma)

    rt_diff = (np.array([rts[iorig] for iorig in resdict['iorig']]) - rt_pred - XRT_shift) / RT_sigma

    ## Will it override the settings provided in args? :k
    prefix = 'DECOY_'
    isdecoy = lambda x: x[0].startswith(prefix)
    isdecoy_key = lambda x: x.startswith(prefix)
    escore = lambda x: -x[1]

    param_grid = {
        'boosting_type': ['gbdt', ],
        'num_leaves': list(range(10, 1000)),
        'learning_rate': list(np.logspace(np.log10(0.001), np.log10(0.3), base = 10, num = 1000)),
        'min_child_samples': list(range(1, 1000, 5)),
        'reg_alpha': list(np.linspace(0, 1)),
        'reg_lambda': list(np.linspace(0, 1)),
        'colsample_bytree': list(np.linspace(0.01, 1, 100)),
        'subsample': list(np.linspace(0.01, 1, 100)),
        'is_unbalance': [True, False],
        'metric': ['rmse', ],
        'verbose': [-1, ],
        'num_threads': [args['nproc'], ],
        # 'use_quantized_grad': [True, ]
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
            'Is',
            # 'Scans',
            'proteins',
            'peptide',
            'md',
            'qpreds',
            'decoy1',
            'decoy2',
            'top_25_targets',
            'G',
        }

        for feature in feature_columns:
            if feature in banned_features:
                columns_to_remove.append(feature)
        feature_columns = feature_columns.drop(columns_to_remove)
        return feature_columns

    def objective_pfms(df, hyperparameters, iteration, threshold=0):
        """Objective function for grid and random search. Returns
        the cross validation score from a set of hyperparameters."""

        all_res = []

        for group_val in range(3):

            mask = df['G'] == group_val
            test_df = df[mask]
            test_ids = set(test_df['ids'])
            train_df = df[(~mask) & (df['ids'].apply(lambda x: x not in test_ids))]



            feature_columns = get_features_pfms(df)
            model = get_cat_model_final_pfms(train_df[~train_df['decoy2']], hyperparameters, feature_columns)

            df.loc[mask, 'preds'] = model.predict(get_X_array(df.loc[mask, :], feature_columns))

            test_df = df[mask]

            fpr, tpr, thresholds = metrics.roc_curve(get_Y_array_pfms(test_df[~test_df['decoy2']]), test_df[~test_df['decoy2']]['preds'])
            shr_v = metrics.auc(fpr, tpr)

            all_res.append(shr_v)

            if shr_v < threshold:
                all_res = [0, ]
                break

        shr_v = np.mean(all_res)

        return np.array([shr_v, hyperparameters, iteration, all_res], dtype=object)

    def random_search_pfms(df, param_grid, out_file, max_evals):
        """Random search for hyperparameter optimization.
        Writes result of search to csv file every search iteration."""

        threshold = 0



        # Dataframe for results
        results = pd.DataFrame(columns = ['sharpe', 'params', 'iteration', 'all_res'],
                                    index = list(range(max_evals)))
        for i in range(max_evals):

            # Choose random hyperparameters
            random_params = {k: np.random.RandomState(seed=SEED).choice(v, 1)[0] for k, v in param_grid.items()}

            # Evaluate randomly selected hyperparameters
            eval_results = objective_pfms(df, random_params, i, threshold)
            results.loc[i, :] = eval_results

            threshold = max(threshold, np.mean(eval_results[3]) - 3 * np.std(eval_results[3]))

            # open connection (append option) and write results
            of_connection = open(out_file, 'a')
            writer = csv.writer(of_connection)
            writer.writerow(eval_results)
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
        model = lgb.train(hyperparameters, dtrain, num_boost_round=500, valid_sets=(dvalid,), valid_names=('valid',), verbose_eval=False,
                    early_stopping_rounds=10, evals_result=evals_result)
        return model

    def get_cat_model_final_pfms(df, hyperparameters, feature_columns):
        feature_columns = list(feature_columns)
        train = df
        dtrain = lgb.Dataset(get_X_array(train, feature_columns), get_Y_array_pfms(train), feature_name=feature_columns, free_raw_data=False)
        np.random.seed(SEED)
        model = lgb.train(hyperparameters, dtrain, num_boost_round=100)

        return model


    logger.info('Prepare ML features')

    df1 = pd.DataFrame()
    for k in resdict.keys():
        df1[k] = resdict[k]

    df1['ids'] = ids[df1['iorig'].values]
    df1['Is'] = Is[df1['iorig'].values]
    df1['Scans'] = Scans[df1['iorig'].values]
    df1['Isotopes'] = Isotopes[df1['iorig'].values]
    df1['mzraw'] = mzraw[df1['iorig'].values]
    df1['rt'] = rts[df1['iorig'].values]
    df1['av'] = avraw[df1['iorig'].values]
    df1['ch'] = charges[df1['iorig'].values]
    df1['im'] = imraw[df1['iorig'].values]

    df1['mass_diff'] = mass_diff
    df1['rt_diff'] = rt_diff
    df1['decoy'] = df1['seqs'].apply(lambda x: all(z.startswith(prefix) for z in pept_prot[x]))

    df1['peptide'] = df1['seqs']
    # mass_dict = {}
    # pI_dict = {}
    # charge_dict = {}
    # for pep in set(df1['peptide']):
    #     try:
    #         mass_dict[pep] = mass.fast_mass2(pep)
    #         pI_dict[pep] = electrochem.pI(pep)
    #         charge_dict[pep] = electrochem.charge(pep, pH=7.0)
    #     except:
    #         mass_dict[pep] = 0
    #         pI_dict[pep] = 0
    #         charge_dict[pep] = 0

    df1['plen'] = df1['peptide'].apply(lambda z: len(z))


    # df1['mass'] = df1['peptide'].apply(lambda x: mass_dict[x])
    # df1['pI'] = df1['peptide'].apply(lambda x: pI_dict[x])
    # df1['charge_theor'] = df1['peptide'].apply(lambda x: charge_dict[x])

    for aa in mass.std_aa_mass:
        df1['c_%s' % (aa, )] = df1['peptide'].apply(lambda x: x.count(aa))
    df1['c_DP'] = df1['peptide'].apply(lambda x: x.count('DP'))
    df1['c_KP'] = df1['peptide'].apply(lambda x: x.count('KP'))
    df1['c_RP'] = df1['peptide'].apply(lambda x: x.count('RP'))

    df1['rt_diff_abs'] = df1['rt_diff'].abs()
    df1['rt_diff_abs_pdiff'] = df1['rt_diff_abs'] - df1.groupby('ids')['rt_diff_abs'].transform('median')
    df1['rt_diff_abs_pnorm'] = df1['rt_diff_abs'] / (df1.groupby('ids')['rt_diff_abs'].transform('sum') + 1e-2)

    df1['mass_diff_abs'] = df1['mass_diff'].abs()
    df1['mass_diff_abs_pdiff'] = df1['mass_diff_abs'] - df1.groupby('ids')['mass_diff_abs'].transform('median')
    df1['mass_diff_abs_pnorm'] = df1['mass_diff_abs'] / (df1.groupby('ids')['mass_diff_abs'].transform('sum') + 1e-2)

    df1['id_count'] = df1.groupby('ids')['mass_diff'].transform('count')

    #limited CSD information
    if args['csd'] == 1:
        logging.info('Using limited CSD')
        df1['seq_count'] = df1.groupby('peptide')['mass_diff'].transform('count')
        df1['charge_count'] = df1.groupby('peptide')['ch'].transform('nunique')
        df1['im_count'] = df1.groupby('peptide')['im'].transform('nunique')

    #complete CSD information
    elif args['csd'] == 2:
        logging.info('Using complete CSD')
        #imports and functions: should be global?
        import json
        from keras import models

        os.environ["CUDA_VISIBLE_DEVICES"] = ""

        #constants
        PROTON = 1.00727646677
        TINY = 1e-6

        psi_to_single_ptm = {'[Acetyl]-': 'B',
                             '[Carbamidomethyl]': '',
                             'M[Oxidation]': 'J',
                             'S[Phosphorylation]': 'X',
                             'T[Phosphorylation]': 'O',
                             'Y[Phosphorylation]': 'U'}
        #functions
        def reshapeOneHot(X):
            X = np.dstack(X)
            X = np.swapaxes(X, 1, 2)
            X = np.swapaxes(X, 0, 1)
            return X

        def get_single_ptm_code(psi_sequence):
            sequence = psi_sequence
            for ptm in psi_to_single_ptm:
                sequence = sequence.replace(ptm, psi_to_single_ptm[ptm])
            return sequence

        def one_hot_encode_peptide(psi_sequence, MAX_LENGTH = 40):
            peptide = get_single_ptm_code(psi_sequence)
            if len(peptide) > MAX_LENGTH:
                print('Peptide length is larger than maximal length of ', str(MAX_LENGTH))
                return ['', None]
            else:
                AA_vocabulary = 'KRPTNAQVSGILCMJHFYWEDBXOU'#B: acetyl; J: oxidized Met; O: PhosphoT; X: PhosphoS; U: PhosphoY

                one_hot_peptide = np.zeros((len(peptide), len(AA_vocabulary)))

                for j in range(0, len(peptide)):
                    try:
                        aa = peptide[j]
                        one_hot_peptide[j, AA_vocabulary.index(aa)] = 1
                    except:
                        print(f'WARN: "{aa}" is not in vocabulary; it will be skipped')

                no_front_paddings = int((MAX_LENGTH - len(peptide))/2)
                peptide_front_paddings = np.zeros((no_front_paddings, one_hot_peptide.shape[1]))

                no_back_paddings = MAX_LENGTH - len(peptide) - no_front_paddings
                peptide_back_paddings = np.zeros((no_back_paddings, one_hot_peptide.shape[1]))

                full_one_hot_peptide = np.vstack((peptide_front_paddings, one_hot_peptide, peptide_back_paddings))

                return peptide, full_one_hot_peptide

        def expand_charges(row, min_charge, max_charge, mode=1):
            charge = np.arange(min_charge, max_charge + 1).reshape(-1, 1)
            mz = row['massCalib']/charge + PROTON * mode
            index = [f'{int(row["id"])}_{z[0]}' for z in charge]
            result = pd.DataFrame(mz, columns=['mz'], index=index)
            result['z'] = charge
            result['rt_start'] = row['rtStart']
            result['rt_end'] = row['rtEnd']
            if 'FAIMS' in row.keys():
                result['scan_filter'] = f'ms cv={row["FAIMS"]}'

            return result

        def AUC(x, y):
            return (0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1])).sum()

        def aggregate_csd(g):
            result = {k: v for k,v in zip(g['z'], g['area'])}
            result['pepid'] = g['pepid'].iloc[0]
            return result

        #imports and functions: end

        logger.info('Running charge-state distribution prediction')

        ## extracting feature CSDs

        faims_mode = np.any(df_features['FAIMS'] != 0)
        logger.debug(f'FAIMS mode: {faims_mode}')

        #find all features that has been maped to sequence
        if faims_mode:
            used_features = df_features.loc[df1['ids'].unique(), ['massCalib', 'rtStart', 'rtEnd', 'FAIMS']].copy()
        else:
            used_features = df_features.loc[df1['ids'].unique(), ['massCalib', 'rtStart', 'rtEnd']].copy()
        used_features['id'] = used_features.index.astype(int)

        logger.debug(f'{used_features.shape[0]} features has been mapped to peptide sequence')

        #expand each feature to 1 - 4 charge states
        mz_table = pd.concat(used_features.apply(expand_charges, args=(1, 4), axis=1).tolist())
        #need better way to filter impossible m/z-s
        mzmin = np.floor(mzraw.min())
        mzmax = np.ceil(mzraw.max())
        mz_table = mz_table.query('mz >= @mzmin & mz <= @mzmax').copy()

        logger.debug(f'{mz_table.shape[0]} XICs scheduled for extraction')

        #writing JSON input for TRFP
        mz_table['tolerance'] = 5
        mz_table['tolerance_unit'] = 'ppm'
        mz_table['comment'] = mz_table.index
        trfp_xic_in = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
        with open(trfp_xic_in, 'w') as out:
            json.dump(mz_table.drop(['z'], axis=1).apply(lambda row: {k:v for k,v in zip(row.index, row)}, axis=1).tolist(), out)

        #running TRFP to extract XICs (needs RAW file)
        trfp_xic_out = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())

        trfp_params = [args['trfp'], 'xic', '-i', f'{base_out_name[:-9]}.raw', '-b', trfp_xic_out, '-j', trfp_xic_in]

        if os.name != 'nt' and os.path.splitext(args['trfp'])[1].lower() == '.exe':
            trfp_params.insert(0, 'mono')

        logger.info('Running XIC extraction')

        try:
            subprocess.check_output(trfp_params)
            os.remove(trfp_xic_in)
        except subprocess.CalledProcessError as ex:
            logger.error(f'TRFP execution failed\nOutput:\n{ex.output}')
            raise Exception('TRFP XIC extraction failed')

        #read TRFP output
        with open(trfp_xic_out, 'r') as xic_input:
            chromatograms = json.load(xic_input)

        os.remove(trfp_xic_out)

        #calculate AUCs
        mz_table['area'] = [AUC(np.array(cr['RetentionTimes']), np.array(cr['Intensities'])) for cr in chromatograms['Content']]
        mz_table['pepid'] = mz_table['comment'].str.split('_', expand=True).iloc[:, 0].astype(int)

        #building CSD table
        csd_table = pd.DataFrame(mz_table.groupby('pepid').apply(aggregate_csd).tolist())
        csd_table.rename(columns={k:v for k,v in zip(csd_table.columns, csd_table.columns.astype(str))}, inplace=True)
        csd_table = csd_table.reindex(columns=sorted(csd_table.columns))

        #add CSD information
        used_features = used_features.join(csd_table.set_index('pepid', verify_integrity=True), how='left')

        #map CSD of features for the PFM information
        pfm_csd = used_features.reindex(df1['ids'])[['1', '2', '3', '4']].values

        ## predict CSDs for all sequences using the model

        #loading model
        CSD_model = models.load_model(
            os.path.normpath(os.path.join(os.path.dirname(__file__), 'models', 'CSD_model_LCMSMS.hdf5')),
            compile=False)
        CSD_model.compile()

        #using only unique sequences for prediction
        unique_sequences = pd.DataFrame(df1.drop_duplicates('seqs')['seqs'])
        CSD_X = reshapeOneHot(unique_sequences['seqs'].apply(lambda s: one_hot_encode_peptide(s)[1]).values)

        logger.debug(f'{unique_sequences.shape[0]} unique peptide sequence detected')
        logger.info('Running charge-state distribution prediction model')

        #prediction
        CSD_pred = CSD_model.predict(CSD_X, batch_size=2048)
        unique_sequences[['1', '2', '3', '4']] = CSD_pred[:, :4]
        pfm_csd_pred = unique_sequences.set_index('seqs').reindex(df1['seqs']).values

        #masking impossible m/z-s in prediction
        mask = np.isnan(pfm_csd)
        pfm_csd[mask] = 0
        pfm_csd_pred[mask] = 0

        #normalizing CSDs
        pfm_csd = pfm_csd / (pfm_csd.sum(axis=1) + TINY).reshape(-1, 1)
        pfm_csd_pred = pfm_csd_pred / (pfm_csd_pred.sum(axis=1) + TINY).reshape(-1, 1)

        #spectrum angle feature
        df1['z_angle'] = 1 - 2 * np.arccos((pfm_csd * pfm_csd_pred).sum(axis=1) /\
            (np.linalg.norm(pfm_csd, 2, axis=1) * np.linalg.norm(pfm_csd_pred, 2, axis=1) + TINY)) / np.pi

        #adding average_charge
        pfm_csd = np.concatenate([pfm_csd, (pfm_csd * np.arange(1, 5).reshape(1, 4)).sum(axis=1).reshape(-1, 1)], axis=1)
        pfm_csd_pred = np.concatenate([pfm_csd_pred, (pfm_csd_pred * np.arange(1, 5).reshape(1, 4)).sum(axis=1).reshape(-1, 1)], axis=1)

        #prediction error features for individual intensities and average charge
        z_deltas = pfm_csd - pfm_csd_pred
        df1[[f'z{z}_err' for z in '1234a']] = (z_deltas - z_deltas.mean(axis=0)) / z_deltas.std(axis=0).reshape(1, -1)

        logger.info('Charge-state distribution features added')

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
    logger.info('Stage 3 search: probability of random match for theoretical peptide = %.3f', p)

    prots_spc = dict()
    all_pvals = utils.calc_sf_all(v_arr, n_arr, p)
    for idx, k in enumerate(names_arr):
        prots_spc[k] = all_pvals[idx]

    target_prots_25_fdr = set([x[0] for x in aux.filter(prots_spc.items(), fdr=0.25, key=escore, is_decoy=isdecoy, remove_decoy=False, formula=1, full_output=True, correction=0)])
    df1['proteins'] = df1['seqs'].apply(lambda x: ';'.join(pept_prot[x]))
    df1['decoy2'] = df1['decoy']
    df1['decoy'] = df1['proteins'].apply(lambda x: all(z not in target_prots_25_fdr for z in x.split(';')))
    df1['top_25_targets'] = df1['decoy']


    if len(target_prots_25_fdr) <= 25:
        logger.info('Low number of identified proteins, turning off LightGBM...')
        filtered_prots = sorted(prots_spc.items(), key=lambda x: -x[1])[:25]
        skip_ml = 1
    else:
        skip_ml = 0

    if args['ml'] and not skip_ml:

        logger.info('Start Machine Learning on PFMs...')

        MAX_EVALS = 25

        out_file = os.path.join(tempfile.gettempdir(), os.urandom(24).hex())
        of_connection = open(out_file, 'w')
        writer = csv.writer(of_connection)

        headers = ['auc', 'params', 'iteration', 'all_res']
        writer.writerow(headers)
        of_connection.close()

        all_id_list = list(set(df1[df1['decoy']]['peptide']))
        np.random.RandomState(seed=SEED).shuffle(all_id_list)
        seq_gmap = {}
        for idx, split in enumerate(np.array_split(all_id_list, 3)):
            for id_ftr in split:
                seq_gmap[id_ftr] = idx

        all_id_list = list(set(df1[~df1['decoy']]['peptide']))
        np.random.RandomState(seed=SEED).shuffle(all_id_list)
        for idx, split in enumerate(np.array_split(all_id_list, 3)):
            for id_ftr in split:
                seq_gmap[id_ftr] = idx



        df1['G'] = df1['peptide'].apply(lambda x: seq_gmap[x])

        # train = df1[~df1['decoy']]
        # train_extra = df1[df1['decoy']]
        # train_extra = train_extra.sample(frac=min(1.0, len(train)/len(train_extra)))
        # train = pd.concat([train, train_extra], ignore_index=True).reset_index(drop=True)
        # random_results = random_search_pfms(train, param_grid, out_file, MAX_EVALS)

        random_results = random_search_pfms(df1, param_grid, out_file, MAX_EVALS)

        random_results = pd.read_csv(out_file)
        random_results = random_results[random_results['auc'] != 'auc']
        random_results['params'] = random_results['params'].apply(lambda x: ast.literal_eval(x))
        convert_dict = {'auc': float,
                    }
        random_results = random_results.astype(convert_dict)


        bestparams = random_results.sort_values(by='auc',ascending=False)['params'].values[0]

        bestparams['num_threads'] = args['nproc']



        for group_val in range(3):

            mask = df1['G'] == group_val
            test_df = df1[mask]
            test_ids = set(test_df['ids'])
            train_df = df1[(~mask) & (df1['ids'].apply(lambda x: x not in test_ids))]

            # mask2 = train['G'] != group_val

            # train_df = train[(mask2) & (train['ids'].apply(lambda x: x not in test_ids))]


            feature_columns = list(get_features_pfms(train_df))
            model = get_cat_model_final_pfms(train_df[~train_df['decoy2']], bestparams, feature_columns)

            df1.loc[mask, 'preds'] = rankdata(model.predict(get_X_array(test_df, feature_columns)), method='ordinal') / sum(mask)


    else:
        df1['preds'] = np.power(df1['mass_diff'], 2) + np.power(df1['rt_diff'], 2)




    df1['qpreds'] = pd.qcut(df1['preds'], 50, labels=range(50))

    df1['decoy'] = df1['decoy2']


    df1u = df1.sort_values(by='preds')
    df1u = df1u.drop_duplicates(subset='seqs')

    qval_ok = 0
    for qval_cur in range(50):
        df1ut = df1u[df1u['qpreds'] == qval_cur]
        decoy_ratio = df1ut['decoy'].sum() / len(df1ut)
        if decoy_ratio < ml_correction:
            qval_ok = qval_cur
        else:
            break
    logger.info('%d %% of PFMs were removed from protein scoring after Machine Learning', (100 - (qval_ok+1)*2))

#CHECK merge conflict was from here
    df1u = df1u[df1u['qpreds'] <= qval_ok]#.copy()

    df1u['qpreds'] = pd.qcut(df1u['preds'], 10, labels=range(10))

    qdict = df1u.set_index('seqs').to_dict()['qpreds']

    df1['qpreds'] = df1['seqs'].apply(lambda x: qdict.get(x, 11))

    # df1un = df1u[df1u['qpreds'] <= qval_ok].copy()
    # df1un['qpreds'] = pd.qcut(df1un['preds'], 10, labels=range(10))

    # qdict = df1un.set_index('seqs').to_dict()['qpreds']


    # df1['qpreds'] = df1['seqs'].apply(lambda x: qdict.get(x, 11))

    # # df1['qz'] = df1['seqs'].apply(lambda x: qval_dict.get(qdict.get(x, 11), 1.0))

#to here

    df1.to_csv(base_out_name + '_PFMs_ML.tsv', sep='\t', index=False)

    df1 = df1[df1['qpreds'] <= 10]

    resdict = {}
    resdict['seqs'] = df1['seqs'].values
    resdict['qpreds'] = df1['qpreds'].values
    resdict['ids'] = df1['ids'].values

    mass_diff = resdict['qpreds']
    rt_diff = []
    if skip_ml:
        mass_diff = np.zeros(len(mass_diff))

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
    logger.info('Final stage search: probability of random match for theoretical peptide = %.3f', p)


    prots_spc = dict()
    all_pvals = utils.calc_sf_all(v_arr, n_arr, p)
    for idx, k in enumerate(names_arr):
        prots_spc[k] = all_pvals[idx]

    sf = args['separate_figures']

    top_proteins = final_iteration(resdict, mass_diff, rt_diff, pept_prot, protsN, base_out_name, prefix, isdecoy, isdecoy_key, escore, fdr, args['nproc'], out_log, fname, separate_figures=sf)

    # pept_prot_limited = dict()
    # for k, v in pept_prot.items():
    #     tmp = set(zz for zz in v if zz in top_proteins)
    #     if len(tmp):
    #         pept_prot_limited[k] = tmp

    # resdict3 = get_resdict(pept_prot_limited, acc_l=25, acc_r=25, aa_mass=kwargs['aa_mass'])

    # p1 = set(resdict3['seqs'])

    # if deeplc_path:

    #     pepdict = dict()

    #     seqs_batch = list(p1)

    #     df_for_check = pd.DataFrame({
    #         'seq': seqs_batch,
    #         'modifications': [utils.mods_for_deepLC(seq, aa_to_psi) for seq in seqs_batch],
    #     })

    #     df_for_check['pr'] =  dlc.make_preds(seq_df=df_for_check)

    #     pepdict_batch = df_for_check.set_index('seq')['pr'].to_dict()

    #     pepdict.update(pepdict_batch)


    # else:

    #     qin = list(p1)
    #     qout = []
    #     pepdict = worker_RT(qin, qout, 0, 1, RC, False, False, True)

    # rt_pred = np.array([pepdict[s] for s in resdict3['seqs']])
    # rt_diff = np.array([rts[iorig] for iorig in resdict3['iorig']]) - rt_pred - XRT_shift
    # # e_all = (rt_diff) ** 2 / (RT_sigma ** 2)
    # # r = 9.0
    # # e_ind = e_all <= r
    # # resdict = filter_results(resdict, e_ind)
    # # rt_diff = rt_diff[e_ind]
    # # rt_pred = rt_pred[e_ind]



    # with open(base_out_name + '_PFMs_extended.tsv', 'w') as output:
    #     output.write('sequence\tmass diff\tRT diff\tpeak_id\tIntensity\tIntensitySum\tnScans\tnIsotopes\tproteins\tm/z\tRT\taveragineCorr\tcharge\tion_mobility\n')
    #     for seq, md, rtd, iorig in zip(resdict3['seqs'], resdict3['md'], rt_diff, resdict3['iorig']):
    #         peak_id = ids[iorig]
    #         I = Is[iorig]
    #         Isum = Isums[iorig]
    #         nScans = Scans[iorig]
    #         nIsotopes = Isotopes[iorig]
    #         mzr = mzraw[iorig]
    #         rtr = rts[iorig]
    #         av = avraw[iorig]
    #         ch = charges[iorig]
    #         im = imraw[iorig]
    #         output.write('\t'.join((seq, str(md), str(rtd), str(peak_id), str(I), str(Isum), str(nScans), str(nIsotopes), ';'.join(pept_prot[seq]), str(mzr), str(rtr), str(av), str(ch), str(im))) + '\n')


    # print('!!!', len(set(resdict3['seqs'])))

    logger.info('The search for file %s is finished.', base_out_name)

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
            pep_pid[pep].add(pid)
            pid_pep[pid].add(pep)
            if pep in banned_dict:
                banned_dict[pep] += 1
            else:
                banned_dict[pep] = 1

        banned_pids_total = set()

        if len(p1):
            prots_spc_final = dict()
            prots_spc_copy = False
            prots_spc2 = False
            unstable_prots = set()
            p0 = False

            prev_best_score = 1e6

            names_arr = False
            tmp_spc_new = False
            decoy_set = False

            cnt_target_final = 0
            cnt_decoy_final = 0

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
                    # decoy_set = list(decoy_set)


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


                    top100decoy_score_tmp = {dprot: tmp_spc_new.get(dprot, 0) for dprot in decoy_set}
                    top100decoy_score_tmp_sum = float(sum(top100decoy_score_tmp.values()))

                    prots_spc_basic = False

                prots_spc = tmp_spc_new
                if not prots_spc_copy:
                    prots_spc_copy = deepcopy(prots_spc)

                # for idx, v in enumerate(decoy_set):
                #     if v in unstable_prots:
                #         top100decoy_score_tmp_sum -= top100decoy_score_tmp[idx]
                #         top100decoy_score_tmp[idx] = prots_spc.get(v, 0)
                #         top100decoy_score_tmp_sum += top100decoy_score_tmp[idx]


                for v in decoy_set.intersection(unstable_prots):
                    top100decoy_score_tmp_sum -= top100decoy_score_tmp[v]
                    top100decoy_score_tmp[v] = prots_spc.get(v, 0)
                    top100decoy_score_tmp_sum += top100decoy_score_tmp[v]


                # p = float(sum(top100decoy_score_tmp)) / top100decoy_N
                p = top100decoy_score_tmp_sum / top100decoy_N
                if not p0:
                    p0 = float(p)

                n_change = set(protsN[k] for k in unstable_prots)
                if len(best_match_dict) == 0:
                    n_change = sorted(n_change)
                for n_val in n_change:
                    for k in n_map_dict[n_val]:
                        v = prots_spc[k]
                        if n_val not in best_match_dict or v > prots_spc[best_match_dict[n_val]]:
                            best_match_dict[n_val] = k
                n_arr_small = []
                names_arr_small = []
                v_arr_small = []
                max_k = 0

                for k, v in best_match_dict.items():
                    num_matched = prots_spc[v]
                    if num_matched >= max_k:

                        max_k = num_matched

                        n_arr_small.append(k)
                        names_arr_small.append(v)
                        v_arr_small.append(prots_spc[v])

                prots_spc_basic = dict()
                all_pvals = utils.calc_sf_all(np.array(v_arr_small), n_arr_small, p)

                for idx, k in enumerate(names_arr_small):
                    prots_spc_basic[k] = all_pvals[idx]

                best_prot = utils.keywithmaxval(prots_spc_basic)

                best_score = min(prots_spc_basic[best_prot], prev_best_score)
                prev_best_score = best_score

                unstable_prots = set()
                if best_prot not in prots_spc_final:
                    prots_spc_final[best_prot] = best_score
                    banned_pids = set()
                    for pep in prots_spc2[best_prot]:
                        for pid in pep_pid[pep]:
                            if pid not in banned_pids_total:
                                banned_pids.add(pid)
                    for pid in banned_pids:#banned_pids.difference(banned_pids_total):
                        for pep in pid_pep[pid]:
                            banned_dict[pep] -= 1
                            if banned_dict[pep] == 0:
                                best_prot_val = features_dict[pep][0]
                                tmp_spc_new[best_prot_val] -= 1
                                unstable_prots.add(best_prot_val)
                                # for bprot in pept_prot[pep]:
                                #     if bprot == best_prot_val:
                                #         tmp_spc_new[bprot] -= 1
                                #         unstable_prots.add(bprot)

                    banned_pids_total.update(banned_pids)

                    if best_prot in decoy_set:
                        cnt_decoy_final += 1
                    else:
                        cnt_target_final += 1

                else:

                    v_arr = np.array([prots_spc[k] for k in names_arr])
                    all_pvals = utils.calc_sf_all(v_arr, n_arr, p)
                    for idx, k in enumerate(names_arr):
                        prots_spc_basic[k] = all_pvals[idx]

                    for k, v in prots_spc_basic.items():
                        if k not in prots_spc_final:
                            prots_spc_final[k] = v

                    break

                prots_spc_basic[best_prot] = 0

                try:
                    prot_fdr = cnt_decoy_final / cnt_target_final
                    # prot_fdr = aux.fdr(prots_spc_final.items(), is_decoy=isdecoy)
                except:
                    prot_fdr = 100.0
                if prot_fdr >= 12.5 * fdr:

                    v_arr = np.array([prots_spc[k] for k in names_arr])
                    all_pvals = utils.calc_sf_all(v_arr, n_arr, p)
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
