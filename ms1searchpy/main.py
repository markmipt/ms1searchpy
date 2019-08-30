import os
import cPickle
import utils
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

def worker_RT(qin, qout, shift, step, RC=False, elude_path=False, ns=False, nr=False):



    pepdict = dict()    
    if elude_path:
        outtrain = tempfile.NamedTemporaryFile(suffix='.txt')
        outres = tempfile.NamedTemporaryFile(suffix='.txt')
        outres_name = outres.name
        outres.close()
        for seq, RT in zip(ns, nr):
            outtrain.write(seq + '\t' + str(RT) + '\n')
        outtrain.flush()

        outtest = tempfile.NamedTemporaryFile(suffix='.txt')

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
    qout.put(pepdict)
    qout.put(None)

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
    for i in xrange(start, end):
        peak_id = ids[i]
        I = Is[i]
        massdiff = (m - nmasses[i]) / m * 1e6
        results.append((seqm, massdiff, rts[i], peak_id, I, Scans[i], Isotopes[i], mzraw[i], avraw[i], charges[i]))
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
    nmasses = []
    rts = []
    charges = []
    ids = []
    Is = []
    Scans = []
    Isotopes = []
    mzraw = []
    avraw = []

    min_ch = args['cmin']
    max_ch = args['cmax']

    min_isotopes = args['i']
    min_scans = args['sc']

    print 'Reading spectra ...'
    for m, RT, c, peak_id, I, nScans, nIsotopes, mzr, avr in utils.iterate_spectra(fname, min_ch, max_ch, min_isotopes, min_scans):
        nmasses.append(m)
        rts.append(RT)
        charges.append(c)
        ids.append(peak_id)
        Is.append(I)
        Scans.append(nScans)
        Isotopes.append(nIsotopes)
        mzraw.append(mzr)
        avraw.append(avr)

    print(len(nmasses))

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
    try:
        outpath = args['outpath']
    except:
        outpath = False

    elude_path = args['elude']
    elude_path = elude_path.strip()

    args['enzyme'] = utils.get_enzyme(args['e'])

    ms1results = []
    peps = utils.peptide_gen(args)
    kwargs = prepare_peptide_processor(fname, args)
    func = peptide_processor_iter_isoforms
    print 'Running the search ...'
    for y in utils.multimap(1, func, peps, **kwargs):
        for result in y:
            if len(result):
                ms1results.extend(result)

    prefix = args['prefix']
    protsN, pept_prot = utils.get_prot_pept_map(args)

    resdict = get_results(ms1results)
    del ms1results

    isdecoy = lambda x: x[0].startswith(prefix)
    isdecoy_key = lambda x: x.startswith(prefix)
    escore = lambda x: -x[1]

    # e_ind = resdict['mc'] == 0
    # resdict2 = filter_results(resdict, e_ind)

    e_ind = resdict['Isotopes'] >= 4
    # e_ind = resdict['Isotopes'] >= 1
    resdict2 = filter_results(resdict, e_ind)

    p1 = set(resdict2['seqs'])

    if len(p1):
        prots_spc2 = defaultdict(set)
        for pep, proteins in pept_prot.iteritems():
            if pep in p1:
                for protein in proteins:
                    prots_spc2[protein].add(pep)

        for k in protsN:
            if k not in prots_spc2:
                prots_spc2[k] = set([])
        prots_spc = dict((k, len(v)) for k, v in prots_spc2.iteritems())

        names_arr = np.array(prots_spc.keys())
        v_arr = np.array(prots_spc.values())
        n_arr = np.array([protsN[k] for k in prots_spc])

        top100decoy_score = [prots_spc.get(dprot, 0) for dprot in protsN if isdecoy_key(dprot)]
        top100decoy_N = [val for key, val in protsN.items() if isdecoy_key(key)]
        p = np.mean(top100decoy_score) / np.mean(top100decoy_N)
        print 'p=%s' % (np.mean(top100decoy_score) / np.mean(top100decoy_N))

        prots_spc = dict()
        all_pvals = calc_sf_all(v_arr, n_arr, p)
        for idx, k in enumerate(names_arr):
            prots_spc[k] = all_pvals[idx]

        checked = set()
        for k, v in prots_spc.items():
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
        print 'results for default search: number of identified proteins = %d' % (identified_proteins, )

        print 'Running mass recalibration...'

        true_md = []
        true_isotopes = []
        true_seqs = []
        true_prots = set(x[0] for x in filtered_prots)
        for pep, proteins in pept_prot.iteritems():
            if any(protein in true_prots for protein in proteins):
                true_seqs.append(pep)

        e_ind = np.in1d(resdict2['seqs'], true_seqs)

        true_seqs = resdict2['seqs'][e_ind]
        true_md.extend(resdict2['md'][e_ind])
        true_md = np.array(true_md)
        true_isotopes.extend(resdict2['Isotopes'][e_ind])
        true_isotopes = np.array(true_isotopes)

        e_ind = true_isotopes >= 4
        # e_ind = true_isotopes >= 1
        true_md = true_md[e_ind]
        true_seqs = true_seqs[e_ind]

        mass_left = args['ptol']
        mass_right = args['ptol']

        mass_shift, mass_sigma, covvalue = calibrate_mass(0.001, mass_left, mass_right, true_md)
        if np.isinf(covvalue):
            mass_shift, mass_sigma, covvalue = calibrate_mass(0.01, mass_left, mass_right, true_md)
        print 'Calibrated mass shift: ', mass_shift
        print 'Calibrated mass sigma in ppm: ', mass_sigma
        # mass_sigma= mass_sigma * 1.1

        e_all = abs(resdict['md'] - mass_shift) / (mass_sigma)
        r = 3.0
        e_ind = e_all <= r
        resdict = filter_results(resdict, e_ind)

        zs_all = e_all[e_ind] ** 2


        # e_ind = resdict['mc'] == 0
        # resdict2 = filter_results(resdict, e_ind)

        p1 = set(resdict['seqs'])

        prots_spc2 = defaultdict(set)
        for pep, proteins in pept_prot.iteritems():
            if pep in p1:
                for protein in proteins:
                    prots_spc2[protein].add(pep)

        for k in protsN:
            if k not in prots_spc2:
                prots_spc2[k] = set([])
        prots_spc = dict((k, len(v)) for k, v in prots_spc2.iteritems())

        names_arr = np.array(prots_spc.keys())
        v_arr = np.array(prots_spc.values())
        n_arr = np.array([protsN[k] for k in prots_spc])

        top100decoy_score = [prots_spc.get(dprot, 0) for dprot in protsN if isdecoy_key(dprot)]
        top100decoy_N = [val for key, val in protsN.items() if isdecoy_key(key)]
        p = np.mean(top100decoy_score) / np.mean(top100decoy_N)
        print 'p=%s' % (np.mean(top100decoy_score) / np.mean(top100decoy_N))

        prots_spc = dict()
        all_pvals = calc_sf_all(v_arr, n_arr, p)
        for idx, k in enumerate(names_arr):
            prots_spc[k] = all_pvals[idx]

        checked = set()
        for k, v in prots_spc.items():
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
        print 'results for default search after mass calibration: number of identified proteins = %d' % (identified_proteins, )



        print 'Running RT prediction...'
        true_seqs = []
        true_rt = []
        true_isotopes = []
        true_prots = set(x[0] for x in filtered_prots[:100])#[:5])
        for pep, proteins in pept_prot.iteritems():
            if any(protein in true_prots for protein in proteins):
                true_seqs.append(pep)
        e_ind = np.in1d(resdict2['seqs'], true_seqs)


        true_seqs = resdict2['seqs'][e_ind]
        true_rt.extend(resdict2['rt'][e_ind])
        true_rt = np.array(true_rt)
        true_isotopes.extend(resdict2['Isotopes'][e_ind])
        true_isotopes = np.array(true_isotopes)

        e_all = abs(resdict2['md'][e_ind] - mass_shift) / (mass_sigma)
        zs_all_tmp = e_all ** 2

        # e_ind = true_isotopes >= 1
        e_ind = true_isotopes >= 4
        true_seqs = true_seqs[e_ind]
        true_rt = true_rt[e_ind]
        true_isotopes = true_isotopes[e_ind]
        zs_all_tmp = zs_all_tmp[e_ind]

        e_ind = np.argsort(zs_all_tmp)
        true_seqs = true_seqs[e_ind]
        true_rt = true_rt[e_ind]
        true_isotopes = true_isotopes[e_ind]

        true_seqs = true_seqs[:1000]
        true_rt = true_rt[:1000]
        true_isotopes = true_isotopes[:1000]

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

        RC, outmask = get_RCs2(true_seqs, true_rt)

        if elude_path:
            outtrain = tempfile.NamedTemporaryFile(suffix='.txt')
            outres = tempfile.NamedTemporaryFile(suffix='.txt')
            outres_name = outres.name
            outres.close()
            ns = true_seqs[~outmask]
            nr = true_rt[~outmask]
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

            start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 50
            print(start_width, 'SW')
            XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
            if np.isinf(covvalue):
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(0.1, RT_left, RT_right, rt_diff_tmp)
            if np.isinf(covvalue):
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(1.0, RT_left, RT_right, rt_diff_tmp)
            print 'Calibrated RT shift: ', XRT_shift
            print 'Calibrated RT sigma in ppm: ', XRT_sigma

            aa, bb, RR, ss = aux.linear_regression(RT_pred, train_RT)
        else:
            RC = achrom.get_RCs_vary_lcp(true_seqs[~outmask], true_rt[~outmask])
            RT_pred = np.array([achrom.calculate_RT(s, RC) for s in true_seqs])
            aa, bb, RR, ss = aux.linear_regression(RT_pred, true_rt)

            rt_diff_tmp = RT_pred - true_rt
            RT_left = -min(rt_diff_tmp)
            RT_right = max(rt_diff_tmp)

            start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 50
            print(start_width, 'SW')
            XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
            if np.isinf(covvalue):
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(0.1, RT_left, RT_right, rt_diff_tmp)
            if np.isinf(covvalue):
                XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(1.0, RT_left, RT_right, rt_diff_tmp)
            print 'Calibrated RT shift: ', XRT_shift
            print 'Calibrated RT sigma in ppm: ', XRT_sigma

        print aa, bb, RR, ss


        best_sigma = XRT_sigma#ss
        RT_sigma = XRT_sigma#best_sigma

    else:
        print 'No matches found'

    p1 = set(resdict['seqs'])

    n = args['nproc']

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
            for k, v in item.iteritems():
                pepdict[k] = v

    for p in procs:
        p.join()


    rt_pred = np.array([pepdict[s] for s in resdict['seqs']])
    rt_diff = resdict['rt'] - rt_pred
    e_all = (rt_diff) ** 2 / (RT_sigma ** 2)
    zs_all = zs_all + e_all
    r = 9.0#16.0#9.0
    e_ind = e_all <= r
    resdict = filter_results(resdict, e_ind)
    rt_diff = rt_diff[e_ind]
    zs_all = zs_all[e_ind]
    rt_pred = rt_pred[e_ind]


    if outpath:
        base_out_name = os.path.splitext(os.path.join(outpath, os.path.basename(fname)))[0]
    else:
        base_out_name = os.path.splitext(fname)[0]

    with open(base_out_name + '_protsN.csv', 'w') as output:
        output.write('dbname\ttheor peptides\n')
        for k, v in protsN.items():
            output.write('\t'.join((k, str(v))) + '\n')
   
    with open(base_out_name + '_PFMs.csv', 'w') as output:
        output.write('sequence\tmass diff\tRT diff\tpeak_id\tIntensity\tnScans\tnIsotopes\tproteins\tm/z\tRT\taveragineCorr\tcharge\n')
        for seq, md, rtd, peak_id, I, nScans, nIsotopes, mzr, rtr, av, ch in zip(resdict['seqs'], resdict['md'], rt_diff, resdict['ids'], resdict['Is'], resdict['Scans'], resdict['Isotopes'], resdict['mzraw'], resdict['rt'], resdict['av'], resdict['ch']):
            output.write('\t'.join((seq, str(md), str(rtd), str(peak_id), str(I), str(nScans), str(nIsotopes), ';'.join(pept_prot[seq]), str(mzr), str(rtr), str(av), str(ch))) + '\n')
            
    mass_diff = np.abs(resdict['md'] - mass_shift) / (mass_sigma)
    rt_diff = np.abs(resdict['rt'] - rt_pred) / RT_sigma

    def final_iteration(resdict, mass_diff, rt_diff, protsN, args):
        n = args['nproc']


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
                    for k, v in protsN.iteritems():
                        n_map_dict[v].append(k)

                    decoy_set = set()
                    for k in protsN:
                        if isdecoy_key(k):
                            decoy_set.add(k)
                    decoy_set = list(decoy_set)
                    

                    prots_spc2 = defaultdict(set)
                    for pep, proteins in pept_prot.iteritems():
                        if pep in p1:
                            for protein in proteins:
                                prots_spc2[protein].add(pep)

                    for k in protsN:
                        if k not in prots_spc2:
                            prots_spc2[k] = set([])
                    prots_spc2 = dict(prots_spc2)
                    unstable_prots = set(prots_spc2.keys())

                    top100decoy_N = sum([val for key, val in protsN.items() if isdecoy_key(key)])

                    names_arr = np.array(prots_spc2.keys())
                    n_arr = np.array([protsN[k] for k in names_arr])

                    tmp_spc_new = dict((k, len(v)) for k, v in prots_spc2.iteritems())


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
                for k, v in best_match_dict.iteritems():
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

                    for k, v in prots_spc_basic.iteritems():
                        if k not in prots_spc_final:
                            prots_spc_final[k] = v

                    break

                prot_fdr = aux.fdr(prots_spc_final.items(), is_decoy=isdecoy)
                if prot_fdr >= 12.5 * fdr:

                    v_arr = np.array([prots_spc[k] for k in names_arr])
                    all_pvals = calc_sf_all(v_arr, n_arr, p)
                    for idx, k in enumerate(names_arr):
                        prots_spc_basic[k] = all_pvals[idx]

                    for k, v in prots_spc_basic.iteritems():
                        if k not in prots_spc_final:
                            prots_spc_final[k] = v
                    break

        prots_spc_basic2 = copy(prots_spc_final)
        prots_spc_final = dict()

        if n == 0:
            try:
                n = cpu_count()
            except NotImplementedError:
                n = 1
        
        
        qin = Queue()
        qout = Queue()

        # for mc in [2, 3, 4, 5, ]:
        # for mc in [0, 1]:
            # for mass_koef in [3, 2, 1, 0.5, 0.25, 0.1]:
            #     for rtt_koef in [3, 2, 1, 0.5, 0.25, 0.1]:
            # for mass_koef in [1.0, 0.8, 0.6, 0.4, 0.2]:
            #     for rtt_koef in [1.0, 0.8, 0.6, 0.4, 0.2]:
        for mass_koef in np.arange(1.0, 0.2, -0.33):
            for rtt_koef in np.arange(1.0, 0.2, -0.33):
            # for mass_koef in np.arange(0.5, 0.0, -0.1):
            #     for rtt_koef in np.arange(0.5, 0.0, -0.1):
            # for mass_koef in np.arange(1.0, 0.0, -0.05):
                # for rtt_koef in np.arange(1.0, 0.0, -0.25):
                # for rtt_koef in [3, ]:
                qin.put((mass_koef, rtt_koef))

        # qin.put((1.0, 1.0))
        # qin.put(None)

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
        # worker(qin, qout, mass_diff, rt_diff, resdict, protsN, pept_prot, isdecoy_key, isdecoy, fdr, prots_spc_basic2)
        # prots_spc_final, prots_spc_copy = qout.get()

        for k in prots_spc_final.keys():
            prots_spc_final[k] = np.mean(prots_spc_final[k])

        prots_spc = prots_spc_final
        sortedlist_spc = sorted(prots_spc.iteritems(), key=operator.itemgetter(1))[::-1]
        with open(base_out_name + '_proteins_full.csv', 'w') as output:
            output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\n')
            for x in sortedlist_spc:
                # output.write('\t'.join((x[0], str(x[1]), str(protsN[x[0]]))) + '\n')
                output.write('\t'.join((x[0], str(x[1]), str(prots_spc_copy[x[0]]), str(protsN[x[0]]))) + '\n')

        checked = set()
        for k, v in prots_spc.items():
            if k not in checked:
                if isdecoy_key(k):
                    if prots_spc.get(k.replace(prefix, ''), -1e6) > v:
                        del prots_spc[k]
                        checked.add(k.replace(prefix, ''))
                else:
                    if prots_spc.get(prefix + k, -1e6) > v:
                        del prots_spc[k]
                        checked.add(prefix + k)

        filtered_prots = aux.filter(prots_spc.items(), fdr=fdr, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1, full_output=True)

        identified_proteins = 0

        for x in filtered_prots:
            identified_proteins += 1

        print 'TOP 5 identified proteins:'
        print 'dbname\tscore\tnum matched peptides\tnum theoretical peptides'
        for x in filtered_prots[:5]:
            # print '\t'.join((str(x[0]), str(x[1]), str(protsN[x[0]])))
            print '\t'.join((str(x[0]), str(x[1]), str(int(prots_spc_copy[x[0]])), str(protsN[x[0]])))
        print 'results:%s;number of identified proteins = %d' % (fname, identified_proteins, )
        print 'R=', r
        with open(base_out_name + '_proteins.csv', 'w') as output:
            output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\n')
            for x in filtered_prots:
                # output.write('\t'.join((x[0], str(x[1]), str(protsN[x[0]]))) + '\n')
                output.write('\t'.join((x[0], str(x[1]), str(prots_spc_copy[x[0]]), str(protsN[x[0]]))) + '\n')


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

    final_iteration(resdict, mass_diff, rt_diff, protsN, args)


def worker(qin, qout, mass_diff, rt_diff, resdict, protsN, pept_prot, isdecoy_key, isdecoy, fdr, prots_spc_basic2):

    for item in iter(qin.get, None):
        mass_koef, rtt_koef = item

        m_k = scoreatpercentile(mass_diff, mass_koef * 100)
        e_ind = mass_diff <= m_k
        resdict2 = filter_results(resdict, e_ind)
        rt_diff2 = rt_diff[e_ind]

        r_k = scoreatpercentile(rt_diff2, rtt_koef * 100)
        e_ind = rt_diff2 <= r_k
        resdict2 = filter_results(resdict2, e_ind)

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
                    for k, v in protsN.iteritems():
                        n_map_dict[v].append(k)

                    decoy_set = set()
                    for k in protsN:
                        if isdecoy_key(k):
                            decoy_set.add(k)
                    decoy_set = list(decoy_set)
                    

                    prots_spc2 = defaultdict(set)
                    for pep, proteins in pept_prot.iteritems():
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

                    names_arr = np.array(prots_spc2.keys())
                    n_arr = np.array([protsN[k] for k in names_arr])

                    tmp_spc_new = dict((k, len(v)) for k, v in prots_spc2.iteritems())


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
                for k, v in best_match_dict.iteritems():
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

                    for k, v in prots_spc_basic.iteritems():
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

                    for k, v in prots_spc_basic.iteritems():
                        if k not in prots_spc_final:
                            prots_spc_final[k] = v
                    break

        if mass_koef + rtt_koef >= 1.99:
            item2 = prots_spc_copy
        else:
            item2 = False

        qout.put((prots_spc_final, item2))
    qout.put(None)
