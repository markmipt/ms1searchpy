import os
import cPickle
import utils
import numpy as np
from scipy.stats import binom
from scipy.optimize import curve_fit
from scipy import exp
import operator
from copy import copy
from collections import defaultdict
import re
from pyteomics import parser, mass, fasta, auxiliary as aux, achrom
try:
    from pyteomics import cmass
except ImportError:
    cmass = mass
import subprocess
from sklearn import linear_model

def get_RCs2(sequences, RTs, lcp = -0.21,
            term_aa = False, **kwargs):
    labels = kwargs.get('labels')

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

    # Use least square linear regression.
    model = linear_model.LinearRegression(n_jobs=12)
    model.fit(np.array(composition_array), np.array(RTs))
    
    model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression(n_jobs=12), min_samples=0.3, max_trials=1000)

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


def process_file(fname, settings):
    ftype = fname.rsplit('.', 1)[-1].lower()
    utils.seen_target.clear()
    utils.seen_decoy.clear()
    return process_peptides(fname, settings)


def peptide_processor(peptide, **kwargs):
    seqm = peptide
    m = cmass.fast_mass(seqm, aa_mass=kwargs['aa_mass'])
    acc_l = kwargs['acc_l']
    acc_r = kwargs['acc_r']
    settings = kwargs['settings']
    dm_l = acc_l * m / 1.0e6
    dm_r = acc_r * m / 1.0e6
    start = nmasses.searchsorted(m - dm_l)
    end = nmasses.searchsorted(m + dm_r)
    idx = set(range(start, end))

    # if idx:
    #     RC = kwargs['RC']
    #     RT_sigma = kwargs['RT_sigma']
    #     RT = achrom.calculate_RT(seqm, RC)

    results = []
    for i in idx:
        # RTdiff = RT - rts[i]
        # if abs(RTdiff) <= 3 * RT_sigma:
        peak_id = ids[i]
        massdiff = (m - nmasses[i]) / m * 1e6
        results.append((seqm, massdiff, rts[i], peak_id))
    return results


def prepare_peptide_processor(fname, settings):
    global nmasses
    global rts
    global charges
    global ids
    nmasses = []
    rts = []
    charges = []
    ids = []

    min_ch = settings.getint('search', 'minimum charge')
    max_ch = settings.getint('search', 'maximum charge')
    min_isotopes = settings.getint('search', 'minimum isotopes')

    print 'Reading spectra ...'
    for m, RT, c, peak_id in utils.iterate_spectra(fname, min_ch, max_ch, min_isotopes):
        nmasses.append(m)
        rts.append(RT)
        charges.append(c)
        ids.append(peak_id)

    i = np.argsort(nmasses)
    nmasses = np.array(nmasses)[i]
    rts = np.array(rts)[i]
    charges = np.array(charges)[i]
    ids = np.array(ids)[i]

    fmods = settings.get('modifications', 'fixed')
    aa_mass = mass.std_aa_mass
    if fmods:
        for mod in re.split(r'[,;]\s*', fmods):
            m, aa = parser._split_label(mod)
            aa_mass[aa] += settings.getfloat('modifications', m)

    acc_l = settings.getfloat('search', 'precursor accuracy left')
    acc_r = settings.getfloat('search', 'precursor accuracy right')

    return {'aa_mass': aa_mass, 'acc_l': acc_l, 'acc_r': acc_r, 'settings': settings}


def peptide_processor_iter_isoforms(peptide, **kwargs):
    out = []
    out.append(peptide_processor(peptide, **kwargs))
    return out


def process_peptides(fname, settings):

    def calc_sf_all(v, n, p):
        sf_values = np.log10(1 / binom.sf(v, n, p))
        sf_values[np.isinf(sf_values)] = 1
        return sf_values

    elude_path = settings.get('performance', 'elude_path')
    elude_path = elude_path.strip()

    ms1results = []
    peps = utils.peptide_gen(settings)
    kwargs = prepare_peptide_processor(fname, settings)
    func = peptide_processor_iter_isoforms
    print 'Running the search ...'
    n = settings.getint('performance', 'processes')
    for y in utils.multimap(n, func, peps, **kwargs):
        for result in y:
            if len(result):
                ms1results.extend(result)

    prefix = settings.get('input', 'decoy prefix')
    protsN, pept_prot = utils.get_prot_pept_map(settings)

    seqs_all, md_all, rt_all, ids_all = zip(*ms1results)
    seqs_all = np.array(seqs_all)
    md_all = np.array(md_all)
    rt_all = np.array(rt_all)
    ids_all = np.array(ids_all)
    del ms1results

    isdecoy = lambda x: x[0].startswith(prefix)
    isdecoy_key = lambda x: x.startswith(prefix)
    escore = lambda x: -x[1]

    p1 = set(seqs_all)

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

        prots_spc_copy = copy(prots_spc)
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

        filtered_prots = aux.filter(prots_spc.items(), fdr=0.01, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1,
                                    full_output=True)

        identified_proteins = 0
        identified_proteins_valid = 0

        for x in filtered_prots:
            identified_proteins += 1

        # for x in filtered_prots[:5]:
        #     print x[0], x[1], int(prots_spc_copy[x[0]]), protsN[x[0]]
        print 'results for default search: number of identified proteins = %d' % (identified_proteins, )


        print 'Running mass recalibration...'


        true_md = []
        true_seqs = []
        true_prots = set(x[0] for x in filtered_prots)
        for pep, proteins in pept_prot.iteritems():
            if any(protein in true_prots for protein in proteins):
                true_seqs.append(pep)

        e_ind = np.in1d(seqs_all, true_seqs)

        true_md.extend(md_all[e_ind])

        mass_left = settings.getfloat('search', 'precursor accuracy left')
        mass_right = settings.getfloat('search', 'precursor accuracy right')
        bwidth = 0.01
        bbins = np.arange(-mass_left, mass_right, bwidth)
        H1, b1 = np.histogram(true_md, bins=bbins)
        b1 = b1 + bwidth
        b1 = b1[:-1]

        def noisygaus(x, a, x0, sigma, b):
            return a * exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b

        popt, pcov = curve_fit(noisygaus, b1, H1, p0=[1, np.median(true_md), 1, 1])
        mass_shift, mass_sigma = popt[1], abs(popt[2])
        print 'Calibrated mass shift: ', mass_shift
        print 'Calibrated mass sigma in ppm: ', mass_sigma
        # print np.sqrt(np.diag(pcov))
        # print '\n'

        e_all = abs(md_all - mass_shift) / (mass_sigma)
        r = 3.0
        # r = zz ** 2
        e_ind = e_all <= r
        seqs_all = seqs_all[e_ind]
        md_all = md_all[e_ind]
        rt_all = rt_all[e_ind]
        ids_all = ids_all[e_ind]


        print 'Running RT prediction...'
        true_seqs = []
        true_md = []
        true_rt = []
        true_prots = set(x[0] for x in filtered_prots[:5])
        for pep, proteins in pept_prot.iteritems():
            if any(protein in true_prots for protein in proteins):
                true_seqs.append(pep)
        e_ind = np.in1d(seqs_all, true_seqs)


        true_seqs = seqs_all[e_ind]
        true_rt.extend(rt_all[e_ind])
        true_md.extend(md_all[e_ind])
        true_rt = np.array(true_rt)
        true_md = np.array(true_md)

        # cPickle.dump([true_seqs, true_rt, true_md], open('/home/mark/test002.cpickle', 'w'))

        best_seq = defaultdict(list)
        newseqs = []
        newRTs = []
        for seq, RT in zip(true_seqs, true_rt):
            best_seq[seq].append(RT)
        for k, v in best_seq.items():
            newseqs.append(k)
            newRTs.append(np.mean(v))
        true_seqs = np.array(newseqs)
        true_rt = np.array(newRTs)

        # RC = achrom.get_RCs_vary_lcp(true_seqs, true_rt)
        RC, outmask = get_RCs2(true_seqs, true_rt)

        if elude_path:
            outtrain = open('/home/mark/Elude/train3.txt', 'w')
            ns = true_seqs[~outmask]
            nr = true_rt[~outmask]
            ll = len(ns)
            for seq, RT in zip(ns[:ll], nr[:ll]):
                outtrain.write(seq + '\t' + str(RT) + '\n')
            outtrain.close()

            subprocess.call([elude_path, '-t', '/home/mark/Elude/train3.txt', '-e', '/home/mark/Elude/train3.txt', '-a', '-g', '-o', '/home/mark/Elude/res3.txt'])
            pepdict = dict()
            train_RT = []
            train_seq = []
            for x in open('/home/mark/Elude/res3.txt').readlines()[3:]:
                seq, RT, RTexp = x.strip().split('\t')
                pepdict[seq] = float(RT)
                train_seq.append(seq)
                train_RT.append(float(RTexp))
            train_RT = np.array(train_RT)
            RT_pred = np.array([pepdict[s] for s in train_seq])#np.array([pepdict[s] for s in true_seqs[~outmask]])
            RT_diff = RT_pred - train_RT
            # print aux.linear_regression(RT_pred, train_RT)
            aa, bb, RR, ss = aux.linear_regression(RT_pred, train_RT)
        else:
            RC = achrom.get_RCs_vary_lcp(true_seqs[~outmask], true_rt[~outmask])
            # RC = cPickle.load(open('/home/mark/MS1_2016/confetti/output/c1_RC.pickle'))
            # RC = achrom.get_RCs(true_seqs, true_rt)
            RT_pred = np.array([achrom.calculate_RT(s, RC) for s in true_seqs])
            RT_diff = RT_pred - true_rt
            aa, bb, RR, ss = aux.linear_regression(RT_pred, true_rt)
            # aa, bb, RR, ss = aux.linear_regression(RT_pred, true_rt)
        print aa, bb, RR, ss

        # RC = achrom.get_RCs(true_seqs[~outmask], true_rt[~outmask])
        # # RC = cPickle.load(open('/home/mark/MS1_2016/confetti/output/c1_RC.pickle'))
        # # RC = achrom.get_RCs(true_seqs, true_rt)
        # RT_pred = np.array([achrom.calculate_RT(s, RC) for s in true_seqs])
        # RT_diff = RT_pred - true_rt
        # aa, bb, RR, ss = aux.linear_regression(RT_pred, true_rt)
        # # aa, bb, RR, ss = aux.linear_regression(RT_pred, true_rt)
        # print aa, bb, RR, ss
        best_sigma = ss#aux.linear_regression(RT_pred, true_rt)[3]
        RT_sigma = best_sigma
        # # RT_sigma = 9.6

        # e_ind = [True for i in range(len(true_seqs))]
        # better = []
        # for i in range(len(true_seqs)):
        #     e_ind[i] = False
        #     if i!=0:
        #         e_ind[i-1] = True
        # #     RC = achrom.get_RCs_vary_lcp(seqs[e_ind], RTs[e_ind])
        #     RC = achrom.get_RCs(true_seqs[e_ind], true_rt[e_ind])
        #     RT_pred = np.array([achrom.calculate_RT(s, RC) for s in true_seqs[e_ind]])
        #     if aux.linear_regression(RT_pred, true_rt[e_ind])[3] <= best_sigma:
        #         better.append(i)
        # # print better
        # e_ind[-1] = True
        # for i in better:
        #     e_ind[i] = False
            
        # RC = achrom.get_RCs(true_seqs[e_ind], true_rt[e_ind])
        # RT_pred = np.array([achrom.calculate_RT(s, RC) for s in true_seqs[e_ind]])
        # print aux.linear_regression(RT_pred, true_rt[e_ind])
        # RT_sigma = aux.linear_regression(RT_pred, true_rt[e_ind])[3]




        # RC = achrom.get_RCs_vary_lcp(true_seqs, true_rt)
        # rt_pred = np.array([achrom.calculate_RT(s, RC) for s in true_seqs])
        # print aux.linear_regression(rt_pred, true_rt)

# #####
#         rt_pred = np.array([achrom.calculate_RT(s, RC) for s in seqs_all])
#         rt_diff = rt_all - rt_pred

#         # RT_m = settings.getfloat('search', 'retention time shift')
#         # RT_sigma = settings.getfloat('search', 'retention time sigma')

#         e_all = (rt_diff) ** 2 / (RT_sigma ** 2)
#         r = settings.getfloat('search', 'r threshold') ** 2
#         e_ind = e_all <= r
#         seqs_all = seqs_all[e_ind]
#         md_all = md_all[e_ind]
#         rt_all = rt_all[e_ind]
#         ids_all = ids_all[e_ind]





    else:
        print 'No matches found'

    p1 = set(seqs_all)

    pepdict = dict()    
    if elude_path:
        outtrain = open('/home/mark/Elude/test3.txt', 'w')
        for seq in p1:
            outtrain.write(seq + '\n')
        outtrain.close()

        subprocess.call([elude_path, '-t', '/home/mark/Elude/train3.txt', '-e', '/home/mark/Elude/test3.txt', '-a', '-o', '/home/mark/Elude/res3.txt'])
        for x in open('/home/mark/Elude/res3.txt').readlines()[3:]:
            seq, RT = x.strip().split('\t')
            pepdict[seq] = float(RT)

    else:
        for seq in p1:
            pepdict[seq] = achrom.calculate_RT(s, RC)
        # for rrr in np.arange(2.0, 0.5, -0.1):
        # for rrr in [1.6, ]:
            # rt_pred = np.array([achrom.calculate_RT(s, RC) for s in seqs_all])
    rt_pred = np.array([pepdict[s] for s in seqs_all])
    rt_diff = rt_all - rt_pred
    e_all = (rt_diff) ** 2 / (RT_sigma ** 2)
    r = settings.getfloat('search', 'r threshold') ** 2
    e_ind = e_all <= r
    seqs_all = seqs_all[e_ind]
    md_all = md_all[e_ind]
    rt_all = rt_all[e_ind]
    ids_all = ids_all[e_ind]



    with open(os.path.splitext(fname)[0] + '_PFMs.csv', 'w') as output:
        output.write('sequence\tmass diff\tRT diff\tpeak_id\tproteins\n')
        for seq, md, rtd, peak_id in zip(seqs_all, md_all, rt_all, ids_all):
            output.write('\t'.join((seq, str(md), str(rtd), str(peak_id), ';'.join(pept_prot[seq]))) + '\n')

    p1 = set(seqs_all)
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

        prots_spc_copy = copy(prots_spc)
        top100decoy_score = [prots_spc.get(dprot, 0) for dprot in protsN if isdecoy_key(dprot)]
        top100decoy_N = [val for key, val in protsN.items() if isdecoy_key(key)]
        p = np.mean(top100decoy_score) / np.mean(top100decoy_N)
        print 'p=%s' % (np.mean(top100decoy_score) / np.mean(top100decoy_N))

        prots_spc = dict()
        all_pvals = calc_sf_all(v_arr, n_arr, p)
        for idx, k in enumerate(names_arr):
            prots_spc[k] = all_pvals[idx]

        sortedlist_spc = sorted(prots_spc.iteritems(), key=operator.itemgetter(1))[::-1]
        with open(os.path.splitext(fname)[0] + '_proteins_full.csv', 'w') as output:
            output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\n')
            for x in sortedlist_spc:
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

        filtered_prots = aux.filter(prots_spc.items(), fdr=0.01, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1, full_output=True)

        identified_proteins = 0
        identified_proteins_valid = 0

        for x in filtered_prots:
            identified_proteins += 1

        print 'TOP 5 identified proteins:'
        print 'dbname\tscore\tnum matched peptides\tnum theoretical peptides'
        for x in filtered_prots[:5]:
            print '\t'.join((x[0], x[1], int(prots_spc_copy[x[0]]), protsN[x[0]]))
        print 'results:%s;number of identified proteins = %d' % (fname, identified_proteins, )
        print 'R=', r
        with open(os.path.splitext(fname)[0] + '_proteins.csv', 'w') as output:
            output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\n')
            for x in filtered_prots:
                output.write('\t'.join((x[0], str(x[1]), str(prots_spc_copy[x[0]]), str(protsN[x[0]]))) + '\n')

    else:
        print 'No matches found'
