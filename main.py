import os
import cPickle
import utils
import numpy as np
from scipy.stats import binom
import operator
from copy import copy
from collections import defaultdict
import re
from pyteomics import parser, mass, fasta, auxiliary as aux, achrom
try:
    from pyteomics import cmass
except ImportError:
    cmass = mass


def process_file(fname, settings):
    ftype = fname.rsplit('.', 1)[-1].lower()
    utils.seen_target.clear()
    utils.seen_decoy.clear()
    return process_peptides(fname, settings)

def settings(fname=None, default_name=os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'default.cfg')):
    """Read a configuration file and return a :py:class:`RawConfigParser` object.
    """

    raw_config = utils.CustomRawConfigParser(dict_type=dict, allow_no_value=True)
    if default_name:
        raw_config.read(default_name)
    if fname:
        raw_config.read(fname)
    return raw_config


def peptide_processor(peptide, **kwargs):
    seqm = peptide
    m = cmass.fast_mass(seqm, aa_mass=kwargs['aa_mass'])
    acc_l = kwargs['acc_l']
    acc_r = kwargs['acc_r']
    settings = kwargs['settings']
    dm_l = acc_l * m / 1.0e6
    dm_r = acc_r * m / 1.0e6
    start = nmasses.searchsorted(m - dm_l)
    end   = nmasses.searchsorted(m + dm_r)
    idx = set(range(start, end))

    if idx:
        RC = kwargs['RC']
        RT_sigma = kwargs['RT_sigma']
        RT = achrom.calculate_RT(seqm, RC)

    results = []
    for i in idx:
        RTdiff = RT - rts[i]
        if abs(RTdiff) <= 3 * RT_sigma:
            intensity = Is[i]
            massdiff = (m - nmasses[i]) / m * 1e6
            results.append((seqm, massdiff, RTdiff, intensity))
    return results


def prepare_peptide_processor(fname, settings):
    global nmasses
    global rts
    global Is
    global charges
    nmasses = []
    rts = []
    Is = []
    charges = []

    min_ch = settings.getint('search', 'minimum charge')
    max_ch = settings.getint('search', 'maximum charge')
    min_i = settings.getint('search', 'intensity threshold')
    nprocs = settings.getint('performance', 'processes')

    print 'Reading spectra ...'
    for m, RT, I, c, peak_id, pI in utils.iterate_spectra(fname, min_ch, max_ch, min_i, nprocs):
        nmasses.append(m)
        rts.append(RT)
        Is.append(I)
        charges.append(c)

    i = np.argsort(nmasses)
    nmasses = np.array(nmasses)[i]
    rts = np.array(rts)[i]
    Is = np.array(Is)[i]
    charges = np.array(charges)[i]

    fmods = settings.get('modifications', 'fixed')
    aa_mass = mass.std_aa_mass
    if fmods:
        for mod in re.split(r'[,;]\s*', fmods):
            m, aa = parser._split_label(mod)
            aa_mass[aa] += settings.getfloat('modifications', m)


    acc_l = settings.getfloat('search', 'precursor accuracy left')
    acc_r = settings.getfloat('search', 'precursor accuracy right')

    RC = cPickle.load(open(settings.get('input', 'RC'), 'r'))
    RT_sigma = settings.getfloat('search', 'retention time sigma')

    return {'aa_mass': aa_mass, 'acc_l': acc_l, 'acc_r': acc_r,
            'RC': RC, 'RT_sigma': RT_sigma, 'settings': settings}

def peptide_processor_iter_isoforms(peptide, **kwargs):
    out = []
    out.append(peptide_processor(peptide, **kwargs))
    return out


def process_peptides(fname, settings):
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

    with open(os.path.splitext(fname)[0] + '_PFMs.csv', 'w') as output:
        output.write('sequence\mass diff\tRT diff\tintensity\tproteins\n')
        for seq, md, rtd, intensity in ms1results:
            output.write('\t'.join((seq, str(md), str(rtd), str(intensity), ';'.join(pept_prot[seq]))) + '\n')

    seqs_all, md_all, rt_all = zip(*ms1results)[:3]
    seqs_all = np.array(seqs_all)
    md_all = np.array(md_all)
    rt_all = np.array(rt_all)
    del ms1results

    mass_m = settings.getfloat('search', 'precursor accuracy shift')
    mass_sigma = settings.getfloat('search', 'precursor accuracy sigma')
    RT_m = settings.getfloat('search', 'retention time shift')
    RT_sigma = settings.getfloat('search', 'retention time sigma')

    e_all = (md_all - mass_m) ** 2 / (mass_sigma ** 2) + (rt_all - RT_m) ** 2 / (RT_sigma ** 2)

    protsV = set()
    path_to_valid_fasta = settings.get('input', 'valid proteins')
    if path_to_valid_fasta:
        for prot in fasta.read(path_to_valid_fasta):
            protsV.add(prot[0].split(' ')[0])

    def calc_sf_all(v, n, p):
        sf_values = np.log10(1 / binom.sf(v, n, p))
        sf_values[np.isinf(sf_values)] = 1
        return sf_values

    r = settings.getfloat('search', 'r threshold') ** 2
    # for zzz in np.arange(1.0, 2.0, 0.1):
    #     r = zzz**2
    p1 = set(seqs_all[e_all <= r])

    if len(p1):
        prots_spc2 = defaultdict(set)
        for pep, proteins in pept_prot.iteritems():
            if pep in p1:
                for protein in proteins:
                    prots_spc2[protein].add(pep)

        isdecoy = lambda x: x[0].startswith(prefix)
        isdecoy_key = lambda x: x.startswith(prefix)
        escore = lambda x: -x[1]

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
            if x[0] in protsV:
                identified_proteins_valid += 1
            identified_proteins += 1

        for x in filtered_prots[:5]:
            print x[0], x[1], int(prots_spc_copy[x[0]]), protsN[x[0]]
        print 'results:%s;number of identified proteins = %d;number of valid proteins = %d' % (fname, identified_proteins, identified_proteins_valid)
        print 'R=', r
        with open(os.path.splitext(fname)[0] + '_proteins.csv', 'w') as output:
            output.write('dbname\tscore\tmatched peptides\ttheoretical peptides\n')
            for x in filtered_prots:
                output.write('\t'.join((x[0], str(x[1]), str(prots_spc_copy[x[0]]), str(protsN[x[0]]))) + '\n')

    else:
        print 'No matches found'
