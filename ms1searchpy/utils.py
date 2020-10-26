from pyteomics import fasta, parser
from multiprocessing import Queue, Process, cpu_count
import os
import csv
import subprocess
from scipy.stats import binom
import numpy as np

def recalc_spc(banned_dict, unstable_prots, prots_spc2):
    tmp = dict()
    for k in unstable_prots:
        tmp[k] = sum(banned_dict.get(l, 1) > 0 for l in prots_spc2[k])
    return tmp

def iterate_spectra(fname, min_ch, max_ch, min_isotopes, min_scans):
    if os.path.splitext(fname)[-1].lower() == '.mzml':
        subprocess.call(['biosaur', fname])
        # advpath = '--advParams=' + os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Dinosaur/adv.txt')
        # subprocess.call(['java', '-Djava.awt.headless=true', '-jar', os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Dinosaur/Dinosaur-1.1.3.free.jar'), advpath, '--concurrency=12', fname])
        fname = os.path.splitext(fname)[0] + '.features.tsv'
    with open(fname, 'r') as infile:
        csvreader = csv.reader(infile, delimiter='\t')
        header = next(csvreader)#.next()

        mass_ind = header.index('massCalib')
        # mass_ind = header.index('mass')
        # mass_ind = header.index('mostAbundantMz')
        RT_ind = header.index('rtApex')
        ch_ind = header.index('charge')
        nIsotopes_ind = header.index('nIsotopes')
        Int_ind = header.index('intensityApex')
        nScans_ind = header.index('nScans')

        # mass_ind = header.index('apex_mz')
        # RT_ind = header.index('rt')
        # ch_ind = header.index('charge')
        # nIsotopes_ind = header.index('isotopes_count')
        # Int_ind = header.index('intensity')
        # nScans_ind = header.index('centroid_mz')

        try:
            mz_ind = header.index('mz')
        except:
            mz_ind = -1


        # try:
        #     FAIMS_ind = header.index('FAIMS')
        # except:
        #     FAIMS_ind = -1


        try:
            FAIMS_ind = header.index('FAIMS')
        except:
            FAIMS_ind = -1

        try:
            av_ind = header.index('averagineCorr')
            # av_ind = header.index('intensity_experimental')
        except:
            av_ind = -1


        try:
            mz3_ind = header.index('mz3')
        except:
            mz3_ind = -1

        idx = 0
        for z in csvreader:
            nm = float(z[mass_ind])
            RT = float(z[RT_ind])
            ch = float(z[ch_ind])
            # if ch != 0:
            #     nm = (nm + 1.0073 * ch) / ch
            # else:
            # nm = (nm + 1.0073 * ch) / ch
            # for ch in [1,2,3]:
            # nm = nm * ch - 1.0073 * ch
            nIsotopes = float(z[nIsotopes_ind])# + 1
            nScans = int(z[nScans_ind])
            # nIsotopes = 4
            # nScans = 4
            I = float(z[Int_ind])
            mz = float(z[mz_ind]) if mz_ind >= 0 else 0
            mz3 = float(z[mz3_ind]) if mz3_ind >= 0 else 0
            im = float(z[FAIMS_ind]) if FAIMS_ind >= 0 else 0
            # av = float(z[av_ind]) if av_ind >= 0 else 0
            av = str(z[av_ind])
            idx += 1
            if nIsotopes >= min_isotopes and min_ch <= ch <= max_ch and min_scans <= nScans:
                yield nm, RT, ch, idx, I, nScans, nIsotopes, mz, av, im
                # yield nm, RT, ch, idx, I, nScans, nIsotopes, mz, av, mz3


def peptide_gen(args):
    prefix = args['prefix']
    enzyme = get_enzyme(args['e'])
    mc = args['mc']
    minlen = args['lmin']
    maxlen = args['lmax']
    for prot in prot_gen(args):
        for pep in prot_peptides(prot[1], enzyme, mc, minlen, maxlen, is_decoy=prot[0].startswith(prefix)):
            yield pep

def get_enzyme(enzyme):
    if enzyme in parser.expasy_rules:
        return parser.expasy_rules.get(enzyme)
    else:
        try:
            enzyme = convert_tandem_cleave_rule_to_regexp(enzyme)
            return enzyme
        except:
            return enzyme

def prot_gen(args):
    db = args['d']
    add_decoy = args['ad']
    prefix = args['prefix']

    read = [fasta.read, lambda f: fasta.decoy_db(f, mode='shuffle', prefix=prefix)][add_decoy]
    with read(db) as f:
        for p in f:
            yield p

def prepare_decoy_db(args):
    add_decoy = args['ad']
    if add_decoy:
        prefix = args['prefix']
        db = args['d']
        out1, out2 = os.path.splitext(db)
        out_db = out1 + '_shuffled' + out2
        print(out_db)
        fasta.write_decoy_db(db, open(out_db, 'w'), mode='shuffle', prefix=prefix).close()
        args['d'] = out_db
        args['ad'] = 0
    return args

seen_target = set()
seen_decoy = set()
def prot_peptides(prot_seq, enzyme, mc, minlen, maxlen, is_decoy, dont_use_seen_peptides=False):

    dont_use_fast_valid = parser.fast_valid(prot_seq)
    peptides = parser.cleave(prot_seq, enzyme, mc)
    for pep in peptides:
        plen = len(pep)
        if minlen <= plen <= maxlen:
            forms = []
            if dont_use_fast_valid or pep in seen_target or pep in seen_decoy or parser.fast_valid(pep):
                if plen <= maxlen:
                    forms.append(pep)
            for f in forms:
                if dont_use_seen_peptides:
                    yield f
                else:
                    if f not in seen_target and f not in seen_decoy:
                        if is_decoy:
                            seen_decoy.add(f)
                        else:
                            seen_target.add(f)
                        yield f

def get_prot_pept_map(args):
    seen_target.clear()
    seen_decoy.clear()


    prefix = args['prefix']
    enzyme = get_enzyme(args['e'])
    mc = args['mc']
    minlen = args['lmin']
    maxlen = args['lmax']

    pept_prot = dict()
    protsN = dict()

    for desc, prot in prot_gen(args):
        dbinfo = desc.split(' ')[0]
        for pep in prot_peptides(prot, enzyme, mc, minlen, maxlen, desc.startswith(prefix), dont_use_seen_peptides=True):
            pept_prot.setdefault(pep, []).append(dbinfo)
            protsN.setdefault(dbinfo, set()).add(pep)
    for k, v in protsN.items():
        protsN[k] = len(v)
    return protsN, pept_prot

def convert_tandem_cleave_rule_to_regexp(cleavage_rule):

    def get_sense(c_term_rule, n_term_rule):
        if '{' in c_term_rule:
            return 'N'
        elif '{' in n_term_rule:
            return 'C'
        else:
            if len(c_term_rule) <= len(n_term_rule):
                return 'C'
            else:
                return 'N'

    def get_cut(cut, no_cut):
        aminoacids = set(parser.std_amino_acids)
        cut = ''.join(aminoacids & set(cut))
        if '{' in no_cut:
            no_cut = ''.join(aminoacids & set(no_cut))
            return cut, no_cut
        else:
            no_cut = ''.join(set(parser.std_amino_acids) - set(no_cut))
            return cut, no_cut

    out_rules = []
    for protease in cleavage_rule.split(','):
        protease = protease.replace('X', ''.join(parser.std_amino_acids))
        c_term_rule, n_term_rule = protease.split('|')
        sense = get_sense(c_term_rule, n_term_rule)
        if sense == 'C':
            cut, no_cut = get_cut(c_term_rule, n_term_rule)
        else:
            cut, no_cut = get_cut(n_term_rule, c_term_rule)

        if no_cut:
            if sense == 'C':
                out_rules.append('([%s](?=[^%s]))' % (cut, no_cut))
            else:
                out_rules.append('([^%s](?=[%s]))' % (no_cut, cut))
        else:
            if sense == 'C':
                out_rules.append('([%s])' % (cut, ))
            else:
                out_rules.append('(?=[%s])' % (cut, ))
    return '|'.join(out_rules)

def multimap(n, func, it, **kw):
    for s in it:
        yield func(s, **kw)

def keywithmaxval(d):
     """ a) create a list of the dict's keys and values; 
         b) return the key with the max value"""  
     v=list(d.values())
     k=list(d.keys())
     return k[v.index(max(v))]

def calc_sf_all(v, n, p):
    sf_values = -np.log10(binom.sf(v, n, p))
    sf_values[np.isinf(sf_values)] = 0
    return sf_values


# def multimap(n, func, it, **kw):
#     if n == 0:
#         try:
#             n = cpu_count()
#         except NotImplementedError:
#             n = 1
#     # if n == 1:
#     #     for s in it:
#     #         result = func(s, best_res, **kw)
#     #         if result:
#     #             for x in result:
#     #                 peptide, m, snp_label, res = x

#     #                 for score, spec_t, c, info in res:
#     #                     if -score <= best_res.get(spec_t, 0):
#     #                         best_res_raw[spec_t] = [peptide, m, snp_label, score, spec_t, c, info]
#     #                         best_res[spec_t] = -score   
#     #     return best_res_raw, best_res
      
#     else:

#         qout = Queue()
#         count = 0

#         while True:
#             qin = list(islice(it, 5000000))
#             if not len(qin):
#                 break
# #           print 'Loaded 500000 items. Ending cycle.'
#             procs = []
#             for proc_num in range(n):
#                 p = Process(target=worker, args=(qin, qout, proc_num, n, best_res, best_res_raw))
#                 p.start()
#                 procs.append(p)

#             count = len(qin)

#             for _ in range(n):
#                 for item in iter(qout.get, None):
#                     for k, v in item.items():
#                         if -v[3] <= best_res.get(k, 0):
#                             best_res_raw[k] = v
#                             best_res[k] = -v[3]
#                     # yield item

#             for p in procs:
#                 p.join()

#         return best_res_raw, best_res