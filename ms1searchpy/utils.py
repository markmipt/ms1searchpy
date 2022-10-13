from pyteomics import fasta, parser
import os
from scipy.stats import binom
import numpy as np
import pandas as pd
import random
import itertools
from biosaur2 import main as bio_main
import logging

logger = logging.getLogger(__name__)


def recalc_spc(banned_dict, unstable_prots, prots_spc2):
    tmp = dict()
    for k in unstable_prots:
        tmp[k] = sum(banned_dict.get(l, 1) > 0 for l in prots_spc2[k])
    return tmp

def iterate_spectra(fname, min_ch, max_ch, min_isotopes, min_scans, nproc):
    if os.path.splitext(fname)[-1].lower() == '.mzml':
        args = {
            'file': fname,
            'mini': 1,
            'minmz': 350,
            'maxmz': 1500,
            'pasefmini': 100,
            'htol': 8,
            'itol': 8,
            'paseftol': 0.05,
            'nm': 0,
            'o': '',
            'hvf': 1.3,
            'minlh': 2,
            'pasefminlh': 1,
            'nprocs': nproc,
            'cmin': 1,
            'cmax': 6,
            'dia': False,
            'diahtol': 25,
            'diaminlh': 1,
            'mgf': '',
            'tof': False,
            'profile': False,
            'write_hills': False,
            'debug': False  # actual debug value is set through logging, not here
        }
        bio_main.process_file(args)
        fname = os.path.splitext(fname)[0] + '.features.tsv'

    df_features = pd.read_csv(fname, sep='\t')

    required_columns = [
        'nIsotopes',
        'nScans',
        'charge',
        'massCalib',
        'rtApex',
        'mz',
        ]

    if not all(req_col in df_features.columns for req_col in required_columns):
        logger.error('input feature file have missing columns: %s', ';'.join([req_col for req_col in required_columns if req_col not in df_features.columns]))
        raise Exception('Exception: wrong columns in feature file')
    logger.info('Total number of peptide isotopic clusters: %d', len(df_features))

    if 'id' not in df_features.columns:
        df_features['id'] = df_features.index
    if 'FAIMS' not in df_features.columns:
        df_features['FAIMS'] = 0
    if 'ion_mobility' not in df_features.columns:
        df_features['ion_mobility'] = 0

    # if 'mz_std_1' in df_features.columns:
    #     df_features['mz_diff_ppm_1'] = df_features.apply(lambda x: 1e6 * (x['mz'] - (x['mz_std_1'] - 1.00335 / x['charge'])) / x['mz'], axis=1)
    #     df_features['mz_diff_ppm_2'] = -100
    #     df_features.loc[df_features['intensity_2'] > 0, 'mz_diff_ppm_2'] = df_features.loc[df_features['intensity_2'] > 0, :].apply(lambda x: 1e6 * (x['mz'] - (x['mz_std_2'] - 2 * 1.00335 / x['charge'])) / x['mz'], axis=1)

    #     df_features['I-0-1'] = df_features.apply(lambda x: x['intensityApex'] / x['intensity_1'], axis=1)
    #     df_features['I-0-2'] = -1
    #     df_features.loc[df_features['intensity_2'] > 0, 'I-0-2'] = df_features.loc[df_features['intensity_2'] > 0, :].apply(lambda x: x['intensityApex'] / x['intensity_2'], axis=1)

    # Check unique ids
    if len(df_features['id']) != len(set(df_features['id'])):
        df_features['id'] = df_features.index + 1

    # Remove features with low number of isotopes
    df_features = df_features[df_features['nIsotopes'] >= min_isotopes]

    # Remove features with low number of Scans
    df_features = df_features[df_features['nScans'] >= min_scans]

    # Remove features using min and max charges
    df_features = df_features[df_features['charge'].apply(lambda x: min_ch <= x <= max_ch)]

    return df_features

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
    return convert_tandem_cleave_rule_to_regexp(enzyme)
    # if enzyme in parser.expasy_rules:
    #     return parser.expasy_rules.get(enzyme)
    # else:
    #     try:
    #         enzyme = convert_tandem_cleave_rule_to_regexp(enzyme)
    #         return enzyme
    #     except:
    #         return enzyme

def prot_gen(args):
    db = args['d']

    with fasta.read(db) as f:
        for p in f:
            yield p

def prepare_decoy_db(args):
    add_decoy = args['ad']
    if add_decoy:

        prefix = args['prefix']
        db = args['d']
        out1, out2 = os.path.splitext(db)
        out_db = out1 + '_shuffled' + out2
        logger.info('Creating decoy database: %s', out_db)

        extra_check = False
        if '{' in args['e']:
            extra_check = True
        if extra_check:
            banned_pairs = set()
            banned_aa = set()
            for enzyme_local in args['e'].split(','):
                if '{' in enzyme_local:
                    lpart, rpart = enzyme_local.split('|')
                    for aa_left, aa_right in itertools.product(lpart[1:-1], rpart[1:-1]):
                        banned_aa.add(aa_left)
                        banned_aa.add(aa_right)
                        banned_pairs.add(aa_left+aa_right)

            logger.debug(banned_aa)
            logger.debug(banned_pairs)

        enzyme = get_enzyme(args['e'])
        cleave_rule_custom = enzyme + '|' + '([BXZUO])'
        # cleave_rule_custom = '([RKBXZUO])'
        logger.debug(cleave_rule_custom)

        shuf_map = dict()

        prots = []

        for p in fasta.read(db):
            if not p[0].startswith(prefix):
                target_peptides = [x[1] for x in parser.icleave(p[1], cleave_rule_custom, 0)]

                checked_peptides = set()
                sample_list = []
                for idx, pep in enumerate(target_peptides):

                    if len(pep) > 2:
                        pep_tmp = pep[1:-1]
                        if extra_check:
                            for bp in banned_pairs:
                                if bp in pep_tmp:
                                    pep_tmp = pep_tmp.replace(bp, '')
                                    checked_peptides.add(idx)


                        sample_list.extend(pep_tmp)
                random.shuffle(sample_list)
                idx_for_shuffle = 0

                decoy_peptides = []
                for idx, pep in enumerate(target_peptides):

                    if len(pep) > 2:

                        if pep in shuf_map:
                            tmp_seq = shuf_map[pep]
                        else:
                            if not extra_check or idx not in checked_peptides:
                                tmp_seq = pep[0]
                                for pep_aa in pep[1:-1]:
                                    tmp_seq += sample_list[idx_for_shuffle]
                                    idx_for_shuffle += 1
                                tmp_seq += pep[-1]
                            else:
                                max_l = len(pep)
                                tmp_seq = ''
                                ii = 0
                                while ii < max_l - 1:
                                # for ii in range(max_l-1):
                                    if pep[ii] in banned_aa and pep[ii+1] in banned_aa and pep[ii] + pep[ii+1] in banned_pairs:
                                        tmp_seq += pep[ii] + pep[ii+1]
                                        ii += 1
                                    else:
                                        if ii == 0:
                                            tmp_seq += pep[ii]
                                        else:
                                            tmp_seq += sample_list[idx_for_shuffle]
                                            idx_for_shuffle += 1

                                    ii += 1
                                tmp_seq += pep[max_l-1]

                            shuf_map[pep] = tmp_seq
                    else:
                        tmp_seq = pep

                    decoy_peptides.append(tmp_seq)

                assert len(target_peptides) == len(decoy_peptides)

                prots.append((p[0], ''.join(target_peptides)))
                prots.append(('DECOY_' + p[0], ''.join(decoy_peptides)))

        fasta.write(prots, open(out_db, 'w')).close()
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

    target_prot_count = 0
    decoy_prot_count = 0
    target_peps = set()
    decoy_peps = set()

    for desc, prot in prot_gen(args):
        dbinfo = desc.split(' ')[0]
        if dbinfo.startswith(prefix):
            decoy_prot_count += 1
        else:
            target_prot_count += 1
        for pep in prot_peptides(prot, enzyme, mc, minlen, maxlen, desc.startswith(prefix), dont_use_seen_peptides=True):
            pept_prot.setdefault(pep, []).append(dbinfo)
            protsN.setdefault(dbinfo, set()).add(pep)
    for k, v in protsN.items():
        if k.startswith(prefix):
            decoy_peps.update(v)
        else:
            target_peps.update(v)

        protsN[k] = len(v)

    logger.info('Database information:')
    logger.info('Target/Decoy proteins: %d/%d', target_prot_count, decoy_prot_count)
    logger.info('Target/Decoy peptides: %d/%d', len(target_peps), len(decoy_peps))
    logger.info('Target-Decoy peptide intersection: %.1f %%',
        100 * len(target_peps.intersection(decoy_peps)) / (len(target_peps) + len(decoy_peps)))
    del decoy_peps
    del target_peps
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

def calc_sf_all(v, n, p, prev_best_score=False):
    sf_values = -np.log10(binom.sf(v-1, n, p))
    sf_values[np.isnan(sf_values)] = 0
    sf_values[np.isinf(sf_values)] = (prev_best_score if prev_best_score is not False else max(sf_values[~np.isinf(sf_values)]) * 2)
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
