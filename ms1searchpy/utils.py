from pyteomics import fasta, parser, mass
import os
from scipy.stats import binom
import numpy as np
import pandas as pd
import random
import itertools
from biosaur2 import main as bio_main
import logging
from copy import deepcopy

logger = logging.getLogger(__name__)

# Temporary for pyteomics <= Version 4.5.5 bug
if 'H-' in mass.std_aa_mass:
    del mass.std_aa_mass['H-']
if '-OH' in mass.std_aa_mass:
    del mass.std_aa_mass['-OH']

mods_custom_dict = {
    'Oxidation': 15.994915,
    'Carbamidomethyl': 57.021464,
    'TMT6plex': 229.162932,
}


def get_aa_mass_with_fixed_mods(fmods, fmods_legend):

    if fmods_legend:
        for mod in fmods_legend.split(','):
            psiname, m = mod.split('@')
            mods_custom_dict[psiname] = float(m)

    aa_mass = deepcopy(mass.std_aa_mass)
    aa_to_psi = dict()

    mass_h2o = mass.calculate_mass('H2O')
    for k in list(aa_mass.keys()):
        aa_mass[k] = round(mass.calculate_mass(sequence=k) - mass_h2o, 7)

    if fmods:
        for mod in fmods.split(','):
            psiname, aa = mod.split('@')
            if psiname not in mods_custom_dict:
                logger.error('PSI Name for modification %s is missing in the modification legend' % (psiname, ))
                raise Exception('Exception: missing PSI Name for modification')
            if aa == '[':
                aa_mass['Nterm'] = float(mods_custom_dict[psiname])#float(m)
                aa_to_psi['Nterm'] = psiname
            elif aa == ']':
                aa_mass['Cterm'] = float(mods_custom_dict[psiname])#float(m)
                aa_to_psi['Cterm'] = psiname
            else:
                aa_mass[aa] += float(mods_custom_dict[psiname])#float(m)
                aa_to_psi[aa] = psiname

    logger.debug(aa_mass)

    return aa_mass, aa_to_psi


def mods_for_deepLC(seq, aa_to_psi):
    if 'Nterm' in aa_to_psi:
        mods_list = ['0|%s' % (aa_to_psi['Nterm'], ), ]
    else:
        mods_list = []
    mods_list.extend([str(idx+1)+'|%s' % (aa_to_psi[aa]) for idx, aa in enumerate(seq) if aa in aa_to_psi])
    if 'Cterm' in aa_to_psi:
        mods_list.append(['-1|%s' % (aa_to_psi['Cterm'], ), ])
    return '|'.join(mods_list)

def recalc_spc(banned_dict, unstable_prots, prots_spc2):
    tmp = dict()
    for k in unstable_prots:
        tmp[k] = sum(banned_dict.get(l, 1) > 0 for l in prots_spc2[k])
    return tmp

def iterate_spectra(fname, min_ch, max_ch, min_isotopes, min_scans, nproc, check_unique=True):
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

    if check_unique:
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
        for pep in prot_peptides(prot, enzyme, mc, minlen, maxlen, desc.startswith(prefix), dont_use_seen_peptides=True):
            pept_prot.setdefault(pep, []).append(dbinfo)
            protsN.setdefault(dbinfo, set()).add(pep)
    for k, v in protsN.items():
        if k.startswith(prefix):
            decoy_prot_count += 1
            decoy_peps.update(v)
        else:
            target_prot_count += 1
            target_peps.update(v)

        protsN[k] = len(v)

    logger.info('Database information:')
    logger.info('Target/Decoy proteins: %d/%d', target_prot_count, decoy_prot_count)
    target_peps_number = len(target_peps)
    decoy_peps_number = len(decoy_peps)
    intersection_number = len(target_peps.intersection(decoy_peps)) / (target_peps_number + decoy_peps_number)
    logger.info('Target/Decoy peptides: %d/%d', target_peps_number, decoy_peps_number)
    logger.info('Target-Decoy peptide intersection: %.1f %%',
        100 * intersection_number)
       
    ml_correction = decoy_peps_number * (1 - intersection_number) / target_peps_number * 0.5
    del decoy_peps
    del target_peps
    return protsN, pept_prot, ml_correction


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

