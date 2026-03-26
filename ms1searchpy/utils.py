from pyteomics import fasta, parser, mass, mzml
import os
from scipy.stats import binom
import numpy as np
import pandas as pd
import random
import itertools
from biosaur2 import main as bio_main
import logging
from copy import deepcopy
from collections import defaultdict, Counter
import string
from time import strftime
from importlib.metadata import version
from copy import copy
from lxml import etree

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

        print(mods_custom_dict)

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




def write_pepxml(inputfile, args, df1, pept_prot):

    outpath = args['o']

    if not args['ms2mzml']:
        df1['b_count'] = 0
        df1['y_count'] = 0
        df1['hyperscore'] = 0
        df1['hyperscore3'] = 0
        df1['sf'] = 0

    if 'mc' not in df1.columns:
        df1['mc'] = 0

    df1['matched_ions'] = df1['b_count'] + df1['y_count']

    # ms2utils.set_mod_dict(settings)

    enzyme = args['e']#settings.get('search', 'enzyme')
    search_engine = 'ms1searchpy'
    database = args['d']#settings.get('input', 'database')
    missed_cleavages = args['mc']#settings.getint('search', 'number of missed cleavages')
    fmods = args['fmods']#settings.get('modifications', 'fixed')
    snp = 0
    # snp = settings.getint('search', 'snp')
    # nterm_mass = settings.getfloat('modifications', 'protein nterm cleavage')
    # cterm_mass = settings.getfloat('modifications', 'protein cterm cleavage')
    # tags = get_tags(settings.get('output', 'tags'))

    # nterm_fixed = 0
    # cterm_fixed = 0

    # for mod in re.split(r'[,;]\s*', fmods):
    #     if mod.startswith('-'):
    #         cterm_fixed = settings.getfloat('modifications', 'protein cterm cleavage')
    #     elif mod.endswith('-'):
    #         nterm_fixed = settings.getfloat('modifications', 'protein nterm cleavage')


    filename = inputfile.replace('.mzML', '.pep.xml')#get_outpath(inputfile, settings, 'pep.xml')
    with open(filename, 'wb') as output:
        logger.info('Writing %s ...', filename)
        line1 = b'<?xml version="1.0" encoding="UTF-8"?>\n\
        <?xml-stylesheet type="text/xsl" href="pepXML_std.xsl"?>\n'
        output.write(line1)

        base_name, ftype = os.path.splitext(inputfile)
        ftype = ftype.lower()

        root = etree.Element('msms_pipeline_analysis')
        root.set("date", strftime("%Y:%m:%d:%H:%M:%S"))
        root.set("summary_xml", '')
        root.set("xmlns", 'http://regis-web.systemsbiology.net/pepXML')
        # TODO
        #root.set("xmlns:xsi", 'http://www.w3.org/2001/XMLSchema-instance')
        #root.set("xsi:schemaLocation", 'http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v117.xsd')

        child1 = etree.Element('msms_run_summary')
        child1.set("base_name", base_name)
        child1.set("search_engine", search_engine)
        child1.set("raw_data_type", "raw")  # ?

        if ftype == '.mgf':
            child1.set("raw_data", ".mgf")
        elif ftype == '.mzml':
            child1.set("raw_data", ".mzML")
        else:
            child1.set("raw_data", ".?")
        root.append(child1)

        child2 = etree.Element('sample_enzyme')
        child2.set('name', enzyme)
        child1.append(child2)

        child3 = etree.Element('specificity')
        child3.set("cut", "KR")
        child3.set("no_cut", "P")
        child3.set("sense", "C")

        child2.append(child3)

        child4 = etree.Element('search_summary')
        child4.set('base_name', base_name)
        child4.set('search_engine', search_engine)
        child4.set("search_engine_version", version('ms1searchpy'))
        child4.set('precursor_mass_type', 'monoisotopic')
        child4.set('fragment_mass_type', 'monoisotopic')
        child4.set('search_id', '1')

        # for child_mod in get_child_for_mods(settings.get('modifications', 'fixed'), settings, fixed=True):
        #     child4.append(child_mod)
        # for child_mod in get_child_for_mods(settings.get('modifications', 'variable_original'), settings, fixed=False):
        #     child4.append(child_mod)
        # for child_mod in get_child_for_mods(settings.get('modifications', 'protein_original'), settings, fixed=False, protein=True):
        #     child4.append(child_mod)

        child1.append(child4)

        child5 = etree.Element('search_database')
        child5.set('local_path', database)
        child5.set('type', 'AA')

        child4.append(copy(child5))

        child5 = etree.Element('enzymatic_search_constraint')
        child5.set('enzyme', enzyme)
        child5.set('max_num_internal_cleavages', str(missed_cleavages))
        child5.set('min_number_termini', '2')

        child4.append(copy(child5))

#         results = [x for x in results if x['candidates'].size]
# #       results = list(get_output(results, settings))
        logger.info('Accumulated results: %s', len(df1))
#         pept_prot, prots, pept_neighbors, pept_ntts = build_pept_prot(settings, results)
        _, _, _, pept_neighbors = get_prot_pept_map(args, neighbors=True)
#         if settings.has_option('misc', 'aa_mass'):
#             aa_mass = settings.get('misc', 'aa_mass')
#         else:
#             aa_mass = get_aa_mass(settings)
#         vmods = set()
#         variablemods = settings.get('modifications', 'variable')
#         if variablemods:
#             for k, v in variablemods.items():
#                 for aa in v:
#                     vmods.add(k + aa)
#                     vmods.add(aa + k)

#         leg = {}
#         if settings.has_option('misc', 'legend'):
#             leg = settings.get('misc', 'legend')
#         if settings.has_option('misc', 'plegend'):
#             leg.update(settings.get('misc', 'plegend'))

#         ntermcleavage = settings.getfloat('modifications', 'protein nterm cleavage')
#         ctermcleavage = settings.getfloat('modifications', 'protein cterm cleavage')

        df1['idx'] = df1.index
        if 'best_spectrum_id' not in df1.columns:
            df1['best_spectrum_id'] = df1['idx']
        for idx, neutral_mass, charge_state, RT, comp_voltage, sequence, matched_ions, md, hyperscore_score, hyperscore3_score, sf_score, ms1_intensity, b_sum, y_sum, missed_cleavages, best_spectrum_id in df1[['idx', 'nmasses', 'ch', 'rt', 'im', 'seqs', 'matched_ions', 'md', 'hyperscore', 'hyperscore3', 'sf', 'Is', 'b_count', 'y_count', 'mc', 'best_spectrum_id']].values:
        # for idx, result in enumerate(df1):
            if 1:
                tmp = etree.Element('spectrum_query')
                # spectrum = result['spectrum']
                tmp.set('spectrum', str(best_spectrum_id))
                tmp.set('spectrumNativeID', str(best_spectrum_id))
                # tmp.set('spectrum', get_title(spectrum))
                # tmp.set('spectrumNativeID', get_title(spectrum))
                tmp.set('start_scan', str(idx))  # ???
                tmp.set('end_scan', str(idx))  # ???
                tmp.set('index', str(idx))  # ???

#                 neutral_mass, charge_state, RT, comp_voltage = get_info(spectrum, result, settings, aa_mass)
                tmp.set('precursor_neutral_mass', str(neutral_mass))
                tmp.set('assumed_charge', str(int(charge_state)))
                if RT:
                    tmp.set('retention_time_sec', str(RT*60))
                if comp_voltage:
                    tmp.set('compensation_voltage', str(comp_voltage))

                tmp2 = etree.Element('search_result')
#                 result['candidates'] = result['candidates'][:len(result['e-values'])]

                flag = 1
                i = 0
                if 1:
#                 for i, candidate in enumerate(result['candidates']):
#                     match = candidate[4]['match']
#                     if match is None:
#                         break
                    tmp3 = etree.Element('search_hit')
                    tmp3.set('hit_rank', str(i + 1))
#                     mod_sequence = normalize_mods(str(candidate[1]), settings)
                    # sequence = re.sub(r'[^A-Z]', '', mod_sequence)
                    if sequence not in pept_prot:
                        flag = 0
                        logger.error('Unaccounted sequence! %s (%s)', sequence, mod_sequence)
                        break
                    else:
                        tmp3.set('peptide', sequence)

                        proteins = list(pept_prot[sequence])
                        tmp3.set('protein', proteins[0])
                        protein_descr = ''


                        # tmp3.set('protein', prots[proteins[0]].split(' ', 1)[0] + (('_' + candidate[7]) if snp else ''))
                        # try:
                        #     protein_descr = prots[proteins[0]].split(' ', 1)[1]
                        # except:
                        #     protein_descr = ''

                        neighbors = pept_neighbors.get(sequence, {}).get(proteins[0], ('-', '-'))

                        tmp3.set('peptide_prev_aa', neighbors[0])
                        tmp3.set('peptide_next_aa', neighbors[1])
                        tmp3.set('protein_descr', protein_descr)

                        num_tot_proteins = len(proteins)
                        tmp3.set('num_tot_proteins', str(num_tot_proteins))
                        tmp3.set('num_matched_ions', str(matched_ions))
                        tmp3.set('tot_num_ions', str((len(sequence) - 1) * 2))
                        neutral_mass_theor = neutral_mass / (1 - md / 1e6)
#                         neutral_mass_theor = custom_mass(str(candidate[1]), aa_mass=aa_mass, nterm_mass=nterm_mass, cterm_mass=cterm_mass)
#                         # neutral_mass_theor = cmass.fast_mass(sequence, aa_mass=aa_mass)
                        tmp3.set('calc_neutral_pep_mass', str(neutral_mass_theor))
                        tmp3.set('massdiff', str(neutral_mass * md / 1e6))
                        tmp3.set('num_tol_term', str(2))
                        # tmp3.set('num_tol_term', str(pept_ntts.get(sequence, {}).get(proteins[0], '?')))
                        tmp3.set('num_missed_cleavages', str(missed_cleavages))
                        # tmp3.set('is_rejected', '0')  # ???

                        if num_tot_proteins > 1 and (not snp or 'wild' not in prots[proteins[0]].split(' ', 1)[0]):
                            for prot in proteins[1:]:
                                tmp4 = etree.Element('alternative_protein')
                                tmp4.set('protein', prot)
                                # tmp4.set('protein', prots[prot].split(' ', 1)[0] + (('_' + candidate[7]) if snp else ''))
                                try:
                                    protein_descr = prots[prot].split(' ', 1)[1]
                                except:
                                    protein_descr = ''
                                tmp4.set('protein_descr', protein_descr)
                                neighbors = pept_neighbors.get(sequence, {}).get(prot, ('-', '-'))
                                tmp4.set('peptide_prev_aa', neighbors[0])
                                tmp4.set('peptide_next_aa', neighbors[1])
                                # tmp4.set('num_tol_term', str(pept_ntts.get(sequence, {}).get(prot, '?')))
                                tmp4.set('num_tol_term', str(2))
                                tmp3.append(copy(tmp4))

#                         labels = parser.std_labels + [la[:-1] if la[-1] == '[' else '-' + la[:-2] if la[-1] == ']' else la for la in leg if len(la) > 1]
# #                       logger.debug('Known labels: %s', labels)
#                         try:
#                             aalist = parser.parse(mod_sequence, labels=labels)
#                         except Exception as e:
#                             logger.debug('Problematic sequence: %s\n%s', mod_sequence, e)
#                             aalist = [a[::-1] for a in parser.parse(mod_sequence[::-1], labels=labels)][::-1]
                        tmp4 = etree.Element('modification_info')
                        ntermmod = 0

#                         # if nterm_fixed:
#                         #     tmp4.set('mod_nterm_mass', str(nterm_fixed))
#                         # if cterm_fixed:
#                         #     tmp4.set('mod_cterm_mass', str(cterm_fixed))

                        # for idx, aminoacid in enumerate(aalist):
                        #     if aminoacid in fmods or aminoacid in vmods:
                        #         if aminoacid.endswith('-') and idx == 0:
                        #             ntermmod = 1
                        #             tmp4.set('mod_nterm_mass', str(str(aa_mass.get(aminoacid) + ntermcleavage)))
                        #         elif aminoacid.startswith('-') and idx == len(aalist) - 1:
                        #             tmp4.set('mod_cterm_mass', str(aa_mass.get(aminoacid) + ctermcleavage))
                        #         else:
                        #             tmp5 = etree.Element('mod_aminoacid_mass')
                        #             tmp5.set('position', str(idx + 1 - ntermmod))
                        #             tmp5.set('mass', str(aa_mass.get(aminoacid)))
                        #             tmp4.append(copy(tmp5))
                        # tmp3.append(copy(tmp4))

                        # if 'RNHS' in candidate[4]:

                        tmp4 = etree.Element('search_score')
                        tmp4.set('name', 'expect')
                        tmp4.set('value', str(1./(hyperscore3_score+1e-3)))
                        tmp3.append(copy(tmp4))


                        tmp4 = etree.Element('search_score')
                        tmp4.set('name', 'hyperscore')
                        tmp4.set('value', str(hyperscore_score))
                        tmp3.append(copy(tmp4))


                        tmp4 = etree.Element('search_score')
                        tmp4.set('name', 'hyperscore3')
                        tmp4.set('value', str(hyperscore3_score))
                        tmp3.append(copy(tmp4))


                        tmp4 = etree.Element('search_score')
                        tmp4.set('name', 'sf')
                        tmp4.set('value', str(sf_score))
                        tmp3.append(copy(tmp4))

                            # tmp4 = etree.Element('search_score')
                            # tmp4.set('name', 'expect')
                            # tmp4.set('value', str(1./candidate[4]['RNHS']))
                            # tmp3.append(copy(tmp4))

#                         else:
#                             tmp4 = etree.Element('search_score')
#                             tmp4.set('name', 'hyperscore')
#                             tmp4.set('value', str(candidate[0]))
#                             tmp3.append(copy(tmp4))

#                             tmp4 = etree.Element('search_score')
#                             tmp4.set('name', 'expect')
#                             tmp4.set('value', str(result['e-values'][i]))
#                             tmp3.append(copy(tmp4))

#                         tmp4 = etree.Element('search_score')
#                         tmp4.set('name', 'sumI')
#                         tmp4.set('value', str(candidate[5]))
#                         tmp3.append(copy(tmp4))

#                         tmp4 = etree.Element('search_score')
#                         tmp4.set('name', 'fragmentMT')
#                         tmp4.set('value', str(candidate[6]))
#                         tmp3.append(copy(tmp4))

#                         tmp4 = etree.Element('search_score')
#                         tmp4.set('name', 'nextscore_std')
#                         tmp4.set('value', str(candidate[8]))
#                         tmp3.append(copy(tmp4))

#                         if 'params' in spectrum:
#                             if 'isowidthdiff' in spectrum['params']:
#                                 tmp4 = etree.Element('search_score')
#                                 tmp4.set('name', 'ISOWIDTHDIFF')
#                                 tmp4.set('value', str(spectrum['params'].get('isowidthdiff', 0)))
#                                 tmp3.append(copy(tmp4))

#                             if 'rtwidth' in spectrum['params']:
#                                 tmp4 = etree.Element('search_score')
#                                 tmp4.set('name', 'RTwidth')
#                                 tmp4.set('value', str(spectrum['params'].get('rtwidth', 0)))
#                                 tmp3.append(copy(tmp4))

                        tmp4 = etree.Element('search_score')
                        tmp4.set('name', 'MS1Intensity')
                        tmp4.set('value', str(ms1_intensity))
                        tmp3.append(copy(tmp4))

#                             if 'pif' in spectrum['params']:
#                                 tmp4 = etree.Element('search_score')
#                                 tmp4.set('name', 'PIF')
#                                 tmp4.set('value', str(spectrum['params'].get('pif', -3)))
#                                 tmp3.append(copy(tmp4))

#                             if 'sulfur' in spectrum['params']:
#                                 tmp4 = etree.Element('search_score')
#                                 tmp4.set('name', 'sulfur')
#                                 tmp4.set('value', str(spectrum['params'].get('sulfur', -1)))
#                                 tmp3.append(copy(tmp4))

                            # tmp4 = etree.Element('search_score')
                            # tmp4.set('name', 'ionmobility')
                            # tmp4.set('value', str(im_value))
                            # tmp3.append(copy(tmp4))

#                         # if tags:
#                         #     for tag_label in tags.keys():
#                         #         tmp4 = etree.Element('search_score')
#                         #         tmp4.set('name', 'tag_' + tag_label)
#                         #         tmp4.set('value', str(spectrum.get(tag_label, 0)))
#                         #         tmp3.append(copy(tmp4))


                        tmp4 = etree.Element('search_score')
                        tmp4.set('name', 'matched_b1_ions')
                        tmp4.set('value', str(b_sum))
                        tmp3.append(copy(tmp4))

                        tmp4 = etree.Element('search_score')
                        tmp4.set('name', 'matched_y1_ions')
                        tmp4.set('value', str(y_sum))
                        tmp3.append(copy(tmp4))

                        tmp2.append(copy(tmp3))
                if flag:
                    tmp.append(copy(tmp2))
                    child1.append(copy(tmp))

        s = etree.tostring(root, pretty_print=True)
        output.write(s)




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

def iterate_spectra(fname, min_ch, max_ch, min_isotopes, min_scans, nproc, check_unique=True, systematic_mass_shift=0):
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
            'ivf': 5,
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
    if 'im' not in df_features.columns:
        df_features['im'] = 0

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

    if systematic_mass_shift:
        df_features['massCalib'] = df_features['massCalib'].apply(lambda x: x * (1 - 1e-6 * systematic_mass_shift))
        df_features['mz'] = df_features['mz'].apply(lambda x: x * (1 - 1e-6 * systematic_mass_shift))

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

try:
    from identipy import cparser
except ImportError:
    from identipy import customparser as cparser


def get_peptides2(prot_seq, enzyme, mc, minlen, maxlen, semitryptic=False):
    peptides = cparser._cleave(prot_seq, enzyme, mc)
    for pep, startposition in peptides:
        plen = len(pep)
        if minlen <= plen <= maxlen:
            if not semitryptic:
                yield pep, startposition, plen
            else:
                for i in range(plen-minlen+1):
                    yield pep[i:], startposition + i, plen - i
                for i in range(1, plen-minlen+1, 1):
                    yield pep[:-i], startposition, plen - i


seen_target = set()
seen_decoy = set()
def prot_peptides2(prot_seq, enzyme, mc, minlen, maxlen, is_decoy, dont_use_seen_peptides=False):


    dont_use_fast_valid = parser.fast_valid(prot_seq)
    # peptides = parser.cleave(prot_seq, enzyme, mc)
    # for pep in peptides:
    #     plen = len(pep)
    for pep, startposition, plen in get_peptides2(prot_seq, enzyme, mc, minlen, maxlen):
        if minlen <= plen <= maxlen:
            forms = []
            if dont_use_fast_valid or pep in seen_target or pep in seen_decoy or parser.fast_valid(pep):
                if plen <= maxlen:
                    forms.append(pep)
            for f in forms:
                if dont_use_seen_peptides:
                    yield (f, startposition)
                else:
                    if f not in seen_target and f not in seen_decoy:
                        if is_decoy:
                            seen_decoy.add(f)
                        else:
                            seen_target.add(f)
                        yield (f, startposition)

def get_prot_pept_map(args, neighbors=False):
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

    pept_neighbors = {}



    for desc, prot in prot_gen(args):
        dbinfo = desc.split(' ')[0]

        if neighbors:

            for pep, startposition in prot_peptides2(prot, enzyme, mc, minlen, maxlen, desc.startswith(prefix), dont_use_seen_peptides=True):
                pept_prot.setdefault(pep, set()).add(dbinfo)

                pept_neighbors.setdefault(pep, {})
                pept_neighbors[pep][dbinfo] = (prot[startposition - 1] if startposition != 0 else '-', prot[startposition + len(pep)] if startposition + len(pep) < len(prot) else '-')

                protsN.setdefault(dbinfo, set()).add(pep)
        else:
            for pep in prot_peptides(prot, enzyme, mc, minlen, maxlen, desc.startswith(prefix), dont_use_seen_peptides=True):
                pept_prot.setdefault(pep, set()).add(dbinfo)
                protsN.setdefault(dbinfo, set()).add(pep)
    for k, v in protsN.items():
        if k.startswith(prefix):
            decoy_prot_count += 1
            decoy_peps.update(v)
        else:
            target_prot_count += 1
            target_peps.update(v)

        protsN[k] = len(v)


    if neighbors:
        del decoy_peps
        del target_peps
        return protsN, pept_prot, 0, pept_neighbors

    else:
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



def get_prot_pept_map_semi(args, pept_prot):


    prefix = args['prefix']
    minlen = args['lmin']
    maxlen = args['lmax']

    pept_prot2 = dict()
    pept_prot3 = defaultdict(set)
    protsN2 = dict()
    protsN3 = dict()

    # split_aa = {'F', 'Y', 'W', 'M'}

    for pep, prots in pept_prot.items():
        # for idx, pa in enumerate(pep):
        for idx in range(len(pep)-1):
            # if pa in split_aa:
            if 1:
                pl = pep[:idx+1]
                pr = pep[idx+1:]
                if minlen <= len(pl) <= maxlen:
                    pept_prot2.setdefault(pl, set()).add(pep)
                    pept_prot3[pl].update(prots)
                    protsN2.setdefault(pep, set()).add(pl)
                    for prot in prots:
                        protsN3.setdefault(prot, set()).add(pl)

                if minlen <= len(pr) <= maxlen:
                    pept_prot2.setdefault(pr, set()).add(pep)
                    pept_prot3[pr].update(prots)
                    protsN2.setdefault(pep, set()).add(pr)
                    for prot in prots:
                        protsN3.setdefault(prot, set()).add(pr)
    for k, v in protsN3.items():
        protsN3[k] = len(v)


    return protsN2, pept_prot2, pept_prot3, protsN3


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

def calc_sf_all(v, n, p, prev_best_score=False, p_array=False):
    if p_array:
        p = p.clip(min=0.0001)
    else:
        p = max(0.0001, p)
    sf_values = -np.log10(binom.sf(v-1, n, p))
    sf_values[np.isnan(sf_values)] = 0
    sf_values[np.isinf(sf_values)] = (prev_best_score if prev_best_score is not False else max(sf_values[~np.isinf(sf_values)]) * 2)
    return sf_values
