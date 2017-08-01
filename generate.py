from sys import argv
from os import path
from pyteomics import fasta, parser, mass, achrom, electrochem
from copy import copy
from random import gauss, random, choice
import re
import pickle
import csv
from collections import Counter
import numpy as np
import scipy.stats
import multiprocessing
from time import sleep
from collections import defaultdict
import argparse

pept_prot = defaultdict(list)

def get_enzyme(enzyme):
    if enzyme in parser.expasy_rules:
        return parser.expasy_rules.get(enzyme)
    else:
        try:
            enzyme = convert_tandem_cleave_rule_to_regexp(enzyme)
            return enzyme
        except:
            return enzyme

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

class ScanWriter():
    def __init__(self, outfile, protease, last_scan_time=None):
        # self.ms1_scan_time = settings.getfloat('ms1', 'ms1_scan_time')
        self.protease = protease
        # self.peptide_elution_time = settings.getfloat('ms1', 'peptide_elution_time')
        # self.output_folder = settings.get('options', 'output_path')
        # if not self.output_folder:
        #     self.output_folder = ''
        self.outfile = outfile
        self.cur_time = 0
        self.last_scan_time = last_scan_time
        # self.norm = scipy.stats.norm(0, self.peptide_elution_time / 2 / 3)
        # self.norm_dict = dict()

    def custom_pdf(self, v):
        v = round(v, 2)
        if v not in self.norm_dict:
            self.norm_dict[v] = self.norm.pdf(v)
        return self.norm_dict[v]

    def write(self, peptides, q_out, nprocs, flagname='', peptide_efficiency=1, mc_efficiency=1):
        filename = '%s%speaks' % (flagname + self.protease.replace(' ', ''), path.extsep)
        # scans = np.arange(0, self.last_scan_time, self.ms1_scan_time / 60000)
        with open(self.outfile + '_temp', 'w') as output:
            peakwriter = csv.writer(output, delimiter='\t')
            # peakwriter.writerow(('mz', 'RT', 'Intensity', 'charge', 'sequence', 'noise', 'mc'))
            peakwriter.writerow(('massCalib', 'rtApex', 'charge', 'sequence', 'noise', 'mc', 'nIsotopes'))
            while q_out.qsize() != nprocs or len(peptides) != 0:
                if len(peptides) == 0:
                    sleep(1)
                else:
                    k, v = peptides.popitem()
                    # left_idx = scans.searchsorted(v['RT'] - self.peptide_elution_time / 2, side='left')
                    # right_idx = scans.searchsorted(v['RT'] + self.peptide_elution_time / 2, side='right')
                    # for scan_time in scans[left_idx:right_idx]:
                    peakwriter.writerow((round(v['mz'] * int(v['ch']) - int(v['ch']) * 1.0073, 5), v['RT'], int(v['ch']), k, v['noise'], v['mc'], 5))
        output_full = open(self.outfile + '_temp', 'r')
        output_short = open(self.outfile, 'w')
        sequences = defaultdict(set)
        mc_todel = dict()
        sequences_tostay = set()
        checked = set()
        output_full.readline()
        for peak in output_full:
            t = peak.strip().split('\t')
            if t[4] == 'no':
                sequences[int(t[5])].add(t[3])
        output_full.close()
        output_full = open(self.outfile + '_temp', 'r')
        mc_todel[0] = peptide_efficiency * mc_efficiency[0]
        for i in [1, 2]:
            mc_todel[i] = float(len(sequences[0])) * mc_todel[0] * mc_efficiency[i] / len(sequences[i])
        output_short.write(output_full.readline())
        for peak in output_full:
            t = peak.strip().split('\t')
            seq = t[3]
            if seq not in checked:
                checked.add(seq)
                mc = int(t[5])
                if mc_todel.get(mc, 0) >= random():
                    sequences_tostay.add(seq)
        output_full.close()
        output_full = open(self.outfile + '_temp', 'r')
        output_full.readline()
        for peak in output_full:
            t = peak.strip().split('\t')
            seq = t[3]
            if seq in sequences_tostay:
                output_short.write(peak)
        output_full.close()
        output_short.close()

        print len(sequences[0]), len(sequences[1]), len(sequences[2])
        print mc_todel


class Peptides():
    def __init__(self):
        self.sequences = np.array([])
        self.mzs = np.array([])
        self.RTs = np.array([])
        self.intensities = np.array([])
        self.charges = np.array([])
        self.noise_labels = np.array([])

    def sort_index(self, idx):
        self.sequences = self.sequences[idx]
        self.mzs = self.mzs[idx]
        self.RTs = self.RTs[idx]
        self.intensities = self.intensities[idx]
        self.charges = self.charges[idx]
        self.noise_labels = self.noise_labels[idx]

    def extend(self, peptides):
        self.sequences = np.append(self.sequences, peptides.keys())
        self.mzs = np.append(self.mzs, [v['mz'] for v in peptides.itervalues()])
        self.RTs = np.append(self.RTs, [v['RT'] for v in peptides.itervalues()])
        self.intensities = np.append(self.intensities, [v['intensity'] for v in peptides.itervalues()])
        self.charges = np.append(self.charges, [v['ch'] for v in peptides.itervalues()])
        self.noise_labels = np.append(self.noise_labels, [v['noise'] for v in peptides.itervalues()])


def get_exp_mass(theor_mass, mz_tol):
    return gauss(theor_mass, mz_tol * 1e-6 * theor_mass)


def get_exp_RT(theor_rt, rt_tol):
    return gauss(theor_rt, rt_tol)

def get_peptide_info(peptide, aa_per_charge, charge_min, charge_max, mz_min, mz_max, RT_min, RT_max, custom_aa_mass,
                     RC, RT_tol, protease, peptide_efficiency, mc_efficiency, mz_tol, noise):
    ch = round(float(len(peptide)) / aa_per_charge, 0)
    if charge_min <= ch <= charge_max:
        peptide_mass = mass.fast_mass(peptide, charge=ch, aa_mass=custom_aa_mass)
        if mz_min <= peptide_mass <= mz_max:
            peptide_RT = achrom.calculate_RT(peptide, RC)
            if RT_min <= peptide_RT <= RT_max:
                peptide_RT_exp = get_exp_RT(peptide_RT, RT_tol)
                mc = parser.num_sites(peptide, protease)#get_number_of_mc(peptide, protease)
                if 1:#noise == 'yes' or (mc in mc_efficiency and peptide_efficiency * mc_efficiency[mc] > random()):
                    peptide_mass_exp = get_exp_mass(peptide_mass, mz_tol)
                    # if pI_tol != -1:
                    #     pI = get_exp_pI(electrochem.pI(peptide), pI_tol)
                    # else:
                    #     pI = 0
                    intensity = 100000  # TODO
                    peptide_dict = {'mz': peptide_mass_exp, 'ch': ch,
                                         'mc': mc, 'RT': peptide_RT_exp,
                                         'intensity': intensity, 'noise': noise}
#                    print 'ok'
                    return peptide_dict
                else:
#                    print 'mc'
                    return None
            else:
#                print 'RT'
                return None
        else:
#            print 'mz'
            return None
    else:
#        print 'ch'
        return None


def get_stats(peptides):
    stats = Counter()
    for peptide in peptides:
        stats[len(peptide)] += 1
    return stats


def generate_wrong_peaks(peptide, peptides, dropped_peptides, aa_per_charge, charge_min, charge_max, mz_min, mz_max, RT_min, RT_max, custom_aa_mass,
                     RC, RT_tol, protease, peptide_efficiency, mc_efficiency, mz_tol, noise_peptides_per_true_peptide):
    aminoacids = mass.std_aa_mass.keys()
    l = len(peptide)
    peptides_wrong = dict()
#    for k, v in stats.iteritems():
    i = noise_peptides_per_true_peptide
    cycle_flag = 0
    while i > 0:
        peptide = ''
        for _ in range(l):
            peptide += choice(aminoacids)
        if cycle_flag > 10000:  # TODO
            print 'skipping generation of spectra for %d length' % (l, )
            break
        if peptide not in dropped_peptides and peptide not in peptides:
            peptide_dict = get_peptide_info(peptide, aa_per_charge, charge_min, charge_max, mz_min,
                                        mz_max, RT_min, RT_max, custom_aa_mass,
                                        RC, RT_tol, protease, peptide_efficiency, mc_efficiency, mz_tol, 'yes')
            if peptide_dict:
                peptides_wrong[peptide] = copy(peptide_dict)
                i -= 1
            else:
                dropped_peptides[peptide] = None
        cycle_flag += 1

    return peptides_wrong

banned_aminoacids = ['B', 'X', 'J', 'Z', 'U', 'O']

arparser = argparse.ArgumentParser(
    description='Search proteins using LC-MS spectra',
    epilog='''

Example usage
-------------
$ generate.py output.tsv -d human.fasta
-------------
''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

arparser.add_argument('file', help='output file with in-silico peptide features')
arparser.add_argument('-db', help='path to protein fasta file', required=True)
arparser.add_argument('-ptol', help='in-silico precursor mass sigma in ppm', default=3.3, type=float)
arparser.add_argument('-rttol', help='in-silico retention time sigma in min', default=3.0, type=float)
arparser.add_argument('-gradient', help='gradient time in min', default=60.0, type=float)
arparser.add_argument('-mzmin', help='min mz of peptides', default=300, type=int)
arparser.add_argument('-mzmax', help='max mz of peptides', default=1500, type=int)
arparser.add_argument('-e', help='cleavage rule in quotes!. X!Tandem style for cleavage rules', default='[RK]|{P}')
arparser.add_argument('-enzeff', help='Enzyme efficiency in %% for 0/1/2 missed cleavages', default='100/10/0', type=str)
arparser.add_argument('-cmin', help='min precursor charge', default=2, type=int)
arparser.add_argument('-cmax', help='max precursor charge', default=5, type=int)
arparser.add_argument('-aach', help='amino acids per charge', default=6, type=int)
arparser.add_argument('-noise', help='noise peaks per true one', default=10, type=int)
arparser.add_argument('-fmods', help='fixed modifications. in mass1@aminoacid1,mass2@aminoacid2 format', default='57.021464@C')
args = vars(arparser.parse_args())

mz_min = args['mzmin']
# mz_min = settings.getfloat('ms1', 'mz_min')
mz_max = args['mzmax']
# mz_max = settings.getfloat('ms1', 'mz_max')
RT_min = 0
RT_max = args['gradient']
# charge_min = settings.getfloat('ms1', 'charge_min')
charge_min = args['cmin']
# charge_max = settings.getfloat('ms1', 'charge_max')
charge_max = args['cmax']

fmods = args['fmods']
custom_aa_mass = copy(mass.std_aa_mass)
if fmods:
    for mod in fmods.split(','):
        m, aa = mod.split('@')
        custom_aa_mass[aa] += float(m)

# custom_aa_mass = copy(mass.std_aa_mass)
# fmods = settings.get('modifications', 'fixed')
# if fmods:
#     for mod in re.split(r'[,;]\s*', fmods):
#         m, aa = parser._split_label(mod)
#         custom_aa_mass[aa] += settings.getfloat('modifications', m)

# RC = settings.get('options', 'path_to_RC')
# if not RC:
RC = {'aa': {'-OH': 0.0,
  'A': 5.7483030141800882,
  'C': 3.8412787644679742,
  'C*': 3.8412787644679742,
  'D': 5.2021408783797467,
  'E': 5.2591061803802877,
  'F': 23.706146443823165,
  'G': 2.6896156182856945,
  'H': -7.4839060978249607,
  'H-': 0.0,
  'I': 18.345336546344498,
  'K': -5.4598784163894685,
  'L': 20.160556727493173,
  'M': 14.670500432066241,
  'N': 2.3154434919886659,
  'O': -2.6046335316136755,
  'P': 6.2944053433804692,
  'Q': 2.7862405120985883,
  'R': -4.1894596240485038,
  'S': 3.1739734369088382,
  'T': 4.4448646506590901,
  'U': 3.8412787644679742,
  'V': 12.108599919296113,
  'W': 25.684657074865044,
  'Y': 11.867335791510143,
  'c-': 19.303536395546583,
  'camC': 3.0998874568392214},
 'const': -18.649446141257613,
 'lcp': -0.20000000000000018}
gradkoef = args['gradient'] / 60
for k in RC['aa'].keys():
    RC['aa'][k] *= gradkoef
RC['const'] *= gradkoef
# else:
#     RC = pickle.load(open(RC, 'r'))

aa_per_charge = args['aach']
# aa_per_charge = settings.getfloat('ms1', 'aa_per_charge')
mz_tol = args['ptol']
RT_tol = args['rttol']
# mz_tol = settings.getfloat('tolerances', 'mz_tol')
# RT_tol = settings.getfloat('tolerances', 'RT_tol')
# pI_tol = settings.getfloat('tolerances', 'pI_tol')

# peptide_efficiency = 1.0#settings.getfloat('ms1', 'peptide_efficiency') / 100
mc0, mc1, mc2 = args['enzeff'].split('/')
mc_efficiency = {0: float(mc0) / 100,
                 1: float(mc1) / 100,
                 2: float(mc2) / 100}
noise_peptides_per_true_peptide = args['noise']
# noise_peptides_per_true_peptide = settings.getfloat('ms1', 'noise_peptides_per_true_peptide')

i = 1
protease = get_enzyme(args['e'])
# protease = settings.get('proteases', 'protease%s' % (i, ))
# if protease:
scan_writer = ScanWriter(args['file'], protease, last_scan_time=200)
manager = multiprocessing.Manager()
peptides = manager.dict()
peptides_valid = manager.dict()
dropped_peptides = manager.dict()
path_to_proteins = args['db']#settings.get('options', 'path_to_proteins')

procs = []
nprocs = 12
q = multiprocessing.Queue()
q_out = multiprocessing.Queue()

def get_peptides(q, q_out, peptides, dropped_peptides, aa_per_charge, charge_min, charge_max, mz_min, mz_max, RT_min, RT_max, custom_aa_mass,
             RC, RT_tol, protease, peptide_efficiency, mc_efficiency, mz_tol, noise_peptides_per_true_peptide, pept_prot):
    for peptide, peptide_efficiency in iter(q.get, None):
        if peptide not in peptides and peptide not in dropped_peptides:
            peptide_dict = get_peptide_info(peptide, aa_per_charge, charge_min, charge_max,
                                            mz_min, mz_max, RT_min, RT_max, custom_aa_mass,
                                            RC, RT_tol, protease, peptide_efficiency, mc_efficiency,
                                            mz_tol, noise='no')
            if peptide_dict:
                peptides[peptide] = copy(peptide_dict)
                wrong_peptide_dict = generate_wrong_peaks(peptide, peptides, dropped_peptides, aa_per_charge, charge_min, charge_max, mz_min, mz_max, RT_min, RT_max, custom_aa_mass,
                     RC, RT_tol, protease, peptide_efficiency, mc_efficiency, mz_tol, noise_peptides_per_true_peptide)
                for k, v in wrong_peptide_dict.iteritems():
                    peptides[k] = copy(v)
    q_out.put(None)

if path.splitext(path_to_proteins)[-1] == '.fasta':
    flagname = ''
    for prot in fasta.read(path_to_proteins):
        peptide_efficiency = np.random.normal(loc=70.0, scale=10.0, size=None) / 100
        for peptide in parser.cleave(prot[1], protease, missed_cleavages=2):
            if all(z not in peptide for z in ['B', 'X', 'J', 'Z', 'U', 'O']):
                q.put((peptide, peptide_efficiency))
                # q.put(peptide)
                pept_prot[peptide].append(prot[0].split('|')[1])
# elif path.splitext(path_to_proteins)[-1] == '.peaks':
#     flagname = path.basename(path.splitext(path_to_proteins)[0]) + '_then_'
#     added_sequences = set()
#     for line in open(path_to_proteins, 'r'):
#         if line.split('\t')[5].strip() == 'no':
#             sequence = line.split('\t')[4]
#             if sequence not in added_sequences:
#                 added_sequences.add(sequence)
#                 for peptide in parser.cleave(sequence, protease, missed_cleavages=2):
#                     q.put(peptide)
#                     pept_prot[peptide].append(sequence)

for _ in range(nprocs):
    p = multiprocessing.Process(target=get_peptides, args=(q, q_out, peptides, dropped_peptides, aa_per_charge, charge_min, charge_max, mz_min, mz_max, RT_min, RT_max, custom_aa_mass,
             RC, RT_tol, protease, peptide_efficiency, mc_efficiency, mz_tol, noise_peptides_per_true_peptide, pept_prot))
    procs.append(p)
    p.start()

for _ in procs:
    q.put(None)

scan_writer.write(peptides, q_out, nprocs, flagname, peptide_efficiency, mc_efficiency)
for p in procs:
    p.terminate()
i += 1

