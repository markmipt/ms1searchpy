from pyteomics import fasta, parser
from multiprocessing import Queue, Process, cpu_count
import os
import csv
import subprocess


def iterate_spectra(fname, min_ch, max_ch, min_isotopes, min_scans):
    if os.path.splitext(fname)[-1].lower() == '.mzml':
        subprocess.call(['java', '-Djava.awt.headless=true', '-jar', os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Dinosaur/Dinosaur-1.1.3.free.jar'), '--advParams=Dinosaur/adv.txt', '--concurrency=12', fname])
        fname = os.path.splitext(fname)[0] + '.features.tsv'
    with open(fname, 'rb') as infile:
        csvreader = csv.reader(infile, delimiter='\t')
        header = csvreader.next()
        mass_ind = header.index('massCalib')
        RT_ind = header.index('rtApex')
        ch_ind = header.index('charge')
        nIsotopes_ind = header.index('nIsotopes')
        Int_ind = header.index('intensityApex')
        nScans_ind = header.index('nScans')
        idx = 0
        for z in csvreader:
            nm = float(z[mass_ind])
            RT = float(z[RT_ind])
            ch = float(z[ch_ind])
            nIsotopes = float(z[nIsotopes_ind])
            nScans = int(z[nScans_ind])
            I = float(z[Int_ind])
            idx += 1
            if nIsotopes >= min_isotopes and min_ch <= ch <= max_ch and min_scans <= nScans:
                yield nm, RT, ch, idx, I, nScans


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

seen_target = set()
seen_decoy = set()
def prot_peptides(prot_seq, enzyme, mc, minlen, maxlen, is_decoy, dont_use_seen_peptides=False):

    dont_use_fast_valid = parser.fast_valid(prot_seq)
    peptides = parser.cleave(prot_seq, enzyme, mc)
    for pep in peptides:
        plen = len(pep)
        if minlen <= plen <= maxlen + 2:
            forms = []
            if dont_use_fast_valid or pep in seen_target or pep in seen_decoy or parser.fast_valid(pep):
                if plen <= maxlen:
                    forms.append(pep)
                if prot_seq[0] == 'M' and prot_seq.startswith(pep):
                    if minlen <= plen - 1 <= maxlen:
                        forms.append(pep[1:])
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