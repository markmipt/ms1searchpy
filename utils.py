from ConfigParser import RawConfigParser
import sys
sys.path.append('/usr/share/mmass')
from sys import argv
from mspy import parser_mzml
from time import time
import os
from pyteomics import mass, electrochem as ec, auxiliary as aux, fasta, mgf, mzml, parser
from multiprocessing import Queue, Process, cpu_count

class CustomRawConfigParser(RawConfigParser):
    def get(self, section, option):
        val = RawConfigParser.get(self, section, option)
        if isinstance(val, basestring):
            return val[::-1].split('|', 1)[-1][::-1]
        return val

    def get_choices(self, section, option):
        val = RawConfigParser.get(self, section, option)
        if isinstance(val, basestring) and len(val.split('|')) > 1:
            return val[::-1].split('|', 1)[0][::-1]
        else:
            return ''

def iterate_spectra(fname, min_ch, max_ch, min_i):
    MZML_data = parser_mzml.parseMZML(fname)
    print 'parsing was started'
    MZML_data.load()
    print 'parsing was finished'
    procs = []
    nprocs = 12
    q = Queue()
    q_output = Queue()

    def deiso(q, q_output, min_ch, max_ch):
        for sc in iter(q.get, None):
            peaks = sc.peaklist
            if len(peaks) and peaks[-1].mz - peaks[0].mz >= 600:
                RT = sc.retentionTime / 60
                peaks = sc.peaklist
                peaks.deisotope(maxCharge=max_ch, mzTolerance=0.01, intTolerance=0.5)
                peaks.remisotopes()
                for peak in peaks:
                    if peak.charge and min_ch <= peak.charge <= max_ch and peak.intensity >= min_i:
                        neutral = peak.mz * peak.charge - 1.0073 * peak.charge
                        q_output.put((neutral, RT, peak.intensity, peak.charge))
        q_output.put(None)

    for i in range(nprocs):
        p = Process(target=deiso, args=(q, q_output, min_ch, max_ch))
        procs.append(p)
        p.start()


    idx = 1
    ms1_idx = 1

    while 1:
        sc = MZML_data.scan(idx)
        if sc == False:
            print 'total number of spectra = %d\ntotal number of ms1 spectra = %d' % (idx - 1, ms1_idx - 1)
            break
        if sc.msLevel == 1:
            q.put(sc)
            ms1_idx += 1
        idx += 1

    for p in procs:
        q.put(None)

    del MZML_data

    while nprocs > 0:
        for rec in iter(q_output.get, None):
            yield rec
        nprocs -= 1

    for p in procs:
        p.terminate()


def iterate_spectra_BU(fname, min_ch, max_ch):
    # filename = os.path.splitext(fname)[0]
    MZML_data = parser_mzml.parseMZML(fname)
    print 'parsing was started'
    MZML_data.load()
    print 'parsing was finished'
    # output_csv = open('%s_peaks.csv' % (filename,), 'w')
    idx = 1
    ms1_idx = 1

    while 1:
        sc = MZML_data.scan(idx)
        if sc == False:
            print 'total number of spectra = %d\ntotal number of ms1 spectra = %d' % (idx - 1, ms1_idx - 1)
            break
        if sc.msLevel == 1:
            peaks = sc.peaklist
            if len(peaks) and peaks[-1].mz - peaks[0].mz >= 600:
                if ms1_idx % 1000 == 0:
                    print '%d of ms1 spectra were processed' % (ms1_idx, )

                RT = sc.retentionTime / 60
                peaks = sc.peaklist
                peaks.deisotope(maxCharge=max_ch, mzTolerance=0.01, intTolerance=0.5)
                peaks.remisotopes()
                for peak in peaks:
                    if peak.charge and min_ch <= peak.charge <= max_ch:
                        neutral = peak.mz * peak.charge - 1.0073 * peak.charge
                        # output_csv.write('%.5f\t%.4f\t%d\t%d\n' % (peak.mz, RT, peak.intensity, peak.charge))
                        yield (neutral, RT, peak.intensity, peak.charge)

                ms1_idx += 1
        idx += 1
    # output_csv.close()


def peptide_gen(settings):
    prefix = settings.get('input', 'decoy prefix')
    enzyme = get_enzyme(settings.get('search', 'enzyme'))
    mc = settings.getint('search', 'number of missed cleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')
    for prot in prot_gen(settings):
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

def prot_gen(settings):
    db = settings.get('input', 'database')
    add_decoy = settings.getboolean('input', 'add decoy')
    prefix = settings.get('input', 'decoy prefix')
    mode = settings.get('input', 'decoy method')

    read = [fasta.read, lambda f: fasta.decoy_db(f, mode=mode, prefix=prefix)][add_decoy]
    with read(db) as f:
        for p in f:
            yield p

seen_target = set()
seen_decoy = set()
def prot_peptides(prot_seq, enzyme, mc, minlen, maxlen, is_decoy):

    dont_use_fast_valid = parser.fast_valid(prot_seq)
    peptides = parser.cleave(prot_seq, enzyme, mc)
    for pep in peptides:
        plen = len(pep)
        if minlen <= plen <= maxlen + 2:
            forms = []
            if dont_use_fast_valid or pep in seen_target or pep in seen_decoy or parser.fast_valid(pep):
                # plen = len(pep)
                # if minlen <= plen <= maxlen:
                if plen <= maxlen:
                    forms.append(pep)
                if prot_seq[0] == 'M' and prot_seq.startswith(pep):
                    if minlen <= plen - 1 <= maxlen:
                        forms.append(pep[1:])
                    if minlen <= plen - 2 <= maxlen:
                        forms.append(pep[2:])
            for f in forms:
                if f not in seen_target and f not in seen_decoy:
                    if is_decoy:
                        seen_decoy.add(f)
                    else:
                        seen_target.add(f)
                    yield f

def get_prot_pept_map(settings):
    seen_target.clear()
    seen_decoy.clear()

    prefix = settings.get('input', 'decoy prefix')
    enzyme = get_enzyme(settings.get('search', 'enzyme'))
    mc = settings.getint('search', 'number of missed cleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')

    pept_prot = dict()
    protsN = dict()

    for desc, prot in prot_gen(settings):
        dbinfo = desc.split(' ')[0]
        for pep in prot_peptides(prot, get_enzyme(enzyme), mc, minlen, maxlen, desc.startswith(prefix)):
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
    if n == 0:
        try:
            n = cpu_count()
        except NotImplementedError:
            n = 1
    if n == 1:
        for s in it:
            yield func(s, **kw)
    else:
        def worker(qin, qout):
            for item in iter(qin.get, None):
                result = func(item, **kw)
                qout.put(result)
        qin = Queue()
        qout = Queue()
        count = 0
        while True:
            procs = []
            for _ in range(n):
                p = Process(target=worker, args=(qin, qout))
                p.start()
                procs.append(p)
            for s in it:
                qin.put(s)
                count += 1
                if count > 5000000:
                    print 'Loaded 5000000 items. Ending cycle.'
                    break
            for _ in range(n):
                qin.put(None)

            if not count:
                print 'No items left. Exiting.'
                break

            while count:
                yield qout.get()
                count -= 1

            for p in procs:
                p.join()

            print 'Cycle finished.'