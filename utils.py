from ConfigParser import RawConfigParser
import numpy as np
from pyteomics import fasta, mzml, parser
from multiprocessing import Queue, Process, cpu_count
from scipy.stats import binom

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

def find_local_max_new(RTs, Is, min_I):
    maximums = []
    lv = RTs.searchsorted(RTs - 1.5, side='left')
    rv = RTs.searchsorted(RTs + 1.5, side='right')
    for idx in range(len(RTs)):
        maxI = max(Is[lv[idx]:rv[idx]])
        if maxI >= min_I and Is[idx] == maxI and lv[idx] != rv[idx] - 1 and RTs[rv[idx] - 1] - RTs[lv[idx]] >= 0.2:
            maximums.append(idx)
    return maximums



def find_peaks(q, q_out, min_I):
    for rec in iter(q.get, None):
        v, x, y, mzst_t, chs_t, pis_t, pI_t = rec
        idt = x.argsort()
        x = x[idt]
        y = y[idt]
        mzst_t = mzst_t[idt]
        chs_t = chs_t[idt]
        pis_t = pis_t[idt]
        pI_t = pI_t[idt]
        temp_idx = np.array([True for _ in x], dtype=bool)
        x = x[temp_idx]
        y = y[temp_idx]
        mzst_t = mzst_t[temp_idx]
        chs_t = chs_t[temp_idx]
        pis_t = pis_t[temp_idx]
        pI_t = pI_t[temp_idx]
        temp_idx = find_local_max_new(x, y, min_I=min_I)
        x_new2 = x[temp_idx]
        y_new2 = y[temp_idx]
        mzst_t = mzst_t[temp_idx]
        chs_t = chs_t[temp_idx]
        pis_t = pis_t[temp_idx]
        pI_t = pI_t[temp_idx]
        q_out.put(zip(mzst_t, x_new2, y_new2, chs_t, pis_t, pI_t))
    q_out.put(None)


def iterate_spectra(fname, min_ch, max_ch, min_i, nprocs, mass_acc):
    procs = []
    q = Queue()
    q_output = Queue()

    def custom_deiso(q, q_output, min_ch, max_ch, mass_acc):
        averagine_mass = 111.1254
        averagine_C = 4.9384
        tmplist = list(range(15))
        charges = list(range(min_ch, max_ch + 1, 1))[::-1]
        prec_masses = []
        prec_isotopes = []
        prec_minisotopes = []
        isotopes_int = []
        for i in range(300, 20000, 100):
            int_arr = binom.pmf(tmplist, float(i) / averagine_mass * averagine_C, 0.0107)
            prec_masses.append(i)
            int_arr_norm = int_arr / int_arr.max()
            prec_is = np.where(int_arr_norm >= 0.25)[0]
            isotopes_int.append(int_arr_norm[prec_is])
            prec_minisotopes.append(prec_is.min())
            prec_isotopes.append(prec_is - prec_minisotopes[-1])
        I_err = 1.5
        prec_masses = np.array(prec_masses)
        for sc in iter(q.get, None):
            banned = set()
            RT = sc['scanList']['scan'][0]['scan start time']
            mz = sc['m/z array']
            Intensities = sc['intensity array']
            tmpidx = Intensities>=10
            mz = mz[tmpidx]
            Intensities = Intensities[tmpidx]

            mz_size = mz.size
            for idx, v in enumerate(mz):
                mz_tol = mass_acc * 1e-6 * v
                for ch in charges:
                    pos_ind = prec_masses.searchsorted(v * ch)
                    shifts = prec_isotopes[pos_ind]
                    found_isotopes = [idx, ]
                    min_shift = prec_minisotopes[pos_ind]
                    flag = 1
                    for shift in shifts[1:]:
                        m = v + (1.00335 * shift / ch)
                        j = mz.searchsorted(m)
                        if j == mz_size or m - mz[j - 1] < mz[j] - m:
                            j -= 1
                        if j not in banned and abs(m - mz[j]) <= mz_tol:
                            found_isotopes.append(j)
                        else:
                            flag = 0
                            break
                    if flag:
                        int_approved = [found_isotopes[0]]
                        for jj_ind, jj in enumerate(found_isotopes[1:]):
                            if 1 / I_err <= (Intensities[found_isotopes[jj_ind]] / Intensities[jj]) / (
                                        isotopes_int[pos_ind][jj_ind] / isotopes_int[pos_ind][
                                    jj_ind + 1]) <= I_err:
                                int_approved.append(jj)
                            elif jj_ind == 0:
                                break
                        if len(int_approved) > 1:
                            for jj in int_approved:
                                banned.add(jj)
                            q_output.put((v - (1.00335 * min_shift / ch), RT, Intensities[idx], ch))
                            break
        q_output.put(None)

    for i in range(nprocs):
        p = Process(target=custom_deiso, args=(q, q_output, min_ch, max_ch, mass_acc))
        procs.append(p)
        p.start()


    idx = 1
    ms1_idx = 1

    for sc in mzml.read(fname):
        idx += 1
        if sc['ms level'] == 1:
            q.put(sc)
            ms1_idx += 1
    print 'total number of spectra = %d\ntotal number of ms1 spectra = %d' % (idx - 1, ms1_idx - 1)
    for p in procs:
        q.put(None)

    results = []

    jj = 0
    while jj < nprocs:
        for rec in iter(q_output.get, None):
            results.append(rec)
        jj += 1

    for p in procs:
        p.terminate()
    print 'Deisotoping is finished'

    peak_id = -1
    number_of_peaks = len(results)

    mzst = np.empty(number_of_peaks, dtype=np.dtype('f4'))
    mzs = np.empty(number_of_peaks, dtype=np.dtype('f4'))
    RTs = np.empty(number_of_peaks, dtype=np.dtype('f2'))
    Is = np.empty(number_of_peaks, dtype=int)
    chs = np.empty(number_of_peaks, dtype=np.dtype('i2'))
    pis = np.empty(number_of_peaks, dtype=np.dtype('i4'))
    pIs = np.empty(number_of_peaks, dtype=np.dtype('f4'))

    for rec in results:
        mz = float(rec[0])
        RT = float(rec[1])
        I = float(rec[2])
        cz = int(rec[3])
        try:
            pI = float(rec[6])
        except:
            pI = 7
        neutral = mz * cz - 1.0073 * cz
        peak_id += 1
        mzst[peak_id] = neutral
        mzs[peak_id] = mz
        RTs[peak_id] = RT
        Is[peak_id] = I
        chs[peak_id] = cz
        pis[peak_id] = peak_id
        pIs[peak_id] = pI
    print len(mzs), 'len(mzs)'
    print peak_id, 'peak_id'
    idx = np.argsort(mzs)
    mzs = mzs[idx]
    mzst = mzst[idx]
    RTs = RTs[idx]
    Is = Is[idx]
    chs = chs[idx]
    pis = pis[idx]
    pIs = pIs[idx]
    minv = -1



    procs = []
    q = Queue()
    q_out = Queue()

    for i in range(nprocs):
        p = Process(target=find_peaks, args=(q, q_out, min_i))
        procs.append(p)
        p.start()

    lv = mzs.searchsorted(mzs - 0.01, side='left')
    rv = mzs.searchsorted(mzs + 0.01, side='right')
    for idx in range(len(mzs)):
        if mzs[idx] >= minv:
            minv = mzs[idx] + 0.005
            if mzs[lv[idx]:rv[idx]].size > 5 and RTs[lv[idx]:rv[idx]].max() - RTs[lv[idx]:rv[idx]].min() > 0.25:
                q.put((mzs[idx], RTs[lv[idx]:rv[idx]], Is[lv[idx]:rv[idx]], mzst[lv[idx]:rv[idx]], chs[lv[idx]:rv[idx]], pis[lv[idx]:rv[idx]], pIs[lv[idx]:rv[idx]]))
    for p in procs:
        q.put(None)

    outres = set()
    jj = 0
    while jj < nprocs:
        for rec in iter(q_out.get, None):
            outres.update(rec)
        jj += 1

    for p in procs:
        p.terminate()

    print 'Total number of peaks = %d' % (len(mzs), )
    print 'Total number of unique peaks = %d' % (len(outres), )
    for res in outres:
        yield res


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