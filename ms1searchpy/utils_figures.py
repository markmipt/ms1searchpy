from __future__ import division
from scipy.stats import scoreatpercentile
from pyteomics import mass
from collections import Counter
import os.path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
try:
    import seaborn
    seaborn.set(rc={'axes.facecolor':'#ffffff'})
    seaborn.set_style('whitegrid')
except ImportError:
    pass
import re
import sys
import logging
import pandas as pd
logger = logging.getLogger(__name__)
logging.getLogger('matplotlib.font_manager').disabled = True
redcolor = '#FC6264'
bluecolor = '#70aed1'
greencolor = '#8AA413'
aa_color_1 = '#fca110'
aa_color_2 = '#a41389'

def _get_sf(fig):
    return isinstance(fig, str) if sys.version_info.major == 3 else isinstance(fig, basestring)


def get_basic_distributions(df):
    if 'mz' in df.columns:
        mz_array = df['mz'].values
        rt_exp_array = df['rtApex'].values
        intensity_array = np.log10(df['intensityApex'].values)
        if 'FAIMS' in df.columns:
            faims_array = df['FAIMS'].replace('None', 0).values
        else:
            faims_array = np.zeros(len(df))
    else:
        mz_array = df['m/z'].values
        rt_exp_array = df['RT'].values
        intensity_array = np.log10(df['Intensity'].values)
        faims_array = df['ion_mobility'].values
    charges_array = df['charge'].values
    nScans_array = np.log2(df['nScans'].values)
    nIsotopes_array = df['nIsotopes'].values
    return mz_array, rt_exp_array, charges_array, intensity_array, nScans_array, nIsotopes_array, faims_array


def get_descriptor_array(df, df_f, dname):
    array_t = df[~df['decoy']][dname].values
    array_d = df[df['decoy']][dname].values
    array_v = df_f[dname].values
    return array_t, array_d, array_v


def plot_hist_basic(array_all, fig, subplot_max_x, subplot_i,
        xlabel, ylabel='# of identifications', idtype='features', bin_size_one=False):

    fig.add_subplot(subplot_max_x, 3, subplot_i)
    cbins = get_bins((array_all, ), bin_size_one)
    if idtype == 'features':
        plt.hist(array_all, bins=cbins, color=redcolor, alpha=0.8, edgecolor='#EEEEEE')
    else:
        plt.hist(array_all, bins=cbins, color=greencolor, alpha=0.8, edgecolor='#EEEEEE')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

def plot_basic_figures(df, fig, subplot_max_x, subplot_start, idtype):
    logger.debug('Plotting %s figures...', idtype)
    logger.debug('Getting basic distributions (all)')
    mz_array, rt_exp_array, lengths_array, intensity_array, nScans_array, nIsotopes_array, faims_array = get_basic_distributions(df)
    # logger.debug('Getting basic distributions (filtered)')
    # mz_array_valid, rt_exp_array_valid, lengths_array_valid, intensity_array_valid, nScans_array_valid, nIsotopes_array_valid, faims_array_valid = get_basic_distributions(df_f)

    logger.debug('Plotting %s m/z distributions.', idtype)
    plot_hist_basic(mz_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, precursor m/z' % (idtype, ), idtype=idtype)
    subplot_start += 1
    logger.debug('Plotting %s RT distributions.', idtype)
    plot_hist_basic(rt_exp_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, RT experimental' % (idtype, ), idtype=idtype)
    subplot_start += 1
    logger.debug('Plotting %s charge distributions.', idtype)
    plot_hist_basic(lengths_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, charge' % (idtype, ), idtype=idtype, bin_size_one=True)
    subplot_start += 1
    logger.debug('Plotting %s intensity distributions.', idtype)
    plot_hist_basic(intensity_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, log10(Intensity)' % (idtype, ), idtype=idtype)
    subplot_start += 1
    logger.debug('Plotting %s scans distributions.', idtype)
    plot_hist_basic(nScans_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, log2(nScans)' % (idtype, ), idtype=idtype)
    subplot_start += 1
    logger.debug('Plotting %s nIsotopes distributions.', idtype)
    plot_hist_basic(nIsotopes_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, nIsotopes' % (idtype, ), idtype=idtype, bin_size_one=True)
    subplot_start += 1
    logger.debug('Plotting %s ion mobility distributions.', idtype)
    plot_hist_basic(faims_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, ion mobility' % (idtype, ), idtype=idtype)


def plot_protein_figures(df, df_f, fig, subplot_max_x, subplot_start):
    logger.debug('Plotting protein figures: NSAF')
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='LOG10_NSAF'), fig, subplot_max_x, subplot_start, xlabel='LOG10(NSAF)')
    logger.debug('Plotting protein figures: SQ')
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='sq'), fig, subplot_max_x, subplot_start+1, xlabel='sequence coverage')


def plot_hist_descriptor(inarrays, fig, subplot_max_x, subplot_i, xlabel, ylabel='# of identifications'):
    logger.debug('Plotting descriptor histogram: %s', xlabel)
    separate_figures = _get_sf(fig)
    if separate_figures:
        plt.figure()
    else:
        fig.add_subplot(subplot_max_x, 3, subplot_i)
    array_t, array_d, array_v = inarrays
    cbins, width = get_bins_for_descriptors(inarrays)
    H1, _ = np.histogram(array_d, bins=cbins)
    H2, _ = np.histogram(array_t, bins=cbins)
    H3, _ = np.histogram(array_v, bins=cbins)
    plt.bar(cbins[:-1], H1, width, align='center',color=redcolor, alpha=0.4, edgecolor='#EEEEEE')
    plt.bar(cbins[:-1], H2, width, align='center',color=bluecolor, alpha=0.4, edgecolor='#EEEEEE')
    plt.bar(cbins[:-1], H3, width, align='center',color=greencolor, alpha=1, edgecolor='#EEEEEE')
    cbins = np.append(cbins[0], cbins)
    H1 = np.append(np.append(0, H1), 0)
    H2 = np.append(np.append(0, H2), 0)
    cbins -= width / 2
    plt.step(cbins, H2, where='post', color=bluecolor, alpha=0.8)
    plt.step(cbins, H1, where='post', color=redcolor, alpha=0.8)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    if 'mass shift' in xlabel:
        plt.xlim(-1.5, 1.5)
        plt.xticks([-1, 0, 1], ['NA', 'unmodified', 'modified'])
    elif width == 1.0:
        plt.xticks(np.arange(int(cbins[0]), cbins[-1], 1))
        plt.gcf().canvas.draw()
    if separate_figures:
        plt.savefig(outpath(fig, xlabel, '.png'))
        plt.close()


def plot_legend(fig, subplot_max_x, subplot_start):
    ax = fig.add_subplot(subplot_max_x, 3, subplot_start)
    legend_elements = [Patch(facecolor=greencolor, label='Positive IDs'),
                    Patch(facecolor=bluecolor, label='Targets'),
                    Patch(facecolor=redcolor, label='Decoys')]
    ax.legend(handles=legend_elements, loc='center', prop={'size': 24})
    ax.set_axis_off()


def plot_aa_stats(df_f, df_proteins_f, fig, subplot_max_x, subplot_i):
    separate_figures = _get_sf(fig)
    if separate_figures:
        plt.figure()
    else:
        fig.add_subplot(subplot_max_x, 3, subplot_i)

    # Generate list of 20 standart amino acids
    std_aa_list = list(mass.std_aa_mass.keys())
    std_aa_list.remove('O')
    std_aa_list.remove('U')

    # Count identified amino acids
    aa_exp = Counter()
    for pep in set(df_f['peptide']):
        for aa in pep:
            aa_exp[aa] += 1

    # Count theoretical amino acids
    aa_theor = Counter()
    for prot_seq in set(df_proteins_f['sequence'].values):
        for aa in prot_seq:
            aa_theor[aa] += 1

    aa_exp_sum = sum(aa_exp.values())
    aa_theor_sum = sum(aa_theor.values())
    lbls, vals = [], []
    for aa in sorted(std_aa_list):
        if aa_theor.get(aa, 0):
            lbls.append(aa)
            vals.append((aa_exp.get(aa, 0)/aa_exp_sum) / (aa_theor[aa]/aa_theor_sum))
    # std_val = np.std(vals)
    clrs = [aa_color_1 if abs(x-1) < 0.4 else aa_color_2 for x in vals]
    plt.bar(range(len(vals)), vals, color=clrs)
    plt.xticks(range(len(lbls)), lbls)
    plt.hlines(1.0, range(len(vals))[0]-1, range(len(vals))[-1]+1)
    plt.ylabel('amino acid ID rate')


def calc_max_x_value(df, df_proteins):
    cnt = 14 # number of basic figures
    peptide_columns = set(df.columns)
    features_list = []#['massdiff_ppm', 'RT diff', 'fragmentMT', 'num_missed_cleavages', 'assumed_charge', 'log_score', 'ISOWIDTHDIFF', 'MS1Intensity']
    for feature in features_list:
        if feature in peptide_columns:
            cnt += 1
    # for feature in peptide_columns:
    #     if feature.startswith('mass shift'):
    #         cnt += 1
    # if len(set(df['massdiff_int'])) > 1:
    #     cnt += 1
    # if 'LOG10_NSAF' in df_proteins.columns:
    #     cnt += 3 # add for NSAF, sequence coverage and aa_stats
    return cnt // 3 + (1 if (cnt % 3) else 0)


def plot_descriptors_figures(df, df_f, fig, subplot_max_x, subplot_start):
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='massdiff_ppm'), fig, subplot_max_x, subplot_start, xlabel='precursor mass difference, ppm')
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='RT diff'), fig, subplot_max_x, subplot_start, xlabel='RT difference, min')
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='num_missed_cleavages'), fig, subplot_max_x, subplot_start, xlabel='missed cleavages')
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='assumed_charge'), fig, subplot_max_x, subplot_start, xlabel='precursor charge')
    subplot_start += 1
    if 'fragmentMT' in df.columns:
        plot_hist_descriptor(get_descriptor_array(df, df_f, dname='fragmentMT'), fig, subplot_max_x, subplot_start, xlabel='median fragment error, Da')
        subplot_start += 1
    if 'ISOWIDTHDIFF' in df.columns:
        plot_hist_descriptor(get_descriptor_array(df, df_f, dname='ISOWIDTHDIFF'), fig, subplot_max_x, subplot_start, xlabel='isolation mass error, Da')
        subplot_start += 1
    if 'MS1Intensity' in df.columns:
        plot_hist_descriptor(get_descriptor_array(df, df_f, dname='MS1Intensity'), fig, subplot_max_x, subplot_start, xlabel='MS1 Intensity')
        subplot_start += 1

    if len(set(df['massdiff_int'])) > 1:
        plot_hist_descriptor(get_descriptor_array(df, df_f, dname='massdiff_int'), fig, subplot_max_x, subplot_start, xlabel='isotope mass difference, Da')
        subplot_start += 1
    for df_col in df.columns:
        if df_col.startswith('mass shift'):
            plot_hist_descriptor(get_descriptor_array(df, df_f, dname=df_col), fig, subplot_max_x, subplot_start, xlabel=df_col)
            subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='log_score'), fig, subplot_max_x, subplot_start, xlabel='LOG10(ML score)')
    subplot_start += 1
    separate_figures = _get_sf(fig)
    if not separate_figures:
        plot_legend(fig, subplot_max_x, subplot_start)
    subplot_start += 1


def get_bins(inarrays, bin_size_one=False):
    tmp = np.concatenate(inarrays)
    minv = tmp.min()
    maxv = tmp.max()
    if bin_size_one:
        return np.arange(minv, maxv+1, 1)
    else:
        return np.linspace(minv, maxv+1, num=100)


def get_bins_for_descriptors(inarrays):
    tmp = np.concatenate(inarrays)
    minv = tmp.min()
    maxv = tmp.max()
    if len(set(tmp)) <= 15:
        return np.arange(minv, maxv + 2, 1.0), 1.0
    binsize = False
    for inar in inarrays:
        binsize_tmp = get_fdbinsize(inar)
        if not binsize or binsize > binsize_tmp:
            binsize = binsize_tmp
    # binsize = get_fdbinsize(tmp)
    if binsize < float(maxv - minv) / 300:
        binsize = float(maxv - minv) / 300

    lbin_s = scoreatpercentile(tmp, 1.0)
    lbin = minv
    if lbin_s and abs((lbin - lbin_s) / lbin_s) > 1.0:
        lbin = lbin_s * 1.05
    rbin_s = scoreatpercentile(tmp, 99.0)
    rbin = maxv
    if rbin_s and abs((rbin - rbin_s) / rbin_s) > 1.0:
        rbin = rbin_s * 1.05
    rbin += 1.5 * binsize
    logger.debug('get_bins_for_descriptors: lbin = %s, rbin = %s, binsize = %s', lbin, rbin, binsize)
    return np.arange(lbin, rbin + binsize, binsize), binsize


def get_fdbinsize(data_list):
    """Calculates the Freedman-Diaconis bin size for
    a data set for use in making a histogram
    Arguments:
    data_list:  1D Data set
    Returns:
    optimal_bin_size:  F-D bin size
    """
    if not isinstance(data_list, np.ndarray):
        data_list = np.array(data_list)
    isnan = np.isnan(data_list)
    logger.debug('Number of nans: %s', isnan.sum())
    data_list = np.sort(data_list[~isnan])
    upperquartile = scoreatpercentile(data_list, 75)
    lowerquartile = scoreatpercentile(data_list, 25)
    iqr = upperquartile - lowerquartile
    logger.debug('IQR: %s, data size: %s', iqr, len(data_list))
    optimal_bin_size = 2. * iqr / len(data_list) ** (1. / 3.)
    logger.debug('Calculated optimal bin size: %s.', optimal_bin_size)
    MINBIN = 1e-8
    if optimal_bin_size < MINBIN:
        logger.debug('Increasing bin size to %s.', MINBIN)
        return MINBIN
    return optimal_bin_size


def normalize_fname(s):
    return re.sub(r'[<>:\|/?*]', '', s)


def outpath(outfolder, s, ext='.png'):
    return os.path.join(outfolder, normalize_fname(s) + ext)


def plot_outfigures(df, df_peptides, df_peptides_f, base_out_name, df_proteins, df_proteins_f):
    fig = plt.figure(figsize=(16, 12))
    dpi = fig.get_dpi()
    fig.set_size_inches(3000.0/dpi, 3000.0/dpi)
    subplot_max_x = calc_max_x_value(df, df_proteins)
    descriptor_start_index = 15
    logger.debug('Plotting feature figures...')
    plot_basic_figures(df, fig, subplot_max_x, 1, 'features')
    logger.debug('Plotting peptide figures...')
    plot_basic_figures(df_peptides_f, fig, subplot_max_x, 8, 'peptides')
    # if 'LOG10_NSAF' in df_proteins.columns:
    #     logger.debug('Plotting protein figures...')
    #     try:
    #         plot_protein_figures(df_proteins, df_proteins_f, fig, subplot_max_x, 4)
    #         logger.debug('Plotting AA stats figures...')
    #         plot_aa_stats(df_f, df_proteins_f, fig, subplot_max_x, 6)
    #         descriptor_start_index += 3
    #     except ZeroDivisionError:
    #         logger.warning('There was an error in protein stats calculation. No protein figures will be plotted.')
    # logger.debug('Plotting descriptor figures...')
    # plot_descriptors_figures(df, df_f, fig, subplot_max_x, descriptor_start_index)
    plt.grid(color='#EEEEEE')
    plt.tight_layout()
    plt.savefig(base_out_name + '.png')
    plt.close()
    logger.info('Figures saved.')
