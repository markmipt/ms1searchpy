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
redcolor = '#FC6264'
bluecolor = '#70aed1'
greencolor = '#8AA413'
aa_color_1 = '#fca110'
aa_color_2 = '#a41389'


def get_basic_distributions(df):
    if 'mz' in df.columns:
        mz_array = df['mz'].values
        rt_exp_array = df['rtApex'].values
        intensity_array = np.log10(df['intensityApex'].values)
        if 'FAIMS' in df.columns:
            faims_array = df['FAIMS'].replace('None', 0).values
        else:
            faims_array = np.zeros(len(df))
        mass_diff_array = []
        RT_diff_array = []
    else:
        mz_array = df['m/z'].values
        rt_exp_array = df['RT'].values
        intensity_array = np.log10(df['Intensity'].values)
        faims_array = df['ion_mobility'].values
        mass_diff_array = df['mass diff'].values
        RT_diff_array = df['RT diff'].values
    charges_array = df['charge'].values
    nScans_array = df['nScans'].values
    nIsotopes_array = df['nIsotopes'].values
    return mz_array, rt_exp_array, charges_array, intensity_array, nScans_array, nIsotopes_array, faims_array, mass_diff_array, RT_diff_array


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
    mz_array, rt_exp_array, lengths_array, intensity_array, nScans_array, nIsotopes_array, faims_array, mass_diff_array, RT_diff_array = get_basic_distributions(df)

    plot_hist_basic(mz_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, precursor m/z' % (idtype, ), idtype=idtype)
    subplot_start += 1
    plot_hist_basic(rt_exp_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, RT experimental' % (idtype, ), idtype=idtype)
    subplot_start += 1
    plot_hist_basic(lengths_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, charge' % (idtype, ), idtype=idtype, bin_size_one=True)
    subplot_start += 1
    plot_hist_basic(intensity_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, log10(Intensity)' % (idtype, ), idtype=idtype)
    subplot_start += 1
    plot_hist_basic(nScans_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, nScans' % (idtype, ), idtype=idtype, bin_size_one=True)
    subplot_start += 1
    plot_hist_basic(nIsotopes_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, nIsotopes' % (idtype, ), idtype=idtype, bin_size_one=True)
    subplot_start += 1
    plot_hist_basic(faims_array, fig, subplot_max_x, subplot_i=subplot_start,
                     xlabel='%s, ion mobility' % (idtype, ), idtype=idtype)
    if idtype=='peptides':
        subplot_start += 1
        plot_hist_basic(mass_diff_array, fig, subplot_max_x, subplot_i=subplot_start,
                        xlabel='%s, mass diff ppm' % (idtype, ), idtype=idtype)
        subplot_start += 1
        plot_hist_basic(RT_diff_array, fig, subplot_max_x, subplot_i=subplot_start,
                        xlabel='%s, RT diff min' % (idtype, ), idtype=idtype)


def plot_protein_figures(df, df_f, fig, subplot_max_x, subplot_start):
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='sq'), fig, subplot_max_x, subplot_start, xlabel='proteins, sequence coverage', ylabel='# of identifications', only_true=True)
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='matched peptides'), fig, subplot_max_x, subplot_start, xlabel='proteins, matched peptides', ylabel='# of identifications', only_true=True, bin_size_one=True)
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='corrected sq'), fig, subplot_max_x, subplot_start, xlabel='proteins, corrected sequence coverage', ylabel='# of identifications', only_true=True)
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='corrected matched peptides'), fig, subplot_max_x, subplot_start, xlabel='proteins, corrected matched peptides', ylabel='# of identifications', only_true=True, bin_size_one=True)
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='score'), fig, subplot_max_x, subplot_start, xlabel='proteins, score', ylabel='# of identifications', only_true=False)
    subplot_start += 1
    return subplot_start


def plot_hist_descriptor(inarrays, fig, subplot_max_x, subplot_i, xlabel, ylabel='# of identifications', only_true=False, bin_size_one=False):
    fig.add_subplot(subplot_max_x, 3, subplot_i)
    array_t, array_d, array_v = inarrays
    if xlabel == 'proteins, score':
        logscale=True
    else:
        logscale=False
    if only_true:
        cbins, width = get_bins_for_descriptors([array_v, ], bin_size_one=bin_size_one)
        H3, _ = np.histogram(array_v, bins=cbins)
        plt.bar(cbins[:-1], H3, width, align='center',color=greencolor, alpha=1, edgecolor='#EEEEEE', log=logscale)
    else:
        cbins, width = get_bins_for_descriptors(inarrays)
        H1, _ = np.histogram(array_d, bins=cbins)
        H2, _ = np.histogram(array_t, bins=cbins)
        H3, _ = np.histogram(array_v, bins=cbins)
        plt.bar(cbins[:-1], H1, width, align='center',color=redcolor, alpha=0.4, edgecolor='#EEEEEE', log=logscale)
        plt.bar(cbins[:-1], H2, width, align='center',color=bluecolor, alpha=0.4, edgecolor='#EEEEEE', log=logscale)
        plt.bar(cbins[:-1], H3, width, align='center',color=greencolor, alpha=1, edgecolor='#EEEEEE', log=logscale)
        cbins = np.append(cbins[0], cbins)
        H1 = np.append(np.append(0, H1), 0)
        H2 = np.append(np.append(0, H2), 0)
        cbins -= width / 2
        plt.step(cbins, H2, where='post', color=bluecolor, alpha=0.8)
        plt.step(cbins, H1, where='post', color=redcolor, alpha=0.8)



    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    if width == 1.0:
        plt.xticks(np.arange(int(cbins[0]), cbins[-1], 1))
        plt.gcf().canvas.draw()


def plot_legend(fig, subplot_max_x, subplot_start):
    ax = fig.add_subplot(subplot_max_x, 3, subplot_start)
    legend_elements = [Patch(facecolor=greencolor, label='Positive IDs'),
                    Patch(facecolor=bluecolor, label='Targets'),
                    Patch(facecolor=redcolor, label='Decoys')]
    ax.legend(handles=legend_elements, loc='center', prop={'size': 24})
    ax.set_axis_off()


def plot_aa_stats(df_f, df, fig, subplot_max_x, subplot_i):
    fig.add_subplot(subplot_max_x, 3, subplot_i)

    # Generate list of 20 standart amino acids
    std_aa_list = list(mass.std_aa_mass.keys())
    std_aa_list.remove('O')
    std_aa_list.remove('U')

    # Count identified amino acids
    aa_exp = Counter()
    for pep in set(df_f['sequence']):
        for aa in pep:
            aa_exp[aa] += 1

    # Count theoretical amino acids
    aa_theor = Counter()
    for pep in set(df[df['decoy']]['sequence']):
        for aa in pep:
            aa_theor[aa] += 1

    aa_exp_sum = sum(aa_exp.values())
    aa_theor_sum = sum(aa_theor.values())
    lbls, vals = [], []
    for aa in sorted(std_aa_list):
        if aa_theor.get(aa, 0):
            lbls.append(aa)
            vals.append((aa_exp.get(aa, 0)/aa_exp_sum) / (aa_theor[aa]/aa_theor_sum))
    clrs = [aa_color_1 if abs(x-1) < 0.4 else aa_color_2 for x in vals]
    plt.bar(range(len(vals)), vals, color=clrs)
    plt.xticks(range(len(lbls)), lbls)
    plt.hlines(1.0, range(len(vals))[0]-1, range(len(vals))[-1]+1)
    plt.ylabel('amino acid ID rate')


def calc_max_x_value(df, df_proteins):
    cnt = 23 # number of basic figures
    peptide_columns = set(df.columns)
    features_list = []
    for feature in features_list:
        if feature in peptide_columns:
            cnt += 1
    return cnt // 3 + (1 if (cnt % 3) else 0)


def plot_descriptors_figures(df, df_f, fig, subplot_max_x, subplot_start):
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='mass diff'), fig, subplot_max_x, subplot_start, xlabel='precursor mass difference, ppm')
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='RT diff'), fig, subplot_max_x, subplot_start, xlabel='RT difference, min')
    subplot_start += 1
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


def get_bins_for_descriptors(inarrays, bin_size_one=False):
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

    if bin_size_one:
        binsize = int(binsize)
        if binsize == 0:
            binsize = 1

    lbin_s = scoreatpercentile(tmp, 1.0)
    lbin = minv
    if lbin_s and abs((lbin - lbin_s) / lbin_s) > 1.0:
        lbin = lbin_s * 1.05
    rbin_s = scoreatpercentile(tmp, 99.0)
    rbin = maxv
    if rbin_s and abs((rbin - rbin_s) / rbin_s) > 1.0:
        rbin = rbin_s * 1.05
    rbin += 1.5 * binsize
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
    data_list = np.sort(data_list[~isnan])
    upperquartile = scoreatpercentile(data_list, 75)
    lowerquartile = scoreatpercentile(data_list, 25)
    iqr = upperquartile - lowerquartile
    optimal_bin_size = 2. * iqr / len(data_list) ** (1. / 3.)
    MINBIN = 1e-8
    if optimal_bin_size < MINBIN:
        return MINBIN
    return optimal_bin_size


def plot_outfigures(df, df_peptides, df_peptides_f, base_out_name, df_proteins, df_proteins_f):
    fig = plt.figure(figsize=(16, 12))
    dpi = fig.get_dpi()
    fig.set_size_inches(3000.0/dpi, 3000.0/dpi)
    subplot_max_x = calc_max_x_value(df, df_proteins)
    descriptor_start_index = 20
    plot_basic_figures(df, fig, subplot_max_x, 1, 'features')
    plot_basic_figures(df_peptides_f, fig, subplot_max_x, 8, 'peptides')
    df_proteins['sq'] = df_proteins['matched peptides'] / df_proteins['theoretical peptides'] * 100
    df_proteins_f['sq'] = df_proteins_f['matched peptides'] / df_proteins_f['theoretical peptides'] * 100

    p_decoy = df_proteins[df_proteins['decoy']]['matched peptides'].sum() / df_proteins[df_proteins['decoy']]['theoretical peptides'].sum()
    df_proteins['corrected matched peptides'] = df_proteins['matched peptides'] - p_decoy * df_proteins['theoretical peptides']
    df_proteins_f['corrected matched peptides'] = df_proteins_f['matched peptides'] - p_decoy * df_proteins_f['theoretical peptides']
    df_proteins['corrected sq'] = df_proteins['corrected matched peptides'] / df_proteins['theoretical peptides'] * 100
    df_proteins_f['corrected sq'] = df_proteins_f['corrected matched peptides'] / df_proteins_f['theoretical peptides'] * 100

    subplot_current = plot_protein_figures(df_proteins, df_proteins_f, fig, subplot_max_x, 17)
    plot_aa_stats(df_peptides_f, df_peptides, fig, subplot_max_x, subplot_current)
    plt.grid(color='#EEEEEE')
    plt.tight_layout()
    plt.savefig(base_out_name + '.png')
    plt.close()
