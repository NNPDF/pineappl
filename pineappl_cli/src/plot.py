#!/usr/bin/env python3

import math
import matplotlib.pyplot as plt
import numpy as np
import pickle

def percent_diff(a, b):
    return (a / b - 1.0) * 100.0

def ylimits(axis):
    # extract the y limits *not* considering margins
    margins = axis.margins()
    axis.margins(y=0.0)
    min, max = axis.get_ylim()
    axis.margins(y=margins[1])

    inc = 1.0

    if (max - min) > 100.0:
        min = -50.0
        max = 50.0
        inc = 25.0
    elif (max - min) > 30.5:
        inc = 10.0
    elif (max - min) > 20.5:
        inc = 5.0
    elif (max - min) > 10.5:
        inc = 2.0
    elif (max - min) < 3.0:
        inc = 0.5

    min = math.floor(min / inc) * inc
    max = math.ceil(max / inc) * inc

    return [min, max, inc]

def plot_int(axis, **kwargs):
    xmin = np.array([])
    xmax = np.array([])
    x = np.array([])
    y = np.array([])

    for index, i in enumerate(kwargs['pdf_results']):
        label, ycentral, ymin, ymax = i
        x = np.append(x, ycentral[:-1])
        xmin = np.append(xmin, ymin[:-1])
        xmax = np.append(xmax, ymax[:-1])
        y = np.append(y, label)

        # draw one- and two-sigma bands
        if label == 'CENTRAL-PDF':
            axis.axvspan(xmin[-1], xmax[-1], alpha=0.3, color='royalblue', linewidth=0)
            # TODO: this is only correct for MC PDF uncertainties
            axis.axvspan(x[-1] - 2.0 * (x[-1] - xmin[-1]), x[-1] + 2.0 * (xmax[-1] - x[-1]), alpha=0.1, color='royalblue', linewidth=0)

    axis.errorbar(x, y, xerr=(x - xmin, xmax - x), fmt='.', capsize=3, markersize=5, linewidth=1.5)
    axis.margins(x=0.1, y=0.1)

def plot_abs(axis, **kwargs):
    x = kwargs['x']
    y = kwargs['y']
    ymin = kwargs['ymin']
    ymax = kwargs['ymax']
    ylog = kwargs['ylog']
    ylabel = kwargs['ylabel']
    slice_label = kwargs['slice_label']

    axis.set_yscale('log' if ylog else 'linear')
    axis.step(x, y, 'royalblue', linewidth=1.0, where='post', label=slice_label)
    axis.fill_between(x, ymin, ymax, alpha=0.4, color='royalblue', linewidth=0.5, step='post')
    axis.set_ylabel(ylabel)

    if slice_label != '':
        axis.legend()

def plot_ratio_pdf(axis, **kwargs):
    x = kwargs['x']
    ylog = kwargs['ylog']
    ylabel = kwargs['ylabel']
    slice_label = kwargs['slice_label']
    pdf_uncertainties = kwargs['pdf_results']
    channels = kwargs['channels']

    axis.set_ylabel('Ratio to ' + pdf_uncertainties[0][0])

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for index, i in enumerate(pdf_uncertainties):
        label, y, ymin, ymax = i
        y = y / pdf_uncertainties[0][1]
        ymin = ymin / pdf_uncertainties[0][1]
        ymax = ymax / pdf_uncertainties[0][1]

        axis.step(x, y, color=colors[index], linewidth=1.0, where='post')
        axis.fill_between(x, ymin, ymax, alpha=0.4, color=colors[index], label=label, linewidth=0.5, step='post')

    axis.legend(bbox_to_anchor=(0,-0.24,1,0.2), loc='upper left', mode='expand', borderaxespad=0, ncol=min(4, len(pdf_uncertainties) + 2))

    if slice_label != '':
        t = axis.text(0.98, 0.98, slice_label, horizontalalignment='right', verticalalignment='top', transform=axis.transAxes, fontsize='x-small')
        t.set_bbox({{ 'alpha': 0.7, 'boxstyle': 'square, pad=0.0', 'edgecolor': 'white', 'facecolor': 'white' }})

def plot_abs_pdfs(axis, **kwargs):
    x = kwargs['x']
    ylog = kwargs['ylog']
    ylabel = kwargs['ylabel']
    slice_label = kwargs['slice_label']
    pdf_uncertainties = kwargs['pdf_results']
    channels = kwargs['channels']

    axis.set_yscale('log' if ylog else 'linear')
    axis.set_ylabel(ylabel)

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for index, i in enumerate(pdf_uncertainties):
        label, y, ymin, ymax = i
        axis.step(x, y, color=colors[index], linewidth=1.0, where='post')
        axis.fill_between(x, ymin, ymax, alpha=0.4, color=colors[index], label=label, linewidth=0.5, step='post')

    linestyles = ['--', ':']
    for index, i in enumerate(channels):
        if index >= len(linestyles):
            break

        label, y = i
        axis.step(x, y, color=colors[0], label=label, linestyle=linestyles[index], linewidth=1.0, where='post')

    axis.legend(bbox_to_anchor=(0,-0.24,1,0.2), loc='upper left', mode='expand', borderaxespad=0, ncol=min(4, len(pdf_uncertainties) + 2))

    if slice_label != '':
        t = axis.text(0.98, 0.98, slice_label, horizontalalignment='right', verticalalignment='top', transform=axis.transAxes, fontsize='x-small')
        t.set_bbox({{ 'alpha': 0.7, 'boxstyle': 'square, pad=0.0', 'edgecolor': 'white', 'facecolor': 'white' }})

def plot_rel_ewonoff(axis, **kwargs):
    x = kwargs['x']
    y = percent_diff(kwargs['y'], kwargs['qcd_y'])
    qcd_y = percent_diff(kwargs['qcd_y'], kwargs['qcd_y'])
    qcd_ymin = percent_diff(kwargs['qcd_min'], kwargs['qcd_y'])
    qcd_ymax = percent_diff(kwargs['qcd_max'], kwargs['qcd_y'])
    ymin = percent_diff(kwargs['ymin'], kwargs['qcd_y'])
    ymax = percent_diff(kwargs['ymax'], kwargs['qcd_y'])
    pdf_min = abs(percent_diff(kwargs['pdf_results'][0][2], kwargs['pdf_results'][0][1]))[:-1]
    pdf_max = abs(percent_diff(kwargs['pdf_results'][0][3], kwargs['pdf_results'][0][1]))[:-1]
    mid = kwargs['mid']

    axis.step(x, qcd_y, 'red', label='NLO QCD', linewidth=1.0, where='post')
    #axis.fill_between(x, qcd_ymin, qcd_ymax, alpha=0.4, color='red', label='7-p.\ scale var.', linewidth=0.5, step='post')
    axis.step(x, y, 'royalblue', label='NLO QCD+EW', linewidth=1.0, where='post')
    axis.fill_between(x, ymin, ymax, alpha=0.4, color='royalblue', label='7-p.\ scale var.', linewidth=0.5, step='post')
    axis.errorbar(mid, y[:-1], yerr=(pdf_min, pdf_max), color='royalblue', label='PDF uncertainty', fmt='.', capsize=1, markersize=0, linewidth=1)
    axis.set_ylabel('NLO EW on/off [\si{{\percent}}]')
    axis.legend(bbox_to_anchor=(0,1.03,1,0.2), loc='lower left', mode='expand', borderaxespad=0, ncol=4)

def plot_rel_pdfunc(axis, **kwargs):
    x = kwargs['x']
    pdf_uncertainties = kwargs['pdf_results']
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    #ymins = np.asmatrix([(ymin / y - 1.0) * 100 for label, y, ymin, ymax in pdf_uncertainties])
    #ymaxs = np.asmatrix([(ymax / y - 1.0) * 100 for label, y, ymin, ymax in pdf_uncertainties])

    for index, i in enumerate(pdf_uncertainties):
        label, y, ymin, ymax = i
        ymin = percent_diff(ymin, y)
        ymax = percent_diff(ymax, y)
        axis.step(x, ymax, color=colors[index], label=label, linewidth=1, where='post')
        axis.step(x, ymin, color=colors[index], linewidth=1, where='post')

    #axis.legend(fontsize='xx-small') #rel_pdfunc
    axis.set_ylabel('PDF uncertainty [\si{{\percent}}]')

    this_ylim = ylimits(axis)

    if False:#SAVE-YLIM-PDFUNC
        with open('ylim-pdfunc', 'wb') as f:
            pickle.dump(this_ylim, f)

    if False:#LOAD-YLIM-PDFUNC
        resave = False

        with open('ylim-pdfunc', 'rb') as f:
            ylim = pickle.load(f)

        if ylim[0] < this_ylim[0]:
            this_ylim[0] = ylim[0]
            resave = True

        if ylim[1] > this_ylim[1]:
            this_ylim[1] = ylim[1]
            resave = True

        if ylim[2] > this_ylim[2]:
            this_ylim[2] = ylim[2]
            resave = True

        if resave:
            with open('ylim-pdfunc', 'wb') as f:
                pickle.dump(this_ylim, f)

    axis.set_yticks(np.arange(this_ylim[0], this_ylim[1] + this_ylim[2], this_ylim[2]))
    space = 0.05 * (this_ylim[1] - this_ylim[0])
    axis.set_ylim((this_ylim[0] - space, this_ylim[1] + space))

def plot_rel_pdfpull(axis, **kwargs):
    central_y = kwargs['pdf_results'][0][1]
    central_ymin = kwargs['pdf_results'][0][2]
    central_ymax = kwargs['pdf_results'][0][3]
    pdf_uncertainties = kwargs['pdf_results']
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    x = kwargs['x']
    y = kwargs['y']

    for index, i in enumerate(pdf_uncertainties):
        label, y, ymin, ymax = i
        diff = y - central_y
        yerr = np.where(diff > 0.0, y - ymin, ymax - y)
        cerr = np.where(diff > 0.0, central_ymax - central_y, central_y - central_ymin)
        pull = diff / np.sqrt(np.power(yerr, 2) + np.power(cerr, 2))

        #axis.fill_between(x, pull, pull_avg, alpha=0.4, color=colors[index], label='sym.\ pull', linewidth=0.5, step='post', zorder=2 * index)
        axis.step(x, pull, color=colors[index], label=label, linewidth=1, where='post', zorder=2 * index + 1)

    axis.legend(bbox_to_anchor=(0,1.03,1,0.2), loc='lower left', mode='expand', borderaxespad=0, ncol=len(pdf_uncertainties), fontsize='x-small', frameon=False, borderpad=0) #rel_pdfpull
    axis.set_ylabel('Pull [$\sigma$]')
    #axis.set_title('Comparison with ' + pdf_uncertainties[0][0], fontdict={{'fontsize': 9}}, loc='left')

    this_ylim = ylimits(axis)

    if False:#SAVE-YLIM-PDFPULL
        with open('ylim-pdfpull', 'wb') as f:
            pickle.dump(this_ylim, f)

    if False:#LOAD-YLIM-PDFPULL
        resave = False

        with open('ylim-pdfpull', 'rb') as f:
            ylim = pickle.load(f)

        if ylim[0] < this_ylim[0]:
            this_ylim[0] = ylim[0]
            resave = True

        if ylim[1] > this_ylim[1]:
            this_ylim[1] = ylim[1]
            resave = True

        if ylim[2] > this_ylim[2]:
            this_ylim[2] = ylim[2]
            resave = True

        if resave:
            with open('ylim-pdfpull', 'wb') as f:
                pickle.dump(this_ylim, f)

    axis.set_yticks(np.arange(this_ylim[0], this_ylim[1] + this_ylim[2], this_ylim[2]))
    space = 0.05 * (this_ylim[1] - this_ylim[0])
    axis.set_ylim((this_ylim[0] - space, this_ylim[1] + space))

def main():
    panels = [
        {inte}plot_int,
        {nint}plot_abs,
        {nint}plot_rel_ewonoff,
        {pdfs}plot_abs_pdfs,
        {pdfs}plot_ratio_pdf,
        {pdfs}plot_rel_pdfunc,
        {pdfs}plot_rel_pdfpull,
    ]

    plt.rc('axes', axisbelow=True, grid=True, labelsize='small')
    {nint}plt.rc('figure', figsize=(6.4,len(panels)*2.4))
    {inte}plt.rc('figure', figsize=(4.2,2.6))
    plt.rc('figure.constrained_layout', hspace=0.0, use=True, wspace=0.0)
    plt.rc('font', family='serif', size=14.0)
    plt.rc('grid', linestyle='dotted')
    plt.rc('legend', borderpad=0.0, fontsize='x-small', frameon=False)
    plt.rc('pdf', compression=0)
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{{siunitx}}\usepackage{{lmodern}}')
    plt.rc('xtick', bottom=True, top=True)
    plt.rc('xtick', direction='in')
    plt.rc('xtick.major', width=0.5)
    plt.rc('xtick.minor', bottom=True, top=True, width=0.5)
    plt.rc('ytick', direction='in')
    plt.rc('ytick', left=True, right=True)
    plt.rc('ytick.major', width=0.5)
    plt.rc('ytick.minor', visible=True, width=0.5)

    xaxis = '{xaxis}'
    xunit = '{xunit}'
    xlabel = r'{xlabel}'
    ylabel = r'{ylabel}'
    ylog = {ylog}
    description = r'{description}'

    data_slices = data()

    for index, dict in enumerate(data_slices):
        dict['xlabel'] = xlabel
        dict['ylabel'] = ylabel
        dict['ylog'] = ylog

        figure, axes = plt.subplots(len(panels), 1, sharex=True, squeeze=False)

        if len(dict['x']) != 2 and xunit != '':
            axes[0, 0].set_xscale('log')

        axes[ 0, 0].set_title(description)
        axes[-1, 0].set_xlabel(xlabel)

        for plot, axis in zip(panels, axes[:, 0]):
            plot(axis, **dict)

        name = '{output}' if len(data_slices) == 1 else '{output}-{{}}'.format(index)
        figure.savefig(name + '.pdf')
        plt.close(figure)

def data():
    return {data}

def metadata():
    return {{
{metadata}
    }}

if __name__ == '__main__':
    main()
