use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use itertools::Itertools;
use lhapdf::{Pdf, PdfSet};
use pineappl::bin::BinInfo;
use pineappl::subgrid::Subgrid;
use rayon::prelude::*;
use std::path::{Path, PathBuf};

/// Creates a matplotlib script plotting the contents of the grid.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id(s) or name of the PDF set(s).
    #[clap(min_values = 1, validator = helpers::validate_pdfset)]
    pdfsets: Vec<String>,
    /// Set the number of scale variations.
    #[clap(default_value = "7", long, possible_values = &["1", "3", "7", "9"], short)]
    scales: usize,
    /// Show the pull for a specific grid three-dimensionally.
    #[clap(
        conflicts_with = "scales",
        long = "subgrid-pull",
        number_of_values = 3,
        use_delimiter = true,
        value_names = &["ORDER", "BIN", "LUMI"]
    )]
    subgrid_pull: Vec<String>,
}

fn pdfset_label(pdfset: &str) -> &str {
    pdfset.rsplit_once('=').map_or(pdfset, |(_, label)| label)
}

fn pdfset_name(pdfset: &str) -> &str {
    pdfset.rsplit_once('=').map_or(pdfset, |(name, _)| name)
}

fn map_format_join(slice: &[f64]) -> String {
    slice.iter().map(|x| format!("{}", x)).join(", ")
}

fn map_format_e_join(slice: &[f64]) -> String {
    slice.iter().map(|x| format!("{:e}", x)).join(", ")
}

fn format_pdf_results(pdf_uncertainties: &[Vec<Vec<f64>>], pdfsets: &[String]) -> String {
    let mut result = String::new();

    for (values, pdfset) in pdf_uncertainties.iter().zip(pdfsets.iter()) {
        result.push_str(&format!(
            "        (
            '{}',
            np.array([{}]),
            np.array([{}]),
            np.array([{}]),
        ),\n",
            pdfset_label(pdfset).replace('_', "\\_"),
            map_format_e_join(&values[0]),
            map_format_e_join(&values[1]),
            map_format_e_join(&values[2]),
        ));
    }

    result
}

fn format_metadata(metadata: &[(&String, &String)]) -> String {
    let mut result = String::new();

    for (key, value) in metadata {
        // skip multi-line entries
        if value.contains('\n') {
            continue;
        }

        result.push_str(&format!(
            "        '{}': r'{}',\n",
            key,
            if *key == "description" {
                value.replace("\u{2013}", "--").replace("\u{2014}", "---")
            } else if key.ends_with("_unit") {
                value
                    .replace("GeV", r#"\giga\electronvolt"#)
                    .replace("/", r#"\per"#)
                    .replace("pb", r#"\pico\barn"#)
            } else {
                (*value).clone()
            }
        ));
    }

    result
}

fn format_script(
    bin_info: &BinInfo,
    output: &str,
    left: &[f64],
    right: &[f64],
    min: &[f64],
    max: &[f64],
    qcd_central: &[f64],
    qcd_min: &[f64],
    qcd_max: &[f64],
    slices: &[(usize, usize)],
    slice_labels: &[String],
    pdf_uncertainties: &[Vec<Vec<f64>>],
    pdfsets: &[String],
    metadata: &[(&String, &String)],
) {
    println!("#!/usr/bin/env python3

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

    if (max - min) > 10.5:
        inc = 2.0
    elif (max - min) < 3.0:
        inc = 0.5

    min = math.floor(min / inc) * inc
    max = math.ceil(max / inc) * inc

    return [min, max, inc]

def plot_int(axis, **kwargs):
    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()
    axis.set_axisbelow(True)
    axis.grid(linestyle='dotted')

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

    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()
    axis.set_yscale('log' if ylog else 'linear')
    axis.set_axisbelow(True)
    axis.grid(linestyle='dotted')
    axis.step(x, y, 'royalblue', linewidth=1.0, where='post', label=slice_label)
    axis.fill_between(x, ymin, ymax, alpha=0.4, color='royalblue', linewidth=0.5, step='post')
    axis.set_ylabel(ylabel)

    if slice_label != '':
        axis.legend(fontsize='xx-small', frameon=False)

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

    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()
    axis.set_axisbelow(True)
    axis.grid(linestyle='dotted')
    axis.step(x, qcd_y, 'red', label='NLO QCD', linewidth=1.0, where='post')
    #axis.fill_between(x, qcd_ymin, qcd_ymax, alpha=0.4, color='red', label='7-p.\\ scale var.', linewidth=0.5, step='post')
    axis.step(x, y, 'royalblue', label='NLO QCD+EW', linewidth=1.0, where='post')
    axis.fill_between(x, ymin, ymax, alpha=0.4, color='royalblue', label='7-p.\\ scale var.', linewidth=0.5, step='post')
    axis.errorbar(mid, y[:-1], yerr=(pdf_min, pdf_max), color='royalblue', label='PDF uncertainty', fmt='.', capsize=1, markersize=0, linewidth=1)
    axis.set_ylabel('NLO EW on/off [\\si{{\\percent}}]')
    axis.legend(bbox_to_anchor=(0,1.03,1,0.2), loc='lower left', mode='expand', borderaxespad=0, ncol=4, fontsize='x-small', frameon=False, borderpad=0)

def plot_rel_pdfunc(axis, **kwargs):
    x = kwargs['x']
    pdf_uncertainties = kwargs['pdf_results']
    colors = ['royalblue', 'brown', 'darkorange', 'darkgreen', 'purple', 'tan']

    #ymins = np.asmatrix([(ymin / y - 1.0) * 100 for label, y, ymin, ymax in pdf_uncertainties])
    #ymaxs = np.asmatrix([(ymax / y - 1.0) * 100 for label, y, ymin, ymax in pdf_uncertainties])

    axis.set_axisbelow(True)
    axis.grid(linestyle='dotted')
    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()

    for index, i in enumerate(pdf_uncertainties):
        label, y, ymin, ymax = i
        ymin = (ymin / y - 1.0) * 100.0
        ymax = (ymax / y - 1.0) * 100.0
        axis.step(x, ymax, color=colors[index], label=label, linewidth=1, where='post')
        axis.step(x, ymin, color=colors[index], linewidth=1, where='post')

    #axis.legend(fontsize='xx-small') #rel_pdfunc
    axis.set_ylabel('PDF uncertainty [\\si{{\\percent}}]')

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
    colors = ['royalblue', 'brown', 'darkorange', 'darkgreen', 'purple', 'tan']
    x = kwargs['x']
    y = kwargs['y']

    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()
    axis.set_axisbelow(True)
    axis.grid(linestyle='dotted')

    for index, i in enumerate(pdf_uncertainties):
        label, y, ymin, ymax = i
        diff = y - central_y
        yerr = np.where(diff > 0.0, y - ymin, ymax - y)
        #pull_avg = (y - central_y) / np.sqrt(np.power(0.5 * (ymax - ymin), 2) + np.power(0.5 * (central_ymax - central_ymin), 2))
        pull = (y - central_y) / np.sqrt(np.power(yerr, 2) + np.power(0.5 * (central_ymax - central_ymin), 2))

        #axis.fill_between(x, pull, pull_avg, alpha=0.4, color=colors[index], label='sym.\\ pull', linewidth=0.5, step='post', zorder=2 * index)
        axis.step(x, pull, color=colors[index], label=label, linewidth=1, where='post', zorder=2 * index + 1)

    axis.legend(bbox_to_anchor=(0,1.03,1,0.2), loc='lower left', mode='expand', borderaxespad=0, ncol=len(pdf_uncertainties), fontsize='x-small', frameon=False, borderpad=0) #rel_pdfpull
    axis.set_ylabel('Pull [$\\sigma$]')
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
        plot_abs,
        plot_rel_ewonoff,
    ]

    data_slices = data()

    if len(data_slices[0]['pdf_results']) > 1:
        panels.extend([
            plot_rel_pdfunc,
            plot_rel_pdfpull,
        ])

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\\usepackage{{siunitx}}\\usepackage{{lmodern}}')
    plt.rc('font', family='serif', size=14.0)
    plt.rc('axes', labelsize='small')
    plt.rc('pdf', compression=0)

    xaxis = '{xaxis}'
    xunit = metadata().get(xaxis + '_unit', '')
    xlabel = metadata()[xaxis + '_label_tex'] + (r' [\\si{{' + xunit + r'}}]' if xunit != '' else '')
    ylabel = metadata()['y_label_tex'] + r' [\\si{{' + metadata()['y_unit'] + r'}}]'
    ylog = xunit != ''
    description = metadata()['description']

    if len(data_slices[0]['x']) == 2:
        panels = [ plot_int ]
        xlabel = ylabel
        plt.rc('figure', figsize=(4.2,2.6))
    else:
        plt.rc('figure', figsize=(6.4,len(panels)*2.4))

    for index, dict in enumerate(data_slices):
        dict['xlabel'] = xlabel
        dict['ylabel'] = ylabel
        dict['ylog'] = ylog

        figure, axes = plt.subplots(len(panels), 1, constrained_layout=True, sharex=True, squeeze=False)
        figure.set_constrained_layout_pads(hspace=0, wspace=0)

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
    left = np.array([{left}])
    right = np.array([{right}])
    min = np.array([{min}])
    max = np.array([{max}])
    qcd_central = np.array([{qcd_central}])
    qcd_min = np.array([{qcd_min}])
    qcd_max = np.array([{qcd_max}])
    slices = {slices}
    slice_labels = {slice_labels}
    pdf_results = [
{pdf_results}    ]

    return [{{
        'mid': 0.5 * (left[slice[0]:slice[1]] + right[slice[0]:slice[1]]),
        'pdf_results': [(
            res[0],
            np.append(res[1][slice[0]:slice[1]], res[1][slice[1]-1]),
            np.append(res[2][slice[0]:slice[1]], res[2][slice[1]-1]),
            np.append(res[3][slice[0]:slice[1]], res[3][slice[1]-1])
            ) for res in pdf_results],
        'qcd_max': np.append(qcd_max[slice[0]:slice[1]], qcd_max[slice[1]-1]),
        'qcd_min': np.append(qcd_min[slice[0]:slice[1]], qcd_min[slice[1]-1]),
        'qcd_y': np.append(qcd_central[slice[0]:slice[1]], qcd_central[slice[1]-1]),
        'x': np.append(left[slice[0]:slice[1]], right[slice[1]-1]),
        'y': np.append(pdf_results[0][1][slice[0]:slice[1]], pdf_results[0][1][slice[1]-1]),
        'ymax': np.append(max[slice[0]:slice[1]], max[slice[1]-1]),
        'ymin': np.append(min[slice[0]:slice[1]], min[slice[1]-1]),
        'slice_label': slice_labels[index],
    }} for (index, slice) in enumerate(slices)]

def metadata():
    return {{
{metadata}    }}

if __name__ == '__main__':
    main()",
        xaxis=format!("x{}", bin_info.dimensions()),
        output=output,
        left=map_format_join(left),
        right=map_format_join(right),
        min=map_format_e_join(min),
        max=map_format_e_join(max),
        qcd_central=map_format_e_join(qcd_central),
        qcd_min=map_format_e_join(qcd_min),
        qcd_max=map_format_e_join(qcd_max),
        slices=format!("{:?}", slices),
        slice_labels=format!("[{}]", slice_labels.iter().map(|string| format!("r'{}'", string)).join(", ")),
        pdf_results=format_pdf_results(pdf_uncertainties, pdfsets),
        metadata=format_metadata(metadata),
    );
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        if self.subgrid_pull.is_empty() {
            let grid = helpers::read_grid(&self.input)?;
            let lhapdf_name = pdfset_name(&self.pdfsets[0]);
            let pdf = lhapdf_name.parse().map_or_else(
                |_| Pdf::with_setname_and_member(lhapdf_name, 0),
                Pdf::with_lhaid,
            );

            let results = helpers::convolute(&grid, &pdf, &[], &[], &[], self.scales);

            let qcd_results = {
                let mut orders = grid.orders().to_vec();
                orders.sort();
                let orders = orders;

                let qcd_orders: Vec<_> = orders
                    .iter()
                    .group_by(|order| order.alphas + order.alpha)
                    .into_iter()
                    .map(|mut group| {
                        let order = group.1.next().unwrap();
                        (order.alphas, order.alpha)
                    })
                    .collect();

                helpers::convolute(&grid, &pdf, &qcd_orders, &[], &[], self.scales)
            };

            let bin_info = grid.bin_info();

            let pdf_uncertainties: Vec<Vec<Vec<f64>>> = self
                .pdfsets
                .par_iter()
                .map(|pdfset| {
                    let lhapdf_name = pdfset_name(pdfset);
                    let set = PdfSet::new(&lhapdf_name.parse().map_or_else(
                        |_| lhapdf_name.to_string(),
                        |lhaid| lhapdf::lookup_pdf(lhaid).unwrap().0,
                    ));

                    let pdf_results: Vec<_> = set
                        .mk_pdfs()
                        .into_par_iter()
                        .flat_map(|pdf| helpers::convolute(&grid, &pdf, &[], &[], &[], 1))
                        .collect();

                    let mut central = Vec::with_capacity(bin_info.bins());
                    let mut min = Vec::with_capacity(bin_info.bins());
                    let mut max = Vec::with_capacity(bin_info.bins());

                    for bin in 0..bin_info.bins() {
                        let values: Vec<_> = pdf_results
                            .iter()
                            .skip(bin)
                            .step_by(bin_info.bins())
                            .copied()
                            .collect();

                        let uncertainty = set.uncertainty(&values, helpers::ONE_SIGMA, false);
                        central.push(uncertainty.central);
                        min.push(uncertainty.central - uncertainty.errminus);
                        max.push(uncertainty.central + uncertainty.errplus);
                    }

                    vec![central, min, max]
                })
                .collect();

            let left_limits: Vec<_> = (0..bin_info.dimensions())
                .map(|i| bin_info.left(i))
                .collect();
            let right_limits: Vec<_> = (0..bin_info.dimensions())
                .map(|i| bin_info.right(i))
                .collect();

            let min: Vec<_> = results
                .chunks_exact(self.scales)
                .map(|variations| {
                    variations
                        .iter()
                        .min_by(|a, b| a.partial_cmp(b).unwrap())
                        .copied()
                        .unwrap()
                })
                .collect();
            let max: Vec<_> = results
                .chunks_exact(self.scales)
                .map(|variations| {
                    variations
                        .iter()
                        .max_by(|a, b| a.partial_cmp(b).unwrap())
                        .copied()
                        .unwrap()
                })
                .collect();

            let qcd_central: Vec<_> = qcd_results.iter().step_by(self.scales).copied().collect();
            let qcd_min: Vec<_> = qcd_results
                .chunks_exact(self.scales)
                .map(|variations| {
                    variations
                        .iter()
                        .min_by(|a, b| a.partial_cmp(b).unwrap())
                        .copied()
                        .unwrap()
                })
                .collect();
            let qcd_max: Vec<_> = qcd_results
                .chunks_exact(self.scales)
                .map(|variations| {
                    variations
                        .iter()
                        .max_by(|a, b| a.partial_cmp(b).unwrap())
                        .copied()
                        .unwrap()
                })
                .collect();

            let slices = bin_info.slices();
            let slice_labels: Vec<_> = slices
                .iter()
                .map(|&(begin, end)| {
                    (0..bin_info.dimensions() - 1)
                        .map(|d| {
                            format!(
                                "${} < {} < {}$",
                                bin_info.left(d)[begin],
                                grid.key_values()
                                    .and_then(|map| map
                                        .get(&format!("x{}_label_tex", d + 1))
                                        .cloned())
                                    .unwrap_or_else(|| format!("x{}", d + 1))
                                    .replace("$", ""),
                                bin_info.right(d)[end - 1]
                            )
                        })
                        .join(r#"\\"#)
                })
                .collect();

            let mut key_values = grid.key_values().cloned().unwrap_or_default();
            key_values.entry("description".to_string()).or_default();
            key_values.entry("x1_label_tex".to_string()).or_default();
            key_values.entry("x1_unit".to_string()).or_default();
            key_values.entry("y_label_tex".to_string()).or_default();
            key_values.entry("y_unit".to_string()).or_default();

            let mut vector: Vec<_> = key_values.iter().collect();
            vector.sort();
            let vector = vector;

            let mut output = self.input.clone();

            // remove ".lz4" and ".pineappl" extension
            match output.extension() {
                Some(x) if x == "lz4" => {
                    output = Path::new(output.file_stem().unwrap()).to_path_buf();
                }
                _ => {}
            }
            match output.extension() {
                Some(x) if x == "pineappl" => {
                    output = Path::new(output.file_stem().unwrap()).to_path_buf();
                }
                _ => {}
            }

            format_script(
                &bin_info,
                output.to_str().unwrap(),
                left_limits.last().unwrap(),
                right_limits.last().unwrap(),
                &min,
                &max,
                &qcd_central,
                &qcd_min,
                &qcd_max,
                &slices,
                &slice_labels,
                &pdf_uncertainties,
                &self.pdfsets,
                &vector,
            );
        } else {
            let (pdfset1, pdfset2) = self.pdfsets.iter().collect_tuple().unwrap();
            let (order, bin, lumi) = self
                .subgrid_pull
                .iter()
                .map(|num| num.parse::<usize>().unwrap())
                .collect_tuple()
                .unwrap();

            let cl = helpers::ONE_SIGMA;
            let grid = helpers::read_grid(&self.input)?;

            let set1 = PdfSet::new(&pdfset1.parse().map_or_else(
                |_| pdfset1.to_string(),
                |lhaid| lhapdf::lookup_pdf(lhaid).unwrap().0,
            ));
            let set2 = PdfSet::new(&pdfset2.parse().map_or_else(
                |_| pdfset2.to_string(),
                |lhaid| lhapdf::lookup_pdf(lhaid).unwrap().0,
            ));
            let pdfset1 = set1.mk_pdfs();
            let pdfset2 = set2.mk_pdfs();

            let values1: Vec<f64> = pdfset1
                .par_iter()
                .flat_map(|pdf| helpers::convolute(&grid, pdf, &[], &[bin], &[], 1))
                .collect();
            let values2: Vec<f64> = pdfset2
                .par_iter()
                .flat_map(|pdf| helpers::convolute(&grid, pdf, &[], &[bin], &[], 1))
                .collect();

            let uncertainty1 = set1.uncertainty(&values1, cl, false);
            let uncertainty2 = set2.uncertainty(&values2, cl, false);

            let full_res1 = {
                let central: Vec<f64> = pdfset1
                    .iter()
                    .flat_map(|pdf| helpers::convolute(&grid, pdf, &[], &[], &[], 1))
                    .collect();
                set1.uncertainty(&central, cl, false).central
            };
            let full_res2 = {
                let central: Vec<f64> = pdfset2
                    .iter()
                    .flat_map(|pdf| helpers::convolute(&grid, pdf, &[], &[], &[], 1))
                    .collect();
                set1.uncertainty(&central, cl, false).central
            };

            let res1 = helpers::convolute_subgrid(&grid, &pdfset1[0], order, bin, lumi);
            let res2 = helpers::convolute_subgrid(&grid, &pdfset2[0], order, bin, lumi);

            let denominator = {
                // use the uncertainties in the direction in which the respective results differ
                let unc1 = if full_res1 > full_res2 {
                    uncertainty1.errminus
                } else {
                    uncertainty1.errplus
                };
                let unc2 = if full_res2 > full_res1 {
                    uncertainty2.errminus
                } else {
                    uncertainty2.errplus
                };

                unc1.hypot(unc2)
            };
            let pull = (res2 - res1) / denominator;

            let subgrid = grid.subgrid(order, bin, lumi);
            //let q2 = subgrid.q2_grid();
            let x1 = subgrid.x1_grid();
            let x2 = subgrid.x2_grid();

            let mut x1_vals = vec![];
            let mut x2_vals = vec![];
            let mut vals = vec![];

            for ((_, ix1, ix2), value) in pull
                .indexed_iter()
                .filter(|((_, _, _), value)| **value != 0.0)
            {
                x1_vals.push(x1[ix1]);
                x2_vals.push(x2[ix2]);
                vals.push(*value);
            }

            println!(
                "#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from math import fabs, log10
from scipy.interpolate import griddata

x1 = np.array([{}])
x2 = np.array([{}])
z = np.array([{}])
x = 0.5 * np.log(x1 / x2)
y = np.sqrt(x1 * x2)

nrap = 50
nmass = 50

sym_min = -max(fabs(np.min(x)), fabs(np.max(x)))
sym_max =  max(fabs(np.min(x)), fabs(np.max(x)))

xi = np.linspace(sym_min, sym_max, (nrap // 2) * 2 + 1)
yi = np.logspace(log10(np.min(y)), log10(np.max(y)), nmass)
zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='linear', rescale=True)

#print(xi.shape)
#print(yi.shape)
#print(zi.shape)

# mask impossible kinematic values
for iy, ix in np.ndindex(zi.shape):
    #print(ix, iy)
    x1v = yi[iy] * np.exp(xi[ix])
    x2v = yi[iy] / np.exp(xi[ix])

    #print('y = {{}} m/s = {{}} -> x1 = {{}} x2 = {{}}'.format(xi[ix], yi[iy], x1v, x2v))

    if x1v > 1.0 or x2v > 1.0:
        zi[iy, ix] = np.nan

figure, axes = plt.subplots(1, 2, constrained_layout=True)
figure.set_size_inches(10, 5)

mesh = axes[0].pcolormesh(xi, yi, zi, shading='nearest', linewidth=0, snap=True)
axes[0].scatter(x, y, marker='*', s=5)
axes[0].set_yscale('log')
axes[0].set_xlabel(r'$y = 1/2 \\log (x_1/x_2)$')
axes[0].set_ylabel(r'$M/\\sqrt{{s}} = \\sqrt{{x_1 x_2}}$')
#axes[0].set_aspect('equal', share=True)

x1i = np.logspace(log10(np.min(x1)), log10(np.max(x1)), 50)
x2i = np.logspace(log10(np.min(x2)), log10(np.max(x2)), 50)
z12i = griddata((x1, x2), z, (x1i[None, :], x2i[:, None]), method='linear', rescale=True)

mesh = axes[1].pcolormesh(x1i, x2i, z12i, shading='nearest', linewidth=0, snap=True)
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[1].scatter(x1, x2, marker='*', s=5)
axes[1].set_aspect('equal', share=True)
axes[1].set_xlabel(r'$x_1$')
axes[1].set_ylabel(r'$x_2$')

figure.colorbar(mesh, ax=axes, extend='min')
figure.savefig('plot.pdf')",
                map_format_e_join(&x1_vals),
                map_format_e_join(&x2_vals),
                map_format_e_join(&vals)
            );
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-plot 
Creates a matplotlib script plotting the contents of the grid

USAGE:
    pineappl plot [OPTIONS] <INPUT> [--] [PDFSETS]...

ARGS:
    <INPUT>         Path to the input grid
    <PDFSETS>...    LHAPDF id(s) or name of the PDF set(s)

OPTIONS:
    -h, --help
            Print help information

    -s, --scales <SCALES>
            Set the number of scale variations [default: 7] [possible values: 1, 3, 7, 9]

        --subgrid-pull <ORDER> <BIN> <LUMI>
            Show the pull for a specific grid three-dimensionally
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["plot", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }
}
