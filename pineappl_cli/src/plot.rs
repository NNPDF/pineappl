use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use itertools::Itertools;
use lhapdf::{Pdf, PdfSet};
use pineappl::bin::BinInfo;
use pineappl::subgrid::Subgrid;
use rayon::{prelude::*, ThreadPoolBuilder};
use std::path::{Path, PathBuf};

/// Creates a matplotlib script plotting the contents of the grid.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id(s) or name of the PDF set(s).
    #[clap(required = true, validator = helpers::validate_pdfset)]
    pdfsets: Vec<String>,
    /// Set the number of scale variations.
    #[clap(default_value = "7", long, possible_values = &["1", "3", "7", "9"], short)]
    scales: usize,
    /// Show the pull for a specific grid three-dimensionally.
    #[clap(
        conflicts_with = "scales",
        long = "subgrid-pull",
        number_of_values = 3,
        use_value_delimiter = true,
        value_names = &["ORDER", "BIN", "LUMI"]
    )]
    subgrid_pull: Vec<String>,
    /// Number of threads to utilize.
    #[clap(default_value = &helpers::NUM_CPUS_STRING, long)]
    threads: usize,
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
                value.replace('\u{2013}', "--").replace('\u{2014}', "---")
            } else if key.ends_with("_unit") {
                value
                    .replace("GeV", r#"\giga\electronvolt"#)
                    .replace('/', r#"\per"#)
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
        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

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
                                    .replace('$', ""),
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
                .map(|pdf| helpers::convolute(&grid, pdf, &[], &[bin], &[], 1)[0])
                .collect();
            let values2: Vec<f64> = pdfset2
                .par_iter()
                .map(|pdf| helpers::convolute(&grid, pdf, &[], &[bin], &[], 1)[0])
                .collect();

            let uncertainty1 = set1.uncertainty(&values1, cl, false);
            let uncertainty2 = set2.uncertainty(&values2, cl, false);

            let denominator = {
                // use the uncertainties in the direction in which the respective results differ
                let unc1 = if uncertainty1.central > uncertainty2.central {
                    uncertainty1.errminus
                } else {
                    uncertainty1.errplus
                };
                let unc2 = if uncertainty2.central > uncertainty1.central {
                    uncertainty2.errminus
                } else {
                    uncertainty2.errplus
                };

                unc1.hypot(unc2)
            };

            let res1 = helpers::convolute_subgrid(&grid, &pdfset1[0], order, bin, lumi);
            let res2 = helpers::convolute_subgrid(&grid, &pdfset2[0], order, bin, lumi);

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
    pineappl plot [OPTIONS] <INPUT> <PDFSETS>...

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

        --threads <THREADS>
            Number of threads to utilize";

    const DEFAULT_STR: &str = r#"#!/usr/bin/env python3

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
    #axis.fill_between(x, qcd_ymin, qcd_ymax, alpha=0.4, color='red', label='7-p.\ scale var.', linewidth=0.5, step='post')
    axis.step(x, y, 'royalblue', label='NLO QCD+EW', linewidth=1.0, where='post')
    axis.fill_between(x, ymin, ymax, alpha=0.4, color='royalblue', label='7-p.\ scale var.', linewidth=0.5, step='post')
    axis.errorbar(mid, y[:-1], yerr=(pdf_min, pdf_max), color='royalblue', label='PDF uncertainty', fmt='.', capsize=1, markersize=0, linewidth=1)
    axis.set_ylabel('NLO EW on/off [\si{\percent}]')
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
    axis.set_ylabel('PDF uncertainty [\si{\percent}]')

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

        #axis.fill_between(x, pull, pull_avg, alpha=0.4, color=colors[index], label='sym.\ pull', linewidth=0.5, step='post', zorder=2 * index)
        axis.step(x, pull, color=colors[index], label=label, linewidth=1, where='post', zorder=2 * index + 1)

    axis.legend(bbox_to_anchor=(0,1.03,1,0.2), loc='lower left', mode='expand', borderaxespad=0, ncol=len(pdf_uncertainties), fontsize='x-small', frameon=False, borderpad=0) #rel_pdfpull
    axis.set_ylabel('Pull [$\sigma$]')
    #axis.set_title('Comparison with ' + pdf_uncertainties[0][0], fontdict={'fontsize': 9}, loc='left')

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
    plt.rc('text.latex', preamble=r'\usepackage{siunitx}\usepackage{lmodern}')
    plt.rc('font', family='serif', size=14.0)
    plt.rc('axes', labelsize='small')
    plt.rc('pdf', compression=0)

    xaxis = 'x1'
    xunit = metadata().get(xaxis + '_unit', '')
    xlabel = metadata()[xaxis + '_label_tex'] + (r' [\si{' + xunit + r'}]' if xunit != '' else '')
    ylabel = metadata()['y_label_tex'] + r' [\si{' + metadata()['y_unit'] + r'}]'
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

        name = 'LHCB_WP_7TEV' if len(data_slices) == 1 else 'LHCB_WP_7TEV-{}'.format(index)
        figure.savefig(name + '.pdf')
        plt.close(figure)

def data():
    left = np.array([2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4])
    right = np.array([2.25, 2.5, 2.75, 3, 3.25, 3.5, 4, 4.5])
    min = np.array([3.611467938709435e2, 3.321433579837943e2, 2.8866119381356157e2, 2.3342680720650836e2, 1.7416313884856686e2, 1.1835555363031739e2, 5.575251801078587e1, 1.3296050881006622e1])
    max = np.array([3.854601120849203e2, 3.548755078840096e2, 3.086037749726691e2, 2.4965490074407626e2, 1.862753392297644e2, 1.2657016115164228e2, 5.956747324807059e1, 1.4165115452781908e1])
    qcd_central = np.array([3.7918224378311083e2, 3.4849529946144844e2, 3.026228651137125e2, 2.4443922829761934e2, 1.8224461208925558e2, 1.2371575021693724e2, 5.8295642942236164e1, 1.3881516194271834e1])
    qcd_min = np.array([3.6496989850717705e2, 3.3543085426531985e2, 2.913455677770412e2, 2.354063202175256e2, 1.755922255294208e2, 1.1926044515809292e2, 5.6248135554764474e1, 1.341968550885382e1])
    qcd_max = np.array([3.893279351359303e2, 3.5803804940350403e2, 3.110439147670205e2, 2.5132876834886176e2, 1.8741679143549374e2, 1.2724140168537167e2, 5.9945545852637885e1, 1.4257123070357537e1])
    slices = [(0, 8)]
    slice_labels = [r'']
    pdf_results = [
        (
            'NNPDF31\_nlo\_as\_0118\_luxqed',
            np.array([3.752886757854729e2, 3.452136455706134e2, 3.0000102441981136e2, 2.4255656111049774e2, 1.8091117793992933e2, 1.228909400732191e2, 5.783713680623624e1, 1.3765721979348262e1]),
            np.array([3.710103567171106e2, 3.4121833365098774e2, 2.9645138819637395e2, 2.396048492382022e2, 1.7861975403791104e2, 1.2123259241866388e2, 5.696677077443632e1, 1.3385529796543375e1]),
            np.array([3.795669948538352e2, 3.49208957490239e2, 3.0355066064324876e2, 2.4550827298279327e2, 1.8320260184194763e2, 1.2454928772777431e2, 5.8707502838036156e1, 1.4145914162153149e1]),
        ),
        (
            'NNPDF40\_nnlo\_as\_01180',
            np.array([3.9213747082652213e2, 3.6039983590424896e2, 3.125269437891755e2, 2.5175697995505837e2, 1.869644665720543e2, 1.2647235045796468e2, 5.9461315083507415e1, 1.4510916784706136e1]),
            np.array([3.901681605939014e2, 3.585546978495585e2, 3.1086141855465786e2, 2.5032094730650383e2, 1.8577848001034008e2, 1.255194354996348e2, 5.877813851984764e1, 1.4032527431051715e1]),
            np.array([3.9410678105914286e2, 3.622449739589394e2, 3.1419246902369315e2, 2.531930126036129e2, 1.8815045313376854e2, 1.2742526541629456e2, 6.014449164716719e1, 1.4989306138360558e1]),
        ),
    ]

    return [{
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
    } for (index, slice) in enumerate(slices)]

def metadata():
    return {
        'arxiv': r'1505.07024',
        'description': r'LHCb differential W-boson production cross section at 7 TeV',
        'hepdata': r'10.17182/hepdata.2114.v1/t4',
        'initial_state_1': r'2212',
        'initial_state_2': r'2212',
        'lumi_id_types': r'pdg_mc_ids',
        'mg5amc_repo': r'',
        'mg5amc_revno': r'',
        'nnpdf_id': r'LHCBWZMU7TEV',
        'pineappl_gitversion': r'v0.4.1-36-gdbdb5d0',
        'runcard_gitversion': r'7b42083',
        'x1_label': r'etal',
        'x1_label_tex': r'$\eta_{\bar{\ell}}$',
        'x1_unit': r'',
        'y_label': r'disg/detal',
        'y_label_tex': r'$\frac{\mathrm{d}\sigma}{\mathrm{d}\eta_{\bar{\ell}}}$',
        'y_unit': r'\pico\barn',
    }

if __name__ == '__main__':
    main()
"#;

    const SUBGRID_PULL_STR: &str = r#"#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from math import fabs, log10
from scipy.interpolate import griddata

x1 = np.array([1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 9.309440808717544e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.956242522922756e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.648139482473823e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 6.01472197967335e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.798989029610255e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.651137041582823e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 8.228122126204893e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 9.699159574043398e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 6.496206194633799e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 8.314068836488144e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4, 1.5745605600841445e-4])
x2 = np.array([1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1.5745605600841445e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 2.3878782918561914e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 5.487795323670796e-4, 3.6205449638139736e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 9.699159574043398e-3, 6.496206194633799e-3, 4.328500638820811e-3, 2.8738675812817515e-3, 1.9034634022867384e-3, 1.2586797144272762e-3, 8.314068836488144e-4, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1.4375068581090129e-2, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 3.0521584007828916e-2, 2.108918668378717e-2, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1, 1.0914375746330703e-1, 8.228122126204893e-2, 6.0480028754447364e-2, 4.341491741702269e-2, 1e0, 9.309440808717544e-1, 8.627839323906108e-1, 7.956242522922756e-1, 7.295868442414312e-1, 6.648139482473823e-1, 6.01472197967335e-1, 5.397572337880445e-1, 4.798989029610255e-1, 4.221667753589648e-1, 3.668753186482242e-1, 3.1438740076927585e-1, 2.651137041582823e-1, 2.195041265003886e-1, 1.7802566042569432e-1, 1.4112080644440345e-1])
z = np.array([-5.018838166148215e-51, -2.957309577866648e-32, 3.7712482974483165e-35, 1.8835659342276553e-33, -2.6480798423054084e-32, -3.1934464086652085e-33, 6.988099665262425e-34, 1.5874498022781137e-35, -6.830398282501106e-35, 4.63166576102029e-34, -1.3288199414883992e-34, 4.683744751357632e-35, -3.73151549811543e-31, 4.879100472095931e-30, -5.490485640640213e-30, -8.087810698631949e-29, -9.461086123482861e-28, -1.7510572614725296e-27, 8.657569974806858e-28, -2.7213960430795284e-26, -7.300314866856712e-26, -1.637063586949902e-25, -1.1210925757903354e-24, -1.6844106778617083e-24, -2.3690406507085948e-24, -3.1907249004494596e-24, -4.450116479843497e-24, -8.820577461317653e-24, -7.318022752874339e-24, 4.414742231630382e-25, 8.0554178961985e-27, -5.141071931119487e-34, -7.173112450563573e-16, 1.7698075118816883e-18, 5.079717312538366e-16, -1.0689970332079113e-14, -1.0571995688532556e-15, 4.196979329379824e-17, -4.6192947048495896e-17, 5.308754236665094e-16, -1.9104927598976504e-16, -5.317383026263638e-15, 3.868395576470808e-13, 3.4945005824253814e-14, 1.13215036088257e-15, -1.105446358106347e-15, 4.3544936619052796e-13, -3.327890886698413e-12, -1.3936018474479922e-11, 1.395178257699206e-11, -1.6951043070556916e-10, -2.444547573373743e-10, 1.2613267510980372e-11, -4.246141657557616e-9, -6.596189350990619e-9, -1.1401024579212645e-8, -8.570130913322604e-8, -1.1330879427883359e-7, -1.8306772375497517e-7, -1.9013471347234604e-7, -9.265513685345548e-8, -4.735729651595605e-7, -2.849675372553423e-7, 9.230571135385536e-8, -7.153825190444441e-9, 1.115711244522852e-34, 2.6341776680866336e-16, 1.3462568849313863e-18, -1.055646163254658e-16, 2.2057022783684553e-15, 2.1860410495444098e-16, -1.6246841022248858e-16, 2.977122876379794e-16, -3.343555290035438e-15, -6.627963399741571e-15, 3.058981787624574e-14, -2.5456626787118678e-12, -2.3559259895689177e-13, 3.56510396594493e-15, 2.9737882466143795e-14, 2.313045275892573e-13, 3.992401768854218e-12, -3.668076936731356e-11, -1.4797707418953525e-10, -2.486866746470346e-10, -2.0654258480637168e-10, -1.0938345299339416e-8, -2.2723271708102823e-8, -9.209653300561791e-8, -1.41945555497321e-7, -3.6466022225172464e-7, -6.487155121584557e-7, -1.0711607300641675e-6, -1.3377825006813592e-6, -3.0231484361232626e-6, -3.863153828414478e-6, -4.509685862353067e-6, -3.9245184380590893e-7, 6.256278359465537e-8, 6.86836354918982e-34, -5.72537087730884e-17, -1.7026167101158498e-17, 2.0704447581869987e-17, -4.099117796794332e-16, -3.899186161775195e-17, 1.2285884779194495e-17, 2.9114574006628753e-15, -3.1693363304091675e-14, -7.339384425181058e-15, 5.785786710978605e-14, -3.5737560077565655e-12, -3.89745040334193e-13, -1.3800121350747045e-13, -3.444900034772229e-13, 1.004676154129716e-12, 7.891561517586924e-12, -1.3070865079779677e-10, -6.544492815758342e-10, -3.898316615690138e-10, 5.639206518747429e-10, -3.45531636237901e-8, -5.5263183823032175e-8, -1.5730986465741266e-7, -3.3354872644088843e-7, -8.922926629303293e-7, -1.534843254740636e-6, -2.0437574757755745e-6, -3.6041927144479273e-6, -5.237728199291294e-6, -1.0080840150622417e-5, -1.5394814606262736e-5, -2.9338255187276347e-7, 1.250853833694477e-7, 1.0299714248696377e-33, -6.745131926174216e-18, -2.3702440067122692e-17, 1.3518454276883178e-18, -6.736406432663886e-17, 6.873014132223288e-16, 1.4624134035483066e-16, -2.472098671262197e-15, 1.1049976277748616e-13, 2.1016689754941557e-16, -1.4055346817271921e-13, 1.059666657527821e-12, -5.438586672314008e-12, -4.050244278527531e-11, 3.777836263633641e-12, -3.1478579088525217e-10, -1.6728215283102646e-9, -1.1640442386973236e-8, -3.8556767721321255e-8, -1.4369520461525108e-7, -2.677075016532712e-7, -5.581948968640844e-7, -1.656667776717956e-6, -2.9051154623212537e-6, -3.2501525283035664e-6, -7.05158946006225e-6, -1.6217340175329973e-5, -2.3840461589033345e-5, -1.7590605828739317e-5, 4.747676625672604e-7, 5.580059752639781e-8, -1.1924890231084554e-34, 7.599886595105224e-19, 2.693232771553793e-18, -1.5394021791067784e-19, 1.7787946989721384e-17, -2.2556836091407168e-14, 4.143288945565264e-14, -5.308398657795211e-13, -1.8973376870357858e-11, -2.1007696721209317e-11, 2.458777000772375e-10, -3.029696951699479e-9, -2.735655750248933e-9, -8.189932203928698e-9, -5.434804827063471e-8, -1.918519892368016e-7, -3.0341201331129486e-7, -9.915582976074323e-7, -2.4759215049548077e-6, -4.918373666201417e-6, -8.780688623415049e-6, -1.4140030322103317e-5, -2.672996433545485e-5, -4.6418488869378784e-5, -4.2330365232827586e-5, -3.838781372491623e-7, 3.279206331978769e-7, -9.862969937484105e-12, 3.528630069114042e-11, 5.5562107518374434e-12, 2.7378816089701576e-12, -1.2255403728627612e-15, -1.938528347338268e-14, -1.0778273774056785e-12, -1.795803983294172e-11, -5.907450068925247e-11, -1.2787351244131467e-10, -1.150606634203518e-9, -2.4029484277722083e-9, -4.074345763069802e-9, -4.4223337414772774e-8, -1.1447605343747708e-7, -3.2264401427021975e-7, -1.014680272635428e-6, -2.7049841116195318e-6, -4.512733024337498e-6, -8.722352806047037e-6, -1.6203565951919712e-5, -2.8922906151958995e-5, -4.3903499028187855e-5, -7.41028506978888e-5, -3.6391489151334335e-7, 4.603515526722768e-7, -8.49243471259974e-14, 4.180014327525392e-12, 2.0398948430830268e-13, -1.884377740696879e-14, 6.111659744110149e-11, -2.0069653470907366e-10, 1.1474497851101905e-10, -2.7571604979250395e-11, -1.299436599216529e-11, 1.1374509481712948e-12, -4.777343219032587e-12, -4.708809741347557e-12, -1.599068280112853e-11, 2.92026835532801e-11, -3.604696604858736e-11, -1.4655357392187374e-9, -8.625924454071327e-9, -2.0535119073828162e-8, -9.235428579649278e-8, -2.4793122718658895e-7, -6.457850423470007e-7, -1.4712474230698272e-6, -3.1644537114773347e-6, -6.616705277242722e-6, -1.2351922096651926e-5, -2.5023982637325794e-5, -5.745585350860222e-5, -7.011656981690364e-5, 4.98970673573542e-7, 4.037067620733834e-7, 5.490213103270969e-13, -2.6988469937507485e-11, -1.315183220379259e-12, 1.209653563369192e-13, 3.061603772721206e-10, -9.246883480767586e-10, 1.2259243326343894e-9, -2.442252179182171e-10, 1.8301605858971215e-10, 2.5003291161721384e-10, 2.1653469431937904e-12, 2.911338746743036e-11, -1.8538204804562524e-10, -6.896075443712793e-10, -2.7640462991727877e-11, -1.2031707859907635e-11, 1.6459852801935045e-9, 3.235225899747997e-9, 2.709028948252355e-8, 1.606388516008557e-7, 6.001056748311089e-7, 1.7780935427202307e-6, 3.911041973186779e-6, 6.881663636050533e-6, 1.1521344706392011e-5, 2.041520526646101e-5, 5.991689883180766e-5, 1.8903338506479895e-5, -4.381214629311404e-6, 2.0716608735657412e-7, 5.794151089818041e-12, -2.8439970126777896e-10, -1.3847392551116896e-11, 1.2718931141994175e-12, -1.3809163349881791e-11, 4.0091735724840464e-11, -6.794442230798617e-11, -3.599538100464055e-11, 2.8777381140184004e-10, 4.842666521841323e-10, -1.1353093191745048e-11, 5.314198417792098e-14, 8.638322336539479e-12, 3.595234127752099e-10, -2.439429946950415e-10, 4.413563204785561e-9, 1.5971428032499683e-8, 2.9266365002048813e-8, 1.6016875023858168e-7, 6.376079794890488e-7, 2.3002179416213522e-6, 6.36129785201614e-6, 1.4636924982002493e-5, 2.995994609732531e-5, 5.7087537528332525e-5, 1.5271311866518194e-4, 8.256110869729453e-5, 1.5140468307739028e-3, 2.6161433771059974e-4, -3.307293875477704e-5, -1.5333300323468706e-13, 7.516462526190007e-12, 3.6578125312641803e-13, -1.789483310948729e-12, 1.004025714851595e-11, 1.74979316473824e-12, -1.3947335043961694e-13, 4.7717660987955436e-12, -3.009564922026557e-11, -5.396274268999562e-11, -1.3718318232040564e-13, 8.868040523353385e-12, -8.350194272401248e-12, 1.024123511350953e-10, 5.619619120093371e-9, 1.008818313233683e-8, 2.725838524697115e-8, 3.577288373170981e-8, 2.4693981056613815e-7, 1.2224387983398187e-6, 4.300847707062365e-6, 1.2995426831954321e-5, 3.14320374068159e-5, 6.756882822637709e-5, 1.4735914121720146e-4, -1.4227889616146903e-4, 3.408209133917789e-3, 2.627137215378902e-2, -1.4314591921201202e-4, -7.714803767318045e-5, 1.0414175577642259e-11, -5.946064343814605e-11, -1.0336547667034882e-11, 8.206006941862092e-13, -1.5430799398198898e-13, 1.8392314348443868e-12, 1.739360761169683e-11, -7.378307588795373e-11, -1.9238819515576855e-10, 3.7118290797487615e-10, 2.934287717559674e-8, 6.612367287926759e-9, 2.6861864437180547e-8, 1.0134510212526812e-8, 4.0185755998280525e-7, 1.994524476068378e-6, 6.88791193305066e-6, 2.2043088073391227e-5, 6.258041032369119e-5, 1.4652709103246689e-4, 4.0378754677683894e-4, -5.684727040595183e-3, 5.0529944132180545e-2, 7.106298437914617e-2, -5.501745822042332e-3, -4.5378212174684806e-6, -2.0537347495773122e-12, 1.2159506585528716e-10, -6.371220501920143e-10, -1.3310051354091922e-10, 1.710189948658285e-10, 2.5104070856295715e-11, 3.6466830928278204e-11, -1.3312630698843395e-9, -3.703317177416344e-10, 3.969002633763495e-11, -9.298255061265103e-11, 5.597408055080082e-10, 6.761058806840334e-10, 2.079658466534954e-9, 9.051142300328965e-9, -5.075831255076704e-10, 6.4379251207863494e-9, 1.1495596720266407e-7, 4.4109665176271884e-7, 2.372457148044651e-6, 9.692649498262611e-6, 3.38665519040332e-5, 1.0214555057966928e-4, 2.702387591421342e-4, 7.397697880776805e-4, -8.050438287137841e-3, 1.753614447514083e-1, 4.8912587034112026e-2, -5.7352341856121945e-3, -1.1702722085316354e-6, -3.4814041033806267e-28, 5.6765909792036995e-12, 2.5338585513770496e-12, -3.092158187141543e-11, 3.3707645176994523e-10, 4.0372689745589033e-10, -1.1072861377151162e-10, 3.437496369128016e-10, -1.4684883297425113e-9, -1.0968575905436102e-10, -3.801566842596621e-10, 1.2453772753703908e-8, 3.4522638250098883e-9, -3.7764716476908304e-10, -6.70542634522906e-10, 3.3255407903493127e-9, 1.0093456645106954e-8, 2.4597208805481817e-8, 2.905706117926122e-8, -2.8831058977925754e-9, -1.1392887989090226e-8, 1.214150681847414e-7, 1.64175922931943e-7, 1.8742101594230223e-6, 1.1979424930334026e-5, 4.494491859497816e-5, 1.5288605597992388e-4, 6.14068856789725e-4, -1.1592011154264234e-2, 8.678486818376237e-2, 2.3874526970410406e-1, -5.8664629999846596e-3, -7.228664949219039e-4, -3.4180807689368596e-7, 1.2893739433831964e-26, -1.9205076329074955e-10, 3.773115314183441e-11, 2.6370856872525703e-10, -1.812022244988107e-9, -3.885036682171987e-9, 6.549049857570029e-9, -6.984555042024442e-10, -8.644994015683317e-10, -9.301544967901853e-10, -5.281545669521937e-11, 6.882908113240923e-9, 1.2021427428312268e-9, -3.482051057327829e-9, -5.294275842926427e-10, -1.1159192536645922e-12, -2.421785028963835e-9, 4.024401814628089e-8, 3.960239475764473e-9, 1.1286753341373527e-8, 3.824421263062287e-8, -4.3713311198113124e-8, 6.836252117883252e-7, 9.126127194084807e-7, 1.0974677458598583e-5, 4.9547269375081254e-5, 2.0176284999758493e-4, 2.2616056738072616e-4, -1.3948119408863221e-2, 2.8704580895746556e-1, 1.0769532936169836e-1, -9.953886073133961e-3, -1.3047759888394538e-4, -8.741684226199378e-9, 2.425172168064874e-26, -1.7252212331371797e-10, -1.0326746536571977e-9, 2.908365491168203e-10, -4.46270127423874e-9, -1.5751659441890632e-9, -9.275139336059565e-8, -2.1772927658501554e-8, -4.5273749585176915e-9, -9.10777511465344e-11, -4.630803922104757e-10, 1.9788778063917307e-10, 4.820684730040187e-9, 2.05175662636411e-8, 2.867281861489526e-9, -5.583072784387118e-10, 1.068902474701102e-8, 3.1892626837463883e-8, 4.009336859894829e-8, 4.5617821822428815e-8, 1.061216540835803e-7, 4.042074296371066e-8, 1.8642516250602268e-6, 1.9128145429322353e-6, 1.0567296997776811e-5, 4.526994127601919e-5, 5.674172416290609e-4, -1.6674793163085577e-2, 1.283742375455506e-1, 2.674221660392871e-1, -4.960117426237542e-3, -9.800712725845723e-5, -8.168034564762207e-5, 2.821285090668156e-9, 6.779658355276525e-27, -9.438369990001063e-11, -1.2209010626638226e-9, -1.673325173211601e-9, 4.2912031641701014e-10, -5.430234273425716e-10, -1.5775142588744487e-8, -7.023113036356527e-9, 1.342648017642789e-10, -8.807929046809012e-9, 1.8013058336110607e-8, 3.2477488044631515e-9, -1.264525305989082e-8, 2.6013610747979785e-8, 4.368355129201483e-9, 3.060973949182115e-10, 7.617349088578155e-9, -1.5122153360957876e-8, -1.4357332301822088e-8, -1.399286001898863e-9, 3.6017910498947884e-8, 1.9016852413840154e-7, 6.745363858928725e-8, 2.298281165063641e-6, 1.2191817140714789e-5, 7.848848863476519e-5, -3.2500323694927068e-3, 1.2185070612013564e-2, 2.4499642852652578e-1, 5.484506557793055e-2, -4.050818340171986e-3, 2.2125283388998895e-4, -3.736544246696656e-5, 6.9046257172061805e-12, -3.0387941780147523e-26, 1.040605265323561e-12, -4.331569513672491e-10, -2.4014557635535628e-9, -2.533771011228578e-9, -7.860502107444659e-9, -6.751716784963645e-8, -2.14610265556215e-8, -3.7195899536523965e-9, -6.427862159106522e-9, 1.4799125126440326e-8, 1.0873805730591043e-8, 9.96894966171686e-8, 1.4633298086745423e-8, -6.810885407672439e-9, -8.580109674615359e-9, 3.93724350032013e-8, 4.176551789045155e-7, 3.317057552079629e-8, -6.57641200013288e-8, 4.288105240608641e-8, -4.9125888531049295e-9, 2.6119192721616453e-7, 1.209540843806445e-6, 2.6198871915724277e-6, -1.528895755487261e-4, -7.067107298838879e-3, 1.0203324125891498e-1, 9.052411714537419e-2, -6.02789353685679e-3, 1.5116762363648806e-3, -4.427238396961415e-6, -9.671708566728603e-6, 7.980879227358106e-26, -6.127115806307832e-10, -3.6647919183050775e-10, -9.499235144024349e-9, -1.9700148480924034e-8, -1.1524110596077857e-8, -4.210388640361268e-8, -2.835227585400416e-8, -3.578975286233115e-8, -2.7992823233799e-9, 2.486653420527589e-8, 2.7595932994763675e-8, 3.9261927193909053e-7, 1.5416557970995013e-7, 1.1019025127987583e-7, 3.527632743253891e-8, 1.0203099322721236e-7, 5.188576249537999e-7, 6.63560257139379e-7, 2.748131389580547e-7, 4.5126312710014856e-7, 3.3898538339342057e-7, -1.0356726575853071e-8, 3.03761909348861e-7, 1.160905033556471e-5, -2.2208483102115744e-3, 2.0109635020953878e-2, 5.2994589818360076e-2, -7.583597934065052e-4, 9.247392798912698e-4, 3.809465714776637e-4, -3.385198347435545e-5, -2.6362246308963196e-8, 1.3550145512818229e-25, -5.598628334014186e-10, -3.7202465293929606e-9, -2.9653881744885464e-9, -5.488369580004962e-8, -7.422648357070086e-8, -2.0861393711065022e-8, -5.5387498613078865e-8, -1.0551034551799497e-8, -7.450696434920689e-10, 7.569990909449514e-8, 1.4956625908202595e-7, 2.3498257539271555e-7, 1.0445545108803347e-7, 8.691653401430152e-8, 5.5590163929177963e-8, 1.6770918987998078e-8, -7.193543628303412e-8, 6.549314883881083e-7, 7.626075930140659e-7, -4.901117266424936e-7, 4.726980813166694e-7, -5.6684641323912826e-8, 9.858660267230046e-7, -2.9154099808690985e-4, 2.0913057488105248e-3, 1.563625300290741e-2, 1.3832929538366578e-3, -2.0701075562106844e-4, 4.2112610923227055e-4, 4.2487512953488115e-6, -4.339323602926813e-6, 7.085567834885643e-8, 9.493871557028016e-27, -1.3641228194963804e-9, -1.101747165322307e-8, -1.0304826413604919e-8, -6.49771250226081e-8, -1.2691718545129315e-7, -3.334729804620414e-8, -3.897382808110507e-8, -5.267880451961027e-9, 4.272174025343819e-8, 1.386011701706142e-7, 2.895024471752082e-7, 4.5094003798627835e-7, 4.232519995200862e-7, 1.7060844441992956e-7, 1.5126132028903135e-7, 5.128232300321026e-7, 2.1439385277593617e-6, 1.2114479306022904e-6, 3.941618389730869e-6, 9.364405168959014e-6, 4.7260475803630445e-7, -1.1655407124147913e-7, -2.137017323585693e-5, 3.378198119517787e-5, 2.8140571401639868e-3, 7.119851800250123e-4, -4.1060450179175015e-4, 1.2024312028185835e-4, -8.173928692811483e-6, -5.61896574861069e-6, 5.030566366798519e-7, 8.351987457464607e-26, -5.200391495551506e-10, -4.3625100987468115e-9, -7.690431702834511e-8, -9.941434293749476e-8, -1.1787703067229141e-7, -8.033972826470539e-8, -1.4914139669180736e-8, 7.388017681434658e-8, 1.3383926036038555e-7, 2.56148974875495e-7, 5.785778623741579e-7, 9.644980873757942e-7, 8.796719723637977e-7, 5.565389000513189e-7, 9.275309690527607e-7, 2.6489553178270213e-6, 4.313819445402903e-6, 3.658195605260452e-6, 1.2004221654707706e-5, 3.062640042967593e-5, 2.0129558833929646e-7, -2.727151505844094e-7, -9.268795657007592e-6, 1.0006300261428845e-4, 3.196637830074756e-5, -1.6036690135698143e-4, 5.56165129735267e-7, -1.0402357392347506e-5, -2.665755058130579e-6, 2.2746274090093512e-7, 1.8155672025404362e-9, 3.506130947785025e-26, 1.5639937421804317e-10, -4.499421657725999e-9, -8.234314306556144e-8, -1.7733784999061594e-7, -2.2452892579838102e-7, -8.766204389170522e-8, 2.9917440459821115e-9, 2.553828128252003e-7, 2.6085978241279777e-7, 5.982131571288384e-7, 1.123024069059967e-6, 1.4294591732699922e-6, 1.777487029773055e-6, 1.2138087940019315e-6, 1.3894369419533295e-6, 3.9572827971404105e-6, 1.3250052147694549e-5, 3.159164562257103e-5, 5.4190417766319634e-5, 6.889491966072254e-5, 3.8765331900020094e-5, -3.9571812161004975e-6, 1.0467347429158892e-6, -9.349763109638392e-6, -2.3001630481965364e-5, -4.5509197917765596e-7, -3.1709699831985273e-7, -3.714363918789627e-8, 7.675620570187417e-9, 2.312403592780134e-25, -6.006130259119563e-10, -8.765847970924451e-9, -7.007225168899025e-8, -1.2661245676428014e-7, -3.7322627887121623e-7, -1.0766637991464257e-7, 7.254814686987587e-8, 3.798464752503007e-7, 1.0018334053883591e-6, 9.815791173527973e-7, 2.044403582886256e-6, 2.8600082796402275e-6, 2.6000782108874545e-6, 2.2576576733838936e-6, 3.3681576576640885e-6, 1.14049232577241e-5, 3.319658903573068e-5, 1.1272565872268767e-4, 3.690897052798232e-4, -1.893677323230782e-3, -1.0683490241891246e-3, 1.241407591974753e-4, 1.5261211212566063e-25, -1.4800178784465605e-10, -1.8961251225215176e-8, -5.3413702462724724e-8, -2.6251284475419256e-7, -1.8254204419875617e-7, -3.944939168673051e-8, 1.8510603059758894e-7, 3.5051641436710144e-7, 1.0508858095086446e-6, 1.7108838327995665e-6, 2.8288668480811353e-6, 3.9408686946477515e-6, 4.852272840268788e-6, 4.22165746795227e-6, 7.3718449787581036e-6, 3.0017589011384806e-5, 1.0515163511761869e-4, 5.504761351582114e-4, -7.374902574584406e-3, 1.4651285339763799e-2, 1.2301024404972057e-2, -1.3074457777774858e-3, 1.8373433479988432e-25, 5.900069073048575e-12, -3.99435231184818e-9, -9.385362475471775e-8, -2.1820843240986026e-7, -2.275876737768987e-7, -6.305650666571046e-8, 2.638044260377187e-7, 8.770663857136225e-7, 1.024334614704583e-6, 2.8669775705234265e-6, 4.6088522312105005e-6, 6.304064671345452e-6, 7.192417910054078e-6, 7.283780616011439e-6, 1.55959129567143e-5, 9.993000183781306e-5, 1.1060542687472752e-4, -1.1897350038911837e-2, 6.490388512795407e-2, 1.615706533007414e-1, 1.3655228208353045e-2, -2.5305047861917955e-3, 1.232854848963344e-25, 1.661185262169757e-10, -1.0854356243495676e-9, -4.625610942946871e-8, -1.137902649554727e-8, -5.110481936685874e-7, -2.8785173819181185e-8, 1.0262720044911914e-7, 9.352842709538371e-7, 1.3899856737596215e-6, 2.8143794317664387e-6, 6.019446269282003e-6, 8.497900912551764e-6, 1.0212256868446952e-5, 1.2910441980836956e-5, 6.282571252559875e-5, -5.514946166603693e-4, -5.945625307632748e-3, 1.0077557071015437e-1, 2.2127678524946448e-1, 4.211362302917454e-2, -5.8762389632420736e-3, 7.591801476029807e-5, 1.072420489699318e-25, -3.1497114922105896e-10, -1.7766186246036212e-8, -2.0790088088924018e-8, -3.639023053878176e-9, -4.249961287515429e-7, 3.855915003770035e-9, 2.5869681259776215e-7, 5.398728849502273e-7, 1.8180000297638846e-6, 2.8920292819533257e-6, 3.849458957380584e-6, 1.0014577881525236e-5, 1.6876756670027697e-5, 2.510317266976481e-5, -4.453734998203431e-4, 1.7889549776332542e-3, 7.951733030305724e-2, 1.585016671848634e-1, 2.0294272387070933e-2, -6.212533597852395e-3, 2.8614866021035986e-4, -2.208375105452173e-6, 1.2088095844837676e-25, 1.39243527388719e-10, -3.026041985915312e-9, -9.810092036540814e-9, -9.777616048905039e-8, -1.3801630123897203e-8, 1.52643898830262e-8, 4.106012695271169e-7, 1.1820176640024301e-6, 7.260086829744045e-7, 2.052956656878686e-6, 5.372364159436626e-6, 1.396912240046814e-5, -4.500710811825104e-5, -1.400297680044393e-4, 2.235747529302771e-3, 2.7318118340945625e-2, 4.6317951786437224e-2, -2.261273889509117e-3, -2.229773746687163e-3, 7.125311257881112e-4, -4.6257804360500875e-5, 3.652124360851531e-7, 4.5031778982373847e-26, 1.2337155069073723e-10, 1.668135813833283e-9, -2.8803767246819858e-8, -1.1250868341829361e-7, -9.771900050631459e-10, 3.574669637420336e-8, 2.2255406589084164e-7, 2.889840889818428e-7, 8.831418781359681e-7, 2.6547840716910873e-6, -1.176048235540906e-5, -4.5233566175574616e-5, 2.5629904718869694e-4, 1.2467428103969949e-3, 3.6474554376674257e-3, 3.648122619110954e-3, -3.4866706516359235e-3, 5.354860611156527e-5, 5.765129482389453e-4, 8.453482484885872e-7, -5.323258917912528e-6, 3.9847483703493954e-8, 9.991540290762291e-27, -1.998757745698704e-10, -6.287649894260476e-9, -1.5015569111231075e-8, -1.656636914126563e-8, -8.803663464774603e-8, 1.2354887782786379e-8, 1.8096313470227688e-7, 4.536374335316951e-7, -1.7316792791266625e-7, 5.22131116848862e-6, 1.0016975849445151e-4, 4.508805253969142e-4, 7.290583352119486e-4, 3.330834951601644e-4, -1.6019569906885537e-4, -4.194069560774249e-4, 1.2069355299807102e-4, 1.1062077548610451e-4, -6.453907547099542e-6, -9.085199519296877e-6, 7.969387544930081e-7, -2.7606584377712776e-28, -9.281183139572747e-13, -8.942239862567355e-10, 8.009353161259802e-10, -2.0037122290748796e-8, -5.355807394217178e-8, 6.103909994039976e-9, 1.9385841171513067e-7, -5.2850673295596135e-8, 3.773051385624605e-6, 7.2584942123928e-5, 1.6133938548074849e-4, 1.0575581418998726e-4, -3.196769377689813e-5, -3.356650671219028e-5, -5.080231026949053e-6, 1.6891418081389211e-6, -4.796943942280705e-6, -1.4136800964637854e-5, -4.427611104680495e-6, 4.64869887141185e-7, 1.1344907896086568e-9, 4.466568854731604e-28, 1.6467610951593003e-11, 1.0603315424421535e-10, -7.645478327915404e-11, 1.0501070713763011e-9, 4.590166923911272e-9, -7.735313010743573e-10, -6.2889866349320606e-9, -2.247986159380717e-8, 1.0059125331328974e-6, -5.294412343592494e-7, -1.3450255808935717e-5, -1.2626269912942293e-5, -1.2661532105981344e-6, -1.6238548053646774e-7, -2.996029411148495e-7, -7.756500904588474e-7, -7.740243087458115e-7, -5.7530654152804694e-8, 1.6620502778272613e-8, -4.5559059584310166e-29, -1.4916858086587347e-12, 6.526293536685834e-13, -2.1420717973282973e-12, 4.146542040268893e-11, 9.725079374349085e-11, -6.0605964541827355e-12, 1.7458359715440917e-10, 4.246390289302025e-10, -1.2401987090276176e-7, -1.9318754542517384e-7, -5.4777467051991385e-9, -2.7642115102801804e-9, -3.6034449927654325e-9, -9.339449387359381e-11, 1.9694148654994408e-11])
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

    #print('y = {} m/s = {} -> x1 = {} x2 = {}'.format(xi[ix], yi[iy], x1v, x2v))

    if x1v > 1.0 or x2v > 1.0:
        zi[iy, ix] = np.nan

figure, axes = plt.subplots(1, 2, constrained_layout=True)
figure.set_size_inches(10, 5)

mesh = axes[0].pcolormesh(xi, yi, zi, shading='nearest', linewidth=0, snap=True)
axes[0].scatter(x, y, marker='*', s=5)
axes[0].set_yscale('log')
axes[0].set_xlabel(r'$y = 1/2 \log (x_1/x_2)$')
axes[0].set_ylabel(r'$M/\sqrt{s} = \sqrt{x_1 x_2}$')
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
figure.savefig('plot.pdf')
"#;

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["plot", "--help"])
            .assert()
            .success()
            .stdout(format!("{} [default: {}]\n", HELP_STR, num_cpus::get()));
    }

    #[test]
    fn default() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "plot",
                "--threads=1",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
                "NNPDF40_nnlo_as_01180",
            ])
            .assert()
            .success()
            .stdout(DEFAULT_STR);
    }

    #[test]
    fn subgrid_pull() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "plot",
                "--subgrid-pull=0,0,0",
                "--threads=1",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
                "NNPDF40_nnlo_as_01180",
            ])
            .assert()
            .success()
            .stdout(SUBGRID_PULL_STR);
    }
}
