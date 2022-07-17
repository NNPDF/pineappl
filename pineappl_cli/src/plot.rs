use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use itertools::Itertools;
use ndarray::Axis;
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
    #[clap(default_value_t = num_cpus::get(), long)]
    threads: usize,
    /// Forces negative PDF values to zero.
    #[clap(long = "force-positive")]
    force_positive: bool,
}

fn map_format_join(slice: &[f64]) -> String {
    slice.iter().map(|x| format!("{}", x)).join(", ")
}

fn map_format_e_join(slice: &[f64]) -> String {
    slice.iter().map(|x| format!("{:.7e}", x)).join(", ")
}

fn format_pdf_results(pdf_uncertainties: &[Vec<Vec<f64>>], pdfsets: &[String]) -> String {
    pdf_uncertainties
        .iter()
        .zip(pdfsets.iter())
        .map(|(values, pdfset)| {
            format!(
                "        (
            '{}',
            np.array([{}]),
            np.array([{}]),
            np.array([{}]),
        ),",
                helpers::pdf_label(pdfset).replace('_', r#"\_"#),
                map_format_e_join(&values[0]),
                map_format_e_join(&values[1]),
                map_format_e_join(&values[2]),
            )
        })
        .join("\n")
}

fn format_metadata(metadata: &[(&String, &String)]) -> String {
    metadata
        .iter()
        .filter_map(|(key, value)| {
            if value.contains('\n') {
                // skip multi-line entries
                None
            } else {
                Some(format!(
                    "        '{}': r'{}',",
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
                ))
            }
        })
        .join("\n")
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

        if self.subgrid_pull.is_empty() {
            let grid = helpers::read_grid(&self.input)?;
            let mut pdf = helpers::create_pdf(&self.pdfsets[0])?;

            let results = helpers::convolute(
                &grid,
                &mut pdf,
                &[],
                &[],
                &[],
                self.scales,
                false,
                self.force_positive,
            );

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

                helpers::convolute(
                    &grid,
                    &mut pdf,
                    &qcd_orders,
                    &[],
                    &[],
                    self.scales,
                    false,
                    self.force_positive,
                )
            };

            let bin_info = grid.bin_info();

            let pdf_uncertainties: Vec<Vec<Vec<f64>>> = self
                .pdfsets
                .par_iter()
                .map(|pdfset| {
                    let (set, member) = helpers::create_pdfset(pdfset).unwrap();

                    let pdf_results: Vec<_> = set
                        .mk_pdfs()
                        .into_par_iter()
                        .flat_map(|mut pdf| {
                            helpers::convolute(
                                &grid,
                                &mut pdf,
                                &[],
                                &[],
                                &[],
                                1,
                                false,
                                self.force_positive,
                            )
                        })
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

                        let uncertainty =
                            set.uncertainty(&values, lhapdf::CL_1_SIGMA, false).unwrap();
                        central.push(if let Some(member) = member {
                            values[member]
                        } else {
                            uncertainty.central
                        });
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

            print!(
                include_str!("plot.py"),
                xaxis = format!("x{}", bin_info.dimensions()),
                output = output.to_str().unwrap(),
                left = map_format_join(left_limits.last().unwrap()),
                right = map_format_join(right_limits.last().unwrap()),
                min = map_format_e_join(&min),
                max = map_format_e_join(&max),
                qcd_central = map_format_e_join(&qcd_central),
                qcd_min = map_format_e_join(&qcd_min),
                qcd_max = map_format_e_join(&qcd_max),
                slices = format!("{:?}", slices),
                slice_labels = format!(
                    "[{}]",
                    slice_labels
                        .iter()
                        .map(|string| format!("r'{}'", string))
                        .join(", ")
                ),
                pdf_results = format_pdf_results(&pdf_uncertainties, &self.pdfsets),
                metadata = format_metadata(&vector),
            );
        } else {
            let (pdfset1, pdfset2) = self.pdfsets.iter().collect_tuple().unwrap();
            let (order, bin, lumi) = self
                .subgrid_pull
                .iter()
                .map(|num| num.parse::<usize>().unwrap())
                .collect_tuple()
                .unwrap();

            let cl = lhapdf::CL_1_SIGMA;
            let grid = helpers::read_grid(&self.input)?;

            let (set1, member1) = helpers::create_pdfset(pdfset1)?;
            let (set2, member2) = helpers::create_pdfset(pdfset2)?;
            let mut pdfset1 = set1.mk_pdfs();
            let mut pdfset2 = set2.mk_pdfs();

            let values1: Vec<f64> = pdfset1
                .par_iter_mut()
                .map(|pdf| {
                    helpers::convolute(&grid, pdf, &[], &[bin], &[], 1, false, self.force_positive)
                        [0]
                })
                .collect();
            let values2: Vec<f64> = pdfset2
                .par_iter_mut()
                .map(|pdf| {
                    helpers::convolute(&grid, pdf, &[], &[bin], &[], 1, false, self.force_positive)
                        [0]
                })
                .collect();

            let uncertainty1 = set1.uncertainty(&values1, cl, false)?;
            let uncertainty2 = set2.uncertainty(&values2, cl, false)?;
            let central1 = if let Some(member1) = member1 {
                values1[member1]
            } else {
                uncertainty1.central
            };
            let central2 = if let Some(member2) = member2 {
                values2[member2]
            } else {
                uncertainty2.central
            };

            let denominator = {
                // use the uncertainties in the direction in which the respective results differ
                let unc1 = if central1 > central2 {
                    uncertainty1.errminus
                } else {
                    uncertainty1.errplus
                };
                let unc2 = if central2 > central1 {
                    uncertainty2.errminus
                } else {
                    uncertainty2.errplus
                };

                unc1.hypot(unc2)
            };

            let res1 = helpers::convolute_subgrid(&grid, &mut pdfset1[0], order, bin, lumi)
                .sum_axis(Axis(0));
            let res2 = helpers::convolute_subgrid(&grid, &mut pdfset2[0], order, bin, lumi)
                .sum_axis(Axis(0));

            let subgrid = grid.subgrid(order, bin, lumi);
            //let q2 = subgrid.q2_grid();
            let x1 = subgrid.x1_grid();
            let x2 = subgrid.x2_grid();

            let mut x1_vals = vec![];
            let mut x2_vals = vec![];
            let mut vals = vec![];

            for (((ix1, ix2), &one), &two) in res1.indexed_iter().zip(res2.iter()) {
                if one == 0.0 {
                    assert_eq!(two, 0.0);
                    continue;
                }

                x1_vals.push(x1[ix1]);
                x2_vals.push(x2[ix2]);
                vals.push((two - one) / denominator);
            }

            print!(
                include_str!("subgrid-pull-plot.py"),
                x1 = map_format_e_join(&x1_vals),
                x2 = map_format_e_join(&x2_vals),
                z = map_format_e_join(&vals)
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
        --force-positive
            Forces negative PDF values to zero

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

def plot_abs_pdfs(axis, **kwargs):
    x = kwargs['x']
    ylog = kwargs['ylog']
    ylabel = kwargs['ylabel']
    slice_label = kwargs['slice_label']
    pdf_uncertainties = kwargs['pdf_results']

    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()
    axis.set_yscale('log' if ylog else 'linear')
    axis.set_axisbelow(True)
    axis.grid(linestyle='dotted')
    axis.set_ylabel(ylabel)

    colors = ['royalblue', 'brown', 'darkorange', 'darkgreen', 'purple', 'tan']
    for index, i in enumerate(pdf_uncertainties):
        label, y, ymin, ymax = i
        axis.step(x, y, color=colors[index], linewidth=1.0, where='post')
        axis.fill_between(x, ymin, ymax, alpha=0.4, color=colors[index], label=label, linewidth=0.5, step='post')

    axis.legend(bbox_to_anchor=(0,-0.24,1,0.2), loc='upper left', mode='expand', borderaxespad=0, ncol=len(pdf_uncertainties), fontsize='x-small', frameon=False, borderpad=0)

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
        ymin = percent_diff(ymin, y)
        ymax = percent_diff(ymax, y)
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
        #plot_abs_pdfs,
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
    min = np.array([3.6114679e2, 3.3214336e2, 2.8866119e2, 2.3342681e2, 1.7416314e2, 1.1835555e2, 5.5752518e1, 1.3296051e1])
    max = np.array([3.8546011e2, 3.5487551e2, 3.0860377e2, 2.4965490e2, 1.8627534e2, 1.2657016e2, 5.9567473e1, 1.4165115e1])
    qcd_central = np.array([3.7918224e2, 3.4849530e2, 3.0262287e2, 2.4443923e2, 1.8224461e2, 1.2371575e2, 5.8295643e1, 1.3881516e1])
    qcd_min = np.array([3.6496990e2, 3.3543085e2, 2.9134557e2, 2.3540632e2, 1.7559223e2, 1.1926045e2, 5.6248136e1, 1.3419686e1])
    qcd_max = np.array([3.8932794e2, 3.5803805e2, 3.1104391e2, 2.5132877e2, 1.8741679e2, 1.2724140e2, 5.9945546e1, 1.4257123e1])
    slices = [(0, 8)]
    slice_labels = [r'']
    pdf_results = [
        (
            'NNPDF31\_nlo\_as\_0118\_luxqed',
            np.array([3.7528868e2, 3.4521365e2, 3.0000102e2, 2.4255656e2, 1.8091118e2, 1.2289094e2, 5.7837137e1, 1.3765722e1]),
            np.array([3.7101036e2, 3.4121833e2, 2.9645139e2, 2.3960485e2, 1.7861975e2, 1.2123259e2, 5.6966771e1, 1.3385530e1]),
            np.array([3.7956699e2, 3.4920896e2, 3.0355066e2, 2.4550827e2, 1.8320260e2, 1.2454929e2, 5.8707503e1, 1.4145914e1]),
        ),
        (
            'NNPDF4.0',
            np.array([3.9213747e2, 3.6039984e2, 3.1252694e2, 2.5175698e2, 1.8696447e2, 1.2647235e2, 5.9461315e1, 1.4510917e1]),
            np.array([3.9016816e2, 3.5855470e2, 3.1086142e2, 2.5032095e2, 1.8577848e2, 1.2551944e2, 5.8778139e1, 1.4032527e1]),
            np.array([3.9410678e2, 3.6224497e2, 3.1419247e2, 2.5319301e2, 1.8815045e2, 1.2742527e2, 6.0144492e1, 1.4989306e1]),
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

x1 = np.array([1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 1.0000000e0, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 9.3094408e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 8.6278393e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.9562425e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 7.2958684e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.6481395e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 6.0147220e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 5.3975723e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.7989890e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 4.2216678e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.6687532e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 3.1438740e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.6511370e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 2.1950413e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.7802566e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.4112081e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 1.0914376e-1, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 8.2281221e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 6.0480029e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 4.3414917e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 3.0521584e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 2.1089187e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 1.4375069e-2, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 9.6991596e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 6.4962062e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 4.3285006e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 2.8738676e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.9034634e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 1.2586797e-3, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 8.3140688e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 5.4877953e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 3.6205450e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 2.3878783e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4, 1.5745606e-4])
x2 = np.array([1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.5745606e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 2.3878783e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 5.4877953e-4, 3.6205450e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 9.6991596e-3, 6.4962062e-3, 4.3285006e-3, 2.8738676e-3, 1.9034634e-3, 1.2586797e-3, 8.3140688e-4, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.4375069e-2, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 3.0521584e-2, 2.1089187e-2, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1, 1.0914376e-1, 8.2281221e-2, 6.0480029e-2, 4.3414917e-2, 1.0000000e0, 9.3094408e-1, 8.6278393e-1, 7.9562425e-1, 7.2958684e-1, 6.6481395e-1, 6.0147220e-1, 5.3975723e-1, 4.7989890e-1, 4.2216678e-1, 3.6687532e-1, 3.1438740e-1, 2.6511370e-1, 2.1950413e-1, 1.7802566e-1, 1.4112081e-1])
z = np.array([-5.0188382e-51, -2.9573096e-32, 3.7712483e-35, 1.8835659e-33, -2.6480798e-32, -3.1934464e-33, 6.9880997e-34, 1.5874498e-35, -6.8303983e-35, 4.6316658e-34, -1.3288199e-34, 4.6837448e-35, -3.7315155e-31, 4.8791005e-30, -5.4904856e-30, -8.0878107e-29, -9.4610861e-28, -1.7510573e-27, 8.6575700e-28, -2.7213960e-26, -7.3003149e-26, -1.6370636e-25, -1.1210926e-24, -1.6844107e-24, -2.3690407e-24, -3.1907249e-24, -4.4501165e-24, -8.8205775e-24, -7.3180228e-24, 4.4147422e-25, 8.0554179e-27, -5.1410719e-34, -7.1731125e-16, 1.7698075e-18, 5.0797173e-16, -1.0689970e-14, -1.0571996e-15, 4.1969793e-17, -4.6192947e-17, 5.3087542e-16, -1.9104928e-16, -5.3173830e-15, 3.8683956e-13, 3.4945006e-14, 1.1321504e-15, -1.1054464e-15, 4.3544937e-13, -3.3278909e-12, -1.3936018e-11, 1.3951783e-11, -1.6951043e-10, -2.4445476e-10, 1.2613268e-11, -4.2461417e-9, -6.5961894e-9, -1.1401025e-8, -8.5701309e-8, -1.1330879e-7, -1.8306772e-7, -1.9013471e-7, -9.2655137e-8, -4.7357297e-7, -2.8496754e-7, 9.2305711e-8, -7.1538252e-9, 1.1157112e-34, 2.6341777e-16, 1.3462569e-18, -1.0556462e-16, 2.2057023e-15, 2.1860410e-16, -1.6246841e-16, 2.9771229e-16, -3.3435553e-15, -6.6279634e-15, 3.0589818e-14, -2.5456627e-12, -2.3559260e-13, 3.5651040e-15, 2.9737882e-14, 2.3130453e-13, 3.9924018e-12, -3.6680769e-11, -1.4797707e-10, -2.4868667e-10, -2.0654258e-10, -1.0938345e-8, -2.2723272e-8, -9.2096533e-8, -1.4194556e-7, -3.6466022e-7, -6.4871551e-7, -1.0711607e-6, -1.3377825e-6, -3.0231484e-6, -3.8631538e-6, -4.5096859e-6, -3.9245184e-7, 6.2562784e-8, 6.8683635e-34, -5.7253709e-17, -1.7026167e-17, 2.0704448e-17, -4.0991178e-16, -3.8991862e-17, 1.2285885e-17, 2.9114574e-15, -3.1693363e-14, -7.3393844e-15, 5.7857867e-14, -3.5737560e-12, -3.8974504e-13, -1.3800121e-13, -3.4449000e-13, 1.0046762e-12, 7.8915615e-12, -1.3070865e-10, -6.5444928e-10, -3.8983166e-10, 5.6392065e-10, -3.4553164e-8, -5.5263184e-8, -1.5730986e-7, -3.3354873e-7, -8.9229266e-7, -1.5348433e-6, -2.0437575e-6, -3.6041927e-6, -5.2377282e-6, -1.0080840e-5, -1.5394815e-5, -2.9338255e-7, 1.2508538e-7, 1.0299714e-33, -6.7451319e-18, -2.3702440e-17, 1.3518454e-18, -6.7364064e-17, 6.8730141e-16, 1.4624134e-16, -2.4720987e-15, 1.1049976e-13, 2.1016690e-16, -1.4055347e-13, 1.0596667e-12, -5.4385867e-12, -4.0502443e-11, 3.7778363e-12, -3.1478579e-10, -1.6728215e-9, -1.1640442e-8, -3.8556768e-8, -1.4369520e-7, -2.6770750e-7, -5.5819490e-7, -1.6566678e-6, -2.9051155e-6, -3.2501525e-6, -7.0515895e-6, -1.6217340e-5, -2.3840462e-5, -1.7590606e-5, 4.7476766e-7, 5.5800598e-8, -1.1924890e-34, 7.5998866e-19, 2.6932328e-18, -1.5394022e-19, 1.7787947e-17, -2.2556836e-14, 4.1432889e-14, -5.3083987e-13, -1.8973377e-11, -2.1007697e-11, 2.4587770e-10, -3.0296970e-9, -2.7356558e-9, -8.1899322e-9, -5.4348048e-8, -1.9185199e-7, -3.0341201e-7, -9.9155830e-7, -2.4759215e-6, -4.9183737e-6, -8.7806886e-6, -1.4140030e-5, -2.6729964e-5, -4.6418489e-5, -4.2330365e-5, -3.8387814e-7, 3.2792063e-7, -9.8629699e-12, 3.5286301e-11, 5.5562108e-12, 2.7378816e-12, -1.2255404e-15, -1.9385283e-14, -1.0778274e-12, -1.7958040e-11, -5.9074501e-11, -1.2787351e-10, -1.1506066e-9, -2.4029484e-9, -4.0743458e-9, -4.4223337e-8, -1.1447605e-7, -3.2264401e-7, -1.0146803e-6, -2.7049841e-6, -4.5127330e-6, -8.7223528e-6, -1.6203566e-5, -2.8922906e-5, -4.3903499e-5, -7.4102851e-5, -3.6391489e-7, 4.6035155e-7, -8.4924347e-14, 4.1800143e-12, 2.0398948e-13, -1.8843777e-14, 6.1116597e-11, -2.0069653e-10, 1.1474498e-10, -2.7571605e-11, -1.2994366e-11, 1.1374509e-12, -4.7773432e-12, -4.7088097e-12, -1.5990683e-11, 2.9202684e-11, -3.6046966e-11, -1.4655357e-9, -8.6259245e-9, -2.0535119e-8, -9.2354286e-8, -2.4793123e-7, -6.4578504e-7, -1.4712474e-6, -3.1644537e-6, -6.6167053e-6, -1.2351922e-5, -2.5023983e-5, -5.7455854e-5, -7.0116570e-5, 4.9897067e-7, 4.0370676e-7, 5.4902131e-13, -2.6988470e-11, -1.3151832e-12, 1.2096536e-13, 3.0616038e-10, -9.2468835e-10, 1.2259243e-9, -2.4422522e-10, 1.8301606e-10, 2.5003291e-10, 2.1653469e-12, 2.9113387e-11, -1.8538205e-10, -6.8960754e-10, -2.7640463e-11, -1.2031708e-11, 1.6459853e-9, 3.2352259e-9, 2.7090289e-8, 1.6063885e-7, 6.0010567e-7, 1.7780935e-6, 3.9110420e-6, 6.8816636e-6, 1.1521345e-5, 2.0415205e-5, 5.9916899e-5, 1.8903339e-5, -4.3812146e-6, 2.0716609e-7, 5.7941511e-12, -2.8439970e-10, -1.3847393e-11, 1.2718931e-12, -1.3809163e-11, 4.0091736e-11, -6.7944422e-11, -3.5995381e-11, 2.8777381e-10, 4.8426665e-10, -1.1353093e-11, 5.3141984e-14, 8.6383223e-12, 3.5952341e-10, -2.4394299e-10, 4.4135632e-9, 1.5971428e-8, 2.9266365e-8, 1.6016875e-7, 6.3760798e-7, 2.3002179e-6, 6.3612979e-6, 1.4636925e-5, 2.9959946e-5, 5.7087538e-5, 1.5271312e-4, 8.2561109e-5, 1.5140468e-3, 2.6161434e-4, -3.3072939e-5, -1.5333300e-13, 7.5164625e-12, 3.6578125e-13, -1.7894833e-12, 1.0040257e-11, 1.7497932e-12, -1.3947335e-13, 4.7717661e-12, -3.0095649e-11, -5.3962743e-11, -1.3718318e-13, 8.8680405e-12, -8.3501943e-12, 1.0241235e-10, 5.6196191e-9, 1.0088183e-8, 2.7258385e-8, 3.5772884e-8, 2.4693981e-7, 1.2224388e-6, 4.3008477e-6, 1.2995427e-5, 3.1432037e-5, 6.7568828e-5, 1.4735914e-4, -1.4227890e-4, 3.4082091e-3, 2.6271372e-2, -1.4314592e-4, -7.7148038e-5, 1.0414176e-11, -5.9460643e-11, -1.0336548e-11, 8.2060069e-13, -1.5430799e-13, 1.8392314e-12, 1.7393608e-11, -7.3783076e-11, -1.9238820e-10, 3.7118291e-10, 2.9342877e-8, 6.6123673e-9, 2.6861864e-8, 1.0134510e-8, 4.0185756e-7, 1.9945245e-6, 6.8879119e-6, 2.2043088e-5, 6.2580410e-5, 1.4652709e-4, 4.0378755e-4, -5.6847270e-3, 5.0529944e-2, 7.1062984e-2, -5.5017458e-3, -4.5378212e-6, -2.0537347e-12, 1.2159507e-10, -6.3712205e-10, -1.3310051e-10, 1.7101899e-10, 2.5104071e-11, 3.6466831e-11, -1.3312631e-9, -3.7033172e-10, 3.9690026e-11, -9.2982551e-11, 5.5974081e-10, 6.7610588e-10, 2.0796585e-9, 9.0511423e-9, -5.0758313e-10, 6.4379251e-9, 1.1495597e-7, 4.4109665e-7, 2.3724571e-6, 9.6926495e-6, 3.3866552e-5, 1.0214555e-4, 2.7023876e-4, 7.3976979e-4, -8.0504383e-3, 1.7536144e-1, 4.8912587e-2, -5.7352342e-3, -1.1702722e-6, -3.4814041e-28, 5.6765910e-12, 2.5338586e-12, -3.0921582e-11, 3.3707645e-10, 4.0372690e-10, -1.1072861e-10, 3.4374964e-10, -1.4684883e-9, -1.0968576e-10, -3.8015668e-10, 1.2453773e-8, 3.4522638e-9, -3.7764716e-10, -6.7054263e-10, 3.3255408e-9, 1.0093457e-8, 2.4597209e-8, 2.9057061e-8, -2.8831059e-9, -1.1392888e-8, 1.2141507e-7, 1.6417592e-7, 1.8742102e-6, 1.1979425e-5, 4.4944919e-5, 1.5288606e-4, 6.1406886e-4, -1.1592011e-2, 8.6784868e-2, 2.3874527e-1, -5.8664630e-3, -7.2286649e-4, -3.4180808e-7, 1.2893739e-26, -1.9205076e-10, 3.7731153e-11, 2.6370857e-10, -1.8120222e-9, -3.8850367e-9, 6.5490499e-9, -6.9845550e-10, -8.6449940e-10, -9.3015450e-10, -5.2815457e-11, 6.8829081e-9, 1.2021427e-9, -3.4820511e-9, -5.2942758e-10, -1.1159193e-12, -2.4217850e-9, 4.0244018e-8, 3.9602395e-9, 1.1286753e-8, 3.8244213e-8, -4.3713311e-8, 6.8362521e-7, 9.1261272e-7, 1.0974677e-5, 4.9547269e-5, 2.0176285e-4, 2.2616057e-4, -1.3948119e-2, 2.8704581e-1, 1.0769533e-1, -9.9538861e-3, -1.3047760e-4, -8.7416842e-9, 2.4251722e-26, -1.7252212e-10, -1.0326747e-9, 2.9083655e-10, -4.4627013e-9, -1.5751659e-9, -9.2751393e-8, -2.1772928e-8, -4.5273750e-9, -9.1077751e-11, -4.6308039e-10, 1.9788778e-10, 4.8206847e-9, 2.0517566e-8, 2.8672819e-9, -5.5830728e-10, 1.0689025e-8, 3.1892627e-8, 4.0093369e-8, 4.5617822e-8, 1.0612165e-7, 4.0420743e-8, 1.8642516e-6, 1.9128145e-6, 1.0567297e-5, 4.5269941e-5, 5.6741724e-4, -1.6674793e-2, 1.2837424e-1, 2.6742217e-1, -4.9601174e-3, -9.8007127e-5, -8.1680346e-5, 2.8212851e-9, 6.7796584e-27, -9.4383700e-11, -1.2209011e-9, -1.6733252e-9, 4.2912032e-10, -5.4302343e-10, -1.5775143e-8, -7.0231130e-9, 1.3426480e-10, -8.8079290e-9, 1.8013058e-8, 3.2477488e-9, -1.2645253e-8, 2.6013611e-8, 4.3683551e-9, 3.0609739e-10, 7.6173491e-9, -1.5122153e-8, -1.4357332e-8, -1.3992860e-9, 3.6017910e-8, 1.9016852e-7, 6.7453639e-8, 2.2982812e-6, 1.2191817e-5, 7.8488489e-5, -3.2500324e-3, 1.2185071e-2, 2.4499643e-1, 5.4845066e-2, -4.0508183e-3, 2.2125283e-4, -3.7365442e-5, 6.9046257e-12, -3.0387942e-26, 1.0406053e-12, -4.3315695e-10, -2.4014558e-9, -2.5337710e-9, -7.8605021e-9, -6.7517168e-8, -2.1461027e-8, -3.7195900e-9, -6.4278622e-9, 1.4799125e-8, 1.0873806e-8, 9.9689497e-8, 1.4633298e-8, -6.8108854e-9, -8.5801097e-9, 3.9372435e-8, 4.1765518e-7, 3.3170576e-8, -6.5764120e-8, 4.2881052e-8, -4.9125889e-9, 2.6119193e-7, 1.2095408e-6, 2.6198872e-6, -1.5288958e-4, -7.0671073e-3, 1.0203324e-1, 9.0524117e-2, -6.0278935e-3, 1.5116762e-3, -4.4272384e-6, -9.6717086e-6, 7.9808792e-26, -6.1271158e-10, -3.6647919e-10, -9.4992351e-9, -1.9700148e-8, -1.1524111e-8, -4.2103886e-8, -2.8352276e-8, -3.5789753e-8, -2.7992823e-9, 2.4866534e-8, 2.7595933e-8, 3.9261927e-7, 1.5416558e-7, 1.1019025e-7, 3.5276327e-8, 1.0203099e-7, 5.1885762e-7, 6.6356026e-7, 2.7481314e-7, 4.5126313e-7, 3.3898538e-7, -1.0356727e-8, 3.0376191e-7, 1.1609050e-5, -2.2208483e-3, 2.0109635e-2, 5.2994590e-2, -7.5835979e-4, 9.2473928e-4, 3.8094657e-4, -3.3851983e-5, -2.6362246e-8, 1.3550146e-25, -5.5986283e-10, -3.7202465e-9, -2.9653882e-9, -5.4883696e-8, -7.4226484e-8, -2.0861394e-8, -5.5387499e-8, -1.0551035e-8, -7.4506964e-10, 7.5699909e-8, 1.4956626e-7, 2.3498258e-7, 1.0445545e-7, 8.6916534e-8, 5.5590164e-8, 1.6770919e-8, -7.1935436e-8, 6.5493149e-7, 7.6260759e-7, -4.9011173e-7, 4.7269808e-7, -5.6684641e-8, 9.8586603e-7, -2.9154100e-4, 2.0913057e-3, 1.5636253e-2, 1.3832930e-3, -2.0701076e-4, 4.2112611e-4, 4.2487513e-6, -4.3393236e-6, 7.0855678e-8, 9.4938716e-27, -1.3641228e-9, -1.1017472e-8, -1.0304826e-8, -6.4977125e-8, -1.2691719e-7, -3.3347298e-8, -3.8973828e-8, -5.2678805e-9, 4.2721740e-8, 1.3860117e-7, 2.8950245e-7, 4.5094004e-7, 4.2325200e-7, 1.7060844e-7, 1.5126132e-7, 5.1282323e-7, 2.1439385e-6, 1.2114479e-6, 3.9416184e-6, 9.3644052e-6, 4.7260476e-7, -1.1655407e-7, -2.1370173e-5, 3.3781981e-5, 2.8140571e-3, 7.1198518e-4, -4.1060450e-4, 1.2024312e-4, -8.1739287e-6, -5.6189657e-6, 5.0305664e-7, 8.3519875e-26, -5.2003915e-10, -4.3625101e-9, -7.6904317e-8, -9.9414343e-8, -1.1787703e-7, -8.0339728e-8, -1.4914140e-8, 7.3880177e-8, 1.3383926e-7, 2.5614897e-7, 5.7857786e-7, 9.6449809e-7, 8.7967197e-7, 5.5653890e-7, 9.2753097e-7, 2.6489553e-6, 4.3138194e-6, 3.6581956e-6, 1.2004222e-5, 3.0626400e-5, 2.0129559e-7, -2.7271515e-7, -9.2687957e-6, 1.0006300e-4, 3.1966378e-5, -1.6036690e-4, 5.5616513e-7, -1.0402357e-5, -2.6657551e-6, 2.2746274e-7, 1.8155672e-9, 3.5061309e-26, 1.5639937e-10, -4.4994217e-9, -8.2343143e-8, -1.7733785e-7, -2.2452893e-7, -8.7662044e-8, 2.9917440e-9, 2.5538281e-7, 2.6085978e-7, 5.9821316e-7, 1.1230241e-6, 1.4294592e-6, 1.7774870e-6, 1.2138088e-6, 1.3894369e-6, 3.9572828e-6, 1.3250052e-5, 3.1591646e-5, 5.4190418e-5, 6.8894920e-5, 3.8765332e-5, -3.9571812e-6, 1.0467347e-6, -9.3497631e-6, -2.3001630e-5, -4.5509198e-7, -3.1709700e-7, -3.7143639e-8, 7.6756206e-9, 2.3124036e-25, -6.0061303e-10, -8.7658480e-9, -7.0072252e-8, -1.2661246e-7, -3.7322628e-7, -1.0766638e-7, 7.2548147e-8, 3.7984648e-7, 1.0018334e-6, 9.8157912e-7, 2.0444036e-6, 2.8600083e-6, 2.6000782e-6, 2.2576577e-6, 3.3681577e-6, 1.1404923e-5, 3.3196589e-5, 1.1272566e-4, 3.6908971e-4, -1.8936773e-3, -1.0683490e-3, 1.2414076e-4, 1.5261211e-25, -1.4800179e-10, -1.8961251e-8, -5.3413702e-8, -2.6251284e-7, -1.8254204e-7, -3.9449392e-8, 1.8510603e-7, 3.5051641e-7, 1.0508858e-6, 1.7108838e-6, 2.8288668e-6, 3.9408687e-6, 4.8522728e-6, 4.2216575e-6, 7.3718450e-6, 3.0017589e-5, 1.0515164e-4, 5.5047614e-4, -7.3749026e-3, 1.4651285e-2, 1.2301024e-2, -1.3074458e-3, 1.8373433e-25, 5.9000691e-12, -3.9943523e-9, -9.3853625e-8, -2.1820843e-7, -2.2758767e-7, -6.3056507e-8, 2.6380443e-7, 8.7706639e-7, 1.0243346e-6, 2.8669776e-6, 4.6088522e-6, 6.3040647e-6, 7.1924179e-6, 7.2837806e-6, 1.5595913e-5, 9.9930002e-5, 1.1060543e-4, -1.1897350e-2, 6.4903885e-2, 1.6157065e-1, 1.3655228e-2, -2.5305048e-3, 1.2328548e-25, 1.6611853e-10, -1.0854356e-9, -4.6256109e-8, -1.1379026e-8, -5.1104819e-7, -2.8785174e-8, 1.0262720e-7, 9.3528427e-7, 1.3899857e-6, 2.8143794e-6, 6.0194463e-6, 8.4979009e-6, 1.0212257e-5, 1.2910442e-5, 6.2825713e-5, -5.5149462e-4, -5.9456253e-3, 1.0077557e-1, 2.2127679e-1, 4.2113623e-2, -5.8762390e-3, 7.5918015e-5, 1.0724205e-25, -3.1497115e-10, -1.7766186e-8, -2.0790088e-8, -3.6390231e-9, -4.2499613e-7, 3.8559150e-9, 2.5869681e-7, 5.3987288e-7, 1.8180000e-6, 2.8920293e-6, 3.8494590e-6, 1.0014578e-5, 1.6876757e-5, 2.5103173e-5, -4.4537350e-4, 1.7889550e-3, 7.9517330e-2, 1.5850167e-1, 2.0294272e-2, -6.2125336e-3, 2.8614866e-4, -2.2083751e-6, 1.2088096e-25, 1.3924353e-10, -3.0260420e-9, -9.8100920e-9, -9.7776160e-8, -1.3801630e-8, 1.5264390e-8, 4.1060127e-7, 1.1820177e-6, 7.2600868e-7, 2.0529567e-6, 5.3723642e-6, 1.3969122e-5, -4.5007108e-5, -1.4002977e-4, 2.2357475e-3, 2.7318118e-2, 4.6317952e-2, -2.2612739e-3, -2.2297737e-3, 7.1253113e-4, -4.6257804e-5, 3.6521244e-7, 4.5031779e-26, 1.2337155e-10, 1.6681358e-9, -2.8803767e-8, -1.1250868e-7, -9.7719001e-10, 3.5746696e-8, 2.2255407e-7, 2.8898409e-7, 8.8314188e-7, 2.6547841e-6, -1.1760482e-5, -4.5233566e-5, 2.5629905e-4, 1.2467428e-3, 3.6474554e-3, 3.6481226e-3, -3.4866707e-3, 5.3548606e-5, 5.7651295e-4, 8.4534825e-7, -5.3232589e-6, 3.9847484e-8, 9.9915403e-27, -1.9987577e-10, -6.2876499e-9, -1.5015569e-8, -1.6566369e-8, -8.8036635e-8, 1.2354888e-8, 1.8096313e-7, 4.5363743e-7, -1.7316793e-7, 5.2213112e-6, 1.0016976e-4, 4.5088053e-4, 7.2905834e-4, 3.3308350e-4, -1.6019570e-4, -4.1940696e-4, 1.2069355e-4, 1.1062078e-4, -6.4539075e-6, -9.0851995e-6, 7.9693875e-7, -2.7606584e-28, -9.2811831e-13, -8.9422399e-10, 8.0093532e-10, -2.0037122e-8, -5.3558074e-8, 6.1039100e-9, 1.9385841e-7, -5.2850673e-8, 3.7730514e-6, 7.2584942e-5, 1.6133939e-4, 1.0575581e-4, -3.1967694e-5, -3.3566507e-5, -5.0802310e-6, 1.6891418e-6, -4.7969439e-6, -1.4136801e-5, -4.4276111e-6, 4.6486989e-7, 1.1344908e-9, 4.4665689e-28, 1.6467611e-11, 1.0603315e-10, -7.6454783e-11, 1.0501071e-9, 4.5901669e-9, -7.7353130e-10, -6.2889866e-9, -2.2479862e-8, 1.0059125e-6, -5.2944123e-7, -1.3450256e-5, -1.2626270e-5, -1.2661532e-6, -1.6238548e-7, -2.9960294e-7, -7.7565009e-7, -7.7402431e-7, -5.7530654e-8, 1.6620503e-8, -4.5559060e-29, -1.4916858e-12, 6.5262935e-13, -2.1420718e-12, 4.1465420e-11, 9.7250794e-11, -6.0605965e-12, 1.7458360e-10, 4.2463903e-10, -1.2401987e-7, -1.9318755e-7, -5.4777467e-9, -2.7642115e-9, -3.6034450e-9, -9.3394494e-11, 1.9694149e-11])
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
                "NNPDF40_nnlo_as_01180=NNPDF4.0",
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
