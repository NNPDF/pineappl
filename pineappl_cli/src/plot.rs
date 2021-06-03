use super::helpers;
use anyhow::Result;
use itertools::Itertools;
use lhapdf::{Pdf, PdfSet};
use pineappl::bin::BinInfo;
use rayon::prelude::*;
use std::path::Path;

fn map_format_join(slice: &[f64]) -> String {
    slice.iter().map(|x| format!("{}", x)).join(", ")
}

fn map_format_e_join(slice: &[f64]) -> String {
    slice.iter().map(|x| format!("{:e}", x)).join(", ")
}

fn format_pdf_results(pdf_uncertainties: &[Vec<Vec<f64>>], pdfsets: &[&str]) -> String {
    let mut result = String::new();

    for (values, pdfset) in pdf_uncertainties.iter().zip(pdfsets.iter()) {
        result.push_str(&format!(
            "        (
            '{}',
            np.array([{}]),
            np.array([{}]),
            np.array([{}]),
        ),\n",
            pdfset.replace('_', "\\_"),
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
    pdfsets: &[&str],
    metadata: &[(&String, &String)],
) {
    println!("#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def percent_diff(a, b):
    return (a / b - 1.0) * 100.0

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
    axis.legend(fontsize='xx-small', frameon=False)

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

    axis.legend(fontsize='xx-small', frameon=False, ncol=2)
    minmax = axis.get_ylim()
    axis.set_yticks(np.arange(np.rint(minmax[0]), np.rint(minmax[1]) + 1.0, 1.0))
    axis.set_ylabel('PDF uncertainty [\\si{{\\percent}}]')

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

    axis.legend(fontsize='xx-small', frameon=False, ncol=2)
    minmax = axis.get_ylim()
    axis.set_yticks(np.arange(np.rint(minmax[0]), np.rint(minmax[1]) + 1.0, 1.0))
    axis.set_ylabel('Pull [$\\sigma$]')
    #axis.set_title('Comparison with ' + pdf_uncertainties[0][0], fontdict={{'fontsize': 9}}, loc='left')

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

pub fn subcommand(input: &str, pdfsets: &[&str], scales: usize) -> Result<()> {
    let grid = helpers::read_grid(input)?;
    let pdf = pdfsets[0].parse().map_or_else(
        |_| Pdf::with_setname_and_member(pdfsets[0], 0),
        Pdf::with_lhaid,
    );

    let results = helpers::convolute(&grid, &pdf, &[], &[], &[], scales);

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

        helpers::convolute(&grid, &pdf, &qcd_orders, &[], &[], scales)
    };

    let bin_info = grid.bin_info();

    let pdf_uncertainties: Vec<Vec<Vec<f64>>> = pdfsets
        .par_iter()
        .map(|pdfset| {
            let set = PdfSet::new(&pdfset.parse().map_or_else(
                |_| (*pdfset).to_string(),
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
                    .cloned()
                    .collect();

                let uncertainty = set.uncertainty(&values, 68.268949213708581, false);
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
        .chunks_exact(scales)
        .map(|variations| {
            variations
                .iter()
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .copied()
                .unwrap()
        })
        .collect();
    let max: Vec<_> = results
        .chunks_exact(scales)
        .map(|variations| {
            variations
                .iter()
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .copied()
                .unwrap()
        })
        .collect();

    let qcd_central: Vec<_> = qcd_results.iter().step_by(scales).copied().collect();
    let qcd_min: Vec<_> = qcd_results
        .chunks_exact(scales)
        .map(|variations| {
            variations
                .iter()
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .copied()
                .unwrap()
        })
        .collect();
    let qcd_max: Vec<_> = qcd_results
        .chunks_exact(scales)
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
                            .map(|map| map.get(&format!("x{}_label_tex", d + 1)).cloned())
                            .flatten()
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

    let mut output = Path::new(input);

    // remove ".lz4" and ".pineappl" extension
    match output.extension() {
        Some(x) if x == "lz4" => output = Path::new(output.file_stem().unwrap()),
        _ => {}
    }
    match output.extension() {
        Some(x) if x == "pineappl" => output = Path::new(output.file_stem().unwrap()),
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
        &pdfsets,
        &vector,
    );

    Ok(())
}
