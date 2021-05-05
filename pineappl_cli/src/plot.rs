use super::helpers;
use anyhow::Result;
use itertools::Itertools;
use lhapdf::{Pdf, PdfSet};
use rayon::prelude::*;

fn map_format_join(slice: &[f64]) -> String {
    slice.iter().map(|x| format!("{}", x)).join(", ")
}

fn map_format_e_join(slice: &[f64]) -> String {
    slice.iter().map(|x| format!("{:e}", x)).join(", ")
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

    println!("#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def percent_diff(a, b):
    return (a / b - 1.0) * 100.0

def plot_abs(axis, **kwargs):
    x = kwargs['x']
    y = kwargs['y']
    ymin = kwargs['ymin']
    ymax = kwargs['ymax']
    ylog = kwargs['ylog']
    ylabel = kwargs['ylabel']

    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()
    axis.set_yscale('log' if ylog else 'linear')
    axis.set_axisbelow(True)
    axis.grid(linestyle='dotted')
    axis.step(x, y, 'royalblue', linewidth=1.0, where='post')
    axis.fill_between(x, ymin, ymax, alpha=0.4, color='royalblue', linewidth=0.5, step='post')
    axis.set_ylabel(ylabel)

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
    plt.rc('figure', figsize=(6.4,len(panels)*2.4))
    plt.rc('font', family='serif', size=14.0)
    plt.rc('axes', labelsize='small')
    plt.rc('pdf', compression=0)

    xaxis = '{}'
    xunit = metadata().get(xaxis + '_unit', '')
    xlabel = metadata()[xaxis + '_label_tex'] + (r' [\\si{{' + xunit + r'}}]' if xunit != '' else '')
    ylabel = metadata()['y_label_tex'] + r' [\\si{{' + metadata()['y_unit'] + r'}}]'
    ylog = xunit != ''
    description = metadata()['description']

    for index, dict in enumerate(data_slices):
        dict['xlabel'] = xlabel
        dict['ylabel'] = ylabel
        dict['ylog'] = ylog
        figure, axes = plt.subplots(len(panels), 1, sharex=True)
        figure.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.6, rect=(0.06,0.03,1.0,0.975))

        axes[0].set_title(description)
        axes[-1].set_xlabel(xlabel)

        if xunit != '':
            axes[0].set_xscale('log')

        for plot, axis in zip(panels, axes):
            plot(axis, **dict)

        name = 'output' if len(data_slices) == 1 else 'output-{{}}'.format(index)
        figure.savefig(name + '.pdf')

def data():
    left = np.array([{}])
    right = np.array([{}])
    min = np.array([{}])
    max = np.array([{}])
    qcd_central = np.array([{}])
    qcd_min = np.array([{}])
    qcd_max = np.array([{}])
    slices = {:?}
    pdf_results = [",
        format!("x{}", bin_info.dimensions()),
        map_format_join(left_limits.last().unwrap()),
        map_format_join(right_limits.last().unwrap()),
        map_format_e_join(&min),
        map_format_e_join(&max),
        map_format_e_join(&qcd_central),
        map_format_e_join(&qcd_min),
        map_format_e_join(&qcd_max),
        slices,
    );

    for (values, pdfset) in pdf_uncertainties.iter().zip(pdfsets.iter()) {
        println!(
            "        (
            '{}',
            np.array([{}]),
            np.array([{}]),
            np.array([{}]),
        ),",
            pdfset.replace('_', "\\_"),
            map_format_e_join(&values[0]),
            map_format_e_join(&values[1]),
            map_format_e_join(&values[2]),
        );
    }

    println!(
        "    ]

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
    }} for slice in slices]

def metadata():
    return {{"
    );

    let mut key_values = grid.key_values().cloned().unwrap_or_default();
    key_values.entry("description".to_string()).or_default();
    key_values.entry("x1_label_tex".to_string()).or_default();
    key_values.entry("x1_unit".to_string()).or_default();
    key_values.entry("y_label_tex".to_string()).or_default();
    key_values.entry("y_unit".to_string()).or_default();

    let mut vector: Vec<_> = key_values.iter().collect();
    vector.sort();
    let vector = vector;

    for (key, value) in &vector {
        // skip multi-line entries
        if value.contains('\n') {
            continue;
        }

        println!(
            "        '{}': r'{}',",
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
        );
    }

    println!("    }}");
    println!();
    println!("if __name__ == '__main__':");
    println!("    main()");

    Ok(())
}
