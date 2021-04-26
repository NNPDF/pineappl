use super::helpers;
use anyhow::Result;
use itertools::Itertools;
use lhapdf::{Pdf, PdfSet};
use rayon::prelude::*;

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
                .unwrap()
        })
        .collect();
    let max: Vec<_> = results
        .chunks_exact(scales)
        .map(|variations| {
            variations
                .iter()
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap()
        })
        .collect();

    let qcd_central: Vec<_> = qcd_results.iter().step_by(scales).collect();
    let qcd_min: Vec<_> = qcd_results
        .chunks_exact(scales)
        .map(|variations| {
            variations
                .iter()
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap()
        })
        .collect();
    let qcd_max: Vec<_> = qcd_results
        .chunks_exact(scales)
        .map(|variations| {
            variations
                .iter()
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap()
        })
        .collect();

    // the following implementation only works for 1D distributions
    assert_eq!(left_limits.len(), 1);
    assert_eq!(right_limits.len(), 1);

    print!("#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import ScaledTranslation
from matplotlib.backends.backend_pdf import PdfPages

def plot_abs(axis):
    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()

    x = data()['left']
    y = data()['pdf_results'][0][1]
    ymin = data()['min']
    ymax = data()['max']
    x1_unit = metadata().get('x1_unit', '')
    ylabel = metadata()['y_label_tex'] + r' [\\si{{' + metadata()['y_unit'] + r'}}]'
    description = metadata()['description']

    if x1_unit != '':
        axis.set_yscale('log')
    axis.set_axisbelow(True)
    axis.grid(linestyle='dotted')
    axis.step(x, y, 'royalblue', linewidth=1.0, where='post')
    axis.fill_between(x, ymin, ymax, alpha=0.4, color='royalblue', linewidth=0.5, step='post')
    axis.set_ylabel(ylabel)
    axis.set_title(description)

def plot_rel_ewonoff(axis):
    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()

    qcd_y = (data()['qcd_central'] / data()['qcd_central'] - 1.0) * 100.0
    qcd_ymin = (data()['qcd_min'] / data()['qcd_central'] - 1.0) * 100.0
    qcd_ymax = (data()['qcd_max'] / data()['qcd_central'] - 1.0) * 100.0
    x = data()['left']
    y = (data()['pdf_results'][0][1] / data()['qcd_central'] - 1.0) * 100.0
    ymin = (data()['min'] / data()['qcd_central'] - 1.0) * 100.0
    ymax = (data()['max'] / data()['qcd_central'] - 1.0) * 100.0
    pdf_min = (abs(data()['pdf_results'][0][2] / data()['pdf_results'][0][1] - 1.0) * 100.0)[:-1]
    pdf_max = (abs(data()['pdf_results'][0][3] / data()['pdf_results'][0][1] - 1.0) * 100.0)[:-1]
    mid = 0.5 * (data()['left'][:-1] + data()['right'])

    axis.set_axisbelow(True)
    axis.grid(linestyle='dotted')
    axis.step(x, qcd_y, 'red', label='NLO QCD', linewidth=1.0, where='post')
    #axis.fill_between(x, qcd_ymin, qcd_ymax, alpha=0.4, color='red', label='7-p.\\ scale var.', linewidth=0.5, step='post')
    axis.step(x, y, 'royalblue', label='NLO QCD+EW', linewidth=1.0, where='post')
    axis.fill_between(x, ymin, ymax, alpha=0.4, color='royalblue', label='7-p.\\ scale var.', linewidth=0.5, step='post')
    axis.errorbar(mid, y[:-1], yerr=(pdf_min, pdf_max), color='royalblue', label='PDF uncertainty', fmt='.', capsize=1, markersize=0, linewidth=1)
    axis.set_ylabel('NLO EW on/off [\\si{{\\percent}}]')
    axis.legend(fontsize='xx-small', frameon=False)

def plot_rel_pdfunc(axis):
    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()

    x = data()['left']
    pdf_uncertainties = data()['pdf_results']

    colors = ['royalblue', 'brown', 'darkorange', 'darkgreen', 'purple', 'tan']

    axis.set_axisbelow(True)
    axis.grid(linestyle='dotted')

    #ymins = np.asmatrix([(ymin / y - 1.0) * 100 for label, y, ymin, ymax in pdf_uncertainties])
    #ymaxs = np.asmatrix([(ymax / y - 1.0) * 100 for label, y, ymin, ymax in pdf_uncertainties])

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

def plot_rel_pdfpull(axis):
    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()

    x = data()['left']
    pdf_uncertainties = data()['pdf_results']

    colors = ['royalblue', 'brown', 'darkorange', 'darkgreen', 'purple', 'tan']

    axis.set_axisbelow(True)
    axis.grid(linestyle='dotted')

    central_y = pdf_uncertainties[0][1]
    central_ymin = pdf_uncertainties[0][2]
    central_ymax = pdf_uncertainties[0][3]

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

def plot_rel_pdfdiff(axis):
    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)
    axis.minorticks_on()

    x = data()['left']
    pdf_uncertainties = data()['pdf_results']

    colors = ['royalblue', 'brown', 'darkorange', 'darkgreen', 'purple', 'tan']

    axis.grid(linestyle='dotted')

    central_y = pdf_uncertainties[0][1]
    central_ymin = pdf_uncertainties[0][2]
    central_ymax = pdf_uncertainties[0][3]

    for index, i in enumerate(pdf_uncertainties):
        label, y, ymin, ymax = i
        pull_max = (y - central_y) / np.sqrt(np.power(np.minimum(ymax - y, y - ymin), 2) + np.power(np.minimum(central_ymax - central_y, central_y - central_ymin), 2))
        pull_min = (y - central_y) / np.sqrt(np.power(np.maximum(ymax - y, y - ymin), 2) + np.power(np.maximum(central_ymax - central_y, central_y - central_ymin), 2))
        diff = (y / central_y - 1.0) * 100.0

        #axis.fill_between(x, pull_min, pull_max, alpha=0.4, color=colors[index], label=label, linewidth=0.5, step='post')
        axis.step(x, diff, color=colors[index], label=label, linewidth=1, where='post')

    axis.legend(fontsize='xx-small', frameon=False, ncol=2)
    minmax = axis.get_ylim()
    axis.set_yticks(np.arange(np.rint(minmax[0]), np.rint(minmax[1]) + 1.0, 2.0))
    axis.set_ylabel('Difference [\\si{{\\percent}}]')
    #axis.set_title('Comparison with ' + pdf_uncertainties[0][0], fontdict={{'fontsize': 9}}, loc='left')

def main():
    panels = [
        plot_abs,
        plot_rel_ewonoff,
        plot_rel_pdfunc,
        #plot_rel_pdfdiff,
        plot_rel_pdfpull,
    ]

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\\usepackage{{siunitx}}')
    plt.rc('figure', figsize=(6.4,len(panels)*2.4))
    plt.rc('font', family='serif', size=14.0)
    plt.rc('axes', labelsize='small')
    plt.rc('pdf', compression=0)

    x1_unit = metadata().get('x1_unit', '')
    xlabel = metadata()['x1_label_tex'] + (r' [\\si{{' + x1_unit + r'}}]' if x1_unit != '' else '')

    with PdfPages('output.pdf') as pp:
        figure, axis = plt.subplots(len(panels), 1, sharex=True)
        figure.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.6, rect=(0.0475,0.03,1.01,0.975))

        if x1_unit != '':
            axis[0].set_xscale('log')

        for index, plot in enumerate(panels):
            plot(axis[index])

        axis[-1].set_xlabel(xlabel)
        figure.savefig(pp, format='pdf')
        plt.close()

def data():
    left = np.array([");

    left_limits[0]
        .iter()
        .chain(right_limits[0].last().iter().copied())
        .for_each(|x| print!("{}, ", x));

    println!("])");
    print!("    right = np.array([");

    right_limits[0].iter().for_each(|x| print!("{}, ", x));

    println!("])");
    print!("    min = np.array([");

    min.iter()
        .chain(min.last().iter().copied())
        .for_each(|x| print!("{:e}, ", x));

    println!("])");
    print!("    max = np.array([");

    max.iter()
        .chain(max.last().iter().copied())
        .for_each(|x| print!("{:e}, ", x));

    println!("])");
    print!("    qcd_central = np.array([");

    qcd_central
        .iter()
        .chain(qcd_central.last().iter().copied())
        .for_each(|x| print!("{:e}, ", x));

    println!("])");
    print!("    qcd_min = np.array([");

    qcd_min
        .iter()
        .chain(qcd_min.last().iter().copied())
        .for_each(|x| print!("{:e}, ", x));

    println!("])");
    print!("    qcd_max = np.array([");

    qcd_max
        .iter()
        .chain(qcd_max.last().iter().copied())
        .for_each(|x| print!("{:e}, ", x));

    println!("])");
    println!("    pdf_results = [");

    for (values, pdfset) in pdf_uncertainties.iter().zip(pdfsets.iter()) {
        println!("        (");
        println!("            '{}',", pdfset.replace('_', "\\_"));
        print!("            np.array([");

        values[0]
            .iter()
            .chain(values[0].last().iter().copied())
            .for_each(|x| print!("{:e}, ", x));

        println!("]),");
        print!("            np.array([");

        values[1]
            .iter()
            .chain(values[1].last().iter().copied())
            .for_each(|x| print!("{:e}, ", x));

        println!("]),");
        print!("            np.array([");

        values[2]
            .iter()
            .chain(values[2].last().iter().copied())
            .for_each(|x| print!("{:e}, ", x));

        println!("]),");
        println!("        ),");
    }

    println!("    ]");

    println!();
    println!("    return {{ 'left': left,");
    println!("             'right': right,");
    println!("             'min': min,");
    println!("             'max': max,");
    println!("             'qcd_central': qcd_central,");
    println!("             'qcd_min': qcd_min,");
    println!("             'qcd_max': qcd_max,");
    println!("             'pdf_results': pdf_results,");
    println!("    }}");
    println!();

    println!("def metadata():");
    println!("    return {{");

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

        match key.as_str() {
            "description" => println!(
                "        '{}': r'{}',",
                key,
                value.replace("\u{2013}", "--").replace("\u{2014}", "---")
            ),
            "x1_unit" | "x2_unit" | "x3_unit" | "y_unit" => println!(
                "        '{}': r'{}',",
                key,
                value
                    .replace("GeV", r#"\giga\electronvolt"#)
                    .replace("/", r#"\per"#)
                    .replace("pb", r#"\pico\barn"#)
            ),
            _ => println!("        '{}': r'{}',", key, value),
        }
    }

    println!("    }}");
    println!();
    println!("if __name__ == '__main__':");
    println!("    main()");

    Ok(())
}
