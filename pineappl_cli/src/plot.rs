use itertools::Itertools;
use lhapdf::{Pdf, PdfSet};
use pineappl::grid::Grid;
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;

use super::helpers;

pub fn subcommand(
    input: &str,
    pdfset: &str,
    other_pdfsets: &[&str],
    scales: usize,
) -> Result<(), Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;
    let pdf = pdfset
        .parse()
        .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);
    let scales_vector = vec![
        (1.0, 1.0),
        (2.0, 2.0),
        (0.5, 0.5),
        (2.0, 1.0),
        (1.0, 2.0),
        (0.5, 1.0),
        (1.0, 0.5),
        (2.0, 0.5),
        (0.5, 2.0),
    ];

    let results = helpers::convolute(&grid, &pdf, &[], &[], &[], &scales_vector[0..scales]);

    let qcd_results = {
        let mut orders = grid.orders().to_vec();
        orders.sort();
        let orders = orders;

        let qcd_orders: Vec<_> = orders
            .iter()
            .group_by(|order| order.alphas + order.alpha)
            .into_iter()
            .map(|mut group| group.1.next().unwrap())
            .collect();
        let qcd_order_mask: Vec<_> = grid
            .orders()
            .iter()
            .map(|order| {
                qcd_orders
                    .iter()
                    .any(|o| (order.alphas == o.alphas) && (order.alpha == o.alpha))
            })
            .collect();

        helpers::convolute(
            &grid,
            &pdf,
            &qcd_order_mask,
            &[],
            &[],
            &scales_vector[0..scales],
        )
    };

    let bin_info = grid.bin_info();

    let other_results: Vec<Vec<Vec<f64>>> = other_pdfsets
        .iter()
        .map(|pdfset| {
            let mut results = vec![];

            results.push(helpers::convolute(
                &grid,
                &pdfset
                    .parse()
                    .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid),
                &[],
                &[],
                &[],
                &[(1.0, 1.0)],
            ));

            let set = PdfSet::new(&pdfset.parse().map_or_else(
                |_| pdfset.to_string(),
                |lhaid| lhapdf::lookup_pdf(lhaid).unwrap().0,
            ));

            let pdf_results: Vec<_> = set
                .mk_pdfs()
                .into_par_iter()
                .flat_map(|pdf| helpers::convolute(&grid, &pdf, &[], &[], &[], &[(1.0, 1.0)]))
                .collect();

            let mut min = vec![];
            let mut max = vec![];

            for bin in 0..bin_info.bins() {
                let values: Vec<_> = pdf_results
                    .iter()
                    .skip(bin)
                    .step_by(bin_info.bins())
                    .cloned()
                    .collect();
                let uncertainty = set.uncertainty(&values, 68.268949213708581, false);
                min.push(uncertainty.central - uncertainty.errminus);
                max.push(uncertainty.central + uncertainty.errplus);
            }

            results.push(min);
            results.push(max);

            results
        })
        .collect();

    let set = PdfSet::new(&pdfset.parse().map_or_else(
        |_| pdfset.to_string(),
        |lhaid| lhapdf::lookup_pdf(lhaid).unwrap().0,
    ));
    let pdfs = set.mk_pdfs();

    let pdf_results: Vec<_> = pdfs
        .into_par_iter()
        //.iter()
        .flat_map(|pdf| helpers::convolute(&grid, &pdf, &[], &[], &[], &[(1.0, 1.0)]))
        .collect();

    let left_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.left(i))
        .collect();
    let right_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.right(i))
        .collect();

    let central: Vec<_> = results.iter().step_by(scales).collect();
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

    let mut pdf_unc_min = Vec::new();
    let mut pdf_unc_max = Vec::new();

    for bin in 0..bin_info.bins() {
        let values: Vec<_> = pdf_results
            .iter()
            .skip(bin)
            .step_by(bin_info.bins())
            .cloned()
            .collect();
        let uncertainty = set.uncertainty(&values, 68.268949213708581, false);
        pdf_unc_min.push(uncertainty.central - uncertainty.errminus);
        pdf_unc_max.push(uncertainty.central + uncertainty.errplus);
    }

    // the following implementation only works for 1D distributions
    assert_eq!(left_limits.len(), 1);
    assert_eq!(right_limits.len(), 1);

    println!("#!/usr/bin/env python3");
    println!("");

    println!("import matplotlib.pyplot as plt");
    println!("import numpy as np");
    println!("from matplotlib.transforms import ScaledTranslation");
    println!("from matplotlib.backends.backend_pdf import PdfPages");
    println!("");

    println!("def data():");
    println!("    left = np.array([");

    for value in left_limits[0]
        .iter()
        .chain(right_limits[0].last().iter().copied())
    {
        println!("        {},", value);
    }

    println!("    ])");
    println!("    right = np.array([");

    for value in &right_limits[0] {
        println!("        {},", value);
    }

    println!("    ])");
    println!("    central = np.array([");

    for value in central.iter().chain(central.last().iter().copied()) {
        println!("        {:e},", value);
    }

    println!("    ])");
    println!("    min = np.array([");

    for value in min.iter().chain(min.last().iter().copied()) {
        println!("        {:e},", value);
    }

    println!("    ])");
    println!("    max = np.array([");

    for value in max.iter().chain(max.last().iter().copied()) {
        println!("        {:e},", value);
    }

    println!("    ])");
    println!("    qcd_central = np.array([");

    for value in qcd_central.iter().chain(qcd_central.last().iter().copied()) {
        println!("        {:e},", value);
    }

    println!("    ])");
    println!("    qcd_min = np.array([");

    for value in qcd_min.iter().chain(qcd_min.last().iter().copied()) {
        println!("        {:e},", value);
    }

    println!("    ])");
    println!("    qcd_max = np.array([");

    for value in qcd_max.iter().chain(qcd_max.last().iter().copied()) {
        println!("        {:e},", value);
    }

    println!("    ])");
    println!("    pdf_unc_min = np.array([");

    for value in pdf_unc_min.iter().chain(pdf_unc_min.last().iter().copied()) {
        println!("        {:e},", value);
    }

    println!("    ])");
    println!("    pdf_unc_max = np.array([");

    for value in pdf_unc_max.iter().chain(pdf_unc_max.last().iter().copied()) {
        println!("        {:e},", value);
    }

    println!("    ])");
    println!("    other_results = [");

    for (values, pdfset) in other_results.iter().zip(other_pdfsets.iter()) {
        println!("        (");
        println!("            '{}',", pdfset.replace('_', "\\_"));
        println!("            np.array([");

        for value in &values[0] {
            println!("                {:e},", value);
        }

        println!("            ]),");
        println!("            np.array([");

        for value in &values[1] {
            println!("                {:e},", value);
        }

        println!("            ]),");
        println!("            np.array([");

        for value in &values[2] {
            println!("                {:e},", value);
        }

        println!("            ]),");
        println!("        ),");
    }

    println!("    ]");

    println!("");
    println!("    return {{ 'pdfset': '{}',", pdfset);
    println!("             'left': left,");
    println!("             'right': right,");
    println!("             'central': central,");
    println!("             'min': min,");
    println!("             'max': max,");
    println!("             'qcd_central': qcd_central,");
    println!("             'qcd_min': qcd_min,");
    println!("             'qcd_max': qcd_max,");
    println!("             'pdf_unc_min': pdf_unc_min,");
    println!("             'pdf_unc_max': pdf_unc_max,");
    println!("             'other_results': other_results,");
    println!("    }}");
    println!("");

    println!("def metadata():");
    println!("    return {{");

    for (key, value) in grid.key_values().unwrap_or(&HashMap::new()) {
        // skip multi-line entries
        if value.contains('\n') {
            continue;
        }

        match key.as_str() {
            "description" => println!(
                "        '{}': r'{}',",
                key,
                value.replace("–", "--").replace("—", "---")
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
    println!("");

    println!("def plot_abs(axis):");
    println!("    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)");
    println!("    axis.minorticks_on()");
    println!("");
    println!("    x = data()['left']");
    println!("    y = data()['central']");
    println!("    ymin = data()['min']");
    println!("    ymax = data()['max']");
    println!(
        "    ylabel = metadata()['y_label_tex'] + r' [\\si{{' + metadata()['y_unit'] + r'}}]'"
    );
    println!("    description = metadata()['description']");
    println!("");
    println!("    axis.set_yscale('log')");
    println!("    axis.set_xscale('log')");
    println!("    axis.step(x, y, 'royalblue', linewidth=1.0, where='post', zorder=-1.0)");
    println!("    axis.fill_between(x, ymin, ymax, alpha=0.4, color='royalblue', linewidth=0, step='post', zorder=-1.0)");
    println!("    axis.set_ylabel(ylabel)");
    println!("    axis.set_title(description)");
    println!("");

    println!("def plot_rel_ewonoff(axis):");
    println!("    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)");
    println!("    axis.minorticks_on()");
    println!("");
    println!("    qcd_y = (data()['qcd_central'] / data()['qcd_central'] - 1.0) * 100.0");
    println!("    qcd_ymin = (data()['qcd_min'] / data()['qcd_central'] - 1.0) * 100.0");
    println!("    qcd_ymax = (data()['qcd_max'] / data()['qcd_central'] - 1.0) * 100.0");
    println!("    x = data()['left']");
    println!("    y = (data()['central'] / data()['qcd_central'] - 1.0) * 100.0");
    println!("    ymin = (data()['min'] / data()['qcd_central'] - 1.0) * 100.0");
    println!("    ymax = (data()['max'] / data()['qcd_central'] - 1.0) * 100.0");
    println!("");
    println!("    axis.step(x, qcd_y, 'red', linewidth=1.0, where='post', zorder=-1.0)");
    println!("    axis.fill_between(x, qcd_ymin, qcd_ymax, alpha=0.4, color='red', label='NLO EW off', linewidth=0, step='post', zorder=-1.0)");
    println!("    axis.step(x, y, 'royalblue', linewidth=1.0, where='post', zorder=-1.0)");
    println!("    axis.fill_between(x, ymin, ymax, alpha=0.4, color='royalblue', label='NLO EW on', linewidth=0, step='post', zorder=-1.0)");
    println!("    axis.set_ylabel('NLO EW on/off [\\si{{\\percent}}]')");
    println!("    axis.legend(fontsize='xx-small')");
    println!("");

    println!("def plot_rel_pdfcomp(axis):");
    println!("    axis.tick_params(axis='both', left=True, right=True, top=True, bottom=True, which='both', direction='in', width=0.5, zorder=10.0)");
    println!("    axis.minorticks_on()");
    println!("");
    println!("    x = 0.5 * (data()['left'][:-1] + data()['right'])");
    println!("    y = (data()['central'] / data()['central'] - 1.0)[:-1]");
    println!("    ymin = (abs(data()['pdf_unc_min'] / data()['central'] - 1.0) * 100.0)[:-1]");
    println!("    ymax = (abs(data()['pdf_unc_max'] / data()['central'] - 1.0) * 100.0)[:-1]");
    println!("    other = data()['other_results']");
    println!("    pdfset = data()['pdfset']");
    println!("    xlabel = metadata()['x1_label_tex'] + (r' [\\si{{' + metadata()['x1_unit'] + r'}}]' if metadata()['x1_unit'] != '' else '')");
    println!("");
    println!("    translation = axis.transData");
    println!("    width = 0.06");
    println!("    translation += ScaledTranslation(-0.5 * width, 0, axis.figure.dpi_scale_trans)");
    println!("");
    println!("    colors = ['royalblue', 'brown', 'darkorange', 'darkgreen', 'purple']");
    println!("    fmt = ['o', '^', 's', 'D', 'x']");
    println!("    axis.errorbar(x, y, yerr=(ymin, ymax), color=colors[0], label=pdfset, fmt=fmt[0], capsize=1, markersize=2, linewidth=1, zorder=-1.0, transform=translation)");
    println!("");
    println!("    for index, i in enumerate(other):");
    println!("        label, y, ymin, ymax = i");
    println!("        ymin = abs(ymin / y - 1.0) * 100.0");
    println!("        ymax = abs(ymax / y - 1.0) * 100.0");
    println!("        y = (y / (data()['central'])[:-1] - 1.0) * 100.0");
    println!("");
    println!("        translation += ScaledTranslation(width / len(other), 0, axis.figure.dpi_scale_trans)");
    println!("        axis.errorbar(x, y, yerr=(ymin, ymax), color=colors[index+1], label=label, fmt=fmt[index+1], capsize=1, markersize=2, linewidth=1, zorder=-1.0, transform=translation)");
    println!("");
    println!("    axis.legend(fontsize='xx-small')");
    println!("    axis.set_xlabel(xlabel)");
    println!("    axis.set_ylabel('PDF uncertainty [\\si{{\\percent}}]')");
    println!("");

    println!("plt.rc('text', usetex=True)");
    println!("plt.rc('text.latex', preamble=r'\\usepackage{{siunitx}}')");
    println!("plt.rc('figure', figsize=(6.4,8.4))");
    println!("plt.rc('font', family='serif', size=14.0)");
    println!("plt.rc('axes', labelsize='large')");
    println!("plt.rc('pdf', compression=0)");
    println!("");
    println!("with PdfPages('output.pdf') as pp:");
    println!("    figure, axis = plt.subplots(3, 1, sharex=True, gridspec_kw={{'height_ratios' : [1.5, 1, 1]}})");
    println!(
        "    figure.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.6, rect=(0.0475,0.03,1.01,0.975))"
    );
    println!("    plot_abs(axis[0])");
    println!("    plot_rel_ewonoff(axis[1])");
    println!("    plot_rel_pdfcomp(axis[2])");
    println!("    figure.savefig(pp, format='pdf')");
    println!("    plt.close()");

    Ok(())
}
