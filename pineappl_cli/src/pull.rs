use super::helpers;
use anyhow::Result;
use lhapdf::PdfSet;
use prettytable::{cell, Row, Table};
use rayon::{prelude::*, ThreadPoolBuilder};

pub fn subcommand(
    input: &str,
    pdfset1: &str,
    pdfset2: &str,
    cl: f64,
    limit: usize,
    threads: usize,
) -> Result<Table> {
    let grid = helpers::read_grid(input)?;

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

    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    let results1: Vec<f64> = pdfset1
        .par_iter()
        .flat_map(|pdf| helpers::convolute(&grid, pdf, &[], &[], &[], 1))
        .collect();
    let results2: Vec<f64> = pdfset2
        .par_iter()
        .flat_map(|pdf| helpers::convolute(&grid, pdf, &[], &[], &[], 1))
        .collect();

    let bin_info = grid.bin_info();
    let left_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.left(i))
        .collect();
    let right_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.right(i))
        .collect();

    let labels = helpers::labels(&grid);
    let (_, x_labels) = labels.split_last().unwrap();
    let mut title = Row::empty();
    title.add_cell(cell!(c->"bin"));
    for x_label in x_labels {
        let mut cell = cell!(c->&x_label);
        cell.set_hspan(2);
        title.add_cell(cell);
    }
    title.add_cell(cell!(c->"total"));
    for _ in 0..limit {
        title.add_cell(cell!(c->"lumi"));
        title.add_cell(cell!(c->"pull"));
    }

    let mut table = helpers::create_table();
    table.set_titles(title);

    for bin in 0..bin_info.bins() {
        let values1: Vec<_> = results1
            .iter()
            .skip(bin)
            .step_by(bin_info.bins())
            .copied()
            .collect();
        let values2: Vec<_> = results2
            .iter()
            .skip(bin)
            .step_by(bin_info.bins())
            .copied()
            .collect();
        let uncertainty1 = set1.uncertainty(&values1, cl, false);
        let uncertainty2 = set2.uncertainty(&values2, cl, false);

        let lumi_results1: Vec<_> = (0..grid.lumi().len())
            .map(|lumi| {
                let mut lumi_mask = vec![false; grid.lumi().len()];
                lumi_mask[lumi] = true;
                let central: Vec<f64> = pdfset1
                    .iter()
                    .flat_map(|pdf| helpers::convolute(&grid, pdf, &[], &[], &lumi_mask, 1))
                    .collect();
                set1.uncertainty(&central, cl, false).central
            })
            .collect();
        let lumi_results2: Vec<_> = (0..grid.lumi().len())
            .map(|lumi| {
                let mut lumi_mask = vec![false; grid.lumi().len()];
                lumi_mask[lumi] = true;
                let central: Vec<f64> = pdfset2
                    .iter()
                    .flat_map(|pdf| helpers::convolute(&grid, pdf, &[], &[], &lumi_mask, 1))
                    .collect();
                set1.uncertainty(&central, cl, false).central
            })
            .collect();

        let mut pull_tuples: Vec<_> = lumi_results2
            .iter()
            .zip(lumi_results1.iter())
            .map(|(res2, res1)| {
                // use the uncertainties in the direction in which the respective results differ
                let unc1 = if res1 > res2 {
                    uncertainty1.errminus
                } else {
                    uncertainty1.errplus
                };
                let unc2 = if res2 > res1 {
                    uncertainty2.errminus
                } else {
                    uncertainty2.errplus
                };
                (res2 - res1) / unc1.hypot(unc2)
            })
            .enumerate()
            .collect();

        let total = pull_tuples
            .iter()
            .fold(0.0, |value, (_, pull)| value + pull);

        let row = table.add_empty_row();

        row.add_cell(cell!(r->&format!("{}", bin)));
        for (left, right) in left_limits.iter().zip(right_limits.iter()) {
            row.add_cell(cell!(r->&format!("{}", left[bin])));
            row.add_cell(cell!(r->&format!("{}", right[bin])));
        }

        row.add_cell(cell!(r->&format!("{:.3}", total)));

        // sort using the absolute value in descending order
        pull_tuples.sort_unstable_by(|(_, pull_left), (_, pull_right)| {
            pull_right.abs().partial_cmp(&pull_left.abs()).unwrap()
        });

        for (lumi, pull) in pull_tuples.iter().take(limit) {
            row.add_cell(cell!(r->&format!("#{}", lumi)));
            row.add_cell(cell!(r->&format!("{:.3}", pull)));
        }
    }

    Ok(table)
}
