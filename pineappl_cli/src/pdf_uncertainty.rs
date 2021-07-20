use super::helpers;
use anyhow::Result;
use lhapdf::PdfSet;
use prettytable::{cell, Row, Table};
use rayon::{prelude::*, ThreadPoolBuilder};

pub fn subcommand(
    input: &str,
    pdfset: &str,
    cl: f64,
    threads: usize,
    orders: &[(u32, u32)],
    integrated: bool,
) -> Result<Table> {
    let grid = helpers::read_grid(input)?;
    let set = PdfSet::new(&pdfset.parse().map_or_else(
        |_| pdfset.to_string(),
        |lhaid| lhapdf::lookup_pdf(lhaid).unwrap().0,
    ));
    let pdfs = set.mk_pdfs();

    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    let results: Vec<f64> = pdfs
        .into_par_iter()
        .flat_map(|pdf| helpers::convolute(&grid, &pdf, orders, &[], &[], 1))
        .collect();

    let bin_info = grid.bin_info();
    let left_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.left(i))
        .collect();
    let right_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.right(i))
        .collect();
    let normalizations = bin_info.normalizations();

    let labels = helpers::labels(&grid);
    let (y_label, x_labels) = labels.split_last().unwrap();
    let mut title = Row::empty();
    title.add_cell(cell!(c->"bin"));
    for x_label in x_labels {
        let mut cell = cell!(c->&x_label);
        cell.set_hspan(2);
        title.add_cell(cell);
    }
    title.add_cell(cell!(c->if integrated { "integ" } else { y_label }));
    title.add_cell(cell!(c->"neg unc"));
    title.add_cell(cell!(c->"pos unc"));

    let mut table = helpers::create_table();
    table.set_titles(title);

    for bin in 0..bin_info.bins() {
        let values: Vec<_> = results
            .iter()
            .skip(bin)
            .step_by(bin_info.bins())
            .copied()
            .collect();
        let uncertainty = set.uncertainty(&values, cl, false);

        let row = table.add_empty_row();

        row.add_cell(cell!(r->&format!("{}", bin)));
        for (left, right) in left_limits.iter().zip(right_limits.iter()) {
            row.add_cell(cell!(r->&format!("{}", left[bin])));
            row.add_cell(cell!(r->&format!("{}", right[bin])));
        }
        row.add_cell(cell!(r->&format!("{:.7e}", if integrated { uncertainty.central * normalizations[bin] } else { uncertainty.central })));
        row.add_cell(
            cell!(r->&format!("{:.2}%", (-uncertainty.errminus / uncertainty.central) * 100.0)),
        );
        row.add_cell(
            cell!(r->&format!("{:.2}%", (uncertainty.errplus / uncertainty.central) * 100.0)),
        );
    }

    Ok(table)
}
