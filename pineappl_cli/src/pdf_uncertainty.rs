use lhapdf::PdfSet;
use pineappl::grid::Grid;
use prettytable::{cell, row, Table};
use rayon::{prelude::*, ThreadPoolBuilder};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;

use super::helpers::create_table;

pub(crate) fn subcommand(
    input: &str,
    pdfset: &str,
    cl: f64,
    threads: usize,
    orders: &[(u32, u32)],
) -> Result<Table, Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;
    let set = PdfSet::new(&pdfset.parse().map_or_else(
        |_| pdfset.to_string(),
        |lhaid| lhapdf::lookup_pdf(lhaid).unwrap().0,
    ));
    let pdfs = set.mk_pdfs();

    let orders: Vec<_> = grid
        .orders()
        .iter()
        .map(|order| {
            orders.is_empty()
                || orders
                    .iter()
                    .any(|other| (order.alphas == other.0) && (order.alpha == other.1))
        })
        .collect();

    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    let results: Vec<f64> = pdfs
        .into_par_iter()
        .flat_map(|pdf| {
            grid.convolute(
                &|id, x, q2| pdf.xfx_q2(id, x, q2),
                &|id, x, q2| pdf.xfx_q2(id, x, q2),
                &|q2| pdf.alphas_q2(q2),
                &orders,
                &[],
                &[],
                &[(1.0, 1.0)],
            )
        })
        .collect();

    let bin_sizes = grid.bin_limits().bin_sizes();
    let bin_limits = grid.bin_limits().limits();

    let mut table = create_table();
    table.set_titles(row![c => "bin", "xmin", "xmax", "diff", "integ", "neg unc", "pos unc"]);

    for bin in 0..bin_sizes.len() {
        let values: Vec<_> = results
            .iter()
            .skip(bin)
            .step_by(bin_sizes.len())
            .cloned()
            .collect();
        let uncertainty = set.uncertainty(&values, cl, false);

        table.add_row(row![r =>
            &format!("{}", bin),
            &format!("{}", bin_limits[bin]),
            &format!("{}", bin_limits[bin + 1]),
            &format!("{:.7e}", uncertainty.central),
            &format!("{:.7e}", uncertainty.central * bin_sizes[bin]),
            &format!("{:.2}%", (-uncertainty.errminus / uncertainty.central) * 100.0),
            &format!("{:.2}%", (uncertainty.errplus / uncertainty.central) * 100.0),
        ]);
    }

    Ok(table)
}
