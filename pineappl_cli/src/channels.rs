use lhapdf::Pdf;
use pineappl::grid::Grid;
use prettytable::{cell, row, Table};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;

use super::helpers::create_table;

pub fn subcommand(
    input: &str,
    pdfset: &str,
    limit: usize,
    orders: &[(u32, u32)],
    absolute: bool,
    lumis: &[usize],
) -> Result<Table, Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;
    let pdf = pdfset
        .parse()
        .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);
    let limit = if lumis.is_empty() { limit } else { usize::MAX };
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

    let results: Vec<_> = (0..grid.lumi().len())
        .map(|lumi| {
            let mut lumi_mask = vec![false; grid.lumi().len()];
            lumi_mask[lumi] = true;
            grid.convolute(
                &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
                &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
                &|q2| pdf.alphas_q2(q2),
                &orders,
                &[],
                &lumi_mask,
                &[(1.0, 1.0)],
            )
        })
        .collect();

    let bin_limits = grid.bin_limits().limits();

    let mut table = create_table();
    table.set_titles(row![c => "bin", "xmin", "xmax", "lumi", "size"]);

    // TODO: add more titles

    for bin in 0..grid.bin_limits().bins() {
        let row = table.add_empty_row();

        row.add_cell(cell!(r->&format!("{}", bin)));
        row.add_cell(cell!(r->&format!("{}", bin_limits[bin])));
        row.add_cell(cell!(r->&format!("{}", bin_limits[bin + 1])));

        if absolute {
            let mut values: Vec<_> = results
                .iter()
                .enumerate()
                .map(|(lumi, vec)| (lumi, vec[bin]))
                .collect();

            // sort using the absolute value in descending order
            values.sort_unstable_by(|(_, left), (_, right)| {
                right.abs().partial_cmp(&left.abs()).unwrap()
            });

            for (lumi, value) in values
                .iter()
                .filter(|(lumi, _)| lumis.is_empty() || lumis.iter().any(|l| l == lumi))
                .take(limit)
            {
                row.add_cell(cell!(r->&format!("#{}", lumi)));
                row.add_cell(cell!(r->&format!("{:.7e}", value)));
            }
        } else {
            let sum: f64 = results.iter().map(|vec| vec[bin]).sum();
            let mut percentages: Vec<_> = results
                .iter()
                .enumerate()
                .map(|(lumi, vec)| (lumi, vec[bin] / sum * 100.0))
                .collect();

            // sort using the absolute value in descending order
            percentages.sort_unstable_by(|(_, left), (_, right)| {
                right.abs().partial_cmp(&left.abs()).unwrap()
            });

            for (lumi, percentage) in percentages
                .iter()
                .filter(|(lumi, _)| lumis.is_empty() || lumis.iter().any(|l| l == lumi))
                .take(limit)
            {
                row.add_cell(cell!(r->&format!("#{}", lumi)));
                row.add_cell(cell!(r->&format!("{:.2}%", percentage)));
            }
        }
    }

    Ok(table)
}
