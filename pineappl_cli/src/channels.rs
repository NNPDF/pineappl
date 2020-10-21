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
    let limit = if lumis.is_empty() {
        grid.lumi().len().min(limit)
    } else {
        lumis.iter().filter(|lumi| **lumi < lumis.len()).count()
    };
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

    let bin_info = grid.bin_info();
    let left_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.left(i))
        .collect();
    let right_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.right(i))
        .collect();

    let mut title_row = row![];
    title_row.add_cell(cell!(c->"bin"));
    for i in 0..bin_info.dimensions() {
        let mut cell = cell!(c->&format!("x{}", i + 1));
        cell.set_hspan(2);
        title_row.add_cell(cell);
    }
    for _ in 0..limit {
        title_row.add_cell(cell!(c->"lumi"));
        title_row.add_cell(cell!(c->"size"));
    }

    let mut table = create_table();
    table.set_titles(title_row);

    for bin in 0..bin_info.bins() {
        let row = table.add_empty_row();

        row.add_cell(cell!(r->&format!("{}", bin)));

        for (left, right) in left_limits.iter().zip(right_limits.iter()) {
            row.add_cell(cell!(r->&format!("{}", left[bin])));
            row.add_cell(cell!(r->&format!("{}", right[bin])));
        }

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
