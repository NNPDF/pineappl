use lhapdf::Pdf;
use pineappl::grid::Grid;
use prettytable::{cell, Row, Table};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;

use super::helpers::create_table;

pub fn subcommand(
    input: &str,
    pdfset: &str,
    absolute: bool,
    normalize: &[(u32, u32)],
) -> Result<Table, Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;
    let pdf = pdfset
        .parse()
        .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);

    let grid_orders = grid.orders();
    let results = grid.convolute(
        &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
        &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
        &|q2| pdf.alphas_q2(q2),
        &[],
        &[],
        &[],
        &[(1.0, 1.0)],
    );

    let order_results: Vec<Vec<f64>> = (0..grid_orders.len())
        .map(|order| {
            let mut order_mask = vec![false; grid_orders.len()];
            order_mask[order] = true;
            grid.convolute(
                &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
                &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
                &|q2| pdf.alphas_q2(q2),
                &order_mask,
                &[],
                &[],
                &[(1.0, 1.0)],
            )
        })
        .collect();

    let mut sorted_grid_orders: Vec<_> = grid_orders
        .iter()
        .filter(|order| (order.logxir == 0) && (order.logxif == 0))
        .collect();
    sorted_grid_orders.sort();

    let unsorted_indices: Vec<_> = sorted_grid_orders
        .iter()
        .map(|sorted| {
            grid_orders
                .iter()
                .position(|unsorted| unsorted == *sorted)
                .unwrap()
        })
        .collect();
    let lo_power = {
        let order = sorted_grid_orders.first().unwrap();
        order.alphas + order.alpha
    };

    let bin_info = grid.bin_info();
    let left_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.left(i))
        .collect();
    let right_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.right(i))
        .collect();

    let mut title = Row::empty();
    title.add_cell(cell!(c->"bin"));
    for i in 0..bin_info.dimensions() {
        let mut cell = cell!(c->&format!("x{}", i + 1));
        cell.set_hspan(2);
        title.add_cell(cell);
    }
    title.add_cell(cell!(c->"diff"));

    for order in &sorted_grid_orders {
        title.add_cell(cell!(c->&format!("O(as^{} a^{})", order.alphas, order.alpha)));
    }

    let mut table = create_table();
    table.set_titles(title);

    for (bin, value) in results.iter().enumerate() {
        let row = table.add_empty_row();

        row.add_cell(cell!(r->&format!("{}", bin)));
        for (left, right) in left_limits.iter().zip(right_limits.iter()) {
            row.add_cell(cell!(r->&format!("{}", left[bin])));
            row.add_cell(cell!(r->&format!("{}", right[bin])));
        }
        row.add_cell(cell!(r->&format!("{:.7e}", value)));

        let mut normalization = 0.0;

        // calculate the sum of all leading orders
        for (index, order) in sorted_grid_orders.iter().enumerate() {
            if (normalize.is_empty() && ((order.alphas + order.alpha) == lo_power))
                || (normalize.iter().any(|o| *o == (order.alphas, order.alpha)))
            {
                normalization += order_results[unsorted_indices[index]][bin];
            }
        }

        // print each order normalized to the sum of all leading orders
        for index in 0..sorted_grid_orders.len() {
            let result = order_results[unsorted_indices[index]][bin];

            if absolute {
                row.add_cell(cell!(r->&format!("{:.7e}", result)));
            } else {
                row.add_cell(cell!(r->&format!("{:.2}%", result / normalization * 100.0)));
            }
        }
    }

    Ok(table)
}
