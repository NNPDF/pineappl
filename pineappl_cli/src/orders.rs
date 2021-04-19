use super::helpers;
use anyhow::Result;
use lhapdf::Pdf;
use prettytable::{cell, Row, Table};

pub fn subcommand(
    input: &str,
    pdfset: &str,
    absolute: bool,
    normalize: &[(u32, u32)],
    integrated: bool,
) -> Result<Table> {
    let grid = helpers::read_grid(input)?;
    let pdf = pdfset
        .parse()
        .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);

    let mut orders: Vec<_> = grid
        .orders()
        .iter()
        .filter(|order| (order.logxir == 0) && (order.logxif == 0))
        .collect();
    orders.sort();
    let orders = orders;

    let results: Vec<Vec<f64>> = orders
        .iter()
        .map(|order| helpers::convolute(&grid, &pdf, &[(order.alphas, order.alpha)], &[], &[], 1))
        .collect();

    let lo_power = {
        let order = orders.first().unwrap();
        order.alphas + order.alpha
    };

    let bin_info = grid.bin_info();
    let left_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.left(i))
        .collect();
    let right_limits: Vec<_> = (0..bin_info.dimensions())
        .map(|i| bin_info.right(i))
        .collect();
    let normalizations = bin_info.normalizations();

    let mut title = Row::empty();
    title.add_cell(cell!(c->"bin"));
    for i in 0..bin_info.dimensions() {
        let mut cell = cell!(c->&format!("x{}", i + 1));
        cell.set_hspan(2);
        title.add_cell(cell);
    }
    title.add_cell(cell!(c->if integrated { "integ" } else { "diff" }));

    for order in &orders {
        title.add_cell(cell!(c->&format!("O(as^{} a^{})", order.alphas, order.alpha)));
    }

    let mut table = helpers::create_table();
    table.set_titles(title);

    for bin in 0..bin_info.bins() {
        let row = table.add_empty_row();
        let bin_norm = if integrated { normalizations[bin] } else { 1.0 };

        row.add_cell(cell!(r->&format!("{}", bin)));
        for (left, right) in left_limits.iter().zip(right_limits.iter()) {
            row.add_cell(cell!(r->&format!("{}", left[bin])));
            row.add_cell(cell!(r->&format!("{}", right[bin])));
        }
        row.add_cell(cell!(r->&format!("{:.7e}",
            bin_norm * results.iter().fold(0.0, |value, results| value + results[bin]))));

        let mut normalization = 0.0;

        // calculate the sum of all leading orders
        for (index, order) in orders.iter().enumerate() {
            if (normalize.is_empty() && ((order.alphas + order.alpha) == lo_power))
                || (normalize.iter().any(|o| *o == (order.alphas, order.alpha)))
            {
                normalization += results[index][bin];
            }
        }

        // print each order normalized to the sum of all leading orders
        for result in results.iter().map(|vec| vec[bin]) {
            if absolute {
                row.add_cell(cell!(r->&format!("{:.7e}", result * bin_norm)));
            } else {
                row.add_cell(cell!(r->&format!("{:.2}%", result / normalization * 100.0)));
            }
        }
    }

    Ok(table)
}
