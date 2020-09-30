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
    other_pdfsets: &[&str],
    show_bins: &[usize],
    scales: usize,
    orders: &[(u32, u32)],
    absolute: bool,
) -> Result<Table, Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;
    let show_bins = if show_bins.is_empty() {
        (0..grid.bin_limits().bins()).collect()
    } else {
        show_bins.to_vec()
    };
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

    let results = grid.convolute(
        &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
        &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
        &|q2| pdf.alphas_q2(q2),
        &orders,
        &show_bins,
        &[],
        &scales_vector[0..scales],
    );

    let other_results: Vec<f64> = other_pdfsets
        .iter()
        .flat_map(|pdfset| {
            let pdf = pdfset
                .parse()
                .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);
            grid.convolute(
                &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
                &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
                &|q2| pdf.alphas_q2(q2),
                &[],
                &show_bins,
                &[],
                &[(1.0, 1.0)],
            )
        })
        .collect();

    let bin_sizes = grid.bin_limits().bin_sizes();
    let bin_limits = grid.bin_limits().limits();

    let mut table = create_table();
    let mut title = Row::empty();
    title.add_cell(cell!(c->"bin"));
    title.add_cell(cell!(c->"xmin"));
    title.add_cell(cell!(c->"xmax"));
    title.add_cell(cell!(c->"diff"));
    title.add_cell(cell!(c->"integ"));

    if absolute {
        for scale in &scales_vector[0..scales] {
            title.add_cell(cell!(c->&format!("({},{})", scale.0, scale.1)));
        }
    } else {
        title.add_cell(cell!(c->"neg unc"));
        title.add_cell(cell!(c->"pos unc"));
    }

    for other in other_pdfsets.iter() {
        let mut cell = cell!(c->other);
        cell.set_hspan(2);
        title.add_cell(cell);
    }

    table.set_titles(title);

    for (bin, values) in results.chunks_exact(scales).enumerate() {
        let min_value = values
            .iter()
            .min_by(|left, right| left.partial_cmp(right).unwrap())
            .unwrap();
        let max_value = values
            .iter()
            .max_by(|left, right| left.partial_cmp(right).unwrap())
            .unwrap();

        let row = table.add_empty_row();

        row.add_cell(cell!(r->&format!("{}", show_bins[bin])));
        row.add_cell(cell!(r->&format!("{}", bin_limits[show_bins[bin]])));
        row.add_cell(cell!(r->&format!("{}", bin_limits[show_bins[bin] + 1])));
        row.add_cell(cell!(r->&format!("{:.7e}", values[0])));
        row.add_cell(cell!(r->&format!("{:.7e}", values[0] * bin_sizes[show_bins[bin]])));

        if absolute {
            for value in values.iter() {
                row.add_cell(cell!(r->&format!("{:.7e}", value * bin_sizes[show_bins[bin]])));
            }
        } else {
            row.add_cell(cell!(r->&format!("{:.2}%", (min_value / values[0] - 1.0) * 100.0)));
            row.add_cell(cell!(r->&format!("{:.2}%", (max_value / values[0] - 1.0) * 100.0)));
        }

        for other in other_results.iter().skip(bin).step_by(show_bins.len()) {
            row.add_cell(cell!(r->&format!("{:.7e}", other)));
            row.add_cell(cell!(r->&format!("{:.2}%", (other / values[0] - 1.0) * 100.0)));
        }
    }

    Ok(table)
}
