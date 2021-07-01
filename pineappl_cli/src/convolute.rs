use super::helpers;
use anyhow::Result;
use lhapdf::Pdf;
use prettytable::{cell, Row, Table};

pub fn subcommand(
    input: &str,
    pdfset: &str,
    other_pdfsets: &[&str],
    show_bins: &[usize],
    scales: usize,
    orders: &[(u32, u32)],
    absolute: bool,
    integrated: bool,
) -> Result<Table> {
    let grid = helpers::read_grid(input)?;
    let pdf = pdfset
        .parse()
        .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);

    let results = helpers::convolute(&grid, &pdf, orders, show_bins, &[], scales);

    let other_results: Vec<f64> = other_pdfsets
        .iter()
        .flat_map(|pdfset| {
            let pdf = pdfset
                .parse()
                .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);
            helpers::convolute(&grid, &pdf, &[], show_bins, &[], 1)
        })
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

    if absolute {
        for scale in &helpers::SCALES_VECTOR[0..scales] {
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

    let mut table = helpers::create_table();
    table.set_titles(title);

    for (index, values) in results.chunks_exact(scales).enumerate() {
        let min_value = values
            .iter()
            .min_by(|left, right| left.partial_cmp(right).unwrap())
            .unwrap();
        let max_value = values
            .iter()
            .max_by(|left, right| left.partial_cmp(right).unwrap())
            .unwrap();
        let bin = if show_bins.is_empty() {
            index
        } else {
            show_bins[index]
        };

        let row = table.add_empty_row();

        row.add_cell(cell!(r->&format!("{}", bin)));
        for (left, right) in left_limits.iter().zip(right_limits.iter()) {
            row.add_cell(cell!(r->&format!("{}", left[bin])));
            row.add_cell(cell!(r->&format!("{}", right[bin])));
        }
        row.add_cell(cell!(r->&format!("{:.7e}", if integrated { values[0] * normalizations[bin] } else { values[0] })));

        if absolute {
            for &value in values.iter() {
                row.add_cell(cell!(r->&format!("{:.7e}", if integrated { value * normalizations[bin] } else { value })));
            }
        } else {
            row.add_cell(cell!(r->&format!("{:.2}%", (min_value / values[0] - 1.0) * 100.0)));
            row.add_cell(cell!(r->&format!("{:.2}%", (max_value / values[0] - 1.0) * 100.0)));
        }

        let bins = if show_bins.is_empty() {
            bin_info.bins()
        } else {
            show_bins.len()
        };

        for &other in other_results.iter().skip(index).step_by(bins) {
            row.add_cell(cell!(r->&format!("{:.7e}", if integrated { other * normalizations[bin] } else { other })));
            row.add_cell(cell!(r->&format!("{:.2}%", (other / values[0] - 1.0) * 100.0)));
        }
    }

    Ok(table)
}
