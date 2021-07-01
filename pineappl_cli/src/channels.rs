use super::helpers;
use anyhow::Result;
use lhapdf::Pdf;
use prettytable::{cell, Row, Table};

pub fn subcommand(
    input: &str,
    pdfset: &str,
    limit: usize,
    orders: &[(u32, u32)],
    absolute: bool,
    lumis: &[usize],
    integrated: bool,
) -> Result<Table> {
    let grid = helpers::read_grid(input)?;
    let pdf = pdfset
        .parse()
        .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);
    let limit = if lumis.is_empty() {
        grid.lumi().len().min(limit)
    } else {
        lumis.iter().filter(|lumi| **lumi < lumis.len()).count()
    };

    let results: Vec<_> = (0..grid.lumi().len())
        .map(|lumi| {
            let mut lumi_mask = vec![false; grid.lumi().len()];
            lumi_mask[lumi] = true;
            helpers::convolute(&grid, &pdf, orders, &[], &lumi_mask, 1)
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
    for _ in 0..limit {
        title.add_cell(cell!(c->"lumi"));
        title.add_cell(
            cell!(c->if absolute { if integrated { "integ" } else { y_label } } else { "size" }),
        );
    }

    let mut table = helpers::create_table();
    table.set_titles(title);

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
                .map(|(lumi, vec)| {
                    (
                        lumi,
                        if integrated {
                            normalizations[bin] * vec[bin]
                        } else {
                            vec[bin]
                        },
                    )
                })
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
