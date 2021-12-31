use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use lhapdf::Pdf;
use prettytable::{cell, Row};
use std::path::PathBuf;

/// Shows the predictions for all bin for each order separately.
#[derive(Parser)]
#[clap(name = "orders")]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id or name of the PDF set.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset: String,
    /// Show absolute numbers of each perturbative order.
    #[clap(long, short)]
    absolute: bool,
    /// Show integrated numbers (without bin widths) instead of differential ones.
    #[clap(long, short)]
    integrated: bool,
    /// Normalize contributions to the specified orders.
    #[clap(
        conflicts_with = "absolute",
        long,
        min_values = 1,
        parse(try_from_str = helpers::parse_order),
        short,
        use_delimiter = true
    )]
    normalize: Vec<(u32, u32)>,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let grid = helpers::read_grid(&self.input)?;
        let pdf = self.pdfset.parse().map_or_else(
            |_| Pdf::with_setname_and_member(&self.pdfset, 0),
            Pdf::with_lhaid,
        );

        let mut orders: Vec<_> = grid
            .orders()
            .iter()
            .filter(|order| (order.logxir == 0) && (order.logxif == 0))
            .collect();
        orders.sort();
        let orders = orders;

        let results: Vec<Vec<f64>> = orders
            .iter()
            .map(|order| {
                helpers::convolute(&grid, &pdf, &[(order.alphas, order.alpha)], &[], &[], 1)
            })
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

        let labels = helpers::labels(&grid);
        let (y_label, x_labels) = labels.split_last().unwrap();
        let mut title = Row::empty();
        title.add_cell(cell!(c->"bin"));
        for x_label in x_labels {
            let mut cell = cell!(c->&x_label);
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        title.add_cell(cell!(c->if self.integrated { "integ" } else { y_label }));

        for order in &orders {
            title.add_cell(cell!(c->&format!("O(as^{} a^{})", order.alphas, order.alpha)));
        }

        let mut table = helpers::create_table();
        table.set_titles(title);

        for bin in 0..bin_info.bins() {
            let row = table.add_empty_row();
            let bin_norm = if self.integrated {
                normalizations[bin]
            } else {
                1.0
            };

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
                if (self.normalize.is_empty() && ((order.alphas + order.alpha) == lo_power))
                    || (self
                        .normalize
                        .iter()
                        .any(|o| *o == (order.alphas, order.alpha)))
                {
                    normalization += results[index][bin];
                }
            }

            // print each order normalized to the sum of all leading orders
            for result in results.iter().map(|vec| vec[bin]) {
                if self.absolute {
                    row.add_cell(cell!(r->&format!("{:.7e}", result * bin_norm)));
                } else {
                    row.add_cell(cell!(r->&format!("{:.2}%", result / normalization * 100.0)));
                }
            }
        }

        table.printstd();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-orders 
Shows the predictions for all bin for each order separately

USAGE:
    pineappl orders [OPTIONS] <INPUT> <PDFSET>

ARGS:
    <INPUT>     Path to the input grid
    <PDFSET>    LHAPDF id or name of the PDF set

OPTIONS:
    -a, --absolute                    Show absolute numbers of each perturbative order
    -h, --help                        Print help information
    -i, --integrated                  Show integrated numbers (without bin widths) instead of
                                      differential ones
    -n, --normalize <NORMALIZE>...    Normalize contributions to the specified orders
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["orders", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }
}
