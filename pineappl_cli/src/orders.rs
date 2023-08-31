use super::helpers::{self, ConvoluteMode};
use super::{GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use prettytable::{cell, Row};
use std::path::PathBuf;
use std::process::ExitCode;

/// Shows the predictions for all bin for each order separately.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id or name of the PDF set.
    #[arg(value_parser = helpers::parse_pdfset)]
    pdfset: String,
    /// Show absolute numbers of each perturbative order.
    #[arg(long, short)]
    absolute: bool,
    /// Show integrated numbers (without bin widths) instead of differential ones.
    #[arg(long, short)]
    integrated: bool,
    /// Normalize contributions to the specified orders.
    #[arg(
        conflicts_with = "absolute",
        long,
        num_args = 1,
        short,
        value_delimiter = ',',
        value_parser = helpers::parse_order
    )]
    normalize: Vec<(u32, u32)>,
    /// Set the number of fractional digits shown for absolute numbers.
    #[arg(default_value_t = 7, long, value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[arg(default_value_t = 2, long, value_name = "REL")]
    digits_rel: usize,
}

impl Subcommand for Opts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        let grid = helpers::read_grid(&self.input)?;
        let mut pdf = helpers::create_pdf(&self.pdfset)?;

        let mut orders: Vec<_> = grid
            .orders()
            .iter()
            .filter(|order| (order.logxir == 0) && (order.logxif == 0))
            .collect();
        orders.sort();
        let orders = orders;

        let limits = helpers::convolute_limits(
            &grid,
            &[],
            if self.integrated {
                ConvoluteMode::Integrated
            } else {
                ConvoluteMode::Normal
            },
        );
        let results: Vec<Vec<f64>> = orders
            .iter()
            .map(|order| {
                helpers::convolute(
                    &grid,
                    &mut pdf,
                    &[(order.alphas, order.alpha)],
                    &[],
                    &[],
                    1,
                    if self.integrated {
                        ConvoluteMode::Integrated
                    } else {
                        ConvoluteMode::Normal
                    },
                    cfg.force_positive,
                )
            })
            .collect();

        let lo_power = {
            let order = orders.first().unwrap();
            order.alphas + order.alpha
        };

        let (x, y_label, y_unit) = helpers::labels_and_units(&grid, self.integrated);
        let mut title = Row::empty();
        title.add_cell(cell!(c->"b"));
        for (x_label, x_unit) in x {
            let mut cell = cell!(c->format!("{x_label}\n[{x_unit}]"));
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        title.add_cell(cell!(c->format!("{y_label}\n[{y_unit}]")));

        for order in &orders {
            title.add_cell(cell!(c->format!("O(as^{} a^{})\n[{}]", order.alphas, order.alpha, if self.absolute { y_unit } else { "%" })));
        }

        let mut table = helpers::create_table();
        table.set_titles(title);

        for (bin, limits) in limits.iter().enumerate() {
            let row = table.add_empty_row();

            row.add_cell(cell!(r->format!("{bin}")));
            for (left, right) in limits {
                row.add_cell(cell!(r->format!("{left}")));
                row.add_cell(cell!(r->format!("{right}")));
            }
            row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs,
            results.iter().fold(0.0, |value, results| value + results[bin]))));

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
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result)));
                } else {
                    row.add_cell(
                        cell!(r->format!("{:.*}", self.digits_rel, result / normalization * 100.0)),
                    );
                }
            }
        }

        table.printstd();

        Ok(ExitCode::SUCCESS)
    }
}
