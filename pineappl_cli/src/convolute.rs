use super::helpers::{self, ConvoluteMode, GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use prettytable::{cell, Row};
use std::ops::RangeInclusive;
use std::path::PathBuf;
use std::process::ExitCode;

/// Convolutes a PineAPPL grid with a PDF set.
#[derive(Parser)]
pub struct Opts {
    /// Path of the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id(s) or name of the PDF set(s).
    #[arg(required = true, value_parser = helpers::parse_pdfset)]
    pdfsets: Vec<String>,
    /// Selects a subset of bins.
    #[arg(
        long,
        short,
        num_args = 1,
        value_delimiter = ',',
        value_parser = helpers::parse_integer_range
    )]
    bins: Vec<RangeInclusive<usize>>,
    /// Show integrated numbers (without bin widths) instead of differential ones.
    #[arg(long, short)]
    integrated: bool,
    /// Select orders manually.
    #[arg(
        long,
        num_args = 1,
        short,
        value_delimiter = ',',
        value_parser = helpers::parse_order
    )]
    orders: Vec<(u32, u32)>,
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
        let mut pdf = helpers::create_pdf(&self.pdfsets[0])?;
        let bins: Vec<_> = self.bins.iter().cloned().flatten().collect();

        let results = helpers::convolute(
            &grid,
            &mut pdf,
            &self.orders,
            &bins,
            &[],
            1,
            if self.integrated {
                ConvoluteMode::Integrated
            } else {
                ConvoluteMode::Normal
            },
            cfg.force_positive,
        );
        let limits = helpers::convolute_limits(
            &grid,
            &bins,
            if self.integrated {
                ConvoluteMode::Integrated
            } else {
                ConvoluteMode::Normal
            },
        );
        let bin_count = limits.len();

        let other_results: Vec<_> = self.pdfsets[1..]
            .iter()
            .flat_map(|pdfset| {
                let mut pdf = helpers::create_pdf(pdfset).unwrap();
                helpers::convolute(
                    &grid,
                    &mut pdf,
                    &self.orders,
                    &bins,
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

        let (x, y_label, y_unit) = helpers::labels_and_units(&grid, self.integrated);
        let mut title = Row::empty();
        title.add_cell(cell!(c->"b"));
        for (x_label, x_unit) in x {
            let mut cell = cell!(c->format!("{x_label}\n[{x_unit}]"));
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        title.add_cell(cell!(c->format!("{y_label}\n[{y_unit}]")));

        for other in self.pdfsets[1..].iter().map(|pdf| helpers::pdf_label(pdf)) {
            let mut cell = cell!(c->format!("{other}\n[{y_unit}] [%]"));
            cell.set_hspan(2);
            title.add_cell(cell);
        }

        let mut table = helpers::create_table();
        table.set_titles(title);

        for (index, (limits, value)) in limits.into_iter().zip(results.iter()).enumerate() {
            let bin = if bins.is_empty() { index } else { bins[index] };

            let row = table.add_empty_row();

            row.add_cell(cell!(r->format!("{bin}")));
            for (left, right) in &limits {
                row.add_cell(cell!(r->format!("{left}")));
                row.add_cell(cell!(r->format!("{right}")));
            }
            row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, value)));

            for &other in other_results.iter().skip(index).step_by(bin_count) {
                row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, other)));
                row.add_cell(
                    cell!(r->format!("{:.*}", self.digits_rel, (other / value - 1.0) * 100.0)),
                );
            }
        }

        table.printstd();

        Ok(ExitCode::SUCCESS)
    }
}
