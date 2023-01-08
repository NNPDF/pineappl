use super::helpers::{self, ConvoluteMode, GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::builder::TypedValueParser;
use clap::{value_parser, Parser, ValueHint};
use prettytable::{cell, Row};
use std::ops::RangeInclusive;
use std::path::PathBuf;
use std::process::ExitCode;

/// Shows the contribution for each partonic channel.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id or name of the PDF set.
    #[arg(value_parser = helpers::parse_pdfset)]
    pdfset: String,
    /// Show absolute numbers of each contribution.
    #[arg(long, short)]
    absolute: bool,
    /// The maximum number of channels displayed.
    #[arg(
        default_value = "10",
        long,
        short,
        // TODO: see https://github.com/clap-rs/clap/issues/4253
        value_parser = value_parser!(u16).range(1..).map(usize::from)
    )]
    limit: usize,
    /// Show integrated numbers (without bin widths) instead of differential ones.
    #[arg(long, requires = "absolute", short)]
    integrated: bool,
    /// Show only the listed channels.
    #[arg(
        conflicts_with = "limit",
        long,
        num_args(1..),
        value_delimiter = ',',
        value_parser = helpers::parse_integer_range
    )]
    lumis: Vec<RangeInclusive<usize>>,
    /// Select orders manually.
    #[arg(
        long,
        num_args(1..),
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
        let mut pdf = helpers::create_pdf(&self.pdfset)?;

        let mut lumis: Vec<_> = self.lumis.iter().cloned().flatten().collect();
        lumis.sort_unstable();
        lumis.dedup();
        let lumis = lumis;

        let limit = grid.lumi().len().min(self.limit);
        let limit = if lumis.is_empty() {
            limit
        } else {
            limit.min(lumis.len())
        };
        let limits = helpers::convolute_limits(
            &grid,
            &[],
            if self.integrated {
                ConvoluteMode::Integrated
            } else {
                ConvoluteMode::Normal
            },
        );
        let results: Vec<_> = (0..grid.lumi().len())
            .map(|lumi| {
                let mut lumi_mask = vec![false; grid.lumi().len()];
                lumi_mask[lumi] = true;
                helpers::convolute(
                    &grid,
                    &mut pdf,
                    &self.orders,
                    &[],
                    &lumi_mask,
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
        for _ in 0..limit {
            title.add_cell(cell!(c->"l"));
            title.add_cell(
                cell!(c->&if self.absolute { format!("{y_label}\n[{y_unit}]") } else { "size\n[%]".to_string() }),
            );
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

            if self.absolute {
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
                    row.add_cell(cell!(r->format!("{lumi}")));
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, value)));
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
                    row.add_cell(cell!(r->format!("{lumi}")));
                    row.add_cell(cell!(r->format!("{:.*}", self.digits_rel, percentage)));
                }
            }
        }

        table.printstd();

        Ok(ExitCode::SUCCESS)
    }
}
