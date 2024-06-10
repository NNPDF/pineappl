use super::helpers::{self, ConvoluteMode};
use super::{GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::builder::TypedValueParser;
use clap::{value_parser, Parser, ValueHint};
use prettytable::{cell, Row};
use std::ops::RangeInclusive;
use std::path::PathBuf;
use std::process::ExitCode;
use std::slice;

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
        default_value_t = 10,
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
        alias = "lumis",
        conflicts_with = "limit",
        long,
        num_args = 1,
        value_delimiter = ',',
        value_parser = helpers::parse_integer_range
    )]
    channels: Vec<RangeInclusive<usize>>,
    /// Select orders manually.
    #[arg(
        long,
        num_args = 1,
        short,
        value_delimiter = ',',
        value_parser = helpers::parse_order
    )]
    orders: Vec<(u32, u32)>,
    /// Do not sort the channels according to their size.
    #[arg(long)]
    dont_sort: bool,
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

        let mut channels: Vec<_> = self.channels.iter().cloned().flatten().collect();
        channels.sort_unstable();
        channels.dedup();
        let channels = channels;

        let limit = grid.channels().len().min(self.limit);
        let limit = if channels.is_empty() {
            limit
        } else {
            limit.min(channels.len())
        };
        let limits = helpers::convolve_limits(
            &grid,
            &[],
            if self.integrated {
                ConvoluteMode::Integrated
            } else {
                ConvoluteMode::Normal
            },
        );
        let results: Vec<_> = (0..grid.channels().len())
            .map(|channel| {
                let mut channel_mask = vec![false; grid.channels().len()];
                channel_mask[channel] = true;
                helpers::convolve(
                    &grid,
                    slice::from_mut(&mut pdf),
                    &self.orders,
                    &[],
                    &channel_mask,
                    1,
                    if self.integrated {
                        ConvoluteMode::Integrated
                    } else {
                        ConvoluteMode::Normal
                    },
                    cfg,
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
            title.add_cell(cell!(c->"c"));
            title.add_cell(
                cell!(c->&if self.absolute { format!("{y_label}\n[{y_unit}]") } else { "size\n[%]".to_owned() }),
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
                    .map(|(channel, vec)| (channel, vec[bin]))
                    .collect();

                if !self.dont_sort {
                    // sort using the absolute value in descending order
                    values.sort_unstable_by(|(_, left), (_, right)| {
                        right.abs().total_cmp(&left.abs())
                    });
                }

                for (channel, value) in values
                    .iter()
                    .filter(|(channel, _)| {
                        channels.is_empty() || channels.iter().any(|c| c == channel)
                    })
                    .take(limit)
                {
                    row.add_cell(cell!(r->format!("{channel}")));
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, value)));
                }
            } else {
                let sum: f64 = results.iter().map(|vec| vec[bin]).sum();
                let mut percentages: Vec<_> = results
                    .iter()
                    .enumerate()
                    .map(|(channel, vec)| (channel, vec[bin] / sum * 100.0))
                    .collect();

                if !self.dont_sort {
                    // sort using the absolute value in descending order
                    percentages.sort_unstable_by(|(_, left), (_, right)| {
                        right.abs().total_cmp(&left.abs())
                    });
                }

                for (channel, percentage) in percentages
                    .iter()
                    .filter(|(channel, _)| {
                        channels.is_empty() || channels.iter().any(|c| c == channel)
                    })
                    .take(limit)
                {
                    row.add_cell(cell!(r->format!("{channel}")));
                    row.add_cell(cell!(r->format!("{:.*}", self.digits_rel, percentage)));
                }
            }
        }

        table.printstd();

        Ok(ExitCode::SUCCESS)
    }
}
