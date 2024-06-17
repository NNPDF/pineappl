use super::helpers::{self, ConvoluteMode, VecConvFun};
use super::{GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::builder::TypedValueParser;
use clap::{value_parser, Parser, ValueHint};
use prettytable::{cell, Row};
use std::path::PathBuf;
use std::process::ExitCode;

/// Perform various analyses with grids.
#[derive(Parser)]
pub struct Opts {
    #[clap(subcommand)]
    subcommand: SubcommandEnum,
}

impl Subcommand for Opts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        self.subcommand.run(cfg)
    }
}

#[derive(Parser)]
enum SubcommandEnum {
    Ckf(CkfOpts),
}

impl Subcommand for SubcommandEnum {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        match self {
            Self::Ckf(opts) => opts.run(cfg),
        }
    }
}

/// Compare K-factors with channel K factors (ckf).
#[derive(Parser)]
pub struct CkfOpts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF ID(s) or name(s) of the PDF(s)/FF(s).
    #[arg(value_parser = helpers::parse_conv_funs)]
    // TODO: it would be better to use `Vec<ConvFun>`, but this consumes all following arguments
    conv_funs: VecConvFun,
    /// Order defining the K factors.
    #[arg(value_parser = helpers::parse_order)]
    order: (u32, u32),
    /// Normalizing orders of the K factors.
    #[arg(value_delimiter = ',', value_parser = helpers::parse_order)]
    orders_den: Vec<(u32, u32)>,
    /// The maximum number of channels displayed.
    #[arg(
        default_value_t = 10,
        long,
        short,
        // TODO: see https://github.com/clap-rs/clap/issues/4253
        value_parser = value_parser!(u16).range(1..).map(usize::from)
    )]
    limit: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[arg(default_value_t = 2, long, value_name = "REL")]
    digits_rel: usize,
}

impl Subcommand for CkfOpts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        let grid = helpers::read_grid(&self.input)?;
        let mut conv_funs = helpers::create_conv_funs(&self.conv_funs.0)?;

        let orders_den = if self.orders_den.is_empty() {
            grid.orders()
                .iter()
                .filter_map(|order| {
                    ((order.alphas != self.order.0) && (order.alpha != self.order.1))
                        .then_some((order.alphas, order.alpha))
                })
                .collect()
        } else {
            self.orders_den.clone()
        };

        let limit = grid.channels().len().min(self.limit);
        let limits = helpers::convolve_limits(&grid, &[], ConvoluteMode::Normal);
        let results: Vec<_> = (0..grid.channels().len())
            .map(|lumi| {
                let mut lumi_mask = vec![false; grid.channels().len()];
                lumi_mask[lumi] = true;
                helpers::convolve(
                    &grid,
                    &mut conv_funs,
                    &[self.order],
                    &[],
                    &lumi_mask,
                    1,
                    ConvoluteMode::Normal,
                    cfg,
                )
            })
            .collect();
        let results_den: Vec<_> = (0..grid.channels().len())
            .map(|lumi| {
                let mut lumi_mask = vec![false; grid.channels().len()];
                lumi_mask[lumi] = true;
                helpers::convolve(
                    &grid,
                    &mut conv_funs,
                    &orders_den,
                    &[],
                    &lumi_mask,
                    1,
                    ConvoluteMode::Normal,
                    cfg,
                )
            })
            .collect();

        let (x, _, _) = helpers::labels_and_units(&grid, false);
        let mut title = Row::empty();
        title.add_cell(cell!(c->"b"));
        for (x_label, x_unit) in x {
            let mut cell = cell!(c->format!("{x_label}\n[{x_unit}]"));
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        title.add_cell(cell!(c->"bin-K"));
        for _ in 0..limit {
            title.add_cell(cell!(c->"c"));
            title.add_cell(cell!(c->"K"));
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

            let mut values: Vec<_> = results
                .iter()
                .zip(results_den.iter())
                .enumerate()
                .map(|(lumi, (vec, vec_den))| (lumi, vec[bin], vec_den[bin]))
                .collect();

            // sort using the absolute value in descending order
            values.sort_unstable_by(|(_, left, left_den), (_, right, right_den)| {
                (right + right_den)
                    .abs()
                    .total_cmp(&(left + left_den).abs())
            });

            let (total, total_den) = values
                .iter()
                .fold((0.0, 0.0), |(nom, den), (_, add_nom, add_den)| {
                    (nom + add_nom, den + add_den)
                });

            row.add_cell(
                cell!(r->format!("{:.*}", self.digits_rel, (total + total_den) / total_den)),
            );

            for (lumi, value, value_den) in values.into_iter().take(limit) {
                row.add_cell(cell!(r->format!("{lumi}")));

                let channel_k = if value_den == 0.0 {
                    if value == 0.0 {
                        0.0
                    } else {
                        f64::INFINITY.copysign(value)
                    }
                } else {
                    (value + value_den) / value_den
                };

                row.add_cell(cell!(r->format!("{:.*}", self.digits_rel, channel_k)));
            }
        }

        table.printstd();

        Ok(ExitCode::SUCCESS)
    }
}
