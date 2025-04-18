use super::helpers::{self, ConvFuns, ConvoluteMode};
use super::{GlobalConfiguration, Subcommand};
use anyhow::{bail, Result};
use clap::{Parser, ValueHint};
use prettytable::{cell, Row};
use std::collections::HashSet;
use std::path::PathBuf;
use std::process::ExitCode;

/// Compares the numerical content of two grids with each other.
#[derive(Parser)]
pub struct Opts {
    /// Path to the first grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input1: PathBuf,
    /// Path to the second grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input2: PathBuf,
    /// LHAPDF ID(s) or name(s) of the PDF(s)/FF(s).
    conv_funs: ConvFuns,
    /// Ignore differences in the orders and sum them.
    #[arg(long)]
    ignore_orders: bool,
    /// Ignore bin limits (but not number of bins).
    #[arg(long)]
    ignore_bin_limits: bool,
    /// Ignore differences in the channel definition.
    #[arg(alias = "ignore-lumis", long)]
    ignore_channels: bool,
    /// Select orders of the first grid.
    #[arg(
        long,
        num_args = 1,
        value_delimiter = ',',
        value_parser = helpers::parse_order
    )]
    orders1: Vec<(u8, u8)>,
    /// Select orders of the second grid.
    #[arg(
        long,
        num_args = 1,
        value_delimiter = ',',
        value_parser = helpers::parse_order
    )]
    orders2: Vec<(u8, u8)>,
    /// Scale all results of the first grid.
    #[arg(long, default_value = "1.0")]
    scale1: f64,
    /// Scale all results of the second grid.
    #[arg(long, default_value = "1.0")]
    scale2: f64,
    /// Set the number of fractional digits shown for absolute numbers.
    #[arg(default_value_t = 7, long, value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[arg(default_value_t = 3, long, value_name = "REL")]
    digits_rel: usize,
}

impl Subcommand for Opts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        let grid1 = helpers::read_grid(&self.input1)?;
        let grid2 = helpers::read_grid(&self.input2)?;

        let orders1: HashSet<_> = grid1
            .orders()
            .iter()
            .filter(|order| {
                (order.logxir == 0)
                    && (order.logxif == 0)
                    && (self.orders1.is_empty()
                        || self.orders1.contains(&(order.alphas, order.alpha)))
            })
            .collect();
        let orders2: HashSet<_> = grid2
            .orders()
            .iter()
            .filter(|order| {
                (order.logxir == 0)
                    && (order.logxif == 0)
                    && (self.orders2.is_empty()
                        || self.orders2.contains(&(order.alphas, order.alpha)))
            })
            .collect();

        let mut diff1: Vec<_> = orders1.difference(&orders2).collect();
        diff1.sort();
        let diff1 = diff1;
        let mut diff2: Vec<_> = orders2.difference(&orders1).collect();
        diff2.sort();
        let diff2 = diff2;

        let mut orders1: Vec<_> = orders1
            .iter()
            .map(|order| (order.alphas, order.alpha))
            .collect();
        orders1.sort_unstable();
        let orders1 = orders1;
        let mut orders2: Vec<_> = orders2
            .iter()
            .map(|order| (order.alphas, order.alpha))
            .collect();
        orders2.sort_unstable();
        let orders2 = orders2;

        if !self.ignore_orders && (!diff1.is_empty() || !diff2.is_empty()) {
            bail!("selected orders differ");
        }

        if self.ignore_bin_limits && (grid1.bwfl().len() != grid2.bwfl().len()) {
            bail!("number of bins differ");
        }

        if !self.ignore_bin_limits && !grid1.bwfl().bins_partial_eq_with_ulps(grid2.bwfl(), 8) {
            bail!("bins limits differ");
        }

        // TODO: use approximate comparison
        if !self.ignore_channels && (grid1.channels() != grid2.channels()) {
            bail!("channels differ");
        }

        let mut conv_funs = helpers::create_conv_funs(&self.conv_funs)?;

        let mut table = helpers::create_table();
        let mut title = Row::empty();
        title.add_cell(cell!(c->"b"));
        for i in 0..grid1.bwfl().dimensions() {
            let mut cell = cell!(c->format!("x{}", i + 1));
            cell.set_hspan(2);
            title.add_cell(cell);
        }

        let limits1 = helpers::convolve_limits(&grid1, &[], ConvoluteMode::Normal);

        if self.ignore_orders {
            let mut cell = cell!(c->"diff");
            cell.set_hspan(3);
            title.add_cell(cell);

            table.set_titles(title);

            let results1 = helpers::convolve(
                &grid1,
                &mut conv_funs,
                &self.conv_funs.conv_types,
                &orders1,
                &[],
                &[],
                1,
                ConvoluteMode::Normal,
                cfg,
            );
            let results2 = helpers::convolve(
                &grid2,
                &mut conv_funs,
                &self.conv_funs.conv_types,
                &orders2,
                &[],
                &[],
                1,
                ConvoluteMode::Normal,
                cfg,
            );

            for (bin, (limits1, (result1, result2))) in limits1
                .iter()
                .zip(results1.iter().zip(results2.iter()))
                .enumerate()
            {
                let row = table.add_empty_row();

                row.add_cell(cell!(r->bin));
                for (left, right) in limits1 {
                    row.add_cell(cell!(r->format!("{left}")));
                    row.add_cell(cell!(r->format!("{right}")));
                }

                let result1 = result1 * self.scale1;
                let result2 = result2 * self.scale2;

                // ALLOW: here we really need an exact comparison
                // TODO: change allow to `expect` if MSRV >= 1.81.0
                #[allow(clippy::float_cmp)]
                let diff = if result1 == result2 {
                    0.0
                } else {
                    result2 / result1 - 1.0
                };

                row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result1)));
                row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result2)));
                row.add_cell(cell!(r->format!("{:.*e}", self.digits_rel, diff)));
            }
        } else {
            let orders = orders1;

            for (alphas, alpha) in &orders {
                let mut cell = cell!(c->format!("O(as^{alphas} a^{alpha})"));
                cell.set_hspan(3);
                title.add_cell(cell);
            }

            table.set_titles(title);

            let order_results1: Vec<Vec<f64>> = orders
                .iter()
                .map(|&order| {
                    helpers::convolve(
                        &grid1,
                        &mut conv_funs,
                        &self.conv_funs.conv_types,
                        &[order],
                        &[],
                        &[],
                        1,
                        ConvoluteMode::Normal,
                        cfg,
                    )
                })
                .collect();
            let order_results2: Vec<Vec<f64>> = orders
                .iter()
                .map(|&order| {
                    helpers::convolve(
                        &grid2,
                        &mut conv_funs,
                        &self.conv_funs.conv_types,
                        &[order],
                        &[],
                        &[],
                        1,
                        ConvoluteMode::Normal,
                        cfg,
                    )
                })
                .collect();

            for (bin, limits1) in limits1.iter().enumerate() {
                let row = table.add_empty_row();

                row.add_cell(cell!(r->bin));
                for (left, right) in limits1 {
                    row.add_cell(cell!(r->format!("{left}")));
                    row.add_cell(cell!(r->format!("{right}")));
                }

                for (result1, result2) in order_results1.iter().zip(order_results2.iter()) {
                    let result1 = result1[bin] * self.scale1;
                    let result2 = result2[bin] * self.scale2;

                    // ALLOW: here we really need an exact comparison
                    // TODO: change allow to `expect` if MSRV >= 1.81.0
                    #[allow(clippy::float_cmp)]
                    let diff = if result1 == result2 {
                        0.0
                    } else {
                        result2 / result1 - 1.0
                    };

                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result1)));
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result2)));
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_rel, diff)));
                }
            }
        }

        table.printstd();

        Ok(ExitCode::SUCCESS)
    }
}
