use super::helpers::{self, ConvoluteMode, GlobalConfiguration, Subcommand};
use anyhow::{bail, Result};
use clap::{Parser, ValueHint};
use prettytable::{cell, Row};
use std::collections::HashSet;
use std::path::PathBuf;

/// Compares the numerical content of two grids with each other.
#[derive(Parser)]
pub struct Opts {
    /// Path to the first grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input1: PathBuf,
    /// Path to the second grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input2: PathBuf,
    /// LHAPDF id or name of the PDF set.
    #[arg(value_parser = helpers::parse_pdfset)]
    pdfset: String,
    /// Ignore differences in the orders and sum them.
    #[arg(alias = "ignore_orders", long = "ignore-orders")]
    ignore_orders: bool,
    /// Ignore bin limits (but not number of bins).
    #[arg(long = "ignore-bin-limits")]
    ignore_bin_limits: bool,
    /// Ignore differences in the luminosity functions.
    #[arg(long = "ignore-lumis")]
    ignore_lumis: bool,
    /// Select orders of the first grid.
    #[arg(
        long,
        num_args(1..),
        value_delimiter = ',',
        value_parser = helpers::parse_order
    )]
    orders1: Vec<(u32, u32)>,
    /// Select orders of the second grid.
    #[arg(
        long,
        num_args(1..),
        value_delimiter = ',',
        value_parser = helpers::parse_order
    )]
    orders2: Vec<(u32, u32)>,
    /// Scale all results of the first grid.
    #[arg(long, default_value = "1.0")]
    scale1: f64,
    /// Scale all results of the second grid.
    #[arg(long, default_value = "1.0")]
    scale2: f64,
    /// Set the number of fractional digits shown for absolute numbers.
    #[arg(default_value_t = 7, long = "digits-abs", value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[arg(default_value_t = 3, long = "digits-rel", value_name = "REL")]
    digits_rel: usize,
}

impl Subcommand for Opts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<u8> {
        let grid1 = helpers::read_grid(&self.input1)?;
        let grid2 = helpers::read_grid(&self.input2)?;

        let orders1: HashSet<_> = grid1
            .orders()
            .iter()
            .filter(|order| {
                (order.logxir == 0)
                    && (order.logxif == 0)
                    && (self.orders1.is_empty()
                        || self
                            .orders1
                            .iter()
                            .any(|&o| (order.alphas, order.alpha) == o))
            })
            .collect();
        let orders2: HashSet<_> = grid2
            .orders()
            .iter()
            .filter(|order| {
                (order.logxir == 0)
                    && (order.logxif == 0)
                    && (self.orders2.is_empty()
                        || self
                            .orders2
                            .iter()
                            .any(|&o| (order.alphas, order.alpha) == o))
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

        if !self.ignore_bin_limits && (grid1.bin_info() != grid2.bin_info()) {
            bail!("bins limits differ");
        }

        if self.ignore_bin_limits && (grid1.bin_info().bins() != grid2.bin_info().bins()) {
            bail!("number of bins differ");
        }

        if !self.ignore_lumis && (grid1.lumi() != grid2.lumi()) {
            bail!("luminosities differ");
        }

        let mut pdf = helpers::create_pdf(&self.pdfset)?;

        let mut table = helpers::create_table();
        let mut title = Row::empty();
        title.add_cell(cell!(c->"b"));
        for i in 0..grid1.bin_info().dimensions() {
            let mut cell = cell!(c->format!("x{}", i + 1));
            cell.set_hspan(2);
            title.add_cell(cell);
        }

        let limits1 = helpers::convolute_limits(&grid1, &[], ConvoluteMode::Normal);

        if self.ignore_orders {
            let mut cell = cell!(c->"diff");
            cell.set_hspan(3);
            title.add_cell(cell);

            table.set_titles(title);

            let results1 = helpers::convolute(
                &grid1,
                &mut pdf,
                &orders1,
                &[],
                &[],
                1,
                ConvoluteMode::Normal,
                cfg.force_positive,
            );
            let results2 = helpers::convolute(
                &grid2,
                &mut pdf,
                &orders2,
                &[],
                &[],
                1,
                ConvoluteMode::Normal,
                cfg.force_positive,
            );

            for (bin, (limits1, (result1, result2))) in limits1
                .iter()
                .zip(results1.iter().zip(results2.iter()))
                .enumerate()
            {
                let row = table.add_empty_row();

                row.add_cell(cell!(r->bin));
                for (left, right) in limits1.iter() {
                    row.add_cell(cell!(r->format!("{left}")));
                    row.add_cell(cell!(r->format!("{right}")));
                }

                let result1 = result1 * self.scale1;
                let result2 = result2 * self.scale2;

                row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result1)));
                row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result2)));
                row.add_cell(cell!(r->format!("{:.*e}", self.digits_rel,
                if result1 == result2 { 0.0 } else { result1 / result2 - 1.0 })));
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
                    helpers::convolute(
                        &grid1,
                        &mut pdf,
                        &[order],
                        &[],
                        &[],
                        1,
                        ConvoluteMode::Normal,
                        cfg.force_positive,
                    )
                })
                .collect();
            let order_results2: Vec<Vec<f64>> = orders
                .iter()
                .map(|&order| {
                    helpers::convolute(
                        &grid2,
                        &mut pdf,
                        &[order],
                        &[],
                        &[],
                        1,
                        ConvoluteMode::Normal,
                        cfg.force_positive,
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
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result1)));
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result2)));
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_rel,
                    if result1 == result2 { 0.0 } else { result1 / result2 - 1.0 })));
                }
            }
        }

        table.printstd();

        Ok(0)
    }
}
