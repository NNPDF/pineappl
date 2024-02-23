use super::helpers::{self, ConvoluteMode};
use super::{GlobalConfiguration, Subcommand};
use anyhow::{anyhow, Result};
use clap::builder::{PossibleValuesParser, TypedValueParser};
use clap::{Parser, ValueHint};
use pineappl::grid::{Grid, Order};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[cfg(feature = "applgrid")]
mod applgrid;

#[cfg(feature = "applgrid")]
fn convert_into_applgrid(
    output: &Path,
    grid: &Grid,
    pdfset: &str,
    member: usize,
    _: usize,
    discard_non_matching_scales: bool,
) -> Result<(&'static str, Vec<f64>, usize, Vec<bool>)> {
    // TODO: check also scale-varied results

    let (mut applgrid, order_mask) =
        applgrid::convert_into_applgrid(grid, output, discard_non_matching_scales)?;
    let results = applgrid::convolute_applgrid(applgrid.pin_mut(), pdfset, member);

    Ok(("APPLgrid", results, 1, order_mask))
}

#[cfg(not(feature = "applgrid"))]
fn convert_into_applgrid(
    _: &Path,
    _: &Grid,
    _: &str,
    _: usize,
    _: usize,
    _: bool,
) -> Result<(&'static str, Vec<f64>, usize, Vec<bool>)> {
    Err(anyhow!(
        "you need to install `pineappl` with feature `applgrid`"
    ))
}

fn convert_into_grid(
    output: &Path,
    grid: &Grid,
    pdfset: &str,
    member: usize,
    scales: usize,
    discard_non_matching_scales: bool,
) -> Result<(&'static str, Vec<f64>, usize, Vec<bool>)> {
    if let Some(extension) = output.extension() {
        if extension == "appl" || extension == "root" {
            return convert_into_applgrid(
                output,
                grid,
                pdfset,
                member,
                scales,
                discard_non_matching_scales,
            );
        }
    }

    Err(anyhow!("could not detect file format"))
}

/// Converts PineAPPL grids to APPLgrid files.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the converted grid.
    #[arg(value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// LHAPDF id or name of the PDF set to check the converted grid with.
    #[arg(value_parser = helpers::parse_pdfset)]
    pdfset: String,
    /// Relative threshold between the table and the converted grid when comparison fails.
    #[arg(default_value = "1e-10", long)]
    accuracy: f64,
    /// Discard non-matching scales that would otherwise lead to panics.
    #[arg(long)]
    discard_non_matching_scales: bool,
    /// Set the number of scale variations to compare with if they are available.
    #[arg(
        default_value_t = 7,
        long,
        short,
        value_parser = PossibleValuesParser::new(["1", "3", "7", "9"]).try_map(|s| s.parse::<usize>())
    )]
    scales: usize,
    /// Set the number of fractional digits shown for absolute numbers.
    #[arg(default_value_t = 7, long, value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[arg(default_value_t = 7, long, value_name = "REL")]
    digits_rel: usize,
}

impl Subcommand for Opts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        use prettytable::{cell, row};

        let grid = helpers::read_grid(&self.input)?;

        // TODO: figure out `member` from `self.pdfset`
        let (grid_type, results, scale_variations, order_mask) = convert_into_grid(
            &self.output,
            &grid,
            &self.pdfset,
            0,
            self.scales,
            self.discard_non_matching_scales,
        )?;

        for Order {
            alphas,
            alpha,
            logxir,
            logxif,
        } in grid
            .orders()
            .iter()
            .zip(order_mask.iter())
            .filter_map(|(order, keep)| (!keep).then_some(order.clone()))
        {
            println!("WARNING: the order O(as^{alphas} a^{alpha} lr^{logxir} lf^{logxif}) isn't supported by {grid_type} and will be skipped.");
        }

        let orders: Vec<_> = grid
            .orders()
            .iter()
            .zip(order_mask)
            .filter_map(
                |(
                    &Order {
                        alphas,
                        alpha,
                        logxir,
                        logxif,
                    },
                    keep,
                )| {
                    (keep && (logxir == 0) && (logxif == 0)).then_some((alphas, alpha))
                },
            )
            .collect();

        let mut different = false;

        if results.is_empty() {
            println!("file was converted, but we cannot check the conversion for this type");
        } else {
            let mut pdf = helpers::create_pdf(&self.pdfset)?;
            let reference_results = helpers::convolute(
                &grid,
                &mut pdf,
                &orders,
                &[],
                &[],
                scale_variations,
                ConvoluteMode::Normal,
                cfg,
            );

            // if both grids don't have the same number of bins there's bug in the program
            assert_eq!(results.len(), reference_results.len());

            let mut table = helpers::create_table();
            let mut titles = row![c => "b", grid_type, "PineAPPL", "rel. diff"];

            if scale_variations > 1 {
                titles.add_cell(cell!(c -> "svmaxreldiff"));
            }

            table.set_titles(titles);

            for (bin, (one, two)) in results
                .chunks_exact(scale_variations)
                .zip(reference_results.chunks_exact(scale_variations))
                .enumerate()
            {
                // catches the case where both results are zero
                let rel_diffs: Vec<_> = one
                    .iter()
                    .zip(two.iter())
                    .map(|(a, b)| if a == b { 0.0 } else { b / a - 1.0 })
                    .collect();

                let max_rel_diff = rel_diffs
                    .iter()
                    .max_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap())
                    .unwrap()
                    .abs();

                if max_rel_diff > self.accuracy {
                    different = true;
                }

                let mut row = row![
                    bin.to_string(),
                    r->format!("{:.*e}", self.digits_abs, one[0]),
                    r->format!("{:.*e}", self.digits_abs, two[0]),
                    r->format!("{:.*e}", self.digits_rel, rel_diffs[0])
                ];

                if scale_variations > 1 {
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_rel, max_rel_diff)));
                }

                table.add_row(row);
            }

            table.printstd();
        }

        if different {
            Err(anyhow!("grids are different"))
        } else {
            Ok(ExitCode::SUCCESS)
        }
    }
}
