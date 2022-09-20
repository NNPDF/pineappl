use super::helpers::{self, ConvoluteMode, Subcommand};
use anyhow::{anyhow, bail, Result};
use clap::{Parser, ValueHint};
use pineappl::grid::Grid;
use std::path::{Path, PathBuf};

#[cfg(feature = "applgrid")]
mod applgrid;
#[cfg(feature = "fastnlo")]
mod fastnlo;
#[cfg(feature = "fktable")]
mod fktable;

#[cfg(feature = "applgrid")]
fn convert_applgrid(
    input: &Path,
    alpha: u32,
    pdfset: &str,
    member: usize,
    dis_pid: i32,
    silence: bool,
) -> Result<(&'static str, Grid, Vec<f64>)> {
    use pineappl_applgrid::ffi;

    let mut grid = ffi::make_grid(input.to_str().unwrap(), silence)?;
    let pgrid = applgrid::convert_applgrid(grid.pin_mut(), alpha, dis_pid)?;
    let results = applgrid::convolute_applgrid(grid.pin_mut(), pdfset, member);

    Ok(("APPLgrid", pgrid, results))
}

#[cfg(not(feature = "applgrid"))]
fn convert_applgrid(
    _: &Path,
    _: u32,
    _: &str,
    _: usize,
    _: i32,
    _: bool,
) -> Result<(&'static str, Grid, Vec<f64>)> {
    Err(anyhow!(
        "you need to install `pineappl` with feature `applgrid`"
    ))
}

#[cfg(feature = "fastnlo")]
fn convert_fastnlo(
    input: &Path,
    alpha: u32,
    pdfset: &str,
    member: usize,
    dis_pid: i32,
    silence: bool,
) -> Result<(&'static str, Grid, Vec<f64>)> {
    use pineappl_fastnlo::ffi::{self, Verbosity};

    if silence {
        ffi::SetGlobalVerbosity(Verbosity::SILENT);
    }

    let mut file = ffi::make_fastnlo_lhapdf_with_name_file_set(
        input.to_str().unwrap(),
        pdfset,
        member.try_into().unwrap(),
        silence,
    );
    let grid = fastnlo::convert_fastnlo_table(&file, alpha, dis_pid)?;
    let results = ffi::GetCrossSection(
        ffi::static_cast_lhapdf_to_reader_mut(file.as_mut().unwrap()),
        false,
    );

    Ok(("fastNLO", grid, results))
}

#[cfg(not(feature = "fastnlo"))]
fn convert_fastnlo(
    _: &Path,
    _: u32,
    _: &str,
    _: usize,
    _: i32,
    _: bool,
) -> Result<(&'static str, Grid, Vec<f64>)> {
    Err(anyhow!(
        "you need to install `pineappl` with feature `fastnlo`"
    ))
}

#[cfg(feature = "fktable")]
fn convert_fktable(input: &Path) -> Result<(&'static str, Grid, Vec<f64>)> {
    let fktable = fktable::convert_fktable(input)?;

    Ok(("fktable", fktable, vec![]))
}

#[cfg(not(feature = "fktable"))]
fn convert_fktable(_: &Path) -> Result<(&'static str, Grid, Vec<f64>)> {
    Err(anyhow!(
        "you need to install `pineappl` with feature `fktable`"
    ))
}

fn convert_grid(
    input: &Path,
    alpha: u32,
    pdfset: &str,
    member: usize,
    dis_pid: i32,
    silence_libraries: bool,
) -> Result<(&'static str, Grid, Vec<f64>)> {
    if let Some(extension) = input.extension() {
        if extension == "tab"
            || (extension == "gz"
                && input
                    .with_extension("")
                    .extension()
                    .map_or(false, |ext| ext == "tab"))
        {
            return convert_fastnlo(input, alpha, pdfset, member, dis_pid, silence_libraries);
        } else if extension == "dat" {
            return convert_fktable(input);
        } else if extension == "appl" || extension == "root" {
            return convert_applgrid(input, alpha, pdfset, member, dis_pid, silence_libraries);
        }
    }

    Err(anyhow!("could not detect file format"))
}

/// Converts fastNLO/FastKernel tables to PineAPPL grids.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the converted grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// LHAPDF id or name of the PDF set to check the converted grid with.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset: String,
    /// LO coupling power in alpha.
    #[clap(default_value_t = 0, long)]
    alpha: u32,
    /// Relative threshold between the table and the converted grid when comparison fails.
    #[clap(default_value = "1e-10", long)]
    accuracy: f64,
    /// Prevents third-party libraries from printing output.
    #[clap(alias = "silence-fastnlo", long = "silence-libraries")]
    silence_libraries: bool,
    /// Set the number of fractional digits shown for absolute numbers.
    #[clap(default_value_t = 7, long = "digits-abs", value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[clap(default_value_t = 7, long = "digits-rel", value_name = "REL")]
    digits_rel: usize,
    /// Do not optimize converted grid.
    #[clap(long = "no-optimize")]
    no_optimize: bool,
    /// Particle ID for the non-hadronic initial states if it cannot be determined from the grid.
    #[clap(long = "dis-pid", default_value_t = 11)]
    dis_pid: i32,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        use prettytable::{cell, row};

        // TODO: figure out `member` from `self.pdfset`
        let (grid_type, mut grid, reference_results) = convert_grid(
            &self.input,
            self.alpha,
            &self.pdfset,
            0,
            self.dis_pid,
            self.silence_libraries,
        )?;

        if !self.no_optimize {
            grid.optimize();
        }

        let mut different = false;

        if reference_results.is_empty() {
            println!("file was converted, but we cannot check the conversion for this type");
        } else {
            let mut pdf = helpers::create_pdf(&self.pdfset)?;
            let results = helpers::convolute(
                &grid,
                &mut pdf,
                &[],
                &[],
                &[],
                1,
                ConvoluteMode::Normal,
                false,
            );

            if results.len() != reference_results.len() {
                bail!(
                    "grids have different number of bins: {} (input) vs {} (output)",
                    reference_results.len(),
                    results.len()
                );
            }

            let mut table = helpers::create_table();
            table.set_titles(row![c => "b", "PineAPPL", grid_type, "rel. diff"]);

            for (bin, (one, two)) in results
                .into_iter()
                .zip(reference_results.into_iter())
                .enumerate()
            {
                // catches the case where both results are zero
                let rel_diff = if one == two { 0.0 } else { two / one - 1.0 };

                if rel_diff.abs() > self.accuracy {
                    different = true;
                }

                table.add_row(row![
                    bin.to_string(),
                    r->format!("{:.*e}", self.digits_abs, one),
                    r->format!("{:.*e}", self.digits_abs, two),
                    r->format!("{:.*e}", self.digits_rel, rel_diff)
                ]);
            }

            table.printstd();
        }

        if different {
            Err(anyhow!("grids are different"))
        } else {
            helpers::write_grid(&self.output, &grid)
        }
    }
}
