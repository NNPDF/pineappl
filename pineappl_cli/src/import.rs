use super::helpers::{self, Subcommand};
use anyhow::{anyhow, Result};
use clap::{Parser, ValueHint};
use pineappl::grid::Grid;
use std::path::{Path, PathBuf};

#[cfg(feature = "fastnlo")]
mod fastnlo;
#[cfg(feature = "fktable")]
mod fktable;

#[cfg(feature = "fastnlo")]
fn convert_fastnlo(
    input: &Path,
    alpha: u32,
    pdfset: &str,
    member: usize,
    silence_fastnlo: bool,
) -> Result<(&'static str, Grid, Vec<f64>)> {
    use pineappl_fastnlo::ffi::{self, Verbosity};
    use std::convert::TryInto;

    if silence_fastnlo {
        ffi::SetGlobalVerbosity(Verbosity::SILENT);
    }

    let mut file = ffi::make_fastnlo_lhapdf_with_name_file_set(
        input.to_str().unwrap(),
        pdfset,
        member.try_into().unwrap(),
        silence_fastnlo,
    );
    let grid = fastnlo::convert_fastnlo_table(&file, alpha.try_into().unwrap())?;

    let results = ffi::GetCrossSection(
        ffi::downcast_lhapdf_to_reader_mut(file.as_mut().unwrap()),
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

#[cfg(not(feature = "fastnlo"))]
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
    silence_fastnlo: bool,
) -> Result<(&'static str, Grid, Vec<f64>)> {
    if let Some(extension) = input.extension() {
        if extension == "tab"
            || (extension == "gz"
                && input
                    .with_extension("")
                    .extension()
                    .map_or(false, |ext| ext == "tab"))
        {
            return convert_fastnlo(input, alpha, pdfset, member, silence_fastnlo);
        } else if extension == "dat" {
            return convert_fktable(input);
        }
    }

    Err(anyhow!("could not detect file format"))
}

/// Converts fastNLO tables to PineAPPL grids.
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
    /// Prevents fastNLO from printing output.
    #[clap(long = "silence-fastnlo")]
    silence_fastnlo: bool,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        use lhapdf::Pdf;
        use prettytable::{cell, row};

        // TODO: figure out `member` from `self.pdfset`
        let (grid_type, grid, reference_results) = convert_grid(
            &self.input,
            self.alpha,
            &self.pdfset,
            0,
            self.silence_fastnlo,
        )?;

        let mut different = false;

        if reference_results.is_empty() {
            println!("can not check conversion for this type");
        } else {
            let pdf = self.pdfset.parse().map_or_else(
                |_| Pdf::with_setname_and_member(&self.pdfset, 0),
                Pdf::with_lhaid,
            );
            let results = helpers::convolute(&grid, &pdf, &[], &[], &[], 1);

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
                    different = false;
                }

                table.add_row(row![
                    bin.to_string(),
                    r->format!("{:.7e}", one),
                    r->format!("{:.7e}", two),
                    r->format!("{:.7e}", rel_diff)
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

#[cfg(test)]
mod tests {
    use assert_cmd::Command;
    use assert_fs::NamedTempFile;

    const HELP_STR: &str = "pineappl-import 
Converts fastNLO tables to PineAPPL grids

USAGE:
    pineappl import [OPTIONS] <INPUT> <OUTPUT> <PDFSET>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path to the converted grid
    <PDFSET>    LHAPDF id or name of the PDF set to check the converted grid with

OPTIONS:
        --accuracy <ACCURACY>    Relative threshold between the table and the converted grid when
                                 comparison fails [default: 1e-10]
        --alpha <ALPHA>          LO coupling power in alpha [default: 0]
    -h, --help                   Print help information
        --silence-fastnlo        Prevents fastNLO from printing output
";

    #[cfg(feature = "fastnlo")]
    const IMPORT_FIX_GRID_STR: &str = "b   PineAPPL     fastNLO      rel. diff
-+------------+------------+--------------
0 2.9158424e-4 2.9158424e-4 -2.9976022e-15
1 2.4657895e-4 2.4657895e-4 -2.8865799e-15
";

    #[cfg(feature = "fastnlo")]
    const IMPORT_FLEX_GRID_STR: &str = "b   PineAPPL     fastNLO      rel. diff
-+------------+------------+--------------
0  8.2754182e1  8.2754182e1 -1.3544721e-14
1  3.6097335e1  3.6097335e1 -6.8833828e-15
2  8.0048746e0  8.0048746e0  5.3290705e-15
3 9.4319392e-1 9.4319392e-1  5.5511151e-15
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["import", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    #[cfg(feature = "fastnlo")]
    fn import_fix_grid() {
        let output = NamedTempFile::new("converted1.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "import",
                "--silence-fastnlo",
                "data/NJetEvents_0-0-2.tab.gz",
                output.path().to_str().unwrap(),
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(IMPORT_FIX_GRID_STR);
    }

    #[test]
    #[cfg(feature = "fastnlo")]
    fn import_flex_grid() {
        let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "import",
                "--silence-fastnlo",
                "data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
                output.path().to_str().unwrap(),
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(IMPORT_FLEX_GRID_STR);
    }
}
