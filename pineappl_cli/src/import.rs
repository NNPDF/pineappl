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
) -> Result<(&'static str, Grid, Vec<f64>)> {
    use pineappl_applgrid::ffi;

    let mut grid = ffi::make_grid(input.to_str().unwrap())?;
    let pgrid = applgrid::convert_applgrid(grid.pin_mut(), alpha)?;
    let results = applgrid::convolute_applgrid(grid.pin_mut(), pdfset, member);

    Ok(("APPLgrid", pgrid, results))
}

#[cfg(not(feature = "applgrid"))]
fn convert_applgrid(_: &Path, _: u32, _: &str, _: usize) -> Result<(&'static str, Grid, Vec<f64>)> {
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
    silence_fastnlo: bool,
) -> Result<(&'static str, Grid, Vec<f64>)> {
    use pineappl_fastnlo::ffi::{self, Verbosity};

    if silence_fastnlo {
        ffi::SetGlobalVerbosity(Verbosity::SILENT);
    }

    let mut file = ffi::make_fastnlo_lhapdf_with_name_file_set(
        input.to_str().unwrap(),
        pdfset,
        member.try_into().unwrap(),
        silence_fastnlo,
    );
    let grid = fastnlo::convert_fastnlo_table(&file, alpha)?;

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
        } else if extension == "appl" || extension == "root" {
            return convert_applgrid(input, alpha, pdfset, member);
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
    /// Prevents fastNLO from printing output.
    #[clap(long = "silence-fastnlo")]
    silence_fastnlo: bool,
    /// Set the number of fractional digits shown for absolute numbers.
    #[clap(default_value_t = 7, long = "digits-abs", value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[clap(default_value_t = 7, long = "digits-rel", value_name = "REL")]
    digits_rel: usize,
    /// Do not optimize converted grid.
    #[clap(long = "no-optimize")]
    no_optimize: bool,
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
            self.silence_fastnlo,
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

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    #[cfg(any(feature = "applgrid", feature = "fastnlo", feature = "fktable"))]
    use assert_fs::NamedTempFile;

    const HELP_STR: &str = "pineappl-import 
Converts fastNLO/FastKernel tables to PineAPPL grids

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
        --digits-abs <ABS>       Set the number of fractional digits shown for absolute numbers
                                 [default: 7]
        --digits-rel <REL>       Set the number of fractional digits shown for relative numbers
                                 [default: 7]
    -h, --help                   Print help information
        --no-optimize            Do not optimize converted grid
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

    #[cfg(feature = "fktable")]
    const IMPORT_DIS_FKTABLE_STR: &str = "b   x1       diff      scale uncertainty
    []        []              [%]       
--+--+--+-------------+--------+--------
 0  0  1   1.8900584e0   -69.81   107.21
 1  1  2   1.4830883e0   -69.63    98.20
 2  2  3   1.1497012e0   -69.68    90.39
 3  3  4  8.7974506e-1   -69.38    82.45
 4  4  5  7.0882550e-1   -69.14    70.54
 5  5  6  5.7345845e-1   -69.02    59.25
 6  6  7  4.7744833e-1   -68.37    44.68
 7  7  8  4.1037984e-1   -67.36    29.06
 8  8  9  4.0362470e-1   -65.97    12.72
 9  9 10  4.2613006e-1   -64.37     0.00
10 10 11  3.7669466e-1   -63.54     0.00
11 11 12  2.9572989e-1   -62.91     0.00
12 12 13  2.0869778e-1   -62.28     0.00
13 13 14  1.2602675e-1   -61.64     0.00
14 14 15  6.4220769e-2   -60.94     0.00
15 15 16  2.5434367e-2   -60.76     0.00
16 16 17  7.6070428e-3   -59.97     0.00
17 17 18  2.1848546e-3   -60.65     0.00
18 18 19  6.2309735e-4   -57.15     0.00
19 19 20 -1.0496472e-4     0.00   -55.42
";

    #[cfg(feature = "fktable")]
    const IMPORT_HADRONIC_FKTABLE_STR: &str = "b x1     diff     scale uncertainty
  []      []             [%]       
-+-+-+-----------+--------+--------
0 0 1 7.7624461e2   -86.97     0.00
";

    #[cfg(feature = "applgrid")]
    const IMPORT_PHOTON_GRID_STR: &str = "b   PineAPPL     APPLgrid     rel. diff
-+------------+------------+--------------
0 5.5621307e-4 5.5621307e-4 -1.5543122e-15
";

    #[cfg(feature = "applgrid")]
    const IMPORT_APPLGRID_STR: &str = "b  PineAPPL    APPLgrid     rel. diff
-+-----------+-----------+--------------
0 2.9884537e6 2.9884537e6 -6.6613381e-16
";

    #[cfg(feature = "applgrid")]
    const IMPORT_NEW_APPLGRID_STR: &str = "b   PineAPPL    APPLgrid     rel. diff
--+-----------+-----------+--------------
0  6.2634897e2 6.2634897e2  1.5543122e-15
1  6.2847078e2 6.2847078e2  4.4408921e-16
2  6.3163323e2 6.3163323e2  6.6613381e-16
3  6.3586556e2 6.3586556e2  4.4408921e-16
4  6.4139163e2 6.4139163e2  1.3322676e-15
5  6.4848088e2 6.4848088e2  2.2204460e-16
6  6.5354150e2 6.5354150e2 -4.2188475e-15
7  6.5377566e2 6.5377566e2  6.6613381e-16
8  6.5094729e2 6.5094729e2  8.8817842e-16
9  6.3588760e2 6.3588760e2 -3.8857806e-15
10 5.9810718e2 5.9810718e2  2.4424907e-15
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

    #[test]
    #[cfg(feature = "fktable")]
    fn import_dis_fktable() {
        let output = NamedTempFile::new("converted3.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "import",
                "data/FK_POSXDQ.dat",
                output.path().to_str().unwrap(),
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout("file was converted, but we cannot check the conversion for this type\n");

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                output.path().to_str().unwrap(),
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(IMPORT_DIS_FKTABLE_STR);
    }

    #[test]
    #[cfg(feature = "fktable")]
    fn import_hadronic_fktable() {
        let output = NamedTempFile::new("converted4.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "import",
                "data/FK_ATLASTTBARTOT13TEV.dat",
                output.path().to_str().unwrap(),
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout("file was converted, but we cannot check the conversion for this type\n");

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                output.path().to_str().unwrap(),
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(IMPORT_HADRONIC_FKTABLE_STR);
    }

    #[test]
    #[cfg(feature = "applgrid")]
    fn import_photon_grid() {
        let output = NamedTempFile::new("converted5.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "import",
                "data/LHCBWZMU7TEV_PI_part1.root",
                output.path().to_str().unwrap(),
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(predicates::str::ends_with(IMPORT_PHOTON_GRID_STR));
    }

    #[test]
    #[cfg(feature = "applgrid")]
    fn import_applgrid() {
        let output = NamedTempFile::new("converted6.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "import",
                "data/ATLASWPT11-Wplus_tot.root",
                output.path().to_str().unwrap(),
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(predicates::str::ends_with(IMPORT_APPLGRID_STR));
    }

    #[test]
    #[cfg(feature = "applgrid")]
    fn import_new_applgrid() {
        let output = NamedTempFile::new("converted7.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "import",
                "data/atlas-atlas-wpm-arxiv-1109.5141-xsec001.appl",
                output.path().to_str().unwrap(),
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(predicates::str::ends_with(IMPORT_NEW_APPLGRID_STR));
    }
}
