use super::helpers::{self, ConvFuns, ConvoluteMode};
use super::{GlobalConfiguration, Subcommand};
use anyhow::{anyhow, Result};
use clap::builder::{PossibleValuesParser, TypedValueParser};
use clap::{Parser, ValueHint};
use lhapdf::Pdf;
use pineappl::grid::Grid;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[cfg(feature = "applgrid")]
mod applgrid;
#[cfg(feature = "fastnlo")]
mod fastnlo;
#[cfg(feature = "fktable")]
mod fktable;

#[cfg(feature = "applgrid")]
fn convert_applgrid(
    input: &Path,
    alpha: u8,
    conv_funs: &mut [Pdf],
    _: usize,
) -> Result<(&'static str, Grid, Vec<f64>, usize)> {
    use pineappl_applgrid::ffi;

    // TODO: check AMCATNLO scale variations

    let mut grid = ffi::make_grid(input.to_str().unwrap())?;
    let pgrid = applgrid::convert_applgrid(grid.pin_mut(), alpha)?;
    let results = applgrid::convolve_applgrid(grid.pin_mut(), conv_funs);

    Ok(("APPLgrid", pgrid, results, 1))
}

#[cfg(not(feature = "applgrid"))]
fn convert_applgrid(
    _: &Path,
    _: u8,
    _: &mut [Pdf],
    _: usize,
) -> Result<(&'static str, Grid, Vec<f64>, usize)> {
    Err(anyhow!(
        "you need to install `pineappl` with feature `applgrid`"
    ))
}

#[cfg(feature = "fastnlo")]
fn convert_fastnlo(
    input: &Path,
    alpha: u8,
    conv_funs: &ConvFuns,
    member: usize,
    scales: usize,
    fnlo_mur: Option<&str>,
    fnlo_muf: Option<&str>,
) -> Result<(&'static str, Grid, Vec<f64>, usize)> {
    use pineappl_fastnlo::ffi;
    use std::ptr;

    // TODO: convert this into an error?
    assert_eq!(conv_funs.lhapdf_names.len(), 1);

    let mut file = ffi::make_fastnlo_lhapdf_with_name_file_set(
        input.to_str().unwrap(),
        &conv_funs.lhapdf_names[0],
        member.try_into().unwrap(),
    );

    {
        let mut reader = ffi::downcast_lhapdf_to_reader_mut(file.as_mut().unwrap());

        if let Some(mur) = fnlo_mur {
            reader.as_mut().SetMuRFunctionalForm(mur.parse()?);
        }
        if let Some(muf) = fnlo_muf {
            reader.as_mut().SetMuFFunctionalForm(muf.parse()?);
        }
    }

    let grid = fastnlo::convert_fastnlo_table(&file, alpha)?;
    let mut reader = ffi::downcast_lhapdf_to_reader_mut(file.as_mut().unwrap());

    // TODO: scale-variation log conversion is only enabled for flex grids
    let scales = if unsafe { reader.GetIsFlexibleScaleTable(ptr::null_mut()) } {
        scales
    } else {
        1
    };

    // fastNLO does not support a fragmentation scale
    let unpermuted_results: Vec<_> = helpers::SCALES_VECTOR_REN_FAC[0..scales]
        .iter()
        .map(|&(xir, xif, _)| {
            if !reader.as_mut().SetScaleFactorsMuRMuF(xir, xif) {
                return None;
            }
            reader.as_mut().CalcCrossSection();
            Some(ffi::GetCrossSection(reader.as_mut(), false))
        })
        .take_while(Option::is_some)
        .map(Option::unwrap)
        .collect();

    assert!(matches!(unpermuted_results.len(), 1 | 3 | 7 | 9));

    let bins = unpermuted_results[0].len();
    let actual_scales = unpermuted_results.len();

    let results: Vec<_> = (0..bins)
        .flat_map(|bin| unpermuted_results.iter().map(move |r| r[bin]))
        .collect();

    Ok(("fastNLO", grid, results, actual_scales))
}

#[cfg(not(feature = "fastnlo"))]
fn convert_fastnlo(
    _: &Path,
    _: u8,
    _: &ConvFuns,
    _: usize,
    _: usize,
    _: Option<&str>,
    _: Option<&str>,
) -> Result<(&'static str, Grid, Vec<f64>, usize)> {
    Err(anyhow!(
        "you need to install `pineappl` with feature `fastnlo`"
    ))
}

#[cfg(feature = "fktable")]
fn convert_fktable(input: &Path) -> Result<(&'static str, Grid, Vec<f64>, usize)> {
    let fktable = fktable::convert_fktable(input)?;

    Ok(("fktable", fktable, Vec::new(), 1))
}

#[cfg(not(feature = "fktable"))]
fn convert_fktable(_: &Path) -> Result<(&'static str, Grid, Vec<f64>, usize)> {
    Err(anyhow!(
        "you need to install `pineappl` with feature `fktable`"
    ))
}

fn convert_grid(
    input: &Path,
    alpha: u8,
    conv_funs: &mut [Pdf],
    fun_names: &ConvFuns,
    member: usize,
    scales: usize,
    fnlo_mur: Option<&str>,
    fnlo_muf: Option<&str>,
) -> Result<(&'static str, Grid, Vec<f64>, usize)> {
    if let Some(extension) = input.extension() {
        if extension == "tab"
            || (extension == "gz"
                && input
                    .with_extension("")
                    .extension()
                    .is_some_and(|ext| ext == "tab"))
        {
            return convert_fastnlo(input, alpha, fun_names, member, scales, fnlo_mur, fnlo_muf);
        } else if extension == "dat" {
            return convert_fktable(input);
        } else if extension == "appl" || extension == "root" {
            return convert_applgrid(input, alpha, conv_funs, scales);
        }
    }

    Err(anyhow!("could not detect file format"))
}

#[cfg(feature = "fastnlo")]
fn fnlo_mu_possible_values() -> Vec<&'static str> {
    vec![
        "kScale1",
        "kScale2",
        "kQuadraticSum",
        "kQuadraticMean",
        "kQuadraticSumOver4",
        "kLinearMean",
        "kLinearSum",
        "kScaleMax",
        "kScaleMin",
        "kProd",
        "kS2plusS1half",
        "kPow4Sum",
        "kWgtAvg",
        "kS2plusS1fourth",
        "kExpProd2",
    ]
}

#[cfg(not(feature = "fastnlo"))]
const fn fnlo_mu_possible_values() -> Vec<&'static str> {
    Vec::new()
}

/// Converts APPLgrid/fastNLO/FastKernel files to PineAPPL grids.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the converted grid.
    #[arg(value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// LHAPDF ID(s) or name of the PDF(s)/FF(s) to check the converted grid with.
    conv_funs: ConvFuns,
    /// LO coupling power in alpha.
    #[arg(default_value_t = 0, long)]
    alpha: u8,
    /// Relative threshold between the table and the converted grid when comparison fails.
    #[arg(default_value = "1e-10", long)]
    accuracy: f64,
    /// Set the number of scale variations to compare with if they are available.
    #[arg(
        default_value_t = 7,
        long,
        short,
        value_parser = PossibleValuesParser::new(["1", "3", "7", "9"]).try_map(|s| s.parse::<usize>())
    )]
    scales: usize,
    /// If importing a fastNLO flexible-scale grid, use the specified functional form for the
    /// renormalization scale.
    #[arg(long, value_parser = PossibleValuesParser::new(fnlo_mu_possible_values()))]
    fnlo_mur: Option<String>,
    /// If importing a fastNLO flexible-scale grid, use the specified functional form for the
    /// factorization scale.
    #[arg(long, value_parser = PossibleValuesParser::new(fnlo_mu_possible_values()))]
    fnlo_muf: Option<String>,
    /// Set the number of fractional digits shown for absolute numbers.
    #[arg(default_value_t = 7, long, value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[arg(default_value_t = 7, long, value_name = "REL")]
    digits_rel: usize,
    /// Do not optimize converted grid.
    #[arg(long)]
    no_optimize: bool,
}

impl Subcommand for Opts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        use prettytable::{cell, row};

        let mut conv_funs = helpers::create_conv_funs(&self.conv_funs)?;

        // TODO: figure out `member` from `self.pdfset`
        let (grid_type, mut grid, reference_results, scale_variations) = convert_grid(
            &self.input,
            self.alpha,
            &mut conv_funs,
            &self.conv_funs,
            0,
            self.scales,
            self.fnlo_mur.as_deref(),
            self.fnlo_muf.as_deref(),
        )?;

        if !self.no_optimize {
            grid.optimize();
        }

        let mut different = false;

        if reference_results.is_empty() {
            println!("file was converted, but we cannot check the conversion for this type");
        } else {
            let results = helpers::convolve(
                &grid,
                &mut conv_funs,
                &self.conv_funs.conv_types,
                &[],
                &[],
                &[],
                scale_variations,
                ConvoluteMode::Normal,
                cfg,
            );

            // if both grids don't have the same number of bins there's bug in the program
            assert_eq!(results.len(), reference_results.len());

            let mut table = helpers::create_table();
            let mut titles = row![c => "b", "PineAPPL", grid_type, "rel. diff"];

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
                    .zip(two)
                    .map(|(&a, &b)| {
                        // ALLOW: here we really need an exact comparison
                        // TODO: change allow to `expect` if MSRV >= 1.81.0
                        #[allow(clippy::float_cmp)]
                        if a == b {
                            0.0
                        } else {
                            b / a - 1.0
                        }
                    })
                    .collect();

                let mut row = row![
                    bin.to_string(),
                    r->format!("{:.*e}", self.digits_abs, one[0]),
                    r->format!("{:.*e}", self.digits_abs, two[0]),
                    r->format!("{:.*e}", self.digits_rel, rel_diffs[0])
                ];

                if rel_diffs[0].abs() > self.accuracy {
                    different = true;
                }

                if scale_variations > 1 {
                    // skip central scale choice
                    let &max_rel_diff = rel_diffs[1..]
                        .iter()
                        .max_by(|a, b| a.abs().total_cmp(&b.abs()))
                        // UNWRAP: in this branch we know there are scale variations
                        .unwrap();

                    if max_rel_diff.abs() > self.accuracy {
                        different = true;
                    }

                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_rel, max_rel_diff)));
                }

                table.add_row(row);
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
