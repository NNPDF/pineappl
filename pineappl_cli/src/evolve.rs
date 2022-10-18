use super::helpers::{self, ConvoluteMode, Subcommand};
use anyhow::{anyhow, Result};
use clap::{Parser, ValueHint};
use pineappl::fk_table::FkTable;
use pineappl::grid::Grid;
use std::path::{Path, PathBuf};

#[cfg(feature = "fktable")]
fn evolve_grid(grid: &Grid, eko: &Path, xir: f64, xif: f64) -> Result<FkTable> {
    use lz4_flex::frame::FrameDecoder;
    use ndarray::{Array1, Array5};
    use ndarray_npy::ReadNpyExt;
    use pineappl::evolution::OperatorInfo;
    use serde::Deserialize;
    use std::fs::File;
    use std::io::BufReader;
    use tar::Archive;

    #[derive(Default, Deserialize)]
    struct Metadata {
        #[serde(rename = "Q2grid")]
        q2_grid: Vec<f64>,
        #[allow(dead_code)]
        eko_version: String,
        inputgrid: Vec<f64>,
        inputpids: Vec<i32>,
        #[allow(dead_code)]
        interpolation_is_log: bool,
        #[allow(dead_code)]
        interpolation_polynomial_degree: usize,
        #[allow(dead_code)]
        interpolation_xgrid: Vec<f64>,
        q2_ref: f64,
        targetgrid: Vec<f64>,
        targetpids: Vec<i32>,
    }

    let mut archive = Archive::new(File::open(eko)?);

    let mut alphas = Vec::new();
    let mut operator = Default::default();
    let mut metadata: Metadata = Default::default();

    for entry in archive.entries()? {
        let file = entry?;
        let path = file.header().path()?;

        if let Some(file_name) = path.file_name() {
            // TODO: get rid of the unwrap
            match file_name.to_str().unwrap() {
                "metadata.yaml" => metadata = serde_yaml::from_reader(file)?,
                "operators.npy.lz4" => {
                    operator = Array5::<f64>::read_npy(FrameDecoder::new(BufReader::new(file)))?
                }
                "alphas.npy.lz4" => {
                    alphas =
                        Array1::<f64>::read_npy(FrameDecoder::new(BufReader::new(file)))?.to_vec()
                }
                _ => {}
            }
        }
    }

    // TODO: handle errors when files in the EKO are not present

    let info = OperatorInfo {
        fac1: metadata.q2_grid.clone(),
        pids0: metadata.inputpids,
        x0: metadata.inputgrid,
        pids1: metadata.targetpids,
        x1: metadata.targetgrid,
        fac0: metadata.q2_ref,
        ren1: metadata.q2_grid, // TODO: check whether this is true in the general case
        alphas,
        xir,
        xif,
        lumi_id_types: "pdg_mc_ids".to_string(), // TODO: determine this from the operator
    };

    Ok(grid.evolve(&operator, &info, &[])?)
}

#[cfg(not(feature = "fktable"))]
fn evolve_grid(_: &Grid, _: &Path, _: f64, _: f64) -> Result<FkTable> {
    Err(anyhow!(
        "you need to install `pineappl` with feature `fktable`"
    ))
}

/// Evolve a grid with an evolution kernel operator to an FK table.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the evolution kernel operator.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    eko: PathBuf,
    /// Path to the converted grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// LHAPDF id or name of the PDF set to check the converted grid with.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset: String,
    /// Relative threshold between the table and the converted grid when comparison fails.
    #[clap(default_value = "1e-3", long)]
    accuracy: f64,
    /// Set the number of fractional digits shown for absolute numbers.
    #[clap(default_value_t = 7, long = "digits-abs", value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[clap(default_value_t = 7, long = "digits-rel", value_name = "REL")]
    digits_rel: usize,
    /// Rescale the renormalization scale with this factor.
    #[clap(default_value_t = 1.0, long)]
    xir: f64,
    /// Rescale the factorization scale with this factor.
    #[clap(default_value_t = 1.0, long)]
    xif: f64,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        use prettytable::{cell, row};

        let grid = helpers::read_grid(&self.input)?;
        let fk_table = evolve_grid(&grid, &self.eko, self.xir, self.xif)?;

        let mut different = false;
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
        let evolved_results = helpers::convolute_scales(
            fk_table.grid(),
            &mut pdf,
            &[],
            &[],
            &[],
            &[(self.xir, self.xif)],
            ConvoluteMode::Normal,
            false,
        );

        // if both grids don't have the same number of bins there's a bug in the program
        assert_eq!(results.len(), evolved_results.len());

        let mut table = helpers::create_table();
        table.set_titles(row![c => "b", "FkTable", "Grid", "rel. diff"]);

        for (bin, (one, two)) in results
            .into_iter()
            .zip(evolved_results.into_iter())
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

        if different {
            Err(anyhow!("grids are different"))
        } else {
            helpers::write_grid(&self.output, fk_table.grid())
        }
    }
}
