use super::helpers::{self, ConvoluteMode, GlobalConfiguration, Subcommand};
use anyhow::{anyhow, Result};
use clap::{Parser, ValueHint};
use lhapdf::Pdf;
use pineappl::fk_table::FkTable;
use pineappl::grid::Grid;
use std::path::{Path, PathBuf};

#[cfg(feature = "fktable")]
fn evolve_grid(
    grid: &Grid,
    eko: &Path,
    pdf: &Pdf,
    orders: &[(u32, u32)],
    xir: f64,
    xif: f64,
) -> Result<FkTable> {
    use lz4_flex::frame::FrameDecoder;
    use ndarray::Array5;
    use ndarray_npy::ReadNpyExt;
    use pineappl::evolution::OperatorInfo;
    use pineappl::subgrid::{Mu2, Subgrid};
    use serde::Deserialize;
    use std::fs::File;
    use std::io::BufReader;
    use tar::Archive;

    #[derive(Default, Deserialize)]
    struct Metadata {
        #[serde(rename = "Q2grid")]
        q2_grid: Vec<f64>,
        inputgrid: Vec<f64>,
        inputpids: Vec<i32>,
        q2_ref: f64,
        targetgrid: Vec<f64>,
        targetpids: Vec<i32>,
    }

    let mut archive = Archive::new(File::open(eko)?);

    let mut operator = Default::default();
    let mut metadata = Metadata::default();

    for entry in archive.entries()? {
        let file = entry?;
        let path = file.header().path()?;

        if let Some(file_name) = path.file_name() {
            // TODO: get rid of the unwrap
            match file_name.to_str().unwrap() {
                "metadata.yaml" => metadata = serde_yaml::from_reader(file)?,
                "operators.npy.lz4" => {
                    operator = Array5::<f64>::read_npy(FrameDecoder::new(BufReader::new(file)))?;
                }
                _ => {}
            }
        }
    }

    // TODO: handle errors when files in the EKO are not present

    // TODO: the following should probably be a method of `Grid`
    let mut ren1: Vec<_> = grid
        .subgrids()
        .iter()
        .flat_map(|subgrid| {
            subgrid
                .mu2_grid()
                .iter()
                .map(|Mu2 { ren, .. }| xir * xir * ren)
                .collect::<Vec<_>>()
        })
        .collect();
    ren1.sort_by(|a, b| a.partial_cmp(b).unwrap());
    ren1.dedup();
    let ren1 = ren1;
    let alphas: Vec<_> = ren1.iter().map(|&mur2| pdf.alphas_q2(mur2)).collect();

    let info = OperatorInfo {
        fac1: metadata.q2_grid.clone(),
        pids0: metadata.inputpids,
        x0: metadata.inputgrid,
        pids1: metadata.targetpids,
        x1: metadata.targetgrid,
        fac0: metadata.q2_ref,
        ren1,
        alphas,
        xir,
        xif,
        lumi_id_types: "pdg_mc_ids".to_string(), // TODO: determine this from the operator
    };

    let orders: Vec<_> = grid
        .orders()
        .iter()
        .map(|order| {
            orders.is_empty()
                || orders
                    .iter()
                    .any(|other| (order.alphas == other.0) && (order.alpha == other.1))
        })
        .collect();

    Ok(grid.evolve(operator.view(), &info, &orders)?)
}

#[cfg(not(feature = "fktable"))]
fn evolve_grid(_: &Grid, _: &Path, _: &Pdf, _: &[(u32, u32)], _: f64, _: f64) -> Result<FkTable> {
    Err(anyhow!(
        "you need to install `pineappl` with feature `fktable`"
    ))
}

/// Evolve a grid with an evolution kernel operator to an FK table.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the evolution kernel operator.
    #[arg(value_hint = ValueHint::FilePath)]
    eko: PathBuf,
    /// Path to the converted grid.
    #[arg(value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// LHAPDF id or name of the PDF set to check the converted grid with.
    #[arg(value_parser = helpers::parse_pdfset)]
    pdfset: String,
    /// Relative threshold between the table and the converted grid when comparison fails.
    #[arg(default_value = "1e-3", long)]
    accuracy: f64,
    /// Set the number of fractional digits shown for absolute numbers.
    #[arg(default_value_t = 7, long = "digits-abs", value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[arg(default_value_t = 7, long = "digits-rel", value_name = "REL")]
    digits_rel: usize,
    /// Select which orders to evolve.
    #[arg(
        long,
        num_args(1..),
        short,
        value_delimiter = ',',
        value_parser = helpers::parse_order
    )]
    orders: Vec<(u32, u32)>,
    /// Rescale the renormalization scale with this factor.
    #[arg(default_value_t = 1.0, long)]
    xir: f64,
    /// Rescale the factorization scale with this factor.
    #[arg(default_value_t = 1.0, long)]
    xif: f64,
}

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<u8> {
        use prettytable::row;

        let grid = helpers::read_grid(&self.input)?;
        let mut pdf = helpers::create_pdf(&self.pdfset)?;
        let results = helpers::convolute_scales(
            &grid,
            &mut pdf,
            &self.orders,
            &[],
            &[],
            &[(self.xir, self.xif)],
            ConvoluteMode::Normal,
            false,
        );

        let fk_table = evolve_grid(&grid, &self.eko, &pdf, &self.orders, self.xir, self.xif)?;
        let evolved_results = helpers::convolute_scales(
            fk_table.grid(),
            &mut pdf,
            &[],
            &[],
            &[],
            &[(1.0, 1.0)],
            ConvoluteMode::Normal,
            false,
        );

        // if both grids don't have the same number of bins there's a bug in the program
        assert_eq!(results.len(), evolved_results.len());

        let mut table = helpers::create_table();
        table.set_titles(row![c => "b", "Grid", "FkTable", "rel. diff"]);

        let mut different = false;

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
