use super::helpers::{self, ConvoluteMode, GlobalConfiguration, Subcommand};
use anyhow::{anyhow, Result};
use clap::{Parser, ValueHint};
use lhapdf::Pdf;
use pineappl::fk_table::FkTable;
use pineappl::grid::Grid;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[cfg(feature = "evolve")]
mod eko {
    use anyhow::{bail, Result};
    use base64::Engine;
    use base64::engine::GeneralPurpose;
    use base64::engine::general_purpose::PAD;
    use base64::alphabet::URL_SAFE;
    use either::Either;
    use lz4_flex::frame::FrameDecoder;
    use ndarray::{Array4, Array5, Axis};
    use ndarray_npy::{NpzReader, ReadNpyExt};
    use pineappl::evolution::OperatorInfo;
    use pineappl::pids;
    use serde::Deserialize;
    use std::fs::File;
    use std::io::BufReader;
    use std::path::Path;
    use tar::Archive;

    #[derive(Deserialize)]
    struct MetadataV0 {
        #[serde(rename = "Q2grid")]
        q2_grid: Vec<f64>,
        inputgrid: Vec<f64>,
        inputpids: Vec<i32>,
        q2_ref: f64,
        targetgrid: Vec<f64>,
        targetpids: Vec<i32>,
    }

    #[derive(Deserialize)]
    struct Rotations {
        #[serde(alias = "_inputgrid")]
        inputgrid: Option<Vec<f64>>,
        #[serde(alias = "_inputpids", with = "either::serde_untagged_optional")]
        inputpids: Option<Either<Vec<Vec<f64>>, Vec<i32>>>,
        #[serde(alias = "_targetgrid")]
        targetgrid: Option<Vec<f64>>,
        #[serde(alias = "_targetpids")]
        targetpids: Option<Vec<i32>>,
        pids: Vec<i32>,
        xgrid: Vec<f64>,
    }

    #[derive(Deserialize)]
    struct MetadataV1 {
        mu20: f64,
        rotations: Rotations,
    }

    #[derive(Deserialize)]
    #[serde(untagged)]
    enum Metadata {
        V0(MetadataV0),
        V1(MetadataV1),
    }

    pub fn read(eko: &Path) -> Result<(OperatorInfo, Array5<f64>)> {
        let mut archive = Archive::new(File::open(eko)?);

        let mut metadata = None;
        let mut operator = None;
        let mut operators = Vec::new();
        let mut fac1 = Vec::new();

        let base64 = GeneralPurpose::new(&URL_SAFE, PAD);

        for entry in archive.entries()? {
            let file = entry?;
            let path = file.header().path()?;

            if let Some(file_name) = path.file_name() {
                // TODO: get rid of the unwrap
                match file_name.to_str().unwrap() {
                    "metadata.yaml" => metadata = Some(serde_yaml::from_reader(file)?),
                    "operators.npy.lz4" => {
                        operator = Some(Array5::<f64>::read_npy(FrameDecoder::new(
                            BufReader::new(file),
                        ))?);
                    }
                    x if x.ends_with(".npz.lz4") => {
                        let name = x.strip_suffix(".npz.lz4").unwrap();
                        let bytes = base64.decode(name.as_bytes())?;
                        let array: [u8; 8] = bytes.as_slice().try_into().unwrap();
                        let muf2 = f64::from_le_bytes(array);

                        let mut reader = BufReader::new(FrameDecoder::new(BufReader::new(file)));
                        let mut buffer = Vec::new();
                        std::io::copy(&mut reader, &mut buffer)?;
                        let mut npz = NpzReader::new(std::io::Cursor::new(buffer))?;
                        let operator: Array4<f64> = npz.by_name("operator.npy")?;

                        fac1.push(muf2);
                        operators.push(operator);
                    }
                    _ => {}
                }
            }
        }

        if !operators.is_empty() {
            let ops: Vec<_> = operators.iter().map(Array4::view).collect();
            operator = Some(ndarray::stack(Axis(0), &ops).unwrap());
        }

        let mut info = match metadata {
            Some(Metadata::V0(metadata)) => OperatorInfo {
                fac1: metadata.q2_grid.clone(),
                pids0: metadata.inputpids,
                x0: metadata.inputgrid,
                pids1: metadata.targetpids,
                x1: metadata.targetgrid,
                fac0: metadata.q2_ref,
                ren1: vec![],
                alphas: vec![],
                xir: 1.0,
                xif: 1.0,
                lumi_id_types: String::new(),
            },
            Some(Metadata::V1(metadata)) => OperatorInfo {
                fac1,
                pids0: metadata.rotations.inputpids.map_or_else(
                    || metadata.rotations.pids.clone(),
                    |either| {
                        either.right_or_else(|basis| {
                            basis
                                .into_iter()
                                .map(|factors| {
                                    let tuples: Vec<_> = metadata
                                        .rotations
                                        .pids
                                        .iter()
                                        .copied()
                                        .zip(factors.into_iter())
                                        .collect();

                                    pids::pdg_mc_ids_to_evol(&tuples).unwrap()
                                })
                                .collect()
                        })
                    },
                ),
                x0: metadata
                    .rotations
                    .inputgrid
                    .unwrap_or_else(|| metadata.rotations.xgrid.clone()),
                pids1: metadata
                    .rotations
                    .targetpids
                    .unwrap_or(metadata.rotations.pids),
                x1: metadata
                    .rotations
                    .targetgrid
                    .unwrap_or(metadata.rotations.xgrid),
                fac0: metadata.mu20,
                ren1: vec![],
                alphas: vec![],
                xir: 1.0,
                xif: 1.0,
                lumi_id_types: String::new(),
            },
            None => bail!("no `metadata.yaml` file found"),
        };

        info.lumi_id_types = pids::determine_lumi_id_types(&info.pids0);

        Ok((info, operator.unwrap()))
    }
}

#[cfg(feature = "evolve")]
fn evolve_grid(
    grid: &Grid,
    eko: &Path,
    pdf: &Pdf,
    orders: &[(u32, u32)],
    xir: f64,
    xif: f64,
) -> Result<FkTable> {
    use pineappl::subgrid::{Mu2, Subgrid};

    let (mut info, operator) = eko::read(eko)?;

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

    info.ren1 = ren1;
    info.alphas = alphas;
    info.xir = xir;
    info.xif = xif;

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

#[cfg(not(feature = "evolve"))]
fn evolve_grid(_: &Grid, _: &Path, _: &Pdf, _: &[(u32, u32)], _: f64, _: f64) -> Result<FkTable> {
    Err(anyhow!(
        "you need to install `pineappl` with feature `evolve`"
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
    #[arg(default_value_t = 7, long, value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[arg(default_value_t = 7, long, value_name = "REL")]
    digits_rel: usize,
    /// Select which orders to evolve.
    #[arg(
        long,
        num_args = 1,
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
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
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
            cfg.force_positive,
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
            cfg.force_positive,
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
