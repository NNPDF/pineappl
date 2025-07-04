use super::helpers::{self, ConvFuns, ConvoluteMode};
use super::{GlobalConfiguration, Subcommand};
use anyhow::{anyhow, Result};
use clap::{Parser, ValueHint};
use lhapdf::Pdf;
use pineappl::fk_table::FkTable;
use pineappl::grid::Grid;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[cfg(feature = "evolve")]
mod eko {
    use anyhow::{anyhow, Result};
    use base64::alphabet::URL_SAFE;
    use base64::engine::general_purpose::PAD;
    use base64::engine::GeneralPurpose;
    use base64::Engine;
    use either::Either;
    use lz4_flex::frame::FrameDecoder;
    use ndarray::iter::AxisIter;
    use ndarray::{Array4, Array5, Axis, CowArray, Ix4};
    use ndarray_npy::{NpzReader, ReadNpyExt};
    use pineappl::convolutions::ConvType;
    use pineappl::evolution::OperatorSliceInfo;
    use pineappl::pids::{self, PidBasis};
    use serde::Deserialize;
    use std::collections::HashMap;
    use std::ffi::OsString;
    use std::fs::File;
    use std::io::{BufReader, Cursor, Read};
    use std::iter::Zip;
    use std::path::Path;
    use std::slice::Iter;
    use tar::{Archive, Entries};

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
        V2(MetadataV2),
        V3(MetadataV3), // v0.15 - v????
    }

    const BASES_V1_DEFAULT_PIDS: [i32; 14] = [22, -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6];

    #[derive(Deserialize)]
    struct OperatorConfigsV1 {
        polarized: bool,
        time_like: bool,
    }

    #[derive(Deserialize)]
    struct OperatorV1 {
        mu0: f64,
        configs: OperatorConfigsV1,
    }

    #[derive(Deserialize)]
    struct OperatorV2 {
        init: Vec<f64>,
        configs: OperatorConfigsV1,
    }

    #[derive(Deserialize)]
    struct OperatorInfoV1 {
        scale: f64,
    }

    #[derive(Deserialize)]
    struct BasesV1 {
        inputgrid: Option<Vec<f64>>,
        inputpids: Option<Vec<Vec<f64>>>,
        targetgrid: Option<Vec<f64>>,
        targetpids: Option<Vec<i32>>,
        xgrid: Vec<f64>,
    }

    #[derive(Deserialize)]
    struct MetadataV2 {
        bases: BasesV1,
    }

    #[derive(Deserialize)]
    struct MetadataV3 {
        xgrid: Vec<f64>,
    }

    pub enum EkoSlices {
        V0 {
            fac1: Vec<f64>,
            info: OperatorSliceInfo,
            operator: Array5<f64>,
        },
        // V1 and V3 are special cases of V2
        V2 {
            fac1: HashMap<OsString, f64>,
            info: OperatorSliceInfo,
            archive: Archive<File>,
        },
    }

    impl EkoSlices {
        /// Read the EKO at `eko_path` and return the contents of the `metadata.yaml` file
        /// deserialized into a [`Metadata`] object.
        fn read_metadata(eko_path: &Path) -> Result<Metadata> {
            for entry in Archive::new(File::open(eko_path)?).entries_with_seek()? {
                let entry = entry?;
                let path = entry.path()?;

                if path.ends_with("metadata.yaml") {
                    return Ok(serde_yaml::from_reader(entry)?);
                }
            }

            Err(anyhow!("no file 'metadata.yaml' in EKO archive found"))
        }

        pub fn new(eko_path: &Path) -> Result<Self> {
            let metadata = Self::read_metadata(eko_path)?;

            match metadata {
                Metadata::V0(v0) => Self::with_v0(v0, eko_path),
                Metadata::V1(v1) => Self::with_v1(v1, eko_path),
                Metadata::V2(v2) => Self::with_v2(v2, eko_path),
                Metadata::V3(v3) => Self::with_v3(v3, eko_path),
            }
        }

        fn with_v0(metadata: MetadataV0, eko_path: &Path) -> Result<Self> {
            let mut operator = None;

            for entry in Archive::new(File::open(eko_path)?).entries_with_seek()? {
                let entry = entry?;
                let path = entry.path()?;

                if path.ends_with("operators.npy.lz4") {
                    operator = Some(Array5::read_npy(FrameDecoder::new(BufReader::new(entry)))?);
                }
            }

            let operator =
                operator.ok_or_else(|| anyhow!("no file 'operator.yaml' in EKO archive found"))?;

            Ok(Self::V0 {
                fac1: metadata.q2_grid,
                info: OperatorSliceInfo {
                    pid_basis: PidBasis::guess(&metadata.inputpids),
                    fac0: metadata.q2_ref,
                    pids0: metadata.inputpids,
                    x0: metadata.inputgrid,
                    fac1: 0.0,
                    pids1: metadata.targetpids,
                    x1: metadata.targetgrid,
                    conv_type: ConvType::UnpolPDF,
                },
                operator,
            })
        }

        fn with_v1(metadata: MetadataV1, eko_path: &Path) -> Result<Self> {
            let mut fac1 = HashMap::new();
            let base64 = GeneralPurpose::new(&URL_SAFE, PAD);

            for entry in Archive::new(File::open(eko_path)?).entries_with_seek()? {
                let entry = entry?;
                let path = entry.path()?;

                if path.starts_with("./operators")
                    && (path.extension().is_some_and(|ext| ext == "lz4"))
                    && (path
                        .with_extension("")
                        .extension()
                        .is_some_and(|ext| ext == "npz"))
                {
                    let Some(file_stem) =
                        path.with_extension("").file_stem().map(ToOwned::to_owned)
                    else {
                        continue;
                    };

                    let bytes = base64.decode(file_stem.to_string_lossy().as_bytes())?;
                    // UNWRAP: we assume that the filenames represent exactly 8 bytes
                    let array: [u8; 8] = bytes.as_slice().try_into().unwrap();
                    let scale = f64::from_le_bytes(array);

                    fac1.insert(file_stem, scale);
                }
            }

            let pids0 = metadata.rotations.inputpids.map_or_else(
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
                                    .zip(factors)
                                    .collect();

                                // UNWRAP: we assume that an evolution basis is specified, if
                                // that's not the case we must make the algorithm more generic
                                pids::pdg_mc_ids_to_evol(&tuples).unwrap()
                            })
                            .collect()
                    })
                },
            );

            Ok(Self::V2 {
                fac1,
                info: OperatorSliceInfo {
                    pid_basis: PidBasis::guess(&pids0),
                    fac0: metadata.mu20,
                    pids0,
                    x0: metadata
                        .rotations
                        .inputgrid
                        .unwrap_or_else(|| metadata.rotations.xgrid.clone()),
                    fac1: 0.0,
                    pids1: metadata
                        .rotations
                        .targetpids
                        .unwrap_or(metadata.rotations.pids),
                    x1: metadata
                        .rotations
                        .targetgrid
                        .unwrap_or(metadata.rotations.xgrid),
                    conv_type: ConvType::UnpolPDF,
                },
                archive: Archive::new(File::open(eko_path)?),
            })
        }

        fn with_v2(metadata: MetadataV2, eko_path: &Path) -> Result<Self> {
            let mut fac1 = HashMap::new();
            let mut operator: Option<OperatorV1> = None;

            for entry in Archive::new(File::open(eko_path)?).entries_with_seek()? {
                let entry = entry?;
                let path = entry.path()?;

                if path.starts_with("./operators")
                    && (path.extension().is_some_and(|ext| ext == "yaml"))
                {
                    let Some(file_stem) = path.file_stem().map(ToOwned::to_owned) else {
                        continue;
                    };

                    let op_info: OperatorInfoV1 = serde_yaml::from_reader(entry)?;
                    fac1.insert(file_stem, op_info.scale);
                } else if path.as_os_str() == "./operator.yaml" {
                    operator = Some(serde_yaml::from_reader(entry)?);
                }
            }

            let operator =
                operator.ok_or_else(|| anyhow!("no file 'operator.yaml' in EKO archive found"))?;

            let pids0 = metadata.bases.inputpids.map_or_else(
                || BASES_V1_DEFAULT_PIDS.to_vec(),
                |basis| {
                    basis
                        .into_iter()
                        .map(|factors| {
                            let tuples: Vec<_> =
                                BASES_V1_DEFAULT_PIDS.iter().copied().zip(factors).collect();

                            // UNWRAP: we assume that an evolution basis is specified, if that's
                            // not the case we must make the algorithm more generic
                            pids::pdg_mc_ids_to_evol(&tuples).unwrap()
                        })
                        .collect()
                },
            );

            Ok(Self::V2 {
                fac1,
                info: OperatorSliceInfo {
                    pid_basis: PidBasis::guess(&pids0),
                    fac0: operator.mu0 * operator.mu0,
                    pids0,
                    x0: metadata
                        .bases
                        .inputgrid
                        .unwrap_or_else(|| metadata.bases.xgrid.clone()),
                    fac1: 0.0,
                    pids1: metadata
                        .bases
                        .targetpids
                        .unwrap_or_else(|| BASES_V1_DEFAULT_PIDS.to_vec()),
                    x1: metadata
                        .bases
                        .targetgrid
                        .unwrap_or_else(|| metadata.bases.xgrid.clone()),
                    conv_type: ConvType::new(
                        operator.configs.polarized,
                        operator.configs.time_like,
                    ),
                },
                archive: Archive::new(File::open(eko_path)?),
            })
        }

        fn with_v3(metadata: MetadataV3, eko_path: &Path) -> Result<Self> {
            let mut fac1 = HashMap::new();
            let mut operator: Option<OperatorV2> = None;

            for entry in Archive::new(File::open(eko_path)?).entries_with_seek()? {
                let entry = entry?;
                let path = entry.path()?;

                if path.starts_with("./operators")
                    && (path.extension().is_some_and(|ext| ext == "yaml"))
                {
                    let Some(file_stem) = path.file_stem().map(ToOwned::to_owned) else {
                        continue;
                    };

                    let op_info: OperatorInfoV1 = serde_yaml::from_reader(entry)?;
                    fac1.insert(file_stem, op_info.scale);
                } else if path.as_os_str() == "./operator.yaml" {
                    operator = Some(serde_yaml::from_reader(entry)?);
                }
            }

            let operator =
                operator.ok_or_else(|| anyhow!("no file 'operator.yaml' in EKO archive found"))?;

            Ok(Self::V2 {
                fac1,
                info: OperatorSliceInfo {
                    // NOTE: Since v0.15, EKOs are always in the flavour basis
                    pid_basis: PidBasis::Pdg,
                    fac0: operator.init[0] * operator.init[0],
                    pids0: BASES_V1_DEFAULT_PIDS.to_vec(),
                    x0: metadata.xgrid.clone(),
                    fac1: 0.0,
                    pids1: BASES_V1_DEFAULT_PIDS.to_vec(),
                    x1: metadata.xgrid,
                    conv_type: ConvType::new(
                        operator.configs.polarized,
                        operator.configs.time_like,
                    ),
                },
                archive: Archive::new(File::open(eko_path)?),
            })
        }

        pub fn iter_mut(&mut self) -> EkoSlicesIter {
            match self {
                Self::V0 {
                    fac1,
                    info,
                    operator,
                } => EkoSlicesIter::V0 {
                    info: info.clone(),
                    iter: fac1.iter().zip(operator.axis_iter(Axis(0))),
                },
                Self::V2 {
                    fac1,
                    info,
                    archive,
                } => {
                    EkoSlicesIter::V2 {
                        fac1: fac1.clone(),
                        info: info.clone(),
                        // UNWRAP: short of changing the return type of this method we can't
                        // propagate the error, so we must panic here
                        entries: archive.entries_with_seek().unwrap(),
                    }
                }
            }
        }
    }

    impl<'a> IntoIterator for &'a mut EkoSlices {
        type Item = Result<(OperatorSliceInfo, CowArray<'a, f64, Ix4>)>;
        type IntoIter = EkoSlicesIter<'a>;

        fn into_iter(self) -> Self::IntoIter {
            self.iter_mut()
        }
    }

    pub enum EkoSlicesIter<'a> {
        V0 {
            info: OperatorSliceInfo,
            iter: Zip<Iter<'a, f64>, AxisIter<'a, f64, Ix4>>,
        },
        V2 {
            fac1: HashMap<OsString, f64>,
            info: OperatorSliceInfo,
            entries: Entries<'a, File>,
        },
    }

    impl<'a> Iterator for EkoSlicesIter<'a> {
        type Item = Result<(OperatorSliceInfo, CowArray<'a, f64, Ix4>)>;

        fn next(&mut self) -> Option<Self::Item> {
            match self {
                Self::V0 { info, iter } => {
                    if let Some((fac1, operator)) = iter.next() {
                        let mut info = info.clone();
                        info.fac1 = *fac1;

                        Some(Ok((info, CowArray::from(operator))))
                    } else {
                        None
                    }
                }
                Self::V2 {
                    fac1,
                    info,
                    entries,
                } => {
                    let fun = || {
                        for entry in entries {
                            let entry = entry?;
                            let path = entry.path()?;

                            // here we're only interested in the operators themselves
                            if path.starts_with("./operators")
                                && (path.extension().is_some_and(|ext| ext == "lz4"))
                                && (path
                                    .with_extension("")
                                    .extension()
                                    .is_some_and(|ext| ext == "npz"))
                            {
                                let Some(file_stem) =
                                    path.with_extension("").file_stem().map(ToOwned::to_owned)
                                else {
                                    continue;
                                };

                                let mut reader = FrameDecoder::new(BufReader::new(entry));
                                let mut buffer = Vec::new();
                                let _ = reader.read_to_end(&mut buffer)?;
                                let mut npz = NpzReader::new(Cursor::new(buffer))?;
                                let operator: Array4<f64> = npz.by_name("operator.npy")?;

                                let mut info = info.clone();
                                info.fac1 = fac1.get(&file_stem).copied().ok_or_else(|| anyhow!("file '{}.yaml' not found, could not determine the operator's factorization scale", file_stem.to_string_lossy()))?;

                                return Ok(Some((info, CowArray::from(operator))));
                            }
                        }

                        Ok(None)
                    };

                    fun().transpose()
                }
            }
        }
    }
}

#[cfg(feature = "evolve")]
fn evolve_grid(
    grid: &Grid,
    ekos: &[&Path],
    use_alphas_from: &Pdf,
    orders: &[(u8, u8)],
    xir: f64,
    xif: f64,
    xia: f64,
) -> Result<FkTable> {
    use eko::EkoSlices;
    use pineappl::evolution::AlphasTable;

    let order_mask: Vec<_> = grid
        .orders()
        .iter()
        .map(|order| {
            orders.is_empty()
                || orders
                    .iter()
                    .any(|other| (order.alphas == other.0) && (order.alpha == other.1))
        })
        .collect();

    let mut eko_slices: Vec<_> = ekos
        .iter()
        .map(|eko| EkoSlices::new(eko))
        .collect::<Result<_, _>>()?;
    let eko_slices: Vec<_> = eko_slices.iter_mut().collect();
    let alphas_table = AlphasTable::from_grid(grid, xir, &|q2| use_alphas_from.alphas_q2(q2));

    Ok(grid.evolve(eko_slices, &order_mask, (xir, xif, xia), &alphas_table)?)
}

#[cfg(not(feature = "evolve"))]
fn evolve_grid(
    _: &Grid,
    _: &[&Path],
    _: &Pdf,
    _: &[(u8, u8)],
    _: f64,
    _: f64,
    _: f64,
) -> Result<FkTable> {
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
    /// Path to the evolution kernel operator(s).
    #[arg(value_name = "EKO1,...")]
    ekos: String,
    /// Path to the converted grid.
    #[arg(value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// LHAPDF ID(s) or name of the PDF(s)/FF(s).
    conv_funs: ConvFuns,
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
    orders: Vec<(u8, u8)>,
    /// Rescale the renormalization scale with this factor.
    #[arg(default_value_t = 1.0, long)]
    xir: f64,
    /// Rescale the factorization scale with this factor.
    #[arg(default_value_t = 1.0, long)]
    xif: f64,
    /// Rescale the fragmentation scale with this factor.
    #[arg(default_value_t = 1.0, long)]
    xia: f64,
}

impl Subcommand for Opts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        use prettytable::row;

        let grid = helpers::read_grid(&self.input)?;
        let mut conv_funs = helpers::create_conv_funs(&self.conv_funs)?;
        let results = helpers::convolve_scales(
            &grid,
            &mut conv_funs,
            &self.conv_funs.conv_types,
            &self.orders,
            &[],
            &[],
            &[(self.xir, self.xif, self.xia)],
            ConvoluteMode::Normal,
            cfg,
        );

        let eko_paths: Vec<_> = self.ekos.split(',').map(Path::new).collect();
        let fk_table = evolve_grid(
            &grid,
            &eko_paths,
            &conv_funs[cfg.use_alphas_from],
            &self.orders,
            self.xir,
            self.xif,
            self.xia,
        )?;

        let evolved_results = helpers::convolve_scales(
            fk_table.grid(),
            &mut conv_funs,
            &self.conv_funs.conv_types,
            &[],
            &[],
            &[],
            &[(1.0, 1.0, 1.0)],
            ConvoluteMode::Normal,
            cfg,
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
            // ALLOW: here we really need an exact comparison
            // TODO: change allow to `expect` if MSRV >= 1.81.0
            #[allow(clippy::float_cmp)]
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
