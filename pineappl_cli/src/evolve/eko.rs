use anyhow::{Result, anyhow};
use base64::Engine as _;
use base64::alphabet::URL_SAFE;
use base64::engine::GeneralPurpose;
use base64::engine::general_purpose::PAD;
use either::Either;
use lz4_flex::frame::FrameDecoder;
use ndarray::iter::AxisIter;
use ndarray::{Array4, Array5, Axis, CowArray, Ix4};
use ndarray_npy::{NpzReader, ReadNpyExt as _};
use pineappl::convolutions::ConvType;
use pineappl::evolution::OperatorSliceInfo;
use pineappl::pids::{self, PidBasis};
use serde::Deserialize;
use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufReader, Cursor, Read as _};
use std::iter::Zip;
use std::path::Path;
use std::slice::Iter;
use tar::{Archive, Entries};

const BASES_V1_DEFAULT_PIDS: [i32; 14] = [22, -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6];

#[derive(Deserialize)]
struct BasesV1 {
    inputgrid: Option<Vec<f64>>,
    inputpids: Option<Vec<Vec<f64>>>,
    targetgrid: Option<Vec<f64>>,
    targetpids: Option<Vec<i32>>,
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
    pub fn iter_mut(&mut self) -> EkoSlicesIter<'_> {
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

    pub fn new(eko_path: &Path) -> Result<Self> {
        let metadata = Self::read_metadata(eko_path)?;

        match metadata {
            Metadata::V0(v0) => Self::with_v0(v0, eko_path),
            Metadata::V1(v1) => Self::with_v1(v1, eko_path),
            Metadata::V2(v2) => Self::with_v2(v2, eko_path),
            Metadata::V3(v3) => Self::with_v3(v3, eko_path),
        }
    }

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
                let Some(file_stem) = path.with_extension("").file_stem().map(ToOwned::to_owned)
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
                conv_type: ConvType::new(operator.configs.polarized, operator.configs.time_like),
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
                conv_type: ConvType::new(operator.configs.polarized, operator.configs.time_like),
            },
            archive: Archive::new(File::open(eko_path)?),
        })
    }
}

impl<'a> IntoIterator for &'a mut EkoSlices {
    type IntoIter = EkoSlicesIter<'a>;
    type Item = Result<(OperatorSliceInfo, CowArray<'a, f64, Ix4>)>;

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

#[derive(Deserialize)]
#[serde(untagged)]
enum Metadata {
    V0(MetadataV0),
    V1(MetadataV1),
    V2(MetadataV2),
    V3(MetadataV3), // v0.15 - v????
}

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
struct MetadataV1 {
    mu20: f64,
    rotations: Rotations,
}

#[derive(Deserialize)]
struct MetadataV2 {
    bases: BasesV1,
}

#[derive(Deserialize)]
struct MetadataV3 {
    xgrid: Vec<f64>,
}

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
