//! Module containing all traits and supporting structures for grids.

use super::bin::{BinInfo, BinLimits, BinRemapper};
use super::boc::{Channel, Order};
use super::convolutions::Convolution;
use super::pids::PidBasis;
use super::subgrid::{SubgridEnum, SubgridParams};
use ndarray::{Array3, ArrayView3};
use serde::Deserialize;
use std::collections::HashMap;
use std::io::{self, BufRead};
use thiserror::Error;

/// This structure represents a position (`x1`, `x2`, `q2`) in a `Subgrid` together with a
/// corresponding `weight`. The type `W` can either be a `f64` or `()`, which is used when multiple
/// weights should be signaled.
#[derive(Deserialize)]
pub struct Ntuple<W> {
    /// Momentum fraction of the first parton.
    pub x1: f64,
    /// Momentum fraction of the second parton.
    pub x2: f64,
    /// Squared scale.
    pub q2: f64,
    /// Weight of this entry.
    pub weight: W,
}

/// Error returned when merging two grids fails.
#[derive(Debug, Error)]
pub enum GridError {
    /// Returned when failed to read a Grid.
    #[error(transparent)]
    ReadFailure(bincode::Error),
    /// Returned while performing IO operations.
    #[error(transparent)]
    IoFailure(io::Error),
}

#[derive(Deserialize)]
struct Mmv1;

#[derive(Deserialize)]
struct Mmv2 {
    remapper: Option<BinRemapper>,
    key_value_db: HashMap<String, String>,
}

#[derive(Deserialize)]
struct Mmv3 {
    remapper: Option<BinRemapper>,
    key_value_db: HashMap<String, String>,
    _subgrid_template: SubgridEnum,
}

// ALLOW: fixing the warning will break the file format
#[allow(clippy::large_enum_variant)]
#[derive(Deserialize)]
enum MoreMembers {
    V1(Mmv1),
    V2(Mmv2),
    V3(Mmv3),
}

/// Main data structure of `PineAPPL`. This structure contains a `Subgrid` for each `LumiEntry`,
/// bin, and coupling order it was created with.
#[derive(Deserialize)]
pub struct Grid {
    subgrids: Array3<SubgridEnum>,
    channels: Vec<Channel>,
    bin_limits: BinLimits,
    orders: Vec<Order>,
    _subgrid_params: SubgridParams,
    more_members: MoreMembers,
}

impl Grid {
    /// Return by which convention the particle IDs are encoded.
    #[must_use]
    pub fn pid_basis(&self) -> PidBasis {
        if let Some(key_values) = self.key_values() {
            if let Some(lumi_id_types) = key_values.get("lumi_id_types") {
                match lumi_id_types.as_str() {
                    "pdg_mc_ids" => return PidBasis::Pdg,
                    "evol" => return PidBasis::Evol,
                    _ => unimplemented!("unknown particle ID convention {lumi_id_types}"),
                }
            }
        }

        // if there's no basis explicitly set we're assuming to use PDG IDs
        PidBasis::Pdg
    }

    /// Construct a `Grid` by deserializing it from `reader`.
    ///
    /// # Errors
    ///
    /// If reading from the compressed or uncompressed stream fails an error is returned.
    ///
    /// # Panics
    ///
    /// Panics if the grid version is not `0`.
    pub fn read_uncompressed(mut reader: impl BufRead) -> Result<Self, GridError> {
        let magic_bytes: [u8; 16] = reader.fill_buf().map_err(GridError::IoFailure)?[0..16]
            .try_into()
            .unwrap_or_else(|_| unreachable!());

        let file_version = if &magic_bytes[0..8] == b"PineAPPL" {
            reader.consume(16);
            u64::from_le_bytes(
                magic_bytes[8..16]
                    .try_into()
                    .unwrap_or_else(|_| unreachable!()),
            )
        } else {
            0
        };

        // should be guarateed not to happen, because `pineappl::grid::Grid::read` only calls this
        // method if the file version matches
        assert_eq!(file_version, 0);

        bincode::deserialize_from(reader).map_err(GridError::ReadFailure)
    }

    /// Return the channels for this `Grid`.
    #[must_use]
    pub fn channels(&self) -> &[Channel] {
        &self.channels
    }

    /// Return a vector containing the type of convolutions performed with this grid.
    ///
    /// # Panics
    ///
    /// Panics if the metadata key--value pairs `convolution_particle_1` and `convolution_type_1`,
    /// or `convolution_particle_2` and `convolution_type_2` are not correctly set.
    #[must_use]
    pub fn convolutions(&self) -> Vec<Convolution> {
        self.key_values().map_or_else(
            // if there isn't any metadata, we assume two unpolarized proton-PDFs are used
            || vec![Convolution::UnpolPDF(2212), Convolution::UnpolPDF(2212)],
            |kv| {
                // the current file format only supports exactly two convolutions
                (1..=2)
                    .map(|index| {
                        // if there are key-value pairs `convolution_particle_1` and
                        // `convolution_type_1` and the same with a higher index, we convert this
                        // metadata into `Convolution`
                        match (
                            kv.get(&format!("convolution_particle_{index}"))
                                .map(|s| s.parse::<i32>()),
                            kv.get(&format!("convolution_type_{index}"))
                                .map(String::as_str),
                        ) {
                            (_, Some("None")) => Convolution::None,
                            (Some(Ok(pid)), Some("UnpolPDF")) => Convolution::UnpolPDF(pid),
                            (Some(Ok(pid)), Some("PolPDF")) => Convolution::PolPDF(pid),
                            (Some(Ok(pid)), Some("UnpolFF")) => Convolution::UnpolFF(pid),
                            (Some(Ok(pid)), Some("PolFF")) => Convolution::PolFF(pid),
                            (None, None) => {
                                // if these key-value pairs are missing use the old metadata
                                match kv
                                    .get(&format!("initial_state_{index}"))
                                    .map(|s| s.parse::<i32>())
                                {
                                    Some(Ok(pid)) => {
                                        // Addresses: https://github.com/NNPDF/pineappl/issues/334
                                        if self.channels().is_empty() && pid == 2212 {
                                            Convolution::UnpolPDF(2212)
                                        } else {
                                            let condition = !self.channels().iter().all(|entry| {
                                                entry.entry().iter().all(|&(a, b, _)|
                                                    match index {
                                                        1 => a,
                                                        2 => b,
                                                        _ => unreachable!()
                                                    } == pid
                                                )
                                            });

                                            if condition {
                                                Convolution::UnpolPDF(pid)
                                            } else {
                                                Convolution::None
                                            }
                                        }
                                    }
                                    None => Convolution::UnpolPDF(2212),
                                    Some(Err(err)) => panic!("metadata 'initial_state_{index}' could not be parsed: {err}"),
                                }
                            }
                            (None, Some(_)) => {
                                panic!("metadata 'convolution_type_{index}' is missing")
                            }
                            (Some(_), None) => {
                                panic!("metadata 'convolution_particle_{index}' is missing")
                            }
                            (Some(Ok(_)), Some(type_)) => {
                                panic!("metadata 'convolution_type_{index} = {type_}' is unknown")
                            }
                            (Some(Err(err)), Some(_)) => panic!(
                                "metadata 'convolution_particle_{index}' could not be parsed: {err}"
                            ),
                        }
                    })
                    .collect()
            },
        )
    }

    /// Returns the subgrid parameters.
    #[must_use]
    pub fn orders(&self) -> &[Order] {
        &self.orders
    }

    /// Return all subgrids as an `ArrayView3`.
    #[must_use]
    pub fn subgrids(&self) -> ArrayView3<'_, SubgridEnum> {
        self.subgrids.view()
    }

    /// Return the currently set remapper, if there is any.
    #[must_use]
    pub const fn remapper(&self) -> Option<&BinRemapper> {
        match &self.more_members {
            MoreMembers::V1(_) => None,
            MoreMembers::V2(mmv2) => mmv2.remapper.as_ref(),
            MoreMembers::V3(mmv3) => mmv3.remapper.as_ref(),
        }
    }

    /// Returns all information about the bins in this grid.
    #[must_use]
    pub const fn bin_info(&self) -> BinInfo<'_> {
        BinInfo::new(&self.bin_limits, self.remapper())
    }

    /// Returns a map with key-value pairs, if there are any stored in this grid.
    #[must_use]
    pub const fn key_values(&self) -> Option<&HashMap<String, String>> {
        match &self.more_members {
            MoreMembers::V3(mmv3) => Some(&mmv3.key_value_db),
            MoreMembers::V2(mmv2) => Some(&mmv2.key_value_db),
            MoreMembers::V1(_) => None,
        }
    }
}
