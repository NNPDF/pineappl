//! Module containing all traits and supporting structures for grids.

use super::bin::{BinInfo, BinLimits, BinRemapper};
use super::boc::{Channel, Order};
use super::convolutions::Convolution;
use super::empty_subgrid::EmptySubgridV1;
use super::pids::{self, PidBasis};
use super::subgrid::{Subgrid, SubgridEnum, SubgridParams};
use bitflags::bitflags;
use git_version::git_version;
use lz4_flex::frame::FrameDecoder;
use ndarray::{Array3, ArrayView3};
use serde::{Deserialize, Serialize, Serializer};
use std::borrow::Cow;
use std::collections::{BTreeMap, HashMap};
use std::io::{self, BufRead, BufReader, Read};
use std::mem;
use std::ops::Range;
use thiserror::Error;

/// This structure represents a position (`x1`, `x2`, `q2`) in a `Subgrid` together with a
/// corresponding `weight`. The type `W` can either be a `f64` or `()`, which is used when multiple
/// weights should be signaled.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
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
    /// Returned when trying to merge two `Grid` objects with incompatible bin limits.
    #[error(transparent)]
    InvalidBinLimits(super::bin::MergeBinError),
    /// Returned if the number of bins in the grid and in the remapper do not agree.
    #[error("the remapper has {remapper_bins} bins, but the grid has {grid_bins}")]
    BinNumberMismatch {
        /// Number of bins in the grid.
        grid_bins: usize,
        /// Number of bins in the remapper.
        remapper_bins: usize,
    },
    /// Returned when it was tried to merge bins that are non-consecutive.
    #[error(transparent)]
    MergeBinError(super::bin::MergeBinError),
    /// Returned when trying to construct a `Grid` using an unknown subgrid type.
    #[error("tried constructing a Grid with unknown Subgrid type `{0}`")]
    UnknownSubgridType(String),
    /// Returned when failed to read a Grid.
    #[error(transparent)]
    ReadFailure(bincode::Error),
    /// Returned when failed to write a Grid.
    #[error(transparent)]
    WriteFailure(bincode::Error),
    /// Returned while performing IO operations.
    #[error(transparent)]
    IoFailure(io::Error),
    /// Returned when trying to read a `PineAPPL` file with file format version that is not
    /// supported.
    #[error("the file version is {file_version}, but supported is only {supported_version}")]
    FileVersionMismatch {
        /// File format version of the file read.
        file_version: u64,
        /// Maximum supported file format version for this library.
        supported_version: u64,
    },
    /// Returned from [`Grid::evolve`] if the evolution failed.
    #[error("failed to evolve grid: {0}")]
    EvolutionFailure(String),
    /// Errors that do no originate from this crate itself.
    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

#[derive(Clone, Deserialize, Serialize)]
struct Mmv1;

#[derive(Clone, Deserialize, Serialize)]
struct Mmv2 {
    remapper: Option<BinRemapper>,
    key_value_db: HashMap<String, String>,
}

fn ordered_map_serialize<S, K: Ord + Serialize, V: Serialize>(
    value: &HashMap<K, V>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let ordered: BTreeMap<_, _> = value.iter().collect();
    ordered.serialize(serializer)
}

#[derive(Clone, Deserialize, Serialize)]
struct Mmv3 {
    remapper: Option<BinRemapper>,
    // order the HashMap before serializing it to make the output stable
    #[serde(serialize_with = "ordered_map_serialize")]
    key_value_db: HashMap<String, String>,
    subgrid_template: SubgridEnum,
}

impl Default for Mmv2 {
    fn default() -> Self {
        Self {
            remapper: None,
            key_value_db: [
                (
                    "pineappl_gitversion".to_owned(),
                    git_version!(
                        args = ["--always", "--dirty", "--long", "--tags"],
                        cargo_prefix = "cargo:",
                        fallback = "unknown"
                    )
                    .to_owned(),
                ),
                // by default we assume there are protons in the initial state
                ("initial_state_1".to_owned(), "2212".to_owned()),
                ("initial_state_2".to_owned(), "2212".to_owned()),
            ]
            .iter()
            .cloned()
            .collect(),
        }
    }
}

// ALLOW: fixing the warning will break the file format
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Deserialize, Serialize)]
enum MoreMembers {
    V1(Mmv1),
    V2(Mmv2),
    V3(Mmv3),
}

impl MoreMembers {
    fn upgrade(&mut self) {
        match self {
            Self::V1(_) => {
                *self = Self::V2(Mmv2::default());
            }
            Self::V2(_) | Self::V3(_) => {}
        }
    }
}

bitflags! {
    /// Bitflags for optimizing a [`Grid`]. See [`Grid::optimize_using`].
    #[derive(Clone, Copy)]
    #[repr(transparent)]
    pub struct GridOptFlags: u32 {
        /// Change the [`Subgrid`] type to optimize storage effeciency.
        const OPTIMIZE_SUBGRID_TYPE = 0b1;
        /// Recognize whether a subgrid was filled with events with a static scale and if this is
        /// the case, optimize it by undoing the interpolation in the scale. This flag requires
        /// [`Self::OPTIMIZE_SUBGRID_TYPE`] to be active.
        const STATIC_SCALE_DETECTION = 0b10;
        /// If two channels differ by transposition of the two initial states and the functions
        /// this grid is convolved with are the same for both initial states, this will merge one
        /// channel into the other, with the correct transpositions.
        const SYMMETRIZE_CHANNELS = 0b100;
        /// Remove all orders ([`Grid::orders`]), which do not contain any non-zero subgrids.
        const STRIP_EMPTY_ORDERS = 0b1000;
        /// Merge the subgrids of channels which have the same definition.
        const MERGE_SAME_CHANNELS = 0b10000;
        /// Remove all channels ([`Grid::channels`]), which do not contain any non-zero subgrids.
        const STRIP_EMPTY_CHANNELS = 0b10_0000;
    }
}

/// Main data structure of `PineAPPL`. This structure contains a `Subgrid` for each `LumiEntry`,
/// bin, and coupling order it was created with.
#[derive(Clone, Deserialize, Serialize)]
pub struct Grid {
    subgrids: Array3<SubgridEnum>,
    channels: Vec<Channel>,
    bin_limits: BinLimits,
    orders: Vec<Order>,
    subgrid_params: SubgridParams,
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

    fn pdg_channels(&self) -> Cow<'_, [Channel]> {
        match self.pid_basis() {
            PidBasis::Evol => self
                .channels
                .iter()
                .map(|entry| Channel::translate(entry, &pids::evol_to_pdg_mc_ids))
                .collect(),
            PidBasis::Pdg => Cow::Borrowed(self.channels()),
        }
    }

    /// Construct a `Grid` by deserializing it from `reader`. Reading is buffered.
    ///
    /// # Errors
    ///
    /// If reading from the compressed or uncompressed stream fails an error is returned.
    pub fn read(reader: impl Read) -> Result<Self, GridError> {
        let mut reader = BufReader::new(reader);
        let buffer = reader.fill_buf().map_err(GridError::IoFailure)?;
        let magic_bytes: [u8; 4] = buffer[0..4].try_into().unwrap_or_else(|_| unreachable!());

        if u32::from_le_bytes(magic_bytes) == 0x18_4D_22_04 {
            Self::read_uncompressed(FrameDecoder::new(reader))
        } else {
            Self::read_uncompressed(reader)
        }
    }

    fn read_uncompressed(mut reader: impl BufRead) -> Result<Self, GridError> {
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

        if file_version != 0 {
            return Err(GridError::FileVersionMismatch {
                file_version,
                supported_version: 0,
            });
        }

        bincode::deserialize_from(reader).map_err(GridError::ReadFailure)
    }

    /// Return the channels for this `Grid`.
    #[must_use]
    pub fn channels(&self) -> &[Channel] {
        &self.channels
    }

    /// Merges the bins for the corresponding range together in a single one.
    ///
    /// # Errors
    ///
    /// When the given bins are non-consecutive, an error is returned.
    pub fn merge_bins(&mut self, bins: Range<usize>) -> Result<(), GridError> {
        self.bin_limits
            .merge_bins(bins.clone())
            .map_err(GridError::MergeBinError)?;

        if let Some(remapper) = self.remapper_mut() {
            remapper
                .merge_bins(bins.clone())
                .map_err(GridError::MergeBinError)?;
        }

        let bin_count = self.bin_info().bins();
        let mut old_subgrids = mem::replace(
            &mut self.subgrids,
            Array3::from_shape_simple_fn(
                (self.orders.len(), bin_count, self.channels.len()),
                || EmptySubgridV1.into(),
            ),
        );

        for ((order, bin, channel), subgrid) in old_subgrids.indexed_iter_mut() {
            if subgrid.is_empty() {
                continue;
            }

            if bins.contains(&bin) {
                let new_subgrid = &mut self.subgrids[[order, bins.start, channel]];

                if new_subgrid.is_empty() {
                    mem::swap(new_subgrid, subgrid);
                } else {
                    new_subgrid.merge(subgrid, false);
                }
            } else {
                let new_bin = if bin > bins.start {
                    bin - (bins.end - bins.start) + 1
                } else {
                    bin
                };

                mem::swap(&mut self.subgrids[[order, new_bin, channel]], subgrid);
            }
        }

        Ok(())
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
                                        let condition = !self.channels().iter().all(|entry| {
                                            entry.entry().iter().all(|&channels| match index {
                                                1 => channels.0 == pid,
                                                2 => channels.1 == pid,
                                                _ => unreachable!(),
                                            })
                                        });

                                        if condition {
                                            Convolution::UnpolPDF(pid)
                                        } else {
                                            Convolution::None
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

    fn remapper_mut(&mut self) -> Option<&mut BinRemapper> {
        match &mut self.more_members {
            MoreMembers::V1(_) => None,
            MoreMembers::V2(mmv2) => mmv2.remapper.as_mut(),
            MoreMembers::V3(mmv3) => mmv3.remapper.as_mut(),
        }
    }

    /// Returns all information about the bins in this grid.
    #[must_use]
    pub const fn bin_info(&self) -> BinInfo<'_> {
        BinInfo::new(&self.bin_limits, self.remapper())
    }

    /// Upgrades the internal data structures to their latest versions.
    pub fn upgrade(&mut self) {
        self.more_members.upgrade();
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
