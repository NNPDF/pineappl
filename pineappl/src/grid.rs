//! Module containing all traits and supporting structures for grids.

use super::bin::{BinInfo, BinLimits, BinRemapper};
use super::boc::{Channel, Order};
use super::convolutions::{Convolution, LumiCache};
use super::empty_subgrid::EmptySubgridV1;
use super::evolution::{self, AlphasTable, EvolveInfo, OperatorInfo, OperatorSliceInfo};
use super::fk_table::FkTable;
use super::import_only_subgrid::ImportOnlySubgridV2;
use super::lagrange_subgrid::{LagrangeSparseSubgridV1, LagrangeSubgridV1, LagrangeSubgridV2};
use super::ntuple_subgrid::NtupleSubgridV1;
use super::pids::{self, PidBasis};
use super::subgrid::{ExtraSubgridParams, Mu2, Subgrid, SubgridEnum, SubgridParams};
use bitflags::bitflags;
use float_cmp::{approx_eq, assert_approx_eq};
use git_version::git_version;
use lz4_flex::frame::{FrameDecoder, FrameEncoder};
use ndarray::{s, Array3, ArrayView3, ArrayView5, ArrayViewMut3, Axis, CowArray, Dimension, Ix4};
use serde::{Deserialize, Serialize, Serializer};
use std::borrow::Cow;
use std::collections::{BTreeMap, HashMap};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::iter;
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

impl Mmv3 {
    fn new(subgrid_template: SubgridEnum) -> Self {
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
                // by default we assume there are unpolarized protons in the initial state
                // do not change these to the new metadata to not break backwards compatibility
                ("initial_state_1".to_owned(), "2212".to_owned()),
                ("initial_state_2".to_owned(), "2212".to_owned()),
            ]
            .iter()
            .cloned()
            .collect(),
            subgrid_template,
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
    /// Constructor.
    #[must_use]
    pub fn new(
        channels: Vec<Channel>,
        orders: Vec<Order>,
        bin_limits: Vec<f64>,
        subgrid_params: SubgridParams,
    ) -> Self {
        Self {
            subgrids: Array3::from_shape_simple_fn(
                (orders.len(), bin_limits.len() - 1, channels.len()),
                || EmptySubgridV1.into(),
            ),
            orders,
            channels,
            bin_limits: BinLimits::new(bin_limits),
            more_members: MoreMembers::V3(Mmv3::new(
                LagrangeSubgridV2::new(&subgrid_params, &ExtraSubgridParams::from(&subgrid_params))
                    .into(),
            )),
            subgrid_params,
        }
    }

    /// Constructor. This function can be used like `new`, but the additional parameter
    /// `subgrid_type` selects the underlying `Subgrid` type. Supported values are:
    /// - `LagrangeSubgrid`
    /// - `LagrangeSparseSubgrid`
    /// - `NtupleSubgrid`
    ///
    /// # Errors
    ///
    /// If `subgrid_type` is none of the values listed above, an error is returned.
    pub fn with_subgrid_type(
        channels: Vec<Channel>,
        orders: Vec<Order>,
        bin_limits: Vec<f64>,
        subgrid_params: SubgridParams,
        extra: ExtraSubgridParams,
        subgrid_type: &str,
    ) -> Result<Self, GridError> {
        let subgrid_template: SubgridEnum = match subgrid_type {
            "LagrangeSubgrid" | "LagrangeSubgridV2" => {
                LagrangeSubgridV2::new(&subgrid_params, &extra).into()
            }
            "LagrangeSubgridV1" => LagrangeSubgridV1::new(&subgrid_params).into(),
            "NtupleSubgrid" => NtupleSubgridV1::new().into(),
            "LagrangeSparseSubgrid" => LagrangeSparseSubgridV1::new(&subgrid_params).into(),
            _ => return Err(GridError::UnknownSubgridType(subgrid_type.to_owned())),
        };

        Ok(Self {
            subgrids: Array3::from_shape_simple_fn(
                (orders.len(), bin_limits.len() - 1, channels.len()),
                || EmptySubgridV1.into(),
            ),
            orders,
            channels,
            bin_limits: BinLimits::new(bin_limits),
            subgrid_params,
            more_members: MoreMembers::V3(Mmv3::new(subgrid_template)),
        })
    }

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

    /// Set the convention by which PIDs of channels are interpreted.
    pub fn set_pid_basis(&mut self, pid_basis: PidBasis) {
        match pid_basis {
            PidBasis::Pdg => self.set_key_value("lumi_id_types", "pdg_mc_ids"),
            PidBasis::Evol => self.set_key_value("lumi_id_types", "evol"),
        }
    }

    fn pdg_channels(&self) -> Cow<[Channel]> {
        match self.pid_basis() {
            PidBasis::Evol => self
                .channels
                .iter()
                .map(|entry| Channel::translate(entry, &pids::evol_to_pdg_mc_ids))
                .collect(),
            PidBasis::Pdg => Cow::Borrowed(self.channels()),
        }
    }

    /// Perform a convolution using the PDFs and strong coupling in `lumi_cache`, and only
    /// selecting only the orders, bins and channels corresponding to `order_mask`, `bin_indices`
    /// and `channel_mask`. A variation of the scales is performed using the factors in `xi`; the
    /// first factor varies the renormalization scale, the second the factorization scale. Note
    /// that for the variation to be trusted all non-zero log-grids must be contained.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn convolve(
        &self,
        lumi_cache: &mut LumiCache,
        order_mask: &[bool],
        bin_indices: &[usize],
        channel_mask: &[bool],
        xi: &[(f64, f64)],
    ) -> Vec<f64> {
        lumi_cache.setup(self, xi).unwrap();

        let bin_indices = if bin_indices.is_empty() {
            (0..self.bin_info().bins()).collect()
        } else {
            bin_indices.to_vec()
        };
        let mut bins = vec![0.0; bin_indices.len() * xi.len()];
        let normalizations = self.bin_info().normalizations();
        let pdg_channels = self.pdg_channels();

        for (xi_index, &(xir, xif)) in xi.iter().enumerate() {
            for ((ord, bin, chan), subgrid) in self.subgrids.indexed_iter() {
                let order = &self.orders[ord];

                if ((order.logxir > 0) && (xir == 1.0)) || ((order.logxif > 0) && (xif == 1.0)) {
                    continue;
                }

                if (!order_mask.is_empty() && !order_mask[ord])
                    || (!channel_mask.is_empty() && !channel_mask[chan])
                {
                    continue;
                }

                let Some(bin_index) = bin_indices.iter().position(|&index| index == bin) else {
                    continue;
                };

                if subgrid.is_empty() {
                    continue;
                }

                let channel = &pdg_channels[chan];
                let mu2_grid = subgrid.mu2_grid();
                let x1_grid = subgrid.x1_grid();
                let x2_grid = subgrid.x2_grid();

                lumi_cache.set_grids(&mu2_grid, &x1_grid, &x2_grid, xir, xif);

                let mut value =
                    subgrid.convolve(&x1_grid, &x2_grid, &mu2_grid, &mut |ix1, ix2, imu2| {
                        let x1 = x1_grid[ix1];
                        let x2 = x2_grid[ix2];
                        let mut lumi = 0.0;

                        for entry in channel.entry() {
                            let xfx1 = lumi_cache.xfx1(entry.0, ix1, imu2);
                            let xfx2 = lumi_cache.xfx2(entry.1, ix2, imu2);
                            lumi += xfx1 * xfx2 * entry.2 / (x1 * x2);
                        }

                        let alphas = lumi_cache.alphas(imu2);

                        lumi *= alphas.powi(order.alphas.try_into().unwrap());
                        lumi
                    });

                if order.logxir > 0 {
                    value *= (xir * xir).ln().powi(order.logxir.try_into().unwrap());
                }

                if order.logxif > 0 {
                    value *= (xif * xif).ln().powi(order.logxif.try_into().unwrap());
                }

                bins[xi_index + xi.len() * bin_index] += value / normalizations[bin];
            }
        }

        bins
    }

    /// Convolutes a single subgrid `(order, bin, channel)` with the PDFs strong coupling given by
    /// `xfx1`, `xfx2` and `alphas`. The convolution result is fully differentially, such that the
    /// axes of the result correspond to the values given by the subgrid `q2`, `x1` and `x2` grid
    /// values.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn convolve_subgrid(
        &self,
        lumi_cache: &mut LumiCache,
        ord: usize,
        bin: usize,
        channel: usize,
        xir: f64,
        xif: f64,
    ) -> Array3<f64> {
        lumi_cache.setup(self, &[(xir, xif)]).unwrap();

        let normalizations = self.bin_info().normalizations();
        let pdg_channels = self.pdg_channels();

        let subgrid = &self.subgrids[[ord, bin, channel]];
        let order = &self.orders[ord];

        let channel = &pdg_channels[channel];
        let mu2_grid = subgrid.mu2_grid();
        let x1_grid = subgrid.x1_grid();
        let x2_grid = subgrid.x2_grid();

        lumi_cache.set_grids(&mu2_grid, &x1_grid, &x2_grid, xir, xif);

        let mut array = Array3::zeros((mu2_grid.len(), x1_grid.len(), x2_grid.len()));

        for ((imu2, ix1, ix2), value) in subgrid.indexed_iter() {
            let x1 = x1_grid[ix1];
            let x2 = x2_grid[ix2];
            let mut lumi = 0.0;

            for entry in channel.entry() {
                let xfx1 = lumi_cache.xfx1(entry.0, ix1, imu2);
                let xfx2 = lumi_cache.xfx2(entry.1, ix2, imu2);
                lumi += xfx1 * xfx2 * entry.2 / (x1 * x2);
            }

            let alphas = lumi_cache.alphas(imu2);

            lumi *= alphas.powi(order.alphas.try_into().unwrap());

            array[[imu2, ix1, ix2]] = lumi * value;
        }

        if order.logxir > 0 {
            array *= (xir * xir).ln().powi(order.logxir.try_into().unwrap());
        }

        if order.logxif > 0 {
            array *= (xif * xif).ln().powi(order.logxif.try_into().unwrap());
        }

        array /= normalizations[bin];
        array
    }

    /// Fills the grid with an ntuple for the given `order`, `observable`, and `channel`.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn fill(&mut self, order: usize, observable: f64, channel: usize, ntuple: &Ntuple<f64>) {
        if let Some(bin) = self.bin_limits.index(observable) {
            let subgrid = &mut self.subgrids[[order, bin, channel]];
            if let SubgridEnum::EmptySubgridV1(_) = subgrid {
                if let MoreMembers::V3(mmv3) = &self.more_members {
                    *subgrid = mmv3.subgrid_template.clone_empty();
                } else {
                    unreachable!();
                }
            }

            subgrid.fill(ntuple);
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

    /// Serializes `self` into `writer`. Writing is buffered.
    ///
    /// # Errors
    ///
    /// If writing fails an error is returned.
    pub fn write(&self, writer: impl Write) -> Result<(), GridError> {
        let mut writer = BufWriter::new(writer);
        let file_header = b"PineAPPL\0\0\0\0\0\0\0\0";

        // first write PineAPPL file header
        writer.write(file_header).map_err(GridError::IoFailure)?;

        // then serialize
        bincode::serialize_into(writer, self).map_err(GridError::WriteFailure)
    }

    /// Serializes `self` into `writer`, using LZ4 compression. Writing is buffered.
    ///
    /// # Errors
    ///
    /// If writing or compression fails an error is returned.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn write_lz4(&self, writer: impl Write) -> Result<(), GridError> {
        let mut encoder = FrameEncoder::new(writer);
        self.write(&mut encoder)?;
        // TODO: get rid of the unwrap call and return the error
        encoder.try_finish().unwrap();

        Ok(())
    }

    /// Fills the grid with events for the parton momentum fractions `x1` and `x2`, the scale `q2`,
    /// and the `order` and `observable`. The events are stored in `weights` and their ordering
    /// corresponds to the ordering of [`Grid::channels`].
    pub fn fill_all(
        &mut self,
        order: usize,
        observable: f64,
        ntuple: &Ntuple<()>,
        weights: &[f64],
    ) {
        for (channel, weight) in weights.iter().enumerate() {
            self.fill(
                order,
                observable,
                channel,
                &Ntuple {
                    x1: ntuple.x1,
                    x2: ntuple.x2,
                    q2: ntuple.q2,
                    weight: *weight,
                },
            );
        }
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

    /// Merges the non-empty `Subgrid`s contained in `other` into `self`.
    ///
    /// # Errors
    ///
    /// If the bin limits of `self` and `other` are different and if the bin limits of `other` can
    /// not be merged with `self` an error is returned.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn merge(&mut self, mut other: Self) -> Result<(), GridError> {
        let mut new_orders: Vec<Order> = Vec::new();
        let mut new_bins = 0;
        let mut new_entries: Vec<Channel> = Vec::new();

        if self.bin_info() != other.bin_info() {
            let lhs_bins = self.bin_info().bins();
            new_bins = other.bin_info().bins();

            let lhs_remapper = self.remapper_mut();
            let rhs_remapper = other.remapper();

            if let Some(lhs) = lhs_remapper {
                if let Some(rhs) = rhs_remapper {
                    lhs.merge(rhs).map_err(GridError::MergeBinError)?;

                    let a = u32::try_from(lhs_bins).unwrap_or_else(|_| unreachable!());
                    let b = u32::try_from(lhs_bins + new_bins).unwrap_or_else(|_| unreachable!());

                    self.bin_limits = BinLimits::new((0..=b).map(f64::from).collect());
                    other.bin_limits = BinLimits::new((a..=b).map(f64::from).collect());
                } else {
                    // Return an error
                    todo!();
                }
            } else if rhs_remapper.is_none() {
                self.bin_limits
                    .merge(&other.bin_limits)
                    .map_err(GridError::InvalidBinLimits)?;
            } else {
                // Return an error
                todo!();
            }
        }

        for ((i, _, k), _) in other
            .subgrids
            .indexed_iter_mut()
            .filter(|((_, _, _), subgrid)| !subgrid.is_empty())
        {
            let other_order = &other.orders[i];
            let other_entry = &other.channels[k];

            if !self
                .orders
                .iter()
                .chain(new_orders.iter())
                .any(|x| x == other_order)
            {
                new_orders.push(other_order.clone());
            }

            if !self
                .channels()
                .iter()
                .chain(new_entries.iter())
                .any(|y| y == other_entry)
            {
                new_entries.push(other_entry.clone());
            }
        }

        if !new_orders.is_empty() || !new_entries.is_empty() || (new_bins != 0) {
            self.increase_shape(&(new_orders.len(), new_bins, new_entries.len()));
        }

        self.orders.append(&mut new_orders);
        self.channels.append(&mut new_entries);

        let bin_indices: Vec<_> = (0..other.bin_info().bins())
            .map(|bin| {
                self.bin_info()
                    .find_bin(&other.bin_info().bin_limits(bin))
                    .unwrap_or_else(|| panic!("failed for {bin}"))
            })
            .collect();

        for ((i, j, k), subgrid) in other
            .subgrids
            .indexed_iter_mut()
            .filter(|((_, _, _), subgrid)| !subgrid.is_empty())
        {
            let other_order = &other.orders[i];
            let other_entry = &other.channels[k];

            let self_i = self.orders.iter().position(|x| x == other_order).unwrap();
            let self_j = bin_indices[j];
            let self_k = self.channels.iter().position(|y| y == other_entry).unwrap();

            if self.subgrids[[self_i, self_j, self_k]].is_empty() {
                mem::swap(&mut self.subgrids[[self_i, self_j, self_k]], subgrid);
            } else {
                self.subgrids[[self_i, self_j, self_k]].merge(&mut *subgrid, false);
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

    /// Set the convolution type for this grid for the corresponding `index`.
    pub fn set_convolution(&mut self, index: usize, convolution: Convolution) {
        // remove outdated metadata
        self.key_values_mut()
            .remove(&format!("initial_state_{}", index + 1));

        let (type_, particle) = match convolution {
            Convolution::UnpolPDF(pid) => ("UnpolPDF".to_owned(), pid.to_string()),
            Convolution::PolPDF(pid) => ("PolPDF".to_owned(), pid.to_string()),
            Convolution::UnpolFF(pid) => ("UnpolFF".to_owned(), pid.to_string()),
            Convolution::PolFF(pid) => ("PolFF".to_owned(), pid.to_string()),
            Convolution::None => ("None".to_owned(), String::new()),
        };

        self.set_key_value(&format!("convolution_type_{}", index + 1), &type_);
        self.set_key_value(&format!("convolution_particle_{}", index + 1), &particle);

        // update the remaining metadata
        for (index, convolution) in self.convolutions().into_iter().enumerate() {
            if self
                .key_values()
                // UNWRAP: we set some key-values before so there must be a storage
                .unwrap_or_else(|| unreachable!())
                .get(&format!("initial_state_{}", index + 1))
                .is_some()
            {
                self.set_convolution(index, convolution);
            }
        }
    }

    fn increase_shape(&mut self, new_dim: &(usize, usize, usize)) {
        let old_dim = self.subgrids.raw_dim().into_pattern();
        let mut new_subgrids = Array3::from_shape_simple_fn(
            (
                old_dim.0 + new_dim.0,
                old_dim.1 + new_dim.1,
                old_dim.2 + new_dim.2,
            ),
            || EmptySubgridV1.into(),
        );

        for ((i, j, k), subgrid) in self.subgrids.indexed_iter_mut() {
            mem::swap(&mut new_subgrids[[i, j, k]], subgrid);
        }

        mem::swap(&mut self.subgrids, &mut new_subgrids);
    }

    /// Scale all subgrids by `factor`.
    pub fn scale(&mut self, factor: f64) {
        self.subgrids
            .iter_mut()
            .for_each(|subgrid| subgrid.scale(factor));
    }

    /// Scales each subgrid by a factor which is the product of the given values `alphas`, `alpha`,
    /// `logxir`, and `logxif`, each raised to the corresponding powers for each subgrid. In
    /// addition, every subgrid is scaled by a factor `global` independently of its order.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn scale_by_order(
        &mut self,
        alphas: f64,
        alpha: f64,
        logxir: f64,
        logxif: f64,
        global: f64,
    ) {
        for ((i, _, _), subgrid) in self.subgrids.indexed_iter_mut() {
            let order = &self.orders[i];
            let factor = global
                * alphas.powi(order.alphas.try_into().unwrap())
                * alpha.powi(order.alpha.try_into().unwrap())
                * logxir.powi(order.logxir.try_into().unwrap())
                * logxif.powi(order.logxif.try_into().unwrap());

            subgrid.scale(factor);
        }
    }

    /// Scales each subgrid by a bin-dependent factor given in `factors`. If a bin does not have a
    /// corresponding entry in `factors` it is not rescaled. If `factors` has more entries than
    /// there are bins the superfluous entries do not have an effect.
    pub fn scale_by_bin(&mut self, factors: &[f64]) {
        for ((_, bin, _), subgrid) in self.subgrids.indexed_iter_mut() {
            if let Some(&factor) = factors.get(bin) {
                subgrid.scale(factor);
            }
        }
    }

    /// Returns the subgrid parameters.
    #[must_use]
    pub fn orders(&self) -> &[Order] {
        &self.orders
    }

    /// Return a mutable reference to the subgrid parameters.
    #[must_use]
    pub fn orders_mut(&mut self) -> &mut [Order] {
        &mut self.orders
    }

    /// Return a mutable reference to the grid's channels.
    pub fn channels_mut(&mut self) -> &mut [Channel] {
        &mut self.channels
    }

    /// Return all subgrids as an `ArrayView3`.
    #[must_use]
    pub fn subgrids(&self) -> ArrayView3<SubgridEnum> {
        self.subgrids.view()
    }

    /// Return all subgrids as an `ArrayViewMut3`.
    #[must_use]
    pub fn subgrids_mut(&mut self) -> ArrayViewMut3<SubgridEnum> {
        self.subgrids.view_mut()
    }

    /// Sets a remapper. A remapper can change the dimensions and limits of each bin in this grid.
    /// This is useful because many Monte Carlo integrators and also `PineAPPL` do not support
    /// multi-dimensional bins. To work around the problem the multi-dimensional bins can be
    /// projected to one-dimensional bins, and the remapper can be used to restore the multi
    /// dimensionality. Furthermore, it allows to normalize each bin separately, and independently
    /// of the bin widths.
    ///
    /// # Errors
    ///
    /// Returns an error if the number of bins in the grid and in the remapper do not agree.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn set_remapper(&mut self, remapper: BinRemapper) -> Result<(), GridError> {
        if remapper.bins() != self.bin_info().bins() {
            return Err(GridError::BinNumberMismatch {
                grid_bins: self.bin_info().bins(),
                remapper_bins: remapper.bins(),
            });
        }

        self.more_members.upgrade();

        match &mut self.more_members {
            MoreMembers::V1(_) => unreachable!(),
            MoreMembers::V2(mmv2) => mmv2.remapper = Some(remapper),
            MoreMembers::V3(mmv3) => mmv3.remapper = Some(remapper),
        }

        Ok(())
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
    pub const fn bin_info(&self) -> BinInfo {
        BinInfo::new(&self.bin_limits, self.remapper())
    }

    /// Calls [`Self::optimize_using`] with all possible optimization options
    /// ([`GridOptFlags::all`]).
    pub fn optimize(&mut self) {
        self.optimize_using(GridOptFlags::all());
    }

    /// Optimizes the internal datastructures for space efficiency. The parameter `flags`
    /// determines which optimizations are applied, see [`GridOptFlags`].
    pub fn optimize_using(&mut self, flags: GridOptFlags) {
        if flags.contains(GridOptFlags::OPTIMIZE_SUBGRID_TYPE) {
            let ssd = flags.contains(GridOptFlags::STATIC_SCALE_DETECTION);
            self.optimize_subgrid_type(ssd);
        }
        if flags.contains(GridOptFlags::SYMMETRIZE_CHANNELS) {
            self.symmetrize_channels();
        }
        if flags.contains(GridOptFlags::STRIP_EMPTY_ORDERS) {
            self.strip_empty_orders();
        }
        if flags.contains(GridOptFlags::MERGE_SAME_CHANNELS) {
            self.merge_same_channels();
        }
        if flags.contains(GridOptFlags::STRIP_EMPTY_CHANNELS) {
            self.strip_empty_channels();
        }
    }

    fn optimize_subgrid_type(&mut self, static_scale_detection: bool) {
        for subgrid in &mut self.subgrids {
            match subgrid {
                // replace empty subgrids of any type with `EmptySubgridV1`
                _ if subgrid.is_empty() => {
                    *subgrid = EmptySubgridV1.into();
                }
                // can't be optimized without losing information
                SubgridEnum::NtupleSubgridV1(_) => continue,
                _ => {
                    // TODO: this requires a `pub(crate)` in `LagrangeSubgridV2`; we should
                    // replace this with a method
                    if !static_scale_detection {
                        if let SubgridEnum::LagrangeSubgridV2(subgrid) = subgrid {
                            // disable static-scale detection
                            subgrid.static_q2 = -1.0;
                        }
                    }

                    let mut new_subgrid = ImportOnlySubgridV2::from(&*subgrid).into();
                    mem::swap(subgrid, &mut new_subgrid);
                }
            }
        }
    }

    /// Try to deduplicate channels by detecting pairs of them that contain the same subgrids. The
    /// numerical equality is tested using a tolerance of `ulps`, given in [units of least
    /// precision](https://docs.rs/float-cmp/latest/float_cmp/index.html#some-explanation).
    pub fn dedup_channels(&mut self, ulps: i64) {
        let mut indices: Vec<usize> = (0..self.channels.len()).collect();

        while let Some(index) = indices.pop() {
            if let Some(other_index) = indices.iter().copied().find(|&other_index| {
                let (mut a, mut b) = self
                    .subgrids
                    .multi_slice_mut((s![.., .., other_index], s![.., .., index]));

                // TODO: use `Iterator::eq_by` once stablizied
                for (lhs, rhs) in a.iter_mut().zip(b.iter_mut()) {
                    let mut it_a = lhs.indexed_iter();
                    let mut it_b = rhs.indexed_iter();

                    loop {
                        let a = it_a.next();
                        let b = it_b.next();

                        match (a, b) {
                            (Some((tuple_a, value_a)), Some((tuple_b, value_b))) => {
                                if tuple_a != tuple_b {
                                    return false;
                                }

                                let u = ulps;
                                if !approx_eq!(f64, value_a, value_b, ulps = u) {
                                    return false;
                                }
                            }
                            (None, None) => break,
                            _ => return false,
                        }
                    }
                }

                true
            }) {
                let old_channel = self.channels.remove(index).entry().to_vec();
                let mut new_channel = self.channels[other_index].entry().to_vec();
                new_channel.extend(old_channel);
                self.channels[other_index] = Channel::new(new_channel);
                self.subgrids.remove_index(Axis(2), index);
            }
        }
    }

    fn merge_same_channels(&mut self) {
        let mut indices: Vec<_> = (0..self.channels.len()).rev().collect();

        // merge channels that are the same
        while let Some(index) = indices.pop() {
            if let Some((other_index, factor)) = indices.iter().find_map(|&i| {
                self.channels[i]
                    .common_factor(&self.channels[index])
                    .map(|factor| (i, factor))
            }) {
                let (mut a, mut b) = self
                    .subgrids
                    .multi_slice_mut((s![.., .., other_index], s![.., .., index]));

                // check if in all cases the limits are compatible with merging
                for (lhs, rhs) in a.iter_mut().zip(b.iter_mut()) {
                    if !rhs.is_empty() {
                        rhs.scale(1.0 / factor);
                        if lhs.is_empty() {
                            // we can't merge into an EmptySubgridV1
                            *lhs = rhs.clone_empty();
                        }
                        lhs.merge(rhs, false);

                        *rhs = EmptySubgridV1.into();
                    }
                }
            }
        }
    }

    fn strip_empty_channels(&mut self) {
        let mut keep_channel_indices = vec![];
        let mut new_channel_entries = vec![];

        // only keep channels that have non-zero factors and for which at least one subgrid is
        // non-empty
        for (channel, entry) in self.channels.iter().enumerate() {
            if !entry.entry().iter().all(|&(_, _, factor)| factor == 0.0)
                && !self
                    .subgrids
                    .slice(s![.., .., channel])
                    .iter()
                    .all(Subgrid::is_empty)
            {
                keep_channel_indices.push(channel);
                new_channel_entries.push(entry.clone());
            }
        }

        // only keep the previously selected subgrids
        let new_subgrids = Array3::from_shape_fn(
            (
                self.orders.len(),
                self.bin_info().bins(),
                keep_channel_indices.len(),
            ),
            |(order, bin, new_channel)| {
                mem::replace(
                    &mut self.subgrids[[order, bin, keep_channel_indices[new_channel]]],
                    EmptySubgridV1.into(),
                )
            },
        );

        self.channels = new_channel_entries;
        self.subgrids = new_subgrids;
    }

    fn strip_empty_orders(&mut self) {
        let mut indices: Vec<_> = (0..self.orders().len()).collect();

        while let Some(index) = indices.pop() {
            if self
                .subgrids
                .slice(s![index, .., ..])
                .iter()
                .all(Subgrid::is_empty)
            {
                self.orders.remove(index);
                self.subgrids.remove_index(Axis(0), index);
            }
        }
    }

    fn symmetrize_channels(&mut self) {
        let convolutions = self.convolutions();
        if convolutions[0] != convolutions[1] {
            return;
        }

        let mut indices: Vec<usize> = (0..self.channels.len()).rev().collect();

        while let Some(index) = indices.pop() {
            let channel_entry = &self.channels[index];

            if *channel_entry == channel_entry.transpose() {
                // check if in all cases the limits are compatible with merging
                self.subgrids
                    .slice_mut(s![.., .., index])
                    .iter_mut()
                    .for_each(|subgrid| {
                        if !subgrid.is_empty() && (subgrid.x1_grid() == subgrid.x2_grid()) {
                            subgrid.symmetrize();
                        }
                    });
            } else if let Some((j, &other_index)) = indices
                .iter()
                .enumerate()
                .find(|(_, i)| self.channels[**i] == channel_entry.transpose())
            {
                indices.remove(j);

                // check if in all cases the limits are compatible with merging
                let (mut a, mut b) = self
                    .subgrids
                    .multi_slice_mut((s![.., .., index], s![.., .., other_index]));

                for (lhs, rhs) in a.iter_mut().zip(b.iter_mut()) {
                    if !rhs.is_empty() {
                        if lhs.is_empty() {
                            // we can't merge into an EmptySubgridV1
                            *lhs = rhs.clone_empty();
                        }

                        lhs.merge(rhs, true);
                        *rhs = EmptySubgridV1.into();
                    }
                }
            }
        }
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

    /// Returns a map with key-value pairs, if there are any stored in this grid.
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
    pub fn key_values_mut(&mut self) -> &mut HashMap<String, String> {
        self.more_members.upgrade();

        match &mut self.more_members {
            MoreMembers::V1(_) => unreachable!(),
            MoreMembers::V2(mmv2) => &mut mmv2.key_value_db,
            MoreMembers::V3(mmv3) => &mut mmv3.key_value_db,
        }
    }

    /// Sets a specific key-value pair in this grid.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn set_key_value(&mut self, key: &str, value: &str) {
        self.key_values_mut()
            .insert(key.to_owned(), value.to_owned());
    }

    /// Returns information for the generation of evolution operators that are being used in
    /// [`Grid::evolve`] with the parameter `order_mask`.
    #[must_use]
    pub fn evolve_info(&self, order_mask: &[bool]) -> EvolveInfo {
        use super::evolution::EVOLVE_INFO_TOL_ULPS;

        let has_pdf1 = self.convolutions()[0] != Convolution::None;
        let has_pdf2 = self.convolutions()[1] != Convolution::None;

        let mut ren1 = Vec::new();
        let mut fac1 = Vec::new();
        let mut x1 = Vec::new();
        let mut pids1 = Vec::new();

        for (channel, subgrid) in self
            .subgrids()
            .indexed_iter()
            .filter_map(|(tuple, subgrid)| {
                (!subgrid.is_empty() && (order_mask.is_empty() || order_mask[tuple.0]))
                    .then_some((tuple.2, subgrid))
            })
        {
            ren1.extend(subgrid.mu2_grid().iter().map(|Mu2 { ren, .. }| *ren));
            ren1.sort_by(f64::total_cmp);
            ren1.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = EVOLVE_INFO_TOL_ULPS));

            fac1.extend(subgrid.mu2_grid().iter().map(|Mu2 { fac, .. }| *fac));
            fac1.sort_by(f64::total_cmp);
            fac1.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = EVOLVE_INFO_TOL_ULPS));

            if has_pdf1 {
                x1.extend(subgrid.x1_grid().iter().copied());
            }
            if has_pdf2 {
                x1.extend(subgrid.x2_grid().iter());
            }

            x1.sort_by(f64::total_cmp);
            x1.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = EVOLVE_INFO_TOL_ULPS));

            if has_pdf1 {
                pids1.extend(self.channels()[channel].entry().iter().map(|(a, _, _)| a));
            }
            if has_pdf2 {
                pids1.extend(self.channels()[channel].entry().iter().map(|(_, b, _)| b));
            }

            pids1.sort_unstable();
            pids1.dedup();
        }

        EvolveInfo {
            fac1,
            pids1,
            x1,
            ren1,
        }
    }

    /// Converts this `Grid` into an [`FkTable`] using an evolution kernel operator (EKO) given as
    /// `operator`. The dimensions and properties of this operator must be described using `info`.
    /// The parameter `order_mask` can be used to include or exclude orders from this operation,
    /// and must correspond to the ordering given by [`Grid::orders`]. Orders that are not given
    /// are enabled, and in particular if `order_mask` is empty all orders are activated.
    ///
    /// # Errors
    ///
    /// Returns a [`GridError::EvolutionFailure`] if either the `operator` or its `info` is
    /// incompatible with this `Grid`.
    #[deprecated(since = "0.7.4", note = "use evolve_with_slice_iter instead")]
    pub fn evolve(
        &self,
        operator: ArrayView5<f64>,
        info: &OperatorInfo,
        order_mask: &[bool],
    ) -> Result<FkTable, GridError> {
        self.evolve_with_slice_iter(
            info.fac1
                .iter()
                .zip(operator.axis_iter(Axis(0)))
                .map(|(&fac1, op)| {
                    Ok::<_, GridError>((
                        OperatorSliceInfo {
                            fac0: info.fac0,
                            pids0: info.pids0.clone(),
                            x0: info.x0.clone(),
                            fac1,
                            pids1: info.pids1.clone(),
                            x1: info.x1.clone(),
                            pid_basis: info.pid_basis,
                        },
                        CowArray::from(op),
                    ))
                }),
            order_mask,
            (info.xir, info.xif),
            &AlphasTable {
                ren1: info.ren1.clone(),
                alphas: info.alphas.clone(),
            },
        )
    }

    // TODO:
    // - try to find a better solution than to require that E must be convertible into
    //   anyhow::Error

    /// Converts this `Grid` into an [`FkTable`] using `slices` that must iterate over a [`Result`]
    /// of tuples of an [`OperatorSliceInfo`] and the corresponding sliced operator. The parameter
    /// `order_mask` can be used to include or exclude orders from this operation, and must
    /// correspond to the ordering given by [`Grid::orders`]. Orders that are not given are
    /// enabled, and in particular if `order_mask` is empty all orders are activated.
    ///
    /// # Errors
    ///
    /// Returns a [`GridError::EvolutionFailure`] if either the `operator` or its `info` is
    /// incompatible with this `Grid`. Returns a [`GridError::Other`] if the iterator from `slices`
    /// return an error.
    pub fn evolve_with_slice_iter<'a, E: Into<anyhow::Error>>(
        &self,
        slices: impl IntoIterator<Item = Result<(OperatorSliceInfo, CowArray<'a, f64, Ix4>), E>>,
        order_mask: &[bool],
        xi: (f64, f64),
        alphas_table: &AlphasTable,
    ) -> Result<FkTable, GridError> {
        use super::evolution::EVOLVE_INFO_TOL_ULPS;

        let mut lhs: Option<Self> = None;
        // Q2 slices we use
        let mut used_op_fac1 = Vec::new();
        // Q2 slices we encounter, but possibly don't use
        let mut op_fac1 = Vec::new();
        // Q2 slices needed by the grid
        let grid_fac1: Vec<_> = self
            .evolve_info(order_mask)
            .fac1
            .into_iter()
            .map(|fac| xi.1 * xi.1 * fac)
            .collect();

        for result in slices {
            let (info, operator) = result.map_err(|err| GridError::Other(err.into()))?;

            op_fac1.push(info.fac1);

            // it's possible that due to small numerical differences we get two slices which are
            // almost the same. We have to skip those in order not to evolve the 'same' slice twice
            if used_op_fac1
                .iter()
                .any(|&fac| approx_eq!(f64, fac, info.fac1, ulps = EVOLVE_INFO_TOL_ULPS))
            {
                continue;
            }

            // skip slices that the grid doesn't use
            if !grid_fac1
                .iter()
                .any(|&fac| approx_eq!(f64, fac, info.fac1, ulps = EVOLVE_INFO_TOL_ULPS))
            {
                continue;
            }

            let op_info_dim = (
                info.pids1.len(),
                info.x1.len(),
                info.pids0.len(),
                info.x0.len(),
            );

            if operator.dim() != op_info_dim {
                return Err(GridError::EvolutionFailure(format!(
                    "operator information {:?} does not match the operator's dimensions: {:?}",
                    op_info_dim,
                    operator.dim(),
                )));
            }

            let view = operator.view();

            let (subgrids, channels) = if self.convolutions()[0] != Convolution::None
                && self.convolutions()[1] != Convolution::None
            {
                evolution::evolve_slice_with_two(self, &view, &info, order_mask, xi, alphas_table)
            } else {
                evolution::evolve_slice_with_one(self, &view, &info, order_mask, xi, alphas_table)
            }?;

            let mut rhs = Self {
                subgrids,
                channels,
                bin_limits: self.bin_limits.clone(),
                orders: vec![Order::new(0, 0, 0, 0)],
                subgrid_params: SubgridParams::default(),
                more_members: self.more_members.clone(),
            };

            // TODO: use a new constructor to set this information
            rhs.set_pid_basis(info.pid_basis);

            if let Some(lhs) = &mut lhs {
                lhs.merge(rhs)?;
            } else {
                lhs = Some(rhs);
            }

            used_op_fac1.push(info.fac1);
        }

        // UNWRAP: if we can't compare two numbers there's a bug
        op_fac1.sort_by(|a, b| a.partial_cmp(b).unwrap_or_else(|| unreachable!()));

        // make sure we've evolved all slices
        if let Some(muf2) = grid_fac1.into_iter().find(|&grid_mu2| {
            !used_op_fac1
                .iter()
                .any(|&eko_mu2| approx_eq!(f64, grid_mu2, eko_mu2, ulps = EVOLVE_INFO_TOL_ULPS))
        }) {
            return Err(GridError::EvolutionFailure(format!(
                "no operator for muf2 = {muf2} found in {op_fac1:?}"
            )));
        }

        // TODO: convert this unwrap into error
        let grid = lhs.unwrap();

        // UNWRAP: merging evolved slices should be a proper FkTable again
        Ok(FkTable::try_from(grid).unwrap_or_else(|_| unreachable!()))
    }

    /// Converts this `Grid` into an [`FkTable`] using `slices` that must iterate over a [`Result`]
    /// of tuples of an [`OperatorSliceInfo`] and the corresponding sliced operator. The parameter
    /// `order_mask` can be used to include or exclude orders from this operation, and must
    /// correspond to the ordering given by [`Grid::orders`]. Orders that are not given are
    /// enabled, and in particular if `order_mask` is empty all orders are activated.
    ///
    /// # Errors
    ///
    /// Returns a [`GridError::EvolutionFailure`] if either the `operator` or its `info` is
    /// incompatible with this `Grid`. Returns a [`GridError::Other`] if the iterator from `slices`
    /// return an error.
    pub fn evolve_with_slice_iter2<'a, E: Into<anyhow::Error>>(
        &self,
        slices_a: impl IntoIterator<Item = Result<(OperatorSliceInfo, CowArray<'a, f64, Ix4>), E>>,
        slices_b: impl IntoIterator<Item = Result<(OperatorSliceInfo, CowArray<'a, f64, Ix4>), E>>,
        order_mask: &[bool],
        xi: (f64, f64),
        alphas_table: &AlphasTable,
    ) -> Result<FkTable, GridError> {
        use super::evolution::EVOLVE_INFO_TOL_ULPS;
        use itertools::izip;

        let mut lhs: Option<Self> = None;
        // Q2 slices we use
        let mut used_op_fac1 = Vec::new();
        // Q2 slices we encounter, but possibly don't use
        let mut op_fac1 = Vec::new();
        // Q2 slices needed by the grid
        let grid_fac1: Vec<_> = self
            .evolve_info(order_mask)
            .fac1
            .into_iter()
            .map(|fac| xi.1 * xi.1 * fac)
            .collect();

        // TODO: simplify the ugly repetition below by offloading some ops into fn
        for (result_a, result_b) in izip!(slices_a, slices_b) {
            // Operate on `slices_a`
            let (info_a, operator_a) = result_a.map_err(|err| GridError::Other(err.into()))?;
            // Operate on `slices_b`
            let (info_b, operator_b) = result_b.map_err(|err| GridError::Other(err.into()))?;

            // TODO: what if the scales of the EKOs don't agree? Is there an ordering problem?
            assert_approx_eq!(f64, info_a.fac1, info_b.fac1, ulps = EVOLVE_INFO_TOL_ULPS);

            // also the PID bases must be the same
            assert_eq!(info_a.pid_basis, info_b.pid_basis);

            op_fac1.push(info_a.fac1);

            // it's possible that due to small numerical differences we get two slices which are
            // almost the same. We have to skip those in order not to evolve the 'same' slice twice
            if used_op_fac1
                .iter()
                .any(|&fac| approx_eq!(f64, fac, info_a.fac1, ulps = EVOLVE_INFO_TOL_ULPS))
            {
                continue;
            }

            // skip slices that the grid doesn't use
            if !grid_fac1
                .iter()
                .any(|&fac| approx_eq!(f64, fac, info_a.fac1, ulps = EVOLVE_INFO_TOL_ULPS))
            {
                continue;
            }

            let op_info_dim_a = (
                info_a.pids1.len(),
                info_a.x1.len(),
                info_a.pids0.len(),
                info_a.x0.len(),
            );

            if operator_a.dim() != op_info_dim_a {
                return Err(GridError::EvolutionFailure(format!(
                    "operator information {:?} does not match the operator's dimensions: {:?}",
                    op_info_dim_a,
                    operator_a.dim(),
                )));
            }

            let op_info_dim_b = (
                info_b.pids1.len(),
                info_b.x1.len(),
                info_b.pids0.len(),
                info_b.x0.len(),
            );

            if operator_b.dim() != op_info_dim_b {
                return Err(GridError::EvolutionFailure(format!(
                    "operator information {:?} does not match the operator's dimensions: {:?}",
                    op_info_dim_b,
                    operator_b.dim(),
                )));
            }

            let views = [operator_a.view(), operator_b.view()];
            let infos = [info_a, info_b];

            let (subgrids, channels) = if self.convolutions()[0] != Convolution::None
                && self.convolutions()[1] != Convolution::None
            {
                evolution::evolve_slice_with_two2(
                    self,
                    &views,
                    &infos,
                    order_mask,
                    xi,
                    alphas_table,
                )
            } else {
                evolution::evolve_slice_with_one(
                    self,
                    &views[0],
                    &infos[1],
                    order_mask,
                    xi,
                    alphas_table,
                )
            }?;

            let mut rhs = Self {
                subgrids,
                channels,
                bin_limits: self.bin_limits.clone(),
                orders: vec![Order::new(0, 0, 0, 0)],
                subgrid_params: SubgridParams::default(),
                more_members: self.more_members.clone(),
            };

            // TODO: use a new constructor to set this information
            rhs.set_pid_basis(infos[0].pid_basis);

            if let Some(lhs) = &mut lhs {
                lhs.merge(rhs)?;
            } else {
                lhs = Some(rhs);
            }

            used_op_fac1.push(infos[0].fac1);
        }

        // UNWRAP: if we can't compare two numbers there's a bug
        op_fac1.sort_by(|a, b| a.partial_cmp(b).unwrap_or_else(|| unreachable!()));

        // make sure we've evolved all slices
        if let Some(muf2) = grid_fac1.into_iter().find(|&grid_mu2| {
            !used_op_fac1
                .iter()
                .any(|&eko_mu2| approx_eq!(f64, grid_mu2, eko_mu2, ulps = EVOLVE_INFO_TOL_ULPS))
        }) {
            return Err(GridError::EvolutionFailure(format!(
                "no operator for muf2 = {muf2} found in {op_fac1:?}"
            )));
        }

        // TODO: convert this unwrap into error
        let grid = lhs.unwrap();

        // UNWRAP: merging evolved slices should be a proper FkTable again
        Ok(FkTable::try_from(grid).unwrap_or_else(|_| unreachable!()))
    }

    /// Deletes bins with the corresponding `bin_indices`. Repeated indices and indices larger or
    /// equal the bin length are ignored.
    pub fn delete_bins(&mut self, bin_indices: &[usize]) {
        let mut bin_indices: Vec<_> = bin_indices
            .iter()
            .copied()
            // ignore indices corresponding to bin that don't exist
            .filter(|&index| index < self.bin_info().bins())
            .collect();

        // sort and remove repeated indices
        bin_indices.sort_unstable();
        bin_indices.dedup();
        let bin_indices = bin_indices;

        let mut bin_ranges: Vec<Range<_>> = Vec::new();

        // convert indices into consecutive ranges
        for &bin_index in &bin_indices {
            match bin_ranges.last_mut() {
                Some(range) if range.end == bin_index => range.end += 1,
                _ => bin_ranges.push(bin_index..(bin_index + 1)),
            }
        }

        let bin_ranges = bin_ranges;
        let mut ranges = bin_ranges.as_slice();
        let old_limits = self.bin_limits.limits();

        // remove the bins from the right first, so as not to invalidate any indices
        if let Some((range, remainder)) = ranges.split_last() {
            if range.end == self.bin_info().bins() {
                self.bin_limits.delete_bins_right(range.end - range.start);
                ranges = remainder;
            }
        }

        // indices on the left aren't affected by removal of bins to their right
        if let Some((range, remainder)) = ranges.split_first() {
            if range.start == 0 {
                self.bin_limits.delete_bins_left(range.end);
                ranges = remainder;
            }
        }

        if !ranges.is_empty() {
            // if there's no remapper we need to store the bin limits in a new remapper
            if self.remapper_mut().is_none() {
                self.set_remapper(
                    BinRemapper::new(
                        old_limits.windows(2).map(|win| win[1] - win[0]).collect(),
                        old_limits.windows(2).map(|win| (win[0], win[1])).collect(),
                    )
                    .unwrap_or_else(|_| unreachable!()),
                )
                .unwrap_or_else(|_| unreachable!());
            }

            // the following should not be needed, but let's set these limits to integer values
            self.bin_limits = BinLimits::new(
                iter::successors(Some(0.0), |x| Some(x + 1.0))
                    .take(old_limits.len() - bin_indices.len())
                    .collect(),
            );
        }

        if let Some(remapper) = self.remapper_mut() {
            remapper.delete_bins(&bin_ranges);
        }

        for &bin_index in bin_indices.iter().rev() {
            self.subgrids.remove_index(Axis(1), bin_index);
        }
    }

    /// Change the particle ID convention.
    pub fn rotate_pid_basis(&mut self, pid_basis: PidBasis) {
        match (self.pid_basis(), pid_basis) {
            (PidBasis::Pdg, PidBasis::Evol) => {
                self.channels = self
                    .channels()
                    .iter()
                    .map(|channel| Channel::translate(channel, &pids::pdg_mc_pids_to_evol))
                    .collect();

                self.set_pid_basis(PidBasis::Evol);
            }
            (PidBasis::Evol, PidBasis::Pdg) => {
                self.channels = self
                    .channels()
                    .iter()
                    .map(|channel| Channel::translate(channel, &pids::evol_to_pdg_mc_ids))
                    .collect();

                self.set_pid_basis(PidBasis::Pdg);
            }
            (PidBasis::Evol, PidBasis::Evol) | (PidBasis::Pdg, PidBasis::Pdg) => {
                // here's nothing to do
            }
        }
    }

    /// Deletes channels with the corresponding `channel_indices`. Repeated indices and indices
    /// larger or equal than the number of channels are ignored.
    pub fn delete_channels(&mut self, channel_indices: &[usize]) {
        let mut channel_indices: Vec<_> = channel_indices
            .iter()
            .copied()
            // ignore indices corresponding to bin that don't exist
            .filter(|&index| index < self.channels().len())
            .collect();

        // sort and remove repeated indices
        channel_indices.sort_unstable();
        channel_indices.dedup();
        channel_indices.reverse();
        let channel_indices = channel_indices;

        for index in channel_indices {
            self.channels.remove(index);
            self.subgrids.remove_index(Axis(2), index);
        }
    }

    pub(crate) fn rewrite_channels(&mut self, add: &[(i32, i32)], del: &[i32]) {
        self.channels = self
            .channels()
            .iter()
            .map(|entry| {
                Channel::new(
                    entry
                        .entry()
                        .iter()
                        .map(|(a, b, f)| {
                            (
                                // if `a` is to be added to another pid replace it with this pid
                                add.iter().fold(
                                    *a,
                                    |id, &(source, target)| if id == source { target } else { id },
                                ),
                                // if `b` is to be added to another pid replace it with this pid
                                add.iter().fold(
                                    *b,
                                    |id, &(source, target)| if id == source { target } else { id },
                                ),
                                // if any of the pids `a` or `b` are to b deleted set the factor to
                                // zero
                                if del.iter().any(|id| id == a || id == b) {
                                    0.0
                                } else {
                                    *f
                                },
                            )
                        })
                        .collect(),
                )
            })
            .collect();
    }

    /// Splits the grid such that each channel contains only a single tuple of PIDs.
    pub fn split_channels(&mut self) {
        let indices: Vec<_> = self
            .channels()
            .iter()
            .enumerate()
            .flat_map(|(index, entry)| iter::repeat(index).take(entry.entry().len()))
            .collect();

        self.subgrids = self.subgrids.select(Axis(2), &indices);
        self.channels = self
            .channels()
            .iter()
            .flat_map(|entry| {
                entry
                    .entry()
                    .iter()
                    .copied()
                    .map(move |entry| Channel::new(vec![entry]))
            })
            .collect();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::channel;
    use std::fs::File;

    #[test]
    fn grid_with_subgrid_type() {
        let subgrid_type = String::from("Idontexist");
        let result = Grid::with_subgrid_type(
            vec![],
            vec![],
            vec![],
            SubgridParams::default(),
            ExtraSubgridParams::default(),
            &subgrid_type,
        );

        matches!(result, Err(GridError::UnknownSubgridType(x)) if x == subgrid_type);
    }

    #[test]
    fn grid_merge_empty_subgrids() {
        let mut grid = Grid::new(
            vec![
                channel![2, 2, 1.0; 4, 4, 1.0],
                channel![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let other = Grid::new(
            vec![
                // differently ordered than `grid`
                channel![1, 1, 1.0; 3, 3, 1.0],
                channel![2, 2, 1.0; 4, 4, 1.0],
            ],
            vec![Order::new(1, 2, 0, 0), Order::new(1, 2, 0, 1)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        // merging with empty subgrids should not change the grid
        grid.merge(other).unwrap();

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);
    }

    #[test]
    fn grid_merge_orders() {
        let mut grid = Grid::new(
            vec![
                channel![2, 2, 1.0; 4, 4, 1.0],
                channel![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            vec![
                channel![2, 2, 1.0; 4, 4, 1.0],
                channel![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![
                Order::new(1, 2, 0, 0),
                Order::new(1, 2, 0, 1),
                Order::new(0, 2, 0, 0),
            ],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        other.fill_all(
            0,
            0.1,
            &Ntuple {
                x1: 0.1,
                x2: 0.2,
                q2: 90.0_f64.powi(2),
                weight: (),
            },
            &[1.0, 2.0],
        );
        other.fill_all(
            1,
            0.1,
            &Ntuple {
                x1: 0.1,
                x2: 0.2,
                q2: 90.0_f64.powi(2),
                weight: (),
            },
            &[1.0, 2.0],
        );

        // merge with four non-empty subgrids
        grid.merge(other).unwrap();

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 3);
    }

    #[test]
    fn grid_merge_channels_entries() {
        let mut grid = Grid::new(
            vec![
                channel![2, 2, 1.0; 4, 4, 1.0],
                channel![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            vec![channel![22, 22, 1.0], channel![2, 2, 1.0; 4, 4, 1.0]],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        // fill the photon-photon entry
        other.fill(
            0,
            0.1,
            0,
            &Ntuple {
                x1: 0.1,
                x2: 0.2,
                q2: 90.0_f64.powi(2),
                weight: 3.0,
            },
        );

        grid.merge(other).unwrap();

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 3);
        assert_eq!(grid.orders().len(), 1);
    }

    #[test]
    fn grid_merge_bins() {
        let mut grid = Grid::new(
            vec![
                channel![2, 2, 1.0; 4, 4, 1.0],
                channel![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.0, 0.25, 0.5],
            SubgridParams::default(),
        );

        assert_eq!(grid.bin_info().bins(), 2);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            vec![
                // channels are differently sorted
                channel![1, 1, 1.0; 3, 3, 1.0],
                channel![2, 2, 1.0; 4, 4, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        other.fill_all(
            0,
            0.1,
            &Ntuple {
                x1: 0.1,
                x2: 0.2,
                q2: 90.0_f64.powi(2),
                weight: (),
            },
            &[2.0, 3.0],
        );

        grid.merge(other).unwrap();

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);
    }

    // TODO: convolve_subgrid, merge_bins, subgrid, set_subgrid

    #[test]
    fn grid_convolutions() {
        let mut grid = Grid::new(
            vec![channel![21, 21, 1.0]],
            vec![Order {
                alphas: 0,
                alpha: 0,
                logxir: 0,
                logxif: 0,
            }],
            vec![0.0, 1.0],
            SubgridParams::default(),
        );

        // by default we assume unpolarized proton PDFs are used
        assert_eq!(
            grid.convolutions(),
            [Convolution::UnpolPDF(2212), Convolution::UnpolPDF(2212)]
        );

        grid.set_convolution(0, Convolution::UnpolPDF(-2212));
        grid.set_convolution(1, Convolution::UnpolPDF(-2212));

        assert_eq!(
            grid.convolutions(),
            [Convolution::UnpolPDF(-2212), Convolution::UnpolPDF(-2212)]
        );
    }

    #[test]
    fn evolve_info() {
        let grid =
            Grid::read(File::open("../test-data/LHCB_WP_7TEV.pineappl.lz4").unwrap()).unwrap();
        let info = grid.evolve_info(&[]);

        assert_eq!(info.fac1.len(), 1);
        assert_approx_eq!(f64, info.fac1[0], 6456.443904000001, ulps = 64);

        assert_eq!(info.pids1, [-3, -1, 2, 4, 21, 22]);

        assert_eq!(info.x1.len(), 50);
        assert_approx_eq!(f64, info.x1[0], 1.9999999999999954e-7, ulps = 64);

        assert_approx_eq!(f64, info.x1[1], 3.034304765867952e-7, ulps = 64);
        assert_approx_eq!(f64, info.x1[2], 4.6035014748963906e-7, ulps = 64);
        assert_approx_eq!(f64, info.x1[3], 6.984208530700364e-7, ulps = 64);
        assert_approx_eq!(f64, info.x1[4], 1.0596094959101024e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[5], 1.607585498470808e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[6], 2.438943292891682e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[7], 3.7002272069854957e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[8], 5.613757716930151e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[9], 8.516806677573355e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[10], 1.292101569074731e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[11], 1.9602505002391748e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[12], 2.97384953722449e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[13], 4.511438394964044e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[14], 6.843744918967896e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[15], 0.00010381172986576898, ulps = 64);
        assert_approx_eq!(f64, info.x1[16], 0.00015745605600841445, ulps = 64);
        assert_approx_eq!(f64, info.x1[17], 0.00023878782918561914, ulps = 64);
        assert_approx_eq!(f64, info.x1[18], 0.00036205449638139736, ulps = 64);
        assert_approx_eq!(f64, info.x1[19], 0.0005487795323670796, ulps = 64);
        assert_approx_eq!(f64, info.x1[20], 0.0008314068836488144, ulps = 64);
        assert_approx_eq!(f64, info.x1[21], 0.0012586797144272762, ulps = 64);
        assert_approx_eq!(f64, info.x1[22], 0.0019034634022867384, ulps = 64);
        assert_approx_eq!(f64, info.x1[23], 0.0028738675812817515, ulps = 64);
        assert_approx_eq!(f64, info.x1[24], 0.004328500638820811, ulps = 64);
        assert_approx_eq!(f64, info.x1[25], 0.006496206194633799, ulps = 64);
        assert_approx_eq!(f64, info.x1[26], 0.009699159574043398, ulps = 64);
        assert_approx_eq!(f64, info.x1[27], 0.014375068581090129, ulps = 64);
        assert_approx_eq!(f64, info.x1[28], 0.02108918668378717, ulps = 64);
        assert_approx_eq!(f64, info.x1[29], 0.030521584007828916, ulps = 64);
        assert_approx_eq!(f64, info.x1[30], 0.04341491741702269, ulps = 64);
        assert_approx_eq!(f64, info.x1[31], 0.060480028754447364, ulps = 64);
        assert_approx_eq!(f64, info.x1[32], 0.08228122126204893, ulps = 64);
        assert_approx_eq!(f64, info.x1[33], 0.10914375746330703, ulps = 64);
        assert_approx_eq!(f64, info.x1[34], 0.14112080644440345, ulps = 64);
        assert_approx_eq!(f64, info.x1[35], 0.17802566042569432, ulps = 64);
        assert_approx_eq!(f64, info.x1[36], 0.2195041265003886, ulps = 64);
        assert_approx_eq!(f64, info.x1[37], 0.2651137041582823, ulps = 64);
        assert_approx_eq!(f64, info.x1[38], 0.31438740076927585, ulps = 64);
        assert_approx_eq!(f64, info.x1[39], 0.3668753186482242, ulps = 64);
        assert_approx_eq!(f64, info.x1[40], 0.4221667753589648, ulps = 64);
        assert_approx_eq!(f64, info.x1[41], 0.4798989029610255, ulps = 64);
        assert_approx_eq!(f64, info.x1[42], 0.5397572337880445, ulps = 64);
        assert_approx_eq!(f64, info.x1[43], 0.601472197967335, ulps = 64);
        assert_approx_eq!(f64, info.x1[44], 0.6648139482473823, ulps = 64);
        assert_approx_eq!(f64, info.x1[45], 0.7295868442414312, ulps = 64);
        assert_approx_eq!(f64, info.x1[46], 0.7956242522922756, ulps = 64);
        assert_approx_eq!(f64, info.x1[47], 0.8627839323906108, ulps = 64);
        assert_approx_eq!(f64, info.x1[48], 0.9309440808717544, ulps = 64);
        assert_approx_eq!(f64, info.x1[49], 1.0, ulps = 64);

        assert_eq!(info.ren1.len(), 1);
        assert_approx_eq!(f64, info.ren1[0], 6456.443904000001, ulps = 64);
    }
}
