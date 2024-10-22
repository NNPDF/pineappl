//! Module containing all traits and supporting structures for grids.

use super::bin::{BinInfo, BinLimits, BinRemapper};
use super::boc::{Channel, Kinematics, Order, Scales};
use super::convolutions::{Conv, ConvType, ConvolutionCache};
use super::empty_subgrid::EmptySubgridV1;
use super::evolution::{self, AlphasTable, EvolveInfo, OperatorSliceInfo};
use super::fk_table::FkTable;
use super::import_subgrid::ImportSubgridV1;
use super::interp_subgrid::InterpSubgridV1;
use super::interpolation::Interp;
use super::pids::PidBasis;
use super::subgrid::{Subgrid, SubgridEnum};
use super::v0;
use bitflags::bitflags;
use float_cmp::{approx_eq, assert_approx_eq};
use git_version::git_version;
use itertools::Itertools;
use lz4_flex::frame::{FrameDecoder, FrameEncoder};
use ndarray::{s, Array3, ArrayD, ArrayView3, ArrayViewMut3, Axis, CowArray, Dimension, Ix4};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::iter;
use std::mem;
use std::ops::Range;
use thiserror::Error;

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
    #[error("file version {file_version} is not supported")]
    FileVersionUnsupported {
        /// File format version of the file read.
        file_version: u64,
    },
    /// Returned from [`Grid::evolve_with_slice_iter`] if the evolution failed.
    #[error("failed to evolve grid: {0}")]
    EvolutionFailure(String),
    /// Errors that do no originate from this crate itself.
    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

#[derive(Clone, Deserialize, Serialize)]
pub(crate) struct Mmv4;

fn default_metadata() -> BTreeMap<String, String> {
    iter::once((
        "pineappl_gitversion".to_owned(),
        git_version!(
            args = ["--always", "--dirty", "--long", "--tags"],
            cargo_prefix = "cargo:",
            fallback = "unknown"
        )
        .to_owned(),
    ))
    .collect()
}

#[derive(Clone, Deserialize, Serialize)]
pub(crate) enum MoreMembers {
    V4(Mmv4),
}

bitflags! {
    /// Bitflags for optimizing a [`Grid`]. See [`Grid::optimize_using`].
    #[derive(Clone, Copy)]
    #[repr(transparent)]
    pub struct GridOptFlags: u32 {
        /// Change the [`Subgrid`] type to optimize storage efficiency.
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
    pub(crate) subgrids: Array3<SubgridEnum>,
    pub(crate) channels: Vec<Channel>,
    pub(crate) bin_limits: BinLimits,
    pub(crate) orders: Vec<Order>,
    pub(crate) metadata: BTreeMap<String, String>,
    pub(crate) convolutions: Vec<Conv>,
    pub(crate) pid_basis: PidBasis,
    pub(crate) more_members: MoreMembers,
    pub(crate) kinematics: Vec<Kinematics>,
    pub(crate) interps: Vec<Interp>,
    pub(crate) remapper: Option<BinRemapper>,
    pub(crate) scales: Scales,
}

impl Grid {
    /// Constructor.
    ///
    /// # Panics
    ///
    /// Panics when the number of PIDs in `channels` is not equal to `convolutions.len()`, or
    /// `interps` and `kinematics` have different lengths or if `kinematics` are not compatible
    /// with `scales`.
    #[must_use]
    pub fn new(
        pid_basis: PidBasis,
        channels: Vec<Channel>,
        orders: Vec<Order>,
        bin_limits: Vec<f64>,
        convolutions: Vec<Conv>,
        interps: Vec<Interp>,
        kinematics: Vec<Kinematics>,
        scales: Scales,
    ) -> Self {
        for (channel_idx, channel) in channels.iter().enumerate() {
            let offending_entry = channel
                .entry()
                .iter()
                .find_map(|(pids, _)| (pids.len() != convolutions.len()).then_some(pids.len()));

            if let Some(pids_len) = offending_entry {
                panic!("channel #{channel_idx} has wrong number of PIDs: expected {}, found {pids_len}", convolutions.len());
            }
        }

        assert_eq!(
            interps.len(),
            kinematics.len(),
            "interps and kinematics have different lengths: {} vs. {}",
            interps.len(),
            kinematics.len(),
        );

        assert!(
            scales.compatible_with(&kinematics),
            "scales and kinematics are not compatible"
        );

        Self {
            subgrids: Array3::from_shape_simple_fn(
                (orders.len(), bin_limits.len() - 1, channels.len()),
                || EmptySubgridV1.into(),
            ),
            orders,
            bin_limits: BinLimits::new(bin_limits),
            metadata: default_metadata(),
            more_members: MoreMembers::V4(Mmv4),
            convolutions,
            pid_basis,
            channels,
            interps,
            kinematics,
            remapper: None,
            scales,
        }
    }

    /// Return the convention by which the channels' PIDs are encoded.
    #[must_use]
    pub const fn pid_basis(&self) -> &PidBasis {
        &self.pid_basis
    }

    /// Set the convention by which PIDs of channels are interpreted.
    pub fn pid_basis_mut(&mut self) -> &mut PidBasis {
        &mut self.pid_basis
    }

    /// TODO
    #[must_use]
    pub fn kinematics(&self) -> &[Kinematics] {
        &self.kinematics
    }

    /// TODO
    #[must_use]
    pub const fn scales(&self) -> &Scales {
        &self.scales
    }

    /// Perform a convolution using the PDFs and strong coupling in `lumi_cache`, and
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
        cache: &mut ConvolutionCache,
        order_mask: &[bool],
        bin_indices: &[usize],
        channel_mask: &[bool],
        xi: &[(f64, f64, f64)],
    ) -> Vec<f64> {
        let mut cache = cache.new_grid_conv_cache(self, xi);

        let bin_indices = if bin_indices.is_empty() {
            (0..self.bin_info().bins()).collect()
        } else {
            bin_indices.to_vec()
        };
        let mut bins = vec![0.0; bin_indices.len() * xi.len()];
        let normalizations = self.bin_info().normalizations();
        let pdg_channels = self.channels_pdg();

        for (xi_index, &xis @ (xir, xif, xia)) in xi.iter().enumerate() {
            for ((ord, bin, chan), subgrid) in self.subgrids.indexed_iter() {
                let order = &self.orders[ord];

                if ((order.logxir > 0) && (xir == 1.0))
                    || ((order.logxif > 0) && (xif == 1.0))
                    || ((order.logxia > 0) && (xia == 1.0))
                {
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
                let mut value = 0.0;

                cache.set_grids(self, subgrid, xis);

                for (idx, v) in subgrid.indexed_iter() {
                    let mut lumi = 0.0;

                    for entry in channel.entry() {
                        // TODO: we assume `idx` to be ordered as scale, x1, x2
                        let fx_prod = cache.as_fx_prod(&entry.0, order.alphas, &idx);
                        lumi += fx_prod * entry.1;
                    }

                    value += lumi * v;
                }

                if order.logxir > 0 {
                    value *= (xir * xir).ln().powi(order.logxir.into());
                }

                if order.logxif > 0 {
                    value *= (xif * xif).ln().powi(order.logxif.into());
                }

                if order.logxia > 0 {
                    value *= (xia * xia).ln().powi(order.logxia.into());
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
        cache: &mut ConvolutionCache,
        ord: usize,
        bin: usize,
        channel: usize,
        xi @ (xir, xif, xia): (f64, f64, f64),
    ) -> ArrayD<f64> {
        let mut cache = cache.new_grid_conv_cache(self, &[(xir, xif, xia)]);

        let normalizations = self.bin_info().normalizations();
        let pdg_channels = self.channels_pdg();

        let subgrid = &self.subgrids[[ord, bin, channel]];
        let order = &self.orders[ord];

        let channel = &pdg_channels[channel];

        cache.set_grids(self, subgrid, xi);

        let node_values: Vec<_> = subgrid.node_values();
        // TODO: generalize this to N dimensions
        assert_eq!(node_values.len(), 3);
        let dim: Vec<_> = node_values.iter().map(Vec::len).collect();

        let mut array = ArrayD::zeros(dim);

        for (idx, value) in subgrid.indexed_iter() {
            assert_eq!(idx.len(), 3);
            let mut lumi = 0.0;

            for entry in channel.entry() {
                debug_assert_eq!(entry.0.len(), 2);
                // TODO: we assume `idx` to be ordered as scale, x1, x2
                let fx_prod = cache.as_fx_prod(&entry.0, order.alphas, &idx);
                lumi += fx_prod * entry.1;
            }

            array[idx.as_slice()] = lumi * value;
        }

        if order.logxir > 0 {
            array *= (xir * xir).ln().powi(order.logxir.into());
        }

        if order.logxif > 0 {
            array *= (xif * xif).ln().powi(order.logxif.into());
        }

        if order.logxia > 0 {
            array *= (xia * xia).ln().powi(order.logxia.into());
        }

        array /= normalizations[bin];
        array
    }

    /// Fills the grid with an ntuple for the given `order`, `observable`, and `channel`. The
    /// parameter `ntuple` must contain the variables specified by the `kinematics` parameter in
    /// the constructor [`Grid::new`] in the same order.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn fill(
        &mut self,
        order: usize,
        observable: f64,
        channel: usize,
        ntuple: &[f64],
        weight: f64,
    ) {
        if let Some(bin) = self.bin_limits.index(observable) {
            let subgrid = &mut self.subgrids[[order, bin, channel]];
            if let SubgridEnum::EmptySubgridV1(_) = subgrid {
                *subgrid = InterpSubgridV1::new(&self.interps).into();
            }

            subgrid.fill(&self.interps, ntuple, weight);
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

        match file_version {
            0 => v0::read_uncompressed_v0(reader),
            1 => bincode::deserialize_from(reader).map_err(GridError::ReadFailure),
            _ => Err(GridError::FileVersionUnsupported { file_version }),
        }
    }

    /// Serializes `self` into `writer`. Writing is buffered.
    ///
    /// # Errors
    ///
    /// If writing fails an error is returned.
    pub fn write(&self, writer: impl Write) -> Result<(), GridError> {
        let mut writer = BufWriter::new(writer);
        let file_header = b"PineAPPL\x01\0\0\0\0\0\0\0";

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

    /// Return the channels for this `Grid`.
    #[must_use]
    pub fn channels(&self) -> &[Channel] {
        &self.channels
    }

    fn channels_pdg(&self) -> Vec<Channel> {
        self.channels()
            .iter()
            .cloned()
            .map(|channel| self.pid_basis().translate(PidBasis::Pdg, channel))
            .collect()
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
                    new_subgrid.merge(subgrid, None);
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
                self.subgrids[[self_i, self_j, self_k]].merge(subgrid, None);
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
    pub fn convolutions(&self) -> &[Conv] {
        &self.convolutions
    }

    /// Return the convolution types.
    pub fn convolutions_mut(&mut self) -> &mut [Conv] {
        &mut self.convolutions
    }

    /// Charge conjugate both the convolution function with index `convolution` and the PIDs in the
    /// channel definition corresponding to it. This leaves the the results returned by
    /// [`Grid::convolve`] invariant.
    pub fn charge_conjugate(&mut self, convolution: usize) {
        let pid_basis = *self.pid_basis();

        for channel in self.channels_mut() {
            *channel = Channel::new(
                channel
                    .entry()
                    .iter()
                    .cloned()
                    .map(|(mut pids, f)| {
                        let (cc_pid, f1) = pid_basis.charge_conjugate(pids[convolution]);
                        pids[convolution] = cc_pid;
                        (pids, f * f1)
                    })
                    .collect(),
            );
        }

        self.convolutions_mut()[convolution] = self.convolutions()[convolution].cc();
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

        for (index, subgrid) in self.subgrids.indexed_iter_mut() {
            mem::swap(&mut new_subgrids[<[usize; 3]>::from(index)], subgrid);
        }

        self.subgrids = new_subgrids;
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
        logxia: f64,
        global: f64,
    ) {
        for ((i, _, _), subgrid) in self.subgrids.indexed_iter_mut() {
            let order = &self.orders[i];
            let factor = global
                * alphas.powi(order.alphas.into())
                * alpha.powi(order.alpha.into())
                * logxir.powi(order.logxir.into())
                * logxif.powi(order.logxif.into())
                * logxia.powi(order.logxia.into());

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

        self.remapper = Some(remapper);

        Ok(())
    }

    /// Return the currently set remapper, if there is any.
    #[must_use]
    pub const fn remapper(&self) -> Option<&BinRemapper> {
        self.remapper.as_ref()
    }

    fn remapper_mut(&mut self) -> Option<&mut BinRemapper> {
        self.remapper.as_mut()
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
                _ => {
                    // TODO: this requires a `pub(crate)` in `InterpSubgridV1`; we should
                    // replace this with a method
                    if !static_scale_detection {
                        if let SubgridEnum::InterpSubgridV1(subgrid) = subgrid {
                            // disable static-scale detection
                            subgrid.static_q2 = -1.0;
                        }
                    }

                    *subgrid = ImportSubgridV1::from(&*subgrid).into();
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
                            *lhs = mem::replace(rhs, EmptySubgridV1.into());
                        } else {
                            lhs.merge(rhs, None);
                            *rhs = EmptySubgridV1.into();
                        }
                    }
                }
            }
        }
    }

    fn strip_empty_channels(&mut self) {
        let mut indices: Vec<_> = (0..self.channels().len()).collect();

        while let Some(index) = indices.pop() {
            if self
                .subgrids
                .slice(s![.., .., index])
                .iter()
                .all(Subgrid::is_empty)
            {
                self.channels.remove(index);
                self.subgrids.remove_index(Axis(2), index);
            }
        }
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
        let pairs: Vec<_> = self
            .convolutions()
            .iter()
            .enumerate()
            .tuple_combinations()
            .filter(|((_, conv_a), (_, conv_b))| conv_a == conv_b)
            .map(|((idx_a, _), (idx_b, _))| (idx_a, idx_b))
            .collect();

        let (idx_a, idx_b) = match *pairs.as_slice() {
            [] => return,
            [pair] => pair,
            _ => panic!("more than two equal convolutions found"),
        };
        let a_subgrid = self
            .kinematics()
            .iter()
            .position(|&kin| kin == Kinematics::X(idx_a))
            // UNWRAP: should be guaranteed by the constructor
            .unwrap();
        let b_subgrid = self
            .kinematics()
            .iter()
            .position(|&kin| kin == Kinematics::X(idx_b))
            // UNWRAP: should be guaranteed by the constructor
            .unwrap();

        let mut indices: Vec<usize> = (0..self.channels.len()).rev().collect();

        while let Some(index) = indices.pop() {
            let channel_entry = &self.channels[index];

            if *channel_entry == channel_entry.transpose(idx_a, idx_b) {
                // check if in all cases the limits are compatible with merging
                self.subgrids
                    .slice_mut(s![.., .., index])
                    .iter_mut()
                    .for_each(|subgrid| {
                        if !subgrid.is_empty()
                            && (subgrid.node_values()[a_subgrid]
                                == subgrid.node_values()[b_subgrid])
                        {
                            subgrid.symmetrize(a_subgrid, b_subgrid);
                        }
                    });
            } else if let Some((j, &other_index)) = indices
                .iter()
                .enumerate()
                .find(|(_, i)| self.channels[**i] == channel_entry.transpose(idx_a, idx_b))
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
                            *lhs = mem::replace(rhs, EmptySubgridV1.into());
                            // transpose `lhs`
                            todo!();
                        } else {
                            lhs.merge(rhs, Some((a_subgrid, b_subgrid)));
                            *rhs = EmptySubgridV1.into();
                        }
                    }
                }
            }
        }
    }

    /// Upgrades the internal data structures to their latest versions.
    pub fn upgrade(&mut self) {}

    /// Return the metadata of this grid.
    #[must_use]
    pub const fn metadata(&self) -> &BTreeMap<String, String> {
        &self.metadata
    }

    /// Return the metadata of this grid.
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
    pub fn metadata_mut(&mut self) -> &mut BTreeMap<String, String> {
        &mut self.metadata
    }

    /// Returns information for the generation of evolution operators that are being used in
    /// [`Grid::convolve`] with the parameter `order_mask`.
    #[must_use]
    pub fn evolve_info(&self, order_mask: &[bool]) -> EvolveInfo {
        use super::evolution::EVOLVE_INFO_TOL_ULPS;

        let mut ren1 = Vec::new();
        let mut fac1 = Vec::new();
        let mut x1 = Vec::new();
        let mut pids1 = Vec::new();

        for (channel, subgrid) in self
            .subgrids()
            .indexed_iter()
            .filter_map(|(tuple, subgrid)| {
                (!subgrid.is_empty() && (order_mask.is_empty() || order_mask[tuple.0]))
                    .then_some((&self.channels()[tuple.2], subgrid))
            })
        {
            // ren1.extend(subgrid.mu2_grid().iter().map(|Mu2 { ren, .. }| *ren));
            ren1.extend(
                self.kinematics()
                    .iter()
                    .zip(subgrid.node_values())
                    .find_map(|(kin, node_values)| {
                        matches!(kin, &Kinematics::Scale(idx) if idx == 0).then_some(node_values)
                    })
                    // TODO: convert this into an error
                    .unwrap()
                    .into_iter(),
            );
            ren1.sort_by(f64::total_cmp);
            ren1.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = EVOLVE_INFO_TOL_ULPS));

            // fac1.extend(subgrid.mu2_grid().iter().map(|Mu2 { fac, .. }| *fac));
            fac1.extend(
                self.kinematics()
                    .iter()
                    .zip(subgrid.node_values())
                    .find_map(|(kin, node_values)| {
                        matches!(kin, &Kinematics::Scale(idx) if idx == 0).then_some(node_values)
                    })
                    // TODO: convert this into an error
                    .unwrap()
                    .into_iter(),
            );
            fac1.sort_by(f64::total_cmp);
            fac1.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = EVOLVE_INFO_TOL_ULPS));

            x1.extend(
                subgrid
                    .node_values()
                    .iter()
                    .zip(self.kinematics())
                    .filter_map(|(nv, kin)| matches!(kin, Kinematics::X(_)).then(|| nv))
                    .flatten(),
            );

            x1.sort_by(f64::total_cmp);
            x1.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = EVOLVE_INFO_TOL_ULPS));

            for (index, _) in self.convolutions().iter().enumerate() {
                pids1.extend(channel.entry().iter().map(|(pids, _)| pids[index]));
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

    // TODO:
    // - try to find a better solution than to require that E must be convertible into
    //   anyhow::Error

    /// Convert this `Grid` into an [`FkTable`] using `slices.len()` evolution operators, which for
    /// each entry must iterate over a [`Result`] of tuples of an [`OperatorSliceInfo`] and the
    /// corresponding sliced operator. The parameter `order_mask` can be used to include or exclude
    /// orders from this operation, and must correspond to the ordering given by [`Grid::orders`].
    /// Orders that are not given are enabled, and in particular if `order_mask` is empty all
    /// orders are activated.
    ///
    /// # Errors
    ///
    /// Returns a [`GridError::EvolutionFailure`] if either the `operator` or its `info` is
    /// incompatible with this `Grid`. Returns a [`GridError::Other`] if the iterator from `slices`
    /// return an error.
    ///
    /// # Panics
    ///
    /// Panics when the operators returned by either slice have different dimensions than promised
    /// by the corresponding [`OperatorSliceInfo`].
    pub fn evolve<
        'a,
        E: Into<anyhow::Error>,
        S: IntoIterator<Item = Result<(OperatorSliceInfo, CowArray<'a, f64, Ix4>), E>>,
    >(
        &self,
        slices: Vec<S>,
        order_mask: &[bool],
        xi: (f64, f64, f64),
        alphas_table: &AlphasTable,
    ) -> Result<FkTable, GridError> {
        struct Iter<T> {
            iters: Vec<T>,
        }

        impl<T: Iterator> Iterator for Iter<T> {
            type Item = Vec<T::Item>;

            fn next(&mut self) -> Option<Self::Item> {
                let v: Vec<_> = self.iters.iter_mut().filter_map(Iterator::next).collect();

                if v.len() == self.iters.len() {
                    Some(v)
                } else {
                    None
                }
            }
        }

        fn zip_n<O, T>(iters: O) -> impl Iterator<Item = Vec<T::Item>>
        where
            O: IntoIterator<Item = T>,
            T: IntoIterator,
        {
            Iter {
                iters: iters.into_iter().map(IntoIterator::into_iter).collect(),
            }
        }

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
            // TODO: also take care of the fragmentation scale
            .map(|fac| xi.1 * xi.1 * fac)
            .collect();

        let mut perm = Vec::new();

        for result in zip_n(slices) {
            let (infos, operators): (Vec<OperatorSliceInfo>, Vec<CowArray<'_, f64, Ix4>>) = result
                .into_iter()
                .map(|res| res.map_err(|err| GridError::Other(err.into())))
                .collect::<Result<_, _>>()?;

            let (info_0, infos_rest) = infos
                .split_first()
                // UNWRAP: TODO
                .unwrap();

            let dim_op_info_0 = (
                info_0.pids1.len(),
                info_0.x1.len(),
                info_0.pids0.len(),
                info_0.x0.len(),
            );

            assert_eq!(
                operators[0].dim(),
                dim_op_info_0,
                "operator information {:?} does not match the operator's dimensions: {:?}",
                dim_op_info_0,
                operators[0].dim(),
            );

            for (index, info) in infos_rest.iter().enumerate() {
                // TODO: what if the scales of the EKOs don't agree? Is there an ordering problem?
                assert_approx_eq!(f64, info_0.fac1, info.fac1, ulps = EVOLVE_INFO_TOL_ULPS);
                assert_eq!(info_0.pid_basis, info.pid_basis);

                let dim_op_info = (
                    info.pids1.len(),
                    info.x1.len(),
                    info.pids0.len(),
                    info.x0.len(),
                );

                assert_eq!(
                    operators[index + 1].dim(),
                    dim_op_info,
                    "operator information {:?} does not match the operator's dimensions: {:?}",
                    dim_op_info,
                    operators[index + 1].dim(),
                );
            }

            if perm.is_empty() {
                let eko_conv_types: Vec<ConvType> =
                    infos.iter().map(|info| info.conv_type).collect();

                perm = self
                    .convolutions()
                    .iter()
                    .enumerate()
                    .map(|(max_idx, conv)| {
                        eko_conv_types
                            .iter()
                            .take(max_idx + 1)
                            .enumerate()
                            .rev()
                            .find_map(|(idx, &eko_conv_type)| {
                                if conv.conv_type() == eko_conv_type {
                                    Some(idx)
                                } else {
                                    None
                                }
                            })
                            // TODO: convert `unwrap` to `Err`
                            .unwrap()
                    })
                    .collect();
            }

            op_fac1.push(info_0.fac1);

            // it's possible that due to small numerical differences we get two slices which are
            // almost the same. We have to skip those in order not to evolve the 'same' slice twice
            if used_op_fac1
                .iter()
                .any(|&fac| approx_eq!(f64, fac, info_0.fac1, ulps = EVOLVE_INFO_TOL_ULPS))
            {
                continue;
            }

            // skip slices that the grid doesn't use
            if !grid_fac1
                .iter()
                .any(|&fac| approx_eq!(f64, fac, info_0.fac1, ulps = EVOLVE_INFO_TOL_ULPS))
            {
                continue;
            }

            let operators: Vec<_> = perm.iter().map(|&idx| operators[idx].view()).collect();
            let infos: Vec<_> = perm.iter().map(|&idx| infos[idx].clone()).collect();

            let (subgrids, channels) = evolution::evolve_slice_with_many(
                self,
                &operators,
                &infos,
                order_mask,
                xi,
                alphas_table,
            )?;

            let rhs = Self {
                subgrids,
                channels,
                bin_limits: self.bin_limits.clone(),
                orders: vec![Order::new(0, 0, 0, 0, 0)],
                interps: self.interps.clone(),
                metadata: self.metadata.clone(),
                convolutions: self.convolutions.clone(),
                pid_basis: infos[0].pid_basis,
                more_members: self.more_members.clone(),
                kinematics: self.kinematics.clone(),
                remapper: self.remapper.clone(),
                // TODO: is this correct?
                scales: self.scales.clone(),
            };

            if let Some(lhs) = &mut lhs {
                lhs.merge(rhs)?;
            } else {
                lhs = Some(rhs);
            }

            used_op_fac1.push(info_0.fac1);
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
        let self_pid_basis = *self.pid_basis();
        for channel in &mut self.channels {
            *channel = self_pid_basis.translate(pid_basis, channel.clone());
        }
        self.pid_basis = pid_basis;
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

    /// Delete orders with the corresponding `order_indices`. Repeated indices and indices larger
    /// or equal than the number of orders are ignored.
    pub fn delete_orders(&mut self, order_indices: &[usize]) {
        let mut order_indices: Vec<_> = order_indices
            .iter()
            .copied()
            // ignore indices corresponding to orders that don't exist
            .filter(|&index| index < self.orders().len())
            .collect();

        // sort and remove repeated indices
        order_indices.sort_unstable();
        order_indices.dedup();
        order_indices.reverse();
        let order_indices = order_indices;

        for index in order_indices {
            self.orders.remove(index);
            self.subgrids.remove_index(Axis(0), index);
        }
    }

    pub(crate) fn rewrite_channels(&mut self, add: &[(i32, i32)], del: &[i32]) {
        // TODO: generalize this method to n convolutions
        assert_eq!(self.convolutions().len(), 2);

        self.channels = self
            .channels()
            .iter()
            .map(|entry| {
                Channel::new(
                    entry
                        .entry()
                        .iter()
                        .map(|(pids, f)| {
                            (
                                vec![
                                    // if `a` is to be added to another pid replace it with this pid
                                    add.iter().fold(pids[0], |id, &(source, target)| {
                                        if id == source {
                                            target
                                        } else {
                                            id
                                        }
                                    }),
                                    // if `b` is to be added to another pid replace it with this pid
                                    add.iter().fold(pids[1], |id, &(source, target)| {
                                        if id == source {
                                            target
                                        } else {
                                            id
                                        }
                                    }),
                                ],
                                // if any of the pids `a` or `b` are to b deleted set the factor to
                                // zero
                                if del.iter().any(|&id| id == pids[0] || id == pids[1]) {
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
                    .cloned()
                    .map(move |entry| Channel::new(vec![entry]))
            })
            .collect();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boc::ScaleFuncForm;
    use crate::channel;
    use crate::convolutions::ConvType;
    use std::fs::File;

    #[test]
    #[should_panic(expected = "channel #0 has wrong number of PIDs: expected 2, found 3")]
    fn grid_new_panic0() {
        let channel = vec![(vec![1, -1, 1], 1.0), (vec![2, -2, 2], 1.0)];

        let _ = Grid::new(
            PidBasis::Pdg,
            vec![Channel::new(channel)],
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![0.0, 1.0],
            vec![
                Conv::new(ConvType::UnpolPDF, 2212),
                Conv::new(ConvType::UnpolPDF, 2212),
            ],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );
    }

    #[test]
    #[should_panic(expected = "interps and kinematics have different lengths: 2 vs. 3")]
    fn grid_new_panic1() {
        let channel = vec![(vec![1, -1], 1.0), (vec![2, -2], 1.0)];

        let _ = Grid::new(
            PidBasis::Pdg,
            vec![Channel::new(channel)],
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![0.0, 1.0],
            vec![
                Conv::new(ConvType::UnpolPDF, 2212),
                Conv::new(ConvType::UnpolPDF, 2212),
            ],
            v0::default_interps(1),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );
    }

    #[test]
    #[should_panic(expected = "scales and kinematics are not compatible")]
    fn grid_new_panic2() {
        let channel = vec![(vec![1, -1], 1.0), (vec![2, -2], 1.0)];

        let _ = Grid::new(
            PidBasis::Pdg,
            vec![Channel::new(channel)],
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![0.0, 1.0],
            vec![
                Conv::new(ConvType::UnpolPDF, 2212),
                Conv::new(ConvType::UnpolPDF, 2212),
            ],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(1),
                frg: ScaleFuncForm::NoScale,
            },
        );
    }

    #[test]
    fn grid_read_file_version_unsupported() {
        assert!(matches!(
            Grid::read(
                &[b'P', b'i', b'n', b'e', b'A', b'P', b'P', b'L', 99, 0, 0, 0, 0, 0, 0, 0][..]
            ),
            Err(GridError::FileVersionUnsupported { file_version: 99 })
        ));
    }

    #[test]
    fn grid_merge_empty_subgrids() {
        let mut grid = Grid::new(
            PidBasis::Pdg,
            vec![
                channel![2, 2, 1.0; 4, 4, 1.0],
                channel![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let other = Grid::new(
            PidBasis::Pdg,
            vec![
                // differently ordered than `grid`
                channel![1, 1, 1.0; 3, 3, 1.0],
                channel![2, 2, 1.0; 4, 4, 1.0],
            ],
            vec![Order::new(1, 2, 0, 0, 0), Order::new(1, 2, 0, 1, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
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
            PidBasis::Pdg,
            vec![
                channel![2, 2, 1.0; 4, 4, 1.0],
                channel![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            PidBasis::Pdg,
            vec![
                channel![2, 2, 1.0; 4, 4, 1.0],
                channel![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![
                Order::new(1, 2, 0, 0, 0),
                Order::new(1, 2, 0, 1, 0),
                Order::new(0, 2, 0, 0, 0),
            ],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        other.fill(0, 0.1, 0, &[90.0_f64.powi(2), 0.1, 0.2], 1.0);
        other.fill(0, 0.1, 1, &[90.0_f64.powi(2), 0.1, 0.2], 2.0);
        other.fill(1, 0.1, 0, &[90.0_f64.powi(2), 0.1, 0.2], 1.0);
        other.fill(1, 0.1, 1, &[90.0_f64.powi(2), 0.1, 0.2], 2.0);

        // merge with four non-empty subgrids
        grid.merge(other).unwrap();

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 3);
    }

    #[test]
    fn grid_merge_channels_entries() {
        let mut grid = Grid::new(
            PidBasis::Pdg,
            vec![
                channel![2, 2, 1.0; 4, 4, 1.0],
                channel![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            PidBasis::Pdg,
            vec![channel![22, 22, 1.0], channel![2, 2, 1.0; 4, 4, 1.0]],
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        // fill the photon-photon entry
        other.fill(0, 0.1, 0, &[90.0_f64.powi(2), 0.1, 0.2], 3.0);

        grid.merge(other).unwrap();

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 3);
        assert_eq!(grid.orders().len(), 1);
    }

    #[test]
    fn grid_merge_bins() {
        let mut grid = Grid::new(
            PidBasis::Pdg,
            vec![
                channel![2, 2, 1.0; 4, 4, 1.0],
                channel![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![0.0, 0.25, 0.5],
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        assert_eq!(grid.bin_info().bins(), 2);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            PidBasis::Pdg,
            vec![
                // channels are differently sorted
                channel![1, 1, 1.0; 3, 3, 1.0],
                channel![2, 2, 1.0; 4, 4, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![0.5, 0.75, 1.0],
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        other.fill(0, 0.1, 0, &[90.0_f64.powi(2), 0.1, 0.2], 2.0);
        other.fill(0, 0.1, 1, &[90.0_f64.powi(2), 0.1, 0.2], 3.0);

        grid.merge(other).unwrap();

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);
    }

    // TODO: convolve_subgrid, merge_bins, subgrid, set_subgrid

    #[test]
    fn grid_convolutions() {
        let mut grid = Grid::new(
            PidBasis::Pdg,
            vec![channel![21, 21, 1.0]],
            vec![Order {
                alphas: 0,
                alpha: 0,
                logxir: 0,
                logxif: 0,
                logxia: 0,
            }],
            vec![0.0, 1.0],
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        // by default we assume unpolarized proton PDFs are used
        assert_eq!(
            grid.convolutions(),
            [
                Conv::new(ConvType::UnpolPDF, 2212),
                Conv::new(ConvType::UnpolPDF, 2212)
            ]
        );

        grid.convolutions_mut()[0] = Conv::new(ConvType::UnpolPDF, -2212);
        grid.convolutions_mut()[1] = Conv::new(ConvType::UnpolPDF, -2212);

        assert_eq!(
            grid.convolutions(),
            [
                Conv::new(ConvType::UnpolPDF, -2212),
                Conv::new(ConvType::UnpolPDF, -2212)
            ]
        );
    }

    #[test]
    fn grid_set_remapper_bin_number_mismatch() {
        let mut grid = Grid::new(
            PidBasis::Pdg,
            vec![
                channel![2, 2, 1.0; 4, 4, 1.0],
                channel![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(2),
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        assert!(matches!(
            grid.set_remapper(BinRemapper::new(vec![1.0], vec![(0.0, 1.0)]).unwrap()),
            Err(GridError::BinNumberMismatch {
                grid_bins: 4,
                remapper_bins: 1
            })
        ));
    }

    #[test]
    fn evolve_info() {
        let grid =
            Grid::read(File::open("../test-data/LHCB_WP_7TEV_opt.pineappl.lz4").unwrap()).unwrap();
        let info = grid.evolve_info(&[]);

        assert_eq!(info.fac1.len(), 1);
        assert_approx_eq!(f64, info.fac1[0], 6456.443904000001, ulps = 64);

        assert_eq!(info.pids1, [-3, -1, 2, 4, 21, 22]);

        assert_eq!(info.x1.len(), 39);
        assert_approx_eq!(f64, info.x1[0], 1.9602505002391748e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[1], 2.97384953722449e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[2], 4.511438394964044e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[3], 6.843744918967896e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[4], 0.00010381172986576898, ulps = 64);
        assert_approx_eq!(f64, info.x1[5], 0.00015745605600841445, ulps = 64);
        assert_approx_eq!(f64, info.x1[6], 0.00023878782918561914, ulps = 64);
        assert_approx_eq!(f64, info.x1[7], 0.00036205449638139736, ulps = 64);
        assert_approx_eq!(f64, info.x1[8], 0.0005487795323670796, ulps = 64);
        assert_approx_eq!(f64, info.x1[9], 0.0008314068836488144, ulps = 64);
        assert_approx_eq!(f64, info.x1[10], 0.0012586797144272762, ulps = 64);
        assert_approx_eq!(f64, info.x1[11], 0.0019034634022867384, ulps = 64);
        assert_approx_eq!(f64, info.x1[12], 0.0028738675812817515, ulps = 64);
        assert_approx_eq!(f64, info.x1[13], 0.004328500638820811, ulps = 64);
        assert_approx_eq!(f64, info.x1[14], 0.006496206194633799, ulps = 64);
        assert_approx_eq!(f64, info.x1[15], 0.009699159574043398, ulps = 64);
        assert_approx_eq!(f64, info.x1[16], 0.014375068581090129, ulps = 64);
        assert_approx_eq!(f64, info.x1[17], 0.02108918668378717, ulps = 64);
        assert_approx_eq!(f64, info.x1[18], 0.030521584007828916, ulps = 64);
        assert_approx_eq!(f64, info.x1[19], 0.04341491741702269, ulps = 64);
        assert_approx_eq!(f64, info.x1[20], 0.060480028754447364, ulps = 64);
        assert_approx_eq!(f64, info.x1[21], 0.08228122126204893, ulps = 64);
        assert_approx_eq!(f64, info.x1[22], 0.10914375746330703, ulps = 64);
        assert_approx_eq!(f64, info.x1[23], 0.14112080644440345, ulps = 64);
        assert_approx_eq!(f64, info.x1[24], 0.17802566042569432, ulps = 64);
        assert_approx_eq!(f64, info.x1[25], 0.2195041265003886, ulps = 64);
        assert_approx_eq!(f64, info.x1[26], 0.2651137041582823, ulps = 64);
        assert_approx_eq!(f64, info.x1[27], 0.31438740076927585, ulps = 64);
        assert_approx_eq!(f64, info.x1[28], 0.3668753186482242, ulps = 64);
        assert_approx_eq!(f64, info.x1[29], 0.4221667753589648, ulps = 64);
        assert_approx_eq!(f64, info.x1[30], 0.4798989029610255, ulps = 64);
        assert_approx_eq!(f64, info.x1[31], 0.5397572337880445, ulps = 64);
        assert_approx_eq!(f64, info.x1[32], 0.601472197967335, ulps = 64);
        assert_approx_eq!(f64, info.x1[33], 0.6648139482473823, ulps = 64);
        assert_approx_eq!(f64, info.x1[34], 0.7295868442414312, ulps = 64);
        assert_approx_eq!(f64, info.x1[35], 0.7956242522922756, ulps = 64);
        assert_approx_eq!(f64, info.x1[36], 0.8627839323906108, ulps = 64);
        assert_approx_eq!(f64, info.x1[37], 0.9309440808717544, ulps = 64);
        assert_approx_eq!(f64, info.x1[38], 1.0, ulps = 64);

        assert_eq!(info.ren1.len(), 1);
        assert_approx_eq!(f64, info.ren1[0], 6456.443904000001, ulps = 64);
    }
}
