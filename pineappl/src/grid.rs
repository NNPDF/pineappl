//! Module containing all traits and supporting structures for grids.

use super::boc::{BinsWithFillLimits, Channel, Kinematics, Order, ScaleFuncForm, Scales};
use super::convolutions::{Conv, ConvType, ConvolutionCache};
use super::empty_subgrid::EmptySubgridV1;
use super::error::{Error, Result};
use super::evolution::{self, AlphasTable, EvolveInfo, OperatorSliceInfo};
use super::fk_table::FkTable;
use super::import_subgrid::ImportSubgridV1;
use super::interp_subgrid::InterpSubgridV1;
use super::interpolation::Interp;
use super::pids::PidBasis;
use super::reference::Reference;
use super::subgrid::{self, Subgrid, SubgridEnum};
use super::v0;
use bitflags::bitflags;
use float_cmp::{approx_eq, assert_approx_eq};
use git_version::git_version;
use itertools::Itertools;
use lz4_flex::frame::{FrameDecoder, FrameEncoder};
use ndarray::{s, Array2, Array3, ArrayView3, ArrayViewMut3, Axis, CowArray, Dimension, Ix4, Zip};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::ops::{Bound, RangeBounds};
use std::{iter, mem};

const BIN_AXIS: Axis = Axis(1);

// const ORDER_AXIS: Axis = Axis(0);
// const CHANNEL_AXIS: Axis = Axis(2);

#[derive(Clone, Deserialize, Serialize)]
struct Mmv4;

#[derive(Clone, Deserialize, Serialize)]
enum MoreMembers {
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
        /// the case, optimize it by undoing the interpolation in the scale.
        const OPTIMIZE_NODES = 0b10;
        /// Deprecated name for [`GridOptFlags::OPTIMIZE_NODES`].
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
    bwfl: BinsWithFillLimits,
    orders: Vec<Order>,
    channels: Vec<Channel>,
    pid_basis: PidBasis,
    convolutions: Vec<Conv>,
    interps: Vec<Interp>,
    kinematics: Vec<Kinematics>,
    scales: Scales,
    metadata: BTreeMap<String, String>,
    more_members: MoreMembers,
    reference: Reference,
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
        bwfl: BinsWithFillLimits,
        orders: Vec<Order>,
        channels: Vec<Channel>,
        pid_basis: PidBasis,
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
                (orders.len(), bwfl.len(), channels.len()),
                || EmptySubgridV1.into(),
            ),
            bwfl,
            orders,
            channels,
            pid_basis,
            convolutions,
            interps,
            kinematics,
            scales,
            metadata: iter::once((
                "pineappl_gitversion".to_owned(),
                git_version!(
                    args = ["--always", "--dirty", "--long", "--tags"],
                    cargo_prefix = "cargo:",
                    fallback = "unknown"
                )
                .to_owned(),
            ))
            .collect(),
            more_members: MoreMembers::V4(Mmv4),
            reference: Reference::default(),
        }
    }

    /// TODO
    pub fn reference(&self) -> &Reference {
        &self.reference
    }

    /// TODO
    pub fn set_reference(&mut self, reference: Reference) {
        // TODO: check that the number of bins and channels is consistent between the grid and
        // `reference`
        self.reference = reference;
    }

    /// Return the convention by which the channels' PIDs are encoded.
    #[must_use]
    pub const fn pid_basis(&self) -> &PidBasis {
        &self.pid_basis
    }

    /// Return a vector containing the interpolation specifications for this grid.
    #[must_use]
    pub fn interpolations(&self) -> &[Interp] {
        &self.interps
    }

    /// Return a vector containing the kinematic specifications for this grid.
    #[must_use]
    pub fn kinematics(&self) -> &[Kinematics] {
        &self.kinematics
    }

    /// Return a vector containg the scale specifications for this grid.
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
            (0..self.bwfl().len()).collect()
        } else {
            bin_indices.to_vec()
        };
        let mut bins = vec![0.0; bin_indices.len() * xi.len()];
        let normalizations = self.bwfl().normalizations();
        let pdg_channels = self.channels_pdg();

        for (xi_index, &xis @ (xir, xif, xia)) in xi.iter().enumerate() {
            for ((ord, bin, chan), subgrid) in self.subgrids.indexed_iter() {
                let order = &self.orders[ord];

                if ((order.logxir != 0) && approx_eq!(f64, xir, 1.0, ulps = 4))
                    || ((order.logxif != 0) && approx_eq!(f64, xif, 1.0, ulps = 4))
                    || ((order.logxia != 0) && approx_eq!(f64, xia, 1.0, ulps = 4))
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

                if order.logxir != 0 {
                    value *= (xir * xir).ln().powi(order.logxir.into());
                }

                if order.logxif != 0 {
                    value *= (xif * xif).ln().powi(order.logxif.into());
                }

                if order.logxia != 0 {
                    value *= (xia * xia).ln().powi(order.logxia.into());
                }

                bins[xi_index + xi.len() * bin_index] += value / normalizations[bin];
            }
        }

        bins
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
        if let Some(bin) = self.bwfl().fill_index(observable) {
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
    pub fn read(reader: impl Read) -> Result<Self> {
        let mut reader = BufReader::new(reader);
        let buffer = reader.fill_buf().map_err(|err| Error::Other(err.into()))?;
        let magic_bytes: [u8; 4] = buffer[0..4].try_into().unwrap_or_else(|_| unreachable!());

        if u32::from_le_bytes(magic_bytes) == 0x18_4D_22_04 {
            Self::read_uncompressed(FrameDecoder::new(reader))
        } else {
            Self::read_uncompressed(reader)
        }
    }

    fn read_uncompressed(mut reader: impl BufRead) -> Result<Self> {
        let magic_bytes: [u8; 16] = reader.fill_buf().map_err(|err| Error::Other(err.into()))?
            [0..16]
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
            1 => bincode::deserialize_from(reader).map_err(|err| Error::Other(err.into())),
            _ => Err(Error::General(format!(
                "file version {file_version} is not supported"
            ))),
        }
    }

    /// Serializes `self` into `writer`. Writing is buffered.
    ///
    /// # Errors
    ///
    /// If writing fails an error is returned.
    pub fn write(&self, writer: impl Write) -> Result<()> {
        let mut writer = BufWriter::new(writer);
        let file_header = b"PineAPPL\x01\0\0\0\0\0\0\0";

        // first write PineAPPL file header
        writer
            .write(file_header)
            .map_err(|err| Error::Other(err.into()))?;

        // then serialize
        bincode::serialize_into(writer, self).map_err(|err| Error::Other(err.into()))
    }

    /// Serializes `self` into `writer`, using LZ4 compression. Writing is buffered.
    ///
    /// # Errors
    ///
    /// If writing or compression fails an error is returned.
    pub fn write_lz4(&self, writer: impl Write) -> Result<()> {
        let mut encoder = FrameEncoder::new(writer);
        self.write(&mut encoder)?;
        encoder
            .try_finish()
            .map_err(|err| Error::Other(err.into()))?;

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

    /// Merge the bins in indices in `range` together in a single one.
    ///
    /// # Errors
    ///
    /// When the given bins are non-consecutive, an error is returned.
    pub fn merge_bins(&mut self, range: impl RangeBounds<usize>) -> Result<()> {
        let range_start = match range.start_bound().cloned() {
            Bound::Included(start) => start,
            Bound::Excluded(start) => start + 1,
            Bound::Unbounded => 0,
        };
        let range_end = match range.end_bound().cloned() {
            Bound::Included(end) => end + 1,
            Bound::Excluded(end) => end,
            Bound::Unbounded => self.bwfl().len(),
        };

        // check if the bins in `range` can be merged - if not return without changing `self`
        self.bwfl = self
            .bwfl()
            .merge(range_start..range_end)
            // TODO: use proper error handling
            .unwrap_or_else(|_| unreachable!());

        let (intermediate, right) = self.subgrids.view().split_at(BIN_AXIS, range_end);
        let (left, merge) = intermediate.split_at(BIN_AXIS, range_start);

        let mut merged: Array2<SubgridEnum> = Array2::from_elem(
            (self.orders().len(), self.channels().len()),
            EmptySubgridV1.into(),
        );

        // merge the corresponding subgrids
        for subview in merge.axis_iter(BIN_AXIS) {
            Zip::from(&mut merged)
                .and(subview)
                .for_each(|lhs, rhs| lhs.merge(rhs, None));
        }
        let merged = merged.insert_axis(BIN_AXIS);

        self.subgrids = ndarray::concatenate(BIN_AXIS, &[left, merged.view(), right])
            // UNWRAP: if this fails there's a bug
            .unwrap_or_else(|_| unreachable!());

        Ok(())
    }

    /// Merge non-empty `Subgrid`s contained in `other` into `self`. Subgrids with the same bin
    /// limits are summed and subgrids with non-overlapping bin limits create new bins in `self`.
    ///
    /// # Errors
    ///
    /// If `self` and `other` in have different convolutions, PID bases, kinematics,
    /// interpolations, or scales an error is returned. If the bin limits of `self` and `other`
    /// are different and if the bin limits of `other` cannot be merged with `self` an error is
    /// returned.
    pub fn merge(&mut self, mut other: Self) -> Result<()> {
        if self.convolutions() != other.convolutions() {
            return Err(Error::General("convolutions do not match".to_owned()));
        }
        if self.pid_basis() != other.pid_basis() {
            return Err(Error::General("PID bases do not match".to_owned()));
        }
        // TODO: relax check if kinematic variables are permutations of each other
        if self.kinematics() != other.kinematics() {
            return Err(Error::General("kinematics do not match".to_owned()));
        }
        // TODO: relax check if subgrid types don't use interpolation
        if self.interpolations() != other.interpolations() {
            return Err(Error::General("interpolations do not match".to_owned()));
        }
        if self.scales() != other.scales() {
            return Err(Error::General("scales do not match".to_owned()));
        }

        let mut new_orders: Vec<Order> = Vec::new();
        let mut new_bins = 0;
        let mut new_entries: Vec<Channel> = Vec::new();

        if !self.bwfl().bins_partial_eq_with_ulps(other.bwfl(), 8) {
            new_bins = other.bwfl().len();

            // TODO: the following just appends bins to self, make this more general
            let lhs_r = self
                .bwfl()
                .fill_limits()
                .last()
                .copied()
                // UNWRAP: `BinsWithFillLimits` should guarantee there's always at least one bin
                .unwrap_or_else(|| unreachable!());
            let rhs_l = other
                .bwfl()
                .fill_limits()
                .first()
                .copied()
                // UNWRAP: `BinsWithFillLimits` should guarantee there's always at least one bin
                .unwrap_or_else(|| unreachable!());
            let new_bwfl = BinsWithFillLimits::new(
                [self.bwfl().bins(), other.bwfl().bins()].concat(),
                self.bwfl()
                    .fill_limits()
                    .iter()
                    .copied()
                    .chain(
                        other
                            .bwfl()
                            .fill_limits()
                            .iter()
                            .skip(1)
                            .map(|&limit| limit + lhs_r - rhs_l),
                    )
                    .collect(),
            )
            // TODO: do proper error handling
            .unwrap_or_else(|_| unreachable!());
            self.bwfl = new_bwfl;
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
            let old_dim = self.subgrids.raw_dim().into_pattern();
            let mut new_subgrids = Array3::from_shape_simple_fn(
                (
                    old_dim.0 + new_orders.len(),
                    old_dim.1 + new_bins,
                    old_dim.2 + new_entries.len(),
                ),
                || EmptySubgridV1.into(),
            );

            for (index, subgrid) in self.subgrids.indexed_iter_mut() {
                mem::swap(&mut new_subgrids[<[usize; 3]>::from(index)], subgrid);
            }

            self.subgrids = new_subgrids;
        }

        self.orders.append(&mut new_orders);
        self.channels.append(&mut new_entries);

        let bin_indices: Vec<_> = other
            .bwfl()
            .bins()
            .iter()
            .map(|bin| {
                self.bwfl()
                    .bins()
                    .iter()
                    .position(|other_bin| bin.partial_eq_with_ulps(other_bin, 8))
                    // UNWRAP: we've inserted the bins above so we must find them
                    .unwrap_or_else(|| unreachable!())
            })
            .collect();

        for ((i, j, k), subgrid) in other
            .subgrids
            .indexed_iter_mut()
            .filter(|((_, _, _), subgrid)| !subgrid.is_empty())
        {
            let other_order = &other.orders[i];
            let other_entry = &other.channels[k];

            let self_i = self
                .orders
                .iter()
                .position(|x| x == other_order)
                // UNWRAP: we added the orders previously so we must find it
                .unwrap_or_else(|| unreachable!());
            let self_j = bin_indices[j];
            let self_k = self
                .channels
                .iter()
                .position(|y| y == other_entry)
                // UNWRAP: we added the channels previously so we must find it
                .unwrap_or_else(|| unreachable!());

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

    /// TODO
    ///
    /// # Errors
    ///
    /// TODO
    pub fn set_bwfl(&mut self, bwfl: BinsWithFillLimits) -> Result<()> {
        let bins = bwfl.len();
        let grid_bins = self.bwfl().len();

        if bins != grid_bins {
            return Err(Error::General(format!(
                "{bins} are given, but the grid has {grid_bins}"
            )));
        }

        self.bwfl = bwfl;

        Ok(())
    }

    /// TODO
    #[must_use]
    pub const fn bwfl(&self) -> &BinsWithFillLimits {
        &self.bwfl
    }

    /// Calls [`Self::optimize_using`] with all possible optimization options
    /// ([`GridOptFlags::all`]).
    pub fn optimize(&mut self) {
        self.optimize_using(GridOptFlags::all());
    }

    /// Optimizes the internal datastructures for space efficiency. The parameter `flags`
    /// determines which optimizations are applied, see [`GridOptFlags`].
    pub fn optimize_using(&mut self, flags: GridOptFlags) {
        if flags.contains(GridOptFlags::OPTIMIZE_NODES) {
            self.optimize_nodes();
        }
        if flags.contains(GridOptFlags::OPTIMIZE_SUBGRID_TYPE) {
            self.optimize_subgrid_type();
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

    fn optimize_nodes(&mut self) {
        for subgrid in &mut self.subgrids {
            subgrid.optimize_nodes();
        }
    }

    fn optimize_subgrid_type(&mut self) {
        for subgrid in &mut self.subgrids {
            match subgrid {
                // replace empty subgrids of any type with `EmptySubgridV1`
                _ if subgrid.is_empty() => {
                    *subgrid = EmptySubgridV1.into();
                }
                _ => {
                    // TODO: check if we should remove this
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
        let mut ren1 = Vec::new();
        let mut fac1 = Vec::new();
        let mut frg1 = Vec::new();
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
            ren1.extend(
                self.scales()
                    .ren
                    .calc(&subgrid.node_values(), self.kinematics())
                    .iter(),
            );
            ren1.sort_by(f64::total_cmp);
            ren1.dedup_by(subgrid::node_value_eq_ref_mut);

            fac1.extend(
                self.scales()
                    .fac
                    .calc(&subgrid.node_values(), self.kinematics())
                    .iter(),
            );
            fac1.sort_by(f64::total_cmp);
            fac1.dedup_by(subgrid::node_value_eq_ref_mut);

            frg1.extend(
                self.scales()
                    .frg
                    .calc(&subgrid.node_values(), self.kinematics())
                    .iter(),
            );
            frg1.sort_by(f64::total_cmp);
            frg1.dedup_by(subgrid::node_value_eq_ref_mut);

            x1.extend(
                subgrid
                    .node_values()
                    .iter()
                    .zip(self.kinematics())
                    .filter(|(_, kin)| matches!(kin, Kinematics::X(_)))
                    .flat_map(|(nv, _)| nv),
            );

            x1.sort_by(f64::total_cmp);
            x1.dedup_by(subgrid::node_value_eq_ref_mut);

            for (index, _) in self.convolutions().iter().enumerate() {
                pids1.extend(channel.entry().iter().map(|(pids, _)| pids[index]));
            }

            pids1.sort_unstable();
            pids1.dedup();
        }

        EvolveInfo {
            fac1,
            frg1,
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
        S: IntoIterator<Item = std::result::Result<(OperatorSliceInfo, CowArray<'a, f64, Ix4>), E>>,
    >(
        &self,
        slices: Vec<S>,
        order_mask: &[bool],
        xi: (f64, f64, f64),
        alphas_table: &AlphasTable,
    ) -> Result<FkTable> {
        struct Iter<T> {
            iters: Vec<T>,
        }

        impl<T: Iterator> Iterator for Iter<T> {
            type Item = Vec<T::Item>;

            fn next(&mut self) -> Option<Self::Item> {
                self.iters.iter_mut().map(Iterator::next).collect()
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
        let mut fac0 = -1.0;
        let mut frg0 = -1.0;

        let mut perm = Vec::new();

        for result in zip_n(slices) {
            let (infos, operators): (Vec<OperatorSliceInfo>, Vec<CowArray<'_, f64, Ix4>>) = result
                .into_iter()
                .map(|res| res.map_err(|err| Error::Other(err.into())))
                .collect::<Result<_>>()?;

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

            if info_0.conv_type.is_pdf() {
                if fac0 < 0.0 {
                    fac0 = info_0.fac0;
                } else {
                    assert_approx_eq!(f64, fac0, info_0.fac0, ulps = 8);
                }
            } else if frg0 < 0.0 {
                frg0 = info_0.fac0;
            } else {
                assert_approx_eq!(f64, frg0, info_0.fac0, ulps = 8);
            }

            for (index, info) in infos_rest.iter().enumerate() {
                // TODO: what if the scales of the EKOs don't agree? Is there an ordering problem?
                assert!(subgrid::node_value_eq(info_0.fac1, info.fac1));

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

                if info.conv_type.is_pdf() {
                    if fac0 < 0.0 {
                        fac0 = info.fac0;
                    } else {
                        assert_approx_eq!(f64, fac0, info.fac0, ulps = 8);
                    }
                } else if frg0 < 0.0 {
                    frg0 = info.fac0;
                } else {
                    assert_approx_eq!(f64, frg0, info.fac0, ulps = 8);
                }
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
                .any(|&fac| subgrid::node_value_eq(fac, info_0.fac1))
            {
                continue;
            }

            // skip slices that the grid doesn't use
            if !grid_fac1
                .iter()
                .any(|&fac| subgrid::node_value_eq(fac, info_0.fac1))
            {
                continue;
            }

            let operators: Vec<_> = perm.iter().map(|&idx| operators[idx].view()).collect();
            let infos: Vec<_> = perm.iter().map(|&idx| infos[idx].clone()).collect();

            let fac = if fac0 < 0.0 {
                ScaleFuncForm::NoScale
            } else {
                ScaleFuncForm::Scale(0)
            };
            let (frg, scale_values) = if frg0 < 0.0 {
                (ScaleFuncForm::NoScale, vec![fac0])
            } else if fac0 < 0.0 || approx_eq!(f64, fac0, frg0, ulps = 8) {
                (ScaleFuncForm::Scale(0), vec![frg0])
            } else {
                (ScaleFuncForm::Scale(1), vec![fac0, frg0])
            };

            let (subgrids, channels) = evolution::evolve_slice(
                self,
                &operators,
                &infos,
                &scale_values,
                order_mask,
                xi,
                alphas_table,
            )?;

            let rhs = Self {
                subgrids,
                bwfl: self.bwfl().clone(),
                orders: vec![Order::new(0, 0, 0, 0, 0)],
                channels,
                pid_basis: infos[0].pid_basis,
                convolutions: self.convolutions.clone(),
                // TODO: the next line is probably wrong for flexible-scale grids
                interps: self.interps.clone(),
                kinematics: iter::once(Kinematics::Scale(0))
                    .chain(if fac0 < 0.0 && frg0 < 0.0 {
                        Some(Kinematics::Scale(1))
                    } else {
                        None
                    })
                    .chain(
                        self.kinematics
                            .iter()
                            .filter(|kin| matches!(kin, Kinematics::X(_)))
                            .copied(),
                    )
                    .collect(),
                scales: Scales {
                    // FK-tables have their renormalization scales burnt in
                    ren: ScaleFuncForm::NoScale,
                    fac,
                    frg,
                },
                metadata: self.metadata.clone(),
                more_members: self.more_members.clone(),
                // TODO: transform the reference result to match the FKTable structure
                reference: self.reference.clone(),
            };

            if let Some(lhs) = &mut lhs {
                lhs.merge(rhs)?;
            } else {
                lhs = Some(rhs);
            }

            used_op_fac1.push(info_0.fac1);
        }

        op_fac1.sort_by(f64::total_cmp);

        // make sure we've evolved all slices
        if let Some(muf2) = grid_fac1.into_iter().find(|&grid_mu2| {
            !used_op_fac1
                .iter()
                .any(|&eko_mu2| subgrid::node_value_eq(grid_mu2, eko_mu2))
        }) {
            return Err(Error::General(format!(
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
            .filter(|&index| index < self.bwfl().len())
            .collect();

        // sort and remove repeated indices
        bin_indices.sort_unstable();
        bin_indices.dedup();
        let bin_indices = bin_indices;

        for &bin_index in bin_indices.iter().rev() {
            self.subgrids.remove_index(Axis(1), bin_index);
            self.bwfl.remove(bin_index);
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
    use crate::interpolation::Map;
    use float_cmp::assert_approx_eq;
    use std::fs::File;

    #[test]
    fn interpolations() {
        let grid = Grid::new(
            BinsWithFillLimits::from_fill_limits([0.0, 1.0].to_vec()).unwrap(),
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![Channel::new(vec![(vec![1, -1], 1.0), (vec![2, -2], 1.0)])],
            PidBasis::Pdg,
            vec![
                Conv::new(ConvType::UnpolPDF, 2212),
                Conv::new(ConvType::UnpolPDF, 2212),
            ],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        let interps = grid.interpolations();
        assert!(matches!(interps[0].map(), Map::ApplGridH0));
        assert!(matches!(interps[1].map(), Map::ApplGridF2));
        assert!(matches!(interps[2].map(), Map::ApplGridF2));
    }

    #[test]
    #[should_panic(expected = "channel #0 has wrong number of PIDs: expected 2, found 3")]
    fn grid_new_panic0() {
        let channel = vec![(vec![1, -1, 1], 1.0), (vec![2, -2, 2], 1.0)];

        let _ = Grid::new(
            BinsWithFillLimits::from_fill_limits([0.0, 1.0].to_vec()).unwrap(),
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![Channel::new(channel)],
            PidBasis::Pdg,
            vec![
                Conv::new(ConvType::UnpolPDF, 2212),
                Conv::new(ConvType::UnpolPDF, 2212),
            ],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
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
            BinsWithFillLimits::from_fill_limits([0.0, 1.0].to_vec()).unwrap(),
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![Channel::new(channel)],
            PidBasis::Pdg,
            vec![
                Conv::new(ConvType::UnpolPDF, 2212),
                Conv::new(ConvType::UnpolPDF, 2212),
            ],
            v0::default_interps(false, 1),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
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
            BinsWithFillLimits::from_fill_limits([0.0, 1.0].to_vec()).unwrap(),
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![Channel::new(channel)],
            PidBasis::Pdg,
            vec![
                Conv::new(ConvType::UnpolPDF, 2212),
                Conv::new(ConvType::UnpolPDF, 2212),
            ],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(1),
                frg: ScaleFuncForm::NoScale,
            },
        );
    }

    #[test]
    fn grid_read_file_version_unsupported() {
        let result = Grid::read(
            &[
                b'P', b'i', b'n', b'e', b'A', b'P', b'P', b'L', 99, 0, 0, 0, 0, 0, 0, 0,
            ][..],
        );

        assert!(
            matches!(result, Err(Error::General(msg)) if msg == "file version 99 is not supported")
        );
    }

    #[test]
    fn grid_merge_empty_subgrids() {
        let mut grid = Grid::new(
            BinsWithFillLimits::from_fill_limits([0.0, 0.25, 0.5, 0.75, 1.0].to_vec()).unwrap(),
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![
                channel![1.0 * (2, 2) + 1.0 * (4, 4)],
                channel![1.0 * (1, 1) + 1.0 * (3, 3)],
            ],
            PidBasis::Pdg,
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        assert_eq!(grid.bwfl().len(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let other = Grid::new(
            BinsWithFillLimits::from_fill_limits([0.0, 0.25, 0.5, 0.75, 1.0].to_vec()).unwrap(),
            vec![Order::new(1, 2, 0, 0, 0), Order::new(1, 2, 0, 1, 0)],
            vec![
                // differently ordered than `grid`
                channel![1.0 * (1, 1) + 1.0 * (3, 3)],
                channel![1.0 * (2, 2) + 1.0 * (4, 4)],
            ],
            PidBasis::Pdg,
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        // merging with empty subgrids should not change the grid
        grid.merge(other).unwrap();

        assert_eq!(grid.bwfl().len(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);
    }

    #[test]
    fn grid_merge_orders() {
        let mut grid = Grid::new(
            BinsWithFillLimits::from_fill_limits([0.0, 0.25, 0.5, 0.75, 1.0].to_vec()).unwrap(),
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![
                channel![1.0 * (2, 2) + 1.0 * (4, 4)],
                channel![1.0 * (1, 1) + 1.0 * (3, 3)],
            ],
            PidBasis::Pdg,
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        assert_eq!(grid.bwfl().len(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            BinsWithFillLimits::from_fill_limits([0.0, 0.25, 0.5, 0.75, 1.0].to_vec()).unwrap(),
            vec![
                Order::new(1, 2, 0, 0, 0),
                Order::new(1, 2, 0, 1, 0),
                Order::new(0, 2, 0, 0, 0),
            ],
            vec![
                channel![1.0 * (2, 2) + 1.0 * (4, 4)],
                channel![1.0 * (1, 1) + 1.0 * (3, 3)],
            ],
            PidBasis::Pdg,
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
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

        assert_eq!(grid.bwfl().len(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 3);
    }

    #[test]
    fn grid_merge_channels_entries() {
        let mut grid = Grid::new(
            BinsWithFillLimits::from_fill_limits([0.0, 0.25, 0.5, 0.75, 1.0].to_vec()).unwrap(),
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![
                channel![1.0 * (2, 2) + 1.0 * (4, 4)],
                channel![1.0 * (1, 1) + 1.0 * (3, 3)],
            ],
            PidBasis::Pdg,
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        assert_eq!(grid.bwfl().len(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            BinsWithFillLimits::from_fill_limits([0.0, 0.25, 0.5, 0.75, 1.0].to_vec()).unwrap(),
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![
                channel![1.0 * (22, 22)],
                channel![1.0 * (2, 2) + 1.0 * (4, 4)],
            ],
            PidBasis::Pdg,
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        // fill the photon-photon entry
        other.fill(0, 0.1, 0, &[90.0_f64.powi(2), 0.1, 0.2], 3.0);

        grid.merge(other).unwrap();

        assert_eq!(grid.bwfl().len(), 4);
        assert_eq!(grid.channels().len(), 3);
        assert_eq!(grid.orders().len(), 1);
    }

    #[test]
    fn grid_merge_bins() {
        let mut grid = Grid::new(
            BinsWithFillLimits::from_fill_limits([0.0, 0.25, 0.5].to_vec()).unwrap(),
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![
                channel![1.0 * (2, 2) + 1.0 * (4, 4)],
                channel![1.0 * (1, 1) + 1.0 * (3, 3)],
            ],
            PidBasis::Pdg,
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        assert_eq!(grid.bwfl().len(), 2);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            BinsWithFillLimits::from_fill_limits([0.5, 0.75, 1.0].to_vec()).unwrap(),
            vec![Order::new(0, 2, 0, 0, 0)],
            vec![
                // channels are differently sorted
                channel![1.0 * (1, 1) + 1.0 * (3, 3)],
                channel![1.0 * (2, 2) + 1.0 * (4, 4)],
            ],
            PidBasis::Pdg,
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        other.fill(0, 0.1, 0, &[90.0_f64.powi(2), 0.1, 0.2], 2.0);
        other.fill(0, 0.1, 1, &[90.0_f64.powi(2), 0.1, 0.2], 3.0);

        grid.merge(other).unwrap();

        assert_eq!(grid.bwfl().len(), 4);
        assert_eq!(grid.channels().len(), 2);
        assert_eq!(grid.orders().len(), 1);
    }

    #[test]
    fn grid_convolutions() {
        let mut grid = Grid::new(
            BinsWithFillLimits::from_fill_limits([0.0, 1.0].to_vec()).unwrap(),
            vec![Order::new(0, 0, 0, 0, 0)],
            vec![channel![1.0 * (21, 21)]],
            PidBasis::Pdg,
            vec![Conv::new(ConvType::UnpolPDF, 2212); 2],
            v0::default_interps(false, 2),
            vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
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
