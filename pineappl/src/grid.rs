//! Module containing all traits and supporting structures for grids.

use super::bin::BinLimits;
use super::lumi::LumiEntry;
use super::ntuple_grid::NtupleSubgrid;
use ndarray::Array3;
use serde::{Deserialize, Serialize};
use std::ops::MulAssign;

/// Coupling powers for each grid.
#[derive(Clone, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub struct Order {
    /// Exponent of the strong coupling.
    pub alphas: u32,
    /// Exponent of the electromagnetic coupling.
    pub alpha: u32,
    /// Exponent of the logarithm of the scale factor of the renomalization scale.
    pub logxir: u32,
    /// Exponent of the logarithm of the scale factor of the factorization scale.
    pub logxif: u32,
}

impl Order {
    /// Compares two vectors of `Order` for equality after sorting them.
    pub fn equal_after_sort(lhs: &Vec<Order>, rhs: &Vec<Order>) -> bool {
        let mut lhs = lhs.clone();
        let mut rhs = rhs.clone();

        lhs.sort();
        rhs.sort();

        lhs == rhs
    }
}

/// Trait each subgrid must implement.
#[typetag::serde(tag = "type")]
pub trait Subgrid {
    /// Convolute the subgrid with a luminosity function
    fn convolute(&self, lumi: &dyn Fn(f64, f64, f64) -> f64) -> f64;

    /// Fills the subgrid with `weight` for the parton momentum fractions `x1` and `x2`, and the
    /// scale `q2`.
    fn fill(&mut self, ntuple: Ntuple<f64>);

    /// Returns true if `fill` was never called for this grid.
    fn is_empty(&self) -> bool;

    /// Scale the subgrid by `factor`.
    fn scale(&mut self, factor: f64);
}

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

impl MulAssign<f64> for Ntuple<f64> {
    fn mul_assign(&mut self, rhs: f64) {
        self.weight *= rhs;
    }
}

/// Main data structure of `PineAPPL`. This structure contains a `Subgrid` for each `LumiEntry`,
/// bin, and coupling order it was created with.
#[derive(Deserialize, Serialize)]
pub struct Grid {
    subgrids: Array3<Box<dyn Subgrid>>,
    lumi: Vec<LumiEntry>,
    bin_limits: BinLimits,
    orders: Vec<Order>,
}

impl Grid {
    /// Constructor.
    #[must_use]
    pub fn new(lumi: Vec<LumiEntry>, orders: Vec<Order>, bin_limits: Vec<f64>) -> Self {
        Self {
            subgrids: Array3::from_shape_simple_fn(
                (orders.len(), bin_limits.len() - 1, lumi.len()),
                || Box::new(NtupleSubgrid::default()) as Box<dyn Subgrid>,
            ),
            orders,
            lumi,
            bin_limits: BinLimits::new(bin_limits),
        }
    }

    /// Returns the bin limits of the observables.
    #[must_use]
    pub const fn bin_limits(&self) -> &BinLimits {
        &self.bin_limits
    }

    /// Performs a convolution of the contained subgrids with the given PDFs, `xfx1` for the first
    /// parton and `xfx2` for the second parton, `alphas` for the evaluation of the strong
    /// coupling. The parameters `order_mask` and `lumi_mask` can be used to selectively enable
    /// perturbative orders and luminosities; they must either be empty (everything enabled) or as
    /// large as the orders and luminosity function, respectively. If the corresponding entry is
    /// `true` the order/luminosity is enable, `false` disables the entry. The tuple `xi` can be
    /// used to independently vary the renormalization (first element) and factorization scale
    /// (second element) from their central value `(1.0, 1.0)`.
    pub fn convolute(
        &self,
        xfx1: &dyn Fn(f64, f64, i32) -> f64,
        xfx2: &dyn Fn(f64, f64, i32) -> f64,
        alphas: &dyn Fn(f64) -> f64,
        order_mask: &[bool],
        lumi_mask: &[bool],
        xi: &(f64, f64),
    ) -> Vec<f64> {
        let mut bins: Vec<f64> = vec![0.0; self.bin_limits.bins()];

        for ((i, j, k), subgrid) in self.subgrids.indexed_iter() {
            let order = &self.orders[i];

            if (!order_mask.is_empty() && !order_mask[i])
                || ((order.logxir > 0) && (xi.0 == 1.0))
                || ((order.logxif > 0) && (xi.1 == 1.0))
                || (!lumi_mask.is_empty() && !lumi_mask[k])
            {
                continue;
            }

            let lumi_entry = &self.lumi[k];

            let mut value = subgrid.convolute(&|x1, x2, q2| {
                let mut lumi = 0.0;

                for entry in lumi_entry.entry() {
                    lumi += xfx1(x1, q2, entry.0) * xfx2(x2, q2, entry.1) * entry.2;
                }

                lumi *= alphas(q2).powi(order.alphas as i32);
                lumi
            });

            if order.logxir > 0 {
                value *= xi.0.ln().powi(order.logxir as i32);
            }

            if order.logxif > 0 {
                value *= xi.1.ln().powi(order.logxif as i32);
            }

            bins[j] += value;
        }

        bins
    }

    /// Fills the grid with an ntuple for the given `order`, `observable`, and `lumi`.
    pub fn fill(&mut self, order: usize, observable: f64, lumi: usize, ntuple: Ntuple<f64>) {
        if let Some(bin) = self.bin_limits.index(observable) {
            self.subgrids[[order, bin, lumi]].fill(ntuple);
        }
    }

    /// Fills the grid with events for the parton momentum fractions `x1` and `x2`, the scale `q2`,
    /// and the `order` and `observable`. The events are stored in `weights` and must be ordered as
    /// the corresponding luminosity function was created.
    pub fn fill_all(&mut self, order: usize, observable: f64, ntuple: Ntuple<()>, weights: &[f64]) {
        if let Some(bin) = self.bin_limits.index(observable) {
            for (lumi, weight) in weights.iter().enumerate() {
                self.subgrids[[order, bin, lumi]].fill(Ntuple {
                    x1: ntuple.x1,
                    x2: ntuple.x2,
                    q2: ntuple.q2,
                    weight: *weight,
                });
            }
        }
    }

    /// Returns the luminosity function.
    #[must_use]
    pub fn lumi(&self) -> &[LumiEntry] {
        &self.lumi
    }

    /// Scale all subgrids by `factor`.
    pub fn scale(&mut self, factor: f64) {
        self.subgrids
            .iter_mut()
            .for_each(|subgrid| subgrid.scale(factor));
    }

    /// Returns the subgrid parameters.
    #[must_use]
    pub fn orders(&self) -> &[Order] {
        &self.orders
    }
}
