//! Optional reference cross sections and uncertainties stored alongside a [`crate::grid::Grid`].

use serde::{Deserialize, Serialize};

/// Absolute reference cross section for a single bin.
#[derive(Clone, Deserialize, Serialize)]
pub enum AbsRefRes {
    /// One value per bin.
    ByBin(f64),
    /// One value per bin and perturbative order.
    ByBinOrder(Vec<f64>),
    /// One value per bin and luminosity channel.
    ByBinChannel(Vec<f64>),
    /// Values for each bin, order, and channel.
    ByBinOrderChannel(Vec<Vec<f64>>),
}

/// Relative reference uncertainty for a single bin.
#[derive(Clone, Deserialize, Serialize)]
pub enum RelRefUnc {
    /// Relative uncertainty per bin.
    ByBin(f64),
    /// Per-order uncertainties combined with `CombOp`.
    ByBinOrder(Vec<f64>, CombOp),
    /// Per-channel uncertainties combined with `CombOp`.
    ByBinChannel(Vec<f64>, CombOp),
    /// Per-order and channel uncertainties combined with `CombOp`.
    ByBinOrderChannel(Vec<Vec<f64>>, CombOp),
}

/// How to combine several relative uncertainties into one.
#[derive(Clone, Deserialize, Serialize)]
pub enum CombOp {
    /// Linear sum of relative terms.
    Sum,
    /// Square root of the sum of squares (in quadrature).
    Quadrature,
}

/// Reference results and convolution-function labels for validation or plotting.
#[derive(Clone, Default, Deserialize, Serialize)]
pub struct Reference {
    ref_res_unc: Vec<(AbsRefRes, RelRefUnc)>,
    ref_conv_fun: Vec<String>,
}
