//! TODO

use serde::{Deserialize, Serialize};

/// Absolute reference result for a single bin.
#[derive(Clone, Deserialize, Serialize)]
pub enum AbsRefRes {
    /// TODO
    ByBin(f64),
    /// TODO
    ByBinOrder(Vec<f64>),
    /// TODO
    ByBinChannel(Vec<f64>),
    /// TODO
    ByBinOrderChannel(Vec<Vec<f64>>),
}

/// Relative reference uncertainty for a single bin.
#[derive(Clone, Deserialize, Serialize)]
pub enum RelRefUnc {
    /// TODO
    ByBin(f64),
    /// TODO
    ByBinOrder(Vec<f64>, CombOp),
    /// TODO
    ByBinChannel(Vec<f64>, CombOp),
    /// TODO
    ByBinOrderChannel(Vec<Vec<f64>>, CombOp),
}

/// TODO
#[derive(Clone, Deserialize, Serialize)]
pub enum CombOp {
    /// TODO
    Sum,
    /// TODO
    Quadrature,
}

/// TODO
#[derive(Clone, Default, Deserialize, Serialize)]
pub struct Reference {
    ref_res_unc: Vec<(AbsRefRes, RelRefUnc)>,
    ref_conv_fun: Vec<String>,
}
