//! TODO

use serde::{Deserialize, Serialize};

#[derive(Clone, Deserialize, Serialize)]
enum RefRes {
    ByBin(f64),
    ByBinOrder(Vec<f64>),
    ByBinChannel(Vec<f64>),
    ByBinOrderChannel(Vec<Vec<f64>>),
}

#[derive(Clone, Deserialize, Serialize)]
enum RefRelUnc {
    ByBin(f64),
    ByBinOrder(Vec<f64>, CombOp),
    ByBinChannel(Vec<f64>, CombOp),
    ByBinOrderChannel(Vec<Vec<f64>>, CombOp),
}

#[derive(Clone, Deserialize, Serialize)]
enum CombOp {
    Sum,
    Quadrature,
}

/// TODO
#[derive(Clone, Default, Deserialize, Serialize)]
pub struct Reference {
    ref_res_unc: Vec<(RefRes, RefRelUnc)>,
    ref_conv_fun: Vec<String>,
}
