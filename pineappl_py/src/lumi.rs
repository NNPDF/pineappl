use pineappl::boc::Channel;

use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::lumi::LumiEntry <lumi/struct.LumiEntry.html>`
///
/// **Usage**: `yadism`, FKTable interface
///
/// Each entry consists of a tuple, which contains, in the following order:
///
/// 1. the PDG id of the first incoming parton
/// 2. the PDG id of the second parton
/// 3. a numerical factor that will multiply the result for this specific combination.
#[pyclass]
#[repr(transparent)]
pub struct PyLumiEntry {
    pub(crate) lumi_entry: Channel,
}

impl PyLumiEntry {
    pub(crate) fn new(lumi_entry: Channel) -> Self {
        Self { lumi_entry }
    }
}

#[pymethods]
impl PyLumiEntry {
    #[new]
    pub fn new_lumi_entry(entry: Vec<(i32, i32, f64)>) -> Self {
        Self::new(Channel::new(entry))
    }

    /// Get list representation.
    ///
    /// **Usage:** FKTable interface
    ///
    /// Returns
    /// -------
    ///     list(tuple(int,int,float)) :
    ///         list representation
    pub fn into_array(&self) -> Vec<(i32, i32, f64)> {
        self.lumi_entry.entry().to_vec()
    }
}
