use pineappl::lumi::LumiEntry;

use pyo3::prelude::*;

/// PyO3 wrapper to [`pineappl::lumi::LumiEntry`]
///
/// **Usage**: `yadism`, FKTable interface
#[pyclass]
#[repr(transparent)]
pub struct PyLumiEntry {
    pub(crate) lumi_entry: LumiEntry,
}

impl PyLumiEntry {
    pub(crate) fn new(lumi_entry: LumiEntry) -> Self {
        Self { lumi_entry }
    }
}

#[pymethods]
impl PyLumiEntry {
    #[new]
    pub fn new_lumi_entry(entry: Vec<(i32, i32, f64)>) -> Self {
        Self::new(LumiEntry::new(entry))
    }

    /// Clone elements to allow Python to access.
    ///
    /// **Usage:** FKTable interface
    pub fn into_array(&self) -> Vec<(i32, i32,f64)> {
        self.lumi_entry.entry().iter().map(|e| e.clone()).collect()
    }
}
