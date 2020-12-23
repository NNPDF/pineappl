use pineappl::lumi::LumiEntry;

use pyo3::prelude::*;

#[pyclass]
#[repr(transparent)]
pub struct PyLumiEntry {
    pub lumi_entry: LumiEntry,
}

impl PyLumiEntry {
    pub(crate) fn new(lumi_entry: LumiEntry) -> Self {
        Self { lumi_entry }
    }
}

#[pymethods]
impl PyLumiEntry {
    #[new]
    pub fn new_lumi_entry(mut entry: Vec<(i32, i32, f64)>) -> Self {
        Self::new(LumiEntry::new(entry))
    }
}
