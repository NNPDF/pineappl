use pineappl::boc::Channel;

use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::boc::Channel <boc/struct.Channel.html>`.
///
/// Each entry consists of a tuple, which contains, in the following order:
///
/// 1. the PDG id of the first incoming parton
/// 2. the PDG id of the second parton
/// 3. a numerical factor that will multiply the result for this specific combination.
#[pyclass(name = "Channel")]
#[repr(transparent)]
pub struct PyChannel {
    pub(crate) entry: Channel,
}

impl PyChannel {
    pub(crate) fn new(entry: Channel) -> Self {
        Self { entry }
    }
}

#[pymethods]
impl PyChannel {
    /// Constructor.
    ///
    /// Parameters
    /// ----------
    /// entry: list(tuple(int, int, float))
    ///     channel configuration
    #[new]
    pub fn new_entry(entry: Vec<(i32, i32, f64)>) -> Self {
        Self::new(Channel::new(entry))
    }

    /// Get list representation.
    ///
    /// Returns
    /// -------
    /// list(tuple(int,int,float)) :
    ///     list representation
    pub fn into_array(&self) -> Vec<(i32, i32, f64)> {
        self.entry.entry().to_vec()
    }
}

#[pymodule]
pub fn channel(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyChannel>()?;
    Ok(())
}
