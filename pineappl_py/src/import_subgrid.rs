//! PyPackedSubgrid* interface.

use super::subgrid::PySubgridEnum;
use ndarray::Dimension;
use numpy::PyReadonlyArrayDyn;
use pineappl::import_subgrid::ImportSubgridV1;
use pineappl::packed_array::PackedArray;
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::import_subgrid::ImportSubgridV1 <import_subgrid/struct.ImportSubgridV1.html>`.
#[pyclass(name = "ImportSubgridV1")]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyImportSubgridV1 {
    pub(crate) import_subgrid: ImportSubgridV1,
}

#[pymethods]
impl PyImportSubgridV1 {
    /// Constructor.
    ///
    /// # Panics
    ///
    /// Parameters
    /// ----------
    /// array : numpy.ndarray(float)
    ///     `N`-dimensional array with all weights
    /// node_values: list(list(float))
    ///     list containing the arrays of energy scales {q1, ..., qn} and momentum fractions
    ///     {x1, ..., xn}.
    #[new]
    #[must_use]
    pub fn new(array: PyReadonlyArrayDyn<f64>, node_values: Vec<Vec<f64>>) -> Self {
        let mut sparse_array: PackedArray<f64> =
            PackedArray::new(node_values.iter().map(Vec::len).collect());

        for (index, value) in array
            .as_array()
            .indexed_iter()
            .filter(|(_, value)| **value != 0.0)
        {
            sparse_array[index.as_array_view().to_slice().unwrap()] = *value;
        }

        Self {
            import_subgrid: ImportSubgridV1::new(sparse_array, node_values),
        }
    }

    /// TODO
    #[must_use]
    pub fn into(&self) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.import_subgrid.clone().into(),
        }
    }
}

/// Register submodule in parent.
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "import_subgrid")?;
    m.setattr(
        pyo3::intern!(m.py(), "__doc__"),
        "Interface for packed subgrid specs.",
    )?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.import_subgrid'] = m"
    );
    m.add_class::<PyImportSubgridV1>()?;
    parent_module.add_submodule(&m)
}
