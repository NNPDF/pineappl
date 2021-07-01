use pineappl::import_only_subgrid::ImportOnlySubgridV1;
use pineappl::sparse_array3::SparseArray3;
use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyImportOnlySubgrid {
    pub import_only_subgrid: ImportOnlySubgridV1,
}

impl PyImportOnlySubgrid {
    pub(crate) fn new(import_only_subgrid: ImportOnlySubgridV1) -> Self {
        Self {
            import_only_subgrid,
        }
    }
}

#[pymethods]
impl PyImportOnlySubgrid {
    #[new]
    pub fn new_import_only_subgrid(
        q2_grid: Vec<f64>,
        x1_grid: Vec<f64>,
        x2_grid: Vec<f64>,
    ) -> Self {
        Self::new(ImportOnlySubgridV1::new(
            SparseArray3::new(q2_grid.len(), x1_grid.len(), x2_grid.len()),
            q2_grid,
            x1_grid,
            x2_grid,
        ))
    }
}
