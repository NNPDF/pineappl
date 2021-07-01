use pineappl::import_only_subgrid::ImportOnlySubgridV1;
use pineappl::subgrid::Subgrid;

use super::subgrid::{PyExtraSubgridParams, PySubgridParams};

use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyImportOnlySubgrid {
    pub import_only_subgrid: ImportOnlySubgridV1,
}

impl PyImportOnlySubgrid {
    pub(crate) fn new(import_only_subgrid: ImportOnlySubgridV1) -> Self {
        Self { import_only_subgrid }
    }
}

#[pymethods]
impl PyImportOnlySubgrid {
    #[new]
    pub fn new_import_only_subgrid(
    ) -> Self {
        // Self::new(ImportOnlySubgridV1::new(
        //     &subgrid_params.subgrid_params,
        //     &extra_params.extra_subgrid_params,
        // ))
    }
}
