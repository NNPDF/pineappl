use ndarray::CowArray;
use pineappl::boc::Order;
use pineappl::convolutions::LumiCache;
use pineappl::evolution::{AlphasTable, OperatorInfo, OperatorSliceInfo};
use pineappl::grid::{Grid, Ntuple};
use pineappl::pids::PidBasis;

use super::bin::PyBinRemapper;
use super::evolution::PyEvolveInfo;
use super::fk_table::PyFkTable;
use super::lumi::PyLumiEntry;
use super::subgrid::{PySubgridEnum, PySubgridParams};

use itertools::izip;
use numpy::{
    IntoPyArray, PyArray1, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray4, PyReadonlyArray5,
};

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyIterator;

/// PyO3 wrapper to :rustdoc:`pineappl::grid::Order <grid/struct.Order.html>`
///
/// **Usage**: `yadism`
#[pyclass]
#[repr(transparent)]
pub struct PyOrder {
    pub(crate) order: Order,
}

// TODO: should probably be in a different module
// TODO: rename to `PidBasis`
#[pyclass]
#[derive(Clone)]
pub enum PyPidBasis {
    Pdg,
    Evol,
}

impl From<PyPidBasis> for PidBasis {
    fn from(basis: PyPidBasis) -> Self {
        match basis {
            PyPidBasis::Pdg => Self::Pdg,
            PyPidBasis::Evol => Self::Evol,
        }
    }
}

// TODO: should probably be in a different module
// TODO: rename to `OperatorSliceInfo`
#[pyclass]
#[derive(Clone)]
pub struct PyOperatorSliceInfo {
    info: OperatorSliceInfo,
}

#[pymethods]
impl PyOperatorSliceInfo {
    #[new]
    pub fn new(
        fac0: f64,
        pids0: Vec<i32>,
        x0: Vec<f64>,
        fac1: f64,
        pids1: Vec<i32>,
        x1: Vec<f64>,
        pid_basis: PyPidBasis,
    ) -> Self {
        Self {
            info: OperatorSliceInfo {
                fac0,
                pids0,
                x0,
                fac1,
                pids1,
                x1,
                pid_basis: pid_basis.into(),
            },
        }
    }
}

impl PyOrder {
    pub(crate) fn new(order: Order) -> Self {
        Self { order }
    }
}

#[pymethods]
impl PyOrder {
    #[new]
    pub fn new_order(alphas: u32, alpha: u32, logxir: u32, logxif: u32) -> Self {
        Self::new(Order::new(alphas, alpha, logxir, logxif))
    }

    /// Tuple representation.
    ///
    /// Returns
    /// -------
    ///     alphas : int
    ///         power of :math:`\alpha_s`
    ///     alpha : int
    ///         power of :math:`\alpha`
    ///     logxir : int
    ///         power of :math:` \ln(\xi_r)`
    ///     logxif : int
    ///         power of :math:` \ln(\xi_f)`
    pub fn as_tuple(&self) -> (u32, u32, u32, u32) {
        (
            self.order.alphas,
            self.order.alpha,
            self.order.logxir,
            self.order.logxif,
        )
    }
    /// Return a mask suitable to pass as the `order_mask` parameter of [`Grid::convolve`]. The
    /// selection of `orders` is controlled using the `max_as` and `max_al` parameters, for
    /// instance setting `max_as = 1` and `max_al = 0` selects the LO QCD only, `max_as = 2` and
    /// `max_al = 0` the NLO QCD; setting `max_as = 3` and `max_al = 2` would select all NLOs, and
    /// the NNLO QCD.
    ///
    /// See `pineappl` crate docs for relevant examples
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(bool)
    ///     boolean array, to be used as orders' mask
    #[staticmethod]
    pub fn create_mask<'py>(
        orders: Vec<PyRef<Self>>,
        max_as: u32,
        max_al: u32,
        logs: bool,
        py: Python<'py>,
    ) -> Bound<'py, PyArray1<bool>> {
        Order::create_mask(
            &orders.iter().map(|o| o.order.clone()).collect::<Vec<_>>(),
            max_as,
            max_al,
            logs,
        )
        .into_pyarray_bound(py)
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::grid::Grid <grid/struct.Grid.html>`
///
/// **Usage**: `yadism`, `pineko`, FKTable interface
#[pyclass]
#[repr(transparent)]
#[derive(Clone)]
pub struct PyGrid {
    pub(crate) grid: Grid,
}

impl PyGrid {
    pub(crate) fn new(grid: Grid) -> Self {
        Self { grid }
    }
}

#[pymethods]
impl PyGrid {
    #[new]
    pub fn new_grid(
        lumi: Vec<PyRef<PyLumiEntry>>,
        orders: Vec<PyRef<PyOrder>>,
        bin_limits: PyReadonlyArray1<f64>,
        subgrid_params: PySubgridParams,
    ) -> Self {
        Self::new(Grid::new(
            lumi.iter().map(|pyl| pyl.lumi_entry.clone()).collect(),
            orders.iter().map(|pyo| pyo.order.clone()).collect(),
            bin_limits.to_vec().unwrap(),
            subgrid_params.subgrid_params,
        ))
    }

    /// Add a point to the grid.
    ///
    /// Parameters
    /// ----------
    ///     x1 : float
    ///         first momentum fraction
    ///     x2 : float
    ///         second momentum fraction
    ///     q2 : float
    ///         process scale
    ///     order : int
    ///         order index
    ///     observable : float
    ///         reference point (to be binned)
    ///     lumi : int
    ///         luminosity index
    ///     weight : float
    ///         cross section weight
    pub fn fill(
        &mut self,
        x1: f64,
        x2: f64,
        q2: f64,
        order: usize,
        observable: f64,
        lumi: usize,
        weight: f64,
    ) {
        self.grid.fill(
            order,
            observable,
            lumi,
            &Ntuple::<f64> { x1, x2, q2, weight },
        );
    }

    /// Add an array to the grid.
    ///
    /// Useful to avoid multiple python calls, leading to performance improvement.
    ///
    /// Parameters
    /// ----------
    ///     ntuples : np.array(float)
    ///         2 dimensional (4, N) array, made of `(x1, x2, q2, weight)` "ntuples"
    ///     order : int
    ///         order index
    ///     observable : float
    ///         reference point (to be binned)
    ///     lumi : int
    ///         luminosity index
    pub fn fill_array(
        &mut self,
        x1s: PyReadonlyArray1<f64>,
        x2s: PyReadonlyArray1<f64>,
        q2s: PyReadonlyArray1<f64>,
        order: usize,
        observables: PyReadonlyArray1<f64>,
        lumi: usize,
        weights: PyReadonlyArray1<f64>,
    ) {
        for (&x1, &x2, &q2, &observable, &weight) in izip!(
            x1s.as_array().iter(),
            x2s.as_array().iter(),
            q2s.as_array().iter(),
            observables.as_array().iter(),
            weights.as_array().iter(),
        ) {
            self.grid.fill(
                order,
                observable,
                lumi,
                &Ntuple::<f64> { x1, x2, q2, weight },
            );
        }
    }

    /// Add a point to the grid for all lumis.
    ///
    /// Parameters
    /// ----------
    ///     x1 : float
    ///         first momentum fraction
    ///     x2 : float
    ///         second momentum fraction
    ///     q2 : float
    ///         process scale
    ///     order : int
    ///         order index
    ///     observable : float
    ///         reference point (to be binned)
    ///     weights : np.array(float)
    ///         cross section weights, one for each lumi
    pub fn fill_all(
        &mut self,
        x1: f64,
        x2: f64,
        q2: f64,
        order: usize,
        observable: f64,
        weights: PyReadonlyArray1<f64>,
    ) {
        self.grid.fill_all(
            order,
            observable,
            &Ntuple::<()> {
                x1,
                x2,
                q2,
                weight: (),
            },
            &weights.to_vec().unwrap(),
        );
    }

    /// Get metadata values stored in the grid.
    ///
    ///
    /// Returns
    /// -------
    ///     dict :
    ///         key, value map
    pub fn key_values(&self) -> HashMap<String, String> {
        self.grid.key_values().unwrap().clone()
    }

    /// Set a metadata key-value pair in the grid.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    ///     key : str
    ///         key
    ///     value : str
    ///         value
    pub fn set_key_value(&mut self, key: &str, value: &str) {
        self.grid.set_key_value(key, value);
    }

    /// Retrieve a subgrid.
    ///
    /// **Usage:** `yadism`
    pub fn subgrid(&self, order: usize, bin: usize, lumi: usize) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.grid.subgrids()[[order, bin, lumi]].clone(),
        }
    }

    /// Set a subgrid.
    ///
    /// **Usage:** `yadism`
    pub fn set_subgrid(&mut self, order: usize, bin: usize, lumi: usize, subgrid: PySubgridEnum) {
        self.grid.subgrids_mut()[[order, bin, lumi]] = subgrid.subgrid_enum;
    }

    /// Set the normalizations.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    ///     remapper: BinRemapper
    ///         Remapper object
    pub fn set_remapper(&mut self, remapper: PyBinRemapper) {
        self.grid.set_remapper(remapper.bin_remapper).unwrap();
    }

    /// Convolute grid with pdf.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     pdg_id : int
    ///         PDG Monte Carlo ID of the hadronic particle `xfx` is the PDF for
    ///     xfx : callable
    ///         lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    ///     alphas : callable
    ///         lhapdf like callable with arguments `Q2` returning :math:`\alpha_s`
    ///     order_mask : numpy.ndarray(bool)
    ///         Mask for selecting specific orders. The value `True` means the corresponding order
    ///         is included. An empty list corresponds to all orders being enabled.
    ///     bin_indices : numpy.ndarray(int)
    ///         A list with the indices of the corresponding bins that should be calculated. An
    ///         empty list means that all orders should be calculated.
    ///     lumi_mask : numpy.ndarray(bool)
    ///         Mask for selecting specific luminosity channels. The value `True` means the
    ///         corresponding channel is included. An empty list corresponds to all channels being
    ///         enabled.
    ///     xi : list((float, float))
    ///         A list with the scale variation factors that should be used to calculate
    ///         scale-varied results. The first entry of a tuple corresponds to the variation of
    ///         the renormalization scale, the second entry to the variation of the factorization
    ///         scale. If only results for the central scale are need the list should contain
    ///         `(1.0, 1.0)`.
    ///
    /// Returns
    /// -------
    ///     numpy.ndarray(float) :
    ///         cross sections for all bins, for each scale-variation tuple (first all bins, then
    ///         the scale variation)
    pub fn convolve_with_one<'py>(
        &self,
        pdg_id: i32,
        xfx: &Bound<'py, PyAny>,
        alphas: &Bound<'py, PyAny>,
        order_mask: PyReadonlyArray1<bool>,
        bin_indices: PyReadonlyArray1<usize>,
        lumi_mask: PyReadonlyArray1<bool>,
        xi: Vec<(f64, f64)>,
        py: Python<'py>,
    ) -> Bound<'py, PyArray1<f64>> {
        let mut xfx = |id, x, q2| xfx.call1((id, x, q2)).unwrap().extract().unwrap();
        // `(q2, )` must have the comma to make it a Rust tuple
        let mut alphas = |q2| alphas.call1((q2,)).unwrap().extract().unwrap();
        let mut lumi_cache = LumiCache::with_one(pdg_id, &mut xfx, &mut alphas);
        self.grid
            .convolve(
                &mut lumi_cache,
                &order_mask.to_vec().unwrap(),
                &bin_indices.to_vec().unwrap(),
                &lumi_mask.to_vec().unwrap(),
                &xi,
            )
            .into_pyarray_bound(py)
    }

    /// Convolute grid with two pdfs.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     pdg_id1 : int
    ///         PDG Monte Carlo ID of the hadronic particle `xfx1` is the PDF for
    ///     xfx1 : callable
    ///         lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    ///     pdg_id2 : int
    ///         PDG Monte Carlo ID of the hadronic particle `xfx2` is the PDF for
    ///     xfx2 : callable
    ///         lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    ///     alphas : callable
    ///         lhapdf like callable with arguments `Q2` returning :math:`\alpha_s`
    ///     order_mask : numpy.ndarray(bool)
    ///         Mask for selecting specific orders. The value `True` means the corresponding order
    ///         is included. An empty list corresponds to all orders being enabled.
    ///     bin_indices : numpy.ndarray(int)
    ///         A list with the indices of the corresponding bins that should be calculated. An
    ///         empty list means that all orders should be calculated.
    ///     lumi_mask : numpy.ndarray(bool)
    ///         Mask for selecting specific luminosity channels. The value `True` means the
    ///         corresponding channel is included. An empty list corresponds to all channels being
    ///         enabled.
    ///     xi : list((float, float))
    ///         A list with the scale variation factors that should be used to calculate
    ///         scale-varied results. The first entry of a tuple corresponds to the variation of
    ///         the renormalization scale, the second entry to the variation of the factorization
    ///         scale. If only results for the central scale are need the list should contain
    ///         `(1.0, 1.0)`.
    ///
    /// Returns
    /// -------
    ///     numpy.ndarray(float) :
    ///         cross sections for all bins, for each scale-variation tuple (first all bins, then
    ///         the scale variation)
    pub fn convolve_with_two<'py>(
        &self,
        pdg_id1: i32,
        xfx1: &Bound<'py, PyAny>,
        pdg_id2: i32,
        xfx2: &Bound<'py, PyAny>,
        alphas: &Bound<'py, PyAny>,
        order_mask: PyReadonlyArray1<bool>,
        bin_indices: PyReadonlyArray1<usize>,
        lumi_mask: PyReadonlyArray1<bool>,
        xi: Vec<(f64, f64)>,
        py: Python<'py>,
    ) -> Bound<'py, PyArray1<f64>> {
        let mut xfx1 = |id, x, q2| xfx1.call1((id, x, q2)).unwrap().extract().unwrap();
        let mut xfx2 = |id, x, q2| xfx2.call1((id, x, q2)).unwrap().extract().unwrap();
        // `(q2, )` must have the comma to make it a Rust tuple
        let mut alphas = |q2| alphas.call1((q2,)).unwrap().extract().unwrap();
        let mut lumi_cache =
            LumiCache::with_two(pdg_id1, &mut xfx1, pdg_id2, &mut xfx2, &mut alphas);
        self.grid
            .convolve(
                &mut lumi_cache,
                &order_mask.to_vec().unwrap(),
                &bin_indices.to_vec().unwrap(),
                &lumi_mask.to_vec().unwrap(),
                &xi,
            )
            .into_pyarray_bound(py)
    }

    /// Convolute with grid with an evolution operator.
    ///
    /// Parameters
    /// ----------
    /// operator : numpy.ndarray(int, rank=5)
    ///     evolution tensor
    /// fac0 : float
    ///     reference scale
    /// pids0 : numpy.ndarray(int)
    ///     sorting of the particles in the tensor for final FkTable
    /// x0 : numpy.ndarray(float)
    ///     final FKTable interpolation grid
    /// fac1 : numpy.ndarray(float)
    ///     list of factorization scales
    /// pids1 : numpy.ndarray(int)
    ///     sorting of the particles in the grid
    /// x1 : numpy.ndarray(float)
    ///     interpolation grid at process level
    /// ren1 : numpy.ndarray(float)
    ///     list of renormalization scales
    /// alphas : numpy.ndarray(float)
    ///     list with :math:`\alpha_s(Q2)` for the process scales
    /// xi : (float, float)
    ///     factorization and renormalization variation
    /// lumi_id_types : str
    ///     type of luminosity identifier
    /// order_mask : numpy.ndarray(bool)
    ///     boolean mask to activate orders
    ///
    /// Returns
    /// -------
    /// PyFkTable :
    ///     produced FK table
    pub fn evolve(
        &self,
        operator: PyReadonlyArray5<f64>,
        fac0: f64,
        pids0: PyReadonlyArray1<i32>,
        x0: PyReadonlyArray1<f64>,
        fac1: PyReadonlyArray1<f64>,
        pids1: PyReadonlyArray1<i32>,
        x1: PyReadonlyArray1<f64>,
        ren1: PyReadonlyArray1<f64>,
        alphas: PyReadonlyArray1<f64>,
        xi: (f64, f64),
        lumi_id_types: String,
        order_mask: PyReadonlyArray1<bool>,
    ) -> PyFkTable {
        let op_info = OperatorInfo {
            fac0: fac0,
            pids0: pids0.to_vec().unwrap(),
            x0: x0.to_vec().unwrap(),
            fac1: fac1.to_vec().unwrap(),
            pids1: pids1.to_vec().unwrap(),
            x1: x1.to_vec().unwrap(),
            ren1: ren1.to_vec().unwrap(),
            alphas: alphas.to_vec().unwrap(),
            xir: xi.0,
            xif: xi.1,
            pid_basis: lumi_id_types.parse().unwrap(),
        };

        let evolved_grid = self
            .grid
            .evolve(
                operator.as_array(),
                &op_info,
                order_mask.as_slice().unwrap(),
            )
            .expect("Nothing returned from evolution.");
        PyFkTable {
            fk_table: evolved_grid,
        }
    }

    /// Convolute with grid with an evolution operator.
    ///
    /// Parameters
    /// ----------
    /// order_mask : numpy.ndarray(bool)
    ///     boolean mask to activate orders
    ///
    /// Returns
    /// -------
    /// PyEvolveInfo :
    ///     produced FK table
    pub fn evolve_info(&self, order_mask: PyReadonlyArray1<bool>) -> PyEvolveInfo {
        PyEvolveInfo {
            evolve_info: self.grid.evolve_info(order_mask.as_slice().unwrap()),
        }
    }

    /// TODO
    ///
    /// Parameters
    /// ----------
    /// slices : TODO
    /// order_mask : TODO
    ///
    /// Returns
    /// -------
    /// TODO
    pub fn evolve_with_slice_iter<'py>(
        &self,
        slices: &Bound<'py, PyIterator>,
        order_mask: PyReadonlyArray1<bool>,
        xi: (f64, f64),
        ren1: Vec<f64>,
        alphas: Vec<f64>,
    ) -> PyResult<PyFkTable> {
        Ok(self
            .grid
            .evolve_with_slice_iter(
                slices.into_iter().map(|slice| {
                    let (info, op) = slice
                        .unwrap()
                        .extract::<(PyOperatorSliceInfo, PyReadonlyArray4<f64>)>()
                        .unwrap();
                    Ok::<_, std::io::Error>((
                        info.info,
                        // TODO: avoid copying
                        CowArray::from(op.as_array().to_owned()),
                    ))
                }),
                // TODO: make `order_mask` a `Vec<f64>`
                &order_mask.to_vec().unwrap(),
                xi,
                &AlphasTable { ren1, alphas },
            )
            .map(|fk_table| PyFkTable { fk_table })
            // TODO: avoid unwrap and convert `Result` into `PyResult`
            .unwrap())
    }

    /// TODO
    ///
    /// Parameters
    /// ----------
    /// slices : TODO
    /// order_mask : TODO
    ///
    /// Returns
    /// -------
    /// TODO
    pub fn evolve_with_slice_iter2<'py>(
        &self,
        slices_a: &Bound<'py, PyIterator>,
        slices_b: &Bound<'py, PyIterator>,
        order_mask: PyReadonlyArray1<bool>,
        xi: (f64, f64),
        ren1: Vec<f64>,
        alphas: Vec<f64>,
    ) -> PyResult<PyFkTable> {
        Ok(self
            .grid
            .evolve_with_slice_iter2(
                slices_a.into_iter().map(|slice| {
                    let (info, op) = slice
                        .unwrap()
                        .extract::<(PyOperatorSliceInfo, PyReadonlyArray4<f64>)>()
                        .unwrap();
                    Ok::<_, std::io::Error>((
                        info.info,
                        // TODO: avoid copying
                        CowArray::from(op.as_array().to_owned()),
                    ))
                }),
                slices_b.into_iter().map(|slice| {
                    let (info, op) = slice
                        .unwrap()
                        .extract::<(PyOperatorSliceInfo, PyReadonlyArray4<f64>)>()
                        .unwrap();
                    Ok::<_, std::io::Error>((
                        info.info,
                        // TODO: avoid copying
                        CowArray::from(op.as_array().to_owned()),
                    ))
                }),
                // TODO: make `order_mask` a `Vec<f64>`
                &order_mask.to_vec().unwrap(),
                xi,
                &AlphasTable { ren1, alphas },
            )
            .map(|fk_table| PyFkTable { fk_table })
            // TODO: avoid unwrap and convert `Result` into `PyResult`
            .unwrap())
    }

    /// Load grid from file.
    ///
    /// **Usage:** `pineko`, FKTable generation
    ///
    /// Parameters
    /// ----------
    ///     path : str
    ///         file path
    ///
    /// Returns
    /// -------
    ///     PyGrid :
    ///         grid
    #[staticmethod]
    pub fn read(path: PathBuf) -> Self {
        Self::new(Grid::read(BufReader::new(File::open(path).unwrap())).unwrap())
    }

    /// Write grid to file.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    ///     path : str
    ///         file path
    pub fn write(&self, path: PathBuf) {
        self.grid.write(File::create(path).unwrap()).unwrap();
    }

    /// Write grid to compressed file.
    ///
    /// Parameters
    /// ----------
    ///     path : str
    ///         file path
    pub fn write_lz4(&self, path: PathBuf) {
        self.grid.write_lz4(File::create(path).unwrap()).unwrap();
    }

    /// Optimize grid content.
    ///
    /// **Usage:** `yadism`
    pub fn optimize(&mut self) {
        self.grid.optimize();
    }

    /// Merge grid with another one
    pub fn merge(&mut self, other: Self) -> PyResult<()> {
        match self.grid.merge(other.grid) {
            Ok(()) => Ok(()),
            Err(x) => Err(PyValueError::new_err(format!("{:?}", x))),
        }
    }

    /// Merge grid with another one, loaded from file
    ///
    /// Note
    /// ----
    /// For a current limitation with the implementation of the bound object `Grid` is not possible
    /// to operate with two `Grid`s in memory, since is not possible to pass a `Grid` by argument
    #[deprecated = "Deprecated in favor of PyGrid::merge"]
    pub fn merge_from_file(&mut self, path: PathBuf) -> PyResult<()> {
        match self.grid.merge(Self::read(path).grid) {
            Ok(()) => Ok(()),
            Err(x) => Err(PyValueError::new_err(format!("{:?}", x))),
        }
    }

    /// Extract the number of dimensions for bins.
    ///
    /// **Usage:** `pineko`
    ///
    /// E.g.: two differential cross-sections will return 2.
    ///
    /// Returns
    /// -------
    ///     int :
    ///         bin dimension
    pub fn bin_dimensions(&self) -> usize {
        self.grid.bin_info().dimensions()
    }

    /// Extract the normalizations for each bin.
    ///
    /// **Usage:** `runcards`
    ///
    /// Returns
    /// -------
    ///     np.ndarray
    ///         bin normalizations
    pub fn bin_normalizations<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.grid.bin_info().normalizations().into_pyarray_bound(py)
    }

    /// Extract the left edges of a specific bin dimension.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     dimension : int
    ///         bin dimension
    ///
    /// Returns
    /// -------
    ///     numpy.ndarray(float) :
    ///         left edges of bins
    pub fn bin_left<'py>(&self, dimension: usize, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.grid.bin_info().left(dimension).into_pyarray_bound(py)
    }

    /// Extract the right edges of a specific bin dimension.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     dimension : int
    ///         bin dimension
    ///
    /// Returns
    /// -------
    ///     numpy.ndarray(float) :
    ///         right edges of bins
    pub fn bin_right<'py>(&self, dimension: usize, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.grid.bin_info().right(dimension).into_pyarray_bound(py)
    }

    /// Return the number of bins.
    ///
    /// Returns
    /// -------
    ///     int :
    ///         Number of bins
    pub fn bins(&self) -> usize {
        self.grid.bin_info().bins()
    }

    /// Extract the available perturbative orders and scale variations.
    ///
    /// Returns
    /// -------
    ///     list(PyOrder) :
    ///         list with perturbative orders and scale variations
    pub fn orders(&self) -> Vec<PyOrder> {
        self.grid
            .orders()
            .iter()
            .map(|order| PyOrder {
                order: order.clone(),
            })
            .collect()
    }

    /// Get luminsosity functions.
    ///
    /// Returns
    /// -------
    ///     list(list(tuple(float,float))) :
    ///         luminosity functions as pid tuples (multiple tuples can bee associated to the same
    ///         contribution)
    pub fn channels(&self) -> Vec<Vec<(i32, i32, f64)>> {
        self.grid
            .channels()
            .iter()
            .map(|entry| entry.entry().to_vec())
            .collect()
    }

    /// Scale all subgrids.
    ///
    /// Parameters
    /// ----------
    /// factor : float
    ///     scalar factor by which scaling
    pub fn scale(&mut self, factor: f64) {
        self.grid.scale(factor);
    }

    /// Scale subgrids bin by bin.
    ///
    /// Parameters
    /// ----------
    /// factors : numpy.ndarray[float]
    ///     bin-dependent factors by which scaling
    pub fn scale_by_bin(&mut self, factors: PyReadonlyArray1<f64>) {
        self.grid.scale_by_bin(&factors.to_vec().unwrap());
    }

    /// Delete bins.
    ///
    /// Repeated bins and those exceeding length are ignored.
    ///
    /// Parameters
    /// ----------
    /// bin_indices : numpy.ndarray[int]
    ///     list of indices of bins to removed
    pub fn delete_bins(&mut self, bin_indices: PyReadonlyArray1<usize>) {
        self.grid.delete_bins(&bin_indices.to_vec().unwrap())
    }
}
