//! Grid interface.

use super::bin::PyBinRemapper;
use super::boc::{PyChannel, PyOrder};
use super::evolution::{PyEvolveInfo, PyOperatorSliceInfo};
use super::fk_table::PyFkTable;
use super::subgrid::{PySubgridEnum, PySubgridParams};
use itertools::izip;
use ndarray::CowArray;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray4};
use pineappl::convolutions::LumiCache;
use pineappl::evolution::AlphasTable;
use pineappl::grid::{Grid, Ntuple};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyIterator;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

/// PyO3 wrapper to :rustdoc:`pineappl::grid::Grid <grid/struct.Grid.html>`.
#[pyclass(name = "Grid")]
#[repr(transparent)]
#[derive(Clone)]
pub struct PyGrid {
    pub(crate) grid: Grid,
}

#[pymethods]
impl PyGrid {
    /// Constructor.
    ///
    /// Parameters
    /// ----------
    /// channels : list(PyChannel)
    ///     channels
    /// orders : list(PyOrder)
    ///     orders
    /// bin_limits : list(float)
    ///     bin configurations
    /// subgrid_params : PySubgridParams
    ///     subgrid parameters
    #[new]
    pub fn new_grid(
        channels: Vec<PyRef<PyChannel>>,
        orders: Vec<PyRef<PyOrder>>,
        bin_limits: Vec<f64>,
        subgrid_params: PySubgridParams,
    ) -> Self {
        Self {
            grid: Grid::new(
                channels.iter().map(|pyc| pyc.entry.clone()).collect(),
                orders.iter().map(|pyo| pyo.order.clone()).collect(),
                bin_limits,
                subgrid_params.subgrid_params,
            ),
        }
    }

    /// Add a point to the grid.
    ///
    /// Parameters
    /// ----------
    /// x1 : float
    ///     first momentum fraction
    /// x2 : float
    ///     second momentum fraction
    /// q2 : float
    ///     process scale
    /// order : int
    ///     order index
    /// observable : float
    ///     reference point (to be binned)
    /// channel : int
    ///     channel index
    /// weight : float
    ///     cross section weight
    pub fn fill(
        &mut self,
        x1: f64,
        x2: f64,
        q2: f64,
        order: usize,
        observable: f64,
        channel: usize,
        weight: f64,
    ) {
        self.grid.fill(
            order,
            observable,
            channel,
            &Ntuple::<f64> { x1, x2, q2, weight },
        );
    }

    /// Add an array to the grid.
    ///
    /// Useful to avoid multiple python calls, leading to performance improvement.
    ///
    /// Parameters
    /// ----------
    /// x1s : np.array(float)
    ///     first momentum fraction of all events
    /// x2s : np.array(float)
    ///     second momentum fraction of all events
    /// x1s : np.array(float)
    ///     process scale of all events
    /// order : int
    ///     order index
    /// observable : float
    ///     reference point (to be binned)
    /// channel : int
    ///     channel index
    /// weights : np.array(float)
    ///     cross section weight for all events
    pub fn fill_array(
        &mut self,
        x1s: Vec<f64>,
        x2s: Vec<f64>,
        q2s: Vec<f64>,
        order: usize,
        observables: Vec<f64>,
        channel: usize,
        weights: Vec<f64>,
    ) {
        for (&x1, &x2, &q2, &observable, &weight) in izip!(
            x1s.iter(),
            x2s.iter(),
            q2s.iter(),
            observables.iter(),
            weights.iter(),
        ) {
            self.grid.fill(
                order,
                observable,
                channel,
                &Ntuple::<f64> { x1, x2, q2, weight },
            );
        }
    }

    /// Add a point to the grid for all channels.
    ///
    /// Parameters
    /// ----------
    /// x1 : float
    ///     first momentum fraction
    /// x2 : float
    ///     second momentum fraction
    /// q2 : float
    ///     process scale
    /// order : int
    ///     order index
    /// observable : float
    ///     reference point (to be binned)
    /// weights : np.array(float)
    ///     cross section weights, one for each channels
    pub fn fill_all(
        &mut self,
        x1: f64,
        x2: f64,
        q2: f64,
        order: usize,
        observable: f64,
        weights: Vec<f64>,
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
            &weights,
        );
    }

    /// Get metadata values stored in the grid.
    ///
    ///
    /// Returns
    /// -------
    /// dict :
    ///     key, value map
    pub fn key_values(&self) -> HashMap<String, String> {
        self.grid.key_values().unwrap().clone()
    }

    /// Set a metadata key-value pair in the grid.
    ///
    /// Parameters
    /// ----------
    /// key : str
    ///     key
    /// value : str
    ///     value
    pub fn set_key_value(&mut self, key: &str, value: &str) {
        self.grid.set_key_value(key, value);
    }

    /// Retrieve a subgrid.
    pub fn subgrid(&self, order: usize, bin: usize, channel: usize) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.grid.subgrids()[[order, bin, channel]].clone(),
        }
    }

    /// Set a subgrid.
    pub fn set_subgrid(
        &mut self,
        order: usize,
        bin: usize,
        channel: usize,
        subgrid: PySubgridEnum,
    ) {
        self.grid.subgrids_mut()[[order, bin, channel]] = subgrid.subgrid_enum;
    }

    /// Set the bin normalizations.
    ///
    /// Parameters
    /// ----------
    /// remapper: BinRemapper
    ///     Remapper object
    pub fn set_remapper(&mut self, remapper: PyBinRemapper) {
        self.grid.set_remapper(remapper.bin_remapper).unwrap();
    }

    /// Convolve with a single distribution.
    ///
    /// Parameters
    /// ----------
    /// pdg_id : int
    ///     PDG Monte Carlo ID of the hadronic particle
    /// xfx : callable
    ///     lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    /// alphas : callable
    ///     lhapdf like callable with arguments `Q2` returning :math:`\alpha_s`
    /// order_mask : numpy.ndarray(bool)
    ///     Mask for selecting specific orders. The value `True` means the corresponding order
    ///     is included. An empty list corresponds to all orders being enabled.
    /// bin_indices : numpy.ndarray(int)
    ///     A list with the indices of the corresponding bins that should be calculated. An
    ///     empty list means that all bins should be calculated.
    /// channel_mask : numpy.ndarray(bool)
    ///     Mask for selecting specific channels. The value `True` means the
    ///     corresponding channel is included. An empty list corresponds to all channels being
    ///     enabled.
    /// xi : list((float, float))
    ///     A list with the scale variation factors that should be used to calculate
    ///     scale-varied results. The first entry of a tuple corresponds to the variation of
    ///     the renormalization scale, the second entry to the variation of the factorization
    ///     scale. If only results for the central scale are need the list should contain
    ///     `(1.0, 1.0)`.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(float) :
    ///     cross sections for all bins, for each scale-variation tuple (first all bins, then
    ///     the scale variation)
    #[pyo3(signature = (pdg_id, xfx, alphas, order_mask = None, bin_indices = None, channel_mask = None, xi = None))]
    pub fn convolve_with_one<'py>(
        &self,
        pdg_id: i32,
        xfx: &Bound<'py, PyAny>,
        alphas: &Bound<'py, PyAny>,
        order_mask: Option<Vec<bool>>,
        bin_indices: Option<Vec<usize>>,
        channel_mask: Option<Vec<bool>>,
        xi: Option<Vec<(f64, f64)>>,
        py: Python<'py>,
    ) -> Bound<'py, PyArray1<f64>> {
        let mut xfx = |id, x, q2| xfx.call1((id, x, q2)).unwrap().extract().unwrap();
        // `(q2, )` must have the comma to make it a Rust tuple
        let mut alphas = |q2| alphas.call1((q2,)).unwrap().extract().unwrap();
        let mut lumi_cache = LumiCache::with_one(pdg_id, &mut xfx, &mut alphas);
        self.grid
            .convolve(
                &mut lumi_cache,
                &order_mask.unwrap_or_default(),
                &bin_indices.unwrap_or_default(),
                &channel_mask.unwrap_or_default(),
                &xi.unwrap_or(vec![(1.0, 1.0)]),
            )
            .into_pyarray_bound(py)
    }

    /// Convolve with two distributions.
    ///
    /// Parameters
    /// ----------
    /// pdg_id1 : int
    ///     PDG Monte Carlo ID of the first hadronic particle
    /// xfx1 : callable
    ///     lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    /// pdg_id2 : int
    ///     PDG Monte Carlo ID of the second hadronic particle
    /// xfx2 : callable
    ///     lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    /// alphas : callable
    ///     lhapdf like callable with arguments `Q2` returning :math:`\alpha_s`
    /// order_mask : numpy.ndarray(bool)
    ///     Mask for selecting specific orders. The value `True` means the corresponding order
    ///     is included. An empty list corresponds to all orders being enabled.
    /// bin_indices : numpy.ndarray(int)
    ///     A list with the indices of the corresponding bins that should be calculated. An
    ///     empty list means that all bins should be calculated.
    /// channel_mask : numpy.ndarray(bool)
    ///     Mask for selecting specific channels. The value `True` means the
    ///     corresponding channel is included. An empty list corresponds to all channels being
    ///     enabled.
    /// xi : list((float, float))
    ///     A list with the scale variation factors that should be used to calculate
    ///     scale-varied results. The first entry of a tuple corresponds to the variation of
    ///     the renormalization scale, the second entry to the variation of the factorization
    ///     scale. If only results for the central scale are need the list should contain
    ///     `(1.0, 1.0)`.
    ///
    /// Returns
    /// -------
    ///     numpy.ndarray(float) :
    ///         cross sections for all bins, for each scale-variation tuple (first all bins, then
    ///         the scale variation)
    #[pyo3(signature = (pdg_id1, xfx1, pdg_id2, xfx2, alphas, order_mask = None, bin_indices = None, channel_mask = None, xi = None))]
    pub fn convolve_with_two<'py>(
        &self,
        pdg_id1: i32,
        xfx1: &Bound<'py, PyAny>,
        pdg_id2: i32,
        xfx2: &Bound<'py, PyAny>,
        alphas: &Bound<'py, PyAny>,
        order_mask: Option<Vec<bool>>,
        bin_indices: Option<Vec<usize>>,
        channel_mask: Option<Vec<bool>>,
        xi: Option<Vec<(f64, f64)>>,
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
                &order_mask.unwrap_or_default(),
                &bin_indices.unwrap_or_default(),
                &channel_mask.unwrap_or_default(),
                &xi.unwrap_or(vec![(1.0, 1.0)]),
            )
            .into_pyarray_bound(py)
    }

    /// Collect information for convolution with an evolution operator.
    ///
    /// Parameters
    /// ----------
    /// order_mask : numpy.ndarray(bool)
    ///     boolean mask to activate orders
    ///
    /// Returns
    /// -------
    /// PyEvolveInfo :
    ///     evolution information
    pub fn evolve_info(&self, order_mask: Vec<bool>) -> PyEvolveInfo {
        PyEvolveInfo {
            evolve_info: self.grid.evolve_info(order_mask.as_slice()),
        }
    }

    /// Convolve with uniform evolution operator slices.
    ///
    /// Parameters
    /// ----------
    /// slices : Iterable
    ///     list of (PyOperatorSliceInfo, 5D array) describing each convolution
    /// order_mask : numpy.ndarray(bool)
    ///     boolean mask to activate orders
    /// xi : (float, float)
    ///     factorization and renormalization variation
    /// ren1 : numpy.ndarray(float)
    ///     list of renormalization scales
    /// alphas : numpy.ndarray(float)
    ///     list with :math:`\alpha_s(Q2)` for the process scales
    ///
    /// Returns
    /// -------
    /// PyFkTable :
    ///     produced FK table
    pub fn evolve_with_slice_iter<'py>(
        &self,
        slices: &Bound<'py, PyIterator>,
        order_mask: Vec<bool>,
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
                &order_mask,
                xi,
                &AlphasTable { ren1, alphas },
            )
            .map(|fk_table| PyFkTable { fk_table })
            // TODO: avoid unwrap and convert `Result` into `PyResult`
            .unwrap())
    }

    /// Convolve with two different evolution operator slices.
    ///
    /// Parameters
    /// ----------
    /// slices_a : Iterable
    ///     list of (PyOperatorSliceInfo, 5D array) describing the first convolution
    /// slices_b : Iterable
    ///     list of (PyOperatorSliceInfo, 5D array) describing the second convolution
    /// order_mask : numpy.ndarray(bool)
    ///     boolean mask to activate orders
    /// xi : (float, float)
    ///     factorization and renormalization variation
    /// ren1 : numpy.ndarray(float)
    ///     list of renormalization scales
    /// alphas : numpy.ndarray(float)
    ///     list with :math:`\alpha_s(Q2)` for the process scales
    ///
    /// Returns
    /// -------
    /// PyFkTable :
    ///     produced FK table
    pub fn evolve_with_slice_iter2<'py>(
        &self,
        slices_a: &Bound<'py, PyIterator>,
        slices_b: &Bound<'py, PyIterator>,
        order_mask: Vec<bool>,
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
                &order_mask,
                xi,
                &AlphasTable { ren1, alphas },
            )
            .map(|fk_table| PyFkTable { fk_table })
            // TODO: avoid unwrap and convert `Result` into `PyResult`
            .unwrap())
    }

    /// Load from file.
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     file path
    ///
    /// Returns
    /// -------
    /// PyGrid :
    ///     grid
    #[staticmethod]
    pub fn read(path: PathBuf) -> Self {
        Self {
            grid: Grid::read(BufReader::new(File::open(path).unwrap())).unwrap(),
        }
    }

    /// Write to file.
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     file path
    pub fn write(&self, path: PathBuf) {
        self.grid.write(File::create(path).unwrap()).unwrap();
    }

    /// Write to compressed file.
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     file path
    pub fn write_lz4(&self, path: PathBuf) {
        self.grid.write_lz4(File::create(path).unwrap()).unwrap();
    }

    /// Optimize content.
    pub fn optimize(&mut self) {
        self.grid.optimize();
    }

    /// Merge with another grid.
    pub fn merge(&mut self, other: Self) -> PyResult<()> {
        match self.grid.merge(other.grid) {
            Ok(()) => Ok(()),
            Err(x) => Err(PyValueError::new_err(format!("{:?}", x))),
        }
    }

    /// Extract the number of dimensions for bins.
    ///
    /// E.g.: two differential cross-sections will return 2.
    ///
    /// Returns
    /// -------
    /// int :
    ///     bin dimension
    pub fn bin_dimensions(&self) -> usize {
        self.grid.bin_info().dimensions()
    }

    /// Extract the normalizations for each bin.
    ///
    /// Returns
    /// -------
    /// np.ndarray
    ///     bin normalizations
    pub fn bin_normalizations<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.grid.bin_info().normalizations().into_pyarray_bound(py)
    }

    /// Extract the left edges of a specific bin dimension.
    ///
    /// Parameters
    /// ----------
    /// dimension : int
    ///     bin dimension
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(float) :
    ///     left edges of bins
    pub fn bin_left<'py>(&self, dimension: usize, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.grid.bin_info().left(dimension).into_pyarray_bound(py)
    }

    /// Extract the right edges of a specific bin dimension.
    ///
    /// Parameters
    /// ----------
    /// dimension : int
    ///     bin dimension
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(float) :
    ///     right edges of bins
    pub fn bin_right<'py>(&self, dimension: usize, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.grid.bin_info().right(dimension).into_pyarray_bound(py)
    }

    /// Return the number of bins.
    ///
    /// Returns
    /// -------
    /// int :
    ///     Number of bins
    pub fn bins(&self) -> usize {
        self.grid.bin_info().bins()
    }

    /// Extract the available perturbative orders and scale variations.
    ///
    /// Returns
    /// -------
    /// list(PyOrder) :
    ///     list with perturbative orders and scale variations
    pub fn orders(&self) -> Vec<PyOrder> {
        self.grid
            .orders()
            .iter()
            .map(|order| PyOrder {
                order: order.clone(),
            })
            .collect()
    }

    /// Extract channels.
    ///
    /// Returns
    /// -------
    /// list(list(tuple(float,float,int))) :
    ///     channels as tuples (pid, pid, factor) (multiple tuples can be associated to the same
    ///     contribution)
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
    ///     scalar factor by which to scale
    pub fn scale(&mut self, factor: f64) {
        self.grid.scale(factor);
    }

    /// Scale subgrids bin by bin.
    ///
    /// Parameters
    /// ----------
    /// factors : numpy.ndarray[float]
    ///     bin-dependent factors by which to scale
    pub fn scale_by_bin(&mut self, factors: Vec<f64>) {
        self.grid.scale_by_bin(&factors);
    }

    /// Delete bins.
    ///
    /// Repeated bins and those exceeding the length are ignored.
    ///
    /// Parameters
    /// ----------
    /// bin_indices : numpy.ndarray[int]
    ///     list of indices of bins to removed
    pub fn delete_bins(&mut self, bin_indices: Vec<usize>) {
        self.grid.delete_bins(&bin_indices)
    }
}

/// Register submodule in parent.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "grid")?;
    m.setattr(pyo3::intern!(m.py(), "__doc__"), "Grid interface.")?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.grid'] = m"
    );
    m.add_class::<PyGrid>()?;
    m.add_class::<PyOrder>()?;
    parent_module.add_submodule(&m)
}
