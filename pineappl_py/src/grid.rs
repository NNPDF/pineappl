//! Grid interface.

use super::bin::PyBinRemapper;
use super::boc::{PyChannel, PyKinematics, PyOrder};
use super::convolutions::PyConv;
use super::evolution::{PyEvolveInfo, PyOperatorSliceInfo};
use super::fk_table::PyFkTable;
use super::interpolation::PyInterp;
use super::pids::PyPidBasis;
use super::subgrid::PySubgridEnum;
use ndarray::CowArray;
use numpy::{IntoPyArray, PyArray1, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray4};
use pineappl::boc::{ScaleFuncForm, Scales};
use pineappl::convolutions::ConvolutionCache;
use pineappl::evolution::AlphasTable;
use pineappl::grid::Grid;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyIterator;
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
    /// Constructor to instantiate a new Grid.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// pid_basis : PidBasis
    ///     choice of basis which can be `Evol` or `Pdg`
    /// channels : list(PyChannel)
    ///     channels
    /// orders : list(PyOrder)
    ///     orders
    /// bin_limits : list(float)
    ///     bin configurations
    /// convolutions : list(PyConv)
    ///     contains the types of convolution
    /// interpolations : list(PyInterp)
    ///     types of interpolations required by each kinematic
    /// kinematics : list(PyKinematics)
    ///     list of kinematics
    #[new]
    #[must_use]
    pub fn new_grid(
        pid_basis: PyPidBasis,
        channels: Vec<PyRef<PyChannel>>,
        orders: Vec<PyRef<PyOrder>>,
        bin_limits: Vec<f64>,
        convolutions: Vec<PyRef<PyConv>>,
        interpolations: Vec<PyRef<PyInterp>>,
        kinematics: Vec<PyRef<PyKinematics>>,
    ) -> Self {
        Self {
            grid: Grid::new(
                pid_basis.into(),
                channels.into_iter().map(|pyc| pyc.entry.clone()).collect(),
                orders.into_iter().map(|pyo| pyo.order.clone()).collect(),
                bin_limits,
                convolutions
                    .into_iter()
                    .map(|pyx| pyx.conv.clone())
                    .collect(),
                interpolations
                    .into_iter()
                    .map(|pyi| pyi.interp.clone())
                    .collect(),
                kinematics.into_iter().map(|pyk| pyk.kinematics).collect(),
                Scales {
                    ren: ScaleFuncForm::Scale(0),
                    fac: ScaleFuncForm::Scale(0),
                    frg: ScaleFuncForm::NoScale,
                },
            ),
        }
    }

    /// Add a point to the grid.
    ///
    /// Parameters
    /// ----------
    /// order : int
    ///     order index
    /// observable : float
    ///     reference point (to be binned)
    /// channel : int
    ///     channel index
    /// ntuple: list(float)
    ///     list containing information on kinematics
    /// weight : float
    ///     cross section weight
    pub fn fill(
        &mut self,
        order: usize,
        observable: f64,
        channel: usize,
        ntuple: Vec<f64>,
        weight: f64,
    ) {
        self.grid.fill(order, observable, channel, &ntuple, weight);
    }

    /// Retrieve a subgrid.
    #[must_use]
    pub fn subgrid(&self, order: usize, bin: usize, channel: usize) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.grid.subgrids()[[order, bin, channel]].clone(),
        }
    }

    /// Set a subgrid.
    ///
    /// Parameters
    /// ----------
    /// order : int
    ///     order index
    /// bin : int
    ///     bin index
    /// channel : int
    ///     channel index
    /// subgrid : PySubgridEnum
    ///     subgrid object
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
    /// # Panics
    ///
    /// TODO
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
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// pdg_conv : PyConv
    ///     contains the types of convolutions and PID
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
    #[must_use]
    #[pyo3(signature = (pdg_conv, xfx, alphas, order_mask = None, bin_indices = None, channel_mask = None, xi = None))]
    pub fn convolve_with_one<'py>(
        &self,
        pdg_conv: PyRef<PyConv>,
        xfx: &Bound<'py, PyAny>,
        alphas: &Bound<'py, PyAny>,
        order_mask: Option<PyReadonlyArray1<bool>>,
        bin_indices: Option<PyReadonlyArray1<usize>>,
        channel_mask: Option<PyReadonlyArray1<bool>>,
        xi: Option<Vec<(f64, f64, f64)>>,
        py: Python<'py>,
    ) -> Bound<'py, PyArray1<f64>> {
        let mut xfx = |id, x, q2| xfx.call1((id, x, q2)).unwrap().extract().unwrap();
        // `(q2, )` must have the comma to make it a Rust tuple
        let mut alphas = |q2| alphas.call1((q2,)).unwrap().extract().unwrap();
        let mut lumi_cache =
            ConvolutionCache::new(vec![pdg_conv.conv.clone()], vec![&mut xfx], &mut alphas);
        self.grid
            .convolve(
                &mut lumi_cache,
                &order_mask.map_or(vec![], |b| b.to_vec().unwrap()),
                &bin_indices.map_or(vec![], |c| c.to_vec().unwrap()),
                &channel_mask.map_or(vec![], |d| d.to_vec().unwrap()),
                &xi.map_or_else(|| vec![(1.0, 1.0, 1.0)], |m| m),
            )
            .into_pyarray_bound(py)
    }

    /// Convolve with two distributions.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// pdg_conv1 : PyConv
    ///     contains the types of convolutions and PID
    /// xfx1 : callable
    ///     lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    /// pdg_conv2 : PyConv
    ///     contains the types of convolutions and PID
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
    #[must_use]
    #[pyo3(signature = (pdg_conv1, xfx1, pdg_conv2, xfx2, alphas, order_mask = None, bin_indices = None, channel_mask = None, xi = None))]
    pub fn convolve_with_two<'py>(
        &self,
        pdg_conv1: PyRef<PyConv>,
        xfx1: &Bound<'py, PyAny>,
        pdg_conv2: PyRef<PyConv>,
        xfx2: &Bound<'py, PyAny>,
        alphas: &Bound<'py, PyAny>,
        order_mask: Option<PyReadonlyArray1<bool>>,
        bin_indices: Option<PyReadonlyArray1<usize>>,
        channel_mask: Option<PyReadonlyArray1<bool>>,
        xi: Option<Vec<(f64, f64, f64)>>,
        py: Python<'py>,
    ) -> Bound<'py, PyArray1<f64>> {
        let mut xfx1 = |id, x, q2| xfx1.call1((id, x, q2)).unwrap().extract().unwrap();
        let mut xfx2 = |id, x, q2| xfx2.call1((id, x, q2)).unwrap().extract().unwrap();
        // `(q2, )` must have the comma to make it a Rust tuple
        let mut alphas = |q2| alphas.call1((q2,)).unwrap().extract().unwrap();
        let mut lumi_cache = ConvolutionCache::new(
            vec![pdg_conv1.conv.clone(), pdg_conv2.conv.clone()],
            vec![&mut xfx1, &mut xfx2],
            &mut alphas,
        );
        self.grid
            .convolve(
                &mut lumi_cache,
                &order_mask.map_or(vec![], |b| b.to_vec().unwrap()),
                &bin_indices.map_or(vec![], |c| c.to_vec().unwrap()),
                &channel_mask.map_or(vec![], |d| d.to_vec().unwrap()),
                &xi.map_or_else(|| vec![(1.0, 1.0, 1.0)], |m| m),
            )
            .into_pyarray_bound(py)
    }

    /// TODO
    // #[pyo3(signature = (pdg_convs, xfxs, alphas, order_mask = None, bin_indices = None, channel_mask = None, xi = None))]
    #[must_use]
    pub fn convolve<'py>(
        &self,
        _pdg_convs: Vec<PyRef<PyConv>>,
        _xfxs: &Bound<'py, PyIterator>,
        _alphas: &Bound<'py, PyAny>,
        _order_mask: Option<PyReadonlyArray1<bool>>,
        _bin_indices: Option<PyReadonlyArray1<usize>>,
        _channel_mask: Option<PyReadonlyArray1<bool>>,
        _xi: Option<Vec<(f64, f64, f64)>>,
        _py: Python<'py>,
    ) -> Bound<'py, PyArray1<f64>> {
        todo!()
    }

    /// Collect information for convolution with an evolution operator.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// order_mask : numpy.ndarray(bool)
    ///     boolean mask to activate orders
    ///
    /// Returns
    /// -------
    /// PyEvolveInfo :
    ///     evolution informations
    #[must_use]
    pub fn evolve_info(&self, order_mask: PyReadonlyArray1<bool>) -> PyEvolveInfo {
        PyEvolveInfo {
            evolve_info: self.grid.evolve_info(order_mask.as_slice().unwrap()),
        }
    }

    /// Evolve grid with as many EKOs as Convolutions.
    ///
    /// TODO: Expose `slices` to be a vector!!!
    ///
    /// # Panics
    /// TODO
    ///
    /// # Errors
    /// TODO
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
    #[allow(clippy::needless_lifetimes)]
    pub fn evolve<'py>(
        &self,
        slices: &Bound<'py, PyIterator>,
        order_mask: PyReadonlyArray1<bool>,
        xi: (f64, f64, f64),
        ren1: Vec<f64>,
        alphas: Vec<f64>,
    ) -> PyResult<PyFkTable> {
        Ok(self
            .grid
            .evolve(
                vec![slices.into_iter().map(|slice| {
                    let (info, op) = slice
                        .unwrap()
                        .extract::<(PyOperatorSliceInfo, PyReadonlyArray4<f64>)>()
                        .unwrap();
                    Ok::<_, std::io::Error>((
                        info.info,
                        // TODO: avoid copying
                        CowArray::from(op.as_array().to_owned()),
                    ))
                })],
                // TODO: make `order_mask` a `Vec<f64>`
                &order_mask.to_vec().unwrap(),
                xi,
                &AlphasTable { ren1, alphas },
            )
            .map(|fk_table| PyFkTable { fk_table })
            // TODO: avoid unwrap and convert `Result` into `PyResult`
            .unwrap())
    }

    /// Evolve grid with one single EKO.
    ///
    /// # Panics
    /// TODO
    ///
    /// # Errors
    /// TODO
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
    #[allow(clippy::needless_lifetimes)]
    pub fn evolve_with_slice_iter<'py>(
        &self,
        slices: &Bound<'py, PyIterator>,
        order_mask: PyReadonlyArray1<bool>,
        xi: (f64, f64, f64),
        ren1: Vec<f64>,
        alphas: Vec<f64>,
    ) -> PyResult<PyFkTable> {
        Ok(self
            .grid
            .evolve(
                vec![slices.into_iter().map(|slice| {
                    let (info, op) = slice
                        .unwrap()
                        .extract::<(PyOperatorSliceInfo, PyReadonlyArray4<f64>)>()
                        .unwrap();
                    Ok::<_, std::io::Error>((
                        info.info,
                        // TODO: avoid copying
                        CowArray::from(op.as_array().to_owned()),
                    ))
                })],
                // TODO: make `order_mask` a `Vec<f64>`
                &order_mask.to_vec().unwrap(),
                xi,
                &AlphasTable { ren1, alphas },
            )
            .map(|fk_table| PyFkTable { fk_table })
            // TODO: avoid unwrap and convert `Result` into `PyResult`
            .unwrap())
    }

    /// Load from file.
    ///
    /// # Panics
    ///
    /// TODO
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
    #[must_use]
    #[staticmethod]
    pub fn read(path: PathBuf) -> Self {
        Self {
            grid: Grid::read(BufReader::new(File::open(path).unwrap())).unwrap(),
        }
    }

    /// Write to file.
    ///
    /// # Panics
    ///
    /// TODO
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
    /// # Panics
    ///
    /// TODO
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
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// # Errors
    ///
    /// TODO
    pub fn merge(&mut self, other: Self) -> PyResult<()> {
        match self.grid.merge(other.grid) {
            Ok(()) => Ok(()),
            Err(x) => Err(PyValueError::new_err(format!("{x:?}"))),
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
    #[must_use]
    pub fn bin_dimensions(&self) -> usize {
        self.grid.bin_info().dimensions()
    }

    /// Extract the normalizations for each bin.
    ///
    /// Returns
    /// -------
    /// np.ndarray
    ///     bin normalizations
    #[must_use]
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
    #[must_use]
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
    #[must_use]
    pub fn bin_right<'py>(&self, dimension: usize, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.grid.bin_info().right(dimension).into_pyarray_bound(py)
    }

    /// Return the number of bins.
    ///
    /// Returns
    /// -------
    /// int :
    ///     Number of bins
    #[must_use]
    pub fn bins(&self) -> usize {
        self.grid.bin_info().bins()
    }

    /// Extract the available perturbative orders and scale variations.
    ///
    /// Returns
    /// -------
    /// list(PyOrder) :
    ///     list with perturbative orders and scale variations
    #[must_use]
    pub fn orders(&self) -> Vec<PyOrder> {
        self.grid
            .orders()
            .iter()
            .map(|order| PyOrder {
                order: order.clone(),
            })
            .collect()
    }

    /// Get the type(s) of convolution(s) for the current Grid.
    ///
    /// Returns
    /// list(PyConv):
    ///     list of convolution type with the corresponding PIDs
    #[must_use]
    pub fn convolutions(&self) -> Vec<PyConv> {
        self.grid
            .convolutions()
            .iter()
            .map(|conv| PyConv { conv: conv.clone() })
            .collect()
    }

    /// Extract channels.
    ///
    /// Returns
    /// -------
    /// list(list(tuple(float,float,int))) :
    ///     channels as tuples (pid, pid, factor) (multiple tuples can be associated to the same
    ///     contribution)
    #[must_use]
    pub fn channels(&self) -> Vec<Vec<(Vec<i32>, f64)>> {
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
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// factors : numpy.ndarray[float]
    ///     bin-dependent factors by which to scale
    pub fn scale_by_bin(&mut self, factors: PyReadonlyArray1<f64>) {
        self.grid.scale_by_bin(&factors.to_vec().unwrap());
    }

    /// Delete bins.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Repeated bins and those exceeding the length are ignored.
    ///
    /// Parameters
    /// ----------
    /// bin_indices : numpy.ndarray[int]
    ///     list of indices of bins to removed
    pub fn delete_bins(&mut self, bin_indices: PyReadonlyArray1<usize>) {
        self.grid.delete_bins(&bin_indices.to_vec().unwrap());
    }
}

/// Register submodule in parent.
/// # Errors
///
/// Raises an error if (sub)module is not found.
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
