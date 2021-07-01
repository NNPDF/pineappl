use lhapdf::Pdf;
use pineappl::grid::{EkoInfo, Grid, Order};

use super::bin::PyBinRemapper;
use super::lumi::PyLumiEntry;
use super::subgrid::{PySubgridEnum, PySubgridParams};

use ndarray::{Array, Ix5};

use std::fs::File;
use std::io::BufReader;

use pyo3::prelude::*;

#[pyclass]
#[repr(transparent)]
pub struct PyOrder {
    pub order: Order,
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
}

#[pyclass]
#[repr(transparent)]
pub struct PyGrid {
    pub grid: Grid,
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
        bin_limits: Vec<f64>,
        subgrid_params: PySubgridParams,
    ) -> Self {
        Self::new(Grid::new(
            lumi.iter().map(|pyl| pyl.lumi_entry.clone()).collect(),
            orders.iter().map(|pyo| pyo.order.clone()).collect(),
            bin_limits,
            subgrid_params.subgrid_params,
        ))
    }

    pub fn set_key_value(&mut self, key: &str, value: &str) {
        self.grid.set_key_value(key, value);
    }

    pub fn bins(&self) -> usize {
        self.grid.bin_info().bins()
    }

    pub fn subgrid(&self, order: usize, bin: usize, lumi: usize) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.grid.subgrid(order, bin, lumi).clone(),
        }
    }

    pub fn set_subgrid(&mut self, order: usize, bin: usize, lumi: usize, subgrid: PySubgridEnum) {
        self.grid
            .set_subgrid(order, bin, lumi, subgrid.subgrid_enum);
    }

    pub fn set_remapper(&mut self, remapper: PyBinRemapper) {
        self.grid.set_remapper(remapper.bin_remapper).unwrap();
    }

    pub fn convolute(&self, pdfset: &str) -> Vec<f64> {
        let pdf_lha = pdfset
            .parse()
            .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);
        let initial_state_1 = self.grid.key_values().map_or(2212, |map| {
            map.get("initial_state_1").unwrap().parse::<i32>().unwrap()
        });
        let initial_state_2 = self.grid.key_values().map_or(2212, |map| {
            map.get("initial_state_2").unwrap().parse::<i32>().unwrap()
        });

        // if the field 'Particle' is missing we assume it's a proton PDF
        let pdf_pdg_id = pdf_lha
            .set()
            .entry("Particle")
            .unwrap_or_else(|| "2212".to_string())
            .parse::<i32>()
            .unwrap();

        let pdf = |id, x, q2| pdf_lha.xfx_q2(id, x, q2);
        let anti_pdf = |id, x, q2| {
            let id = match id {
                -6..=6 | 11 | 13 | -11 | -13 => -id,
                21 | 22 => id,
                _ => unimplemented!(),
            };
            pdf_lha.xfx_q2(id, x, q2)
        };
        let no_pdf = |_, x, _| x;

        let xfx1: Box<dyn Fn(i32, f64, f64) -> f64> = if initial_state_1 == pdf_pdg_id {
            Box::new(&pdf)
        } else if initial_state_1 == -pdf_pdg_id {
            Box::new(&anti_pdf)
        } else {
            match initial_state_1 {
                11 | 13 | -11 | -13 => Box::new(&no_pdf),
                _ => unimplemented!(),
            }
        };
        let xfx2: Box<dyn Fn(i32, f64, f64) -> f64> = if initial_state_2 == pdf_pdg_id {
            Box::new(&pdf)
        } else if initial_state_2 == -pdf_pdg_id {
            Box::new(&anti_pdf)
        } else {
            match initial_state_2 {
                11 | 13 | -11 | -13 => Box::new(&no_pdf),
                _ => unimplemented!(),
            }
        };
        let alphas = |q2| pdf_lha.alphas_q2(q2);

        let orders: Vec<bool> = Vec::new();
        let bins: Vec<usize> = Vec::new();
        let lumis: Vec<bool> = Vec::new();
        let scales = [(1., 1.)];
        self.grid.convolute(
            &xfx1,
            &xfx2,
            &alphas,
            &orders[..],
            &bins[..],
            &lumis[..],
            &scales,
        )
    }

    pub fn eko_info(&self) -> (Vec<f64>, Vec<f64>) {
        let EkoInfo { x_grid, q2_grid } = self.grid.eko_info().unwrap();
        (x_grid, q2_grid)
    }

    pub fn convolute_eko(
        &self,
        q2: f64,
        alphas: Vec<f64>,
        pids: Vec<i32>,
        x_grid: Vec<f64>,
        q2_grid: Vec<f64>,
        operator_flattened: Vec<f64>,
        operator_shape: Vec<usize>,
    ) -> Self {
        let operator = Array::from_shape_vec(operator_shape, operator_flattened).unwrap();
        let evolved_grid = self
            .grid
            .convolute_eko(
                q2,
                &alphas,
                (1., 1.),
                &pids,
                x_grid,
                q2_grid,
                operator.into_dimensionality::<Ix5>().unwrap(),
            )
            .unwrap();
        Self::new(evolved_grid)
    }

    #[staticmethod]
    pub fn read(path: &str) -> Self {
        Self::new(Grid::read(BufReader::new(File::open(path).unwrap())).unwrap())
    }

    pub fn write(&self, path: &str) {
        self.grid.write(File::create(path).unwrap()).unwrap();
    }
}
