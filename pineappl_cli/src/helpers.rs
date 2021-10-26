use anyhow::{Context, Result};
use lhapdf::Pdf;
use ndarray::Array3;
use pineappl::grid::Grid;
use pineappl::lumi::LumiCache;
use prettytable::format::{FormatBuilder, LinePosition, LineSeparator};
use prettytable::Table;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter};

pub const ONE_SIGMA: f64 = 68.268_949_213_708_58;

pub fn read_grid(input: &str) -> Result<Grid> {
    Grid::read(BufReader::new(
        File::open(input).context(format!("unable to open '{}'", input))?,
    ))
    .context(format!("unable to read '{}'", input))
}

pub fn write_grid(output: &str, grid: &Grid) -> Result<()> {
    grid.write(BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)
            .context(format!("unable to write '{}'", output))?,
    ))
}

pub fn create_table() -> Table {
    let mut table = Table::new();
    table.set_format(
        FormatBuilder::new()
            .column_separator(' ')
            .separator(LinePosition::Title, LineSeparator::new('-', '+', ' ', ' '))
            .build(),
    );
    table
}

pub const SCALES_VECTOR: [(f64, f64); 9] = [
    (1.0, 1.0),
    (2.0, 2.0),
    (0.5, 0.5),
    (2.0, 1.0),
    (1.0, 2.0),
    (0.5, 1.0),
    (1.0, 0.5),
    (2.0, 0.5),
    (0.5, 2.0),
];

pub fn labels(grid: &Grid) -> Vec<String> {
    let mut labels = vec![];
    let key_values = grid.key_values().cloned().unwrap_or_default();

    for d in 0..grid.bin_info().dimensions() {
        labels.push(
            key_values
                .get(&format!("x{}_label", d + 1))
                .unwrap_or(&format!("x{}", d))
                .clone(),
        );
    }

    labels.push(
        key_values
            .get("y_label")
            .unwrap_or(&"diff".to_owned())
            .clone(),
    );
    labels
}

pub fn convolute(
    grid: &Grid,
    lhapdf: &Pdf,
    orders: &[(u32, u32)],
    bins: &[usize],
    lumis: &[bool],
    scales: usize,
) -> Vec<f64> {
    let orders: Vec<_> = grid
        .orders()
        .iter()
        .map(|order| {
            orders.is_empty()
                || orders
                    .iter()
                    .any(|other| (order.alphas == other.0) && (order.alpha == other.1))
        })
        .collect();

    let pdf_pdg_id = lhapdf
        .set()
        .entry("Particle")
        .unwrap_or_else(|| "2212".to_string())
        .parse::<i32>()
        .unwrap();
    let mut pdf = |id, x, q2| lhapdf.xfx_q2(id, x, q2);
    let mut alphas = |q2| lhapdf.alphas_q2(q2);
    let mut cache = LumiCache::with_one(pdf_pdg_id, &mut pdf, &mut alphas);

    grid.convolute(&mut cache, &orders, bins, lumis, &SCALES_VECTOR[0..scales])
}

pub fn convolute_subgrid(
    grid: &Grid,
    lhapdf: &Pdf,
    order: usize,
    bin: usize,
    lumi: usize,
) -> Array3<f64> {
    let initial_state_1 = grid.key_values().map_or(2212, |map| {
        map.get("initial_state_1").unwrap().parse::<i32>().unwrap()
    });
    let initial_state_2 = grid.key_values().map_or(2212, |map| {
        map.get("initial_state_2").unwrap().parse::<i32>().unwrap()
    });

    // if the field 'Particle' is missing we assume it's a proton PDF
    let pdf_pdg_id = lhapdf
        .set()
        .entry("Particle")
        .unwrap_or_else(|| "2212".to_string())
        .parse::<i32>()
        .unwrap();

    let pdf = |id, x, q2| lhapdf.xfx_q2(id, x, q2);
    let anti_pdf = |id, x, q2| {
        let id = match id {
            -6..=6 | 11 | 13 | -11 | -13 => -id,
            21 | 22 => id,
            _ => unimplemented!(),
        };
        lhapdf.xfx_q2(id, x, q2)
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
    let alphas = |q2| lhapdf.alphas_q2(q2);

    grid.convolute_subgrid(&xfx1, &xfx2, &alphas, order, bin, lumi, 1.0, 1.0)
}
