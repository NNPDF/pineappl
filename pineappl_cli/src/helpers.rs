use anyhow::{Context, Result};
use lhapdf::Pdf;
use ndarray::Array3;
use pineappl::grid::Grid;
use pineappl::lumi::LumiCache;
use prettytable::format::{FormatBuilder, LinePosition, LineSeparator};
use prettytable::Table;
use std::fs::{File, OpenOptions};
use std::path::Path;

pub const ONE_SIGMA: f64 = 68.268_949_213_708_58;

pub fn read_grid(input: &str) -> Result<Grid> {
    Grid::read(File::open(input).context(format!("unable to open '{}'", input))?)
        .context(format!("unable to read '{}'", input))
}

pub fn write_grid(output: &str, grid: &Grid) -> Result<()> {
    let path = Path::new(output);
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)
        .context(format!("unable to write '{}'", output))?;

    if path.extension().map_or(false, |ext| ext == "lz4") {
        Ok(grid.write_lz4(file)?)
    } else {
        Ok(grid.write(file)?)
    }
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
    // if the field 'Particle' is missing we assume it's a proton PDF
    let pdf_pdg_id = lhapdf
        .set()
        .entry("Particle")
        .unwrap_or_else(|| "2212".to_string())
        .parse::<i32>()
        .unwrap();
    let mut pdf = |id, x, q2| lhapdf.xfx_q2(id, x, q2);
    let mut alphas = |q2| lhapdf.alphas_q2(q2);
    let mut cache = LumiCache::with_one(pdf_pdg_id, &mut pdf, &mut alphas);

    grid.convolute_subgrid(&mut cache, order, bin, lumi, 1.0, 1.0)
}
