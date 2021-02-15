use lhapdf::Pdf;
use pineappl::grid::Grid;
use prettytable::format::{FormatBuilder, LinePosition, LineSeparator};
use prettytable::Table;

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

pub fn convolute(
    grid: &Grid,
    lhapdf: &Pdf,
    orders: &[bool],
    bins: &[usize],
    lumis: &[bool],
    scales: &[(f64, f64)],
) -> Vec<f64> {
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

    grid.convolute(&xfx1, &xfx2, &alphas, orders, bins, lumis, scales)
}
