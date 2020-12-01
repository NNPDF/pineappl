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
    pdf: &Pdf,
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
    // TODO: make sure this is a proton PDF
    let lhapdf_pdf = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let no_pdf = |_, x, _| x;
    let xfx1: Box<dyn Fn(i32, f64, f64) -> f64> = match initial_state_1 {
        2212 => Box::new(&lhapdf_pdf),
        // TODO: handle anti-proton PDF
        11 | 13 | -11 | -13 => Box::new(&no_pdf),
        _ => unimplemented!(),
    };
    let xfx2: Box<dyn Fn(i32, f64, f64) -> f64> = match initial_state_2 {
        2212 => Box::new(&lhapdf_pdf),
        // TODO: handle anti-proton PDF
        11 | 13 | -11 | -13 => Box::new(&no_pdf),
        _ => unimplemented!(),
    };
    let alphas = |q2| pdf.alphas_q2(q2);

    grid.convolute(&xfx1, &xfx2, &alphas, &orders, &bins, &lumis, &scales)
}
