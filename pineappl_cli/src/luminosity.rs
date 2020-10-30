use pineappl::grid::Grid;
use prettytable::{cell, row, Table};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;

use super::helpers::create_table;

pub fn subcommand(input: &str) -> Result<Table, Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;

    let mut table = create_table();
    let mut titles = row![c => "id"];
    for _ in 0..grid.lumi().iter().map(|lumi| lumi.entry().len()).max().unwrap() {
        titles.add_cell(cell!(c->"entry"));
    }
    table.set_titles(titles);

    for (index, entry) in grid.lumi().iter().enumerate() {
        let row = table.add_empty_row();

        row.add_cell(cell!(&format!("{}", index)));

        for (id1, id2, factor) in entry.entry().iter() {
            row.add_cell(cell!(&format!(
                "{} \u{d7} ({:2.}, {:2.})",
                factor, id1, id2
            )));
        }
    }

    Ok(table)
}
