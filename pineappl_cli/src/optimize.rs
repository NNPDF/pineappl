use pineappl::grid::Grid;
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter};

pub fn subcommand(input: &str, output: &str) -> Result<(), Box<dyn Error>> {
    let input = File::open(input)?;
    let mut grid = Grid::read(BufReader::new(input))?;
    grid.optimize();
    let output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
    grid.write(BufWriter::new(output))?;

    Ok(())
}
