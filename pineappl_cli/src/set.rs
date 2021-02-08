use pineappl::grid::Grid;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter};

pub fn subcommand(
    input: &str,
    output: &str,
    entries: Vec<&str>,
    entries_from_file: Vec<&str>,
    deletes: Vec<&str>,
) -> Result<(), Box<dyn Error>> {
    let mut grid = Grid::read(BufReader::new(File::open(input)?))?;

    for key_value in entries.chunks(2) {
        grid.set_key_value(key_value[0], key_value[1]);
    }

    for key_file in entries_from_file.chunks(2) {
        grid.set_key_value(key_file[0], &fs::read_to_string(key_file[1])?);
    }

    for delete in &deletes {
        grid.key_values_mut().remove(*delete);
    }

    grid.write(BufWriter::new(File::create(output)?))?;

    Ok(())
}
