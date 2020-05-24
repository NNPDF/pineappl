#[macro_use]
extern crate clap;

use pineappl::grid::Grid;
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter};

fn merge(output: &str, input0: &str, input_rest: &[&str]) -> Result<(), Box<dyn Error>> {
    let output = File::create(output)?;
    let input0 = File::open(input0)?;
    let input_rest = input_rest
        .iter()
        .map(File::open)
        .collect::<Result<Vec<_>, std::io::Error>>()?;

    let mut grid0 = Grid::read(BufReader::new(input0))?;

    for i in input_rest {
        grid0.merge(Grid::read(BufReader::new(i))?)?;
    }

    grid0.write(BufWriter::new(output))?;

    Ok(())
}

fn main() {
    let matches = clap_app!(pineappl_cli =>
        (version: crate_version!())
        (author: crate_authors!())
        (about: crate_description!())
        (@subcommand merge =>
            (about: "Merges one or more PineAPPL grids together")
            (@arg output: +required "Path of the merged PineAPPL file")
            (@arg input: ... +required "Path(s) of the files that should be merged")
        )
    )
    .get_matches();

    if let Some(matches) = matches.subcommand_matches("merge") {
        let output = matches.value_of("output").unwrap();
        let input: Vec<_> = matches.values_of("input").unwrap().collect();

        merge(output, input.first().unwrap(), &input[1..]).expect("Merging failed");
    }
}
