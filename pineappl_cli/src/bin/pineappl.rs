#[macro_use]
extern crate clap;

use pineappl::grid::Grid;
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter};

fn merge(
    output: &str,
    input0: &str,
    input_rest: &[&str],
    scale: Option<f64>,
    scale_by_order: &[f64],
) -> Result<(), Box<dyn Error>> {
    let output = OpenOptions::new().write(true).create_new(true).open(output)?;
    let input0 = File::open(input0)?;
    let input_rest = input_rest
        .iter()
        .map(File::open)
        .collect::<Result<Vec<_>, std::io::Error>>()?;

    let mut grid0 = Grid::read(BufReader::new(input0))?;

    for i in input_rest {
        grid0.merge(Grid::read(BufReader::new(i))?)?;
    }

    if let Some(scale) = scale {
        grid0.scale(scale);
    } else if !scale_by_order.is_empty() {
        grid0.scale_by_order(
            scale_by_order[0],
            scale_by_order[1],
            scale_by_order[2],
            scale_by_order[3],
            scale_by_order[4],
        );
    }

    grid0.write(BufWriter::new(output))?;

    Ok(())
}

fn main() {
    let matches = clap_app!(pineappl_cli =>
        (version: crate_version!())
        (author: crate_authors!())
        (about: crate_description!())
        (@setting SubcommandRequiredElseHelp)
        (@subcommand merge =>
            (about: "Merges one or more PineAPPL grids together")
            (@arg output: +required "Path of the merged PineAPPL file")
            (@arg input: ... +required "Path(s) of the files that should be merged")
            (@arg scale: -s --scale +takes_value "Scales all grids with the given factor")
            (@arg scale_by_order: -b --scale_by_order +takes_value conflicts_with[scale]
                number_of_values(5) "Scales all grids with order-dependent factors")
        )
    )
    .get_matches();

    if let Some(matches) = matches.subcommand_matches("merge") {
        let str_to_f64 =
            |s| str::parse::<f64>(s).expect("Could not convert string to floating point number");

        let output = matches.value_of("output").unwrap();
        let input: Vec<_> = matches.values_of("input").unwrap().collect();
        let scale = matches.value_of("scale").map(str_to_f64);
        let scale_by_order: Vec<_> = matches
            .values_of("scale_by_order")
            .map(|s| s.map(str_to_f64).collect())
            .unwrap_or(Vec::new());

        merge(
            output,
            input.first().unwrap(),
            &input[1..],
            scale,
            &scale_by_order,
        )
        .expect("Merging failed");
    }
}
