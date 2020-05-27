#[macro_use]
extern crate clap;

use lhapdf::Pdf;
use pineappl::grid::Grid;
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter};

fn parse_integer_list(list: &str) -> Result<Vec<usize>, Box<dyn Error>> {
    let mut integers = Vec::new();

    for s in list.split_terminator(',') {
        if let Some(at) = s.find('-') {
            let (left, right) = s.split_at(at);
            integers.extend(str::parse::<usize>(left)?..=str::parse::<usize>(&right[1..])?);
        } else {
            integers.push(str::parse::<usize>(s)?);
        }
    }

    Ok(integers)
}

fn merge(
    output: &str,
    input0: &str,
    input_rest: &[&str],
    scale: Option<f64>,
    scale_by_order: &[f64],
) -> Result<(), Box<dyn Error>> {
    let output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
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

fn convolute(
    input: &str,
    pdfset: &str,
    other_pdfsets: &[&str],
    show_bins: &[usize],
) -> Result<(), Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;
    let show_bins = if show_bins.is_empty() {
        (0..grid.bin_limits().bins()).collect()
    } else {
        show_bins.to_vec()
    };
    let pdf = Pdf::new(pdfset, 0);

    let results = grid.convolute(
        &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
        &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
        &|q2| pdf.alphas_q2(q2),
        &[],
        &show_bins,
        &[],
        &(1.0, 1.0),
    );

    let other_results: Vec<f64> = other_pdfsets
        .iter()
        .map(|pdfset| {
            let pdf = Pdf::new(pdfset, 0);
            grid.convolute(
                &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
                &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
                &|q2| pdf.alphas_q2(q2),
                &[],
                &show_bins,
                &[],
                &(1.0, 1.0),
            )
        })
        .flatten()
        .collect();

    let bin_sizes = grid.bin_limits().bin_sizes();

    for (bin, value) in results.iter().enumerate() {
        print!(
            "{:<3}  {:>12.7e}  {:>12.7e}",
            show_bins[bin],
            value,
            value * bin_sizes[show_bins[bin]],
        );

        for other in other_results.iter().skip(bin).step_by(show_bins.len()) {
            print!("  {:+6.2}%", (other / value - 1.0) * 100.0);
        }

        println!("");
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = clap_app!(pineappl =>
        (version: crate_version!())
        (author: crate_authors!())
        (about: crate_description!())
        (@setting SubcommandRequiredElseHelp)
        (@subcommand merge =>
            (about: "Merges one or more PineAPPL grids together")
            (@arg output: +required "Path of the merged PineAPPL file")
            (@arg input: ... +required "Path(s) of the files that should be merged")
            (@arg scale: -s --scale +takes_value "Scales all grids with the given factor")
            (@arg scale_by_order: --scale_by_order +takes_value conflicts_with[scale]
                number_of_values(5) "Scales all grids with order-dependent factors")
        )
        (@subcommand convolute =>
            (about: "Convolutes a PineAPPL grid with a PDF set")
            (@arg input: +required "Path of the input grid")
            (@arg pdfset: ... +required "LHAPDF id(s) or name(s) of the PDF set(s)")
            (@arg bins: -b --bins +takes_value "Selects a subset of bins")
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

        return merge(
            output,
            input.first().unwrap(),
            &input[1..],
            scale,
            &scale_by_order,
        );
    } else if let Some(matches) = matches.subcommand_matches("convolute") {
        let input = matches.value_of("input").unwrap();
        let pdfset: Vec<_> = matches.values_of("pdfset").unwrap().collect();
        let mut bins = parse_integer_list(matches.value_of("bins").unwrap_or(""))?;

        return convolute(input, pdfset.first().unwrap(), &pdfset[1..], &bins);
    }

    Ok(())
}
