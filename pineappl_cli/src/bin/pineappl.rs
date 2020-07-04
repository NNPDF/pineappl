#[macro_use]
extern crate clap;

use itertools::Itertools;
use lhapdf::{Pdf, PdfSet};
use pineappl::grid::Grid;
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter};
use std::rc::Rc;

fn validate_pos_non_zero(argument: String) -> Result<(), String> {
    if let Ok(number) = argument.parse::<usize>() {
        if number != 0 {
            return Ok(());
        }
    }

    Err(format!(
        "The value `{}` is not positive and non-zero",
        argument
    ))
}

fn validate_pdfset(argument: String) -> Result<(), String> {
    if let Ok(lhaid) = argument.parse() {
        if lhapdf::lookup_pdf(lhaid).is_some() {
            return Ok(());
        }

        return Err(format!(
            "The PDF set for the LHAPDF ID `{}` was not found",
            argument
        ));
    } else if lhapdf::available_pdf_sets()
        .iter()
        .any(|set| *set == argument)
    {
        return Ok(());
    }

    Err(format!("The PDF set `{}` was not found", argument))
}

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

fn parse_order(order: &str) -> Result<(u32, u32), Box<dyn Error>> {
    let mut alphas = 0;
    let mut alpha = 0;

    let matches: Vec<_> = order.match_indices("a").collect();

    if matches.len() > 2 {
        todo!();
    } else {
        for (index, _) in matches {
            if &order[index..index + 2] == "as" {
                let len = order[index + 2..]
                    .chars()
                    .take_while(|c| c.is_numeric())
                    .count();
                alphas = str::parse::<u32>(&order[index + 2..index + 2 + len])?;
            } else {
                let len = order[index + 1..]
                    .chars()
                    .take_while(|c| c.is_numeric())
                    .count();
                alpha = str::parse::<u32>(&order[index + 1..index + 1 + len])?;
            }
        }
    }

    Ok((alphas, alpha))
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
    scales: usize,
    orders: &[(u32, u32)],
) -> Result<(), Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;
    let show_bins = if show_bins.is_empty() {
        (0..grid.bin_limits().bins()).collect()
    } else {
        show_bins.to_vec()
    };
    let pdf = pdfset
        .parse()
        .map(Pdf::with_lhaid)
        .unwrap_or_else(|_| Pdf::with_setname_and_member(pdfset, 0));
    let scales_vector = vec![
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

    let results = grid.convolute(
        &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
        &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
        &|q2| pdf.alphas_q2(q2),
        &orders,
        &show_bins,
        &[],
        &scales_vector[0..scales],
    );

    let other_results: Vec<f64> = other_pdfsets
        .iter()
        .map(|pdfset| {
            let pdf = pdfset
                .parse()
                .map(Pdf::with_lhaid)
                .unwrap_or_else(|_| Pdf::with_setname_and_member(pdfset, 0));
            grid.convolute(
                &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
                &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
                &|q2| pdf.alphas_q2(q2),
                &[],
                &show_bins,
                &[],
                &[(1.0, 1.0)],
            )
        })
        .flatten()
        .collect();

    let bin_sizes = grid.bin_limits().bin_sizes();

    for (bin, values) in results.chunks_exact(scales).enumerate() {
        let min_value = values
            .iter()
            .min_by(|left, right| left.partial_cmp(right).unwrap())
            .unwrap();
        let max_value = values
            .iter()
            .max_by(|left, right| left.partial_cmp(right).unwrap())
            .unwrap();

        print!(
            "{:<3}  {:>12.7e}  {:>12.7e}  {:+5.2}% {:+5.2}%",
            show_bins[bin],
            values[0],
            values[0] * bin_sizes[show_bins[bin]],
            (min_value / values[0] - 1.0) * 100.0,
            (max_value / values[0] - 1.0) * 100.0,
        );

        for other in other_results.iter().skip(bin).step_by(show_bins.len()) {
            print!("  {:+6.2}%", (other / values[0] - 1.0) * 100.0);
        }

        println!();
    }

    Ok(())
}

fn orders(input: &str, pdfset: &str) -> Result<(), Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;
    let pdf = pdfset
        .parse()
        .map(Pdf::with_lhaid)
        .unwrap_or_else(|_| Pdf::with_setname_and_member(pdfset, 0));

    let grid_orders = grid.orders();
    let results = grid.convolute(
        &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
        &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
        &|q2| pdf.alphas_q2(q2),
        &[],
        &[],
        &[],
        &[(1.0, 1.0)],
    );

    let order_results: Vec<Vec<f64>> = (0..grid_orders.len())
        .map(|order| {
            let mut order_mask = vec![false; grid_orders.len()];
            order_mask[order] = true;
            grid.convolute(
                &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
                &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
                &|q2| pdf.alphas_q2(q2),
                &order_mask,
                &[],
                &[],
                &[(1.0, 1.0)],
            )
        })
        .collect();

    for (index, order) in grid_orders.iter().enumerate() {
        if (order.logxir != 0) || (order.logxif != 0) {
            continue;
        }

        println!("{:<2}: O(as^{} a^{})", index, order.alphas, order.alpha);
    }

    println!();

    for (bin, value) in results.iter().enumerate() {
        print!("{:<3}  {:>12.7e}", bin, value);

        for index in 0..grid_orders.len() {
            if (grid_orders[index].logxir != 0) || (grid_orders[index].logxif != 0) {
                continue;
            }

            let order_value = order_results[index][bin];

            print!(" {:>6.2}%", order_value / order_results[0][bin] * 100.0);
        }

        println!();
    }

    Ok(())
}

fn luminosity(input: &str) -> Result<(), Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;

    for (index, entry) in grid.lumi().iter().enumerate() {
        print!("{:>2}: ", index);
        println!(
            "{}",
            entry
                .entry()
                .iter()
                .map(|(id1, id2, factor)| format!("{} × ({}, {})", factor, id1, id2))
                .join(" + ")
        );
    }

    Ok(())
}

fn channels(input: &str, pdfset: &str, limit: usize) -> Result<(), Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;
    let pdf = pdfset
        .parse()
        .map(Pdf::with_lhaid)
        .unwrap_or_else(|_| Pdf::with_setname_and_member(pdfset, 0));

    let results: Vec<_> = (0..grid.lumi().len())
        .map(|lumi| {
            let mut lumi_mask = vec![false; grid.lumi().len()];
            lumi_mask[lumi] = true;
            grid.convolute(
                &|id, x1, q2| pdf.xfx_q2(id, x1, q2),
                &|id, x2, q2| pdf.xfx_q2(id, x2, q2),
                &|q2| pdf.alphas_q2(q2),
                &[],
                &[],
                &lumi_mask,
                &[(1.0, 1.0)],
            )
        })
        .collect();

    for bin in 0..grid.bin_limits().bins() {
        print!("{:>3}:", bin);

        let sum: f64 = results.iter().map(|vec| vec[bin]).sum();
        let mut percentages: Vec<_> = results
            .iter()
            .enumerate()
            .map(|(lumi, vec)| (lumi, vec[bin] / sum * 100.0))
            .collect();

        // sort using the absolute value in descending order
        percentages.sort_unstable_by(|(_, left), (_, right)| {
            right.abs().partial_cmp(&left.abs()).unwrap()
        });

        for (lumi, percentage) in percentages.iter().take(limit) {
            print!("  #{:<3} {:5.1}%", lumi, percentage);
        }
        println!();
    }

    Ok(())
}

fn pdf_uncertainty(input: &str, pdfset: &str, cl: f64) -> Result<(), Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;
    let set = PdfSet::new(pdfset);
    let pdfs: Vec<_> = set.mk_pdfs().into_iter().map(|pdf| Rc::new(pdf)).collect();

    let xfx_array: Vec<_> = pdfs
        .iter()
        .map(|pdf| pdf.clone())
        .map(|pdf| {
            Box::new(move |id, x, q2| pdf.xfx_q2(id, x, q2)) as Box<dyn Fn(i32, f64, f64) -> f64>
        })
        .collect();
    let alphas_array: Vec<_> = pdfs
        .iter()
        .map(|pdf| pdf.clone())
        .map(|pdf| Box::new(move |q2| pdf.alphas_q2(q2)) as Box<dyn Fn(f64) -> f64>)
        .collect();

    let results = grid.convolute_multiple(
        &xfx_array,
        &xfx_array,
        &alphas_array,
        &[],
        &[],
        &[],
        &[(1.0, 1.0)],
    );

    for (bin, values) in results.iter().enumerate() {
        let uncertainty = set.uncertainty(values, cl, false);

        println!(
            "{:<3}  {:>12.7e}  {:>12.7e}  {:+5.2}% {:+5.2}%",
            bin,
            uncertainty.central,
            uncertainty.central,
            (-uncertainty.errminus / uncertainty.central) * 100.0,
            (uncertainty.errplus / uncertainty.central) * 100.0,
        );
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = clap_app!(pineappl =>
        (author: crate_authors!())
        (about: crate_description!())
        (version: crate_version!())
        (@setting SubcommandRequiredElseHelp)
        (@subcommand channels =>
            (about: "Shows the contribution for each partonic channel")
            (@arg input: +required "Path to the input grid")
            (@arg pdfset: +required validator(validate_pdfset) "LHAPDF id or name of the PDF set")
            (@arg limit: -l --limit default_value("10") validator(validate_pos_non_zero)
                "The maximum number of channels displayed")
        )
        (@subcommand convolute =>
            (about: "Convolutes a PineAPPL grid with a PDF set")
            (@arg input: +required "Path of the input grid")
            (@arg pdfset: ... +required validator(validate_pdfset)
                "LHAPDF id(s) or name of the PDF set(s)")
            (@arg bins: -b --bins +takes_value "Selects a subset of bins")
            (@arg scales: -s --scales default_value("7") possible_values(&["1", "3", "7", "9"])
                "Set the number of scale variations")
            (@arg orders: -o --orders +use_delimiter min_values(1) "Select orders manually")
        )
        (@subcommand luminosity =>
            (about: "Shows the luminosity function")
            (@arg input: +required "Path to the input grid")
        )
        (@subcommand merge =>
            (about: "Merges one or more PineAPPL grids together")
            (@arg output: +required "Path of the merged PineAPPL file")
            (@arg input: ... +required "Path(s) of the files that should be merged")
            (@arg scale: -s --scale +takes_value "Scales all grids with the given factor")
            (@arg scale_by_order: --scale_by_order +takes_value conflicts_with[scale]
                number_of_values(5) "Scales all grids with order-dependent factors")
        )
        (@subcommand orders =>
            (about: "Shows thw predictions for all bin for each order separately")
            (@arg input: +required "Path to the input grid")
            (@arg pdfset: +required validator(validate_pdfset) "LHAPDF id or name of the PDF set")
        )
        (@subcommand pdf_uncertainty =>
            (about: "Calculates PDF uncertainties")
            (@arg cl: --cl default_value("68.268949213708581") "Confidence level in per cent")
            (@arg input: +required "Path to the input grid")
            (@arg pdfset: +required validator(validate_pdfset) "LHAPDF id or name of the PDF set")
        )
    )
    .get_matches();

    if let Some(matches) = matches.subcommand_matches("channels") {
        let input = matches.value_of("input").unwrap();
        let pdfset = matches.value_of("pdfset").unwrap();
        let limit = matches.value_of("limit").unwrap().parse()?;

        return channels(input, pdfset, limit);
    } else if let Some(matches) = matches.subcommand_matches("convolute") {
        let input = matches.value_of("input").unwrap();
        let pdfset: Vec<_> = matches.values_of("pdfset").unwrap().collect();
        let bins = parse_integer_list(matches.value_of("bins").unwrap_or(""))?;
        let scales = matches.value_of("scales").unwrap().parse()?;
        let orders: Vec<_> = matches
            .values_of("orders")
            .map_or(vec![], |values| values.map(parse_order).collect())
            .into_iter()
            .collect::<Result<_, _>>()?;

        return convolute(
            input,
            pdfset.first().unwrap(),
            &pdfset[1..],
            &bins,
            scales,
            &orders,
        );
    } else if let Some(matches) = matches.subcommand_matches("luminosity") {
        let input = matches.value_of("input").unwrap();

        return luminosity(input);
    } else if let Some(matches) = matches.subcommand_matches("merge") {
        let output = matches.value_of("output").unwrap();
        let input: Vec<_> = matches.values_of("input").unwrap().collect();
        let scale = matches
            .value_of("scale")
            .map(str::parse::<f64>)
            .transpose()?;
        let scale_by_order: Vec<_> = matches
            .values_of("scale_by_order")
            .map_or(vec![], |s| s.map(str::parse::<f64>).collect())
            .into_iter()
            .collect::<Result<_, _>>()?;

        return merge(
            output,
            input.first().unwrap(),
            &input[1..],
            scale,
            &scale_by_order,
        );
    } else if let Some(matches) = matches.subcommand_matches("orders") {
        let input = matches.value_of("input").unwrap();
        let pdfset = matches.value_of("pdfset").unwrap();

        return orders(input, pdfset);
    } else if let Some(matches) = matches.subcommand_matches("pdf_uncertainty") {
        let input = matches.value_of("input").unwrap();
        let pdfset = matches.value_of("pdfset").unwrap();
        let cl = matches.value_of("cl").unwrap().parse()?;

        return pdf_uncertainty(input, pdfset, cl);
    }

    Ok(())
}
