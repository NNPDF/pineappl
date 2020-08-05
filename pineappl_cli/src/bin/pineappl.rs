#[macro_use]
extern crate clap;

use itertools::Itertools;
use lhapdf::{Pdf, PdfSet};
use pineappl::grid::Grid;
use prettytable::format::{FormatBuilder, LinePosition, LineSeparator};
use prettytable::{cell, row, Row, Table};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::HashSet;
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter};

fn create_table() -> Table {
    let mut table = Table::new();
    table.set_format(
        FormatBuilder::new()
            .column_separator(' ')
            .separator(LinePosition::Title, LineSeparator::new('-', '+', ' ', ' '))
            .build(),
    );
    table
}

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

    let matches: Vec<_> = order.match_indices('a').collect();

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
    absolute: bool,
) -> Result<Table, Box<dyn Error>> {
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
    let bin_limits = grid.bin_limits().limits();

    let mut table = create_table();
    let mut title = Row::empty();
    title.add_cell(cell!(c->"bin"));
    title.add_cell(cell!(c->"xmin"));
    title.add_cell(cell!(c->"xmax"));
    title.add_cell(cell!(c->"diff"));
    title.add_cell(cell!(c->"integ"));

    if absolute {
        for scale in &scales_vector[0..scales] {
            title.add_cell(cell!(c->&format!("({},{})", scale.0, scale.1)));
        }
    } else {
        title.add_cell(cell!(c->"neg unc"));
        title.add_cell(cell!(c->"pos unc"));
    }

    for other in other_pdfsets.iter() {
        let mut cell = cell!(c->other);
        cell.set_hspan(2);
        title.add_cell(cell);
    }

    table.set_titles(title);

    for (bin, values) in results.chunks_exact(scales).enumerate() {
        let min_value = values
            .iter()
            .min_by(|left, right| left.partial_cmp(right).unwrap())
            .unwrap();
        let max_value = values
            .iter()
            .max_by(|left, right| left.partial_cmp(right).unwrap())
            .unwrap();

        let row = table.add_empty_row();

        row.add_cell(cell!(r->&format!("{}", show_bins[bin])));
        row.add_cell(cell!(r->&format!("{}", bin_limits[show_bins[bin]])));
        row.add_cell(cell!(r->&format!("{}", bin_limits[show_bins[bin] + 1])));
        row.add_cell(cell!(r->&format!("{:.7e}", values[0])));
        row.add_cell(cell!(r->&format!("{:.7e}", values[0] * bin_sizes[show_bins[bin]])));

        if absolute {
            for value in values.iter() {
                row.add_cell(cell!(r->&format!("{:.7e}", value * bin_sizes[show_bins[bin]])));
            }
        } else {
            row.add_cell(cell!(r->&format!("{:.2}%", (min_value / values[0] - 1.0) * 100.0)));
            row.add_cell(cell!(r->&format!("{:.2}%", (max_value / values[0] - 1.0) * 100.0)));
        }

        for other in other_results.iter().skip(bin).step_by(show_bins.len()) {
            row.add_cell(cell!(r->&format!("{:.7e}", other)));
            row.add_cell(cell!(r->&format!("{:.2}%", (other / values[0] - 1.0) * 100.0)));
        }
    }

    Ok(table)
}

fn orders(input: &str, pdfset: &str, absolute: bool) -> Result<Table, Box<dyn Error>> {
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

    let mut sorted_grid_orders: Vec<_> = grid_orders
        .iter()
        .filter(|order| (order.logxir == 0) && (order.logxif == 0))
        .collect();
    sorted_grid_orders.sort();

    let unsorted_indices: Vec<_> = sorted_grid_orders
        .iter()
        .map(|sorted| {
            grid_orders
                .iter()
                .position(|unsorted| unsorted == *sorted)
                .unwrap()
        })
        .collect();
    let lo_power = {
        let order = sorted_grid_orders.first().unwrap();
        order.alphas + order.alpha
    };

    let bin_limits = grid.bin_limits().limits();

    let mut table = create_table();
    let mut title = Row::empty();
    title.add_cell(cell!(c->"bin"));
    title.add_cell(cell!(c->"xmin"));
    title.add_cell(cell!(c->"xmax"));
    title.add_cell(cell!(c->"diff"));

    for order in sorted_grid_orders.iter() {
        title.add_cell(cell!(c->&format!("O(as^{} a^{})", order.alphas, order.alpha)));
    }

    table.set_titles(title);

    for (bin, value) in results.iter().enumerate() {
        let row = table.add_empty_row();

        row.add_cell(cell!(r->&format!("{}", bin)));
        row.add_cell(cell!(r->&format!("{}", bin_limits[bin])));
        row.add_cell(cell!(r->&format!("{}", bin_limits[bin + 1])));
        row.add_cell(cell!(r->&format!("{:.7e}", value)));

        let mut leading_order = 0.0;

        // calculate the sum of all leading orders
        for (index, order) in sorted_grid_orders.iter().enumerate() {
            if (order.alphas + order.alpha) == lo_power {
                leading_order += order_results[unsorted_indices[index]][bin];
            }
        }

        // print each order normalized to the sum of all leading orders
        for index in 0..sorted_grid_orders.len() {
            let result = order_results[unsorted_indices[index]][bin];

            if absolute {
                row.add_cell(cell!(r->&format!("{:.7e}", result)));
            } else {
                row.add_cell(cell!(r->&format!("{:.2}%", result / leading_order * 100.0)));
            }
        }
    }

    Ok(table)
}

fn diff(input1: &str, input2: &str, pdfset: &str) -> Result<Table, Box<dyn Error>> {
    let grid1 = Grid::read(BufReader::new(File::open(input1)?))?;
    let grid2 = Grid::read(BufReader::new(File::open(input2)?))?;
    let pdf = pdfset
        .parse()
        .map(Pdf::with_lhaid)
        .unwrap_or_else(|_| Pdf::with_setname_and_member(pdfset, 0));

    let mut table = create_table();

    if grid1.bin_limits() != grid2.bin_limits() {
        print!("--- Bin limits: ");
        for limit in grid1.bin_limits().limits() {
            print!("{} ", limit);
        }
        println!();
        print!("+++ Bin limits: ");
        for limit in grid2.bin_limits().limits() {
            print!("{} ", limit);
        }
        println!();
    } else {
        let orders1: HashSet<_> = grid1
            .orders()
            .iter()
            .filter(|order| (order.logxir == 0) && (order.logxif == 0))
            .collect();
        let orders2: HashSet<_> = grid2
            .orders()
            .iter()
            .filter(|order| (order.logxir == 0) && (order.logxif == 0))
            .collect();

        let mut diff1 = orders1.difference(&orders2).peekable();
        let mut diff2 = orders2.difference(&orders1).peekable();
        if diff1.peek().is_some() || diff2.peek().is_some() {
            print!("--- Orders: ");
            for order in diff1 {
                if order.alphas == 0 {
                    print!("O(a^{}) ", order.alpha);
                } else if order.alpha == 0 {
                    print!("O(as^{}) ", order.alphas);
                } else {
                    print!("O(as^{} a^{}) ", order.alphas, order.alpha);
                }
            }
            println!();
            print!("+++ Orders: ");
            for order in diff2 {
                if order.alphas == 0 {
                    print!("O(a^{}) ", order.alpha);
                } else if order.alpha == 0 {
                    print!("O(as^{}) ", order.alphas);
                } else {
                    print!("O(as^{} a^{}) ", order.alphas, order.alpha);
                }
            }
            println!();
            println!();
        }

        let mut title = Row::empty();
        title.add_cell(cell!(c->"bin"));
        title.add_cell(cell!(c->"xmin"));
        title.add_cell(cell!(c->"xmax"));

        for order in orders1.intersection(&orders2) {
            let mut cell = cell!(c->&format!("O(as^{} a^{})", order.alphas, order.alpha));
            cell.set_hspan(3);
            title.add_cell(cell);
        }

        table.set_titles(title);

        let order_results1: Vec<Vec<f64>> = orders1
            .intersection(&orders2)
            .map(|order| {
                let mut order_mask = vec![false; grid1.orders().len()];
                order_mask[grid1.orders().iter().position(|o| o == *order).unwrap()] = true;
                grid1.convolute(
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
        let order_results2: Vec<Vec<f64>> = orders1
            .intersection(&orders2)
            .map(|order| {
                let mut order_mask = vec![false; grid2.orders().len()];
                order_mask[grid2.orders().iter().position(|o| o == *order).unwrap()] = true;
                grid2.convolute(
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

        for (bin, limits) in grid1.bin_limits().limits().windows(2).enumerate() {
            let row = table.add_empty_row();

            row.add_cell(cell!(r->bin));
            row.add_cell(cell!(r->limits[0]));
            row.add_cell(cell!(r->limits[1]));

            for (result1, result2) in order_results1.iter().zip(order_results2.iter()) {
                let result1 = result1[bin];
                let result2 = result2[bin];
                row.add_cell(cell!(r->&format!("{:.3e}", result1)));
                row.add_cell(cell!(r->&format!("{:.3e}", result2)));
                row.add_cell(cell!(r->&format!("{:.2}%",
                    if result1 == result2 { 0.0 } else { result1 / result2 - 1.0 })));
            }
        }
    }

    Ok(table)
}

fn info(input: &str, mode: &str) -> Result<(), Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;

    let mut sorted_grid_orders: Vec<_> = grid
        .orders()
        .iter()
        .filter(|order| (order.logxir == 0) && (order.logxif == 0))
        .collect();
    sorted_grid_orders.sort();

    let orders = sorted_grid_orders
        .into_iter()
        .group_by(|order| order.alphas + order.alpha)
        .into_iter()
        .map(|mut iter| match mode {
            "qcd" => iter.1.next().unwrap(),
            "ew" => iter.1.last().unwrap(),
            _ => unreachable!(),
        })
        .map(|order| {
            if order.alphas == 0 {
                format!("a{}", order.alpha)
            } else if order.alpha == 0 {
                format!("as{}", order.alphas)
            } else {
                format!("as{}a{}", order.alphas, order.alpha)
            }
        })
        .collect::<Vec<_>>()
        .join(",");

    println!("{}", orders);

    Ok(())
}

fn luminosity(input: &str) -> Result<Table, Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;

    let mut table = create_table();
    table.set_titles(row![c => "id", "entry"]);

    for (index, entry) in grid.lumi().iter().enumerate() {
        let row = table.add_empty_row();

        row.add_cell(cell!(&format!("{}", index)));

        for (id1, id2, factor) in entry.entry().iter() {
            row.add_cell(cell!(&format!("{} × ({:2.}, {:2.})", factor, id1, id2)));
        }
    }

    Ok(table)
}

fn channels(input: &str, pdfset: &str, limit: usize) -> Result<Table, Box<dyn Error>> {
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

    let bin_limits = grid.bin_limits().limits();

    let mut table = create_table();
    table.set_titles(row![c => "bin", "xmin", "xmax", "lumi", "size"]);

    // TODO: add more titles

    for bin in 0..grid.bin_limits().bins() {
        let row = table.add_empty_row();

        row.add_cell(cell!(r->&format!("{}", bin)));
        row.add_cell(cell!(r->&format!("{}", bin_limits[bin])));
        row.add_cell(cell!(r->&format!("{}", bin_limits[bin + 1])));

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
            row.add_cell(cell!(r->&format!("#{}", lumi)));
            row.add_cell(cell!(r->&format!("{:.2}%", percentage)));
        }
    }

    Ok(table)
}

fn pdf_uncertainty(
    input: &str,
    pdfset: &str,
    cl: f64,
    threads: usize,
    orders: &[(u32, u32)],
) -> Result<Table, Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;
    let set = PdfSet::new(pdfset);
    let pdfs = set.mk_pdfs();

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

    ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();

    let results: Vec<f64> = pdfs
        .into_par_iter()
        .flat_map(|pdf| {
            grid.convolute(
                &|id, x, q2| pdf.xfx_q2(id, x, q2),
                &|id, x, q2| pdf.xfx_q2(id, x, q2),
                &|q2| pdf.alphas_q2(q2),
                &orders,
                &[],
                &[],
                &[(1.0, 1.0)],
            )
        })
        .collect();

    let bin_sizes = grid.bin_limits().bin_sizes();
    let bin_limits = grid.bin_limits().limits();

    let mut table = create_table();
    table.set_titles(row![c => "bin", "xmin", "xmax", "diff", "integ", "neg unc", "pos unc"]);

    for bin in 0..bin_sizes.len() {
        let values: Vec<_> = results
            .iter()
            .skip(bin)
            .step_by(bin_sizes.len())
            .cloned()
            .collect();
        let uncertainty = set.uncertainty(&values, cl, false);

        table.add_row(row![r =>
            &format!("{}", bin),
            &format!("{}", bin_limits[bin]),
            &format!("{}", bin_limits[bin + 1]),
            &format!("{:.7e}", uncertainty.central),
            &format!("{:.7e}", uncertainty.central * bin_sizes[bin]),
            &format!("{:.2}%", (-uncertainty.errminus / uncertainty.central) * 100.0),
            &format!("{:.2}%", (uncertainty.errplus / uncertainty.central) * 100.0),
        ]);
    }

    Ok(table)
}

fn main() -> Result<(), Box<dyn Error>> {
    let num_cpus = num_cpus::get().to_string();
    let matches = clap_app!(pineappl =>
        (author: crate_authors!())
        (about: crate_description!())
        (version: crate_version!())
        (@setting DisableHelpSubcommand)
        (@setting SubcommandRequiredElseHelp)
        (@setting VersionlessSubcommands)
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
            (@arg absolute: -a --absolute "Show absolute numbers of the scale variation")
        )
        (@subcommand diff =>
            (about: "Compares the contents of two grids with each other")
            (@arg input1: +required "Path to the first grid")
            (@arg input2: +required "Path to the second grid")
            (@arg pdfset: +required validator(validate_pdfset)
                "LHAPDF id(s) or name of the PDF set(s)")
        )
        (@subcommand info =>
            (about: "Shows information about the grid")
            (@arg input: +required "Path to the input grid")
            (@group mode +required =>
                (@arg qcd: --qcd "For each order print a list of the largest QCD order")
                (@arg ew: --ew "For each order print a list of the largest EW order")
            )
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
            (about: "Shows the predictions for all bin for each order separately")
            (@arg input: +required "Path to the input grid")
            (@arg pdfset: +required validator(validate_pdfset) "LHAPDF id or name of the PDF set")
            (@arg absolute: -a --absolute "Show absolute numbers of each perturbative order")
        )
        (@subcommand pdf_uncertainty =>
            (about: "Calculates PDF uncertainties")
            (@arg cl: --cl default_value("68.268949213708581") "Confidence level in per cent")
            (@arg input: +required "Path to the input grid")
            (@arg pdfset: +required validator(validate_pdfset) "LHAPDF id or name of the PDF set")
            (@arg threads: --threads default_value(&num_cpus) "Number of threads to utilize")
            (@arg orders: -o --orders +use_delimiter min_values(1) "Select orders manually")
        )
    )
    .get_matches();

    if let Some(matches) = matches.subcommand_matches("channels") {
        let input = matches.value_of("input").unwrap();
        let pdfset = matches.value_of("pdfset").unwrap();
        let limit = matches.value_of("limit").unwrap().parse()?;

        channels(input, pdfset, limit)?.printstd();
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
        let absolute = matches.is_present("absolute");

        convolute(
            input,
            pdfset.first().unwrap(),
            &pdfset[1..],
            &bins,
            scales,
            &orders,
            absolute,
        )?
        .printstd();
    } else if let Some(matches) = matches.subcommand_matches("diff") {
        let input1 = matches.value_of("input1").unwrap();
        let input2 = matches.value_of("input2").unwrap();
        let pdfset = matches.value_of("pdfset").unwrap();

        diff(input1, input2, pdfset)?.printstd();
    } else if let Some(matches) = matches.subcommand_matches("info") {
        let input = matches.value_of("input").unwrap();

        info(
            input,
            if matches.is_present("qcd") {
                "qcd"
            } else if matches.is_present("ew") {
                "ew"
            } else {
                unreachable!()
            },
        )?;
    } else if let Some(matches) = matches.subcommand_matches("luminosity") {
        let input = matches.value_of("input").unwrap();

        luminosity(input)?.printstd();
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
        let absolute = matches.is_present("absolute");

        orders(input, pdfset, absolute)?.printstd();
    } else if let Some(matches) = matches.subcommand_matches("pdf_uncertainty") {
        let input = matches.value_of("input").unwrap();
        let pdfset = matches.value_of("pdfset").unwrap();
        let cl = matches.value_of("cl").unwrap().parse()?;
        let threads = matches.value_of("threads").unwrap().parse()?;
        let orders: Vec<_> = matches
            .values_of("orders")
            .map_or(vec![], |values| values.map(parse_order).collect())
            .into_iter()
            .collect::<Result<_, _>>()?;

        pdf_uncertainty(input, pdfset, cl, threads, &orders)?.printstd();
    }

    Ok(())
}
