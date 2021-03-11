#![warn(clippy::all, clippy::cargo, clippy::nursery, clippy::pedantic)]

mod channels;
mod convolute;
mod diff;
mod helpers;
mod info;
mod luminosity;
mod merge;
mod optimize;
mod orders;
mod pdf_uncertainty;
mod remap;
mod set;
mod subgrids;

use clap::{clap_app, crate_authors, crate_description, crate_version};
use std::error::Error;
use std::str::FromStr;

fn validate_pos_non_zero<T: Default + FromStr + PartialEq>(argument: &str) -> Result<(), String> {
    if let Ok(number) = argument.parse::<T>() {
        if number != T::default() {
            return Ok(());
        }
    }

    Err(format!(
        "The value `{}` is not positive and non-zero",
        argument
    ))
}

fn validate_pdfset(argument: &str) -> Result<(), String> {
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

fn main() -> Result<(), Box<dyn Error>> {
    let num_cpus = num_cpus::get().to_string();
    let matches = clap_app!(pineappl =>
        (author: crate_authors!())
        (about: crate_description!())
        (version: crate_version!())
        (@arg silence_lhapdf: alias("silence_lhapdf") long("silence-lhapdf")
            "Prevents LHAPDF from printing banners")
        (@setting DisableHelpSubcommand)
        (@setting SubcommandRequiredElseHelp)
        (@setting VersionlessSubcommands)
        (@subcommand channels =>
            (about: "Shows the contribution for each partonic channel")
            (@arg input: +required "Path to the input grid")
            (@arg pdfset: +required validator(validate_pdfset) "LHAPDF id or name of the PDF set")
            (@arg limit: -l --limit default_value("10") validator(validate_pos_non_zero::<usize>)
                "The maximum number of channels displayed")
            (@arg orders: -o --orders +use_delimiter min_values(1) "Select orders manually")
            (@arg absolute: -a --absolute "Show absolute numbers of each contribution")
            (@arg lumis: --lumis +takes_value conflicts_with("limit")
                "Show only the listed channels")
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
            (@arg ignore_orders: alias("ignore_orders") long("ignore-orders")
                "Sums over all orders")
        )
        (@subcommand info =>
            (about: "Shows information about the grid")
            (@arg input: +required "Path to the input grid")
            (@group mode +required =>
                (@arg qcd: --qcd "For each order print a list of the largest QCD order")
                (@arg ew: --ew "For each order print a list of the largest EW order")
                (@arg get: --get +takes_value value_name("key") "Gets an internal key-value pair")
                (@arg keys: --keys "Show all keys stored in the grid")
                (@arg show: --show "Shows all key-value pairs stored in the grid")
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
            (@arg scale_by_order: alias("scale_by_order") long("scale-by-order") +takes_value
                conflicts_with[scale] number_of_values(5) value_names(&["alphas", "alpha",
                "logxir", "logxif", "global"]) "Scales all grids with order-dependent factors")
        )
        (@subcommand optimize =>
            (about: "Optimizes the internal data structure to minimize memory usage")
            (@arg input: +required "Path to the input grid")
            (@arg output: +required "Path to the optimized PineAPPL file")
        )
        (@subcommand orders =>
            (about: "Shows the predictions for all bin for each order separately")
            (@arg input: +required "Path to the input grid")
            (@arg pdfset: +required validator(validate_pdfset) "LHAPDF id or name of the PDF set")
            (@arg absolute: -a --absolute "Show absolute numbers of each perturbative order")
            (@arg normalize: -n --normalize +use_delimiter min_values(1) conflicts_with("absolute")
                "Normalize contributions to the specified orders")
        )
        (@subcommand pdf_uncertainty =>
            (about: "Calculates PDF uncertainties")
            (@arg cl: --cl default_value("68.268949213708581") "Confidence level in per cent")
            (@arg input: +required "Path to the input grid")
            (@arg pdfset: +required validator(validate_pdfset) "LHAPDF id or name of the PDF set")
            (@arg threads: --threads default_value(&num_cpus) "Number of threads to utilize")
            (@arg orders: -o --orders +use_delimiter min_values(1) "Select orders manually")
        )
        (@subcommand remap =>
            (about: "Modifies the bin dimensions, widths and normalizations")
            (@arg input: +required "Path to the input grid")
            (@arg output: +required "Path of the modified PineAPPL file")
            (@arg remapping: +required "Remapping string")
            (@arg norm: --norm default_value("1.0") validator(validate_pos_non_zero::<f64>)
                "Normalization factor in addition to the given bin widths")
            (@arg ignore_obs_norm: alias("ignore_obs_norm") long("ignore-obs-norm") +use_delimiter
                "Ignore the given observables for differential normalization")
        )
        (@subcommand set =>
            (about: "Modifies the internal key-value storage")
            (@arg input: +required "Path to the input grid")
            (@arg output: +required "Path of the modified PineAPPL file")
            (@group mode +multiple =>
                (@arg entry: --entry +allow_hyphen_values number_of_values(2)
                    value_names(&["key", "value"]) +multiple "Sets an internal key-value pair")
                (@arg entry_from_file: alias("entry_from_file") long("entry-from-file")
                    number_of_values(2) value_names(&["key", "file"]) +multiple
                    "Sets an internal key-value pair, with value being read from a file")
                (@arg delete: --delete value_name("key") +multiple
                    "Deletes an internal key-value pair")
            )
        )
        (@subcommand subgrids =>
            (about: "Print information about the internal subgrid types")
            (@arg input: +required "Path to the input grid")
        )
    )
    .get_matches();

    if matches.is_present("silence_lhapdf") {
        lhapdf::set_verbosity(0);
    }

    if let Some(matches) = matches.subcommand_matches("channels") {
        let input = matches.value_of("input").unwrap();
        let pdfset = matches.value_of("pdfset").unwrap();
        let limit = matches.value_of("limit").unwrap().parse()?;
        let orders: Vec<_> = matches
            .values_of("orders")
            .map_or(vec![], |values| values.map(parse_order).collect())
            .into_iter()
            .collect::<Result<_, _>>()?;
        let absolute = matches.is_present("absolute");
        let lumis = parse_integer_list(matches.value_of("lumis").unwrap_or(""))?;

        channels::subcommand(input, pdfset, limit, &orders, absolute, &lumis)?.printstd();
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

        convolute::subcommand(
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
        let ignore_orders = matches.is_present("ignore_orders");

        diff::subcommand(input1, input2, pdfset, ignore_orders)?.printstd();
    } else if let Some(matches) = matches.subcommand_matches("info") {
        let input = matches.value_of("input").unwrap();

        if matches.is_present("ew") || matches.is_present("qcd") {
            info::subcommand_qcd_ew(
                input,
                if matches.is_present("ew") {
                    "ew"
                } else {
                    "qcd"
                },
            )?;
        } else if matches.is_present("get") {
            info::subcommand_get(input, matches.value_of("get").unwrap())?;
        } else if matches.is_present("keys") {
            info::subcommand_keys(input)?;
        } else if matches.is_present("show") {
            info::subcommand_show(input)?;
        } else {
            unreachable!();
        }
    } else if let Some(matches) = matches.subcommand_matches("luminosity") {
        let input = matches.value_of("input").unwrap();

        luminosity::subcommand(input)?.printstd();
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

        return merge::subcommand(
            output,
            input.first().unwrap(),
            &input[1..],
            scale,
            &scale_by_order,
        );
    } else if let Some(matches) = matches.subcommand_matches("optimize") {
        let input = matches.value_of("input").unwrap();
        let output = matches.value_of("output").unwrap();

        optimize::subcommand(input, output)?;
    } else if let Some(matches) = matches.subcommand_matches("orders") {
        let input = matches.value_of("input").unwrap();
        let pdfset = matches.value_of("pdfset").unwrap();
        let absolute = matches.is_present("absolute");
        let normalize: Vec<_> = matches
            .values_of("normalize")
            .map_or(vec![], |values| values.map(parse_order).collect())
            .into_iter()
            .collect::<Result<_, _>>()?;

        orders::subcommand(input, pdfset, absolute, &normalize)?.printstd();
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

        pdf_uncertainty::subcommand(input, pdfset, cl, threads, &orders)?.printstd();
    } else if let Some(matches) = matches.subcommand_matches("remap") {
        let input = matches.value_of("input").unwrap();
        let output = matches.value_of("output").unwrap();
        let remapping = matches.value_of("remapping").unwrap();
        let norm = matches.value_of("norm").unwrap().parse()?;
        let ignore_obs_norm: Vec<_> = matches
            .values_of("ignore_obs_norm")
            .map_or(vec![], |values| values.map(str::parse::<usize>).collect())
            .into_iter()
            .collect::<Result<_, _>>()?;

        remap::subcommand(input, output, remapping, norm, &ignore_obs_norm)?;
    } else if let Some(matches) = matches.subcommand_matches("set") {
        let input = matches.value_of("input").unwrap();
        let output = matches.value_of("output").unwrap();

        let entries: Vec<_> = matches.values_of("entry").map_or(vec![], Iterator::collect);
        let entries_from_file: Vec<_> = matches
            .values_of("entry_from_file")
            .map_or(vec![], Iterator::collect);
        let deletes: Vec<_> = matches
            .values_of("delete")
            .map_or(vec![], Iterator::collect);

        set::subcommand(input, output, &entries, &entries_from_file, &deletes)?;
    } else if let Some(matches) = matches.subcommand_matches("subgrids") {
        let input = matches.value_of("input").unwrap();

        subgrids::subcommand(input)?;
    }

    Ok(())
}
