use super::GlobalConfiguration;
use anyhow::{ensure, Context, Result};
use lhapdf::{Pdf, PdfSet};
use ndarray::Array3;
use pineappl::grid::Grid;
use pineappl::lumi::LumiCache;
use prettytable::format::{FormatBuilder, LinePosition, LineSeparator};
use prettytable::Table;
use std::fs::{File, OpenOptions};
use std::iter;
use std::ops::RangeInclusive;
use std::path::Path;
use std::process::ExitCode;

pub fn create_pdf(pdf: &str) -> Result<Pdf> {
    let pdf = pdf.split_once('=').map_or(pdf, |(name, _)| name);

    Ok(pdf
        .parse()
        .map_or_else(|_| Pdf::with_setname_and_nmem(pdf), Pdf::with_lhaid)?)
}

pub fn create_pdfset(pdfset: &str) -> Result<(PdfSet, Option<usize>)> {
    let pdfset = pdfset.split_once('=').map_or(pdfset, |(name, _)| name);
    let (pdfset, member) = pdfset
        .rsplit_once('/')
        .map_or((pdfset, None), |(set, member)| {
            (set, Some(member.parse::<usize>().unwrap()))
        });

    Ok((
        PdfSet::new(&pdfset.parse().map_or_else(
            |_| pdfset.to_owned(),
            |lhaid| lhapdf::lookup_pdf(lhaid).unwrap().0,
        ))?,
        member,
    ))
}

pub fn pdf_label(pdf: &str) -> &str {
    pdf.split_once('=').map_or(pdf, |(_, label)| label)
}

pub fn read_grid(input: &Path) -> Result<Grid> {
    Grid::read(File::open(input).context(format!("unable to open '{}'", input.display()))?)
        .context(format!("unable to read '{}'", input.display()))
}

pub fn write_grid(output: &Path, grid: &Grid) -> Result<ExitCode> {
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)
        .context(format!("unable to write '{}'", output.display()))?;

    if output.extension().map_or(false, |ext| ext == "lz4") {
        grid.write_lz4(file)?;
    } else {
        grid.write(file)?;
    }

    Ok(ExitCode::SUCCESS)
}

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

pub const SCALES_VECTOR: [(f64, f64); 9] = [
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

pub fn labels_and_units(grid: &Grid, integrated: bool) -> (Vec<(String, &str)>, &str, &str) {
    let key_values = grid.key_values();

    (
        (0..grid.bin_info().dimensions())
            .map(|d| {
                (
                    key_values
                        .and_then(|kv| kv.get(&format!("x{}_label", d + 1)).cloned())
                        .unwrap_or_else(|| format!("x{}", d + 1)),
                    key_values
                        .and_then(|kv| kv.get(&format!("x{}_unit", d + 1)))
                        .map_or("", String::as_str),
                )
            })
            .collect(),
        if integrated {
            "integ"
        } else {
            key_values
                .and_then(|kv| kv.get("y_label").map(String::as_str))
                .unwrap_or("diff")
        },
        if integrated {
            "" // TODO: compute the units for the integrated cross section
        } else {
            key_values
                .and_then(|kv| kv.get("y_unit").map(String::as_str))
                .unwrap_or("")
        },
    )
}

#[derive(Clone, Copy)]
pub enum ConvoluteMode {
    Asymmetry,
    Integrated,
    Normal,
}

pub fn convolve_scales(
    grid: &Grid,
    lhapdf: &mut Pdf,
    orders: &[(u32, u32)],
    bins: &[usize],
    lumis: &[bool],
    scales: &[(f64, f64)],
    mode: ConvoluteMode,
    cfg: &GlobalConfiguration,
) -> Vec<f64> {
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

    // if the field 'Particle' is missing we assume it's a proton PDF
    let pdf_pdg_id = lhapdf
        .set()
        .entry("Particle")
        .map_or(Ok(2212), |string| string.parse::<i32>())
        .unwrap();

    if cfg.force_positive {
        lhapdf.set_force_positive(1);
    }

    let x_max = lhapdf.x_max();
    let x_min = lhapdf.x_min();
    let mut pdf = |id, x, q2| {
        if !cfg.allow_extrapolation && (x < x_min || x > x_max) {
            0.0
        } else {
            lhapdf.xfx_q2(id, x, q2)
        }
    };
    let mut alphas = |q2| lhapdf.alphas_q2(q2);
    let mut cache = LumiCache::with_one(pdf_pdg_id, &mut pdf, &mut alphas);
    let mut results = grid.convolute(&mut cache, &orders, bins, lumis, scales);

    match mode {
        ConvoluteMode::Asymmetry => {
            let bin_count = grid.bin_info().bins();

            // calculating the asymmetry for a subset of bins doesn't work
            assert!((bins.is_empty() || (bins.len() == bin_count)) && (bin_count % 2 == 0));

            results
                .iter()
                .skip((bin_count / 2) * scales.len())
                .zip(
                    results
                        .chunks_exact(scales.len())
                        .take(bin_count / 2)
                        .rev()
                        .flatten(),
                )
                .map(|(pos, neg)| (pos - neg) / (pos + neg))
                .collect()
        }
        ConvoluteMode::Integrated => {
            let normalizations = grid.bin_info().normalizations();

            results
                .iter_mut()
                .zip(
                    normalizations
                        .iter()
                        .enumerate()
                        .filter(|(index, _)| (bins.is_empty() || bins.contains(index)))
                        .flat_map(|(_, norm)| iter::repeat(norm).take(scales.len())),
                )
                .for_each(|(value, norm)| *value *= norm);

            results
        }
        ConvoluteMode::Normal => results,
    }
}

pub fn convolve(
    grid: &Grid,
    lhapdf: &mut Pdf,
    orders: &[(u32, u32)],
    bins: &[usize],
    lumis: &[bool],
    scales: usize,
    mode: ConvoluteMode,
    cfg: &GlobalConfiguration,
) -> Vec<f64> {
    convolve_scales(
        grid,
        lhapdf,
        orders,
        bins,
        lumis,
        &SCALES_VECTOR[0..scales],
        mode,
        cfg,
    )
}

pub fn convolve_limits(grid: &Grid, bins: &[usize], mode: ConvoluteMode) -> Vec<Vec<(f64, f64)>> {
    let limits: Vec<_> = grid
        .bin_info()
        .limits()
        .into_iter()
        .enumerate()
        .filter_map(|(index, limits)| (bins.is_empty() || bins.contains(&index)).then_some(limits))
        .collect();

    match mode {
        ConvoluteMode::Asymmetry => limits[limits.len() / 2..].to_vec(),
        ConvoluteMode::Integrated | ConvoluteMode::Normal => limits,
    }
}

pub fn convolve_subgrid(
    grid: &Grid,
    lhapdf: &mut Pdf,
    order: usize,
    bin: usize,
    lumi: usize,
    cfg: &GlobalConfiguration,
) -> Array3<f64> {
    // if the field 'Particle' is missing we assume it's a proton PDF
    let pdf_pdg_id = lhapdf
        .set()
        .entry("Particle")
        .map_or(Ok(2212), |string| string.parse::<i32>())
        .unwrap();

    if cfg.force_positive {
        lhapdf.set_force_positive(1);
    }

    let x_max = lhapdf.x_max();
    let x_min = lhapdf.x_min();
    let mut pdf = |id, x, q2| {
        if !cfg.allow_extrapolation && (x < x_min || x > x_max) {
            0.0
        } else {
            lhapdf.xfx_q2(id, x, q2)
        }
    };
    let mut alphas = |q2| lhapdf.alphas_q2(q2);
    let mut cache = LumiCache::with_one(pdf_pdg_id, &mut pdf, &mut alphas);

    grid.convolute_subgrid(&mut cache, order, bin, lumi, 1.0, 1.0)
}

pub fn parse_pdfset(argument: &str) -> std::result::Result<String, String> {
    // TODO: figure out how to validate `argument` with `managed-lhapdf`
    Ok(argument.to_owned())
}

pub fn parse_integer_range(range: &str) -> Result<RangeInclusive<usize>> {
    if let Some(at) = range.find('-') {
        let (left, right) = range.split_at(at);
        let left = str::parse::<usize>(left).context(format!(
            "unable to parse integer range '{range}'; couldn't convert '{left}'"
        ))?;
        let right = str::parse::<usize>(&right[1..]).context(format!(
            "unable to parse integer range '{range}'; couldn't convert '{right}'"
        ))?;

        Ok(left..=right)
    } else {
        let value =
            str::parse::<usize>(range).context(format!("unable to parse integer '{range}'"))?;

        Ok(value..=value)
    }
}

pub fn parse_order(order: &str) -> Result<(u32, u32)> {
    let mut alphas = 0;
    let mut alpha = 0;

    let matches: Vec<_> = order.match_indices('a').collect();

    ensure!(
        matches.len() <= 2,
        "unable to parse order; too many couplings in '{}'",
        order
    );

    for (index, _) in matches {
        if &order[index..index + 2] == "as" {
            let len = order[index + 2..]
                .chars()
                .take_while(|c| c.is_numeric())
                .count();
            alphas = str::parse::<u32>(&order[index + 2..index + 2 + len])
                .context(format!("unable to parse order '{order}'"))?;
        } else {
            let len = order[index + 1..]
                .chars()
                .take_while(|c| c.is_numeric())
                .count();
            alpha = str::parse::<u32>(&order[index + 1..index + 1 + len])
                .context(format!("unable to parse order '{order}'"))?;
        }
    }

    Ok((alphas, alpha))
}
