use super::GlobalConfiguration;
use anyhow::{anyhow, ensure, Context, Error, Result};
use lhapdf::{Pdf, PdfSet};
use ndarray::Array3;
use pineappl::convolutions::LumiCache;
use pineappl::grid::Grid;
use prettytable::format::{FormatBuilder, LinePosition, LineSeparator};
use prettytable::Table;
use std::fs::{File, OpenOptions};
use std::iter;
use std::ops::RangeInclusive;
use std::path::Path;
use std::process::ExitCode;
use std::str::FromStr;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConvFuns {
    pub lhapdf_names: Vec<String>,
    pub members: Vec<Option<usize>>,
    pub label: String,
}

impl FromStr for ConvFuns {
    type Err = Error;

    fn from_str(arg: &str) -> std::result::Result<Self, Self::Err> {
        let (names, label) = arg.split_once('=').unwrap_or((arg, arg));
        let (lhapdf_names, members) = names
            .split(',')
            .map(|fun| {
                Ok::<_, Error>(if let Some((name, mem)) = fun.split_once('/') {
                    (name.to_owned(), Some(mem.parse()?))
                } else {
                    (fun.to_owned(), None)
                })
            })
            .collect::<Result<Vec<(_, _)>, _>>()?
            .into_iter()
            .unzip();

        Ok(Self {
            lhapdf_names,
            members,
            label: label.to_owned(),
        })
    }
}

pub fn create_conv_funs(funs: &ConvFuns) -> Result<Vec<Pdf>> {
    Ok(funs
        .lhapdf_names
        .iter()
        .zip(&funs.members)
        .map(|(lhapdf_name, member)| {
            lhapdf_name.parse().map_or_else(
                |_| {
                    let member = member.unwrap_or(0);
                    // UNWRAP: we don't support sets with more members than `i32`
                    Pdf::with_setname_and_member(lhapdf_name, member.try_into().unwrap())
                },
                Pdf::with_lhaid,
            )
        })
        .collect::<Result<_, _>>()?)
}

pub fn create_conv_funs_for_set(
    funs: &ConvFuns,
    index_of_set: usize,
) -> Result<(PdfSet, Vec<Vec<Pdf>>)> {
    let setname = &funs.lhapdf_names[index_of_set];
    let set = setname.parse().map_or_else(
        |_| Ok::<_, Error>(PdfSet::new(setname)?),
        |lhaid| {
            Ok(PdfSet::new(
                &lhapdf::lookup_pdf(lhaid)
                    .map(|(set, _)| set)
                    .ok_or_else(|| {
                        anyhow!("no convolution function for LHAID = `{lhaid}` found")
                    })?,
            )?)
        },
    )?;

    let conv_funs = set
        .mk_pdfs()?
        .into_iter()
        .map(|conv_fun| {
            // TODO: do not create objects that are getting overwritten in any case
            let mut conv_funs = create_conv_funs(funs)?;
            conv_funs[index_of_set] = conv_fun;

            Ok::<_, Error>(conv_funs)
        })
        .collect::<Result<_, _>>()?;

    Ok((set, conv_funs))
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
    let metadata = grid.metadata();

    (
        (0..grid.bin_info().dimensions())
            .map(|d| {
                (
                    metadata
                        .get(&format!("x{}_label", d + 1))
                        .cloned()
                        .unwrap_or_else(|| format!("x{}", d + 1)),
                    metadata
                        .get(&format!("x{}_unit", d + 1))
                        .map_or("", String::as_str),
                )
            })
            .collect(),
        if integrated {
            "integ"
        } else {
            metadata.get("y_label").map_or("diff", String::as_str)
        },
        if integrated {
            "" // TODO: compute the units for the integrated cross section
        } else {
            metadata.get("y_unit").map_or("", String::as_str)
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
    conv_funs: &mut [Pdf],
    orders: &[(u32, u32)],
    bins: &[usize],
    channels: &[bool],
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

    if cfg.force_positive {
        for fun in conv_funs.iter_mut() {
            fun.set_force_positive(1);
        }
    }

    let mut results = match conv_funs {
        [fun] => {
            // there's only one convolution function from which we can use the strong coupling
            assert_eq!(cfg.use_alphas_from, 0);

            // if the field 'Particle' is missing we assume it's a proton PDF
            let pdg_id = fun
                .set()
                .entry("Particle")
                .map_or(Ok(2212), |string| string.parse::<i32>())
                .unwrap();

            let x_max = fun.x_max();
            let x_min = fun.x_min();
            let mut alphas = |q2| fun.alphas_q2(q2);
            let mut fun = |id, x, q2| {
                if !cfg.allow_extrapolation && (x < x_min || x > x_max) {
                    0.0
                } else {
                    fun.xfx_q2(id, x, q2)
                }
            };

            let mut cache = LumiCache::with_one(pdg_id, &mut fun, &mut alphas);

            grid.convolve(&mut cache, &orders, bins, channels, scales)
        }
        [fun1, fun2] => {
            let pdg_id1 = fun1
                .set()
                .entry("Particle")
                .map_or(Ok(2212), |string| string.parse::<i32>())
                .unwrap();

            let pdg_id2 = fun2
                .set()
                .entry("Particle")
                .map_or(Ok(2212), |string| string.parse::<i32>())
                .unwrap();

            let x_max1 = fun1.x_max();
            let x_min1 = fun1.x_min();
            let x_max2 = fun2.x_max();
            let x_min2 = fun2.x_min();

            let mut alphas = |q2| match cfg.use_alphas_from {
                0 => fun1.alphas_q2(q2),
                1 => fun2.alphas_q2(q2),
                _ => panic!(
                    "expected `use_alphas_from` to be `0` or `1`, is {}",
                    cfg.use_alphas_from
                ),
            };
            let mut fun1 = |id, x, q2| {
                if !cfg.allow_extrapolation && (x < x_min1 || x > x_max1) {
                    0.0
                } else {
                    fun1.xfx_q2(id, x, q2)
                }
            };

            let mut fun2 = |id, x, q2| {
                if !cfg.allow_extrapolation && (x < x_min2 || x > x_max2) {
                    0.0
                } else {
                    fun2.xfx_q2(id, x, q2)
                }
            };

            let mut cache =
                LumiCache::with_two(pdg_id1, &mut fun1, pdg_id2, &mut fun2, &mut alphas);

            grid.convolve(&mut cache, &orders, bins, channels, scales)
        }
        _ => unimplemented!(),
    };

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
    conv_funs: &mut [Pdf],
    orders: &[(u32, u32)],
    bins: &[usize],
    lumis: &[bool],
    scales: usize,
    mode: ConvoluteMode,
    cfg: &GlobalConfiguration,
) -> Vec<f64> {
    convolve_scales(
        grid,
        conv_funs,
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
    conv_funs: &mut [Pdf],
    order: usize,
    bin: usize,
    lumi: usize,
    cfg: &GlobalConfiguration,
) -> Array3<f64> {
    if cfg.force_positive {
        for fun in conv_funs.iter_mut() {
            fun.set_force_positive(1);
        }
    }

    match conv_funs {
        [fun] => {
            // there's only one convolution function from which we can use the strong coupling
            assert_eq!(cfg.use_alphas_from, 0);

            // if the field 'Particle' is missing we assume it's a proton PDF
            let pdg_id = fun
                .set()
                .entry("Particle")
                .map_or(Ok(2212), |string| string.parse::<i32>())
                .unwrap();

            let x_max = fun.x_max();
            let x_min = fun.x_min();
            let mut alphas = |q2| fun.alphas_q2(q2);
            let mut fun = |id, x, q2| {
                if !cfg.allow_extrapolation && (x < x_min || x > x_max) {
                    0.0
                } else {
                    fun.xfx_q2(id, x, q2)
                }
            };

            let mut cache = LumiCache::with_one(pdg_id, &mut fun, &mut alphas);

            grid.convolve_subgrid(&mut cache, order, bin, lumi, 1.0, 1.0)
        }
        [fun1, fun2] => {
            let pdg_id1 = fun1
                .set()
                .entry("Particle")
                .map_or(Ok(2212), |string| string.parse::<i32>())
                .unwrap();

            let pdg_id2 = fun2
                .set()
                .entry("Particle")
                .map_or(Ok(2212), |string| string.parse::<i32>())
                .unwrap();

            let x_max1 = fun1.x_max();
            let x_min1 = fun1.x_min();
            let x_max2 = fun2.x_max();
            let x_min2 = fun2.x_min();

            let mut alphas = |q2| match cfg.use_alphas_from {
                0 => fun1.alphas_q2(q2),
                1 => fun2.alphas_q2(q2),
                _ => panic!(
                    "expected `use_alphas_from` to be `0` or `1`, is {}",
                    cfg.use_alphas_from
                ),
            };
            let mut fun1 = |id, x, q2| {
                if !cfg.allow_extrapolation && (x < x_min1 || x > x_max1) {
                    0.0
                } else {
                    fun1.xfx_q2(id, x, q2)
                }
            };

            let mut fun2 = |id, x, q2| {
                if !cfg.allow_extrapolation && (x < x_min2 || x > x_max2) {
                    0.0
                } else {
                    fun2.xfx_q2(id, x, q2)
                }
            };

            let mut cache =
                LumiCache::with_two(pdg_id1, &mut fun1, pdg_id2, &mut fun2, &mut alphas);

            grid.convolve_subgrid(&mut cache, order, bin, lumi, 1.0, 1.0)
        }
        _ => unimplemented!(),
    }
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

#[cfg(test)]
mod test {
    use super::ConvFuns;

    #[test]
    fn conv_fun_from_str() {
        assert_eq!(
            "A/2,B/1,C/0,D=X".parse::<ConvFuns>().unwrap(),
            ConvFuns {
                lhapdf_names: vec![
                    "A".to_owned(),
                    "B".to_owned(),
                    "C".to_owned(),
                    "D".to_owned()
                ],
                members: vec![Some(2), Some(1), Some(0), None],
                label: "X".to_owned()
            }
        );
    }
}
