use super::GlobalConfiguration;
use anyhow::{anyhow, bail, ensure, Context, Error, Result};
use itertools::Itertools;
use lhapdf::{Pdf, PdfSet};
use pineappl::boc::{ScaleFuncForm, Scales};
use pineappl::convolutions::{Conv, ConvType, ConvolutionCache};
use pineappl::grid::Grid;
use prettytable::format::{FormatBuilder, LinePosition, LineSeparator};
use prettytable::Table;
use std::fs::{File, OpenOptions};
use std::iter;
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::str::FromStr;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConvFuns {
    pub lhapdf_names: Vec<String>,
    pub members: Vec<Option<usize>>,
    pub conv_types: Vec<ConvType>,
    pub label: String,
}

impl FromStr for ConvFuns {
    type Err = Error;

    fn from_str(arg: &str) -> std::result::Result<Self, Self::Err> {
        let (names, label) = arg.split_once('=').unwrap_or((arg, arg));
        let (lhapdf_names, members, conv_types) = names
            .split(',')
            .map(|fun| {
                let (name, typ) = fun.split_once('+').unwrap_or((fun, ""));
                let (name, mem) = name.split_once('/').map_or((name, None), |(name, mem)| {
                    (
                        name,
                        Some(
                            mem.parse()
                                // TODO: do proper error handling
                                .unwrap(),
                        ),
                    )
                });
                let name = name.to_owned();
                let typ = match typ {
                    "" => ConvType::UnpolPDF,
                    "p" => ConvType::PolPDF,
                    "f" => ConvType::UnpolFF,
                    "pf" | "fp" => ConvType::PolFF,
                    _ => bail!("unknown convolution type '{typ}'"),
                };
                Ok::<_, Error>((name, mem, typ))
            })
            .collect::<Result<Vec<(_, _, _)>, _>>()?
            .into_iter()
            .multiunzip();

        Ok(Self {
            lhapdf_names,
            members,
            conv_types,
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

    if output.extension().is_some_and(|ext| ext == "lz4") {
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

pub const SCALES_VECTOR_REN_FAC: [(f64, f64, f64); 9] = [
    (1.0, 1.0, 1.0),
    (2.0, 2.0, 1.0),
    (0.5, 0.5, 1.0),
    (2.0, 1.0, 1.0),
    (1.0, 2.0, 1.0),
    (0.5, 1.0, 1.0),
    (1.0, 0.5, 1.0),
    (2.0, 0.5, 1.0),
    (0.5, 2.0, 1.0),
];

const SCALES_VECTOR_REN_FRG: [(f64, f64, f64); 9] = [
    (1.0, 1.0, 1.0),
    (2.0, 1.0, 2.0),
    (0.5, 1.0, 0.5),
    (2.0, 1.0, 1.0),
    (1.0, 1.0, 2.0),
    (0.5, 1.0, 1.0),
    (1.0, 1.0, 0.5),
    (2.0, 1.0, 0.5),
    (0.5, 1.0, 2.0),
];

const SCALES_VECTOR_27: [(f64, f64, f64); 27] = [
    (1.0, 1.0, 1.0),
    (2.0, 2.0, 2.0),
    (0.5, 0.5, 0.5),
    (0.5, 0.5, 1.0),
    (0.5, 1.0, 0.5),
    (0.5, 1.0, 1.0),
    (0.5, 1.0, 2.0),
    (1.0, 0.5, 0.5),
    (1.0, 0.5, 1.0),
    (1.0, 1.0, 0.5),
    (1.0, 1.0, 2.0),
    (1.0, 2.0, 1.0),
    (1.0, 2.0, 2.0),
    (2.0, 1.0, 0.5),
    (2.0, 1.0, 1.0),
    (2.0, 1.0, 2.0),
    (2.0, 2.0, 1.0),
    (2.0, 0.5, 0.5),
    (0.5, 2.0, 0.5),
    (1.0, 2.0, 0.5),
    (2.0, 2.0, 0.5),
    (2.0, 0.5, 1.0),
    (0.5, 2.0, 1.0),
    (0.5, 0.5, 2.0),
    (1.0, 0.5, 2.0),
    (2.0, 0.5, 2.0),
    (0.5, 2.0, 2.0),
];

pub fn labels_and_units(grid: &Grid, integrated: bool) -> (Vec<(String, &str)>, &str, &str) {
    let metadata = grid.metadata();

    (
        (0..grid.bwfl().dimensions())
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
    conv_types: &[ConvType],
    orders: &[(u8, u8)],
    bins: &[usize],
    channels: &[bool],
    scales: &[(f64, f64, f64)],
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

    // TODO: promote this to an error
    assert!(
        cfg.use_alphas_from < conv_funs.len(),
        "expected `use_alphas_from` to be an integer within `[0, {})`, but got `{}`",
        conv_funs.len(),
        cfg.use_alphas_from
    );

    let x_min_max: Vec<_> = conv_funs
        .iter_mut()
        .map(|fun| (fun.x_min(), fun.x_max()))
        .collect();
    let mut funs: Vec<_> = conv_funs
        .iter()
        .zip(x_min_max)
        .map(|(fun, (x_min, x_max))| {
            move |id, x, q2| {
                if !cfg.allow_extrapolation && (x < x_min || x > x_max) {
                    0.0
                } else {
                    fun.xfx_q2(id, x, q2)
                }
            }
        })
        .collect();
    let xfx: Vec<_> = funs
        .iter_mut()
        .map(|fun| fun as &mut dyn FnMut(i32, f64, f64) -> f64)
        .collect();
    let mut alphas_funs: Vec<_> = conv_funs
        .iter()
        .map(|fun| move |q2| fun.alphas_q2(q2))
        .collect();
    let convolutions: Vec<_> = conv_funs
        .iter()
        .zip(conv_types)
        .map(|(fun, &conv_type)| {
            let pid = fun
                .set()
                .entry("Particle")
                // if the field 'Particle' is missing we assume it's a proton PDF
                .map_or(Ok(2212), |string| string.parse::<i32>())
                // UNWRAP: if this fails, there's a non-integer string in the LHAPDF info file
                .unwrap();

            Conv::new(conv_type, pid)
        })
        .collect();

    let mut cache = ConvolutionCache::new(convolutions, xfx, &mut alphas_funs[cfg.use_alphas_from]);
    let mut results = grid.convolve(&mut cache, &orders, bins, channels, scales);

    match mode {
        ConvoluteMode::Asymmetry => {
            let bin_count = grid.bwfl().len();

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
            results
                .iter_mut()
                .zip(
                    grid.bwfl()
                        .normalizations()
                        .into_iter()
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

pub fn scales_vector(grid: &Grid, scales: usize) -> &[(f64, f64, f64)] {
    let Scales { fac, frg, .. } = grid.scales();

    match (fac, frg, scales) {
        (_, _, 1) => &SCALES_VECTOR_27[0..1],
        (_, _, 3) => &SCALES_VECTOR_27[0..3],
        (_, ScaleFuncForm::NoScale, 7) => &SCALES_VECTOR_REN_FAC[0..7],
        (_, ScaleFuncForm::NoScale, 9) => &SCALES_VECTOR_REN_FAC[..],
        (ScaleFuncForm::NoScale, _, 7) => &SCALES_VECTOR_REN_FRG[0..7],
        (ScaleFuncForm::NoScale, _, 9) => &SCALES_VECTOR_REN_FRG[..],
        (_, _, 17) => &SCALES_VECTOR_27[0..17],
        (_, _, 27) => &SCALES_VECTOR_27[..],
        _ => unreachable!(),
    }
}

pub fn convolve(
    grid: &Grid,
    conv_funs: &mut [Pdf],
    conv_types: &[ConvType],
    orders: &[(u8, u8)],
    bins: &[usize],
    lumis: &[bool],
    scales: usize,
    mode: ConvoluteMode,
    cfg: &GlobalConfiguration,
) -> Vec<f64> {
    convolve_scales(
        grid,
        conv_funs,
        conv_types,
        orders,
        bins,
        lumis,
        scales_vector(grid, scales),
        mode,
        cfg,
    )
}

pub fn convolve_limits(grid: &Grid, bins: &[usize], mode: ConvoluteMode) -> Vec<Vec<(f64, f64)>> {
    let limits: Vec<_> = grid
        .bwfl()
        .bins()
        .iter()
        .map(|bin| bin.limits().to_vec())
        .enumerate()
        .filter_map(|(index, limits)| (bins.is_empty() || bins.contains(&index)).then_some(limits))
        .collect();

    match mode {
        ConvoluteMode::Asymmetry => limits[limits.len() / 2..].to_vec(),
        ConvoluteMode::Integrated | ConvoluteMode::Normal => limits,
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

pub fn parse_order(order: &str) -> Result<(u8, u8)> {
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
            alphas = str::parse::<u8>(&order[index + 2..index + 2 + len])
                .context(format!("unable to parse order '{order}'"))?;
        } else {
            let len = order[index + 1..]
                .chars()
                .take_while(|c| c.is_numeric())
                .count();
            alpha = str::parse::<u8>(&order[index + 1..index + 1 + len])
                .context(format!("unable to parse order '{order}'"))?;
        }
    }

    Ok((alphas, alpha))
}

pub fn parse_ekos(ekos_str: &str) -> Vec<PathBuf> {
    let eko_names = ekos_str.split(',');
    eko_names.into_iter().map(PathBuf::from).collect()
}

#[cfg(test)]
mod test {
    use super::ConvFuns;
    use pineappl::convolutions::ConvType;

    #[test]
    fn conv_fun_from_str() {
        assert_eq!(
            "A/2+p,B/1+f,C/0+fp,D+pf,E=X".parse::<ConvFuns>().unwrap(),
            ConvFuns {
                lhapdf_names: vec![
                    "A".to_owned(),
                    "B".to_owned(),
                    "C".to_owned(),
                    "D".to_owned(),
                    "E".to_owned()
                ],
                members: vec![Some(2), Some(1), Some(0), None, None],
                conv_types: vec![
                    ConvType::PolPDF,
                    ConvType::UnpolFF,
                    ConvType::PolFF,
                    ConvType::PolFF,
                    ConvType::UnpolPDF
                ],
                label: "X".to_owned()
            }
        );
    }
}
