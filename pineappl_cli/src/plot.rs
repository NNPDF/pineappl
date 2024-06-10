use super::helpers::{self, ConvoluteMode};
use super::{GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::builder::{PossibleValuesParser, TypedValueParser};
use clap::{Parser, ValueHint};
use itertools::Itertools;
use ndarray::Axis;
use pineappl::boc::Channel;
use pineappl::convolutions::Convolution;
use pineappl::subgrid::Subgrid;
use rayon::{prelude::*, ThreadPoolBuilder};
use std::fmt::Write;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::slice;
use std::thread;

/// Creates a matplotlib script plotting the contents of the grid.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id(s) or name of the PDF set(s).
    #[arg(required = true, value_parser = helpers::parse_pdfset)]
    pdfsets: Vec<String>,
    /// Set the number of scale variations.
    #[arg(
        default_value_t = 7,
        long,
        short,
        value_parser = PossibleValuesParser::new(["1", "3", "7", "9"]).try_map(|s| s.parse::<usize>())
    )]
    scales: usize,
    /// Show the pull for a specific grid three-dimensionally.
    #[arg(
        conflicts_with = "scales",
        long,
        num_args = 1,
        value_delimiter = ',',
        value_name = "ORDER,BIN,CHAN"
    )]
    subgrid_pull: Vec<String>,
    /// Plot the asymmetry.
    #[arg(conflicts_with = "subgrid_pull", long)]
    asymmetry: bool,
    /// Number of threads to utilize.
    #[arg(default_value_t = thread::available_parallelism().map_or(1, NonZeroUsize::get), long)]
    threads: usize,
    /// Disable the (time-consuming) calculation of PDF uncertainties.
    #[arg(long)]
    no_pdf_unc: bool,
}

fn map_format_join(slice: &[f64]) -> String {
    slice.iter().map(|x| format!("{x}")).join(", ")
}

fn map_format_e_join(slice: &[f64]) -> String {
    slice.iter().map(|x| format!("{x:.7e}")).join(", ")
}

fn map_format_e_join_repeat_last(slice: &[f64]) -> String {
    slice
        .iter()
        .chain(slice.last())
        .map(|x| format!("{x:.7e}"))
        .join(", ")
}

// TODO: this function should take into account what type the particle IDs are
fn map_format_parton(parton: i32) -> &'static str {
    match parton {
        -6 => r"\bar{\mathrm{t}}",
        -5 => r"\bar{\mathrm{b}}",
        -4 => r"\bar{\mathrm{c}}",
        -3 => r"\bar{\mathrm{s}}",
        -2 => r"\bar{\mathrm{u}}",
        -1 => r"\bar{\mathrm{d}}",
        1 => r"\mathrm{d}",
        2 => r"\mathrm{u}",
        3 => r"\mathrm{s}",
        4 => r"\mathrm{c}",
        5 => r"\mathrm{b}",
        6 => r"\mathrm{t}",
        0 | 21 => r"\mathrm{g}",
        22 => r"\gamma",
        100 => r"\Sigma",
        103 => r"\mathrm{T}_3",
        108 => r"\mathrm{T}_8",
        115 => r"\mathrm{T}_{15}",
        124 => r"\mathrm{T}_{24}",
        135 => r"\mathrm{T}_{35}",
        200 => r"\mathrm{V}",
        203 => r"\mathrm{V}_3",
        208 => r"\mathrm{V}_8",
        215 => r"\mathrm{V}_{15}",
        224 => r"\mathrm{V}_{24}",
        235 => r"\mathrm{V}_{35}",
        _ => unimplemented!("PID = {parton} unknown"),
    }
}

fn map_format_channel(channel: &Channel, has_pdf1: bool, has_pdf2: bool) -> String {
    channel
        .entry()
        .iter()
        .map(|&(a, b, _)| {
            format!(
                "{}{}",
                if has_pdf1 { map_format_parton(a) } else { "" },
                if has_pdf2 { map_format_parton(b) } else { "" }
            )
        })
        .join(" + ")
}

fn map_format_channels(channels: &[(String, Vec<f64>)]) -> String {
    channels
        .iter()
        .map(|(label, bins)| {
            format!(
                "                (r\"${}$\", np.array([{}]))",
                label,
                map_format_e_join_repeat_last(bins)
            )
        })
        .join(",\n")
}

fn format_pdf_results(pdf_uncertainties: &[Vec<Vec<f64>>], pdfsets: &[String]) -> String {
    pdf_uncertainties
        .iter()
        .zip(pdfsets.iter())
        .map(|(values, pdfset)| {
            format!(
                "                (
                    \"{}\",
                    np.array([{}]),
                    np.array([{}]),
                    np.array([{}]),
                ),",
                helpers::pdf_label(pdfset).replace('_', r"\_"),
                map_format_e_join_repeat_last(&values[0]),
                map_format_e_join_repeat_last(&values[1]),
                map_format_e_join_repeat_last(&values[2]),
            )
        })
        .join("\n")
}

fn format_metadata(metadata: &[(&String, &String)]) -> String {
    metadata
        .iter()
        .filter_map(|(key, value)| {
            if value.contains('\n') {
                // skip multi-line entries
                None
            } else {
                Some(format!(
                    "        \"{}\": r\"{}\",",
                    key,
                    if *key == "description" {
                        value.replace('\u{2013}', "--").replace('\u{2014}', "---")
                    } else if key.ends_with("_unit") {
                        value
                            .replace("GeV", r"\giga\electronvolt")
                            .replace('/', r"\per")
                            .replace("pb", r"\pico\barn")
                    } else {
                        (*value).clone()
                    }
                ))
            }
        })
        .join("\n")
}

impl Subcommand for Opts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

        if self.subgrid_pull.is_empty() {
            let mode = if self.asymmetry {
                ConvoluteMode::Asymmetry
            } else {
                ConvoluteMode::Normal
            };

            let grid = helpers::read_grid(&self.input)?;
            let mut pdf = helpers::create_pdf(&self.pdfsets[0])?;
            let slices = grid.bin_info().slices();
            let mut data_string = String::new();

            data_string.push_str("[\n");

            for (slice, label) in slices.iter().zip(slices.iter().map(|&(begin, end)| {
                (0..grid.bin_info().dimensions() - 1)
                    .map(|d| {
                        format!(
                            "$\\SI{{{left}}}{{{unit}}} < {obs} < \\SI{{{right}}}{{{unit}}}$",
                            left = grid.bin_info().left(d)[begin],
                            obs = grid
                                .key_values()
                                .and_then(|map| map.get(&format!("x{}_label_tex", d + 1)).cloned())
                                .unwrap_or_else(|| format!("x{}", d + 1))
                                .replace('$', ""),
                            right = grid.bin_info().right(d)[end - 1],
                            unit = grid
                                .key_values()
                                .and_then(|map| map.get(&format!("x{}_unit", d + 1)).cloned())
                                .unwrap_or_default()
                        )
                    })
                    .join(r"\\")
            })) {
                let bins: Vec<_> = (slice.0..slice.1).collect();

                let results = helpers::convolve(
                    &grid,
                    slice::from_mut(&mut pdf),
                    &[],
                    &bins,
                    &[],
                    self.scales,
                    mode,
                    cfg,
                );

                let qcd_results = {
                    let mut orders = grid.orders().to_vec();
                    orders.sort();
                    let orders = orders;

                    let qcd_orders: Vec<_> = orders
                        .iter()
                        .group_by(|order| order.alphas + order.alpha)
                        .into_iter()
                        .map(|mut group| {
                            let order = group.1.next().unwrap();
                            (order.alphas, order.alpha)
                        })
                        .collect();

                    helpers::convolve(
                        &grid,
                        slice::from_mut(&mut pdf),
                        &qcd_orders,
                        &bins,
                        &[],
                        self.scales,
                        mode,
                        cfg,
                    )
                };

                let bin_limits: Vec<_> = helpers::convolve_limits(&grid, &bins, mode)
                    .into_iter()
                    .map(|limits| limits.last().copied().unwrap())
                    .collect();
                let x: Vec<_> = bin_limits
                    .iter()
                    .map(|(left, _)| left)
                    .chain(bin_limits.last().map(|(_, right)| right))
                    .copied()
                    .collect();
                let mid: Vec<_> = x
                    .windows(2)
                    .map(|limits| 0.5 * (limits[0] + limits[1]))
                    .collect();

                let pdf_uncertainties: Vec<Vec<Vec<_>>> = self
                    .pdfsets
                    .par_iter()
                    .map(|pdfset| {
                        if self.no_pdf_unc {
                            let mut pdf = helpers::create_pdf(pdfset).unwrap();

                            let results = helpers::convolve(
                                &grid,
                                slice::from_mut(&mut pdf),
                                &[],
                                &bins,
                                &[],
                                1,
                                mode,
                                cfg,
                            );

                            Ok(vec![results; 3])
                        } else {
                            let (set, member) = helpers::create_pdfset(pdfset).unwrap();

                            let pdf_results: Vec<_> = set
                                .mk_pdfs()?
                                .into_par_iter()
                                .flat_map(|mut pdf| {
                                    helpers::convolve(
                                        &grid,
                                        slice::from_mut(&mut pdf),
                                        &[],
                                        &bins,
                                        &[],
                                        1,
                                        mode,
                                        cfg,
                                    )
                                })
                                .collect();

                            let bins = mid.len();

                            let mut central = Vec::with_capacity(bins);
                            let mut min = Vec::with_capacity(bins);
                            let mut max = Vec::with_capacity(bins);

                            for bin in 0..bins {
                                let values: Vec<_> = pdf_results
                                    .iter()
                                    .skip(bin)
                                    .step_by(bins)
                                    .copied()
                                    .collect();

                                let uncertainty =
                                    set.uncertainty(&values, lhapdf::CL_1_SIGMA, false).unwrap();
                                central.push(
                                    member.map_or(uncertainty.central, |member| values[member]),
                                );
                                min.push(uncertainty.central - uncertainty.errminus);
                                max.push(uncertainty.central + uncertainty.errplus);
                            }

                            Ok(vec![central, min, max])
                        }
                    })
                    .collect::<Result<_>>()?;

                let central: Vec<_> = results.iter().step_by(self.scales).copied().collect();
                let min: Vec<_> = results
                    .chunks_exact(self.scales)
                    .map(|variations| {
                        variations
                            .iter()
                            .min_by(|a, b| a.total_cmp(b))
                            .copied()
                            .unwrap()
                    })
                    .collect();
                let max: Vec<_> = results
                    .chunks_exact(self.scales)
                    .map(|variations| {
                        variations
                            .iter()
                            .max_by(|a, b| a.total_cmp(b))
                            .copied()
                            .unwrap()
                    })
                    .collect();

                let qcd_central: Vec<_> =
                    qcd_results.iter().step_by(self.scales).copied().collect();
                let qcd_min: Vec<_> = qcd_results
                    .chunks_exact(self.scales)
                    .map(|variations| {
                        variations
                            .iter()
                            .min_by(|a, b| a.total_cmp(b))
                            .copied()
                            .unwrap()
                    })
                    .collect();
                let qcd_max: Vec<_> = qcd_results
                    .chunks_exact(self.scales)
                    .map(|variations| {
                        variations
                            .iter()
                            .max_by(|a, b| a.total_cmp(b))
                            .copied()
                            .unwrap()
                    })
                    .collect();

                let channels = if matches!(mode, ConvoluteMode::Asymmetry) {
                    vec![]
                } else {
                    let mut channels: Vec<_> = (0..grid.channels().len())
                        .map(|channel| {
                            let mut channel_mask = vec![false; grid.channels().len()];
                            channel_mask[channel] = true;
                            (
                                map_format_channel(
                                    &grid.channels()[channel],
                                    grid.convolutions()[0] != Convolution::None,
                                    grid.convolutions()[1] != Convolution::None,
                                ),
                                helpers::convolve(
                                    &grid,
                                    slice::from_mut(&mut pdf),
                                    &[],
                                    &bins,
                                    &channel_mask,
                                    1,
                                    mode,
                                    cfg,
                                ),
                            )
                        })
                        .collect();

                    // sort channels by importance
                    channels.sort_by(|(_, lhs), (_, rhs)| {
                        let lhs = lhs.iter().fold(0.0, |prev, x| prev + x.abs());
                        let rhs = rhs.iter().fold(0.0, |prev, x| prev + x.abs());
                        rhs.total_cmp(&lhs)
                    });
                    channels
                };

                writeln!(
                    &mut data_string,
                    "        {{
            \"slice_label\"    : r\"{slice_label}\",
            \"x\"        : np.array([{x}]),
            \"mid\"      : np.array([{mid}]),
            \"pdf_results\" : [
{pdf_results}
            ],
            \"qcd_y\"    : np.array([{qcd_y}]),
            \"qcd_min\"  : np.array([{qcd_min}]),
            \"qcd_max\"  : np.array([{qcd_max}]),
            \"y\"        : np.array([{y}]),
            \"ymin\"     : np.array([{ymin}]),
            \"ymax\"     : np.array([{ymax}]),
            \"channels\" : [
{channels}
            ],
        }},",
                    slice_label = label,
                    mid = map_format_join(&mid),
                    pdf_results = format_pdf_results(&pdf_uncertainties, &self.pdfsets),
                    qcd_y = map_format_e_join_repeat_last(&qcd_central),
                    qcd_min = map_format_e_join_repeat_last(&qcd_min),
                    qcd_max = map_format_e_join_repeat_last(&qcd_max),
                    x = map_format_join(&x),
                    y = map_format_e_join_repeat_last(&central),
                    ymin = map_format_e_join_repeat_last(&min),
                    ymax = map_format_e_join_repeat_last(&max),
                    channels = map_format_channels(&channels),
                )
                .unwrap_or_else(|_| unreachable!());
            }

            data_string.push_str("    ]");

            // prepare metadata
            let key_values = grid.key_values().cloned().unwrap_or_default();
            let mut vector: Vec<_> = key_values.iter().collect();
            vector.sort();
            let vector = vector;

            let mut output = self.input.clone();

            // remove ".lz4" and ".pineappl" extension
            if let Some(x) = output.extension() {
                if x == "lz4" {
                    output = Path::new(output.file_stem().unwrap()).to_path_buf();
                }
            }
            if let Some(x) = output.extension() {
                if x == "pineappl" {
                    output = Path::new(output.file_stem().unwrap()).to_path_buf();
                }
            }

            let xaxis = format!("x{}", grid.bin_info().dimensions());
            let xunit = key_values
                .get(&format!("{xaxis}_unit"))
                .map_or("", String::as_str);
            let xlabel = format!(
                "{}{}",
                key_values
                    .get(&format!("{xaxis}_label_tex"))
                    .map_or("", String::as_str),
                if xunit.is_empty() {
                    String::new()
                } else {
                    format!(" [\\si{{{xunit}}}]")
                }
            );
            let yunit = key_values.get("y_unit").map_or("", String::as_str);
            let ylabel = format!(
                "{}{}",
                key_values.get("y_label_tex").map_or("", String::as_str),
                if yunit.is_empty() {
                    String::new()
                } else {
                    format!(" [\\si{{{yunit}}}]")
                }
            );
            let xlog = if xunit.is_empty() { "False" } else { "True" };
            let ylog = xlog;
            let title = key_values.get("description").map_or("", String::as_str);
            let bins = grid.bin_info().bins();
            let pdfs = self.pdfsets.len();

            print!(
                include_str!("plot.py"),
                inte = if bins == 1 { "" } else { "# " },
                nint = if bins == 1 { "# " } else { "" },
                pdfs = if pdfs == 1 || bins == 1 { "# " } else { "" },
                xlabel = xlabel,
                ylabel = ylabel,
                xlog = xlog,
                ylog = ylog,
                title = title,
                output = output.to_str().unwrap(),
                data = data_string,
                metadata = format_metadata(&vector),
            );
        } else {
            let (pdfset1, pdfset2) = self.pdfsets.iter().collect_tuple().unwrap();
            let (order, bin, channel) = self
                .subgrid_pull
                .iter()
                .map(|num| num.parse::<usize>().unwrap())
                .collect_tuple()
                .unwrap();

            let cl = lhapdf::CL_1_SIGMA;
            let grid = helpers::read_grid(&self.input)?;

            let (set1, member1) = helpers::create_pdfset(pdfset1)?;
            let (set2, member2) = helpers::create_pdfset(pdfset2)?;
            let mut pdfset1 = set1.mk_pdfs()?;
            let mut pdfset2 = set2.mk_pdfs()?;

            let values1: Vec<_> = pdfset1
                .par_iter_mut()
                .map(|pdf| {
                    let values = helpers::convolve(
                        &grid,
                        slice::from_mut(pdf),
                        &[],
                        &[bin],
                        &[],
                        1,
                        ConvoluteMode::Normal,
                        cfg,
                    );
                    assert_eq!(values.len(), 1);
                    values[0]
                })
                .collect();
            let values2: Vec<_> = pdfset2
                .par_iter_mut()
                .map(|pdf| {
                    let values = helpers::convolve(
                        &grid,
                        slice::from_mut(pdf),
                        &[],
                        &[bin],
                        &[],
                        1,
                        ConvoluteMode::Normal,
                        cfg,
                    );
                    assert_eq!(values.len(), 1);
                    values[0]
                })
                .collect();

            let uncertainty1 = set1.uncertainty(&values1, cl, false)?;
            let uncertainty2 = set2.uncertainty(&values2, cl, false)?;
            let central1 = member1.map_or(uncertainty1.central, |member1| values1[member1]);
            let central2 = member2.map_or(uncertainty2.central, |member2| values2[member2]);

            let denominator = {
                // use the uncertainties in the direction in which the respective results differ
                let unc1 = if central1 > central2 {
                    uncertainty1.errminus
                } else {
                    uncertainty1.errplus
                };
                let unc2 = if central2 > central1 {
                    uncertainty2.errminus
                } else {
                    uncertainty2.errplus
                };

                unc1.hypot(unc2)
            };

            let res1 = helpers::convolve_subgrid(&grid, &mut pdfset1[0], order, bin, channel, cfg)
                .sum_axis(Axis(0));
            let res2 = helpers::convolve_subgrid(&grid, &mut pdfset2[0], order, bin, channel, cfg)
                .sum_axis(Axis(0));

            let subgrid = &grid.subgrids()[[order, bin, channel]];
            //let q2 = subgrid.q2_grid();
            let x1 = subgrid.x1_grid();
            let x2 = subgrid.x2_grid();

            let mut x1_vals = vec![];
            let mut x2_vals = vec![];
            let mut vals = vec![];

            for (((ix1, ix2), &one), &two) in res1.indexed_iter().zip(res2.iter()) {
                if one == 0.0 {
                    assert_eq!(two, 0.0);
                    continue;
                }

                x1_vals.push(x1[ix1]);
                x2_vals.push(x2[ix2]);
                vals.push((two - one) / denominator);
            }

            print!(
                include_str!("subgrid-pull-plot.py"),
                x1 = map_format_e_join(&x1_vals),
                x2 = map_format_e_join(&x2_vals),
                z = map_format_e_join(&vals)
            );
        }

        Ok(ExitCode::SUCCESS)
    }
}
