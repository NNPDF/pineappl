use super::helpers::{self, ConvFuns, ConvoluteMode};
use super::{GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::builder::{PossibleValuesParser, TypedValueParser};
use clap::{Parser, ValueHint};
use float_cmp::assert_approx_eq;
use itertools::Itertools;
use ndarray::Axis;
use pineappl::boc::{Channel, Kinematics};
use pineappl::grid::Grid;
use pineappl::subgrid::Subgrid;
use rayon::{prelude::*, ThreadPoolBuilder};
use std::fmt::Write;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::thread;

const MARKER_CONFIG_BEGIN: &str = "# CLI_INSERT_CONFIG_BEGIN\n";
const MARKER_CONFIG_END: &str = "# CLI_INSERT_CONFIG_END";
const MARKER_DATA_INSERT: &str = "# CLI_INSERT_DATA";

/// Creates a matplotlib script plotting the contents of the grid.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id(s) or name of the PDF set(s).
    #[arg(required = true)]
    conv_funs: Vec<ConvFuns>,
    /// Choose for which convolution function the uncertainty should be calculated.
    #[arg(default_value = "0", long, value_name = "IDX")]
    conv_fun_uncert_from: usize,
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
    no_conv_fun_unc: bool,
}

/// Convert `slice` to (unformatted) Python list.
fn map_format_join(slice: &[f64]) -> String {
    slice.iter().map(|x| format!("{x}")).join(", ")
}

/// Convert `slice` to Python list using `.7e`.
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

/// Convert a channel to a good Python string representation.
fn map_format_channel(
    channel: &Channel,
    has_pdf1: bool,
    has_pdf2: bool,
    pid_basis: PidBasis,
) -> String {
    channel
        .entry()
        .iter()
        .map(|(pids, _)| {
            pids.iter()
                .map(|&pid| grid.pid_basis().to_latex_str(pid))
                .collect::<Vec<_>>()
                .join("")
        })
        .join(" + ")
}

/// Convert channel contributions to Python tuples.
fn map_format_channels(channels: &[(String, Vec<f64>)]) -> String {
    channels
        .iter()
        .map(|(label, bins)| {
            format!(
                "            (r\"${}$\", np.array([{}]))",
                label,
                map_format_e_join_repeat_last(bins)
            )
        })
        .join(",\n")
}

/// Convert PDF results to a Python tuple.
fn format_pdf_results(pdf_uncertainties: &[Vec<Vec<f64>>], conv_funs: &[ConvFuns]) -> String {
    pdf_uncertainties
        .iter()
        .zip(conv_funs.iter().map(|fun| &fun.label))
        .map(|(values, label)| {
            format!(
                "            (
                    r\"{}\",
                    np.array([{}]),
                    np.array([{}]),
                    np.array([{}]),
                ),",
                label.replace('_', r"\_"),
                map_format_e_join_repeat_last(&values[0]),
                map_format_e_join_repeat_last(&values[1]),
                map_format_e_join_repeat_last(&values[2]),
            )
        })
        .join("\n")
}

/// Convert metadata into a Python dict.
fn format_metadata(metadata: &[(&String, &String)]) -> String {
    metadata
        .iter()
        .filter_map(|(key, value)| {
            if value.contains('\n') {
                // skip multi-line entries
                None
            } else {
                Some(format!(
                    "    \"{}\": r\"{}\",",
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

/// Convert `b` into a Python bool literal.
fn map_bool(b: bool) -> String {
    if b {
        "True".to_owned()
    } else {
        "False".to_owned()
    }
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
            let mut conv_funs = helpers::create_conv_funs(&self.conv_funs[0])?;
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
                                .metadata()
                                .get(&format!("x{}_label_tex", d + 1))
                                .cloned()
                                .unwrap_or_else(|| format!("x{}", d + 1))
                                .replace('$', ""),
                            right = grid.bin_info().right(d)[end - 1],
                            unit = grid
                                .metadata()
                                .get(&format!("x{}_unit", d + 1))
                                .cloned()
                                .unwrap_or_default()
                        )
                    })
                    .join(r"\\")
            })) {
                let bins: Vec<_> = (slice.0..slice.1).collect();

                let results = helpers::convolve(
                    &grid,
                    &mut conv_funs,
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
                        &mut conv_funs,
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

                let conv_fun_uncertainties: Vec<Vec<Vec<_>>> = self
                    .conv_funs
                    .par_iter()
                    .map(|conv_funs| {
                        if self.no_conv_fun_unc {
                            let mut conv_funs = helpers::create_conv_funs(conv_funs)?;

                            let results = helpers::convolve(
                                &grid,
                                &mut conv_funs,
                                &[],
                                &bins,
                                &[],
                                1,
                                mode,
                                cfg,
                            );

                            Ok(vec![results; 3])
                        } else {
                            let (set, funs) = helpers::create_conv_funs_for_set(
                                conv_funs,
                                self.conv_fun_uncert_from,
                            )?;

                            let pdf_results: Vec<_> = funs
                                .into_par_iter()
                                .flat_map(|mut funs| {
                                    helpers::convolve(
                                        &grid,
                                        &mut funs,
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
                                    conv_funs.members[self.conv_fun_uncert_from]
                                        .map_or(uncertainty.central, |member| values[member]),
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
                                map_format_channel(&grid.channels()[channel], &grid),
                                helpers::convolve(
                                    &grid,
                                    &mut conv_funs,
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
                    "    {{
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
                    pdf_results = format_pdf_results(&conv_fun_uncertainties, &self.conv_funs),
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

            data_string.push_str("]");

            // prepare metadata
            let metadata = grid.metadata();
            let vector: Vec<_> = metadata.iter().collect();

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
            let xunit = metadata
                .get(&format!("{xaxis}_unit"))
                .map_or("", String::as_str);
            let xlabel = format!(
                "{}{}",
                metadata
                    .get(&format!("{xaxis}_label_tex"))
                    .map_or("", String::as_str),
                if xunit.is_empty() {
                    String::new()
                } else {
                    format!(" [\\si{{{xunit}}}]")
                }
            );
            let yunit = metadata.get("y_unit").map_or("", String::as_str);
            let ylabel = format!(
                "{}{}",
                metadata.get("y_label_tex").map_or("", String::as_str),
                if yunit.is_empty() {
                    String::new()
                } else {
                    format!(" [\\si{{{yunit}}}]")
                }
            );
            let xlog = !xunit.is_empty();
            let ylog = xlog;
            let title = metadata.get("description").map_or("", String::as_str);
            let bins = grid.bin_info().bins();
            let nconvs = self.conv_funs.len();

            let enable_int = bins == 1;
            let enable_abs = !enable_int;
            // TODO: only enable if there are EW corrections
            let enable_rel_ewonoff = enable_abs;
            let enable_abs_pdfs = !(nconvs == 1 || bins == 1);
            let enable_ratio_pdf = enable_abs_pdfs;
            let enable_double_ratio_pdf = enable_abs_pdfs;
            let enable_rel_pdfunc = !(nconvs == 1 || bins == 1 || self.no_conv_fun_unc);
            let enable_rel_pdfpull = enable_rel_pdfunc;

            let config = format!(
                "title = r\"{title}\"
xlabel = r\"{xlabel}\"
ylabel = r\"{ylabel}\"
xlog = {xlog}
ylog = {ylog}
scales = {scales}
plot_panels = {{
    \"int\": {enable_int},
    \"abs\": {enable_abs},
    \"rel_ewonoff\": {enable_rel_ewonoff},
    \"abs_pdfs\": {enable_abs_pdfs},
    \"ratio_pdf\": {enable_ratio_pdf},
    \"double_ratio_pdf\": {enable_double_ratio_pdf},
    \"rel_pdfunc\": {enable_rel_pdfunc},
    \"rel_pdfpull\": {enable_rel_pdfpull},
}}
output = r\"{output}\"",
                enable_int = map_bool(enable_int),
                enable_abs = map_bool(enable_abs),
                enable_rel_ewonoff = map_bool(enable_rel_ewonoff),
                enable_abs_pdfs = map_bool(enable_abs_pdfs),
                enable_ratio_pdf = map_bool(enable_ratio_pdf),
                enable_double_ratio_pdf = map_bool(enable_double_ratio_pdf),
                enable_rel_pdfunc = map_bool(enable_rel_pdfunc),
                enable_rel_pdfpull = map_bool(enable_rel_pdfpull),
                xlabel = xlabel,
                ylabel = ylabel,
                xlog = map_bool(xlog),
                ylog = map_bool(ylog),
                title = title,
                scales = self.scales,
                output = output.to_str().unwrap(),
            );

            let data = format!(
                "data = {data_string}
metadata = {{
{metadata}
}}",
                data_string = data_string,
                metadata = format_metadata(&vector),
            );
            let template = include_str!("plot.py");
            // UNWRAP: if there are no markers the template is wrong
            let config_marker_begin = template.find(MARKER_CONFIG_BEGIN).unwrap();
            let config_marker_end = template.find(MARKER_CONFIG_END).unwrap();
            let data_marker = template.find(MARKER_DATA_INSERT).unwrap();
            // echo template and dynamic content
            print!("{}", template.get(0..config_marker_begin).unwrap());
            print!("{}", config);
            print!(
                "{}",
                template
                    .get((config_marker_end + MARKER_CONFIG_END.len())..data_marker)
                    .unwrap()
            );
            print!("{}", data);
            print!(
                "{}",
                template
                    .get((data_marker + MARKER_DATA_INSERT.len())..)
                    .unwrap()
            );
        } else {
            // TODO: enforce two arguments with clap
            assert_eq!(self.conv_funs.len(), 2);

            let (set1, mut conv_funs1) =
                helpers::create_conv_funs_for_set(&self.conv_funs[0], self.conv_fun_uncert_from)?;
            let (set2, mut conv_funs2) =
                helpers::create_conv_funs_for_set(&self.conv_funs[1], self.conv_fun_uncert_from)?;
            let index @ (order, bin, channel) = self
                .subgrid_pull
                .iter()
                .map(|num| num.parse::<usize>().unwrap())
                .collect_tuple()
                .unwrap();

            let cl = lhapdf::CL_1_SIGMA;
            let grid = helpers::read_grid(&self.input)?;

            // TODO: convert this into an error
            assert_eq!(grid.convolutions().len(), 2);

            let member1 = self.conv_funs[0].members[self.conv_fun_uncert_from];
            let member2 = self.conv_funs[1].members[self.conv_fun_uncert_from];

            let values1: Vec<_> = conv_funs1
                .par_iter_mut()
                .map(|conv_funs| {
                    let values = helpers::convolve(
                        &grid,
                        conv_funs,
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
            let values2: Vec<_> = conv_funs2
                .par_iter_mut()
                .map(|conv_funs| {
                    let values = helpers::convolve(
                        &grid,
                        conv_funs,
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

            // TODO: if no member is given, the zeroth is used, but we should show the averaged
            // result of all members instead
            let res1 = helpers::convolve_subgrid(
                &grid,
                &mut conv_funs1[member1.unwrap_or(0)],
                order,
                bin,
                channel,
                cfg,
            )
            .sum_axis(Axis(0));
            let res2 = helpers::convolve_subgrid(
                &grid,
                &mut conv_funs2[member2.unwrap_or(0)],
                order,
                bin,
                channel,
                cfg,
            )
            .sum_axis(Axis(0));

            let subgrid = &grid.subgrids()[<[usize; 3]>::from(index)];
            //let q2 = subgrid.q2_grid();
            let x1 = grid
                .kinematics()
                .iter()
                .zip(subgrid.node_values())
                .find_map(|(kin, node_values)| {
                    matches!(kin, &Kinematics::X(idx) if idx == 0).then_some(node_values)
                })
                // TODO: convert this into an error
                .unwrap();

            let x2 = grid
                .kinematics()
                .iter()
                .zip(subgrid.node_values())
                .find_map(|(kin, node_values)| {
                    matches!(kin, &Kinematics::X(idx) if idx == 1).then_some(node_values)
                })
                // TODO: convert this into an error
                .unwrap();

            let mut x1_vals = vec![];
            let mut x2_vals = vec![];
            let mut vals = vec![];

            for (((ix1, ix2), &one), &two) in res1.indexed_iter().zip(res2.iter()) {
                if one == 0.0 {
                    assert_approx_eq!(f64, two, 0.0, ulps = 4);
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
