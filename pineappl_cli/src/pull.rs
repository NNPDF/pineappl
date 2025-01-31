use super::helpers::{self, ConvFuns, ConvoluteMode};
use super::{GlobalConfiguration, Subcommand};
use anyhow::{Error, Result};
use clap::{Parser, ValueHint};
use lhapdf::{Pdf, PdfSet};
use prettytable::{cell, Row};
use rayon::{prelude::*, ThreadPoolBuilder};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::process::ExitCode;
use std::thread;

// TODO: do we need the CL parameter?

/// Calculates the pull between two different PDF sets.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF ID(s) or name(s) of the first PDF(s)/FF(s).
    conv_funs1: ConvFuns,
    /// LHAPDF ID(s) or name(s) of the second PDF(s)/FF(s).
    conv_funs2: ConvFuns,
    /// Index of the convolution functions for which the pull should be calculated.
    #[arg(default_value = "0", long, value_name = "IDX")]
    pull_from: usize,
    /// Confidence level in per cent.
    #[arg(default_value_t = lhapdf::CL_1_SIGMA, long)]
    cl: f64,
    /// The maximum number of channels displayed.
    #[arg(default_value_t = 10, long, short)]
    limit: usize,
    /// Select orders manually.
    #[arg(
        long,
        num_args = 1,
        short,
        value_delimiter = ',',
        value_parser = helpers::parse_order
    )]
    orders: Vec<(u8, u8)>,
    /// Number of threads to utilize.
    #[arg(default_value_t = thread::available_parallelism().map_or(1, NonZeroUsize::get), long)]
    threads: usize,
    /// Set the number of digits shown for numerical values.
    #[arg(default_value_t = 3, long)]
    digits: usize,
}

impl Subcommand for Opts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        let grid = helpers::read_grid(&self.input)?;

        let (set1, mut conv_funs1) =
            helpers::create_conv_funs_for_set(&self.conv_funs1, self.pull_from)?;
        let (set2, mut conv_funs2) =
            helpers::create_conv_funs_for_set(&self.conv_funs2, self.pull_from)?;

        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

        let limit = grid.channels().len().min(self.limit);
        let bin_limits = helpers::convolve_limits(&grid, &[], ConvoluteMode::Normal);
        let results1: Vec<_> = conv_funs1
            .par_iter_mut()
            .map(|funs| {
                Ok::<_, Error>(helpers::convolve(
                    &grid,
                    funs,
                    &self.conv_funs1.conv_types,
                    &self.orders,
                    &[],
                    &[],
                    1,
                    ConvoluteMode::Normal,
                    cfg,
                ))
            })
            .collect::<Result<_, _>>()?;
        let results1: Vec<Vec<_>> = (0..results1[0].len())
            .map(|bin| (0..results1.len()).map(|pdf| results1[pdf][bin]).collect())
            .collect();
        let results2: Vec<_> = conv_funs2
            .par_iter_mut()
            .map(|funs| {
                Ok::<_, Error>(helpers::convolve(
                    &grid,
                    funs,
                    &self.conv_funs2.conv_types,
                    &self.orders,
                    &[],
                    &[],
                    1,
                    ConvoluteMode::Normal,
                    cfg,
                ))
            })
            .collect::<Result<_, _>>()?;
        let results2: Vec<Vec<_>> = (0..results2[0].len())
            .map(|bin| (0..results2.len()).map(|pdf| results2[pdf][bin]).collect())
            .collect();

        let mut title = Row::empty();
        title.add_cell(cell!(c->"b"));
        for (x_label, x_unit) in helpers::labels_and_units(&grid, false).0 {
            let mut cell = cell!(c->format!("{x_label}\n[{x_unit}]"));
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        title.add_cell(cell!(c->"total\n[\u{3c3}]"));
        for _ in 0..limit {
            title.add_cell(cell!(c->"c"));
            title.add_cell(cell!(c->"pull\n[\u{3c3}]"));
        }

        let mut table = helpers::create_table();
        table.set_titles(title);

        for (bin, limits) in bin_limits.iter().enumerate() {
            let (total, unc1, unc2) = {
                let values1 = &results1[bin];
                let values2 = &results2[bin];
                let uncertainty1 = set1.uncertainty(values1, self.cl, false)?;
                let uncertainty2 = set2.uncertainty(values2, self.cl, false)?;

                // if requested use the given member instead of the central value
                let diff = self.conv_funs2.members[self.pull_from]
                    .map_or(uncertainty2.central, |member| values2[member])
                    - self.conv_funs1.members[self.pull_from]
                        .map_or(uncertainty1.central, |member| values1[member]);

                // use the uncertainties in the direction in which they point to each other
                let (unc1, unc2) = if diff > 0.0 {
                    (uncertainty2.errminus, uncertainty1.errplus)
                } else {
                    (uncertainty1.errminus, uncertainty2.errplus)
                };
                (diff / unc1.hypot(unc2), unc1, unc2)
            };

            let channel_results =
                |conv_funs: &ConvFuns, pdfset: &mut [Vec<Pdf>], set: &PdfSet| -> Vec<f64> {
                    if let Some(member) = conv_funs.members[self.pull_from] {
                        (0..grid.channels().len())
                            .map(|channel| {
                                let mut channel_mask = vec![false; grid.channels().len()];
                                channel_mask[channel] = true;
                                match helpers::convolve(
                                    &grid,
                                    &mut pdfset[member],
                                    &conv_funs.conv_types,
                                    &self.orders,
                                    &[bin],
                                    &channel_mask,
                                    1,
                                    ConvoluteMode::Normal,
                                    cfg,
                                )
                                .as_slice()
                                {
                                    [value] => *value,
                                    _ => unreachable!(),
                                }
                            })
                            .collect()
                    } else {
                        let results: Vec<_> = pdfset
                            .iter_mut()
                            .flat_map(|fun| {
                                (0..grid.channels().len())
                                    .map(|channel| {
                                        let mut channel_mask = vec![false; grid.channels().len()];
                                        channel_mask[channel] = true;
                                        match helpers::convolve(
                                            &grid,
                                            fun,
                                            &conv_funs.conv_types,
                                            &self.orders,
                                            &[bin],
                                            &channel_mask,
                                            1,
                                            ConvoluteMode::Normal,
                                            cfg,
                                        )
                                        .as_slice()
                                        {
                                            [value] => *value,
                                            _ => unreachable!(),
                                        }
                                    })
                                    .collect::<Vec<_>>()
                            })
                            .collect();

                        (0..grid.channels().len())
                            .map(|channel| {
                                let central: Vec<_> = results
                                    .iter()
                                    .skip(channel)
                                    .step_by(grid.channels().len())
                                    .copied()
                                    .collect();
                                set.uncertainty(&central, self.cl, false).unwrap().central
                            })
                            .collect()
                    }
                };

            let mut pull_tuples = if self.limit == 0 {
                Vec::new()
            } else {
                let channel_results1 = channel_results(&self.conv_funs1, &mut conv_funs1, &set1);
                let channel_results2 = channel_results(&self.conv_funs2, &mut conv_funs2, &set2);

                let pull_tuples: Vec<_> = channel_results2
                    .iter()
                    .zip(channel_results1.iter())
                    .map(|(res2, res1)| (res2 - res1) / unc1.hypot(unc2))
                    .enumerate()
                    .collect();

                pull_tuples
            };

            let row = table.add_empty_row();

            row.add_cell(cell!(r->format!("{bin}")));
            for (left, right) in limits {
                row.add_cell(cell!(r->format!("{left}")));
                row.add_cell(cell!(r->format!("{right}")));
            }

            row.add_cell(cell!(r->format!("{:.*}", self.digits, total)));

            // sort using the absolute value in descending order
            pull_tuples.sort_unstable_by(|(_, pull_left), (_, pull_right)| {
                pull_right.abs().total_cmp(&pull_left.abs())
            });

            for (channel, pull) in pull_tuples.iter().take(self.limit) {
                row.add_cell(cell!(r->format!("{channel}")));
                row.add_cell(cell!(r->format!("{:.*}", self.digits, pull)));
            }
        }

        table.printstd();

        Ok(ExitCode::SUCCESS)
    }
}
