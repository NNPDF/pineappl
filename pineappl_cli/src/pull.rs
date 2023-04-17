use super::helpers::{self, ConvoluteMode, GlobalConfiguration, Subcommand};
use anyhow::Result;
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
    /// LHAPDF id or name of the first PDF set.
    #[arg(value_parser = helpers::parse_pdfset)]
    pdfset1: String,
    /// LHAPDF id or name of the second PDF set.
    #[arg(value_parser = helpers::parse_pdfset)]
    pdfset2: String,
    /// Confidence level in per cent.
    #[arg(default_value_t = lhapdf::CL_1_SIGMA, long)]
    cl: f64,
    /// The maximum number of luminosities displayed.
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
    orders: Vec<(u32, u32)>,
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

        let (set1, member1) = helpers::create_pdfset(&self.pdfset1)?;
        let (set2, member2) = helpers::create_pdfset(&self.pdfset2)?;
        let mut pdfset1 = set1.mk_pdfs();
        let mut pdfset2 = set2.mk_pdfs();

        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

        let limit = grid.lumi().len().min(self.limit);
        let bin_limits = helpers::convolute_limits(&grid, &[], ConvoluteMode::Normal);
        let results1: Vec<_> = pdfset1
            .par_iter_mut()
            .flat_map(|pdf| {
                helpers::convolute(
                    &grid,
                    pdf,
                    &self.orders,
                    &[],
                    &[],
                    1,
                    ConvoluteMode::Normal,
                    cfg.force_positive,
                )
            })
            .collect();
        let results2: Vec<_> = pdfset2
            .par_iter_mut()
            .flat_map(|pdf| {
                helpers::convolute(
                    &grid,
                    pdf,
                    &self.orders,
                    &[],
                    &[],
                    1,
                    ConvoluteMode::Normal,
                    cfg.force_positive,
                )
            })
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
            title.add_cell(cell!(c->"l"));
            title.add_cell(cell!(c->"pull\n[\u{3c3}]"));
        }

        let mut table = helpers::create_table();
        table.set_titles(title);

        for (bin, limits) in bin_limits.iter().enumerate() {
            let (total, unc1, unc2) = {
                let values1: Vec<_> = results1
                    .iter()
                    .skip(bin)
                    .step_by(bin_limits.len())
                    .copied()
                    .collect();
                let values2: Vec<_> = results2
                    .iter()
                    .skip(bin)
                    .step_by(bin_limits.len())
                    .copied()
                    .collect();
                let uncertainty1 = set1.uncertainty(&values1, self.cl, false)?;
                let uncertainty2 = set2.uncertainty(&values2, self.cl, false)?;

                // if requested use the given member instead of the central value
                let diff = member2.map_or(uncertainty2.central, |member| values2[member])
                    - member1.map_or(uncertainty1.central, |member| values1[member]);

                // use the uncertainties in the direction in which they point to each other
                let (unc1, unc2) = if diff > 0.0 {
                    (uncertainty2.errminus, uncertainty1.errplus)
                } else {
                    (uncertainty1.errminus, uncertainty2.errplus)
                };
                (diff / unc1.hypot(unc2), unc1, unc2)
            };

            let lumi_results =
                |member: Option<usize>, pdfset: &mut Vec<Pdf>, set: &PdfSet| -> Vec<f64> {
                    if let Some(member) = member {
                        (0..grid.lumi().len())
                            .map(|lumi| {
                                let mut lumi_mask = vec![false; grid.lumi().len()];
                                lumi_mask[lumi] = true;
                                match helpers::convolute(
                                    &grid,
                                    &mut pdfset[member],
                                    &self.orders,
                                    &[bin],
                                    &lumi_mask,
                                    1,
                                    ConvoluteMode::Normal,
                                    cfg.force_positive,
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
                            .flat_map(|pdf| {
                                (0..grid.lumi().len())
                                    .map(|lumi| {
                                        let mut lumi_mask = vec![false; grid.lumi().len()];
                                        lumi_mask[lumi] = true;
                                        match helpers::convolute(
                                            &grid,
                                            pdf,
                                            &self.orders,
                                            &[bin],
                                            &lumi_mask,
                                            1,
                                            ConvoluteMode::Normal,
                                            cfg.force_positive,
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

                        (0..grid.lumi().len())
                            .map(|lumi| {
                                let central: Vec<_> = results
                                    .iter()
                                    .skip(lumi)
                                    .step_by(grid.lumi().len())
                                    .copied()
                                    .collect();
                                set.uncertainty(&central, self.cl, false).unwrap().central
                            })
                            .collect()
                    }
                };

            let mut pull_tuples = if self.limit == 0 {
                vec![]
            } else {
                let lumi_results1 = lumi_results(member1, &mut pdfset1, &set1);
                let lumi_results2 = lumi_results(member2, &mut pdfset2, &set2);

                let pull_tuples: Vec<_> = lumi_results2
                    .iter()
                    .zip(lumi_results1.iter())
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

            for (lumi, pull) in pull_tuples.iter().take(self.limit) {
                row.add_cell(cell!(r->format!("{lumi}")));
                row.add_cell(cell!(r->format!("{:.*}", self.digits, pull)));
            }
        }

        table.printstd();

        Ok(ExitCode::SUCCESS)
    }
}
