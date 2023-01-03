use super::helpers::{self, ConvoluteMode, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use prettytable::{cell, Row};
use rayon::{prelude::*, ThreadPoolBuilder};
use std::path::PathBuf;

// TODO: do we need the CL parameter?

/// Calculates the pull between two different PDF sets.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id or name of the first PDF set.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset1: String,
    /// LHAPDF id or name of the second PDF set.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset2: String,
    /// Confidence level in per cent.
    #[clap(default_value_t = lhapdf::CL_1_SIGMA, long)]
    cl: f64,
    /// The maximum number of luminosities displayed.
    #[clap(
        default_value = "10",
        long,
        short,
        validator = helpers::validate_pos_non_zero::<usize>
    )]
    limit: usize,
    /// Number of threads to utilize.
    #[clap(default_value_t = num_cpus::get(), long)]
    threads: usize,
    /// Set the number of digits shown for numerical values.
    #[clap(default_value_t = 3, long = "digits")]
    digits: usize,
    /// Forces negative PDF values to zero.
    #[clap(long = "force-positive")]
    force_positive: bool,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<u8> {
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
        let results1: Vec<f64> = pdfset1
            .par_iter_mut()
            .flat_map(|pdf| {
                helpers::convolute(
                    &grid,
                    pdf,
                    &[],
                    &[],
                    &[],
                    1,
                    ConvoluteMode::Normal,
                    self.force_positive,
                )
            })
            .collect();
        let results2: Vec<f64> = pdfset2
            .par_iter_mut()
            .flat_map(|pdf| {
                helpers::convolute(
                    &grid,
                    pdf,
                    &[],
                    &[],
                    &[],
                    1,
                    ConvoluteMode::Normal,
                    self.force_positive,
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

            let lumi_results1: Vec<_> = (0..grid.lumi().len())
                .map(|lumi| {
                    let mut lumi_mask = vec![false; grid.lumi().len()];
                    lumi_mask[lumi] = true;

                    if let Some(member1) = member1 {
                        match helpers::convolute(
                            &grid,
                            &mut pdfset1[member1],
                            &[],
                            &[bin],
                            &lumi_mask,
                            1,
                            ConvoluteMode::Normal,
                            self.force_positive,
                        )
                        .as_slice()
                        {
                            [value] => *value,
                            _ => unreachable!(),
                        }
                    } else {
                        let central: Vec<f64> = pdfset1
                            .par_iter_mut()
                            .map(|pdf| {
                                match helpers::convolute(
                                    &grid,
                                    pdf,
                                    &[],
                                    &[bin],
                                    &lumi_mask,
                                    1,
                                    ConvoluteMode::Normal,
                                    self.force_positive,
                                )
                                .as_slice()
                                {
                                    [value] => *value,
                                    _ => unreachable!(),
                                }
                            })
                            .collect();
                        set1.uncertainty(&central, self.cl, false).unwrap().central
                    }
                })
                .collect();
            let lumi_results2: Vec<_> = (0..grid.lumi().len())
                .map(|lumi| {
                    let mut lumi_mask = vec![false; grid.lumi().len()];
                    lumi_mask[lumi] = true;

                    if let Some(member2) = member2 {
                        match helpers::convolute(
                            &grid,
                            &mut pdfset2[member2],
                            &[],
                            &[bin],
                            &lumi_mask,
                            1,
                            ConvoluteMode::Normal,
                            self.force_positive,
                        )
                        .as_slice()
                        {
                            [value] => *value,
                            _ => unreachable!(),
                        }
                    } else {
                        let central: Vec<f64> = pdfset2
                            .par_iter_mut()
                            .map(|pdf| {
                                match helpers::convolute(
                                    &grid,
                                    pdf,
                                    &[],
                                    &[bin],
                                    &lumi_mask,
                                    1,
                                    ConvoluteMode::Normal,
                                    self.force_positive,
                                )
                                .as_slice()
                                {
                                    [value] => *value,
                                    _ => unreachable!(),
                                }
                            })
                            .collect();
                        set2.uncertainty(&central, self.cl, false).unwrap().central
                    }
                })
                .collect();

            let mut pull_tuples: Vec<_> = lumi_results2
                .iter()
                .zip(lumi_results1.iter())
                .map(|(res2, res1)| {
                    // use the uncertainties in the direction in which the respective results differ
                    let unc1 = if res1 > res2 {
                        uncertainty1.errminus
                    } else {
                        uncertainty1.errplus
                    };
                    let unc2 = if res2 > res1 {
                        uncertainty2.errminus
                    } else {
                        uncertainty2.errplus
                    };
                    (res2 - res1) / unc1.hypot(unc2)
                })
                .enumerate()
                .collect();

            let total = pull_tuples
                .iter()
                .fold(0.0, |value, (_, pull)| value + pull);

            let row = table.add_empty_row();

            row.add_cell(cell!(r->format!("{bin}")));
            for (left, right) in limits {
                row.add_cell(cell!(r->format!("{left}")));
                row.add_cell(cell!(r->format!("{right}")));
            }

            row.add_cell(cell!(r->format!("{:.*}", self.digits, total)));

            // sort using the absolute value in descending order
            pull_tuples.sort_unstable_by(|(_, pull_left), (_, pull_right)| {
                pull_right.abs().partial_cmp(&pull_left.abs()).unwrap()
            });

            for (lumi, pull) in pull_tuples.iter().take(self.limit) {
                row.add_cell(cell!(r->format!("{lumi}")));
                row.add_cell(cell!(r->format!("{:.*}", self.digits, pull)));
            }
        }

        table.printstd();

        Ok(0)
    }
}
