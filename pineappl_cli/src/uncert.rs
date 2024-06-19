use super::helpers::{self, ConvFuns, ConvoluteMode};
use super::{GlobalConfiguration, Subcommand};
use anyhow::{Error, Result};
use clap::builder::{PossibleValuesParser, TypedValueParser};
use clap::{Args, Parser, ValueHint};
use prettytable::{cell, Row};
use rayon::{prelude::*, ThreadPoolBuilder};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::process::ExitCode;
use std::thread;

#[derive(Args)]
#[group(multiple = true, required = true)]
struct Group {
    /// Calculate convolution function uncertainties.
    #[arg(
        default_missing_value = "0",
        num_args = 0..=1,
        long,
        require_equals = true,
        value_delimiter = ',',
        value_name = "IDX"
    )]
    conv_fun: Vec<usize>,
    /// Show absolute numbers of the scale-varied results.
    #[arg(
        default_missing_value = "7",
        num_args = 0..=1,
        long,
        require_equals = true,
        value_name = "SCALES",
        value_parser = PossibleValuesParser::new(["3", "7", "9"]).try_map(|s| s.parse::<u16>())
    )]
    scale_abs: Option<u16>,
    /// Calculate scale uncertainties using the covariance method.
    #[arg(
        default_missing_value = "7",
        num_args = 0..=1,
        long,
        require_equals = true,
        value_name = "SCALES",
        value_parser = PossibleValuesParser::new(["3", "7", "9"]).try_map(|s| s.parse::<u16>())
    )]
    scale_cov: Option<u16>,
    /// Calculate the envelope of results where renormalization and factorization scales varied.
    #[arg(
        default_missing_value = "7",
        num_args = 0..=1,
        long,
        require_equals = true,
        value_name = "SCALES",
        value_parser = PossibleValuesParser::new(["3", "7", "9"]).try_map(|s| s.parse::<u16>())
    )]
    scale_env: Option<u16>,
}

/// Calculate scale and convolution function uncertainties.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF ID(s) or name(s) of the PDF(s)/FF(s).
    conv_funs: ConvFuns,
    #[command(flatten)]
    group: Group,
    /// Confidence level in per cent, for convolution function uncertainties.
    #[arg(default_value_t = lhapdf::CL_1_SIGMA, long)]
    cl: f64,
    /// Show integrated numbers (without bin widths) instead of differential ones.
    #[arg(long, short)]
    integrated: bool,
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
    /// Set the number of fractional digits shown for absolute numbers.
    #[arg(default_value_t = 7, long, value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[arg(default_value_t = 2, long, value_name = "REL")]
    digits_rel: usize,
}

impl Subcommand for Opts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        let grid = helpers::read_grid(&self.input)?;
        let mut conv_funs = helpers::create_conv_funs(&self.conv_funs)?;

        let limits = helpers::convolve_limits(
            &grid,
            &[],
            if self.integrated {
                ConvoluteMode::Integrated
            } else {
                ConvoluteMode::Normal
            },
        );

        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

        let conv_fun_results: Vec<Vec<_>> = self
            .group
            .conv_fun
            .iter()
            .map(|&index| {
                let set = conv_funs[index].set();
                let results: Vec<_> = set
                    .mk_pdfs()?
                    .into_par_iter()
                    .map(|fun| {
                        // TODO: do not create objects that are getting overwritten in any case
                        let mut conv_funs = helpers::create_conv_funs(&self.conv_funs)?;
                        conv_funs[index] = fun;

                        Ok::<_, Error>(helpers::convolve(
                            &grid,
                            &mut conv_funs,
                            &self.orders,
                            &[],
                            &[],
                            1,
                            if self.integrated {
                                ConvoluteMode::Integrated
                            } else {
                                ConvoluteMode::Normal
                            },
                            cfg,
                        ))
                    })
                    .collect::<Result<_, _>>()?;

                // transpose results
                let results: Vec<Vec<_>> = (0..results[0].len())
                    .map(|bin| (0..results.len()).map(|pdf| results[pdf][bin]).collect())
                    .collect();

                results
                    .into_iter()
                    .map(|values| Ok(set.uncertainty(&values, self.cl, false)?))
                    .collect::<Result<_, _>>()
            })
            .collect::<Result<_, Error>>()?;
        let scales_max = self
            .group
            .scale_env
            .iter()
            .chain(self.group.scale_abs.iter())
            .chain(self.group.scale_cov.iter())
            .map(|&x| usize::from(x))
            .max()
            .unwrap_or(1);
        let scale_results = helpers::convolve(
            &grid,
            &mut conv_funs,
            &self.orders,
            &[],
            &[],
            scales_max,
            if self.integrated {
                ConvoluteMode::Integrated
            } else {
                ConvoluteMode::Normal
            },
            cfg,
        );

        let (x, y_label, y_unit) = helpers::labels_and_units(&grid, self.integrated);
        let mut title = Row::empty();
        title.add_cell(cell!(c->"b"));
        for (x_label, x_unit) in x {
            let mut cell = cell!(c->format!("{x_label}\n[{x_unit}]"));
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        title.add_cell(cell!(c->format!("{y_label}\n[{y_unit}]")));

        for &index in &self.group.conv_fun {
            title.add_cell(
                cell!(c->format!("{}", self.conv_funs.lhapdf_names[index])).with_hspan(3),
            );
        }

        if let Some(scales) = self.group.scale_abs {
            for scale in &helpers::SCALES_VECTOR[0..scales.into()] {
                title.add_cell(cell!(c->format!("(r={},f={})\n[{}]", scale.0, scale.1, y_unit)));
            }
        }

        if let Some(scales) = self.group.scale_cov {
            title.add_cell(cell!(c->format!("{}pt scale (cov)\n[%]", scales)).with_hspan(2));
        }

        if let Some(scales) = self.group.scale_env {
            title.add_cell(cell!(c->format!("{}pt-svar (env)\n[%]", scales)).with_hspan(2));
        }

        let mut table = helpers::create_table();
        table.set_titles(title);

        for (bin, (left_right_limits, scale_res)) in limits
            .iter()
            .zip(scale_results.chunks_exact(scales_max))
            .enumerate()
        {
            let row = table.add_empty_row();
            row.add_cell(cell!(r->format!("{bin}")));
            for (left, right) in left_right_limits {
                row.add_cell(cell!(r->format!("{left}")));
                row.add_cell(cell!(r->format!("{right}")));
            }

            row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, scale_res[0])));

            for uncertainty in conv_fun_results.iter().map(|results| &results[bin]) {
                row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, uncertainty.central)));
                row.add_cell(cell!(r->format!("{:.*}", self.digits_rel, -100.0 * uncertainty.errminus / uncertainty.central)));
                row.add_cell(cell!(r->format!("{:.*}", self.digits_rel, 100.0 * uncertainty.errplus / uncertainty.central)));
            }

            if let Some(scales) = self.group.scale_abs {
                for result in scale_res.iter().take(scales.into()) {
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result)));
                }
            }

            if let Some(scales) = self.group.scale_cov {
                let ns = if scales == 3 { 1.0 } else { 2.0 } / f64::from(scales - 1);
                let unc = (ns
                    * scale_res
                        .iter()
                        .take(scales.into())
                        .skip(1)
                        .map(|x| (x - scale_res[0]).powi(2))
                        .sum::<f64>())
                .sqrt();
                let rel_unc = 100.0 * unc / scale_res[0];

                row.add_cell(cell!(r->format!("{:.*}", self.digits_rel, -rel_unc)));
                row.add_cell(cell!(r->format!("{:.*}", self.digits_rel, rel_unc)));
            }

            if let Some(scales) = self.group.scale_env {
                let min_value = scale_res
                    .iter()
                    .take(usize::from(scales))
                    .min_by(|left, right| left.total_cmp(right))
                    .unwrap();
                let max_value = scale_res
                    .iter()
                    .take(usize::from(scales))
                    .max_by(|left, right| left.total_cmp(right))
                    .unwrap();
                let scale_neg = 100.0 * (min_value / scale_res[0] - 1.0);
                let scale_pos = 100.0 * (max_value / scale_res[0] - 1.0);

                row.add_cell(cell!(r->format!("{:.*}", self.digits_rel, scale_neg)));
                row.add_cell(cell!(r->format!("{:.*}", self.digits_rel, scale_pos)));
            }
        }

        table.printstd();

        Ok(ExitCode::SUCCESS)
    }
}
