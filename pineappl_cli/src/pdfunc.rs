use super::helpers::{self, ConvoluteMode, GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use prettytable::{cell, Row};
use rayon::{prelude::*, ThreadPoolBuilder};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::process::ExitCode;
use std::thread;

/// Calculates PDF uncertainties.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id or name of the PDF set.
    #[arg(value_parser = helpers::parse_pdfset)]
    pdfset: String,
    /// Confidence level in per cent.
    #[arg(default_value_t = lhapdf::CL_1_SIGMA, long)]
    cl: f64,
    /// Show integrated numbers (without bin widths) instead of differential ones.
    #[arg(long, short)]
    integrated: bool,
    /// Select orders manually.
    #[arg(
        long,
        num_args(1),
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
        let (set, member) = helpers::create_pdfset(&self.pdfset)?;
        let pdfs = set.mk_pdfs();

        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

        let limits = helpers::convolute_limits(
            &grid,
            &[],
            if self.integrated {
                ConvoluteMode::Integrated
            } else {
                ConvoluteMode::Normal
            },
        );
        let results: Vec<_> = pdfs
            .into_par_iter()
            .flat_map(|mut pdf| {
                helpers::convolute(
                    &grid,
                    &mut pdf,
                    &self.orders,
                    &[],
                    &[],
                    1,
                    if self.integrated {
                        ConvoluteMode::Integrated
                    } else {
                        ConvoluteMode::Normal
                    },
                    cfg.force_positive,
                )
            })
            .collect();

        let (x, y_label, y_unit) = helpers::labels_and_units(&grid, self.integrated);
        let mut title = Row::empty();
        title.add_cell(cell!(c->"b"));
        for (x_label, x_unit) in x {
            let mut cell = cell!(c->format!("{x_label}\n[{x_unit}]"));
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        title.add_cell(cell!(c->format!("{y_label}\n[{y_unit}]")));
        title.add_cell(cell!(c->"PDF uncertainty\n[%]").with_hspan(2));

        let mut table = helpers::create_table();
        table.set_titles(title);

        for (bin, left_right_limits) in limits.iter().enumerate() {
            let values: Vec<_> = results
                .iter()
                .skip(bin)
                .step_by(limits.len())
                .copied()
                .collect();
            let uncertainty = set.uncertainty(&values, self.cl, false)?;

            let row = table.add_empty_row();

            row.add_cell(cell!(r->format!("{bin}")));
            for (left, right) in left_right_limits {
                row.add_cell(cell!(r->format!("{left}")));
                row.add_cell(cell!(r->format!("{right}")));
            }
            row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, member.map_or(uncertainty.central, |member| values[member]))));
            row.add_cell(
                cell!(r->format!("{:.*}", self.digits_rel, (-uncertainty.errminus / uncertainty.central) * 100.0)),
            );
            row.add_cell(
                cell!(r->format!("{:.*}", self.digits_rel, (uncertainty.errplus / uncertainty.central) * 100.0)),
            );
        }

        table.printstd();

        Ok(ExitCode::SUCCESS)
    }
}
