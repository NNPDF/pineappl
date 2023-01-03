use super::helpers::{self, ConvoluteMode, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use prettytable::{cell, Row};
use rayon::{prelude::*, ThreadPoolBuilder};
use std::path::PathBuf;

/// Calculates PDF uncertainties.
#[derive(Parser)]
#[clap(aliases = &["pdf-uncertainty", "pdf_uncertainty"])]
pub struct Opts {
    /// Path to the input grid.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id or name of the PDF set.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset: String,
    /// Confidence level in per cent.
    #[clap(default_value_t = lhapdf::CL_1_SIGMA, long)]
    cl: f64,
    /// Show integrated numbers (without bin widths) instead of differential ones.
    #[clap(long, short)]
    integrated: bool,
    /// Select orders manually.
    #[clap(
        long,
        min_values = 1,
        parse(try_from_str = helpers::parse_order),
        short,
        use_value_delimiter = true
    )]
    orders: Vec<(u32, u32)>,
    /// Number of threads to utilize.
    #[clap(default_value_t = num_cpus::get(), long)]
    threads: usize,
    /// Set the number of fractional digits shown for absolute numbers.
    #[clap(default_value_t = 7, long = "digits-abs", value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[clap(default_value_t = 2, long = "digits-rel", value_name = "REL")]
    digits_rel: usize,
    /// Forces negative PDF values to zero.
    #[clap(long = "force-positive")]
    force_positive: bool,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<u8> {
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
        let results: Vec<f64> = pdfs
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
                    self.force_positive,
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

        Ok(0)
    }
}
