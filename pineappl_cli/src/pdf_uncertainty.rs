use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use lhapdf::PdfSet;
use prettytable::{cell, Row};
use rayon::{prelude::*, ThreadPoolBuilder};
use std::path::PathBuf;

/// Calculates PDF uncertainties.
#[derive(Parser)]
#[clap(name = "pdfunc", aliases = &["pdf-uncertainty", "pdf_uncertainty"])]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id or name of the PDF set.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset: String,
    /// Confidence level in per cent.
    #[clap(default_value = helpers::ONE_SIGMA_STR, long)]
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
        use_delimiter = true
    )]
    orders: Vec<(u32, u32)>,
    /// Number of threads to utilize.
    #[clap(default_value = &helpers::NUM_CPUS_STRING, long)]
    threads: usize,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let grid = helpers::read_grid(&self.input)?;
        let set = PdfSet::new(&self.pdfset.parse().map_or_else(
            |_| self.pdfset.to_string(),
            |lhaid| lhapdf::lookup_pdf(lhaid).unwrap().0,
        ));
        let pdfs = set.mk_pdfs();

        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

        let results: Vec<f64> = pdfs
            .into_par_iter()
            .flat_map(|pdf| helpers::convolute(&grid, &pdf, &self.orders, &[], &[], 1))
            .collect();

        let bin_info = grid.bin_info();
        let left_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.left(i))
            .collect();
        let right_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.right(i))
            .collect();
        let normalizations = bin_info.normalizations();

        let labels = helpers::labels(&grid);
        let (y_label, x_labels) = labels.split_last().unwrap();
        let mut title = Row::empty();
        title.add_cell(cell!(c->"bin"));
        for x_label in x_labels {
            let mut cell = cell!(c->&x_label);
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        title.add_cell(cell!(c->if self.integrated { "integ" } else { y_label }));
        title.add_cell(cell!(c->"neg unc"));
        title.add_cell(cell!(c->"pos unc"));

        let mut table = helpers::create_table();
        table.set_titles(title);

        for bin in 0..bin_info.bins() {
            let values: Vec<_> = results
                .iter()
                .skip(bin)
                .step_by(bin_info.bins())
                .copied()
                .collect();
            let uncertainty = set.uncertainty(&values, self.cl, false);

            let row = table.add_empty_row();

            row.add_cell(cell!(r->&format!("{}", bin)));
            for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                row.add_cell(cell!(r->&format!("{}", left[bin])));
                row.add_cell(cell!(r->&format!("{}", right[bin])));
            }
            row.add_cell(cell!(r->&format!("{:.7e}", if self.integrated { uncertainty.central * normalizations[bin] } else { uncertainty.central })));
            row.add_cell(
                cell!(r->&format!("{:.2}%", (-uncertainty.errminus / uncertainty.central) * 100.0)),
            );
            row.add_cell(
                cell!(r->&format!("{:.2}%", (uncertainty.errplus / uncertainty.central) * 100.0)),
            );
        }

        table.printstd();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-pdf-uncertainty 

Calculates PDF uncertainties

USAGE:
    pineappl pdf-uncertainty [OPTIONS] <INPUT> <PDFSET>

ARGS:
    <INPUT>     Path to the input grid
    <PDFSET>    LHAPDF id or name of the PDF set

OPTIONS:
        --cl <CL>               Confidence level in per cent [default: 68.26894921370858]
    -h, --help                  Print help information
    -i, --integrated            Show integrated numbers (without bin widths) instead of differential
                                ones
    -o, --orders <ORDERS>...    Select orders manually
        --threads <THREADS>     Number of threads to utilize";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["pdf_uncertainty", "--help"])
            .assert()
            .success()
            .stdout(format!("{} [default: {}]\n", HELP_STR, num_cpus::get()));
    }
}
