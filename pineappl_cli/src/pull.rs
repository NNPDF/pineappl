use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use lhapdf::PdfSet;
use prettytable::{cell, Row};
use rayon::{prelude::*, ThreadPoolBuilder};
use std::path::PathBuf;

/// Calculates the pull between two different PDF sets.
#[derive(Parser)]
#[clap(name = "pull")]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id or name of the first PDF set.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset1: String,
    /// LHAPDF id or name of the second PDF set.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset2: String,
    /// Confidence level in per cent.
    #[clap(default_value = helpers::ONE_SIGMA_STR, long)]
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
    #[clap(default_value = &helpers::NUM_CPUS_STRING, long)]
    threads: usize,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let grid = helpers::read_grid(&self.input)?;

        let set1 = PdfSet::new(&self.pdfset1.parse().map_or_else(
            |_| self.pdfset1.to_string(),
            |lhaid| lhapdf::lookup_pdf(lhaid).unwrap().0,
        ));
        let set2 = PdfSet::new(&self.pdfset2.parse().map_or_else(
            |_| self.pdfset2.to_string(),
            |lhaid| lhapdf::lookup_pdf(lhaid).unwrap().0,
        ));
        let pdfset1 = set1.mk_pdfs();
        let pdfset2 = set2.mk_pdfs();

        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

        let results1: Vec<f64> = pdfset1
            .par_iter()
            .flat_map(|pdf| helpers::convolute(&grid, pdf, &[], &[], &[], 1))
            .collect();
        let results2: Vec<f64> = pdfset2
            .par_iter()
            .flat_map(|pdf| helpers::convolute(&grid, pdf, &[], &[], &[], 1))
            .collect();

        let bin_info = grid.bin_info();
        let left_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.left(i))
            .collect();
        let right_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.right(i))
            .collect();

        let labels = helpers::labels(&grid);
        let (_, x_labels) = labels.split_last().unwrap();
        let mut title = Row::empty();
        title.add_cell(cell!(c->"bin"));
        for x_label in x_labels {
            let mut cell = cell!(c->&x_label);
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        title.add_cell(cell!(c->"total"));
        for _ in 0..self.limit {
            title.add_cell(cell!(c->"lumi"));
            title.add_cell(cell!(c->"pull"));
        }

        let mut table = helpers::create_table();
        table.set_titles(title);

        for bin in 0..bin_info.bins() {
            let values1: Vec<_> = results1
                .iter()
                .skip(bin)
                .step_by(bin_info.bins())
                .copied()
                .collect();
            let values2: Vec<_> = results2
                .iter()
                .skip(bin)
                .step_by(bin_info.bins())
                .copied()
                .collect();
            let uncertainty1 = set1.uncertainty(&values1, self.cl, false);
            let uncertainty2 = set2.uncertainty(&values2, self.cl, false);

            let lumi_results1: Vec<_> = (0..grid.lumi().len())
                .map(|lumi| {
                    let mut lumi_mask = vec![false; grid.lumi().len()];
                    lumi_mask[lumi] = true;
                    let central: Vec<f64> = pdfset1
                        .iter()
                        .map(|pdf| helpers::convolute(&grid, pdf, &[], &[bin], &lumi_mask, 1)[0])
                        .collect();
                    set1.uncertainty(&central, self.cl, false).central
                })
                .collect();
            let lumi_results2: Vec<_> = (0..grid.lumi().len())
                .map(|lumi| {
                    let mut lumi_mask = vec![false; grid.lumi().len()];
                    lumi_mask[lumi] = true;
                    let central: Vec<f64> = pdfset2
                        .iter()
                        .map(|pdf| helpers::convolute(&grid, pdf, &[], &[bin], &lumi_mask, 1)[0])
                        .collect();
                    set1.uncertainty(&central, self.cl, false).central
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

            row.add_cell(cell!(r->&format!("{}", bin)));
            for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                row.add_cell(cell!(r->&format!("{}", left[bin])));
                row.add_cell(cell!(r->&format!("{}", right[bin])));
            }

            row.add_cell(cell!(r->&format!("{:.3}", total)));

            // sort using the absolute value in descending order
            pull_tuples.sort_unstable_by(|(_, pull_left), (_, pull_right)| {
                pull_right.abs().partial_cmp(&pull_left.abs()).unwrap()
            });

            for (lumi, pull) in pull_tuples.iter().take(self.limit) {
                row.add_cell(cell!(r->&format!("#{}", lumi)));
                row.add_cell(cell!(r->&format!("{:.3}", pull)));
            }
        }

        table.printstd();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-pull 
Calculates the pull between two different PDF sets

USAGE:
    pineappl pull [OPTIONS] <INPUT> <PDFSET1> <PDFSET2>

ARGS:
    <INPUT>      Path to the input grid
    <PDFSET1>    LHAPDF id or name of the first PDF set
    <PDFSET2>    LHAPDF id or name of the second PDF set

OPTIONS:
        --cl <CL>              Confidence level in per cent [default: 68.26894921370858]
    -h, --help                 Print help information
    -l, --limit <LIMIT>        The maximum number of luminosities displayed [default: 10]
        --threads <THREADS>    Number of threads to utilize";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["pull", "--help"])
            .assert()
            .success()
            .stdout(format!("{} [default: {}]\n", HELP_STR, num_cpus::get()));
    }
}
