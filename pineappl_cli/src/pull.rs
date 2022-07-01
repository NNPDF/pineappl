use super::helpers::{self, Subcommand};
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
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
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
    fn run(&self) -> Result<()> {
        let grid = helpers::read_grid(&self.input)?;

        let set1 = helpers::create_pdfset(&self.pdfset1)?;
        let set2 = helpers::create_pdfset(&self.pdfset2)?;
        let mut pdfset1 = set1.mk_pdfs();
        let mut pdfset2 = set2.mk_pdfs();

        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

        let results1: Vec<f64> = pdfset1
            .par_iter_mut()
            .flat_map(|pdf| {
                helpers::convolute(&grid, pdf, &[], &[], &[], 1, false, self.force_positive)
            })
            .collect();
        let results2: Vec<f64> = pdfset2
            .par_iter_mut()
            .flat_map(|pdf| {
                helpers::convolute(&grid, pdf, &[], &[], &[], 1, false, self.force_positive)
            })
            .collect();

        let bin_info = grid.bin_info();
        let left_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.left(i))
            .collect();
        let right_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.right(i))
            .collect();

        let mut title = Row::empty();
        title.add_cell(cell!(c->"b"));
        for (x_label, x_unit) in helpers::labels_and_units(&grid, false).0 {
            let mut cell = cell!(c->format!("{}\n[{}]", x_label, x_unit));
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        title.add_cell(cell!(c->"total\n[\u{3c3}]"));
        for _ in 0..self.limit {
            title.add_cell(cell!(c->"l"));
            title.add_cell(cell!(c->"pull\n[\u{3c3}]"));
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
            let uncertainty1 = set1.uncertainty(&values1, self.cl, false)?;
            let uncertainty2 = set2.uncertainty(&values2, self.cl, false)?;

            let lumi_results1: Vec<_> = (0..grid.lumi().len())
                .map(|lumi| {
                    let mut lumi_mask = vec![false; grid.lumi().len()];
                    lumi_mask[lumi] = true;
                    let central: Vec<f64> = pdfset1
                        .par_iter_mut()
                        .map(|pdf| {
                            helpers::convolute(
                                &grid,
                                pdf,
                                &[],
                                &[bin],
                                &lumi_mask,
                                1,
                                false,
                                self.force_positive,
                            )[0]
                        })
                        .collect();
                    set1.uncertainty(&central, self.cl, false).unwrap().central
                })
                .collect();
            let lumi_results2: Vec<_> = (0..grid.lumi().len())
                .map(|lumi| {
                    let mut lumi_mask = vec![false; grid.lumi().len()];
                    lumi_mask[lumi] = true;
                    let central: Vec<f64> = pdfset2
                        .par_iter_mut()
                        .map(|pdf| {
                            helpers::convolute(
                                &grid,
                                pdf,
                                &[],
                                &[bin],
                                &lumi_mask,
                                1,
                                false,
                                self.force_positive,
                            )[0]
                        })
                        .collect();
                    set2.uncertainty(&central, self.cl, false).unwrap().central
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

            row.add_cell(cell!(r->format!("{}", bin)));
            for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                row.add_cell(cell!(r->format!("{}", left[bin])));
                row.add_cell(cell!(r->format!("{}", right[bin])));
            }

            row.add_cell(cell!(r->format!("{:.*}", self.digits, total)));

            // sort using the absolute value in descending order
            pull_tuples.sort_unstable_by(|(_, pull_left), (_, pull_right)| {
                pull_right.abs().partial_cmp(&pull_left.abs()).unwrap()
            });

            for (lumi, pull) in pull_tuples.iter().take(self.limit) {
                row.add_cell(cell!(r->format!("{}", lumi)));
                row.add_cell(cell!(r->format!("{:.*}", self.digits, pull)));
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
        --digits <DIGITS>      Set the number of digits shown for numerical values [default: 3]
        --force-positive       Forces negative PDF values to zero
    -h, --help                 Print help information
    -l, --limit <LIMIT>        The maximum number of luminosities displayed [default: 10]
        --threads <THREADS>    Number of threads to utilize";

    const DEFAULT_STR: &str = "b   etal    total l pull  l  pull  l  pull  l  pull  l  pull 
     []      [\u{3c3}]     [\u{3c3}]     [\u{3c3}]      [\u{3c3}]      [\u{3c3}]      [\u{3c3}]  
-+----+----+-----+-+-----+-+------+-+------+-+------+-+------
0    2 2.25 3.577 0 3.762 1 -0.108 3 -0.052 4 -0.016 2 -0.009
1 2.25  2.5 3.451 0 3.630 1 -0.095 3 -0.062 4 -0.016 2 -0.006
2  2.5 2.75 3.195 0 3.338 1 -0.072 3 -0.056 4 -0.010 2 -0.005
3 2.75    3 2.803 0 2.887 1 -0.045 3 -0.024 4 -0.011 2 -0.004
4    3 3.25 2.346 0 2.349 3  0.023 1 -0.013 4 -0.009 2 -0.004
5 3.25  3.5 1.873 0 1.810 3  0.082 2 -0.011 4 -0.007 1 -0.001
6  3.5    4 1.468 0 1.389 3  0.177 1 -0.088 4 -0.007 2 -0.003
7    4  4.5 1.219 0 1.439 1 -0.358 3  0.147 4 -0.006 2 -0.001
";

    const CL_90_STR: &str = "b   etal    total l pull  l  pull  l  pull  l  pull  l  pull 
     []      [\u{3c3}]     [\u{3c3}]     [\u{3c3}]      [\u{3c3}]      [\u{3c3}]      [\u{3c3}]  
-+----+----+-----+-+-----+-+------+-+------+-+------+-+------
0    2 2.25 2.175 0 2.287 1 -0.066 3 -0.031 4 -0.009 2 -0.005
1 2.25  2.5 2.098 0 2.207 1 -0.058 3 -0.038 4 -0.010 2 -0.003
2  2.5 2.75 1.942 0 2.029 1 -0.044 3 -0.034 4 -0.006 2 -0.003
3 2.75    3 1.704 0 1.755 1 -0.027 3 -0.015 4 -0.007 2 -0.003
4    3 3.25 1.426 0 1.428 3  0.014 1 -0.008 4 -0.005 2 -0.003
5 3.25  3.5 1.138 0 1.100 3  0.050 2 -0.007 4 -0.004 1 -0.001
6  3.5    4 0.892 0 0.845 3  0.107 1 -0.053 4 -0.004 2 -0.002
7    4  4.5 0.741 0 0.875 1 -0.218 3  0.089 4 -0.004 2 -0.001
";

    const LIMIT_STR: &str = "b   etal    total l pull 
     []      [\u{3c3}]     [\u{3c3}] 
-+----+----+-----+-+-----
0    2 2.25 3.577 0 3.762
1 2.25  2.5 3.451 0 3.630
2  2.5 2.75 3.195 0 3.338
3 2.75    3 2.803 0 2.887
4    3 3.25 2.346 0 2.349
5 3.25  3.5 1.873 0 1.810
6  3.5    4 1.468 0 1.389
7    4  4.5 1.219 0 1.439
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["pull", "--help"])
            .assert()
            .success()
            .stdout(format!("{} [default: {}]\n", HELP_STR, num_cpus::get()));
    }

    #[test]
    fn default() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "pull",
                "--threads=1",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
                "NNPDF40_nnlo_as_01180",
            ])
            .assert()
            .success()
            .stdout(DEFAULT_STR);
    }

    #[test]
    fn cl_90() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "pull",
                "--cl=90",
                "--threads=1",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
                "NNPDF40_nnlo_as_01180",
            ])
            .assert()
            .success()
            .stdout(CL_90_STR);
    }

    #[test]
    fn limit() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "pull",
                "--limit=1",
                "--threads=1",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
                "NNPDF40_nnlo_as_01180",
            ])
            .assert()
            .success()
            .stdout(LIMIT_STR);
    }
}
