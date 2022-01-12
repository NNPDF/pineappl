use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use lhapdf::Pdf;
use prettytable::{cell, Row};
use std::ops::RangeInclusive;
use std::path::PathBuf;

/// Convolutes a PineAPPL grid with a PDF set.
#[derive(Parser)]
pub struct Opts {
    /// Path of the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id(s) or name of the PDF set(s).
    #[clap(min_values = 1, validator = helpers::validate_pdfset)]
    pdfsets: Vec<String>,
    /// Show absolute numbers of the scale variation.
    #[clap(long, short)]
    absolute: bool,
    /// Selects a subset of bins.
    #[clap(
        long,
        short,
        multiple_values = true,
        parse(try_from_str = helpers::try_parse_integer_range),
        use_delimiter = true
    )]
    bins: Vec<RangeInclusive<usize>>,
    /// Show integrated numbers (without bin widths) instead of differential ones.
    #[clap(long, short)]
    integrated: bool,
    /// Select orders manually.
    #[clap(
        long,
        multiple_values = true,
        parse(try_from_str = helpers::parse_order),
        short,
        use_delimiter = true
    )]
    orders: Vec<(u32, u32)>,
    /// Set the number of scale variations.
    #[clap(default_value = "7", long, possible_values = ["1", "3", "7", "9"], short)]
    scales: usize,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let grid = helpers::read_grid(&self.input)?;
        let pdf = self.pdfsets[0].parse().map_or_else(
            |_| Pdf::with_setname_and_member(&self.pdfsets[0], 0),
            Pdf::with_lhaid,
        );
        let bins: Vec<_> = self.bins.iter().cloned().flatten().collect();

        let results = helpers::convolute(&grid, &pdf, &self.orders, &bins, &[], self.scales);

        let other_results: Vec<f64> = self.pdfsets[1..]
            .iter()
            .flat_map(|pdfset| {
                let pdf = pdfset
                    .parse()
                    .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);
                helpers::convolute(&grid, &pdf, &[], &bins, &[], 1)
            })
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

        if self.absolute {
            for scale in &helpers::SCALES_VECTOR[0..self.scales] {
                title.add_cell(cell!(c->&format!("({},{})", scale.0, scale.1)));
            }
        } else {
            title.add_cell(cell!(c->"scale uncertainty").with_hspan(2));
        }

        for other in self.pdfsets[1..].iter() {
            let mut cell = cell!(c->other);
            cell.set_hspan(2);
            title.add_cell(cell);
        }

        let mut table = helpers::create_table();
        table.set_titles(title);

        for (index, values) in results.chunks_exact(self.scales).enumerate() {
            let min_value = values
                .iter()
                .min_by(|left, right| left.partial_cmp(right).unwrap())
                .unwrap();
            let max_value = values
                .iter()
                .max_by(|left, right| left.partial_cmp(right).unwrap())
                .unwrap();
            let bin = if bins.is_empty() { index } else { bins[index] };

            let row = table.add_empty_row();

            row.add_cell(cell!(r->&format!("{}", bin)));
            for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                row.add_cell(cell!(r->&format!("{}", left[bin])));
                row.add_cell(cell!(r->&format!("{}", right[bin])));
            }
            row.add_cell(cell!(r->&format!("{:.7e}", if self.integrated { values[0] * normalizations[bin] } else { values[0] })));

            if self.absolute {
                for &value in values.iter() {
                    row.add_cell(cell!(r->&format!("{:.7e}", if self.integrated { value * normalizations[bin] } else { value })));
                }
            } else {
                row.add_cell(cell!(r->&format!("{:.2}%", (min_value / values[0] - 1.0) * 100.0)));
                row.add_cell(cell!(r->&format!("{:.2}%", (max_value / values[0] - 1.0) * 100.0)));
            }

            let bins = if bins.is_empty() {
                bin_info.bins()
            } else {
                self.bins.len()
            };

            for &other in other_results.iter().skip(index).step_by(bins) {
                row.add_cell(cell!(r->&format!("{:.7e}", if self.integrated { other * normalizations[bin] } else { other })));
                row.add_cell(cell!(r->&format!("{:.2}%", (other / values[0] - 1.0) * 100.0)));
            }
        }

        table.printstd();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-convolute 
Convolutes a PineAPPL grid with a PDF set

USAGE:
    pineappl convolute [OPTIONS] <INPUT> [--] [PDFSETS]...

ARGS:
    <INPUT>         Path of the input grid
    <PDFSETS>...    LHAPDF id(s) or name of the PDF set(s)

OPTIONS:
    -a, --absolute              Show absolute numbers of the scale variation
    -b, --bins <BINS>...        Selects a subset of bins
    -h, --help                  Print help information
    -i, --integrated            Show integrated numbers (without bin widths) instead of differential
                                ones
    -o, --orders <ORDERS>...    Select orders manually
    -s, --scales <SCALES>       Set the number of scale variations [default: 7] [possible values: 1,
                                3, 7, 9]
";

    const DEFAULT_STR: &str = "bin   etal    disg/detal  scale uncertainty
---+----+----+-----------+--------+--------
  0    2 2.25 3.7527620e2   -3.77%    2.71%
  1 2.25  2.5 3.4521553e2   -3.79%    2.80%
  2  2.5 2.75 3.0001406e2   -3.78%    2.86%
  3 2.75    3 2.4257663e2   -3.77%    2.92%
  4    3 3.25 1.8093343e2   -3.74%    2.95%
  5 3.25  3.5 1.2291115e2   -3.71%    2.98%
  6  3.5    4 5.7851018e1   -3.63%    2.97%
  7    4  4.5 1.3772029e1   -3.46%    2.85%
";

    const DEFAULT_MULTIPLE_PDFS_STR: &str =
        "bin   etal    disg/detal  scale uncertainty NNPDF31_nlo_as_0118_luxqed 
---+----+----+-----------+--------+--------+-------------+-------------
  0    2 2.25 3.7527620e2   -3.77%    2.71%   3.7527620e2         0.00%
  1 2.25  2.5 3.4521553e2   -3.79%    2.80%   3.4521553e2         0.00%
  2  2.5 2.75 3.0001406e2   -3.78%    2.86%   3.0001406e2         0.00%
  3 2.75    3 2.4257663e2   -3.77%    2.92%   2.4257663e2         0.00%
  4    3 3.25 1.8093343e2   -3.74%    2.95%   1.8093343e2         0.00%
  5 3.25  3.5 1.2291115e2   -3.71%    2.98%   1.2291115e2         0.00%
  6  3.5    4 5.7851018e1   -3.63%    2.97%   5.7851018e1         0.00%
  7    4  4.5 1.3772029e1   -3.46%    2.85%   1.3772029e1         0.00%
";

    const ABSOLUTE_STR: &str =
"bin   etal    disg/detal     (1,1)       (2,2)     (0.5,0.5)     (2,1)       (1,2)      (0.5,1)     (1,0.5)  
---+----+----+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------
  0    2 2.25 3.7527620e2 3.7527620e2 3.8169721e2 3.6948620e2 3.7004750e2 3.8546011e2 3.8178087e2 3.6114679e2
  1 2.25  2.5 3.4521553e2 3.4521553e2 3.5114200e2 3.3973093e2 3.4031501e2 3.5487551e2 3.5131193e2 3.3214336e2
  2  2.5 2.75 3.0001406e2 3.0001406e2 3.0511561e2 2.9517901e2 2.9567460e2 3.0860377e2 3.0541248e2 2.8866119e2
  3 2.75    3 2.4257663e2 2.4257663e2 2.4665730e2 2.3861256e2 2.3902145e2 2.4965490e2 2.4699938e2 2.3342681e2
  4    3 3.25 1.8093343e2 1.8093343e2 1.8387616e2 1.7800964e2 1.7821415e2 1.8627534e2 1.8431630e2 1.7416314e2
  5 3.25  3.5 1.2291115e2 1.2291115e2 1.2481060e2 1.2097578e2 1.2099928e2 1.2657016e2 1.2528958e2 1.1835555e2
  6  3.5    4 5.7851018e1 5.7851018e1 5.8647563e1 5.7008512e1 5.6897537e1 5.9567473e1 5.9037178e1 5.5752518e1
  7    4  4.5 1.3772029e1 1.3772029e1 1.3903642e1 1.3622752e1 1.3512675e1 1.4165115e1 1.4094674e1 1.3296051e1
";

    const ABSOLUTE_INTEGRATED_STR: &str =
"bin   etal       integ       (1,1)       (2,2)     (0.5,0.5)     (2,1)       (1,2)      (0.5,1)     (1,0.5)  
---+----+----+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------
  0    2 2.25 9.3819050e1 9.3819050e1 9.5424302e1 9.2371549e1 9.2511874e1 9.6365028e1 9.5445218e1 9.0286698e1
  1 2.25  2.5 8.6303882e1 8.6303882e1 8.7785500e1 8.4932734e1 8.5078752e1 8.8718877e1 8.7827983e1 8.3035839e1
  2  2.5 2.75 7.5003515e1 7.5003515e1 7.6278903e1 7.3794753e1 7.3918649e1 7.7150944e1 7.6353121e1 7.2165298e1
  3 2.75    3 6.0644157e1 6.0644157e1 6.1664326e1 5.9653141e1 5.9755362e1 6.2413725e1 6.1749845e1 5.8356702e1
  4    3 3.25 4.5233358e1 4.5233358e1 4.5969040e1 4.4502410e1 4.4553538e1 4.6568835e1 4.6079075e1 4.3540785e1
  5 3.25  3.5 3.0727788e1 3.0727788e1 3.1202651e1 3.0243945e1 3.0249820e1 3.1642540e1 3.1322395e1 2.9588888e1
  6  3.5    4 2.8925509e1 2.8925509e1 2.9323782e1 2.8504256e1 2.8448768e1 2.9783737e1 2.9518589e1 2.7876259e1
  7    4  4.5 6.8860146e0 6.8860146e0 6.9518210e0 6.8113762e0 6.7563375e0 7.0825577e0 7.0473370e0 6.6480254e0
";

    const BINS_13567_STR: &str = "bin   etal   disg/detal  scale uncertainty
---+----+---+-----------+--------+--------
  1 2.25 2.5 3.4521553e2   -3.79%    2.80%
  3 2.75   3 2.4257663e2   -3.77%    2.92%
  5 3.25 3.5 1.2291115e2   -3.71%    2.98%
  6  3.5   4 5.7851018e1   -3.63%    2.97%
  7    4 4.5 1.3772029e1   -3.46%    2.85%
";

    const INTEGRATED_STR: &str = "bin   etal       integ    scale uncertainty
---+----+----+-----------+--------+--------
  0    2 2.25 9.3819050e1   -3.77%    2.71%
  1 2.25  2.5 8.6303882e1   -3.79%    2.80%
  2  2.5 2.75 7.5003515e1   -3.78%    2.86%
  3 2.75    3 6.0644157e1   -3.77%    2.92%
  4    3 3.25 4.5233358e1   -3.74%    2.95%
  5 3.25  3.5 3.0727788e1   -3.71%    2.98%
  6  3.5    4 2.8925509e1   -3.63%    2.97%
  7    4  4.5 6.8860146e0   -3.46%    2.85%
";

    const INTEGRATED_MULTIPLE_PDFS_STR: &str =
        "bin   etal       integ    scale uncertainty NNPDF31_nlo_as_0118_luxqed 
---+----+----+-----------+--------+--------+-------------+-------------
  0    2 2.25 9.3819050e1   -3.77%    2.71%   9.3819050e1         0.00%
  1 2.25  2.5 8.6303882e1   -3.79%    2.80%   8.6303882e1         0.00%
  2  2.5 2.75 7.5003515e1   -3.78%    2.86%   7.5003515e1         0.00%
  3 2.75    3 6.0644157e1   -3.77%    2.92%   6.0644157e1         0.00%
  4    3 3.25 4.5233358e1   -3.74%    2.95%   4.5233358e1         0.00%
  5 3.25  3.5 3.0727788e1   -3.71%    2.98%   3.0727788e1         0.00%
  6  3.5    4 2.8925509e1   -3.63%    2.97%   2.8925509e1         0.00%
  7    4  4.5 6.8860146e0   -3.46%    2.85%   6.8860146e0         0.00%
";

    const ORDERS_A2_A3_STR: &str = "bin   etal    disg/detal  scale uncertainty
---+----+----+-----------+--------+--------
  0    2 2.25 3.2092052e2   -9.18%    7.92%
  1 2.25  2.5 2.9427151e2   -8.68%    7.41%
  2  2.5 2.75 2.5490261e2   -8.12%    6.84%
  3 2.75    3 2.0561831e2   -7.55%    6.26%
  4    3 3.25 1.5266481e2   -6.97%    5.68%
  5 3.25  3.5 1.0303603e2   -6.38%    5.09%
  6  3.5    4 4.7938981e1   -5.59%    4.31%
  7    4  4.5 1.1075878e1   -4.60%    3.35%
";

    const SCALES_9_STR: &str = "bin   etal    disg/detal  scale uncertainty
---+----+----+-----------+--------+--------
  0    2 2.25 3.7527620e2   -5.55%    3.96%
  1 2.25  2.5 3.4521553e2   -5.55%    4.14%
  2  2.5 2.75 3.0001406e2   -5.53%    4.31%
  3 2.75    3 2.4257663e2   -5.49%    4.46%
  4    3 3.25 1.8093343e2   -5.45%    4.60%
  5 3.25  3.5 1.2291115e2   -5.42%    4.76%
  6  3.5    4 5.7851018e1   -5.37%    4.95%
  7    4  4.5 1.3772029e1   -5.36%    5.22%
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["convolute", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn default() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(DEFAULT_STR);
    }

    #[test]
    fn default_multiple_pdfs() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(DEFAULT_MULTIPLE_PDFS_STR);
    }

    #[test]
    fn absolute() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                "--absolute",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(ABSOLUTE_STR);
    }

    #[test]
    fn absolute_integrated() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                "--absolute",
                "--integrated",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(ABSOLUTE_INTEGRATED_STR);
    }

    #[test]
    fn bins_13567() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                "--bins=1,3,5-7",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(BINS_13567_STR);
    }

    #[test]
    fn integrated() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                "--integrated",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(INTEGRATED_STR);
    }

    #[test]
    fn integrated_multiple_pdfs() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                "--integrated",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(INTEGRATED_MULTIPLE_PDFS_STR);
    }

    #[test]
    fn orders_a2_a3() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                "--orders=a2,a3",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(ORDERS_A2_A3_STR);
    }

    #[test]
    fn scales_9() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                "--scales=9",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(SCALES_9_STR);
    }
}
