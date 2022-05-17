use super::helpers::{self, Subcommand};
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
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
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
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let grid = helpers::read_grid(&self.input)?;
        let set = helpers::create_pdfset(&self.pdfset)?;
        let pdfs = set.mk_pdfs();

        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

        let results: Vec<f64> = pdfs
            .into_par_iter()
            .flat_map(|mut pdf| {
                helpers::convolute(&grid, &mut pdf, &self.orders, &[], &[], 1, self.integrated)
            })
            .collect();

        let bin_info = grid.bin_info();
        let left_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.left(i))
            .collect();
        let right_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.right(i))
            .collect();

        let (x, y_label, y_unit) = helpers::labels_and_units(&grid, self.integrated);
        let mut title = Row::empty();
        title.add_cell(cell!(c->"b"));
        for (x_label, x_unit) in x {
            let mut cell = cell!(c->format!("{}\n[{}]", x_label, x_unit));
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        title.add_cell(cell!(c->format!("{}\n[{}]", y_label, y_unit)));
        title.add_cell(cell!(c->"PDF uncertainty\n[%]").with_hspan(2));

        let mut table = helpers::create_table();
        table.set_titles(title);

        for bin in 0..bin_info.bins() {
            let values: Vec<_> = results
                .iter()
                .skip(bin)
                .step_by(bin_info.bins())
                .copied()
                .collect();
            let uncertainty = set.uncertainty(&values, self.cl, false)?;

            let row = table.add_empty_row();

            row.add_cell(cell!(r->format!("{}", bin)));
            for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                row.add_cell(cell!(r->format!("{}", left[bin])));
                row.add_cell(cell!(r->format!("{}", right[bin])));
            }
            row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, uncertainty.central)));
            row.add_cell(
                cell!(r->format!("{:.*}", self.digits_rel, (-uncertainty.errminus / uncertainty.central) * 100.0)),
            );
            row.add_cell(
                cell!(r->format!("{:.*}", self.digits_rel, (uncertainty.errplus / uncertainty.central) * 100.0)),
            );
        }

        table.printstd();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-pdfunc 
Calculates PDF uncertainties

USAGE:
    pineappl pdfunc [OPTIONS] <INPUT> <PDFSET>

ARGS:
    <INPUT>     Path to the input grid
    <PDFSET>    LHAPDF id or name of the PDF set

OPTIONS:
        --cl <CL>               Confidence level in per cent [default: 68.26894921370858]
        --digits-abs <ABS>      Set the number of fractional digits shown for absolute numbers
                                [default: 7]
        --digits-rel <REL>      Set the number of fractional digits shown for relative numbers
                                [default: 2]
    -h, --help                  Print help information
    -i, --integrated            Show integrated numbers (without bin widths) instead of differential
                                ones
    -o, --orders <ORDERS>...    Select orders manually
        --threads <THREADS>     Number of threads to utilize";

    const DEFAULT_STR: &str = "b   etal    disg/detal  PDF uncertainty
     []        [pb]           [%]      
-+----+----+-----------+-------+-------
0    2 2.25 3.7528868e2   -1.14    1.14
1 2.25  2.5 3.4521365e2   -1.16    1.16
2  2.5 2.75 3.0000102e2   -1.18    1.18
3 2.75    3 2.4255656e2   -1.22    1.22
4    3 3.25 1.8091118e2   -1.27    1.27
5 3.25  3.5 1.2289094e2   -1.35    1.35
6  3.5    4 5.7837137e1   -1.50    1.50
7    4  4.5 1.3765722e1   -2.76    2.76
";

    const CL_90_STR: &str = "b   etal    disg/detal  PDF uncertainty
     []        [pb]           [%]      
-+----+----+-----------+-------+-------
0    2 2.25 3.7528868e2   -1.88    1.88
1 2.25  2.5 3.4521365e2   -1.90    1.90
2  2.5 2.75 3.0000102e2   -1.95    1.95
3 2.75    3 2.4255656e2   -2.00    2.00
4    3 3.25 1.8091118e2   -2.08    2.08
5 3.25  3.5 1.2289094e2   -2.22    2.22
6  3.5    4 5.7837137e1   -2.48    2.48
7    4  4.5 1.3765722e1   -4.54    4.54
";

    const INTEGRATED_STR: &str = "b   etal       integ    PDF uncertainty
     []         []            [%]      
-+----+----+-----------+-------+-------
0    2 2.25 9.3822169e1   -1.14    1.14
1 2.25  2.5 8.6303411e1   -1.16    1.16
2  2.5 2.75 7.5000256e1   -1.18    1.18
3 2.75    3 6.0639140e1   -1.22    1.22
4    3 3.25 4.5227794e1   -1.27    1.27
5 3.25  3.5 3.0722735e1   -1.35    1.35
6  3.5    4 2.8918568e1   -1.50    1.50
7    4  4.5 6.8828610e0   -2.76    2.76
";

    const ORDERS_A2_AS1A2_STR: &str = "b   etal    disg/detal  PDF uncertainty
     []        [pb]           [%]      
-+----+----+-----------+-------+-------
0    2 2.25 3.7919477e2   -1.14    1.14
1 2.25  2.5 3.4849336e2   -1.16    1.16
2  2.5 2.75 3.0260975e2   -1.18    1.18
3 2.75    3 2.4441905e2   -1.22    1.22
4    3 3.25 1.8222226e2   -1.26    1.26
5 3.25  3.5 1.2369548e2   -1.35    1.35
6  3.5    4 5.8281739e1   -1.50    1.50
7    4  4.5 1.3875186e1   -2.77    2.77
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["pdfunc", "--help"])
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
                "pdfunc",
                "--threads=1",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
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
                "pdfunc",
                "--cl=90",
                "--threads=1",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(CL_90_STR);
    }

    #[test]
    fn integrated() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "pdfunc",
                "--integrated",
                "--threads=1",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(INTEGRATED_STR);
    }

    #[test]
    fn orders_a2_as1a2() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "pdfunc",
                "--orders=a2,as1a2",
                "--threads=1",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(ORDERS_A2_AS1A2_STR);
    }
}
