use super::helpers::{self, Subcommand};
use anyhow::{bail, Result};
use clap::{Parser, ValueHint};
use prettytable::{cell, Row};
use std::collections::HashSet;
use std::path::PathBuf;

/// Compares the numerical content of two grids with each other.
#[derive(Parser)]
pub struct Opts {
    /// Path to the first grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input1: PathBuf,
    /// Path to the second grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input2: PathBuf,
    /// LHAPDF id or name of the PDF set.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset: String,
    /// Ignore differences in the orders and sum them.
    #[clap(alias = "ignore_orders", long = "ignore-orders")]
    ignore_orders: bool,
    /// Ignore bin limits (but not number of bins).
    #[clap(long = "ignore-bin-limits")]
    ignore_bin_limits: bool,
    /// Ignore differences in the luminosity functions.
    #[clap(long = "ignore-lumis")]
    ignore_lumis: bool,
    /// Select orders of the first grid.
    #[clap(
        long,
        multiple_values = true,
        parse(try_from_str = helpers::parse_order),
        use_value_delimiter = true
    )]
    orders1: Vec<(u32, u32)>,
    /// Select orders of the second grid.
    #[clap(
        long,
        multiple_values = true,
        parse(try_from_str = helpers::parse_order),
        use_value_delimiter = true
    )]
    orders2: Vec<(u32, u32)>,
    /// Scale all results of the first grid.
    #[clap(long, default_value = "1.0")]
    scale1: f64,
    /// Scale all results of the second grid.
    #[clap(long, default_value = "1.0")]
    scale2: f64,
    /// Set the number of fractional digits shown for absolute numbers.
    #[clap(default_value_t = 7, long = "digits-abs", value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[clap(default_value_t = 3, long = "digits-rel", value_name = "REL")]
    digits_rel: usize,
    /// Forces negative PDF values to zero.
    #[clap(long = "force-positive")]
    force_positive: bool,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let grid1 = helpers::read_grid(&self.input1)?;
        let grid2 = helpers::read_grid(&self.input2)?;

        let orders1: HashSet<_> = grid1
            .orders()
            .iter()
            .filter(|order| {
                (order.logxir == 0)
                    && (order.logxif == 0)
                    && (self.orders1.is_empty()
                        || self
                            .orders1
                            .iter()
                            .any(|&o| (order.alphas, order.alpha) == o))
            })
            .collect();
        let orders2: HashSet<_> = grid2
            .orders()
            .iter()
            .filter(|order| {
                (order.logxir == 0)
                    && (order.logxif == 0)
                    && (self.orders2.is_empty()
                        || self
                            .orders2
                            .iter()
                            .any(|&o| (order.alphas, order.alpha) == o))
            })
            .collect();

        let mut diff1: Vec<_> = orders1.difference(&orders2).collect();
        diff1.sort();
        let diff1 = diff1;
        let mut diff2: Vec<_> = orders2.difference(&orders1).collect();
        diff2.sort();
        let diff2 = diff2;

        let mut orders1: Vec<_> = orders1
            .iter()
            .map(|order| (order.alphas, order.alpha))
            .collect();
        orders1.sort_unstable();
        let orders1 = orders1;
        let mut orders2: Vec<_> = orders2
            .iter()
            .map(|order| (order.alphas, order.alpha))
            .collect();
        orders2.sort_unstable();
        let orders2 = orders2;

        if !self.ignore_orders && (!diff1.is_empty() || !diff2.is_empty()) {
            bail!("selected orders differ");
        }

        if !self.ignore_bin_limits && (grid1.bin_info() != grid2.bin_info()) {
            bail!("bins limits differ");
        }

        if grid1.bin_info().bins() != grid2.bin_info().bins() {
            bail!("number of bins differ");
        }

        if (grid1.lumi() != grid2.lumi()) && !self.ignore_lumis {
            bail!("luminosities differ");
        }

        let mut pdf = helpers::create_pdf(&self.pdfset)?;

        let mut table = helpers::create_table();
        let bin_info = grid1.bin_info();

        let left_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.left(i))
            .collect();
        let right_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.right(i))
            .collect();

        let mut title = Row::empty();
        title.add_cell(cell!(c->"b"));
        for i in 0..bin_info.dimensions() {
            let mut cell = cell!(c->format!("x{}", i + 1));
            cell.set_hspan(2);
            title.add_cell(cell);
        }

        if self.ignore_orders {
            let mut cell = cell!(c->"diff");
            cell.set_hspan(3);
            title.add_cell(cell);

            table.set_titles(title);

            let results1 = helpers::convolute(
                &grid1,
                &mut pdf,
                &orders1,
                &[],
                &[],
                1,
                false,
                self.force_positive,
            );
            let results2 = helpers::convolute(
                &grid2,
                &mut pdf,
                &orders2,
                &[],
                &[],
                1,
                false,
                self.force_positive,
            );

            for (bin, (result1, result2)) in results1.iter().zip(results2.iter()).enumerate() {
                let row = table.add_empty_row();

                row.add_cell(cell!(r->bin));
                for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                    row.add_cell(cell!(r->format!("{}", left[bin])));
                    row.add_cell(cell!(r->format!("{}", right[bin])));
                }

                let result1 = result1 * self.scale1;
                let result2 = result2 * self.scale2;

                row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result1)));
                row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result2)));
                row.add_cell(cell!(r->format!("{:.*e}", self.digits_rel,
                if result1 == result2 { 0.0 } else { result1 / result2 - 1.0 })));
            }
        } else {
            let orders = orders1;

            for (alphas, alpha) in &orders {
                let mut cell = cell!(c->format!("O(as^{} a^{})", alphas, alpha));
                cell.set_hspan(3);
                title.add_cell(cell);
            }

            table.set_titles(title);

            let order_results1: Vec<Vec<f64>> = orders
                .iter()
                .map(|&order| {
                    helpers::convolute(
                        &grid1,
                        &mut pdf,
                        &[order],
                        &[],
                        &[],
                        1,
                        false,
                        self.force_positive,
                    )
                })
                .collect();
            let order_results2: Vec<Vec<f64>> = orders
                .iter()
                .map(|&order| {
                    helpers::convolute(
                        &grid2,
                        &mut pdf,
                        &[order],
                        &[],
                        &[],
                        1,
                        false,
                        self.force_positive,
                    )
                })
                .collect();

            for bin in 0..bin_info.bins() {
                let row = table.add_empty_row();

                row.add_cell(cell!(r->bin));
                for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                    row.add_cell(cell!(r->format!("{}", left[bin])));
                    row.add_cell(cell!(r->format!("{}", right[bin])));
                }

                for (result1, result2) in order_results1.iter().zip(order_results2.iter()) {
                    let result1 = result1[bin] * self.scale1;
                    let result2 = result2[bin] * self.scale2;
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result1)));
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_abs, result2)));
                    row.add_cell(cell!(r->format!("{:.*e}", self.digits_rel,
                    if result1 == result2 { 0.0 } else { result1 / result2 - 1.0 })));
                }
            }
        }

        table.printstd();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-diff 
Compares the numerical content of two grids with each other

USAGE:
    pineappl diff [OPTIONS] <INPUT1> <INPUT2> <PDFSET>

ARGS:
    <INPUT1>    Path to the first grid
    <INPUT2>    Path to the second grid
    <PDFSET>    LHAPDF id or name of the PDF set

OPTIONS:
        --digits-abs <ABS>        Set the number of fractional digits shown for absolute numbers
                                  [default: 7]
        --digits-rel <REL>        Set the number of fractional digits shown for relative numbers
                                  [default: 3]
    -h, --help                    Print help information
        --ignore-bin-limits       Ignore bin limits (but not number of bins)
        --ignore-lumis            Ignore differences in the luminosity functions
        --ignore-orders           Ignore differences in the orders and sum them
        --orders1 <ORDERS1>...    Select orders of the first grid
        --orders2 <ORDERS2>...    Select orders of the second grid
        --scale1 <SCALE1>         Scale all results of the first grid [default: 1.0]
        --scale2 <SCALE2>         Scale all results of the second grid [default: 1.0]
";

    const ORDERS1_A2_ORDERS2_A2_STR: &str = "b    x1               O(as^0 a^2)          
-+----+----+-----------+-----------+-------
0    2 2.25 3.2482657e2 3.2482657e2 0.000e0
1 2.25  2.5 2.9755128e2 2.9755128e2 0.000e0
2  2.5 2.75 2.5751142e2 2.5751142e2 0.000e0
3 2.75    3 2.0748091e2 2.0748091e2 0.000e0
4    3 3.25 1.5397599e2 1.5397599e2 0.000e0
5 3.25  3.5 1.0384063e2 1.0384063e2 0.000e0
6  3.5    4 4.8383606e1 4.8383606e1 0.000e0
7    4  4.5 1.1185365e1 1.1185365e1 0.000e0
";

    const ORDERS1_A2_A2AS1_ORDERS2_A2_A2AS1_STR: &str =
        "b    x1               O(as^0 a^2)                     O(as^1 a^2)          
-+----+----+-----------+-----------+-------+-----------+-----------+-------
0    2 2.25 3.2482657e2 3.2482657e2 0.000e0 5.4355679e1 5.4355679e1 0.000e0
1 2.25  2.5 2.9755128e2 2.9755128e2 0.000e0 5.0944018e1 5.0944018e1 0.000e0
2  2.5 2.75 2.5751142e2 2.5751142e2 0.000e0 4.5111446e1 4.5111446e1 0.000e0
3 2.75    3 2.0748091e2 2.0748091e2 0.000e0 3.6958317e1 3.6958317e1 0.000e0
4    3 3.25 1.5397599e2 1.5397599e2 0.000e0 2.8268620e1 2.8268620e1 0.000e0
5 3.25  3.5 1.0384063e2 1.0384063e2 0.000e0 1.9875123e1 1.9875123e1 0.000e0
6  3.5    4 4.8383606e1 4.8383606e1 0.000e0 9.9120372e0 9.9120372e0 0.000e0
7    4  4.5 1.1185365e1 1.1185365e1 0.000e0 2.6961509e0 2.6961509e0 0.000e0
";

    const ORDERS1_A2_A2AS1_IGNORE_ORDERS_STR: &str = "b    x1                   diff              
-+----+----+-----------+-----------+--------
0    2 2.25 3.7918224e2 3.7527620e2 1.041e-2
1 2.25  2.5 3.4849530e2 3.4521553e2 9.501e-3
2  2.5 2.75 3.0262287e2 3.0001406e2 8.696e-3
3 2.75    3 2.4443923e2 2.4257663e2 7.678e-3
4    3 3.25 1.8224461e2 1.8093343e2 7.247e-3
5 3.25  3.5 1.2371575e2 1.2291115e2 6.546e-3
6  3.5    4 5.8295643e1 5.7851018e1 7.686e-3
7    4  4.5 1.3881516e1 1.3772029e1 7.950e-3
";

    const SCALE2_2_STR: &str = "b    x1                O(as^0 a^2)                         O(as^0 a^3)                         O(as^1 a^2)           
-+----+----+-----------+-----------+---------+-------------+-------------+---------+-----------+-----------+---------
0    2 2.25 3.2482657e2 6.4965313e2 -5.000e-1  -3.9060418e0  -7.8120836e0 -5.000e-1 5.4355679e1 1.0871136e2 -5.000e-1
1 2.25  2.5 2.9755128e2 5.9510256e2 -5.000e-1  -3.2797697e0  -6.5595394e0 -5.000e-1 5.0944018e1 1.0188804e2 -5.000e-1
2  2.5 2.75 2.5751142e2 5.1502284e2 -5.000e-1  -2.6088069e0  -5.2176138e0 -5.000e-1 4.5111446e1 9.0222892e1 -5.000e-1
3 2.75    3 2.0748091e2 4.1496182e2 -5.000e-1  -1.8626008e0  -3.7252015e0 -5.000e-1 3.6958317e1 7.3916633e1 -5.000e-1
4    3 3.25 1.5397599e2 3.0795198e2 -5.000e-1  -1.3111794e0  -2.6223588e0 -5.000e-1 2.8268620e1 5.6537240e1 -5.000e-1
5 3.25  3.5 1.0384063e2 2.0768125e2 -5.000e-1 -8.0459807e-1  -1.6091961e0 -5.000e-1 1.9875123e1 3.9750247e1 -5.000e-1
6  3.5    4 4.8383606e1 9.6767212e1 -5.000e-1 -4.4462513e-1 -8.8925027e-1 -5.000e-1 9.9120372e0 1.9824074e1 -5.000e-1
7    4  4.5 1.1185365e1 2.2370731e1 -5.000e-1 -1.0948700e-1 -2.1897400e-1 -5.000e-1 2.6961509e0 5.3923018e0 -5.000e-1
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["diff", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn orders1_a2_orders2_a2() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "diff",
                "--orders1=a2",
                "--orders2=a2",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(ORDERS1_A2_ORDERS2_A2_STR);
    }

    #[test]
    fn orders1_a2_a2as1_orders2_a2_a2as1() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "diff",
                "--orders1=a2,a2as1",
                "--orders2=a2,a2as1",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(ORDERS1_A2_A2AS1_ORDERS2_A2_A2AS1_STR);
    }

    #[test]
    fn orders1_a2_a2as1_ignore_orders() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "diff",
                "--orders1=a2,a2as1",
                "--ignore-orders",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(ORDERS1_A2_A2AS1_IGNORE_ORDERS_STR);
    }

    #[test]
    fn scale2_2() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "diff",
                "--scale2=2",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(SCALE2_2_STR);
    }
}
