use assert_cmd::Command;
use std::num::NonZeroUsize;
use std::thread;

const HELP_STR: &str = "Calculates PDF uncertainties

Usage: pineappl pdfunc [OPTIONS] <INPUT> <PDFSET>

Arguments:
  <INPUT>   Path to the input grid
  <PDFSET>  LHAPDF id or name of the PDF set

Options:
      --cl <CL>            Confidence level in per cent [default: 68.26894921370858]
  -i, --integrated         Show integrated numbers (without bin widths) instead of differential ones
  -o, --orders <ORDERS>    Select orders manually
      --threads <THREADS>  Number of threads to utilize [default: {}]
      --digits-abs <ABS>   Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>   Set the number of fractional digits shown for relative numbers [default: 2]
  -h, --help               Print help
";

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

const REPLICA0_STR: &str = "b   etal    disg/detal  PDF uncertainty
     []        [pb]           [%]      
-+----+----+-----------+-------+-------
0    2 2.25 3.7527620e2   -1.14    1.14
1 2.25  2.5 3.4521553e2   -1.16    1.16
2  2.5 2.75 3.0001406e2   -1.18    1.18
3 2.75    3 2.4257663e2   -1.22    1.22
4    3 3.25 1.8093343e2   -1.27    1.27
5 3.25  3.5 1.2291115e2   -1.35    1.35
6  3.5    4 5.7851018e1   -1.50    1.50
7    4  4.5 1.3772029e1   -2.76    2.76
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
        .args(["pdfunc", "--help"])
        .assert()
        .success()
        .stdout(
            HELP_STR.replace(
                "{}",
                &thread::available_parallelism()
                    .map_or(1, NonZeroUsize::get)
                    .to_string(),
            ),
        );
}

#[test]
fn default() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
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
        .args([
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
        .args([
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
fn replica0() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--silence-lhapdf",
            "pdfunc",
            "--threads=1",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed/0",
        ])
        .assert()
        .success()
        .stdout(REPLICA0_STR);
}

#[test]
fn orders_a2_as1a2() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
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
