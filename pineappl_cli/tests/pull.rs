use assert_cmd::Command;
use std::num::NonZeroUsize;
use std::thread;

const HELP_STR: &str = "Calculates the pull between two different PDF sets

Usage: pineappl pull [OPTIONS] <INPUT> <PDFSET1> <PDFSET2>

Arguments:
  <INPUT>    Path to the input grid
  <PDFSET1>  LHAPDF id or name of the first PDF set
  <PDFSET2>  LHAPDF id or name of the second PDF set

Options:
      --cl <CL>            Confidence level in per cent [default: 68.26894921370858]
  -l, --limit <LIMIT>      The maximum number of channels displayed [default: 10]
  -o, --orders <ORDERS>    Select orders manually
      --threads <THREADS>  Number of threads to utilize [default: {}]
      --digits <DIGITS>    Set the number of digits shown for numerical values [default: 3]
  -h, --help               Print help
";

const DEFAULT_STR: &str = "b   etal    total c pull  c  pull  c  pull  c  pull  c  pull 
     []      [\u{3c3}]     [\u{3c3}]     [\u{3c3}]      [\u{3c3}]      [\u{3c3}]      [\u{3c3}]  
-+----+----+-----+-+-----+-+------+-+------+-+------+-+------
0    2 2.25 3.578 0 3.765 1 -0.108 3 -0.052 4 -0.018 2 -0.009
1 2.25  2.5 3.444 0 3.627 1 -0.095 3 -0.062 4 -0.016 2 -0.010
2  2.5 2.75 3.190 0 3.340 1 -0.074 3 -0.055 4 -0.015 2 -0.006
3 2.75    3 2.795 0 2.883 1 -0.044 3 -0.025 4 -0.011 2 -0.009
4    3 3.25 2.340 0 2.348 3  0.022 1 -0.013 4 -0.010 2 -0.007
5 3.25  3.5 1.874 0 1.808 3  0.084 4 -0.010 2 -0.005 1 -0.002
6  3.5    4 1.466 0 1.390 3  0.176 1 -0.088 2 -0.006 4 -0.006
7    4  4.5 1.224 0 1.435 1 -0.353 3  0.147 4 -0.003 2 -0.002
";

const ORDERS_STR: &str = "b   etal    total c pull  c pull  c pull  c pull  c pull 
     []      [\u{3c3}]     [\u{3c3}]     [\u{3c3}]     [\u{3c3}]     [\u{3c3}]     [\u{3c3}] 
-+----+----+-----+-+-----+-+-----+-+-----+-+-----+-+-----
0    2 2.25 3.631 0 3.631 1 0.000 2 0.000 3 0.000 4 0.000
1 2.25  2.5 3.475 0 3.475 1 0.000 2 0.000 3 0.000 4 0.000
2  2.5 2.75 3.164 0 3.164 1 0.000 2 0.000 3 0.000 4 0.000
3 2.75    3 2.701 0 2.701 1 0.000 2 0.000 3 0.000 4 0.000
4    3 3.25 2.162 0 2.162 1 0.000 2 0.000 3 0.000 4 0.000
5 3.25  3.5 1.632 0 1.632 1 0.000 2 0.000 3 0.000 4 0.000
6  3.5    4 1.240 0 1.240 1 0.000 2 0.000 3 0.000 4 0.000
7    4  4.5 1.202 0 1.202 1 0.000 2 0.000 3 0.000 4 0.000
";

const CL_90_STR: &str = "b   etal    total c pull  c  pull  c  pull  c  pull  c  pull 
     []      [\u{3c3}]     [\u{3c3}]     [\u{3c3}]      [\u{3c3}]      [\u{3c3}]      [\u{3c3}]  
-+----+----+-----+-+-----+-+------+-+------+-+------+-+------
0    2 2.25 2.175 0 2.289 1 -0.065 3 -0.031 4 -0.011 2 -0.006
1 2.25  2.5 2.094 0 2.205 1 -0.058 3 -0.038 4 -0.010 2 -0.006
2  2.5 2.75 1.939 0 2.031 1 -0.045 3 -0.034 4 -0.009 2 -0.004
3 2.75    3 1.699 0 1.753 1 -0.027 3 -0.015 4 -0.007 2 -0.005
4    3 3.25 1.423 0 1.427 3  0.013 1 -0.008 4 -0.006 2 -0.004
5 3.25  3.5 1.140 0 1.099 3  0.051 4 -0.006 2 -0.003 1 -0.001
6  3.5    4 0.891 0 0.845 3  0.107 1 -0.053 2 -0.004 4 -0.003
7    4  4.5 0.744 0 0.872 1 -0.215 3  0.089 4 -0.002 2 -0.001
";

const LIMIT_STR: &str = "b   etal    total c pull 
     []      [\u{3c3}]     [\u{3c3}] 
-+----+----+-----+-+-----
0    2 2.25 3.578 0 3.765
1 2.25  2.5 3.444 0 3.627
2  2.5 2.75 3.190 0 3.340
3 2.75    3 2.795 0 2.883
4    3 3.25 2.340 0 2.348
5 3.25  3.5 1.874 0 1.808
6  3.5    4 1.466 0 1.390
7    4  4.5 1.224 0 1.435
";

const REPLICA0_STR: &str = "b   etal    total c pull  c  pull  c  pull  c  pull  c  pull 
     []      [\u{3c3}]     [\u{3c3}]     [\u{3c3}]      [\u{3c3}]      [\u{3c3}]      [\u{3c3}]  
-+----+----+-----+-+-----+-+------+-+------+-+------+-+------
0    2 2.25 3.583 0 3.770 1 -0.108 3 -0.052 4 -0.018 2 -0.009
1 2.25  2.5 3.445 0 3.628 1 -0.095 3 -0.062 4 -0.016 2 -0.010
2  2.5 2.75 3.188 0 3.338 1 -0.074 3 -0.055 4 -0.015 2 -0.006
3 2.75    3 2.790 0 2.879 1 -0.044 3 -0.025 4 -0.011 2 -0.009
4    3 3.25 2.333 0 2.340 3  0.022 1 -0.013 4 -0.010 2 -0.007
5 3.25  3.5 1.865 0 1.799 3  0.084 4 -0.010 2 -0.005 1 -0.003
6  3.5    4 1.455 0 1.379 3  0.176 1 -0.088 2 -0.006 4 -0.006
7    4  4.5 1.215 0 1.426 1 -0.353 3  0.147 4 -0.003 2 -0.002
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["pull", "--help"])
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
            "pull",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
            "NNPDF40_nnlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}

#[test]
fn orders() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "pull",
            "--orders=a2",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed/0",
            "NNPDF40_nnlo_as_01180/0",
        ])
        .assert()
        .success()
        .stdout(ORDERS_STR);
}

#[test]
fn cl_90() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "pull",
            "--cl=90",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .args([
            "pull",
            "--limit=1",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
            "NNPDF40_nnlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(LIMIT_STR);
}

#[test]
fn replica0() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "pull",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed/0",
            "NNPDF40_nnlo_as_01180/0",
        ])
        .assert()
        .success()
        .stdout(REPLICA0_STR);
}
