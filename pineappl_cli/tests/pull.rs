use assert_cmd::Command;
use std::num::NonZeroUsize;
use std::thread;

const HELP_STR: &str = "Calculates the pull between two different PDF sets

Usage: pineappl pull [OPTIONS] <INPUT> <CONV_FUNS1> <CONV_FUNS2>

Arguments:
  <INPUT>       Path to the input grid
  <CONV_FUNS1>  LHAPDF ID(s) or name(s) of the first PDF(s)/FF(s)
  <CONV_FUNS2>  LHAPDF ID(s) or name(s) of the second PDF(s)/FF(s)

Options:
      --pull-from <IDX>    Index of the convolution functions for which the pull should be calculated [default: 0]
      --cl <CL>            Confidence level in per cent [default: 68.26894921370858]
  -l, --limit <LIMIT>      The maximum number of channels displayed [default: 10]
  -o, --orders <ORDERS>    Select orders manually
      --threads <THREADS>  Number of threads to utilize [default: {}]
      --digits <DIGITS>    Set the number of digits shown for numerical values [default: 3]
  -h, --help               Print help
";

// last two columns are zero because the PDFs don't have a photon
const DEFAULT_STR: &str = "b   etal    total c pull  c pull  c  pull  c pull  c pull 
     []      [\u{3c3}]     [\u{3c3}]     [\u{3c3}]     [\u{3c3}]      [\u{3c3}]     [\u{3c3}] 
-+----+----+-----+-+-----+-+-----+-+------+-+-----+-+-----
0    2 2.25 0.636 0 0.525 3 0.062 1  0.049 2 0.000 4 0.000
1 2.25  2.5 0.695 0 0.571 3 0.069 1  0.054 2 0.000 4 0.000
2  2.5 2.75 0.709 0 0.567 3 0.080 1  0.062 2 0.000 4 0.000
3 2.75    3 0.690 0 0.527 3 0.090 1  0.072 2 0.000 4 0.000
4    3 3.25 0.636 0 0.450 3 0.100 1  0.086 2 0.000 4 0.000
5 3.25  3.5 0.556 0 0.353 3 0.102 1  0.101 2 0.000 4 0.000
6  3.5    4 0.434 0 0.242 1 0.118 3  0.073 2 0.000 4 0.000
7    4  4.5 0.370 0 0.304 1 0.094 3 -0.028 2 0.000 4 0.000
";

// last two columns are zero because the PDFs don't have a photon
const CL_90_STR: &str = "b   etal    total c pull  c  pull  c  pull  c pull  c pull 
     []      [\u{3c3}]     [\u{3c3}]     [\u{3c3}]      [\u{3c3}]      [\u{3c3}]     [\u{3c3}] 
-+----+----+-----+-+-----+-+------+-+------+-+-----+-+-----
0    2 2.25 0.203 0 0.191 1  0.036 3 -0.024 2 0.000 4 0.000
1 2.25  2.5 0.222 0 0.211 1  0.038 3 -0.027 2 0.000 4 0.000
2  2.5 2.75 0.222 0 0.210 1  0.043 3 -0.031 2 0.000 4 0.000
3 2.75    3 0.212 0 0.196 1  0.050 3 -0.034 2 0.000 4 0.000
4    3 3.25 0.199 0 0.176 1  0.058 3 -0.035 2 0.000 4 0.000
5 3.25  3.5 0.198 0 0.173 1  0.062 3 -0.037 2 0.000 4 0.000
6  3.5    4 0.277 0 0.270 1  0.044 3 -0.038 2 0.000 4 0.000
7    4  4.5 0.603 0 0.698 1 -0.063 3 -0.032 2 0.000 4 0.000
";

const LIMIT_STR: &str = "b   etal    total  c  pull 
     []      [\u{3c3}]      [\u{3c3}]  
-+----+----+------+-+------
0    2 2.25 -0.504 0 -0.504
1 2.25  2.5 -0.535 0 -0.535
2  2.5 2.75 -0.526 0 -0.526
3 2.75    3 -0.483 0 -0.483
4    3 3.25 -0.406 0 -0.406
5 3.25  3.5 -0.312 0 -0.312
6  3.5    4 -0.211 0 -0.211
7    4  4.5 -0.270 0 -0.270
";

const REPLICA0_STR: &str = "b   etal      total   c   pull    c    pull    c    pull    c    pull    c    pull   
     []        [\u{3c3}]         [\u{3c3}]         [\u{3c3}]          [\u{3c3}]          [\u{3c3}]          [\u{3c3}]    
-+----+----+---------+-+---------+-+----------+-+----------+-+----------+-+----------
0    2 2.25 0.8254315 0 0.8205578 1  0.0381204 3 -0.0228547 4 -0.0068566 2 -0.0035354
1 2.25  2.5 0.9013633 0 0.8981709 1  0.0412455 3 -0.0280520 4 -0.0063410 2 -0.0036601
2  2.5 2.75 0.9374630 0 0.9268791 1  0.0496495 3 -0.0308432 4 -0.0058393 2 -0.0023832
3 2.75    3 0.9270558 0 0.9018028 1  0.0625365 3 -0.0289913 4 -0.0046141 2 -0.0036782
4    3 3.25 0.8818159 0 0.8334133 1  0.0791614 3 -0.0236047 4 -0.0042131 2 -0.0029411
5 3.25  3.5 0.8095801 0 0.7409452 1  0.0917985 3 -0.0160410 4 -0.0048014 2 -0.0023212
6  3.5    4 0.6695089 0 0.5991452 1  0.0806321 3 -0.0038535 2 -0.0032678 4 -0.0031471
7    4  4.5 0.1659039 0 0.1773824 3 -0.0075349 4 -0.0025330 2 -0.0014339 1  0.0000233
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
            "ABMP16als118_5_nnlo",
            "CT18NNLO",
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
            "pull",
            "--cl=90",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "MSHT20nnlo_as118",
            "CT18NNLO",
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
            // since we only show one column, the higher orders aren't that relevant
            "--orders=a2",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "CT18NNLO",
            "ABMP16als118_5_nnlo",
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
            "--digits=7",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed/0",
            "CT18NNLO",
        ])
        .assert()
        .success()
        .stdout(REPLICA0_STR);
}
