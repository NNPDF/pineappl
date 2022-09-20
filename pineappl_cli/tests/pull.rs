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

const REPLICA0_STR: &str = "b   etal    total l pull  l  pull  l  pull  l  pull  l  pull 
     []      [\u{3c3}]     [\u{3c3}]     [\u{3c3}]      [\u{3c3}]      [\u{3c3}]      [\u{3c3}]  
-+----+----+-----+-+-----+-+------+-+------+-+------+-+------
0    2 2.25 3.582 0 3.766 1 -0.108 3 -0.052 4 -0.016 2 -0.009
1 2.25  2.5 3.452 0 3.631 1 -0.095 3 -0.062 4 -0.016 2 -0.006
2  2.5 2.75 3.193 0 3.336 1 -0.073 3 -0.056 4 -0.010 2 -0.005
3 2.75    3 2.798 0 2.882 1 -0.045 3 -0.024 4 -0.011 2 -0.004
4    3 3.25 2.339 0 2.342 3  0.023 1 -0.013 4 -0.009 2 -0.004
5 3.25  3.5 1.863 0 1.801 3  0.082 2 -0.011 4 -0.007 1 -0.002
6  3.5    4 1.457 0 1.378 3  0.177 1 -0.088 4 -0.007 2 -0.003
7    4  4.5 1.211 0 1.430 1 -0.358 3  0.147 4 -0.006 2 -0.001
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

#[test]
fn replica0() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "pull",
            "--threads=1",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed/0",
            "NNPDF40_nnlo_as_01180/0",
        ])
        .assert()
        .success()
        .stdout(REPLICA0_STR);
}
