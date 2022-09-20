use assert_cmd::Command;
use assert_fs::NamedTempFile;

const HELP_STR: &str = "pineappl-ops 
A collection of various modifying operations on grids

USAGE:
    pineappl ops [OPTIONS] <INPUT> <OUTPUT>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path of the modified PineAPPL file

OPTIONS:
        --cc1                               Charge conjugate the first initial state
        --cc2                               Charge conjugate the second initial state
    -h, --help                              Print help information
        --scale-by-bin <SCALE_BY_BIN>...    Scale each bin with a different factor
";

const DEFAULT_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 3.7527620e2    -3.77     2.71
1 2.25  2.5 3.4521553e2    -3.79     2.80
2  2.5 2.75 3.0001406e2    -3.78     2.86
3 2.75    3 2.4257663e2    -3.77     2.92
4    3 3.25 1.8093343e2    -3.74     2.95
5 3.25  3.5 1.2291115e2    -3.71     2.98
6  3.5    4 5.7851018e1    -3.63     2.97
7    4  4.5 1.3772029e1    -3.46     2.85
";

const SCALE_BY_BIN_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 3.7527620e2    -3.77     2.71
1 2.25  2.5 6.9043106e2    -3.79     2.80
2  2.5 2.75 9.0004217e2    -3.78     2.86
3 2.75    3 9.7030651e2    -3.77     2.92
4    3 3.25 9.0466716e2    -3.74     2.95
5 3.25  3.5 7.3746691e2    -3.71     2.98
6  3.5    4 4.0495712e2    -3.63     2.97
7    4  4.5 1.1017623e2    -3.46     2.85
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["ops", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn cc1() {
    let output = NamedTempFile::new("cc1.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "ops",
            "--cc1",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}

#[test]
fn cc2() {
    let output = NamedTempFile::new("cc2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "ops",
            "--cc2",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}

#[test]
fn scale_by_bin() {
    let output = NamedTempFile::new("scale_by_bin.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "ops",
            "--scale-by-bin=1,2,3,4,5,6,7,8",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_BY_BIN_STR);
}
