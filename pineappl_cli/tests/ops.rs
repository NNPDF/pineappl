use assert_cmd::Command;
use assert_fs::NamedTempFile;

const HELP_STR: &str = "A collection of various modifying operations on grids

Usage: pineappl ops [OPTIONS] <INPUT> <OUTPUT>

Arguments:
  <INPUT>   Path to the input grid
  <OUTPUT>  Path of the modified PineAPPL file

Options:
      --cc1                                  Charge conjugate the first initial state
      --cc2                                  Charge conjugate the second initial state
      --delete-bins <BIN1-BIN2,...>          Delete bins with the specified indices
  -s, --scale <SCALE>                        Scales all grids with the given factor
      --scale-by-bin <BIN1,BIN2,...>         Scale each bin with a different factor
      --scale-by-order <AS,AL,LR,LF,GLOBAL>  Scales all grids with order-dependent factors
      --upgrade                              Convert the file format to the most recent version
  -h, --help                                 Print help
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

const DELETE_BINS_02_57_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0 2.75    3 2.4257663e2    -3.77     2.92
1    3 3.25 1.8093343e2    -3.74     2.95
";

const DELETE_BINS_25_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 3.7527620e2    -3.77     2.71
1 2.25  2.5 3.4521553e2    -3.79     2.80
2  3.5    4 5.7851018e1    -3.63     2.97
3    4  4.5 1.3772029e1    -3.46     2.85
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

const SCALE_BY_ORDER_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 2.1481594e2    -4.53     3.53
1 2.25  2.5 1.9807977e2    -4.50     3.55
2  2.5 2.75 1.7256275e2    -4.43     3.54
3 2.75    3 1.3976747e2    -4.36     3.52
4    3 3.25 1.0460103e2    -4.26     3.49
5 3.25  3.5 7.1393138e1    -4.17     3.45
6  3.5    4 3.3881527e1    -4.01     3.50
7    4  4.5 8.2340900e0    -3.70     3.92
";
#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["ops", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn cc1() {
    let output = NamedTempFile::new("cc1.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
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
        .args([
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
        .args([
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
        .args([
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
fn delete_bins_02_57() {
    let output = NamedTempFile::new("deleted.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "ops",
            "--delete-bins=0-2,5-7",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--silence-lhapdf",
            "convolute",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DELETE_BINS_02_57_STR);
}

#[test]
fn delete_bins_25() {
    let output = NamedTempFile::new("deleted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "ops",
            "--delete-bins=2-5",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--silence-lhapdf",
            "convolute",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DELETE_BINS_25_STR);
}

#[test]
fn scale_by_bin() {
    let output = NamedTempFile::new("scale_by_bin.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
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
        .args([
            "--silence-lhapdf",
            "convolute",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_BY_BIN_STR);
}

#[test]
fn scale_by_order() {
    let output = NamedTempFile::new("merged.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "ops",
            "--scale-by-order=2,1,0.5,0.5,0.5",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--silence-lhapdf",
            "convolute",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_BY_ORDER_STR);
}

#[test]
fn upgrade() {
    let output = NamedTempFile::new("upgraded.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "ops",
            "--upgrade",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");
}
