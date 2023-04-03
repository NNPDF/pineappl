use assert_cmd::Command;
use assert_fs::NamedTempFile;

const HELP_STR: &str = "Sums two or more bins of a grid together

Usage: pineappl sum <--integrated|--bins <BINS>> <INPUT> <OUTPUT>

Arguments:
  <INPUT>   Path to the input grid
  <OUTPUT>  Path to the modified PineAPPL file

Options:
      --integrated   Sums all bins into a single bin
  -b, --bins <BINS>  Merge specific bins together
  -h, --help         Print help
";

const BINS_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 3.7527620e2    -3.77     2.71
1 2.25  2.5 3.4521553e2    -3.79     2.80
2  2.5 2.75 3.0001406e2    -3.78     2.86
3 2.75    3 2.4257663e2    -3.77     2.92
4    3 3.25 1.8093343e2    -3.74     2.95
5 3.25  3.5 1.2291115e2    -3.71     2.98
6  3.5  4.5 3.5811524e1    -3.59     2.95
";

const INTEGRATED_STR: &str = "b x1     diff     scale uncertainty
  []      []             [%]       
-+-+-+-----------+--------+--------
0 0 1 4.2754327e2    -3.75     2.85
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["sum", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn bins() {
    let output = NamedTempFile::new("bins.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "sum",
            "--bins=6-7",
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
        .stdout(BINS_STR);
}

#[test]
fn integrated() {
    let output = NamedTempFile::new("integrated.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "sum",
            "--integrated",
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
        .stdout(INTEGRATED_STR);
}
