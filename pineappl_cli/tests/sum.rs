use assert_cmd::Command;
use assert_fs::NamedTempFile;

const HELP_STR: &str = "Sums two or more bins of a grid together

Usage: pineappl sum <--integrated> <INPUT> <OUTPUT>

Arguments:
  <INPUT>   Path to the input grid
  <OUTPUT>  Path to the modified PineAPPL file

Options:
      --integrated  Sums all bins into a single bin
  -h, --help        Print help
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
