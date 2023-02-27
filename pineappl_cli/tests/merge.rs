use assert_cmd::Command;
use assert_fs::NamedTempFile;

const HELP_STR: &str = "Merges one or more PineAPPL grids together

Usage: pineappl merge <OUTPUT> <INPUT>...

Arguments:
  <OUTPUT>    Path of the merged PineAPPL file
  <INPUT>...  Path(s) of the files that should be merged

Options:
  -h, --help  Print help
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

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["merge", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn default() {
    let output = NamedTempFile::new("merged.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "merge",
            output.path().to_str().unwrap(),
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
