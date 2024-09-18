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

const DEFAULT_STR: &str = "b   etal    dsig/detal
     []        [pb]
-+----+----+-----------
0    2 2.25 1.5769498e3
1 2.25  2.5 1.4412311e3
2  2.5 2.75 1.2505273e3
3 2.75    3 1.0077177e3
4    3 3.25 7.4800340e2
5 3.25  3.5 5.0601781e2
6  3.5    4 2.3818927e2
7    4  4.5 5.8009213e1
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
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF40_nnlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}
