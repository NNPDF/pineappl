use assert_cmd::Command;

#[cfg(feature = "applgrid")]
use assert_fs::NamedTempFile;

const HELP_STR: &str = "Converts PineAPPL grids to APPLgrid files

Usage: pineappl export [OPTIONS] <INPUT> <OUTPUT> <PDFSET>

Arguments:
  <INPUT>   Path to the input grid
  <OUTPUT>  Path to the converted grid
  <PDFSET>  LHAPDF id or name of the PDF set to check the converted grid with

Options:
      --accuracy <ACCURACY>          Relative threshold between the table and the converted grid when comparison fails [default: 1e-10]
      --discard-non-matching-scales  Discard non-matching scales that would otherwise lead to panics
  -s, --scales <SCALES>              Set the number of scale variations to compare with if they are available [default: 7] [possible values: 1, 3, 7, 9]
      --silence-libraries            Prevents third-party libraries from printing output
      --digits-abs <ABS>             Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>             Set the number of fractional digits shown for relative numbers [default: 7]
  -h, --help                         Print help information
";

#[cfg(feature = "applgrid")]
const EXPORT_APPLGRID_STR: &str =
    "WARNING: the order O(as^1 a^2 lr^1 lf^0) isn't supported by APPLgrid and will be skipped.
WARNING: the order O(as^1 a^2 lr^0 lf^1) isn't supported by APPLgrid and will be skipped.
WARNING: the order O(as^0 a^3 lr^0 lf^0) isn't supported by APPLgrid and will be skipped.
WARNING: the order O(as^0 a^3 lr^1 lf^0) isn't supported by APPLgrid and will be skipped.
WARNING: the order O(as^0 a^3 lr^0 lf^1) isn't supported by APPLgrid and will be skipped.
b    APPLgrid     PineAPPL     rel. diff
--+------------+------------+--------------
0   7.9566291e0  7.9566291e0 -2.5535130e-15
1   2.3289219e1  2.3289219e1 -1.9984014e-15
2   3.7442697e1  3.7442697e1 -2.1094237e-15
3   5.0087316e1  5.0087316e1 -3.1086245e-15
4   6.0873237e1  6.0873237e1 -2.8865799e-15
5   6.8944378e1  6.8944378e1 -3.8857806e-15
6   7.4277783e1  7.4277783e1 -2.8865799e-15
7   7.6356931e1  7.6356931e1 -3.7747583e-15
8   7.5009607e1  7.5009607e1 -1.5543122e-15
9   7.0045787e1  7.0045787e1 -1.2212453e-15
10  6.0009803e1  6.0009803e1 -7.7715612e-16
11  4.6770515e1  4.6770515e1  4.4408921e-16
12  3.3569217e1  3.3569217e1  1.5543122e-15
13  2.1820341e1  2.1820341e1  1.3322676e-15
14  1.2542026e1  1.2542026e1  2.2204460e-16
15  6.0879666e0  6.0879666e0 -1.3322676e-15
16  1.5789361e0  1.5789361e0 -1.5543122e-15
17 7.4959880e-2 7.4959880e-2 -1.3322676e-15
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["export", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
#[ignore]
#[cfg(feature = "applgrid")]
fn export_applgrid() {
    let output = NamedTempFile::new("converted.appl").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "export",
            "--silence-libraries",
            "../test-data/LHCB_DY_8TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(EXPORT_APPLGRID_STR);
}
