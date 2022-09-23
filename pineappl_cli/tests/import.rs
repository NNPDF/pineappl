use assert_cmd::Command;

#[cfg(any(feature = "applgrid", feature = "fastnlo", feature = "fktable"))]
use assert_fs::NamedTempFile;

const HELP_STR: &str = "pineappl-import 
Converts APPLgrid/fastNLO/FastKernel files to PineAPPL grids

USAGE:
    pineappl import [OPTIONS] <INPUT> <OUTPUT> <PDFSET>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path to the converted grid
    <PDFSET>    LHAPDF id or name of the PDF set to check the converted grid with

OPTIONS:
        --accuracy <ACCURACY>    Relative threshold between the table and the converted grid when
                                 comparison fails [default: 1e-10]
        --alpha <ALPHA>          LO coupling power in alpha [default: 0]
        --digits-abs <ABS>       Set the number of fractional digits shown for absolute numbers
                                 [default: 7]
        --digits-rel <REL>       Set the number of fractional digits shown for relative numbers
                                 [default: 7]
        --dis-pid <DIS_PID>      Particle ID for the non-hadronic initial states if it cannot be
                                 determined from the grid [default: 11]
    -h, --help                   Print help information
        --no-optimize            Do not optimize converted grid
        --silence-libraries      Prevents third-party libraries from printing output
";

#[cfg(feature = "fastnlo")]
const IMPORT_FIX_GRID_STR: &str = "b   PineAPPL     fastNLO      rel. diff
-+------------+------------+--------------
0 2.9158424e-4 2.9158424e-4 -2.9976022e-15
1 2.4657895e-4 2.4657895e-4 -2.8865799e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_STR: &str = "b   PineAPPL     fastNLO      rel. diff
-+------------+------------+--------------
0  8.2754182e1  8.2754182e1 -1.3544721e-14
1  3.6097335e1  3.6097335e1 -6.8833828e-15
2  8.0048746e0  8.0048746e0  5.3290705e-15
3 9.4319392e-1 9.4319392e-1  5.5511151e-15
";

#[cfg(feature = "fktable")]
const IMPORT_DIS_FKTABLE_STR: &str = "b   x1       diff      scale uncertainty
    []        []              [%]       
--+--+--+-------------+--------+--------
 0  0  1   1.8900584e0   -69.81   107.21
 1  1  2   1.4830883e0   -69.63    98.20
 2  2  3   1.1497012e0   -69.68    90.39
 3  3  4  8.7974506e-1   -69.38    82.45
 4  4  5  7.0882550e-1   -69.14    70.54
 5  5  6  5.7345845e-1   -69.02    59.25
 6  6  7  4.7744833e-1   -68.37    44.68
 7  7  8  4.1037984e-1   -67.36    29.06
 8  8  9  4.0362470e-1   -65.97    12.72
 9  9 10  4.2613006e-1   -64.37     0.00
10 10 11  3.7669466e-1   -63.54     0.00
11 11 12  2.9572989e-1   -62.91     0.00
12 12 13  2.0869778e-1   -62.28     0.00
13 13 14  1.2602675e-1   -61.64     0.00
14 14 15  6.4220769e-2   -60.94     0.00
15 15 16  2.5434367e-2   -60.76     0.00
16 16 17  7.6070428e-3   -59.97     0.00
17 17 18  2.1848546e-3   -60.65     0.00
18 18 19  6.2309735e-4   -57.15     0.00
19 19 20 -1.0496472e-4     0.00   -55.42
";

#[cfg(feature = "fktable")]
const IMPORT_HADRONIC_FKTABLE_STR: &str = "b x1     diff     scale uncertainty
  []      []             [%]       
-+-+-+-----------+--------+--------
0 0 1 7.7624461e2   -86.97     0.00
";

#[cfg(feature = "applgrid")]
const IMPORT_PHOTON_GRID_STR: &str = "b   PineAPPL     APPLgrid     rel. diff
-+------------+------------+--------------
0 5.5621307e-4 5.5621307e-4 -1.5543122e-15
";

#[cfg(feature = "applgrid")]
const IMPORT_APPLGRID_STR: &str = "b  PineAPPL    APPLgrid     rel. diff
-+-----------+-----------+--------------
0 2.9884537e6 2.9884537e6 -6.6613381e-16
";

#[cfg(feature = "applgrid")]
const IMPORT_NEW_APPLGRID_STR: &str = "b   PineAPPL    APPLgrid     rel. diff
--+-----------+-----------+--------------
0  6.2634897e2 6.2634897e2  1.5543122e-15
1  6.2847078e2 6.2847078e2  4.4408921e-16
2  6.3163323e2 6.3163323e2  6.6613381e-16
3  6.3586556e2 6.3586556e2  4.4408921e-16
4  6.4139163e2 6.4139163e2  1.3322676e-15
5  6.4848088e2 6.4848088e2  2.2204460e-16
6  6.5354150e2 6.5354150e2 -4.2188475e-15
7  6.5377566e2 6.5377566e2  6.6613381e-16
8  6.5094729e2 6.5094729e2  8.8817842e-16
9  6.3588760e2 6.3588760e2 -3.8857806e-15
10 5.9810718e2 5.9810718e2  2.4424907e-15
";

const IMPORT_FILE_FORMAT_FAILURE_STR: &str = "Error: could not detect file format
";

#[cfg(feature = "fastnlo")]
const IMPORT_GRID_COMPARISON_FAILURE_STR: &str = "Error: grids are different
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["import", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_fix_grid() {
    let output = NamedTempFile::new("converted1.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "import",
            "--silence-libraries",
            "data/NJetEvents_0-0-2.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(IMPORT_FIX_GRID_STR);
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "import",
            "--silence-libraries",
            "data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(IMPORT_FLEX_GRID_STR);
}

#[test]
#[cfg(feature = "fktable")]
fn import_dis_fktable() {
    let output = NamedTempFile::new("converted3.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "import",
            "data/FK_POSXDQ.dat",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout("file was converted, but we cannot check the conversion for this type\n");

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
        .stdout(IMPORT_DIS_FKTABLE_STR);
}

#[test]
#[cfg(feature = "fktable")]
fn import_hadronic_fktable() {
    let output = NamedTempFile::new("converted4.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "import",
            "data/FK_ATLASTTBARTOT13TEV.dat",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout("file was converted, but we cannot check the conversion for this type\n");

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
        .stdout(IMPORT_HADRONIC_FKTABLE_STR);
}

#[test]
#[cfg(feature = "applgrid")]
fn import_photon_grid() {
    let output = NamedTempFile::new("converted5.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "import",
            "--silence-libraries",
            "data/LHCBWZMU7TEV_PI_part1.root",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(IMPORT_PHOTON_GRID_STR);
}

#[test]
#[cfg(feature = "applgrid")]
fn import_applgrid() {
    let output = NamedTempFile::new("converted6.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "import",
            "--silence-libraries",
            "data/ATLASWPT11-Wplus_tot.root",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(IMPORT_APPLGRID_STR);
}

#[test]
#[cfg(feature = "applgrid")]
fn import_new_applgrid() {
    let output = NamedTempFile::new("converted7.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "import",
            "--silence-libraries",
            "data/atlas-atlas-wpm-arxiv-1109.5141-xsec001.appl",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(IMPORT_NEW_APPLGRID_STR);
}

#[test]
fn import_file_format_failure() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "import",
            "file-doesnt.exist",
            "no.output",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(IMPORT_FILE_FORMAT_FAILURE_STR);
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_grid_comparison_failure() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "import",
            "--accuracy=0",
            "--silence-libraries",
            "data/NJetEvents_0-0-2.tab.gz",
            "this-file-wont-be-writen",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(IMPORT_GRID_COMPARISON_FAILURE_STR)
        .stdout(IMPORT_FIX_GRID_STR);
}
