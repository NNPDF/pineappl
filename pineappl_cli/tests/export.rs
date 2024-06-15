use assert_cmd::Command;

#[cfg(feature = "applgrid")]
use assert_fs::NamedTempFile;

const HELP_STR: &str = "Converts PineAPPL grids to APPLgrid files

Usage: pineappl export [OPTIONS] <INPUT> <OUTPUT> <CONV_FUNS>...

Arguments:
  <INPUT>         Path to the input grid
  <OUTPUT>        Path to the converted grid
  <CONV_FUNS>...  LHAPDF ID(s) or name of the PDF(s)/FF(s) to check the converted grid with

Options:
      --accuracy <ACCURACY>          Relative threshold between the table and the converted grid when comparison fails [default: 1e-10]
      --discard-non-matching-scales  Discard non-matching scales that would otherwise lead to panics
  -s, --scales <SCALES>              Set the number of scale variations to compare with if they are available [default: 7] [possible values: 1, 3, 7, 9]
      --digits-abs <ABS>             Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>             Set the number of fractional digits shown for relative numbers [default: 7]
  -h, --help                         Print help
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

#[cfg(feature = "applgrid")]
const EXPORT_DIS_APPLGRID_STR: &str =
    "WARNING: the order O(as^1 a^0 lr^0 lf^1) isn't supported by APPLgrid and will be skipped.
b   APPLgrid    PineAPPL     rel. diff
--+-----------+-----------+--------------
0  2.8829972e0 2.8829972e0  3.3306691e-15
1  9.4531032e0 9.4531032e0 -3.3306691e-15
2  1.3796888e1 1.3796888e1 -1.5543122e-15
3  1.1500842e1 1.1500842e1 -7.7715612e-16
4  1.0201793e1 1.0201793e1    0.0000000e0
5  1.8236948e1 1.8236948e1 -2.2204460e-16
6  1.7932872e1 1.7932872e1  4.4408921e-16
7  1.1772337e1 1.1772337e1  1.1102230e-15
8  1.8758110e1 1.8758110e1  6.6613381e-16
9  1.1491805e1 1.1491805e1 -7.7715612e-16
10 2.2083762e1 2.2083762e1  8.8817842e-16
11 1.3902297e1 1.3902297e1  6.6613381e-16
12 1.2194279e1 1.2194279e1  6.6613381e-16
13 2.2536654e1 2.2536654e1  1.3322676e-15
14 1.5937161e1 1.5937161e1  1.3322676e-15
15 9.4835177e0 9.4835177e0 -3.3306691e-16
16 1.2293034e1 1.2293034e1    0.0000000e0
17 1.5390536e1 1.5390536e1  1.7763568e-15
18 1.3361386e1 1.3361386e1  2.2204460e-16
19 1.2423806e1 1.2423806e1  4.4408921e-16
20 1.5550940e1 1.5550940e1  2.4424907e-15
21 1.7110366e1 1.7110366e1  1.7763568e-15
22 1.4208475e1 1.4208475e1  1.1102230e-15
23 1.0023924e1 1.0023924e1 -8.8817842e-16
24 4.0031117e0 4.0031117e0 -4.4408921e-16
25 1.0531971e1 1.0531971e1 -1.2212453e-15
26 1.3879598e1 1.3879598e1  2.2204460e-16
27 1.7044239e1 1.7044239e1  1.7763568e-15
28 1.0076863e1 1.0076863e1 -1.4432899e-15
29 1.1685934e1 1.1685934e1 -5.5511151e-16
30 1.4765962e1 1.4765962e1 -2.2204460e-16
31 1.5383672e1 1.5383672e1 -6.6613381e-16
32 4.5167719e0 4.5167719e0    0.0000000e0
33 1.1469329e1 1.1469329e1 -1.4432899e-15
34 4.6285601e0 4.6285601e0    0.0000000e0
35 1.5054457e1 1.5054457e1  1.1102230e-15
36 4.3142208e0 4.3142208e0 -2.2204460e-16
37 1.1179505e1 1.1179505e1 -1.2212453e-15
38 1.1405936e1 1.1405936e1 -6.6613381e-16
39 4.5331796e0 4.5331796e0  4.4408921e-16
40 4.3137869e0 4.3137869e0    0.0000000e0
41 1.1074106e1 1.1074106e1 -1.1102230e-15
42 4.1007590e0 4.1007590e0    0.0000000e0
43 4.2016757e0 4.2016757e0    0.0000000e0
44 3.9803998e0 3.9803998e0  2.2204460e-16
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["export", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
#[cfg(feature = "applgrid")]
fn export_applgrid() {
    let output = NamedTempFile::new("converted.appl").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "export",
            "../test-data/LHCB_DY_8TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(EXPORT_APPLGRID_STR));
}

#[test]
#[cfg(feature = "applgrid")]
fn export_dis_applgrid() {
    let output1 = NamedTempFile::new("remapped.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--remap=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45",
            "../test-data/NUTEV_CC_NU_FE_SIGMARED.pineappl.lz4",
            output1.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    let output2 = NamedTempFile::new("converted-dis.appl").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "export",
            output1.path().to_str().unwrap(),
            output2.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(EXPORT_DIS_APPLGRID_STR));
}
