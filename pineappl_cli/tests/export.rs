#![allow(missing_docs)]

use assert_cmd::Command;

#[cfg(feature = "applgrid")]
use assert_fs::NamedTempFile;

const HELP_STR: &str = "Converts PineAPPL grids to APPLgrid files

Usage: pineappl export [OPTIONS] <INPUT> <OUTPUT> <CONV_FUNS>

Arguments:
  <INPUT>      Path to the input grid
  <OUTPUT>     Path to the converted grid
  <CONV_FUNS>  LHAPDF ID(s) or name of the PDF(s)/FF(s) to check the converted grid with

Options:
      --accuracy <ACCURACY>          Relative threshold between the table and the converted grid when comparison fails [default: 1e-10]
      --discard-non-matching-values  Discard non-matching scales and momentum fractions that would otherwise fail the export
  -s, --scales <SCALES>              Set the number of scale variations to compare with if they are available [default: 7] [possible values: 1, 3, 7, 9]
      --digits-abs <ABS>             Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>             Set the number of fractional digits shown for relative numbers [default: 7]
  -h, --help                         Print help
";

#[cfg(feature = "applgrid")]
const EXPORT_APPLGRID_STR: &str =
    "WARNING: the order O(as^1 a^2 lr^1 lf^0 la^0) isn't supported by APPLgrid and will be skipped.
WARNING: the order O(as^1 a^2 lr^0 lf^1 la^0) isn't supported by APPLgrid and will be skipped.
WARNING: the order O(as^0 a^3 lr^0 lf^0 la^0) isn't supported by APPLgrid and will be skipped.
WARNING: the order O(as^0 a^3 lr^1 lf^0 la^0) isn't supported by APPLgrid and will be skipped.
WARNING: the order O(as^0 a^3 lr^0 lf^1 la^0) isn't supported by APPLgrid and will be skipped.
b    APPLgrid     PineAPPL     rel. diff
--+------------+------------+--------------
0   7.9566291e0  7.9566291e0 -2.5535130e-15
1   2.3289219e1  2.3289219e1 -1.7763568e-15
2   3.7442697e1  3.7442697e1 -2.1094237e-15
3   5.0087316e1  5.0087316e1 -2.9976022e-15
4   6.0873237e1  6.0873237e1 -2.5535130e-15
5   6.8944378e1  6.8944378e1 -3.5527137e-15
6   7.4277783e1  7.4277783e1 -2.8865799e-15
7   7.6356931e1  7.6356931e1 -3.3306691e-15
8   7.5009607e1  7.5009607e1 -1.8873791e-15
9   7.0045787e1  7.0045787e1 -9.9920072e-16
10  6.0009803e1  6.0009803e1 -7.7715612e-16
11  4.6770515e1  4.6770515e1  4.4408921e-16
12  3.3569217e1  3.3569217e1  1.5543122e-15
13  2.1820341e1  2.1820341e1  1.1102230e-15
14  1.2542026e1  1.2542026e1  2.2204460e-16
15  6.0879666e0  6.0879666e0 -1.3322676e-15
16  1.5789361e0  1.5789361e0 -1.5543122e-15
17 7.4959880e-2 7.4959880e-2 -1.1102230e-15
";

#[cfg(feature = "applgrid")]
const EXPORT_DIS_APPLGRID_STR: &str =
    "WARNING: the order O(as^1 a^0 lr^0 lf^1 la^0) isn't supported by APPLgrid and will be skipped.
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

#[cfg(feature = "applgrid")]
const EXPORT_NNLO_AK5_PTJ_STR: &str = "b    APPLgrid     PineAPPL     rel. diff
--+------------+------------+--------------
0   7.1475199e8  7.1475253e8   7.6053554e-7
1   2.5536645e8  2.5536667e8   8.6645764e-7
2   8.0742409e7  8.0742409e7  1.2478907e-13
3   4.0260433e7  4.0260445e7   2.8169024e-7
4   1.9839883e7  1.9839887e7   2.1808465e-7
5   1.0119027e7  1.0119028e7   1.3544486e-7
6   5.5702376e6  5.5702379e6   5.8673922e-8
7   3.1322190e6  3.1322207e6   5.3759875e-7
8   1.9447742e6  1.9447743e6   4.2597541e-8
9   1.1701564e6  1.1701565e6   7.4993014e-8
10  7.2929160e5  7.2929198e5   5.1385162e-7
11  4.7885350e5  4.7885359e5   1.8458927e-7
12  3.0583779e5  3.0583781e5   7.0650601e-8
13  2.0537268e5  2.0537272e5   1.9724071e-7
14  1.3581347e5  1.3581349e5   1.4894264e-7
15  9.3481731e4  9.3481745e4   1.5065258e-7
16  6.1789306e4  6.1789317e4   1.7921741e-7
17  4.2234799e4  4.2234800e4   3.7378332e-9
18  2.9446997e4  2.9446998e4   2.7827810e-8
19  2.0711158e4  2.0711158e4   3.6083760e-9
20  1.4399288e4  1.4399288e4  -1.6959262e-9
21  1.0204248e4  1.0204248e4  2.1560531e-13
22  7.0541335e3  7.0541335e3  2.0916602e-13
23  5.1274795e3  5.1274798e3   5.3440854e-8
24  3.5785021e3  3.5785021e3 -7.3925543e-10
25  2.5497806e3  2.5497806e3  2.0650148e-13
26  1.8525134e3  1.8525134e3 -7.1513917e-10
27  1.3124163e3  1.3124164e3   6.5466644e-8
28  9.3353298e2  9.3353296e2  -2.4834690e-8
29  6.8182565e2  6.8182566e2   1.9250868e-8
30  4.8152109e2  4.8152109e2  2.0583535e-13
31  3.4599110e2  3.4599110e2  2.0228264e-13
32  2.4853801e2  2.4853801e2  1.9961810e-13
33  1.8153858e2  1.8153858e2  1.9295676e-13
34  1.2893001e2  1.2893001e2  1.8895996e-13
35  9.2863897e1  9.2863897e1 -1.2712831e-11
36  6.8219858e1  6.8219863e1   6.7323209e-8
37  4.8499048e1  4.8499049e1   2.6932878e-9
38  3.5260803e1  3.5260804e1   3.1629175e-8
39  2.4875207e1  2.4875207e1   1.5327475e-9
40  1.8248343e1  1.8248344e1   9.2649788e-8
41  1.2848351e1  1.2848351e1  -3.4248875e-8
42  9.3637433e0  9.3637433e0   1.0366541e-9
43  6.5025357e0  6.5025357e0   2.9029490e-9
44  4.7299678e0  4.7299683e0   1.1183387e-7
45  3.3169392e0  3.3169393e0   1.6485298e-8
46  2.3366424e0  2.3366426e0   7.2538949e-8
47  1.6375345e0  1.6375347e0   1.3000311e-7
48  1.1395783e0  1.1395784e0   8.0374843e-8
49 7.7573740e-1 7.7573745e-1   6.3753925e-8
50 5.3385646e-1 5.3385655e-1   1.7113346e-7
51 3.5431304e-1 3.5431287e-1  -4.7859888e-7
52 2.3984658e-1 2.3984660e-1   7.3858889e-8
53 1.5742653e-1 1.5742660e-1   4.7452613e-7
54 1.0381680e-1 1.0381681e-1   1.3036277e-7
55 6.7367194e-2 6.7367194e-2   4.7924602e-9
56 4.1754029e-2 4.1754049e-2   4.8244591e-7
57 2.5167475e-2 2.5167484e-2   3.3871366e-7
58 1.5430966e-2 1.5430971e-2   3.0639168e-7
59 8.9218263e-3 8.9218286e-3   2.6013356e-7
60 5.0265483e-3 5.0265486e-3   6.8485890e-8
61 2.9066861e-3 2.9066853e-3  -2.9518261e-7
62 1.4089133e-3 1.4089139e-3   4.5456056e-7
63 6.7833493e-4 6.7833993e-4   7.3775731e-6
64 1.0497365e-4 1.0497395e-4   2.8590177e-6
";

#[cfg(feature = "applgrid")]
const EXPORT_NNLO_AK5_PTJ_NO_DISCARD_FAILS_STR: &str = "Error: factorization scale muf2 = 33634.5450347851 not found in APPLgrid; try exporting with `--discard-non-matching-values`
";

#[test]
#[cfg(feature = "applgrid")]
fn export_nnlo_ak5_ptj_discard_non_matching_values() {
    let output = NamedTempFile::new("nnlo.ak5_ptj.appl").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "export",
            "--discard-non-matching-values",
            "--accuracy=1e-5",
            "../test-data/nnlo.ak5_ptj.pineappl.lz4",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(EXPORT_NNLO_AK5_PTJ_STR));
}

#[test]
#[cfg(feature = "applgrid")]
fn export_nnlo_ak5_ptj_discard_non_matching_scales_alias() {
    let output = NamedTempFile::new("nnlo.ak5_ptj.appl").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "export",
            "--discard-non-matching-scales",
            "--accuracy=1e-5",
            "../test-data/nnlo.ak5_ptj.pineappl.lz4",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(EXPORT_NNLO_AK5_PTJ_STR));
}

#[test]
#[cfg(feature = "applgrid")]
fn export_nnlo_ak5_ptj_discard_non_matching_momentum_alias() {
    let output = NamedTempFile::new("nnlo.ak5_ptj.appl").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "export",
            "--discard-non-matching-momentum",
            "--accuracy=1e-5",
            "../test-data/nnlo.ak5_ptj.pineappl.lz4",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(EXPORT_NNLO_AK5_PTJ_STR));
}

#[test]
#[cfg(feature = "applgrid")]
fn export_nnlo_ak5_ptj_no_discard_fails() {
    let output = NamedTempFile::new("converted.appl").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "export",
            "../test-data/nnlo.ak5_ptj.pineappl.lz4",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(EXPORT_NNLO_AK5_PTJ_NO_DISCARD_FAILS_STR);
}

#[test]
#[cfg(feature = "applgrid")]
fn export_dis_applgrid() {
    let output1 = NamedTempFile::new("remapped.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--set-bins=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45",
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
