use assert_cmd::Command;

#[cfg(any(feature = "applgrid", feature = "fastnlo", feature = "fktable"))]
use assert_fs::NamedTempFile;

#[cfg(feature = "fastnlo")]
const HELP_STR: &str = "Converts APPLgrid/fastNLO/FastKernel files to PineAPPL grids

Usage: pineappl import [OPTIONS] <INPUT> <OUTPUT> <CONV_FUNS>

Arguments:
  <INPUT>      Path to the input grid
  <OUTPUT>     Path to the converted grid
  <CONV_FUNS>  LHAPDF ID(s) or name of the PDF(s)/FF(s) to check the converted grid with

Options:
      --alpha <ALPHA>        LO coupling power in alpha [default: 0]
      --accuracy <ACCURACY>  Relative threshold between the table and the converted grid when comparison fails [default: 1e-10]
  -s, --scales <SCALES>      Set the number of scale variations to compare with if they are available [default: 7] [possible values: 1, 3, 7, 9]
      --fnlo-mur <FNLO_MUR>  If importing a fastNLO flexible-scale grid, use the specified functional form for the renormalization scale [possible values: kScale1, kScale2, kQuadraticSum, kQuadraticMean, kQuadraticSumOver4, kLinearMean, kLinearSum, kScaleMax, kScaleMin, kProd, kS2plusS1half, kPow4Sum, kWgtAvg, kS2plusS1fourth, kExpProd2]
      --fnlo-muf <FNLO_MUF>  If importing a fastNLO flexible-scale grid, use the specified functional form for the factorization scale [possible values: kScale1, kScale2, kQuadraticSum, kQuadraticMean, kQuadraticSumOver4, kLinearMean, kLinearSum, kScaleMax, kScaleMin, kProd, kS2plusS1half, kPow4Sum, kWgtAvg, kS2plusS1fourth, kExpProd2]
      --digits-abs <ABS>     Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>     Set the number of fractional digits shown for relative numbers [default: 7]
      --no-optimize          Do not optimize converted grid
      --dis-pid <DIS_PID>    Particle ID for the non-hadronic initial states if it cannot be determined from the grid [default: 11]
  -h, --help                 Print help
";

#[cfg(not(feature = "fastnlo"))]
const HELP_STR: &str = "Converts APPLgrid/fastNLO/FastKernel files to PineAPPL grids

Usage: pineappl import [OPTIONS] <INPUT> <OUTPUT> <CONV_FUNS>

Arguments:
  <INPUT>      Path to the input grid
  <OUTPUT>     Path to the converted grid
  <CONV_FUNS>  LHAPDF ID(s) or name of the PDF(s)/FF(s) to check the converted grid with

Options:
      --alpha <ALPHA>        LO coupling power in alpha [default: 0]
      --accuracy <ACCURACY>  Relative threshold between the table and the converted grid when comparison fails [default: 1e-10]
  -s, --scales <SCALES>      Set the number of scale variations to compare with if they are available [default: 7] [possible values: 1, 3, 7, 9]
      --fnlo-mur <FNLO_MUR>  If importing a fastNLO flexible-scale grid, use the specified functional form for the renormalization scale
      --fnlo-muf <FNLO_MUF>  If importing a fastNLO flexible-scale grid, use the specified functional form for the factorization scale
      --digits-abs <ABS>     Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>     Set the number of fractional digits shown for relative numbers [default: 7]
      --no-optimize          Do not optimize converted grid
      --dis-pid <DIS_PID>    Particle ID for the non-hadronic initial states if it cannot be determined from the grid [default: 11]
  -h, --help                 Print help
";

#[cfg(feature = "fastnlo")]
const IMPORT_FIX_GRID_STR: &str = "b   PineAPPL     fastNLO      rel. diff
-+------------+------------+--------------
0 2.9158424e-4 2.9158424e-4 -2.9976022e-15
1 2.4657895e-4 2.4657895e-4 -2.8865799e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.2754182e1  8.2754182e1 -1.3544721e-14 1.3544721e-14
1  3.6097335e1  3.6097335e1 -6.8833828e-15 8.8817842e-15
2  8.0048746e0  8.0048746e0  5.3290705e-15 6.8833828e-15
3 9.4319392e-1 9.4319392e-1  5.5511151e-15 5.5511151e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_SCALE_1_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.1965747e1  8.1965747e1  1.5543122e-15 7.6605389e-15
1  3.6115068e1  3.6115068e1 -3.1086245e-15 1.4321877e-14
2  8.1057136e0  8.1057136e0  8.8817842e-16 5.7731597e-15
3 9.5444782e-1 9.5444782e-1  5.5511151e-15 5.5511151e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_SCALE_2_STR: &str = "b   PineAPPL     fastNLO      rel. diff   svmaxreldiff
-+------------+------------+-------------+-------------
0  8.3815533e1  8.3815533e1 4.8849813e-15 4.8849813e-15
1  3.6084994e1  3.6084994e1 2.6645353e-15 7.7715612e-15
2  7.8842272e0  7.8842272e0 1.9984014e-15 4.3298698e-15
3 9.1960866e-1 9.1960866e-1 3.1086245e-15 5.3290705e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_QUADRATIC_SUM_STR: &str =
    "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.1098571e1  8.1098571e1 -4.7739590e-15 7.3274720e-15
1  3.5222658e1  3.5222658e1  1.1102230e-15 6.6613381e-15
2  7.7939468e0  7.7939468e0  1.9984014e-15 4.5519144e-15
3 9.1540624e-1 9.1540624e-1 -5.9952043e-15 7.7715612e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_QUADRATIC_MEAN_STR: &str =
    "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.2712488e1  8.2712488e1  2.2204460e-16 1.0214052e-14
1  3.6091182e1  3.6091182e1 -7.7715612e-16 5.9952043e-15
2  7.9809031e0  7.9809031e0 -6.9944051e-15 9.5479180e-15
3 9.3467326e-1 9.3467326e-1  6.6613381e-16 2.4424907e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_5_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.4122343e1  8.4122343e1  8.6597396e-15 8.6597396e-15
1  3.6813708e1  3.6813708e1  5.9952043e-15 7.3274720e-15
2  8.1178188e0  8.1178188e0 -1.1102230e-15 1.3322676e-14
3 9.5090947e-1 9.5090947e-1  6.6613381e-15 6.6613381e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_6_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.2853942e1  8.2853942e1    0.0000000e0 5.1070259e-15
1  3.6103118e1  3.6103118e1  7.1054274e-15 9.2148511e-15
2  8.0161351e0  8.0161351e0 -2.6645353e-15 1.1324275e-14
3 9.4536395e-1 9.4536395e-1  5.1070259e-15 5.8841820e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_7_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  7.9163980e1  7.9163980e1  3.7747583e-15 1.2323476e-14
1  3.4313126e1  3.4313126e1 -4.7739590e-15 9.2148511e-15
2  7.7006079e0  7.7006079e0 -3.2196468e-15 1.1546319e-14
3 9.2392932e-1 9.2392932e-1 -4.3298698e-15 5.5511151e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_8_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.1965747e1  8.1965747e1  1.5543122e-15 7.6605389e-15
1  3.5968167e1  3.5968167e1    0.0000000e0 7.7715612e-15
2  7.9289155e0  7.9289155e0  1.3322676e-15 1.0436096e-14
3 9.3523838e-1 9.3523838e-1 -2.6645353e-15 9.5479180e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_9_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.3760665e1  8.3760665e1  3.5527137e-15 9.5479180e-15
1  3.6243722e1  3.6243722e1 -7.4384943e-15 1.7430501e-14
2  8.1057136e0  8.1057136e0  8.8817842e-16 5.7731597e-15
3 9.5444782e-1 9.5444782e-1  5.5511151e-15 5.5511151e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_10_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  6.9429217e1  6.9429217e1 -2.6645353e-15 5.1070259e-15
1  2.9273448e1  2.9273448e1 -8.8817842e-16 5.7731597e-15
2  6.6031456e0  6.6031456e0  2.6645353e-15 5.5511151e-15
3 8.2741590e-1 8.2741590e-1 -8.7707619e-15 8.7707619e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_11_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.2113143e1  8.2113143e1 -2.9976022e-15 6.4392935e-15
1  3.5603233e1  3.5603233e1 -2.2204460e-16 1.2212453e-14
2  7.8899185e0  7.8899185e0 -4.8849813e-15 8.4376950e-15
3 9.3402696e-1 9.3402696e-1 -2.3314684e-15 6.8833828e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_12_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.1720563e1  8.1720563e1 -8.5487173e-15 1.0658141e-14
1  3.5668000e1  3.5668000e1  3.9968029e-15 9.4368957e-15
2  7.9130511e0  7.9130511e0 -1.6653345e-15 6.4392935e-15
3 9.3503500e-1 9.3503500e-1 -2.2204460e-16 5.1070259e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_13_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.2403226e1  8.2403226e1 -4.6629367e-15 5.5511151e-15
1  3.6074907e1  3.6074907e1 -4.4408921e-16 5.8841820e-15
2  7.9668487e0  7.9668487e0 -6.6613381e-16 1.0880186e-14
3 9.3711914e-1 9.3711914e-1 -3.3306691e-15 8.6597396e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_14_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  8.2850540e1  8.2850540e1  6.8833828e-15 7.1054274e-15
1  3.5828674e1  3.5828674e1  2.6645353e-15 1.0103030e-14
2  7.9087501e0  7.9087501e0 -8.1046281e-15 8.1046281e-15
3 9.3462321e-1 9.3462321e-1  4.4408921e-16 8.2156504e-15
";

#[cfg(feature = "fastnlo")]
const IMPORT_FLEX_GRID_15_STR: &str = "b   PineAPPL     fastNLO      rel. diff    svmaxreldiff
-+------------+------------+--------------+-------------
0  6.6997861e1  6.6997861e1  5.7731597e-15 1.1879386e-14
1  2.6049196e1  2.6049196e1 -7.7715612e-16 1.3100632e-14
2  5.2022797e0  5.2022797e0  8.4376950e-15 8.4376950e-15
3 7.2427500e-1 7.2427500e-1  1.9984014e-15 9.5479180e-15
";

#[cfg(feature = "fktable")]
const IMPORT_DIS_FKTABLE_STR: &str = "b   x1       diff     
    []        []      
--+--+--+-------------
 0  0  1   1.8900584e0
 1  1  2   1.4830883e0
 2  2  3   1.1497012e0
 3  3  4  8.7974506e-1
 4  4  5  7.0882550e-1
 5  5  6  5.7345845e-1
 6  6  7  4.7744833e-1
 7  7  8  4.1037984e-1
 8  8  9  4.0362470e-1
 9  9 10  4.2613006e-1
10 10 11  3.7669466e-1
11 11 12  2.9572989e-1
12 12 13  2.0869778e-1
13 13 14  1.2602675e-1
14 14 15  6.4220769e-2
15 15 16  2.5434367e-2
16 16 17  7.6070428e-3
17 17 18  2.1848546e-3
18 18 19  6.2309735e-4
19 19 20 -1.0496472e-4
";

#[cfg(feature = "fktable")]
const IMPORT_HADRONIC_FKTABLE_STR: &str = "b x1     diff    
  []      []     
-+-+-+-----------
0 0 1 7.7624461e2
";

#[cfg(feature = "applgrid")]
const IMPORT_PHOTON_GRID_STR: &str = "b   PineAPPL     APPLgrid    rel. diff
-+------------+------------+-----------
0 5.5621307e-4 5.5621307e-4 0.0000000e0
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
1  6.2847078e2 6.2847078e2    0.0000000e0
2  6.3163323e2 6.3163323e2  2.2204460e-16
3  6.3586556e2 6.3586556e2  2.2204460e-16
4  6.4139163e2 6.4139163e2  1.7763568e-15
5  6.4848088e2 6.4848088e2 -2.6645353e-15
6  6.5354150e2 6.5354150e2 -3.6637360e-15
7  6.5377566e2 6.5377566e2 -1.7763568e-15
8  6.5094729e2 6.5094729e2  1.7763568e-15
9  6.3588760e2 6.3588760e2  2.2204460e-15
10 5.9810718e2 5.9810718e2  2.6645353e-15
";

const IMPORT_FILE_FORMAT_FAILURE_STR: &str = "Error: could not detect file format
";

#[cfg(feature = "fastnlo")]
const IMPORT_GRID_COMPARISON_FAILURE_STR: &str = "Error: grids are different
";

#[cfg(feature = "applgrid")]
const IMPORT_DIS_APPLGRID_STR: &str = "b   PineAPPL     APPLgrid     rel. diff
-+------------+------------+--------------
0 9.3514881e-2 9.3514881e-2 -3.3306691e-16
1 3.9993061e-2 3.9993061e-2  2.2204460e-16
2 1.3593440e-2 1.3593440e-2 -2.2204460e-16
3 2.0825199e-3 2.0825199e-3 -4.4408921e-16
";

#[cfg(feature = "fastnlo")]
const IMPORT_DOUBLE_HADRONIC_FASTNLO_STR: &str =
    "b    PineAPPL     fastNLO      rel. diff    svmaxreldiff
--+------------+------------+--------------+-------------
0   9.6382069e5  9.6382069e5  4.4408921e-16 8.3266727e-15
1   3.7342594e5  3.7342594e5  1.7985613e-14 1.9095836e-14
2   1.4195038e5  1.4195038e5 -1.0880186e-14 2.2648550e-14
3   5.7043791e4  5.7043791e4  4.2188475e-15 7.9936058e-15
4   2.3327746e4  2.3327746e4  8.4376950e-15 1.2101431e-14
5   1.0495603e4  1.0495603e4  1.3100632e-14 1.7985613e-14
6   4.8153483e3  4.8153483e3 -1.6098234e-14 2.9753977e-14
7   2.2957587e3  2.2957587e3  4.8849813e-15 3.0642155e-14
8   1.1142545e3  1.1142545e3 -2.4424907e-15 1.5765167e-14
9   5.3699925e2  5.3699925e2 -6.5503158e-15 1.8429702e-14
10  2.5460314e2  2.5460314e2 -7.6605389e-15 1.3544721e-14
11  1.1847638e2  1.1847638e2  1.0658141e-14 1.2656542e-14
12  5.7567355e1  5.7567355e1 -2.9976022e-15 9.2148511e-15
13  2.7189719e1  2.7189719e1  1.1102230e-15 1.5543122e-14
14  1.2791922e1  1.2791922e1 -6.9944051e-15 1.2656542e-14
15  5.8346996e0  5.8346996e0  2.8865799e-15 1.4988011e-14
16  2.6521590e0  2.6521590e0  7.3274720e-15 1.4765966e-14
17  1.1726035e0  1.1726035e0  1.3100632e-14 1.4432899e-14
18 4.8823596e-1 4.8823596e-1  8.6597396e-15 1.3655743e-14
19 1.9564964e-1 1.9564964e-1 -4.4408921e-15 1.1102230e-14
20 2.0326950e-2 2.0326950e-2  6.6613381e-15 1.3211654e-14
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["import", "--help"])
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
        .args([
            "import",
            "../test-data/NJetEvents_0-0-2.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FIX_GRID_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_scale_1() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kScale1",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_SCALE_1_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_scale_2() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kScale2",
            "--fnlo-muf=kScale2",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_SCALE_2_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_quadratic_sum() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kQuadraticSum",
            "--fnlo-muf=kQuadraticSum",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(
            IMPORT_FLEX_GRID_QUADRATIC_SUM_STR,
        ));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_quadratic_mean() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kQuadraticMean",
            "--fnlo-muf=kQuadraticMean",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(
            IMPORT_FLEX_GRID_QUADRATIC_MEAN_STR,
        ));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_5() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kQuadraticSumOver4",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_5_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_6() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kLinearMean",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_6_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_7() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kLinearSum",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_7_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_8() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kScaleMax",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_8_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_9() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kScaleMin",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_9_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_10() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kProd",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_10_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_11() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kS2plusS1half",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_11_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_12() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kPow4Sum",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_12_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_13() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kWgtAvg",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_13_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_14() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kS2plusS1fourth",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_14_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_flex_grid_15() {
    let output = NamedTempFile::new("converted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-incjets-fnlo-arxiv-0706.3722-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
            "--fnlo-mur=kExpProd2",
            "--fnlo-muf=kScale1",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_FLEX_GRID_15_STR));
}

#[test]
#[cfg(feature = "fktable")]
fn import_dis_fktable() {
    use pineappl::fk_table::FkTable;
    use pineappl::grid::Grid;
    use std::fs::File;

    let output = NamedTempFile::new("converted3.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/FK_POSXDQ.dat",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout("file was converted, but we cannot check the conversion for this type\n");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(IMPORT_DIS_FKTABLE_STR);

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["read", "--fktable", output.path().to_str().unwrap()])
        .assert()
        .success()
        .stdout("yes\n");

    let fk_table =
        FkTable::try_from(Grid::read(File::open(output.path()).unwrap()).unwrap()).unwrap();

    // TODO: this should ideally be a unit test, but we need an FK table that we don't convert

    assert_eq!(fk_table.muf2(), 1.65 * 1.65);
    assert_eq!(
        fk_table.x_grid(),
        [
            4.949999999999766e-7,
            6.07738163568489e-7,
            7.461527647610767e-7,
            9.160916270024875e-7,
            1.1247343892713325e-6,
            1.3808958061045167e-6,
            1.6953981117803826e-6,
            2.0815281713594856e-6,
            2.5555987186183014e-6,
            3.1376373712233846e-6,
            3.852232500999415e-6,
            4.729571647582583e-6,
            5.806716273568501e-6,
            7.1291666068230954e-6,
            8.752782514758356e-6,
            1.0746141311783372e-5,
            1.3193431732093381e-5,
            1.6198005757901242e-5,
            1.988673749378987e-5,
            2.4415371932960747e-5,
            2.9975087617735706e-5,
            3.680054747936805e-5,
            4.517977350729615e-5,
            5.54662556790657e-5,
            6.809379655659362e-5,
            8.35947033863248e-5,
            0.0001026220732293567,
            0.0001259770779412776,
            0.0001546423495893497,
            0.00018982279841601318,
            0.00023299547017757838,
            0.00028597037276041366,
            0.0003509645774833949,
            0.00043069233015845184,
            0.0005284743881553908,
            0.0006483703223240056,
            0.0007953380646144586,
            0.0009754255030559555,
            0.0011959993572648567,
            0.0014660168019741714,
            0.0017963451801748858,
            0.002200134425171399,
            0.0026932451647119623,
            0.0032947324788041603,
            0.004027380393488671,
            0.00491827481991599,
            0.005999392230815202,
            0.007308167572260819,
            0.008887987953627191,
            0.010788539757904647,
            0.013065918705455737,
            0.015782399751004605,
            0.0190057630326846,
            0.022808090849411308,
            0.0272639948527007,
            0.032448303490771324,
            0.038433330036855784,
            0.04528593412657552,
            0.053064660208759794,
            0.061817260572283074,
            0.07157887493476887,
            0.08237104713583011,
            0.09420163447031454,
            0.10706553869960851,
            0.12094608945898123,
            0.1358168578945569,
            0.1516436724678353,
            0.16838663922883584,
            0.18600201927146398,
            0.20444387120690388,
            0.22366541552584715,
            0.2436201152370417,
            0.26426249207236646,
            0.2855487113993947,
            0.307436974454643,
            0.32988775640071194,
            0.35286392532847166,
            0.3763307724088937,
            0.4002559780692545,
            0.42460953399658175,
            0.4493636362927147,
            0.47449256134271817,
            0.4999725329101106,
            0.5257815865664291,
            0.5518994357134863,
            0.5783073420529784,
            0.6049879923226599,
            0.6319253823550638,
            0.6591047089777237,
            0.6865122698903083,
            0.7141353713969553,
            0.741962243702359,
            0.7699819633745686,
            0.798184382520784,
            0.8265600641917732,
            0.8551002235319206,
            0.8837966741980419,
            0.9126417795942889,
            0.9416284084927907,
            0.9707498946430192
        ]
    );

    let table = fk_table.table();

    assert_eq!(table.dim(), (20, 9, 100, 1));
    assert_eq!(
        table
            .indexed_iter()
            .filter(|(_, value)| **value != 0.0)
            .take(10)
            .collect::<Vec<_>>(),
        [
            ((0, 0, 0, 0), &4.506605409085538e-8),
            ((0, 0, 1, 0), &1.8561090273141668e-8),
            ((0, 0, 2, 0), &-3.3821015317570252e-9),
            ((0, 0, 3, 0), &1.980084314325426e-9),
            ((0, 0, 4, 0), &2.187815687938248e-9),
            ((0, 0, 5, 0), &1.3280152778522626e-9),
            ((0, 0, 6, 0), &1.3848470515483116e-9),
            ((0, 0, 7, 0), &1.5145898293299224e-9),
            ((0, 0, 8, 0), &1.6942313031679552e-9),
            ((0, 0, 9, 0), &1.9734220063025288e-9),
        ]
    );
}

#[test]
#[cfg(feature = "fktable")]
fn import_hadronic_fktable() {
    use float_cmp::assert_approx_eq;
    use lhapdf::Pdf;
    use pineappl::convolutions::Convolution;
    use pineappl::convolutions::LumiCache;
    use pineappl::fk_table::{FkAssumptions, FkTable};
    use pineappl::grid::Grid;
    use std::fs::File;

    let output = NamedTempFile::new("converted4.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/FK_ATLASTTBARTOT13TEV.dat",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout("file was converted, but we cannot check the conversion for this type\n");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(IMPORT_HADRONIC_FKTABLE_STR);

    // TODO: this should ideally be a unit test, but we need an FK table that we don't convert

    let file = File::open(output.path()).unwrap();
    let grid = Grid::read(file).unwrap();

    let pdf = Pdf::with_setname_and_member("NNPDF31_nlo_as_0118_luxqed", 0).unwrap();
    let mut xfx = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut alphas = |_| 0.0;
    let mut lumi_cache = LumiCache::with_one(2212, &mut xfx, &mut alphas);
    let results = grid.convolve(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0)]);

    let mut fk_table = FkTable::try_from(grid).unwrap();
    let table = fk_table.table();

    assert_eq!(table.dim(), (1, 45, 30, 30));
    assert_eq!(
        table
            .indexed_iter()
            .filter(|(_, value)| **value != 0.0)
            .take(10)
            .collect::<Vec<_>>(),
        [
            ((0, 0, 0, 22), &1.735017465956091e-13),
            ((0, 0, 0, 23), &2.122735969332162e-12),
            ((0, 0, 0, 24), &-1.377087610955672e-12),
            ((0, 0, 0, 25), &3.348327389530187e-10),
            ((0, 0, 0, 26), &1.7760262237114066e-9),
            ((0, 0, 0, 27), &4.061604348778655e-9),
            ((0, 0, 0, 28), &6.702574861447763e-9),
            ((0, 0, 0, 29), &8.930615742623519e-9),
            ((0, 0, 1, 19), &1.5892116835024906e-12),
            ((0, 0, 1, 20), &8.160271984903452e-12),
        ]
    );

    assert_eq!(fk_table.bins(), 1);
    assert_eq!(fk_table.bin_normalizations(), [1.0]);
    assert_eq!(fk_table.bin_dimensions(), 1);
    assert_eq!(fk_table.bin_left(0), [0.0]);
    assert_eq!(fk_table.bin_right(0), [1.0]);
    assert_eq!(
        fk_table.grid().convolutions(),
        [Convolution::UnpolPDF(2212), Convolution::UnpolPDF(2212)]
    );
    let channels = fk_table.channels();
    assert_eq!(
        channels,
        [
            [100, 100],
            [100, 21],
            [100, 200],
            [100, 203],
            [100, 208],
            [100, 215],
            [100, 103],
            [100, 108],
            [100, 115],
            [21, 21],
            [21, 200],
            [21, 203],
            [21, 208],
            [21, 215],
            [21, 103],
            [21, 108],
            [21, 115],
            [200, 200],
            [200, 203],
            [200, 208],
            [200, 215],
            [200, 103],
            [200, 108],
            [200, 115],
            [203, 203],
            [203, 208],
            [203, 215],
            [203, 103],
            [203, 108],
            [203, 115],
            [208, 208],
            [208, 215],
            [208, 103],
            [208, 108],
            [208, 115],
            [215, 215],
            [215, 103],
            [215, 108],
            [215, 115],
            [103, 103],
            [103, 108],
            [103, 115],
            [108, 108],
            [108, 115],
            [115, 115]
        ]
    );
    assert_eq!(fk_table.muf2(), 1.65 * 1.65);
    assert_eq!(
        fk_table.x_grid(),
        [
            0.000289610509361025,
            0.00046359018311242675,
            0.0007416226869940695,
            0.0011852237616640315,
            0.0018911851231829696,
            0.0030101731988996802,
            0.004772784166476982,
            0.007522794084020634,
            0.011752513921971337,
            0.018125783312196084,
            0.02746253640499583,
            0.0406567654515178,
            0.058525420074312585,
            0.08163793076418926,
            0.11020879289284256,
            0.14410230325000434,
            0.1829266515794375,
            0.22615523325198245,
            0.27322658065490807,
            0.32360666435374563,
            0.3768187300969618,
            0.43245213927164394,
            0.4901600611528503,
            0.549652367583778,
            0.6106872012258386,
            0.6730628689129965,
            0.736610712155086,
            0.8011891117676506,
            0.8666785643244361,
            0.9329776879738404
        ]
    );

    assert_eq!(results, fk_table.convolve(&mut lumi_cache, &[], &[]));

    fk_table.optimize(FkAssumptions::Nf6Ind);
    assert_eq!(fk_table.channels(), channels);
    assert_approx_eq!(
        f64,
        results[0],
        fk_table.convolve(&mut lumi_cache, &[], &[])[0],
        ulps = 4
    );
    fk_table.optimize(FkAssumptions::Nf6Sym);
    assert_eq!(fk_table.channels(), channels);
    assert_approx_eq!(
        f64,
        results[0],
        fk_table.convolve(&mut lumi_cache, &[], &[])[0],
        ulps = 4
    );
    fk_table.optimize(FkAssumptions::Nf5Ind);
    assert_eq!(fk_table.channels(), channels);
    assert_approx_eq!(
        f64,
        results[0],
        fk_table.convolve(&mut lumi_cache, &[], &[])[0]
    );
    fk_table.optimize(FkAssumptions::Nf5Sym);
    assert_eq!(fk_table.channels(), channels);
    assert_approx_eq!(
        f64,
        results[0],
        fk_table.convolve(&mut lumi_cache, &[], &[])[0]
    );
    fk_table.optimize(FkAssumptions::Nf4Ind);
    assert_eq!(fk_table.channels(), channels);
    assert_approx_eq!(
        f64,
        results[0],
        fk_table.convolve(&mut lumi_cache, &[], &[])[0]
    );

    fk_table.optimize(FkAssumptions::Nf4Sym);
    assert_eq!(
        fk_table.channels(),
        [
            [100, 100],
            [100, 21],
            [100, 203],
            [100, 208],
            [100, 200],
            [100, 103],
            [100, 108],
            [100, 115],
            [21, 21],
            [21, 203],
            [21, 208],
            [21, 200],
            [21, 103],
            [21, 108],
            [21, 115],
            [200, 203],
            [200, 208],
            [203, 203],
            [203, 208],
            [203, 103],
            [203, 108],
            [203, 115],
            [208, 208],
            [208, 103],
            [208, 108],
            [208, 115],
            [200, 200],
            [200, 103],
            [200, 108],
            [200, 115],
            [103, 103],
            [103, 108],
            [103, 115],
            [108, 108],
            [108, 115],
            [115, 115]
        ]
    );
    fk_table.optimize(FkAssumptions::Nf3Ind);
    assert_eq!(
        fk_table.channels(),
        [
            [100, 21],
            [100, 203],
            [100, 208],
            [100, 200],
            [100, 103],
            [100, 108],
            [21, 21],
            [21, 203],
            [21, 208],
            [21, 200],
            [21, 103],
            [21, 108],
            [200, 203],
            [200, 208],
            [203, 203],
            [203, 208],
            [203, 103],
            [203, 108],
            [208, 208],
            [208, 103],
            [208, 108],
            [200, 200],
            [200, 103],
            [200, 108],
            [103, 103],
            [103, 108],
            [108, 108],
            [100, 100]
        ]
    );
    fk_table.optimize(FkAssumptions::Nf3Sym);
    assert_eq!(
        fk_table.channels(),
        [
            [100, 21],
            [100, 203],
            [100, 200],
            [100, 103],
            [100, 108],
            [21, 21],
            [21, 203],
            [21, 200],
            [21, 103],
            [21, 108],
            [200, 203],
            [203, 203],
            [203, 103],
            [203, 108],
            [200, 200],
            [200, 103],
            [200, 108],
            [103, 103],
            [103, 108],
            [108, 108],
            [100, 100],
        ]
    );
}

#[test]
#[cfg(feature = "applgrid")]
fn import_photon_grid() {
    let output = NamedTempFile::new("converted5.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/LHCBWZMU7TEV_PI_part1.appl",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_PHOTON_GRID_STR));
}

#[test]
#[cfg(feature = "applgrid")]
fn import_applgrid() {
    let output = NamedTempFile::new("converted6.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/ATLASWPT11-Wplus_tot.appl",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_APPLGRID_STR));
}

#[test]
#[cfg(feature = "applgrid")]
fn import_new_applgrid() {
    let output = NamedTempFile::new("converted7.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/atlas-atlas-wpm-arxiv-1109.5141-xsec001.appl",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_NEW_APPLGRID_STR));
}

#[test]
fn import_file_format_failure() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
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
        .args([
            "import",
            "--accuracy=0",
            "../test-data/NJetEvents_0-0-2.tab.gz",
            "this-file-wont-be-writen",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(IMPORT_GRID_COMPARISON_FAILURE_STR)
        .stdout(predicates::str::ends_with(IMPORT_FIX_GRID_STR));
}

#[test]
#[cfg(feature = "applgrid")]
fn import_dis_applgrid() {
    let output = NamedTempFile::new("converted8.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-h1-dijets-appl-arxiv-0010054-xsec000.appl",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(IMPORT_DIS_APPLGRID_STR));
}

#[test]
#[cfg(feature = "fastnlo")]
fn import_double_hadronic_fastnlo() {
    let output = NamedTempFile::new("converted9.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "import",
            "../test-data/applfast-atlas-dijets-fnlo-arxiv-1312.3524-xsec000.tab.gz",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(predicates::str::ends_with(
            IMPORT_DOUBLE_HADRONIC_FASTNLO_STR,
        ));
}
