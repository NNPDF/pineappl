use assert_cmd::Command;
use assert_fs::NamedTempFile;

const HELP_STR: &str = "pineappl-evolve 
Evolve a grid with an evolution kernel operator to an FK table

USAGE:
    pineappl evolve [OPTIONS] <INPUT> <EKO> <OUTPUT> <PDFSET>

ARGS:
    <INPUT>     Path to the input grid
    <EKO>       Path to the evolution kernel operator
    <OUTPUT>    Path to the converted grid
    <PDFSET>    LHAPDF id or name of the PDF set to check the converted grid with

OPTIONS:
        --accuracy <ACCURACY>    Relative threshold between the table and the converted grid when
                                 comparison fails [default: 1e-3]
        --digits-abs <ABS>       Set the number of fractional digits shown for absolute numbers
                                 [default: 7]
        --digits-rel <REL>       Set the number of fractional digits shown for relative numbers
                                 [default: 7]
    -h, --help                   Print help information
";

const E906NLO_BIN_00_STR: &str = "b   FkTable        Grid       rel. diff
-+------------+------------+-------------
0 1.0659807e-1 1.0657910e-1 -1.7800240e-4
1 1.0659807e-1 1.0657910e-1 -1.7800240e-4
2 1.0659807e-1 1.0657910e-1 -1.7800240e-4
3 1.0659807e-1 1.0657910e-1 -1.7800240e-4
4  3.2698655e0  3.2710283e0  3.5561489e-4
5  1.6039253e0  1.6047363e0  5.0563376e-4
";

const LHCB_WP_7TEV_STR: &str = "b   FkTable      Grid       rel. diff
-+-----------+-----------+-------------
0 7.7911994e2 7.7891206e2 -2.6680685e-4
1 7.1170872e2 7.1152118e2 -2.6351282e-4
2 6.1766223e2 6.1750067e2 -2.6157398e-4
3 4.9817054e2 4.9804168e2 -2.5866843e-4
4 3.7040220e2 3.7030899e2 -2.5164930e-4
5 2.5123713e2 2.5117697e2 -2.3945309e-4
6 1.1883712e2 1.1881217e2 -2.0989153e-4
7 2.9013827e1 2.9010212e1 -1.2459216e-4
";

const NUTEV_CC_NU_FE_SIGMARED_STR: &str = "b    FkTable      Grid       rel. diff
--+-----------+-----------+-------------
0  8.2920022e0 1.1046839e1  3.3222819e-1
1  1.3975037e1 1.4232414e1  1.8416931e-2
2  1.9422915e1 1.8434111e1 -5.0909145e-2
3  1.4891673e1 1.4867720e1 -1.6084417e-3
4  1.4086222e1 1.4088115e1  1.3438079e-4
5  2.2998409e1 2.2998443e1  1.4697729e-6
6  2.1899598e1 2.1899091e1 -2.3138103e-5
7  1.5578822e1 1.5578334e1 -3.1319781e-5
8  2.2226765e1 2.2225956e1 -3.6375906e-5
9  1.3669291e1 1.3669038e1 -1.8446253e-5
10 2.5946596e1 2.5945504e1 -4.2110358e-5
11 1.8363999e1 1.8363361e1 -3.4704937e-5
12 1.5819602e1 1.5819077e1 -3.3204532e-5
13 2.5959140e1 2.5957783e1 -5.2285210e-5
14 2.0457793e1 2.0456964e1 -4.0496689e-5
15 1.1393262e1 1.1392626e1 -5.5781822e-5
16 1.4106814e1 1.4106786e1 -1.9959280e-6
17 1.9464293e1 1.9463520e1 -3.9713070e-5
18 1.5721645e1 1.5721550e1 -6.0470703e-6
19 1.4033170e1 1.4033158e1 -8.5170032e-7
20 1.9366211e1 1.9365482e1 -3.7645378e-5
21 2.1174542e1 2.1173688e1 -4.0315154e-5
22 1.6606739e1 1.6606635e1 -6.2333155e-6
23 1.1354320e1 1.1353660e1 -5.8066443e-5
24 5.5456589e0 5.5455822e0 -1.3826910e-5
25 1.2257404e1 1.2256637e1 -6.2534396e-5
26 1.5838419e1 1.5838379e1 -2.5004885e-6
27 2.0844501e1 2.0844671e1  8.1437423e-6
28 1.1127050e1 1.1126496e1 -4.9723938e-5
29 1.3257033e1 1.3256687e1 -2.6046580e-5
30 1.6809659e1 1.6810591e1  5.5429049e-5
31 1.7616888e1 1.7618088e1  6.8090881e-5
32 5.8908612e0 5.8906267e0 -3.9802439e-5
33 1.2668418e1 1.2667676e1 -5.8532882e-5
34 6.1718100e0 6.1715800e0 -3.7259928e-5
35 1.6993748e1 1.6993216e1 -3.1306297e-5
36 5.5089321e0 5.5087361e0 -3.5579336e-5
37 1.2163126e1 1.2162352e1 -6.3616512e-5
38 1.2575328e1 1.2574501e1 -6.5776146e-5
39 5.9912886e0 5.9910701e0 -3.6479256e-5
40 5.5396399e0 5.5394441e0 -3.5343558e-5
41 1.2035714e1 1.2034895e1 -6.7972844e-5
42 5.1605061e0 5.1603299e0 -3.4143438e-5
43 5.3541916e0 5.3540055e0 -3.4757362e-5
44 4.9727490e0 4.9725844e0 -3.3105968e-5
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["evolve", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
#[ignore]
fn lhcb_wp_7tev() {
    let output = NamedTempFile::new("fktable1.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "evolve",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(LHCB_WP_7TEV_STR);
}

#[test]
#[ignore]
fn e906nlo_bin_00() {
    let output = NamedTempFile::new("fktable2.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "evolve",
            "../test-data/E906nlo_bin_00.pineappl.lz4",
            "../test-data/E906nlo_bin_00.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(E906NLO_BIN_00_STR);
}

#[test]
#[ignore]
fn nutev_cc_nu_fe_sigmared() {
    let output = NamedTempFile::new("fktable3.lz4").unwrap();

    // TODO: find out the reason why this evolution fails
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "evolve",
            "../test-data/NUTEV_CC_NU_FE_SIGMARED.pineappl.lz4",
            "../test-data/NUTEV_CC_NU_FE_SIGMARED.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nlo_as_01180",
        ])
        .assert()
        .failure()
        .stderr("Error: grids are different\n")
        .stdout(NUTEV_CC_NU_FE_SIGMARED_STR);
}
