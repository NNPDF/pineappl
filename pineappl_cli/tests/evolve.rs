#![allow(missing_docs)]
#![cfg(feature = "evolve")]

use assert_cmd::Command;
use assert_fs::NamedTempFile;

const HELP_STR: &str = "Evolve a grid with an evolution kernel operator to an FK table

Usage: pineappl evolve [OPTIONS] <INPUT> <EKO1,...> <OUTPUT> <CONV_FUNS>

Arguments:
  <INPUT>      Path to the input grid
  <EKO1,...>   Path to the evolution kernel operator(s)
  <OUTPUT>     Path to the converted grid
  <CONV_FUNS>  LHAPDF ID(s) or name of the PDF(s)/FF(s)

Options:
      --accuracy <ACCURACY>  Relative threshold between the table and the converted grid when comparison fails [default: 1e-3]
      --digits-abs <ABS>     Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>     Set the number of fractional digits shown for relative numbers [default: 7]
  -o, --orders <ORDERS>      Select which orders to evolve
      --xir <XIR>            Rescale the renormalization scale with this factor [default: 1]
      --xif <XIF>            Rescale the factorization scale with this factor [default: 1]
      --xia <XIA>            Rescale the fragmentation scale with this factor [default: 1]
  -h, --help                 Print help
";

const E906NLO_BIN_00_STR: &str = "b     Grid       FkTable      rel. diff
-+------------+------------+-------------
0 1.0659807e-1 1.0657904e-1 -1.7851986e-4
1  3.2698655e0  3.2711890e0  4.0477586e-4
2  1.6039253e0  1.6047566e0  5.1825508e-4
";

const LHCB_DY_8TEV_STR: &str = "b      Grid       FkTable      rel. diff
--+------------+------------+-------------
0   8.1932148e0  8.1911178e0 -2.5593854e-4
1   2.3943610e1  2.3937840e1 -2.4100020e-4
2   3.8438471e1  3.8429305e1 -2.3845440e-4
3   5.1375221e1  5.1362447e1 -2.4865338e-4
4   6.2453684e1  6.2437402e1 -2.6071449e-4
5   7.0822177e1  7.0803399e1 -2.6514173e-4
6   7.6444585e1  7.6424290e1 -2.6547832e-4
7   7.8730588e1  7.8709827e1 -2.6370710e-4
8   7.7456336e1  7.7436387e1 -2.5756344e-4
9   7.2375785e1  7.2357933e1 -2.4665727e-4
10  6.1963102e1  6.1948582e1 -2.3434293e-4
11  4.8188526e1  4.8177411e1 -2.3065387e-4
12  3.4441638e1  3.4433119e1 -2.4734868e-4
13  2.2221332e1  2.2215664e1 -2.5507029e-4
14  1.2601951e1  1.2599333e1 -2.0769978e-4
15  5.9418647e0  5.9412615e0 -1.0151683e-4
16  1.3823472e0  1.3825567e0  1.5158285e-4
17 5.0414246e-2 5.0460460e-2  9.1667268e-4
";

const LHCB_WP_7TEV_STR: &str = "b    Grid       FkTable     rel. diff
-+-----------+-----------+-------------
0 7.8752127e2 7.8731123e2 -2.6670347e-4
1 7.1872113e2 7.1853181e2 -2.6340965e-4
2 6.2322357e2 6.2306064e2 -2.6144415e-4
3 5.0216763e2 5.0203783e2 -2.5848913e-4
4 3.7314506e2 3.7305123e2 -2.5144024e-4
5 2.5302044e2 2.5295990e2 -2.3927343e-4
6 1.1971046e2 1.1968534e2 -2.0983475e-4
7 2.9272102e1 2.9268451e1 -1.2474967e-4
";

const LHCB_WP_7TEV_V2_STR: &str = "b        Grid             FkTable         rel. diff
-+------------------+------------------+-------------
0 7.87521267980686e2 7.87310643809286e2 -2.6745204e-4
1 7.18721130803477e2 7.18531231478480e2 -2.6421837e-4
2 6.23223573918485e2 6.23060099284591e2 -2.6230496e-4
3 5.02167629888729e2 5.02037373633690e2 -2.5938800e-4
4 3.73145056990031e2 3.73050898328477e2 -2.5233796e-4
5 2.53020442272921e2 2.52959682618899e2 -2.4013733e-4
6 1.19710459847744e2 1.19685254122495e2 -2.1055575e-4
7 2.92721022139301e1 2.92684433666511e1 -1.2499435e-4
";

const LHCB_WP_7TEV_V2_XIR2_STR: &str = "b        Grid             FkTable         rel. diff
-+------------------+------------------+-------------
0 7.76348332927370e2 7.76140378165194e2 -2.6786270e-4
1 7.08661998751250e2 7.08474448397818e2 -2.6465417e-4
2 6.14275560249818e2 6.14114173745311e2 -2.6272656e-4
3 4.94828199827837e2 4.94699640811431e2 -2.5980536e-4
4 3.67562574493549e2 3.67469675694897e2 -2.5274281e-4
5 2.49126427018341e2 2.49066510299154e2 -2.4050728e-4
6 1.17762540400323e2 1.17737720394934e2 -2.1076316e-4
7 2.87498912976683e1 2.87462994796563e1 -1.2493327e-4
";

const LHCB_WP_7TEV_V2_XIF_2_STR: &str = "b        Grid             FkTable         rel. diff
-+------------------+------------------+-------------
0 8.09024497135338e2 8.08801090895792e2 -2.7614274e-4
1 7.38692425698934e2 7.38491131004839e2 -2.7250136e-4
2 6.41024959047782e2 6.40851780258714e2 -2.7015920e-4
3 5.16685638376539e2 5.16547861676678e2 -2.6665479e-4
4 3.84050669911243e2 3.83951276776197e2 -2.5880214e-4
5 2.60476971252294e2 2.60412959132739e2 -2.4574963e-4
6 1.23243647450223e2 1.23217157841843e2 -2.1493691e-4
7 3.01346299826566e1 3.01308723713458e1 -1.2469412e-4
";

const LHCB_WP_7TEV_V2_XIF_2_ERROR_STR: &str =
    "Error: no operator for fac1 = '25825.775616000003' found
";

const LHCB_WP_7TEV_OPTIMIZED_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.8731123e2
1 2.25  2.5 7.1853181e2
2  2.5 2.75 6.2306064e2
3 2.75    3 5.0203783e2
4    3 3.25 3.7305123e2
5 3.25  3.5 2.5295990e2
6  3.5    4 1.1968534e2
7    4  4.5 2.9268451e1
";

const NUTEV_CC_NU_FE_SIGMARED_STR: &str = "b     Grid       FkTable     rel. diff
--+-----------+-----------+-------------
0  8.2920022e0 1.0954648e1  3.2111014e-1
1  1.3975037e1 1.4236502e1  1.8709416e-2
2  1.9422915e1 1.8431736e1 -5.1031421e-2
3  1.4891673e1 1.4867911e1 -1.5956063e-3
4  1.4086222e1 1.4087867e1  1.1676254e-4
5  2.2998409e1 2.2998414e1  2.3434412e-7
6  2.1899598e1 2.1899098e1 -2.2841392e-5
7  1.5578822e1 1.5578327e1 -3.1774716e-5
8  2.2226765e1 2.2225948e1 -3.6746158e-5
9  1.3669291e1 1.3669040e1 -1.8302270e-5
10 2.5946596e1 2.5945509e1 -4.1912903e-5
11 1.8363999e1 1.8363357e1 -3.4940966e-5
12 1.5819602e1 1.5819075e1 -3.3330771e-5
13 2.5959140e1 2.5957771e1 -5.2747753e-5
14 2.0457793e1 2.0456973e1 -4.0087031e-5
15 1.1393262e1 1.1392625e1 -5.5841791e-5
16 1.4106814e1 1.4106781e1 -2.3548459e-6
17 1.9464293e1 1.9463530e1 -3.9217573e-5
18 1.5721645e1 1.5721553e1 -5.8661712e-6
19 1.4033170e1 1.4033163e1 -4.8582905e-7
20 1.9366211e1 1.9365473e1 -3.8125440e-5
21 2.1174542e1 2.1173681e1 -4.0655919e-5
22 1.6606739e1 1.6606637e1 -6.1144617e-6
23 1.1354320e1 1.1353662e1 -5.7937744e-5
24 5.5456589e0 5.5455826e0 -1.3762742e-5
25 1.2257404e1 1.2256635e1 -6.2679973e-5
26 1.5838419e1 1.5838375e1 -2.8058203e-6
27 2.0844501e1 2.0844372e1 -6.2280076e-6
28 1.1127050e1 1.1126384e1 -5.9865346e-5
29 1.3257033e1 1.3256528e1 -3.8071417e-5
30 1.6809659e1 1.6810166e1  3.0140970e-5
31 1.7616888e1 1.7617649e1  4.3196882e-5
32 5.8908612e0 5.8906261e0 -3.9897551e-5
33 1.2668418e1 1.2667673e1 -5.8804601e-5
34 6.1718100e0 6.1715801e0 -3.7253769e-5
35 1.6993748e1 1.6993214e1 -3.1436426e-5
36 5.5089321e0 5.5087362e0 -3.5567214e-5
37 1.2163126e1 1.2162354e1 -6.3452698e-5
38 1.2575328e1 1.2574503e1 -6.5657683e-5
39 5.9912886e0 5.9910700e0 -3.6491168e-5
40 5.5396399e0 5.5394441e0 -3.5344028e-5
41 1.2035714e1 1.2034895e1 -6.7998362e-5
42 5.1605061e0 5.1603299e0 -3.4144995e-5
43 5.3541916e0 5.3540055e0 -3.4756819e-5
44 4.9727490e0 4.9725844e0 -3.3104317e-5
";

const CMS_TTB_8TEV_2D_TTM_TRAP_TOT_STR: &str = "b    Grid       FkTable     rel. diff
-+-----------+-----------+-------------
0 2.0680644e2 2.0666857e2 -6.6663644e-4
";

const STAR_WMWP_510GEV_WM_AL_POL: &str = "b    Grid       FkTable     rel. diff
-+-----------+-----------+-------------
0 3.2222870e2 3.2226654e2  1.1745654e-4
1 1.8038157e3 1.8037829e3 -1.8192479e-5
2 3.4767572e3 3.4762728e3 -1.3933339e-4
3 4.3157563e3 4.3154783e3 -6.4409623e-5
4 3.6443947e3 3.6443481e3 -1.2807044e-5
5 5.8386697e2 5.8336795e2 -8.5468266e-4
";

const ZEUS_2JET_STR: &str = "b     Grid       FkTable      rel. diff
-+------------+------------+-------------
0 8.8139165e-2 8.8034383e-2 -1.1888247e-3
1 3.3383261e-2 3.3326664e-2 -1.6953697e-3
2 3.6247796e-3 3.6162230e-3 -2.3605729e-3
";

const LHCB_WP_8TEV_STR: &str = "b    Grid       FkTable     rel. diff
-+-----------+-----------+-------------
0 8.8660824e2 8.8745467e2  9.5468156e-4
1 8.3324869e2 8.3388816e2  7.6744152e-4
2 7.4379285e2 7.4420759e2  5.5761143e-4
3 6.2114832e2 6.2135970e2  3.4030039e-4
4 4.8212545e2 4.8218796e2  1.2966015e-4
5 3.4357834e2 3.4355392e2 -7.1080989e-5
6 1.7271792e2 1.7266488e2 -3.0707061e-4
7 4.6738298e1 4.6715819e1 -4.8096830e-4
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["evolve", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn lhcb_wp_7tev() {
    let output = NamedTempFile::new("fktable1a.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "evolve",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nlo_as_01180",
            "--orders=a2,as1a2",
        ])
        .assert()
        .success()
        .stdout(LHCB_WP_7TEV_STR);

    let optimized = NamedTempFile::new("fktable1b.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--optimize-fk-table",
            "Nf4Sym",
            output.path().to_str().unwrap(),
            optimized.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            optimized.path().to_str().unwrap(),
            "NNPDF40_nlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(LHCB_WP_7TEV_OPTIMIZED_STR);
}

#[test]
fn lhcb_wp_7tev_v2() {
    let output = NamedTempFile::new("fktable2a.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "evolve",
            "--digits-abs=14",
            "--digits-rel=7",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV_v2.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nlo_as_01180",
            "--orders=a2,as1a2",
        ])
        .assert()
        .success()
        .stdout(LHCB_WP_7TEV_V2_STR);
}

#[test]
fn lhcb_wp_7tev_v2_xir_2() {
    let output = NamedTempFile::new("fktable2b.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "evolve",
            "--digits-abs=14",
            "--digits-rel=7",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV_v2.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nlo_as_01180",
            "--orders=a2,as1a2",
            "--xir=2",
        ])
        .assert()
        .success()
        .stdout(LHCB_WP_7TEV_V2_XIR2_STR);
}

#[test]
fn lhcb_wp_7tev_v2_xif_2() {
    let output = NamedTempFile::new("fktable2c.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "evolve",
            "--digits-abs=14",
            "--digits-rel=7",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV_v2_xif_2.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nlo_as_01180",
            "--orders=a2,as1a2",
            "--xif=2",
        ])
        .assert()
        .success()
        .stdout(LHCB_WP_7TEV_V2_XIF_2_STR);
}

#[test]
fn lhcb_wp_7tev_v2_xif_2_error() {
    let output = NamedTempFile::new("fktable2c.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "evolve",
            "--digits-abs=16",
            "--digits-rel=16",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV_v2.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nlo_as_01180",
            "--orders=a2,as1a2",
            "--xif=2",
        ])
        .assert()
        .failure()
        .stderr(LHCB_WP_7TEV_V2_XIF_2_ERROR_STR)
        .stdout("");
}

#[test]
fn e906nlo_bin_00() {
    let input = NamedTempFile::new("E906nlo_bin_00_unique_bin_limits.pineappl.lz4").unwrap();
    let output = NamedTempFile::new("fktable2.lz4").unwrap();

    // we need to delete bins with the same bin limits for `Grid::merge` to work properly
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--delete-bins=1-3",
            "../test-data/E906nlo_bin_00.pineappl.lz4",
            input.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "evolve",
            input.path().to_str().unwrap(),
            "../test-data/E906nlo_bin_00.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(E906NLO_BIN_00_STR);
}

#[test]
fn nutev_cc_nu_fe_sigmared() {
    let output = NamedTempFile::new("fktable3.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
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

#[test]
fn lhcb_dy_8tev() {
    let output = NamedTempFile::new("fktable4.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "evolve",
            "../test-data/LHCB_DY_8TEV.pineappl.lz4",
            "../test-data/LHCB_DY_8TEV.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nlo_as_01180",
            "--orders=a2,as1a2",
        ])
        .assert()
        .success()
        .stdout(LHCB_DY_8TEV_STR);
}

#[test]
fn cms_ttb_8tev_2d_ttm_trap_tot() {
    let output = NamedTempFile::new("fktable5.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "evolve",
            "--orders=as2,as3,as4",
            "--xir=2",
            "../test-data/CMS_TTB_8TEV_2D_TTM_TRAP_TOT-opt.pineappl.lz4",
            "../test-data/CMS_TTB_8TEV_2D_TTM_TRAP_TOT.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nnlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(CMS_TTB_8TEV_2D_TTM_TRAP_TOT_STR);
}

#[test]
fn star_wmwp_510gev_wm_al_pol() {
    let output = NamedTempFile::new("fktable6.lz4").unwrap();

    // Grid is (PolPDF, UnpolPDF) but EKOs are ordered as {UnpolEko, PolEko} to
    // check that order doesn't matter.
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "evolve",
            "../test-data/STAR_WMWP_510GEV_WM-AL-POL.pineappl.lz4",
            "../test-data/STAR_WMWP_510GEV_WM-AL-POL_UnpolPDF.tar,../test-data/STAR_WMWP_510GEV_WM-AL-POL_PolPDF.tar",
            output.path().to_str().unwrap(),
            "240608-tr-pol-nlo-100+p,NNPDF40_nlo_pch_as_01180",
        ])
        .assert()
        .success()
        .stdout(STAR_WMWP_510GEV_WM_AL_POL);
}

#[test]
fn zeus_2jet() {
    let output = NamedTempFile::new("fktable7.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "evolve",
            "--accuracy=3e-3",
            "../test-data/ZEUS_2JET_319GEV_374PB-1_DIF_ETQ2_BIN6.pineappl.lz4",
            "../test-data/ZEUS_2JET_319GEV_374PB-1_DIF_ETQ2_BIN6.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nnlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(ZEUS_2JET_STR);
}

#[test]
fn lhcb_wp_8tev() {
    let output = NamedTempFile::new("fktable8.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "evolve",
            "../test-data/LHCB_WP_8TEV.pineappl.lz4",
            "../test-data/LHCB_WP_8TEV.tar",
            output.path().to_str().unwrap(),
            "NNPDF40_nnlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(LHCB_WP_8TEV_STR);
}
