use assert_cmd::Command;
use assert_fs::{fixture::FileWriteStr, NamedTempFile};

const HELP_STR: &str = "Write a grid modified by various operations

Usage: pineappl write [OPTIONS] <INPUT> <OUTPUT>

Arguments:
  <INPUT>   Path to the input grid
  <OUTPUT>  Path of the modified PineAPPL file

Options:
      --cc1[=<ENABLE>]                Charge conjugate the first initial state [possible values: true, false]
      --cc2[=<ENABLE>]                Charge conjugate the second initial state [possible values: true, false]
      --dedup-channels[=<ULPS>]       Deduplicate channels assuming numbers differing by ULPS are the same
      --delete-bins <BIN1-BIN2,...>   Delete bins with the specified indices
      --delete-key <KEY>              Delete an internal key-value pair
      --merge-bins <BIN1-BIN2,...>    Merge specific bins together
      --optimize[=<ENABLE>]           Optimize internal data structure to minimize memory and disk usage [possible values: true, false]
      --optimize-fk-table <OPTIMI>    Optimize internal data structure of an FkTable to minimize memory and disk usage [possible values: Nf6Ind, Nf6Sym, Nf5Ind, Nf5Sym, Nf4Ind, Nf4Sym, Nf3Ind, Nf3Sym]
      --remap <REMAPPING>             Modify the bin dimensions and widths
      --remap-norm <NORM>             Modify the bin normalizations with a common factor
      --remap-norm-ignore <DIM1,...>  Modify the bin normalizations by multiplying with the bin lengths for the given dimensions
      --rewrite-channel <IDX> <CHAN>  Rewrite the definition of the channel with index IDX
  -s, --scale <SCALE>                 Scales all grids with the given factor
      --scale-by-bin <BIN1,BIN2,...>  Scale each bin with a different factor
      --scale-by-order <AS,AL,LR,LF>  Scales all grids with order-dependent factors
      --set-key-value <KEY> <VALUE>   Set an internal key-value pair
      --set-key-file <KEY> <FILE>     Set an internal key-value pair, with value being read from a file
      --split-lumi[=<ENABLE>]         Split the grid such that the luminosity function contains only a single combination per channel [possible values: true, false]
      --upgrade[=<ENABLE>]            Convert the file format to the most recent version [possible values: true, false]
  -h, --help                          Print help
";

const CHANNEL_STR: &str = "l    entry        entry
-+------------+------------
0 1 × ( 2, -1) 1 × ( 4, -3)
1 1 × (21, -3) 1 × (21, -1)
2 1 × (22, -3) 1 × (22, -1)
3 1 × ( 2, 21) 1 × ( 4, 21)
4 1 × ( 2, 22) 1 × ( 4, 22)
";

const DEDUP_CHANNEL_DIFF_STR: &str = "b    x1               O(as^0 a^2)                       O(as^0 a^3)                       O(as^1 a^2)          
-+----+----+-----------+-----------+-------+-------------+-------------+-------+-----------+-----------+-------
0    2 2.25 6.5070305e2 6.5070305e2 0.000e0  -7.8692484e0  -7.8692484e0 0.000e0 1.1175729e2 1.1175729e2 0.000e0
1 2.25  2.5 5.9601236e2 5.9601236e2 0.000e0  -6.5623495e0  -6.5623495e0 0.000e0 1.0083341e2 1.0083341e2 0.000e0
2  2.5 2.75 5.1561247e2 5.1561247e2 0.000e0  -5.2348261e0  -5.2348261e0 0.000e0 8.9874343e1 8.9874343e1 0.000e0
3 2.75    3 4.1534629e2 4.1534629e2 0.000e0  -3.7590420e0  -3.7590420e0 0.000e0 7.3935106e1 7.3935106e1 0.000e0
4    3 3.25 3.0812719e2 3.0812719e2 0.000e0  -2.5871885e0  -2.5871885e0 0.000e0 5.6414554e1 5.6414554e1 0.000e0
5 3.25  3.5 2.0807482e2 2.0807482e2 0.000e0  -1.6762487e0  -1.6762487e0 0.000e0 3.9468336e1 3.9468336e1 0.000e0
6  3.5    4 9.6856769e1 9.6856769e1 0.000e0 -8.1027456e-1 -8.1027456e-1 0.000e0 1.9822014e1 1.9822014e1 0.000e0
7    4  4.5 2.2383492e1 2.2383492e1 0.000e0 -2.2022770e-1 -2.2022770e-1 0.000e0 5.3540011e0 5.3540011e0 0.000e0
";

const DEFAULT_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5459110e2
1 2.25  2.5 6.9028342e2
2  2.5 2.75 6.0025198e2
3 2.75    3 4.8552235e2
4    3 3.25 3.6195456e2
5 3.25  3.5 2.4586691e2
6  3.5    4 1.1586851e2
7    4  4.5 2.7517266e1
";

const DELETE_BINS_02_57_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0 2.75    3 4.8552235e2
1    3 3.25 3.6195456e2
";

const DELETE_BINS_25_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5459110e2
1 2.25  2.5 6.9028342e2
2  3.5    4 1.1586851e2
3    4  4.5 2.7517266e1
";

const KEY_VALUE_STR: &str = r"arxiv: 1505.07024
description: LHCb differential W-boson production cross section at 7 TeV
hepdata: 10.17182/hepdata.2114.v1/t4
initial_state_1: 2212
initial_state_2: 2212
key: value
lumi_id_types: pdg_mc_ids
mg5amc_repo: http://bazaar.launchpad.net/~maddevelopers/mg5amcnlo/3.1.2/
mg5amc_revno: 983
multiline: one
two
three
four
nnpdf_id: LHCBWZMU7TEV
pineappl_gitversion: v0.4.1-114-gdce19e0
results: ----------------------------------------------------------------------
   PineAPPL         MC        sigma      central         min      max 
                              1/100   sigma   1/1000   1/1000   1/1000
----------------------------------------------------------------------
 3.772955e+02  3.772821e+02   0.165   0.022   0.0357   0.0392   0.0313
 3.451417e+02  3.451342e+02   0.179   0.012   0.0217   0.0251   0.0172
 3.001260e+02  3.001231e+02   0.029   0.033   0.0096   0.0104   0.0076
 2.427612e+02  2.427624e+02   0.024   0.021   0.0049   0.0046   0.0060
 1.809773e+02  1.809799e+02   0.023   0.062   0.0143   0.0134   0.0154
 1.229334e+02  1.229354e+02   0.028   0.056   0.0157   0.0120   0.0200
 1.158685e+02  1.158603e+02   0.029   0.245   0.0708   0.0859   0.0514
 2.751727e+01  2.749798e+01   0.074   0.944   0.7014   0.7554   0.6281

runcard_gitversion: 82de4ad
x1_label: etal
x1_label_tex: $\eta_{\bar{\ell}}$
x1_unit: 
y_label: dsig/detal
y_label_tex: $\frac{\mathrm{d}\sigma}{\mathrm{d}\eta_{\bar{\ell}}}$
y_unit: pb
";

const MERGE_BINS_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5459110e2
1 2.25  2.5 6.9028342e2
2  2.5 2.75 6.0025198e2
3 2.75    3 4.8552235e2
4    3 3.25 3.6195456e2
5 3.25  3.5 2.4586691e2
6  3.5  4.5 7.1692887e1
";

const REMAP_STR: &str = "b etal  x2  x3  dsig/detal 
   []   []  []     [pb]    
-+--+--+-+-+-+-+-----------
0  0  1 0 2 1 2 3.7729555e1
1  0  1 0 2 2 3 3.4514171e1
2  0  1 0 2 3 4 3.0012599e1
3  0  1 0 2 4 5 2.4276118e1
4  0  1 2 4 1 2 1.8097728e1
5  1  2 0 2 8 9 1.2293345e1
6  1  2 2 4 3 4 1.1586851e1
7  1  2 2 4 4 5 2.7517266e0
";

const REMAP_NO_REMAPPER_STR: &str = "Error: grid does not have a remapper
";

const REWRITE_CHANNELS_CONVOLUTE_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5534392e2
1 2.25  2.5 6.9342538e2
2  2.5 2.75 6.0526279e2
3 2.75    3 4.9140786e2
4    3 3.25 3.6782869e2
5 3.25  3.5 2.5085041e2
6  3.5    4 1.1874486e2
7    4  4.5 2.8214633e1
";

const REWRITE_CHANNELS_LUMIS_STR: &str = "l              entry                            entry                       entry                       entry                        entry                  entry
-+--------------------------------+-------------------------------+-----------------------+--------------------------------+-----------------------+---------------------
0 0.0000128881 × ( 2, -5)          0.050940490000000005 × ( 2, -3) 0.9490461561 × ( 2, -1) 0.0017222500000000003 × ( 4, -5) 0.9473907556 × ( 4, -3) 0.05089536 × ( 4, -1)
1 0.0017351381000000003 × (-5, 21) 0.9983312456 × (-3, 21)         0.9999415161 × (-1, 21)                                                          
2 1 × (22, -3)                     1 × (22, -1)                                                                                                     
3 0.9999995342 × ( 2, 21)          1.0000083656 × ( 4, 21)                                                                                          
4 1 × ( 2, 22)                     1 × ( 4, 22)                                                                                                     
";

const SCALE_BY_BIN_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5459110e2
1 2.25  2.5 1.3805668e3
2  2.5 2.75 1.8007559e3
3 2.75    3 1.9420894e3
4    3 3.25 1.8097728e3
5 3.25  3.5 1.4752015e3
6  3.5    4 8.1107956e2
7    4  4.5 2.2013813e2
";

const SCALE_BY_ORDER_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 4.3317419e2
1 2.25  2.5 3.9555841e2
2  2.5 2.75 3.4506316e2
3 2.75    3 2.7972873e2
4    3 3.25 2.0918456e2
5 3.25  3.5 1.4266762e2
6  3.5    4 6.7845261e1
7    4  4.5 1.6435633e1
";

const SPLIT_LUMI_STR: &str = "l    entry
-+------------
0 1 × ( 2, -1)
1 1 × ( 4, -3)
2 1 × (21, -3)
3 1 × (21, -1)
4 1 × (22, -3)
5 1 × (22, -1)
6 1 × ( 2, 21)
7 1 × ( 4, 21)
8 1 × ( 2, 22)
9 1 × ( 4, 22)
";

const MULTIPLE_ARGUMENTS_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2  2.5 7.5454524e2
1  2.5    3 5.6456123e2
2    3 3.25 3.7400170e2
3 3.25  3.5 2.5300890e2
4  3.5    4 1.1909464e2
5    4  4.5 2.9004607e1
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["write", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn cc1() {
    let output = NamedTempFile::new("cc1.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--cc1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .stdout(DEFAULT_STR);
}

#[test]
fn cc2() {
    let output = NamedTempFile::new("cc2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--cc2",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .stdout(DEFAULT_STR);
}

#[test]
fn delete_bins_02_57() {
    let output = NamedTempFile::new("deleted.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--delete-bins=0-2,5-7",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .stdout(DELETE_BINS_02_57_STR);
}

#[test]
fn delete_bins_25() {
    let output = NamedTempFile::new("deleted2.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--delete-bins=2-5",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .stdout(DELETE_BINS_25_STR);
}

#[test]
fn key_value() {
    let output = NamedTempFile::new("set.pineappl.lz4").unwrap();
    let file = NamedTempFile::new("file").unwrap();

    file.write_str("one\ntwo\nthree\nfour").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--delete-key=runcard",
            "--set-key-value",
            "key",
            "value",
            "--set-key-file",
            "multiline",
            file.path().to_str().unwrap(),
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["read", "--show", output.path().to_str().unwrap()])
        .assert()
        .success()
        .stdout(KEY_VALUE_STR);
}

#[test]
fn merge_bins() {
    let output = NamedTempFile::new("bins.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--merge-bins=6-7",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .stdout(MERGE_BINS_STR);
}

#[test]
fn optimize() {
    // use `.pineappl` extension without `.lz4` to test `Grid::write` without compresssion
    let output = NamedTempFile::new("optimized.pineappl").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--optimize",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");
}

#[test]
fn remap() {
    let output = NamedTempFile::new("remapped.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--remap=0,1,2;0,2,4;1,2,3,4,5|:3|5:1,2,3,4,5,8,9|2:2",
            "--remap-norm-ignore=1",
            "--remap-norm=5",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .stdout(REMAP_STR);
}

#[test]
fn remap_norm_no_remapper() {
    let output = NamedTempFile::new("remapped.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--remap-norm=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(REMAP_NO_REMAPPER_STR);
}

#[test]
fn remap_norm_ignore_no_remapper() {
    let output = NamedTempFile::new("remapped.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--remap-norm-ignore=0",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(REMAP_NO_REMAPPER_STR);
}

#[test]
fn scale_by_bin() {
    let output = NamedTempFile::new("scale_by_bin.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--scale-by-bin=1,2,3,4,5,6,7,8",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .stdout(SCALE_BY_BIN_STR);
}

#[test]
fn scale_by_order() {
    let output = NamedTempFile::new("merged.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--scale-by-order=2,1,0.5,0.5",
            "--scale=0.5",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .stdout(SCALE_BY_ORDER_STR);
}

#[test]
fn split_lumi() {
    let output = NamedTempFile::new("split-lumi.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--split-lumi",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .stdout(DEFAULT_STR);

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--silence-lhapdf",
            "read",
            "--lumis",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(SPLIT_LUMI_STR);
}

#[test]
fn dedup_channels() {
    let output = NamedTempFile::new("dedup-channels.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--split-lumi",
            "--dedup-channels",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--silence-lhapdf",
            "diff",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEDUP_CHANNEL_DIFF_STR);

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--silence-lhapdf",
            "read",
            "--lumis",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(CHANNEL_STR);
}

#[test]
fn upgrade() {
    let output = NamedTempFile::new("upgraded.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--upgrade",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");
}

#[test]
fn multiple_arguments() {
    let output = NamedTempFile::new("multiple.merge.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "--merge-bins=0-1,2-3",
            "--scale=2",
            "--merge-bins=0-0",
            "--scale=0.5",
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
            "NNPDF40_nnlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(MULTIPLE_ARGUMENTS_STR);
}

#[test]
fn rewrite_channels() {
    let output = NamedTempFile::new("ckm_channels.pineappl.lz4").unwrap();

    // 0 1 × ( 2, -1) 1 × ( 4, -3)
    // 1 1 × (21, -3) 1 × (21, -1)
    // 2 1 × (22, -3) 1 × (22, -1)
    // 3 1 × ( 2, 21) 1 × ( 4, 21)
    // 4 1 × ( 2, 22) 1 × ( 4, 22)

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--rewrite-channel", "0", "0.9490461561 * ( 2, -1) + 0.050940490000000005 * (2, -3) + 0.0000128881 * (2, -5) + 0.05089536 * (4, -1) + 0.9473907556 * (4, -3) + 0.0017222500000000003 * (4, -5)",
            "--rewrite-channel", "1", "0.9999415161 * (-1, 21) + 0.9983312456 * (-3, 21) + 0.0017351381000000003 * (-5, 21)",
            "--rewrite-channel", "3", "0.9999995342 * ( 2, 21) + 1.0000083656 * ( 4, 21)",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["read", "--lumis", output.path().to_str().unwrap()])
        .assert()
        .success()
        .stdout(REWRITE_CHANNELS_LUMIS_STR);

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
        .stdout(REWRITE_CHANNELS_CONVOLUTE_STR);
}
