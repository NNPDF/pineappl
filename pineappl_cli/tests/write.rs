use assert_cmd::Command;
use assert_fs::{fixture::FileWriteStr, NamedTempFile};

const HELP_STR: &str = "Write a grid modified by various operations

Usage: pineappl write [OPTIONS] <INPUT> <OUTPUT>

Arguments:
  <INPUT>   Path to the input grid
  <OUTPUT>  Path of the modified PineAPPL file

Options:
      --cc1                           Charge conjugate the first initial state
      --cc2                           Charge conjugate the second initial state
      --delete-bins <BIN1-BIN2,...>   Delete bins with the specified indices
      --delete-key <KEY>              Delete an internal key-value pair
      --merge-bins <BIN1-BIN2>        Merge specific bins together
      --optimize                      Optimize internal data structure to minimize memory and disk usage
      --optimize-fk-table <OPTIMI>    Optimize internal data structure of an FkTable to minimize memory and disk usage [possible values: Nf6Ind, Nf6Sym, Nf5Ind, Nf5Sym, Nf4Ind, Nf4Sym, Nf3Ind, Nf3Sym]
      --remap <REMAPPING>             Modify the bin dimensions and widths
      --remap-norm <NORM>             Modify the bin normalizations with a common factor
      --remap-norm-ignore <DIMS>      Modify the bin normalizations by multiplying with the bin lengths for the given dimensions
  -s, --scale <SCALE>                 Scales all grids with the given factor
      --scale-by-bin <BIN1,BIN2,...>  Scale each bin with a different factor
      --scale-by-order <AS,AL,LR,LF>  Scales all grids with order-dependent factors
      --set-key-value <KEY> <VALUE>   Set an internal key-value pair
      --set-key-file <KEY> <FILE>     Set an internal key-value pair, with value being read from a file
      --upgrade                       Convert the file format to the most recent version
  -h, --help                          Print help
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

const DELETE_BINS_02_57_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0 2.75    3 2.4257663e2    -3.77     2.92
1    3 3.25 1.8093343e2    -3.74     2.95
";

const DELETE_BINS_25_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 3.7527620e2    -3.77     2.71
1 2.25  2.5 3.4521553e2    -3.79     2.80
2  3.5    4 5.7851018e1    -3.63     2.97
3    4  4.5 1.3772029e1    -3.46     2.85
";

const KEY_VALUE_STR: &str = r#"arxiv: 1505.07024
description: LHCb differential W-boson production cross section at 7 TeV
hepdata: 10.17182/hepdata.2114.v1/t4
initial_state_1: 2212
initial_state_2: 2212
key: value
lumi_id_types: pdg_mc_ids
mg5amc_repo: 
mg5amc_revno: 
multiline: one
two
three
four
nnpdf_id: LHCBWZMU7TEV
pineappl_gitversion: v0.4.1-36-gdbdb5d0
results: ----------------------------------------------------------------------
   PineAPPL         MC        sigma      central         min      max 
                              1/100   sigma   1/1000   1/1000   1/1000
----------------------------------------------------------------------
 1.876381e+02  1.876313e+02   0.082   0.044   0.0360   0.0396   0.0317
 1.726078e+02  1.726041e+02   0.082   0.026   0.0211   0.0246   0.0166
 1.500070e+02  1.500056e+02   0.051   0.019   0.0095   0.0103   0.0075
 1.212883e+02  1.212890e+02   0.047   0.012   0.0056   0.0052   0.0068
 9.046672e+01  9.046795e+01   0.057   0.024   0.0136   0.0127   0.0146
 6.145558e+01  6.145650e+01   0.063   0.024   0.0151   0.0116   0.0193
 5.785102e+01  5.784687e+01   0.075   0.095   0.0717   0.0874   0.0518
 1.377203e+01  1.376219e+01   0.119   0.599   0.7153   0.7646   0.6465

runcard_gitversion: 7b42083
x1_label: etal
x1_label_tex: $\eta_{\bar{\ell}}$
x1_unit: 
y_label: disg/detal
y_label_tex: $\frac{\mathrm{d}\sigma}{\mathrm{d}\eta_{\bar{\ell}}}$
y_unit: pb
"#;

const MERGE_BINS_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 3.7527620e2    -3.77     2.71
1 2.25  2.5 3.4521553e2    -3.79     2.80
2  2.5 2.75 3.0001406e2    -3.78     2.86
3 2.75    3 2.4257663e2    -3.77     2.92
4    3 3.25 1.8093343e2    -3.74     2.95
5 3.25  3.5 1.2291115e2    -3.71     2.98
6  3.5  4.5 3.5811524e1    -3.59     2.95
";

const REMAP_STR: &str = "b etal  x2  x3  disg/detal  scale uncertainty
   []   []  []     [pb]            [%]       
-+--+--+-+-+-+-+-----------+--------+--------
0  0  1 0 2 1 2 1.8763810e1    -3.77     2.71
1  0  1 0 2 2 3 1.7260776e1    -3.79     2.80
2  0  1 0 2 3 4 1.5000703e1    -3.78     2.86
3  0  1 0 2 4 5 1.2128831e1    -3.77     2.92
4  0  1 2 4 1 2 9.0466716e0    -3.74     2.95
5  1  2 0 2 8 9 6.1455576e0    -3.71     2.98
6  1  2 2 4 3 4 5.7851018e0    -3.63     2.97
7  1  2 2 4 4 5 1.3772029e0    -3.46     2.85
";

const SCALE_BY_BIN_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 3.7527620e2    -3.77     2.71
1 2.25  2.5 6.9043106e2    -3.79     2.80
2  2.5 2.75 9.0004217e2    -3.78     2.86
3 2.75    3 9.7030651e2    -3.77     2.92
4    3 3.25 9.0466716e2    -3.74     2.95
5 3.25  3.5 7.3746691e2    -3.71     2.98
6  3.5    4 4.0495712e2    -3.63     2.97
7    4  4.5 1.1017623e2    -3.46     2.85
";

const SCALE_BY_ORDER_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 2.1481594e2    -4.53     3.53
1 2.25  2.5 1.9807977e2    -4.50     3.55
2  2.5 2.75 1.7256275e2    -4.43     3.54
3 2.75    3 1.3976747e2    -4.36     3.52
4    3 3.25 1.0460103e2    -4.26     3.49
5 3.25  3.5 7.1393138e1    -4.17     3.45
6  3.5    4 3.3881527e1    -4.01     3.50
7    4  4.5 8.2340900e0    -3.70     3.92
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
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["info", "--show", output.path().to_str().unwrap()])
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
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
fn scale_by_bin() {
    let output = NamedTempFile::new("scale_by_bin.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--scale-by-bin=1,2,3,4,5,6,7,8",
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
fn upgrade() {
    let output = NamedTempFile::new("upgraded.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "write",
            "--upgrade",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");
}
