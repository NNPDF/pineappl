use assert_cmd::Command;
use assert_fs::{fixture::FileWriteStr, NamedTempFile};

const HELP_STR: &str = "pineappl-set 
Modifies the internal key-value storage

USAGE:
    pineappl set [OPTIONS] <INPUT> <OUTPUT>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path of the modified PineAPPL file

OPTIONS:
        --delete <KEY>
            Deletes an internal key-value pair

        --entry <KEY> <VALUE>
            Sets an internal key-value pair

        --entry-from-file <KEY> <FILE>
            Sets an internal key-value pair, with value being read from a file

    -h, --help
            Print help information
";

const DEFAULT_STR: &str = r#"arxiv: 1505.07024
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

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["set", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn default() {
    let output = NamedTempFile::new("set.pineappl.lz4").unwrap();
    let file = NamedTempFile::new("file").unwrap();

    file.write_str("one\ntwo\nthree\nfour").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "set",
            "--delete=runcard",
            "--entry",
            "key",
            "value",
            "--entry-from-file",
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
        .args(&["info", "--show", output.path().to_str().unwrap()])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}
