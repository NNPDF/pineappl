use assert_cmd::Command;

const HELP_STR: &str = "pineappl-analyze 
Perform various analyses with grids

USAGE:
    pineappl analyze <SUBCOMMAND>

OPTIONS:
    -h, --help    Print help information

SUBCOMMANDS:
    ckf    Compare K-factors with channel K factors (ckf)
";

const CKF_HELP_STR: &str = "pineappl-analyze-ckf 
Compare K-factors with channel K factors (ckf)

USAGE:
    pineappl analyze ckf [OPTIONS] <INPUT> <PDFSET> <ORDER> [ORDERS_DEN]...

ARGS:
    <INPUT>            Path to the input grid
    <PDFSET>           LHAPDF id or name of the PDF set
    <ORDER>            Order defining the K factors
    <ORDERS_DEN>...    Normalizing orders of the K factors

OPTIONS:
        --digits-rel <REL>    Set the number of fractional digits shown for relative numbers
                              [default: 2]
        --force-positive      Forces negative PDF values to zero
    -h, --help                Print help information
    -l, --limit <LIMIT>       The maximum number of channels displayed [default: 10]
";

const CKF_STR: &str = "b   etal    bin-K l  K   l  K   l  K   l  K   l  K  
     []                                             
-+----+----+-----+-+----+-+----+-+----+-+----+-+----
0    2 2.25  1.17 0 1.30 3 -inf 1 -inf 2 0.00 4 0.00
1 2.25  2.5  1.17 0 1.31 3 -inf 1 -inf 2 0.00 4 0.00
2  2.5 2.75  1.18 0 1.33 3 -inf 1 -inf 2 0.00 4 0.00
3 2.75    3  1.18 0 1.34 3 -inf 1 -inf 2 0.00 4 0.00
4    3 3.25  1.18 0 1.36 3 -inf 1 -inf 2 0.00 4 0.00
5 3.25  3.5  1.19 0 1.38 3 -inf 1 -inf 2 0.00 4 0.00
6  3.5    4  1.20 0 1.40 3 -inf 1 -inf 2 0.00 4 0.00
7    4  4.5  1.24 0 1.44 3 -inf 1 -inf 2 0.00 4 0.00
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["analyze", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn ckf_help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["analyze", "ckf", "--help"])
        .assert()
        .success()
        .stdout(CKF_HELP_STR);
}

#[test]
fn ckf() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "analyze",
            "ckf",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
            "a2as1",
            "a2",
        ])
        .assert()
        .success()
        .stdout(CKF_STR);
}
