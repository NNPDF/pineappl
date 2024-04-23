use assert_cmd::Command;

const HELP_STR: &str = "Perform various analyses with grids

Usage: pineappl analyze <COMMAND>

Commands:
  ckf  Compare K-factors with channel K factors (ckf)

Options:
  -h, --help  Print help
";

const CKF_HELP_STR: &str = "Compare K-factors with channel K factors (ckf)

Usage: pineappl analyze ckf [OPTIONS] <INPUT> <PDFSET> <ORDER> [ORDERS_DEN]...

Arguments:
  <INPUT>          Path to the input grid
  <PDFSET>         LHAPDF id or name of the PDF set
  <ORDER>          Order defining the K factors
  [ORDERS_DEN]...  Normalizing orders of the K factors

Options:
  -l, --limit <LIMIT>     The maximum number of channels displayed [default: 10]
      --digits-rel <REL>  Set the number of fractional digits shown for relative numbers [default: 2]
  -h, --help              Print help
";

const CKF_STR: &str = "b   etal    bin-K l  K   l  K   l  K   l  K   l  K  
     []                                             
-+----+----+-----+-+----+-+----+-+----+-+----+-+----
0    2 2.25  1.17 0 1.30 3 -inf 1 -inf 2 0.00 4 0.00
1 2.25  2.5  1.17 0 1.31 3 -inf 1 -inf 2 0.00 4 0.00
2  2.5 2.75  1.17 0 1.33 3 -inf 1 -inf 2 0.00 4 0.00
3 2.75    3  1.18 0 1.34 3 -inf 1 -inf 2 0.00 4 0.00
4    3 3.25  1.18 0 1.36 3 -inf 1 -inf 2 0.00 4 0.00
5 3.25  3.5  1.19 0 1.37 3 -inf 1 -inf 2 0.00 4 0.00
6  3.5    4  1.20 0 1.40 3 -inf 1 -inf 2 0.00 4 0.00
7    4  4.5  1.24 0 1.43 3 -inf 1 -inf 2 0.00 4 0.00
";

// TODO: understand these factors
const CKF_WITH_DEFAULT_DENOMINATOR_STR: &str =
    "b   etal    bin-K  l   K    l  K   l  K   l  K   l  K  
     []                                                
-+----+----+------+-+------+-+----+-+----+-+----+-+----
0    2 2.25 -13.20 0 -23.29 3 -inf 1 -inf 4 1.00 2 1.00
1 2.25  2.5 -14.37 0 -26.28 3 -inf 1 -inf 4 1.00 2 1.00
2  2.5 2.75 -16.17 0 -30.26 3 -inf 1 -inf 4 1.00 2 1.00
3 2.75    3 -18.67 0 -35.49 3 -inf 1 -inf 4 1.00 2 1.00
4    3 3.25 -20.81 0 -40.24 3 -inf 1 -inf 4 1.00 2 1.00
5 3.25  3.5 -22.55 0 -43.96 3 -inf 1 -inf 4 1.00 2 1.00
6  3.5    4 -23.46 0 -45.26 3 -inf 1 -inf 2 1.00 4 1.00
7    4  4.5 -23.31 0 -42.13 3 -inf 1 -inf 4 1.00 2 1.00
";

const CKF_WITH_BAD_LIMIT_STR: &str =
    "error: invalid value '0' for '--limit <LIMIT>': 0 is not in 1..=65535

For more information, try '--help'.
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["analyze", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn ckf_help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["analyze", "ckf", "--help"])
        .assert()
        .success()
        .stdout(CKF_HELP_STR);
}

#[test]
fn ckf() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "analyze",
            "ckf",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
            "a2as1",
            "a2",
        ])
        .assert()
        .success()
        .stdout(CKF_STR);
}

#[test]
fn ckf_with_default_denominator() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "analyze",
            "ckf",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
            "a2as1",
        ])
        .assert()
        .success()
        .stdout(CKF_WITH_DEFAULT_DENOMINATOR_STR);
}

#[test]
fn ckf_with_bad_limit() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "analyze",
            "ckf",
            "--limit=0",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
            "a2as1",
        ])
        .assert()
        .failure()
        .stderr(CKF_WITH_BAD_LIMIT_STR);
}
