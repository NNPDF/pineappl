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
