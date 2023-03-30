use assert_cmd::Command;
use assert_fs::NamedTempFile;

const HELP_STR: &str = "Converts the file format to the most recent version

Usage: pineappl upgrade <INPUT> <OUTPUT>

Arguments:
  <INPUT>   Path to the input grid
  <OUTPUT>  Path to the upgraded PineAPPL file

Options:
  -h, --help  Print help
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["upgrade", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn upgrade() {
    let output = NamedTempFile::new("upgraded.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "upgrade",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");
}
