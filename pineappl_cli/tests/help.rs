use assert_cmd::Command;

const HELP_STR: &str = "Display a manpage for selected subcommands

Usage: pineappl help [SUBCOMMAND]...

Arguments:
  [SUBCOMMAND]...  Name of the (chain of) subcommand(s) to read the manpage of

Options:
  -h, --help  Print help
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["help", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}
