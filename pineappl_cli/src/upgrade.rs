use super::helpers;
use anyhow::Result;
use clap::{Parser, ValueHint};
use std::path::PathBuf;

/// Converts the file format to the most recent version.
#[derive(Parser)]
#[clap(name = "upgrade")]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the upgraded PineAPPL file.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
}

impl Opts {
    pub fn subcommand(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;
        grid.upgrade();
        helpers::write_grid(&self.output, &grid)
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-upgrade 

Converts the file format to the most recent version

USAGE:
    pineappl upgrade <INPUT> <OUTPUT>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path to the upgraded PineAPPL file

OPTIONS:
    -h, --help    Print help information
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
}
