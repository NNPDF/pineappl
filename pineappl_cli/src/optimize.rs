use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use std::path::PathBuf;

/// Optimizes the internal data structure to minimize memory usage.
#[derive(Parser)]
#[clap(name = "optimize")]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the optimized PineAPPL file.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;
        grid.optimize();
        helpers::write_grid(&self.output, &grid)
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;
    use assert_fs::NamedTempFile;

    const HELP_STR: &str = "pineappl-optimize 
Optimizes the internal data structure to minimize memory usage

USAGE:
    pineappl optimize <INPUT> <OUTPUT>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path to the optimized PineAPPL file

OPTIONS:
    -h, --help    Print help information
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["optimize", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn default() {
        let output = NamedTempFile::new("optimized.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "optimize",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                output.path().to_str().unwrap(),
            ])
            .assert()
            .success()
            .stdout("");
    }
}
