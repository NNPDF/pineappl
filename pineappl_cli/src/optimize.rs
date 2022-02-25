use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use pineappl::fk_table::{FkAssumptions, FkTable};
use std::convert::TryFrom;
use std::path::PathBuf;

/// Optimizes the internal data structure to minimize memory usage.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the optimized PineAPPL file.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
    #[clap(
        long = "fk-table",
        parse(try_from_str = FkAssumptions::try_from),
        possible_values = &[
            "Nf6Ind",
            "Nf6Sym",
            "Nf5Ind",
            "Nf5Sym",
            "Nf4Ind",
            "Nf4Sym",
            "Nf3Ind",
            "Nf3Sym",
        ],
        value_name = "ASSUMPTIONS",
    )]
    fk_table: Option<FkAssumptions>,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;

        if let Some(assumptions) = self.fk_table {
            let mut fk_table = FkTable::try_from(grid)?;
            fk_table.optimize(assumptions);
            helpers::write_grid(&self.output, fk_table.grid())
        } else {
            grid.optimize();
            helpers::write_grid(&self.output, &grid)
        }
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;
    use assert_fs::NamedTempFile;

    const HELP_STR: &str = "pineappl-optimize 
Optimizes the internal data structure to minimize memory usage

USAGE:
    pineappl optimize [OPTIONS] <INPUT> <OUTPUT>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path to the optimized PineAPPL file

OPTIONS:
        --fk-table <ASSUMPTIONS>    [possible values: Nf6Ind, Nf6Sym, Nf5Ind, Nf5Sym, Nf4Ind,
                                    Nf4Sym, Nf3Ind, Nf3Sym]
    -h, --help                      Print help information
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
