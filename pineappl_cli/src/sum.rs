use super::helpers::{self, Subcommand};
use anyhow::{bail, Result};
use clap::{ArgGroup, Parser, ValueHint};
use pineappl::bin::BinRemapper;
use std::path::PathBuf;

/// Sums two or more bins of a grid together.
#[derive(Parser)]
#[clap(group = ArgGroup::new("mode").required(true), name = "sum")]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the modified PineAPPL file.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// Sums all bins into a single bin.
    #[clap(long, group = "mode")]
    integrated: bool,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        if self.integrated {
            let mut grid = helpers::read_grid(&self.input)?;

            if grid.merge_bins(0..grid.bin_info().bins()).is_err() {
                bail!("TODO");
            }
            grid.set_remapper(
                BinRemapper::new(vec![1.0], vec![(0.0, 1.0)]).unwrap_or_else(|_| unreachable!()),
            )?;

            let dimensions = grid.bin_info().dimensions();
            let key_values = grid.key_values_mut();
            for dim in 0..dimensions {
                key_values.remove(&format!("x{}_label", dim + 1));
                key_values.remove(&format!("x{}_label_tex", dim + 1));
                key_values.remove(&format!("x{}_unit", dim + 1));
            }
            key_values.remove("y_label");
            key_values.remove("y_label_tex");
            key_values.remove("y_unit");

            helpers::write_grid(&self.output, &grid)
        } else {
            unreachable!();
        }
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;
    use assert_fs::NamedTempFile;

    const HELP_STR: &str = "pineappl-sum 

Sums two or more bins of a grid together

USAGE:
    pineappl sum <--integrated> <INPUT> <OUTPUT>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path to the modified PineAPPL file

OPTIONS:
    -h, --help          Print help information
        --integrated    Sums all bins into a single bin
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["sum", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn integrated() {
        let output = NamedTempFile::new("summed.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "sum",
                "--integrated",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                output.path().to_str().unwrap(),
            ])
            .assert()
            .success()
            .stdout("");
    }
}
