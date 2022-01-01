use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use std::fs;
use std::path::PathBuf;

/// Modifies the internal key-value storage.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path of the modified PineAPPL file.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// Deletes an internal key-value pair.
    #[clap(long, multiple_occurrences = true, value_name = "KEY")]
    delete: Vec<String>,
    /// Sets an internal key-value pair.
    #[clap(
        long,
        multiple_occurrences = true,
        number_of_values = 2,
        value_names = &["KEY", "VALUE"]
    )]
    entry: Vec<String>,
    /// Sets an internal key-value pair, with value being read from a file.
    #[clap(
        alias = "entry_from_file",
        long = "entry-from-file",
        multiple_occurrences = true,
        number_of_values = 2,
        value_names = &["KEY", "FILE"]
    )]
    entry_from_file: Vec<String>,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;

        for key_value in self.entry.chunks(2) {
            grid.set_key_value(&key_value[0], &key_value[1]);
        }

        for key_file in self.entry_from_file.chunks(2) {
            grid.set_key_value(&key_file[0], &fs::read_to_string(&key_file[1])?);
        }

        for delete in &self.delete {
            grid.key_values_mut().remove(delete);
        }

        helpers::write_grid(&self.output, &grid)
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

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

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["set", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }
}
