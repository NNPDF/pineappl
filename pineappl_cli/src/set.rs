use super::helpers;
use anyhow::Result;
use clap::Parser;
use std::fs;

/// Modifies the internal key-value storage.
#[derive(Parser)]
#[clap(name = "set")]
pub struct Opts {
    /// Path to the input grid.
    input: String,
    /// Path of the modified PineAPPL file.
    output: String,
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

impl Opts {
    pub fn subcommand(&self) -> Result<()> {
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
