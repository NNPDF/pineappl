use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use std::path::PathBuf;

/// Converts the file format to the most recent version.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the upgraded PineAPPL file.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    output: PathBuf,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<u8> {
        let mut grid = helpers::read_grid(&self.input)?;
        grid.upgrade();
        helpers::write_grid(&self.output, &grid)
    }
}
