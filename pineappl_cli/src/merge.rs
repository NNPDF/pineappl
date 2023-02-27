use super::helpers::{self, GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use std::path::PathBuf;
use std::process::ExitCode;

/// Merges one or more PineAPPL grids together.
#[derive(Parser)]
pub struct Opts {
    /// Path of the merged PineAPPL file.
    #[arg(value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// Path(s) of the files that should be merged.
    #[arg(required = true, value_hint = ValueHint::FilePath)]
    input: Vec<PathBuf>,
}

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<ExitCode> {
        let (input0, input_rest) = self.input.split_first().unwrap();
        let mut grid0 = helpers::read_grid(input0)?;

        for i in input_rest {
            grid0.merge(helpers::read_grid(i)?)?;
        }

        helpers::write_grid(&self.output, &grid0)
    }
}
