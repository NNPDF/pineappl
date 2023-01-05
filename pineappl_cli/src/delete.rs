use super::helpers::{self, GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::{ArgGroup, Parser, ValueHint};
use std::ops::RangeInclusive;
use std::path::PathBuf;
use std::process::ExitCode;

/// Deletes parts from a PineAPPL grid.
#[derive(Parser)]
#[command(group = ArgGroup::new("mode").multiple(true).required(true))]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the modified PineAPPL file.
    #[arg(value_hint = ValueHint::FilePath)]
    output: PathBuf,
    #[arg(
        group = "mode",
        long,
        num_args(1..),
        value_delimiter = ',',
        value_parser = helpers::parse_integer_range
    )]
    /// Indices of bins that should be deleted.
    bins: Vec<RangeInclusive<usize>>,
}

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<ExitCode> {
        let mut grid = helpers::read_grid(&self.input)?;
        let bins: Vec<usize> = self
            .bins
            .iter()
            .flat_map(|range| range.clone().collect::<Vec<_>>())
            .collect();

        grid.delete_bins(&bins);

        helpers::write_grid(&self.output, &grid)
    }
}
