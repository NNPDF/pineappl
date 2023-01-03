use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{ArgGroup, Parser, ValueHint};
use std::ops::RangeInclusive;
use std::path::PathBuf;

/// Deletes parts from a PineAPPL grid.
#[derive(Parser)]
#[clap(group = ArgGroup::new("mode").multiple(true).required(true))]
pub struct Opts {
    /// Path to the input grid.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the modified PineAPPL file.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    output: PathBuf,
    #[clap(
        group = "mode",
        long,
        multiple_values = true,
        use_value_delimiter = true,
        value_parser = helpers::parse_integer_range
    )]
    /// Indices of bins that should be deleted.
    bins: Vec<RangeInclusive<usize>>,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<u8> {
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
