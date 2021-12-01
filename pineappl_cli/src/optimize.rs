use super::helpers;
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

impl Opts {
    pub fn subcommand(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;
        grid.optimize();
        helpers::write_grid(&self.output, &grid)
    }
}
