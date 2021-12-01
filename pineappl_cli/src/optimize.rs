use super::helpers;
use anyhow::Result;
use clap::Parser;

/// Optimizes the internal data structure to minimize memory usage.
#[derive(Parser)]
#[clap(name = "optimize")]
pub struct Opts {
    /// Path to the input grid.
    input: String,
    /// Path to the optimized PineAPPL file.
    output: String,
}

impl Opts {
    pub fn subcommand(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;
        grid.optimize();
        helpers::write_grid(&self.output, &grid)
    }
}
