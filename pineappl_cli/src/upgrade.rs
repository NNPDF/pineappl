use super::helpers;
use anyhow::Result;
use clap::Parser;

/// Converts the file format to the most recent version.
#[derive(Parser)]
#[clap(name = "upgrade")]
pub struct Opts {
    /// Path to the input grid.
    input: String,
    /// Path to the upgraded PineAPPL file.
    output: String,
}

impl Opts {
    pub fn subcommand(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;
        grid.upgrade();
        helpers::write_grid(&self.output, &grid)
    }
}
