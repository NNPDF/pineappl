use super::helpers;
use anyhow::Result;

pub fn subcommand(input: &str, output: &str) -> Result<()> {
    let mut grid = helpers::read_grid(input)?;
    grid.upgrade();
    helpers::write_grid(output, &grid)
}
