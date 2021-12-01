use super::helpers;
use anyhow::{bail, Result};
use clap::{ArgGroup, Parser};
use pineappl::bin::BinRemapper;

/// Sums two or more bins of a grid together.
#[derive(Parser)]
#[clap(group = ArgGroup::new("mode").required(true), name = "sum")]
pub struct Opts {
    /// Path to the input grid.
    input: String,
    /// Path to the modified PineAPPL file.
    output: String,
    /// Sums all bins into a single bin.
    #[clap(long, group = "mode")]
    integrated: bool,
}

impl Opts {
    pub fn subcommand(&self) -> Result<()> {
        if self.integrated {
            let mut grid = helpers::read_grid(&self.input)?;

            if grid.merge_bins(0..grid.bin_info().bins()).is_err() {
                bail!("TODO");
            }
            grid.set_remapper(
                BinRemapper::new(vec![1.0], vec![(0.0, 1.0)]).unwrap_or_else(|_| unreachable!()),
            )?;

            let dimensions = grid.bin_info().dimensions();
            let key_values = grid.key_values_mut();
            for dim in 0..dimensions {
                key_values.remove(&format!("x{}_label", dim + 1));
                key_values.remove(&format!("x{}_label_tex", dim + 1));
                key_values.remove(&format!("x{}_unit", dim + 1));
            }
            key_values.remove("y_label");
            key_values.remove("y_label_tex");
            key_values.remove("y_unit");

            helpers::write_grid(&self.output, &grid)
        } else {
            unreachable!();
        }
    }
}
