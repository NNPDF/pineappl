use super::helpers;
use anyhow::{bail, Result};
use pineappl::bin::BinRemapper;

pub fn subcommand_integrated(input: &str, output: &str) -> Result<()> {
    let mut grid = helpers::read_grid(input)?;

    if grid.merge_bins(0..grid.bin_info().bins()).is_err() {
        bail!("TODO");
    }
    grid.set_remapper(
        BinRemapper::new(vec![1.0], vec![(0.0, 1.0)]).unwrap_or_else(|_| unreachable!()),
    )?;

    let dimensions = grid.bin_info().dimensions();
    let key_values = grid.key_values_mut();
    for dim in 0..dimensions {
        key_values.remove(&format!("x{}_label", dim));
        key_values.remove(&format!("x{}_label_tex", dim));
        key_values.remove(&format!("x{}_unit", dim));
    }
    key_vales.remove("y_label");
    key_values.remove("y_label_tex");
    key_values.remove("y_unit");

    helpers::write_grid(output, &grid)
}
