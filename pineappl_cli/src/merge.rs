use super::helpers;
use anyhow::Result;

pub fn subcommand(
    output: &str,
    input0: &str,
    input_rest: &[&str],
    scale: Option<f64>,
    scale_by_order: &[f64],
) -> Result<()> {
    let mut grid0 = helpers::read_grid(input0)?;

    for i in input_rest {
        grid0.merge(helpers::read_grid(i)?)?;
    }

    if let Some(scale) = scale {
        grid0.scale(scale);
    } else if !scale_by_order.is_empty() {
        grid0.scale_by_order(
            scale_by_order[0],
            scale_by_order[1],
            scale_by_order[2],
            scale_by_order[3],
            scale_by_order[4],
        );
    }

    helpers::write_grid(output, &grid0)
}
