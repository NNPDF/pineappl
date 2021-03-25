use super::helpers;
use anyhow::Result;
use pineappl::subgrid::SubgridEnum;

pub fn subcommand(input: &str) -> Result<()> {
    let grid = helpers::read_grid(input)?;

    for order in 0..grid.orders().len() {
        for bin in 0..grid.bin_info().bins() {
            for lumi in 0..grid.lumi().len() {
                println!(
                    "({}, {}, {}) = {}",
                    order,
                    bin,
                    lumi,
                    match grid.subgrid(order, bin, lumi) {
                        SubgridEnum::LagrangeSubgridV1(_) => "LagrangeSubgridV1",
                        SubgridEnum::NtupleSubgridV1(_) => "NtupleSubgridV1",
                        SubgridEnum::LagrangeSparseSubgridV1(_) => "LagrangeSparseSubgridV1",
                        SubgridEnum::LagrangeSubgridV2(_) => "LagrangeSubgridV2",
                        SubgridEnum::ReadOnlySparseSubgridV1(_) => "ReadOnlySparseSubgridV1",
                    }
                );
            }
        }
    }

    Ok(())
}
