use pineappl::grid::Grid;
use pineappl::subgrid::SubgridEnum;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;

pub fn subcommand(input: &str) -> Result<(), Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;

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
