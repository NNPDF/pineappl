use super::helpers;
use anyhow::Result;
use pineappl::subgrid::SubgridEnum;
use prettytable::{cell, row, Table};

pub fn subcommand(input: &str) -> Result<Table> {
    let grid = helpers::read_grid(input)?;
    let mut table = helpers::create_table();

    let titles = row![c => "order", "bin", "lumi", "type"];
    table.set_titles(titles);

    for order in 0..grid.orders().len() {
        for bin in 0..grid.bin_info().bins() {
            for lumi in 0..grid.lumi().len() {
                let row = table.add_empty_row();
                row.add_cell(cell!(l->&format!("{}", order)));
                row.add_cell(cell!(l->&format!("{}", bin)));
                row.add_cell(cell!(l->&format!("{}", lumi)));
                row.add_cell(cell!(l->
                    match grid.subgrid(order, bin, lumi) {
                        SubgridEnum::LagrangeSubgridV1(_) => "LagrangeSubgridV1",
                        SubgridEnum::NtupleSubgridV1(_) => "NtupleSubgridV1",
                        SubgridEnum::LagrangeSparseSubgridV1(_) => "LagrangeSparseSubgridV1",
                        SubgridEnum::LagrangeSubgridV2(_) => "LagrangeSubgridV2",
                        SubgridEnum::ImportOnlySubgridV1(_) => "ImportOnlySubgridV1",
                        SubgridEnum::ImportOnlySubgridV2(_) => "ImportOnlySubgridV2",
                        SubgridEnum::EmptySubgridV1(_) => "EmptySubgridV1",
                    }
                ));
            }
        }
    }

    Ok(table)
}
