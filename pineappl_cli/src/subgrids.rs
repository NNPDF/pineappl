use super::helpers;
use anyhow::Result;
use clap::Parser;
use pineappl::subgrid::{Subgrid, SubgridEnum};
use prettytable::{cell, row, Table};

/// Print information about the internal subgrid types.
#[derive(Parser)]
#[clap(name = "subgrids")]
pub struct Opts {
    /// Path to the input grid.
    input: String,
    /// Show empty subgrids.
    #[clap(long = "show-empty")]
    show_empty: bool,
}

impl Opts {
    pub fn subcommand(&self) -> Result<Table> {
        let grid = helpers::read_grid(&self.input)?;
        let mut table = helpers::create_table();

        let titles = row![c => "order", "bin", "lumi", "type"];
        table.set_titles(titles);

        for ((order, bin, lumi), subgrid) in grid.subgrids().indexed_iter() {
            if !self.show_empty && subgrid.is_empty() {
                continue;
            }

            let row = table.add_empty_row();
            row.add_cell(cell!(l->&format!("{}", order)));
            row.add_cell(cell!(l->&format!("{}", bin)));
            row.add_cell(cell!(l->&format!("{}", lumi)));
            row.add_cell(cell!(l->
                match subgrid {
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

        Ok(table)
    }
}
