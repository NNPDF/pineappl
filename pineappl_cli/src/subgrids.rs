use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use pineappl::subgrid::{Subgrid, SubgridEnum};
use prettytable::{cell, row};
use std::path::PathBuf;

/// Print information about the internal subgrid types.
#[derive(Parser)]
#[clap(name = "subgrids")]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Show empty subgrids.
    #[clap(long = "show-empty")]
    show_empty: bool,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
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

        table.printstd();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-subgrids 

Print information about the internal subgrid types

USAGE:
    pineappl subgrids [OPTIONS] <INPUT>

ARGS:
    <INPUT>    Path to the input grid

OPTIONS:
    -h, --help          Print help information
        --show-empty    Show empty subgrids
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["subgrids", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }
}
