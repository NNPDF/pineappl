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

    const DEFAULT_STR: &str = "order bin lumi        type
-----+---+----+-------------------
0     0   0    ImportOnlySubgridV1
0     1   0    ImportOnlySubgridV1
0     2   0    ImportOnlySubgridV1
0     3   0    ImportOnlySubgridV1
0     4   0    ImportOnlySubgridV1
0     5   0    ImportOnlySubgridV1
0     6   0    ImportOnlySubgridV1
0     7   0    ImportOnlySubgridV1
1     0   0    ImportOnlySubgridV1
1     0   1    ImportOnlySubgridV1
1     0   3    ImportOnlySubgridV1
1     1   0    ImportOnlySubgridV1
1     1   1    ImportOnlySubgridV1
1     1   3    ImportOnlySubgridV1
1     2   0    ImportOnlySubgridV1
1     2   1    ImportOnlySubgridV1
1     2   3    ImportOnlySubgridV1
1     3   0    ImportOnlySubgridV1
1     3   1    ImportOnlySubgridV1
1     3   3    ImportOnlySubgridV1
1     4   0    ImportOnlySubgridV1
1     4   1    ImportOnlySubgridV1
1     4   3    ImportOnlySubgridV1
1     5   0    ImportOnlySubgridV1
1     5   1    ImportOnlySubgridV1
1     5   3    ImportOnlySubgridV1
1     6   0    ImportOnlySubgridV1
1     6   1    ImportOnlySubgridV1
1     6   3    ImportOnlySubgridV1
1     7   0    ImportOnlySubgridV1
1     7   1    ImportOnlySubgridV1
1     7   3    ImportOnlySubgridV1
3     0   0    ImportOnlySubgridV1
3     0   1    ImportOnlySubgridV1
3     0   3    ImportOnlySubgridV1
3     1   0    ImportOnlySubgridV1
3     1   1    ImportOnlySubgridV1
3     1   3    ImportOnlySubgridV1
3     2   0    ImportOnlySubgridV1
3     2   1    ImportOnlySubgridV1
3     2   3    ImportOnlySubgridV1
3     3   0    ImportOnlySubgridV1
3     3   1    ImportOnlySubgridV1
3     3   3    ImportOnlySubgridV1
3     4   0    ImportOnlySubgridV1
3     4   1    ImportOnlySubgridV1
3     4   3    ImportOnlySubgridV1
3     5   0    ImportOnlySubgridV1
3     5   1    ImportOnlySubgridV1
3     5   3    ImportOnlySubgridV1
3     6   0    ImportOnlySubgridV1
3     6   1    ImportOnlySubgridV1
3     6   3    ImportOnlySubgridV1
3     7   0    ImportOnlySubgridV1
3     7   1    ImportOnlySubgridV1
3     7   3    ImportOnlySubgridV1
4     0   0    ImportOnlySubgridV1
4     0   2    ImportOnlySubgridV1
4     0   4    ImportOnlySubgridV1
4     1   0    ImportOnlySubgridV1
4     1   2    ImportOnlySubgridV1
4     1   4    ImportOnlySubgridV1
4     2   0    ImportOnlySubgridV1
4     2   2    ImportOnlySubgridV1
4     2   4    ImportOnlySubgridV1
4     3   0    ImportOnlySubgridV1
4     3   2    ImportOnlySubgridV1
4     3   4    ImportOnlySubgridV1
4     4   0    ImportOnlySubgridV1
4     4   2    ImportOnlySubgridV1
4     4   4    ImportOnlySubgridV1
4     5   0    ImportOnlySubgridV1
4     5   2    ImportOnlySubgridV1
4     5   4    ImportOnlySubgridV1
4     6   0    ImportOnlySubgridV1
4     6   2    ImportOnlySubgridV1
4     6   4    ImportOnlySubgridV1
4     7   0    ImportOnlySubgridV1
4     7   2    ImportOnlySubgridV1
4     7   4    ImportOnlySubgridV1
6     0   0    ImportOnlySubgridV1
6     0   2    ImportOnlySubgridV1
6     0   4    ImportOnlySubgridV1
6     1   0    ImportOnlySubgridV1
6     1   2    ImportOnlySubgridV1
6     1   4    ImportOnlySubgridV1
6     2   0    ImportOnlySubgridV1
6     2   2    ImportOnlySubgridV1
6     2   4    ImportOnlySubgridV1
6     3   0    ImportOnlySubgridV1
6     3   2    ImportOnlySubgridV1
6     3   4    ImportOnlySubgridV1
6     4   0    ImportOnlySubgridV1
6     4   2    ImportOnlySubgridV1
6     4   4    ImportOnlySubgridV1
6     5   0    ImportOnlySubgridV1
6     5   2    ImportOnlySubgridV1
6     5   4    ImportOnlySubgridV1
6     6   0    ImportOnlySubgridV1
6     6   2    ImportOnlySubgridV1
6     6   4    ImportOnlySubgridV1
6     7   0    ImportOnlySubgridV1
6     7   2    ImportOnlySubgridV1
6     7   4    ImportOnlySubgridV1
";

    const SHOW_EMPTY_STR: &str = "order bin lumi        type
-----+---+----+-------------------
0     0   0    ImportOnlySubgridV1
0     0   1    EmptySubgridV1
0     0   2    EmptySubgridV1
0     0   3    EmptySubgridV1
0     0   4    EmptySubgridV1
0     1   0    ImportOnlySubgridV1
0     1   1    EmptySubgridV1
0     1   2    EmptySubgridV1
0     1   3    EmptySubgridV1
0     1   4    EmptySubgridV1
0     2   0    ImportOnlySubgridV1
0     2   1    EmptySubgridV1
0     2   2    EmptySubgridV1
0     2   3    EmptySubgridV1
0     2   4    EmptySubgridV1
0     3   0    ImportOnlySubgridV1
0     3   1    EmptySubgridV1
0     3   2    EmptySubgridV1
0     3   3    EmptySubgridV1
0     3   4    EmptySubgridV1
0     4   0    ImportOnlySubgridV1
0     4   1    EmptySubgridV1
0     4   2    EmptySubgridV1
0     4   3    EmptySubgridV1
0     4   4    EmptySubgridV1
0     5   0    ImportOnlySubgridV1
0     5   1    EmptySubgridV1
0     5   2    EmptySubgridV1
0     5   3    EmptySubgridV1
0     5   4    EmptySubgridV1
0     6   0    ImportOnlySubgridV1
0     6   1    EmptySubgridV1
0     6   2    EmptySubgridV1
0     6   3    EmptySubgridV1
0     6   4    EmptySubgridV1
0     7   0    ImportOnlySubgridV1
0     7   1    EmptySubgridV1
0     7   2    EmptySubgridV1
0     7   3    EmptySubgridV1
0     7   4    EmptySubgridV1
1     0   0    ImportOnlySubgridV1
1     0   1    ImportOnlySubgridV1
1     0   2    EmptySubgridV1
1     0   3    ImportOnlySubgridV1
1     0   4    EmptySubgridV1
1     1   0    ImportOnlySubgridV1
1     1   1    ImportOnlySubgridV1
1     1   2    EmptySubgridV1
1     1   3    ImportOnlySubgridV1
1     1   4    EmptySubgridV1
1     2   0    ImportOnlySubgridV1
1     2   1    ImportOnlySubgridV1
1     2   2    EmptySubgridV1
1     2   3    ImportOnlySubgridV1
1     2   4    EmptySubgridV1
1     3   0    ImportOnlySubgridV1
1     3   1    ImportOnlySubgridV1
1     3   2    EmptySubgridV1
1     3   3    ImportOnlySubgridV1
1     3   4    EmptySubgridV1
1     4   0    ImportOnlySubgridV1
1     4   1    ImportOnlySubgridV1
1     4   2    EmptySubgridV1
1     4   3    ImportOnlySubgridV1
1     4   4    EmptySubgridV1
1     5   0    ImportOnlySubgridV1
1     5   1    ImportOnlySubgridV1
1     5   2    EmptySubgridV1
1     5   3    ImportOnlySubgridV1
1     5   4    EmptySubgridV1
1     6   0    ImportOnlySubgridV1
1     6   1    ImportOnlySubgridV1
1     6   2    EmptySubgridV1
1     6   3    ImportOnlySubgridV1
1     6   4    EmptySubgridV1
1     7   0    ImportOnlySubgridV1
1     7   1    ImportOnlySubgridV1
1     7   2    EmptySubgridV1
1     7   3    ImportOnlySubgridV1
1     7   4    EmptySubgridV1
2     0   0    EmptySubgridV1
2     0   1    EmptySubgridV1
2     0   2    EmptySubgridV1
2     0   3    EmptySubgridV1
2     0   4    EmptySubgridV1
2     1   0    EmptySubgridV1
2     1   1    EmptySubgridV1
2     1   2    EmptySubgridV1
2     1   3    EmptySubgridV1
2     1   4    EmptySubgridV1
2     2   0    EmptySubgridV1
2     2   1    EmptySubgridV1
2     2   2    EmptySubgridV1
2     2   3    EmptySubgridV1
2     2   4    EmptySubgridV1
2     3   0    EmptySubgridV1
2     3   1    EmptySubgridV1
2     3   2    EmptySubgridV1
2     3   3    EmptySubgridV1
2     3   4    EmptySubgridV1
2     4   0    EmptySubgridV1
2     4   1    EmptySubgridV1
2     4   2    EmptySubgridV1
2     4   3    EmptySubgridV1
2     4   4    EmptySubgridV1
2     5   0    EmptySubgridV1
2     5   1    EmptySubgridV1
2     5   2    EmptySubgridV1
2     5   3    EmptySubgridV1
2     5   4    EmptySubgridV1
2     6   0    EmptySubgridV1
2     6   1    EmptySubgridV1
2     6   2    EmptySubgridV1
2     6   3    EmptySubgridV1
2     6   4    EmptySubgridV1
2     7   0    EmptySubgridV1
2     7   1    EmptySubgridV1
2     7   2    EmptySubgridV1
2     7   3    EmptySubgridV1
2     7   4    EmptySubgridV1
3     0   0    ImportOnlySubgridV1
3     0   1    ImportOnlySubgridV1
3     0   2    EmptySubgridV1
3     0   3    ImportOnlySubgridV1
3     0   4    EmptySubgridV1
3     1   0    ImportOnlySubgridV1
3     1   1    ImportOnlySubgridV1
3     1   2    EmptySubgridV1
3     1   3    ImportOnlySubgridV1
3     1   4    EmptySubgridV1
3     2   0    ImportOnlySubgridV1
3     2   1    ImportOnlySubgridV1
3     2   2    EmptySubgridV1
3     2   3    ImportOnlySubgridV1
3     2   4    EmptySubgridV1
3     3   0    ImportOnlySubgridV1
3     3   1    ImportOnlySubgridV1
3     3   2    EmptySubgridV1
3     3   3    ImportOnlySubgridV1
3     3   4    EmptySubgridV1
3     4   0    ImportOnlySubgridV1
3     4   1    ImportOnlySubgridV1
3     4   2    EmptySubgridV1
3     4   3    ImportOnlySubgridV1
3     4   4    EmptySubgridV1
3     5   0    ImportOnlySubgridV1
3     5   1    ImportOnlySubgridV1
3     5   2    EmptySubgridV1
3     5   3    ImportOnlySubgridV1
3     5   4    EmptySubgridV1
3     6   0    ImportOnlySubgridV1
3     6   1    ImportOnlySubgridV1
3     6   2    EmptySubgridV1
3     6   3    ImportOnlySubgridV1
3     6   4    EmptySubgridV1
3     7   0    ImportOnlySubgridV1
3     7   1    ImportOnlySubgridV1
3     7   2    EmptySubgridV1
3     7   3    ImportOnlySubgridV1
3     7   4    EmptySubgridV1
4     0   0    ImportOnlySubgridV1
4     0   1    EmptySubgridV1
4     0   2    ImportOnlySubgridV1
4     0   3    EmptySubgridV1
4     0   4    ImportOnlySubgridV1
4     1   0    ImportOnlySubgridV1
4     1   1    EmptySubgridV1
4     1   2    ImportOnlySubgridV1
4     1   3    EmptySubgridV1
4     1   4    ImportOnlySubgridV1
4     2   0    ImportOnlySubgridV1
4     2   1    EmptySubgridV1
4     2   2    ImportOnlySubgridV1
4     2   3    EmptySubgridV1
4     2   4    ImportOnlySubgridV1
4     3   0    ImportOnlySubgridV1
4     3   1    EmptySubgridV1
4     3   2    ImportOnlySubgridV1
4     3   3    EmptySubgridV1
4     3   4    ImportOnlySubgridV1
4     4   0    ImportOnlySubgridV1
4     4   1    EmptySubgridV1
4     4   2    ImportOnlySubgridV1
4     4   3    EmptySubgridV1
4     4   4    ImportOnlySubgridV1
4     5   0    ImportOnlySubgridV1
4     5   1    EmptySubgridV1
4     5   2    ImportOnlySubgridV1
4     5   3    EmptySubgridV1
4     5   4    ImportOnlySubgridV1
4     6   0    ImportOnlySubgridV1
4     6   1    EmptySubgridV1
4     6   2    ImportOnlySubgridV1
4     6   3    EmptySubgridV1
4     6   4    ImportOnlySubgridV1
4     7   0    ImportOnlySubgridV1
4     7   1    EmptySubgridV1
4     7   2    ImportOnlySubgridV1
4     7   3    EmptySubgridV1
4     7   4    ImportOnlySubgridV1
5     0   0    EmptySubgridV1
5     0   1    EmptySubgridV1
5     0   2    EmptySubgridV1
5     0   3    EmptySubgridV1
5     0   4    EmptySubgridV1
5     1   0    EmptySubgridV1
5     1   1    EmptySubgridV1
5     1   2    EmptySubgridV1
5     1   3    EmptySubgridV1
5     1   4    EmptySubgridV1
5     2   0    EmptySubgridV1
5     2   1    EmptySubgridV1
5     2   2    EmptySubgridV1
5     2   3    EmptySubgridV1
5     2   4    EmptySubgridV1
5     3   0    EmptySubgridV1
5     3   1    EmptySubgridV1
5     3   2    EmptySubgridV1
5     3   3    EmptySubgridV1
5     3   4    EmptySubgridV1
5     4   0    EmptySubgridV1
5     4   1    EmptySubgridV1
5     4   2    EmptySubgridV1
5     4   3    EmptySubgridV1
5     4   4    EmptySubgridV1
5     5   0    EmptySubgridV1
5     5   1    EmptySubgridV1
5     5   2    EmptySubgridV1
5     5   3    EmptySubgridV1
5     5   4    EmptySubgridV1
5     6   0    EmptySubgridV1
5     6   1    EmptySubgridV1
5     6   2    EmptySubgridV1
5     6   3    EmptySubgridV1
5     6   4    EmptySubgridV1
5     7   0    EmptySubgridV1
5     7   1    EmptySubgridV1
5     7   2    EmptySubgridV1
5     7   3    EmptySubgridV1
5     7   4    EmptySubgridV1
6     0   0    ImportOnlySubgridV1
6     0   1    EmptySubgridV1
6     0   2    ImportOnlySubgridV1
6     0   3    EmptySubgridV1
6     0   4    ImportOnlySubgridV1
6     1   0    ImportOnlySubgridV1
6     1   1    EmptySubgridV1
6     1   2    ImportOnlySubgridV1
6     1   3    EmptySubgridV1
6     1   4    ImportOnlySubgridV1
6     2   0    ImportOnlySubgridV1
6     2   1    EmptySubgridV1
6     2   2    ImportOnlySubgridV1
6     2   3    EmptySubgridV1
6     2   4    ImportOnlySubgridV1
6     3   0    ImportOnlySubgridV1
6     3   1    EmptySubgridV1
6     3   2    ImportOnlySubgridV1
6     3   3    EmptySubgridV1
6     3   4    ImportOnlySubgridV1
6     4   0    ImportOnlySubgridV1
6     4   1    EmptySubgridV1
6     4   2    ImportOnlySubgridV1
6     4   3    EmptySubgridV1
6     4   4    ImportOnlySubgridV1
6     5   0    ImportOnlySubgridV1
6     5   1    EmptySubgridV1
6     5   2    ImportOnlySubgridV1
6     5   3    EmptySubgridV1
6     5   4    ImportOnlySubgridV1
6     6   0    ImportOnlySubgridV1
6     6   1    EmptySubgridV1
6     6   2    ImportOnlySubgridV1
6     6   3    EmptySubgridV1
6     6   4    ImportOnlySubgridV1
6     7   0    ImportOnlySubgridV1
6     7   1    EmptySubgridV1
6     7   2    ImportOnlySubgridV1
6     7   3    EmptySubgridV1
6     7   4    ImportOnlySubgridV1
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

    #[test]
    fn default() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["subgrids", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(DEFAULT_STR);
    }

    #[test]
    fn show_empty() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["subgrids", "--show-empty", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(SHOW_EMPTY_STR);
    }
}
