use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use prettytable::{cell, row};
use std::path::PathBuf;

/// Shows the luminosity function.
#[derive(Parser)]
#[clap(aliases = &["lumis", "luminosities", "luminosity"])]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let grid = helpers::read_grid(&self.input)?;

        let mut table = helpers::create_table();
        let mut titles = row![c => "id"];
        for _ in 0..grid
            .lumi()
            .iter()
            .map(|lumi| lumi.entry().len())
            .max()
            .unwrap()
        {
            titles.add_cell(cell!(c->"entry"));
        }
        table.set_titles(titles);

        for (index, entry) in grid.lumi().iter().enumerate() {
            let row = table.add_empty_row();

            row.add_cell(cell!(&format!("{}", index)));

            for (id1, id2, factor) in entry.entry().iter() {
                row.add_cell(cell!(&format!(
                    "{} \u{d7} ({:2.}, {:2.})",
                    factor, id1, id2
                )));
            }
        }

        table.printstd();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-lumi 
Shows the luminosity function

USAGE:
    pineappl lumi <INPUT>

ARGS:
    <INPUT>    Path to the input grid

OPTIONS:
    -h, --help    Print help information
";

    const DEFAULT_STR: &str = "id    entry        entry
--+------------+------------
0  1 × ( 2, -1) 1 × ( 4, -3)
1  1 × ( 0, -3) 1 × ( 0, -1)
2  1 × (22, -3) 1 × (22, -1)
3  1 × ( 2,  0) 1 × ( 4,  0)
4  1 × ( 2, 22) 1 × ( 4, 22)
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["luminosity", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn default() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["luminosity", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(DEFAULT_STR);
    }
}
