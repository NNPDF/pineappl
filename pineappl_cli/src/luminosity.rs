use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use prettytable::{cell, row};
use std::path::PathBuf;

/// Shows the luminosity function.
#[derive(Parser)]
#[clap(name = "lumi", aliases = &["luminosities", "luminosity"])]
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

    const HELP_STR: &str = "pineappl-luminosity 

Shows the luminosity function

USAGE:
    pineappl luminosity <INPUT>

ARGS:
    <INPUT>    Path to the input grid

OPTIONS:
    -h, --help    Print help information
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
}
