use super::helpers::{self, GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::builder::{PossibleValuesParser, TypedValueParser};
use clap::{Parser, ValueHint};
use pineappl::fk_table::{FkAssumptions, FkTable};
use std::path::PathBuf;
use std::process::ExitCode;

/// Optimizes the internal data structure to minimize memory usage.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the optimized PineAPPL file.
    #[arg(value_hint = ValueHint::FilePath)]
    output: PathBuf,
    #[arg(
        long = "fk-table",
        value_name = "ASSUMPTIONS",
        value_parser = PossibleValuesParser::new(["Nf6Ind", "Nf6Sym", "Nf5Ind", "Nf5Sym", "Nf4Ind", "Nf4Sym", "Nf3Ind", "Nf3Sym"]).try_map(|s| s.parse::<FkAssumptions>())
    )]
    fk_table: Option<FkAssumptions>,
}

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<ExitCode> {
        let mut grid = helpers::read_grid(&self.input)?;

        if let Some(assumptions) = self.fk_table {
            let mut fk_table = FkTable::try_from(grid)?;
            fk_table.optimize(assumptions);
            helpers::write_grid(&self.output, fk_table.grid())
        } else {
            grid.optimize();
            helpers::write_grid(&self.output, &grid)
        }
    }
}
