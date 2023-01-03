use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use pineappl::fk_table::{FkAssumptions, FkTable};
use std::convert::TryFrom;
use std::path::PathBuf;

/// Optimizes the internal data structure to minimize memory usage.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the optimized PineAPPL file.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    output: PathBuf,
    #[clap(
        long = "fk-table",
        parse(try_from_str = FkAssumptions::try_from),
        possible_values = &[
            "Nf6Ind",
            "Nf6Sym",
            "Nf5Ind",
            "Nf5Sym",
            "Nf4Ind",
            "Nf4Sym",
            "Nf3Ind",
            "Nf3Sym",
        ],
        value_name = "ASSUMPTIONS",
    )]
    fk_table: Option<FkAssumptions>,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<u8> {
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
