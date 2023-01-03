use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueEnum, ValueHint};
use pineappl::fk_table::{FkAssumptions, FkTable};
use std::convert::TryFrom;
use std::path::PathBuf;

#[derive(Clone, Copy, ValueEnum)]
#[clap(rename_all = "PascalCase")]
enum FkAssum {
    Nf6Ind,
    Nf6Sym,
    Nf5Ind,
    Nf5Sym,
    Nf4Ind,
    Nf4Sym,
    Nf3Ind,
    Nf3Sym,
}

impl From<FkAssum> for FkAssumptions {
    fn from(assumptions: FkAssum) -> FkAssumptions {
        match assumptions {
            FkAssum::Nf6Ind => Self::Nf6Ind,
            FkAssum::Nf6Sym => Self::Nf6Sym,
            FkAssum::Nf5Ind => Self::Nf5Ind,
            FkAssum::Nf5Sym => Self::Nf5Sym,
            FkAssum::Nf4Ind => Self::Nf4Ind,
            FkAssum::Nf4Sym => Self::Nf4Sym,
            FkAssum::Nf3Ind => Self::Nf3Ind,
            FkAssum::Nf3Sym => Self::Nf3Sym,
        }
    }
}

/// Optimizes the internal data structure to minimize memory usage.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the optimized PineAPPL file.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    output: PathBuf,
    #[clap(long = "fk-table", value_enum, value_name = "ASSUMPTIONS")]
    fk_table: Option<FkAssum>,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<u8> {
        let mut grid = helpers::read_grid(&self.input)?;

        if let Some(assumptions) = self.fk_table {
            let mut fk_table = FkTable::try_from(grid)?;
            fk_table.optimize(assumptions.into());
            helpers::write_grid(&self.output, fk_table.grid())
        } else {
            grid.optimize();
            helpers::write_grid(&self.output, &grid)
        }
    }
}
