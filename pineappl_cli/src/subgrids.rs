use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{ArgGroup, Parser, ValueHint};
use pineappl::subgrid::Mu2;
use pineappl::subgrid::{Subgrid, SubgridEnum};
use prettytable::{cell, row};
use std::path::PathBuf;

/// Print information about the internal subgrid types.
#[derive(Parser)]
#[clap(group = ArgGroup::new("show").multiple(true).required(true))]
pub struct Opts {
    /// Path to the input grid.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Show empty subgrids.
    #[clap(long = "show-empty")]
    show_empty: bool,
    /// Show the subgrid type.
    #[clap(group = "show", long = "type")]
    type_: bool,
    /// Show the renormalization grid values.
    #[clap(group = "show", long)]
    mur: bool,
    /// Show the squared renormalization grid values.
    #[clap(group = "show", long)]
    mur2: bool,
    /// Show the factorization grid values.
    #[clap(group = "show", long)]
    muf: bool,
    /// Show the squared factorization grid values.
    #[clap(group = "show", long)]
    muf2: bool,
    /// Show the x1 grid values.
    #[clap(group = "show", long)]
    x1: bool,
    /// Show the x2 grid values.
    #[clap(group = "show", long)]
    x2: bool,
    /// Show grid statistics (figures are the number of entries).
    #[clap(group = "show", long)]
    stats: bool,
    /// Set the number of digits shown for numerical values.
    #[clap(default_value_t = 3, long = "digits")]
    digits: usize,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<u8> {
        let grid = helpers::read_grid(&self.input)?;
        let mut table = helpers::create_table();
        let mut titles = row![c => "o", "b", "l"];

        if self.type_ {
            titles.add_cell(cell!(c->"type"));
        }
        if self.mur {
            titles.add_cell(cell!(c->"mur"));
        }
        if self.mur2 {
            titles.add_cell(cell!(c->"mur2"));
        }
        if self.muf {
            titles.add_cell(cell!(c->"muf"));
        }
        if self.muf2 {
            titles.add_cell(cell!(c->"muf2"));
        }
        if self.x1 {
            titles.add_cell(cell!(c->"x1"));
        }
        if self.x2 {
            titles.add_cell(cell!(c->"x2"));
        }
        if self.stats {
            titles.add_cell(cell!(c->"total"));
            titles.add_cell(cell!(c->"allocated"));
            titles.add_cell(cell!(c->"zeros"));
            titles.add_cell(cell!(c->"overhead"));
        }
        table.set_titles(titles);

        for ((order, bin, lumi), subgrid) in grid.subgrids().indexed_iter() {
            if !self.show_empty && subgrid.is_empty() {
                continue;
            }

            let row = table.add_empty_row();

            row.add_cell(cell!(l->format!("{order}")));
            row.add_cell(cell!(l->format!("{bin}")));
            row.add_cell(cell!(l->format!("{lumi}")));

            if self.type_ {
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
            if self.mur {
                let values: Vec<_> = subgrid
                    .mu2_grid()
                    .iter()
                    .map(|Mu2 { ren, fac: _ }| format!("{:.*}", self.digits, ren.sqrt()))
                    .collect();

                row.add_cell(cell!(l->values.join(", ")));
            }
            if self.mur2 {
                let values: Vec<_> = subgrid
                    .mu2_grid()
                    .iter()
                    .map(|Mu2 { ren, fac: _ }| format!("{:.*}", self.digits, ren))
                    .collect();

                row.add_cell(cell!(l->values.join(", ")));
            }
            if self.muf {
                let values: Vec<_> = subgrid
                    .mu2_grid()
                    .iter()
                    .map(|Mu2 { ren: _, fac }| format!("{:.*}", self.digits, fac.sqrt()))
                    .collect();

                row.add_cell(cell!(l->values.join(", ")));
            }
            if self.muf2 {
                let values: Vec<_> = subgrid
                    .mu2_grid()
                    .iter()
                    .map(|Mu2 { ren: _, fac }| format!("{:.*}", self.digits, fac))
                    .collect();

                row.add_cell(cell!(l->values.join(", ")));
            }
            if self.x1 {
                let values: Vec<_> = subgrid
                    .x1_grid()
                    .iter()
                    .map(|x| format!("{:.*e}", self.digits, x))
                    .collect();

                row.add_cell(cell!(l->values.join(", ")));
            }
            if self.x2 {
                let values: Vec<_> = subgrid
                    .x2_grid()
                    .iter()
                    .map(|x| format!("{:.*e}", self.digits, x))
                    .collect();

                row.add_cell(cell!(l->values.join(", ")));
            }
            if self.stats {
                let stats = subgrid.stats();
                row.add_cell(cell!(r->stats.total.to_string()));
                row.add_cell(cell!(r->stats.allocated.to_string()));
                row.add_cell(cell!(r->stats.zeros.to_string()));
                row.add_cell(cell!(r->stats.overhead.to_string()));
            }
        }

        table.printstd();

        Ok(0)
    }
}
