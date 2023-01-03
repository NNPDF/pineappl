use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{ArgGroup, Parser, ValueHint};
use itertools::Itertools;
use pineappl::fk_table::FkTable;
use pineappl::grid::Order;
use prettytable::{cell, row, Row};
use std::path::PathBuf;

/// Shows information about orders (o), bins (b), or luminosities (l) of a grid.
#[derive(Parser)]
#[clap(group = ArgGroup::new("mode").required(true))]
pub struct Opts {
    /// Path to the input grid.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Show the orders of a grid, stripping zero powers.
    #[clap(group = "mode", long, short)]
    orders: bool,
    /// Show the orders of a grid, replacing zero powers with spaces.
    #[clap(group = "mode", long)]
    orders_spaces: bool,
    /// Show the orders of a grid, including zero powers.
    #[clap(group = "mode", long)]
    orders_long: bool,
    /// Show the bins of a grid.
    #[clap(group = "mode", long, short)]
    bins: bool,
    /// Show the luminsities a grid.
    #[clap(group = "mode", long, short)]
    lumis: bool,
    /// Check if input is an FK table.
    #[clap(group = "mode", long)]
    fktable: bool,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<u8> {
        let grid = helpers::read_grid(&self.input)?;

        let mut table = helpers::create_table();

        if self.bins {
            let mut titles = Row::empty();
            titles.add_cell(cell!(c->"b"));

            for (x_label, _) in helpers::labels_and_units(&grid, false).0 {
                let mut cell = cell!(c->x_label);
                cell.set_hspan(2);
                titles.add_cell(cell);
            }
            titles.add_cell(cell!(c->"norm"));

            table.set_titles(titles);

            let left_limits: Vec<_> = (0..grid.bin_info().dimensions())
                .map(|i| grid.bin_info().left(i))
                .collect();
            let right_limits: Vec<_> = (0..grid.bin_info().dimensions())
                .map(|i| grid.bin_info().right(i))
                .collect();
            let normalizations = grid.bin_info().normalizations();

            for bin in 0..grid.bin_info().bins() {
                let row = table.add_empty_row();
                row.add_cell(cell!(bin.to_string()));

                for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                    row.add_cell(cell!(r->format!("{}", left[bin])));
                    row.add_cell(cell!(r->format!("{}", right[bin])));
                }

                row.add_cell(cell!(r->format!("{}", normalizations[bin])));
            }
        } else if self.fktable {
            if let Err(err) = FkTable::try_from(grid) {
                println!("no\n{err}");
                return Ok(1);
            }

            println!("yes");
            return Ok(0);
        } else if self.lumis {
            let mut titles = row![c => "l"];
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

                row.add_cell(cell!(format!("{index}")));

                for (id1, id2, factor) in entry.entry().iter() {
                    row.add_cell(cell!(format!("{factor} \u{d7} ({id1:2}, {id2:2})")));
                }
            }
        } else {
            table.set_titles(row![c => "o", "order"]);

            for (index, order) in grid.orders().iter().enumerate() {
                let row = table.add_empty_row();

                let Order {
                    alphas,
                    alpha,
                    logxir,
                    logxif,
                } = order;

                let order_string = [alphas, alpha, logxir, logxif]
                    .iter()
                    .zip(["as^", "a^", "lr^", "lf^"].iter())
                    .filter_map(|(num, string)| {
                        if **num == 0 && self.orders {
                            None
                        } else if **num == 0 && self.orders_spaces {
                            Some(" ".repeat(string.len() + 1))
                        } else {
                            let mut result = (*string).to_string();
                            result.push_str(&num.to_string());
                            Some(result)
                        }
                    })
                    .join(" ");

                row.add_cell(cell!(index.to_string()));
                row.add_cell(cell!(format!("O({order_string})")));
            }
        }

        table.printstd();

        Ok(0)
    }
}
