use super::helpers;
use super::{GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::{Args, Parser, ValueHint};
use itertools::Itertools;
use pineappl::boc::Order;
use pineappl::fk_table::FkTable;
use prettytable::{cell, row, Row};
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Args)]
#[group(multiple = false, required = true)]
struct Group {
    /// Show the bins of a grid.
    #[arg(long, short)]
    bins: bool,
    /// Show the channel definition of a grid.
    #[arg(alias = "lumis", long)]
    channels: bool,
    /// For each order print a list of the largest EW order.
    #[arg(long)]
    ew: bool,
    /// Check if input is an FK table.
    #[arg(long)]
    fk_table: bool,
    /// Return the (squared) factorization scale of the FK-table.
    #[arg(long)]
    fk_table_fac0: bool,
    /// Return the (squared) fragmentation scale of the FK-table.
    #[arg(long)]
    fk_table_frg0: bool,
    /// Gets an internal key-value pair.
    #[arg(long, num_args = 1, value_name = "KEY")]
    get: Option<String>,
    /// Show all keys stored in the grid.
    #[arg(long)]
    keys: bool,
    /// Show the orders of a grid, stripping zero powers.
    #[arg(long, short)]
    orders: bool,
    /// Show the orders of a grid, including zero powers.
    #[arg(long)]
    orders_long: bool,
    /// Show the orders of a grid, replacing zero powers with spaces.
    #[arg(long)]
    orders_spaces: bool,
    /// For each order print a list of the largest QCD order.
    #[arg(long)]
    qcd: bool,
    /// Shows all key-value pairs stored in the grid.
    #[arg(long)]
    show: bool,
}

/// Read out information of a grid.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    #[command(flatten)]
    group: Group,
}

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<ExitCode> {
        let grid = helpers::read_grid(&self.input)?;

        let mut table = helpers::create_table();

        if self.group.bins {
            let mut titles = Row::empty();
            titles.add_cell(cell!(c->"b"));

            for (x_label, _) in helpers::labels_and_units(&grid, false).0 {
                let mut cell = cell!(c->x_label);
                cell.set_hspan(2);
                titles.add_cell(cell);
            }
            titles.add_cell(cell!(c->"norm"));

            table.set_titles(titles);

            for (bin_index, bin) in grid.bwfl().bins().iter().enumerate() {
                let row = table.add_empty_row();
                row.add_cell(cell!(bin_index.to_string()));

                for (left, right) in bin.limits() {
                    row.add_cell(cell!(r->left.to_string()));
                    row.add_cell(cell!(r->right.to_string()));
                }

                row.add_cell(cell!(r->bin.normalization().to_string()));
            }
        } else if self.group.channels {
            let mut titles = row![c => "c"];

            // if there are no channels print at least one column
            for _ in 0..grid
                .channels()
                .iter()
                .map(|channel| channel.entry().len())
                .max()
                .unwrap_or(1)
            {
                titles.add_cell(cell!(c->"entry"));
            }
            table.set_titles(titles);

            for (index, channel) in grid.channels().iter().enumerate() {
                let row = table.add_empty_row();

                row.add_cell(cell!(format!("{index}")));

                for (pids, factor) in channel.entry() {
                    row.add_cell(cell!(format!(
                        "{factor} \u{d7} ({})",
                        pids.iter().map(|pid| format!("{pid:2}")).join(", ")
                    )));
                }
            }
        } else if self.group.ew || self.group.qcd {
            let mut sorted_grid_orders: Vec<_> = grid
                .orders()
                .iter()
                .filter(|order| (order.logxir == 0) && (order.logxif == 0))
                .collect();
            sorted_grid_orders.sort();

            let orders = sorted_grid_orders
                .into_iter()
                .group_by(|order| order.alphas + order.alpha)
                .into_iter()
                .map(|mut iter| {
                    if self.group.qcd {
                        iter.1.next().unwrap()
                    } else {
                        iter.1.last().unwrap()
                    }
                })
                .map(|order| {
                    if order.alphas == 0 {
                        format!("a{}", order.alpha)
                    } else if order.alpha == 0 {
                        format!("as{}", order.alphas)
                    } else {
                        format!("as{}a{}", order.alphas, order.alpha)
                    }
                })
                .collect::<Vec<_>>()
                .join(",");

            println!("{orders}");
        } else if self.group.fk_table {
            if let Err(err) = FkTable::try_from(grid) {
                println!("no\n{err}");
                return Ok(ExitCode::FAILURE);
            }

            println!("yes");
            return Ok(ExitCode::SUCCESS);
        } else if self.group.fk_table_fac0 {
            let fk_table = FkTable::try_from(grid)?;

            if let Some(fac0) = fk_table.fac0() {
                println!("{fac0}");
            } else {
                println!("None");
            }

            return Ok(ExitCode::SUCCESS);
        } else if self.group.fk_table_frg0 {
            let fk_table = FkTable::try_from(grid)?;

            if let Some(frg0) = fk_table.frg0() {
                println!("{frg0}");
            } else {
                println!("None");
            }

            return Ok(ExitCode::SUCCESS);
        } else if let Some(key) = &self.group.get {
            if let Some(value) = grid.metadata().get(key) {
                println!("{value}");
            }
        } else if self.group.keys {
            for key in grid.metadata().keys() {
                println!("{key}");
            }
        } else if self.group.orders || self.group.orders_spaces || self.group.orders_long {
            table.set_titles(row![c => "o", "order"]);

            for (index, order) in grid.orders().iter().enumerate() {
                let row = table.add_empty_row();

                let Order {
                    alphas,
                    alpha,
                    logxir,
                    logxif,
                    logxia,
                } = order;

                let order_string = [alphas, alpha, logxir, logxif, logxia]
                    .iter()
                    .zip(["as^", "a^", "lr^", "lf^", "la^"].iter())
                    .filter_map(|(num, string)| {
                        if **num == 0 && self.group.orders {
                            None
                        } else if **num == 0 && self.group.orders_spaces {
                            Some(" ".repeat(string.len() + 1))
                        } else {
                            let mut result = (*string).to_owned();
                            result.push_str(&num.to_string());
                            Some(result)
                        }
                    })
                    .join(" ");

                row.add_cell(cell!(index.to_string()));
                row.add_cell(cell!(format!("O({order_string})")));
            }
        } else if self.group.show {
            for (key, value) in grid.metadata() {
                println!("{key}: {value}");
            }
        }

        table.printstd();

        Ok(ExitCode::SUCCESS)
    }
}
