use super::helpers::{self, GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::{ArgGroup, Parser, ValueHint};
use itertools::Itertools;
use std::path::PathBuf;
use std::process::ExitCode;

/// Shows information about the grid.
#[derive(Parser)]
#[command(group = ArgGroup::new("mode").required(true))]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// For each order print a list of the largest EW order.
    #[arg(group = "mode", long)]
    ew: bool,
    /// Gets an internal key-value pair.
    #[arg(group = "mode", long, num_args = 1, value_name = "key")]
    get: Option<String>,
    /// Show all keys stored in the grid.
    #[arg(group = "mode", long)]
    keys: bool,
    /// For each order print a list of the largest QCD order.
    #[arg(group = "mode", long)]
    qcd: bool,
    /// Shows all key-value pairs stored in the grid.
    #[arg(group = "mode", long)]
    show: bool,
}

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<ExitCode> {
        let mut grid = helpers::read_grid(&self.input)?;

        if self.ew || self.qcd {
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
                    if self.qcd {
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
        } else if let Some(ref key) = self.get {
            grid.upgrade();

            grid.key_values().map_or_else(
                || unreachable!(),
                |key_values| {
                    if let Some(value) = key_values.get(key) {
                        println!("{value}");
                    }
                },
            );
        } else if self.keys {
            grid.upgrade();

            grid.key_values().map_or_else(
                || unreachable!(),
                |key_values| {
                    let mut vector = key_values.iter().collect::<Vec<_>>();
                    vector.sort();

                    for (key, _) in &vector {
                        println!("{key}");
                    }
                },
            );
        } else if self.show {
            grid.upgrade();

            grid.key_values().map_or_else(
                || unreachable!(),
                |key_values| {
                    let mut vector = key_values.iter().collect::<Vec<_>>();
                    vector.sort();

                    for (key, value) in &vector {
                        println!("{key}: {value}");
                    }
                },
            );
        }

        Ok(ExitCode::SUCCESS)
    }
}
