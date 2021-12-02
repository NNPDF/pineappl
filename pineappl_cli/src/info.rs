use super::helpers;
use anyhow::Result;
use clap::{ArgGroup, Parser, ValueHint};
use itertools::Itertools;
use std::path::PathBuf;

/// Shows information about the grid.
#[derive(Parser)]
#[clap(group = ArgGroup::new("mode").required(true), name = "info")]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// For each order print a list of the largest EW order.
    #[clap(group = "mode", long)]
    ew: bool,
    /// Gets an internal key-value pair.
    #[clap(group = "mode", long, takes_value = true, value_name = "key")]
    get: Option<String>,
    /// Show all keys stored in the grid.
    #[clap(group = "mode", long)]
    keys: bool,
    /// For each order print a list of the largest QCD order.
    #[clap(group = "mode", long)]
    qcd: bool,
    /// Shows all key-value pairs stored in the grid.
    #[clap(group = "mode", long)]
    show: bool,
}

impl Opts {
    pub fn subcommand(&self) -> Result<()> {
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

            println!("{}", orders);
        } else if let Some(ref key) = self.get {
            grid.upgrade();

            if let Some(key_values) = grid.key_values() {
                if let Some(value) = key_values.get(key) {
                    println!("{}", value);
                }
            } else {
                unreachable!();
            }
        } else if self.keys {
            grid.upgrade();

            if let Some(key_values) = grid.key_values() {
                let mut vector = key_values.iter().collect::<Vec<_>>();
                vector.sort();

                for (key, _) in &vector {
                    println!("{}", key);
                }
            } else {
                unreachable!();
            }
        } else if self.show {
            grid.upgrade();

            if let Some(key_values) = grid.key_values() {
                let mut vector = key_values.iter().collect::<Vec<_>>();
                vector.sort();

                for (key, value) in &vector {
                    println!("{}: {}", key, value);
                }
            } else {
                unreachable!();
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-info 

Shows information about the grid

USAGE:
    pineappl info <--ew|--get <key>|--keys|--qcd|--show> <INPUT>

ARGS:
    <INPUT>    Path to the input grid

OPTIONS:
        --ew           For each order print a list of the largest EW order
        --get <key>    Gets an internal key-value pair
    -h, --help         Print help information
        --keys         Show all keys stored in the grid
        --qcd          For each order print a list of the largest QCD order
        --show         Shows all key-value pairs stored in the grid
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["info", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }
}
