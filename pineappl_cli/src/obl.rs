use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{ArgGroup, Parser, ValueHint};
use itertools::Itertools;
use pineappl::grid::Order;
use prettytable::{cell, row, Row};
use std::path::PathBuf;

/// Shows information about orders (o), bins (b), or luminosities (l) of a grid.
#[derive(Parser)]
#[clap(group = ArgGroup::new("mode").required(true))]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
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
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let grid = helpers::read_grid(&self.input)?;

        let mut table = helpers::create_table();

        if self.bins {
            let mut titles = Row::empty();
            titles.add_cell(cell!(c->"b"));

            for x_label in helpers::labels(&grid, false).split_last().unwrap().1 {
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

                row.add_cell(cell!(format!("{}", index)));

                for (id1, id2, factor) in entry.entry().iter() {
                    row.add_cell(cell!(format!("{} \u{d7} ({:2.}, {:2.})", factor, id1, id2)));
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
                row.add_cell(cell!(format!("O({})", order_string)));
            }
        }

        table.printstd();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-obl 
Shows information about orders (o), bins (b), or luminosities (l) of a grid

USAGE:
    pineappl obl <--orders|--orders-spaces|--orders-long|--bins|--lumis> <INPUT>

ARGS:
    <INPUT>    Path to the input grid

OPTIONS:
    -b, --bins             Show the bins of a grid
    -h, --help             Print help information
    -l, --lumis            Show the luminsities a grid
    -o, --orders           Show the orders of a grid, stripping zero powers
        --orders-long      Show the orders of a grid, including zero powers
        --orders-spaces    Show the orders of a grid, replacing zero powers with spaces
";

    const BINS_STR: &str = "b   etal    norm
-+----+----+----
0    2 2.25 0.25
1 2.25  2.5 0.25
2  2.5 2.75 0.25
3 2.75    3 0.25
4    3 3.25 0.25
5 3.25  3.5 0.25
6  3.5    4  0.5
7    4  4.5  0.5
";

    const LUMIS_STR: &str = "l    entry        entry
-+------------+------------
0 1 × ( 2, -1) 1 × ( 4, -3)
1 1 × ( 0, -3) 1 × ( 0, -1)
2 1 × (22, -3) 1 × (22, -1)
3 1 × ( 2,  0) 1 × ( 4,  0)
4 1 × ( 2, 22) 1 × ( 4, 22)
";

    const ORDERS_STR: &str = "o      order
-+----------------
0 O(a^2)
1 O(as^1 a^2)
2 O(as^1 a^2 lr^1)
3 O(as^1 a^2 lf^1)
4 O(a^3)
5 O(a^3 lr^1)
6 O(a^3 lf^1)
";

    const ORDERS_LONG_STR: &str = "o         order
-+---------------------
0 O(as^0 a^2 lr^0 lf^0)
1 O(as^1 a^2 lr^0 lf^0)
2 O(as^1 a^2 lr^1 lf^0)
3 O(as^1 a^2 lr^0 lf^1)
4 O(as^0 a^3 lr^0 lf^0)
5 O(as^0 a^3 lr^1 lf^0)
6 O(as^0 a^3 lr^0 lf^1)
";

    const ORDERS_SPACES_STR: &str = "o         order
-+---------------------
0 O(     a^2          )
1 O(as^1 a^2          )
2 O(as^1 a^2 lr^1     )
3 O(as^1 a^2      lf^1)
4 O(     a^3          )
5 O(     a^3 lr^1     )
6 O(     a^3      lf^1)
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["obl", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn bins() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["obl", "--bins", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(BINS_STR);
    }

    #[test]
    fn lumis() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["obl", "--lumis", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(LUMIS_STR);
    }

    #[test]
    fn orders() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["obl", "--orders", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(ORDERS_STR);
    }

    #[test]
    fn orders_long() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["obl", "--orders-long", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(ORDERS_LONG_STR);
    }

    #[test]
    fn orders_spaces() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["obl", "--orders-spaces", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(ORDERS_SPACES_STR);
    }
}
