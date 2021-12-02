use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use std::path::PathBuf;

/// Merges one or more PineAPPL grids together.
#[derive(Parser)]
#[clap(name = "merge")]
pub struct Opts {
    /// Path of the merged PineAPPL file.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// Path(s) of the files that should be merged.
    #[clap(min_values = 1, parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: Vec<PathBuf>,
    /// Scales all grids with the given factor.
    #[clap(long, short)]
    scale: Option<f64>,
    /// Scales all grids with order-dependent factors.
    #[clap(
        alias = "scale_by_order",
        conflicts_with = "scale",
        long = "scale-by-order",
        value_names = &["ALPHAS", "ALPHA", "LOGXIR", "LOGXIF", "GLOBAL"]
    )]
    scale_by_order: Vec<f64>,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let (input0, input_rest) = self.input.split_first().unwrap();
        let mut grid0 = helpers::read_grid(input0)?;

        for i in input_rest {
            grid0.merge(helpers::read_grid(i)?)?;
        }

        if let Some(scale) = self.scale {
            grid0.scale(scale);
        } else if !self.scale_by_order.is_empty() {
            grid0.scale_by_order(
                self.scale_by_order[0],
                self.scale_by_order[1],
                self.scale_by_order[2],
                self.scale_by_order[3],
                self.scale_by_order[4],
            );
        }

        helpers::write_grid(&self.output, &grid0)
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-merge 

Merges one or more PineAPPL grids together

USAGE:
    pineappl merge [OPTIONS] <OUTPUT> [--] [INPUT]...

ARGS:
    <OUTPUT>      Path of the merged PineAPPL file
    <INPUT>...    Path(s) of the files that should be merged

OPTIONS:
    -h, --help
            Print help information

    -s, --scale <SCALE>
            Scales all grids with the given factor

        --scale-by-order <ALPHAS> <ALPHA> <LOGXIR> <LOGXIF> <GLOBAL>
            Scales all grids with order-dependent factors
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["merge", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }
}
