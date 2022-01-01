use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use std::path::PathBuf;

/// Merges one or more PineAPPL grids together.
#[derive(Parser)]
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

    const DEFAULT_STR: &str = "bin   etal    disg/detal  neg unc pos unc
---+----+----+-----------+-------+-------
  0    2 2.25 3.7527620e2  -3.77%   2.71%
  1 2.25  2.5 3.4521553e2  -3.79%   2.80%
  2  2.5 2.75 3.0001406e2  -3.78%   2.86%
  3 2.75    3 2.4257663e2  -3.77%   2.92%
  4    3 3.25 1.8093343e2  -3.74%   2.95%
  5 3.25  3.5 1.2291115e2  -3.71%   2.98%
  6  3.5    4 5.7851018e1  -3.63%   2.97%
  7    4  4.5 1.3772029e1  -3.46%   2.85%
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

    #[test]
    fn default() {
        let output = NamedTempFile::new("merged.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "merge",
                "--scale=0.5",
                output.path().to_str().unwrap(),
                "data/LHCB_WP_7TEV.pineappl.lz4",
                "data/LHCB_WP_7TEV.pineappl.lz4",
            ])
            .assert()
            .success()
            .stdout("");

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                output.path().to_str().unwrap(),
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout(DEFAULT_STR);
    }
}
