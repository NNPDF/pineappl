use super::helpers::{self, Subcommand};
use anyhow::{bail, Result};
use clap::{ArgGroup, Parser, ValueHint};
use pineappl::bin::BinRemapper;
use std::ops::RangeInclusive;
use std::path::PathBuf;

/// Sums two or more bins of a grid together.
#[derive(Parser)]
#[clap(group = ArgGroup::new("mode").required(true))]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the modified PineAPPL file.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// Sums all bins into a single bin.
    #[clap(long, group = "mode")]
    integrated: bool,
    /// Merge specific bins together.
    #[clap(
        long,
        group = "mode",
        parse(try_from_str = helpers::try_parse_integer_range),
        short,
        use_value_delimiter = true,
    )]
    bins: Option<RangeInclusive<usize>>,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;

        if self.integrated {
            if grid.merge_bins(0..grid.bin_info().bins()).is_err() {
                bail!("TODO");
            }
            grid.set_remapper(
                BinRemapper::new(vec![1.0], vec![(0.0, 1.0)]).unwrap_or_else(|_| unreachable!()),
            )?;

            let dimensions = grid.bin_info().dimensions();
            let key_values = grid.key_values_mut();
            for dim in 0..dimensions {
                key_values.remove(&format!("x{}_label", dim + 1));
                key_values.remove(&format!("x{}_label_tex", dim + 1));
                key_values.remove(&format!("x{}_unit", dim + 1));
            }
            key_values.remove("y_label");
            key_values.remove("y_label_tex");
            key_values.remove("y_unit");
        } else if let Some(range) = self.bins.as_ref() {
            if grid.merge_bins(*range.start()..(range.end() + 1)).is_err() {
                bail!("TODO");
            }
        } else {
            unreachable!();
        }

        helpers::write_grid(&self.output, &grid)
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;
    use assert_fs::NamedTempFile;

    const HELP_STR: &str = "pineappl-sum 
Sums two or more bins of a grid together

USAGE:
    pineappl sum <--integrated|--bins <BINS>> <INPUT> <OUTPUT>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path to the modified PineAPPL file

OPTIONS:
    -b, --bins <BINS>    Merge specific bins together
    -h, --help           Print help information
        --integrated     Sums all bins into a single bin
";

    const BINS_STR: &str = "bin   etal    disg/detal  scale uncertainty
---+----+----+-----------+--------+--------
  0    2 2.25 3.7527620e2   -3.77%    2.71%
  1 2.25  2.5 3.4521553e2   -3.79%    2.80%
  2  2.5 2.75 3.0001406e2   -3.78%    2.86%
  3 2.75    3 2.4257663e2   -3.77%    2.92%
  4    3 3.25 1.8093343e2   -3.74%    2.95%
  5 3.25  3.5 1.2291115e2   -3.71%    2.98%
  6  3.5  4.5 3.5811524e1   -3.59%    2.95%
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["sum", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn bins() {
        let output = NamedTempFile::new("bins.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "sum",
                "--bins=6-7",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                output.path().to_str().unwrap(),
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
            .stdout(BINS_STR);
    }

    #[test]
    fn integrated() {
        let output = NamedTempFile::new("summed.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "sum",
                "--integrated",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                output.path().to_str().unwrap(),
            ])
            .assert()
            .success()
            .stdout("");
    }
}
