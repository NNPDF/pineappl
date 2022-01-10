use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{ArgGroup, Parser, ValueHint};
use std::ops::RangeInclusive;
use std::path::PathBuf;

/// Print information about the internal subgrid types.
#[derive(Parser)]
#[clap(group = ArgGroup::new("mode").multiple(true).required(true))]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the modified PineAPPL file.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
    #[clap(
        group = "mode",
        long,
        multiple_values = true,
        parse(try_from_str = helpers::try_parse_integer_range),
        use_delimiter = true
    )]
    /// Indices of bins that should be deleted.
    bins: Vec<RangeInclusive<usize>>,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;
        let bins: Vec<usize> = self
            .bins
            .iter()
            .flat_map(|range| range.clone().collect::<Vec<_>>())
            .collect();

        grid.delete_bins(&bins);

        helpers::write_grid(&self.output, &grid)
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;
    use assert_fs::NamedTempFile;

    const HELP_STR: &str = "pineappl-delete 
Print information about the internal subgrid types

USAGE:
    pineappl delete <--bins <BINS>...> <INPUT> <OUTPUT>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path to the modified PineAPPL file

OPTIONS:
        --bins <BINS>...    Indices of bins that should be deleted
    -h, --help              Print help information
";

    const BINS_02_57_STR: &str = "bin   etal    disg/detal  neg unc pos unc
---+----+----+-----------+-------+-------
  0 2.75    3 2.4257663e2  -3.77%   2.92%
  1    3 3.25 1.8093343e2  -3.74%   2.95%
";

    const BINS_25_STR: &str = "bin   etal    disg/detal  neg unc pos unc
---+----+----+-----------+-------+-------
  0    2 2.25 3.7527620e2  -3.77%   2.71%
  1 2.25  2.5 3.4521553e2  -3.79%   2.80%
  2  3.5    4 5.7851018e1  -3.63%   2.97%
  3    4  4.5 1.3772029e1  -3.46%   2.85%
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["delete", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn bins_02_57() {
        let output = NamedTempFile::new("deleted.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "delete",
                "--bins=0-2,5-7",
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
            .stdout(BINS_02_57_STR);
    }

    #[test]
    fn bins_25() {
        let output = NamedTempFile::new("deleted2.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "delete",
                "--bins=2-5",
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
            .stdout(BINS_25_STR);
    }
}
