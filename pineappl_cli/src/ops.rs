use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use pineappl::lumi::LumiEntry;
use pineappl::pids;
use std::ops::Deref;
use std::path::PathBuf;

/// A collection of various modifying operations on grids.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path of the modified PineAPPL file.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// Charge conjugate the first initial state.
    #[clap(long)]
    cc1: bool,
    /// Charge conjugate the second initial state.
    #[clap(long)]
    cc2: bool,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;

        if self.cc1 || self.cc2 {
            let lumi_id_types = grid.key_values().map_or("pdg_mc_ids", |kv| {
                kv.get("lumi_id_types").map_or("pdg_mc_ids", Deref::deref)
            });
            let lumis = grid
                .lumi()
                .iter()
                .map(|entry| {
                    LumiEntry::new(
                        entry
                            .entry()
                            .iter()
                            .map(|&(a, b, f)| {
                                let (ap, f1) = if self.cc1 {
                                    pids::charge_conjugate(lumi_id_types, a)
                                } else {
                                    (a, 1.0)
                                };
                                let (bp, f2) = if self.cc2 {
                                    pids::charge_conjugate(lumi_id_types, b)
                                } else {
                                    (b, 1.0)
                                };
                                (ap, bp, f * f1 * f2)
                            })
                            .collect(),
                    )
                })
                .collect();

            let mut initial_state_1: i32 = grid
                .key_values()
                .map_or("2212", |kv| &kv["initial_state_1"])
                .parse()?;
            let mut initial_state_2: i32 = grid
                .key_values()
                .map_or("2212", |kv| &kv["initial_state_2"])
                .parse()?;

            if self.cc1 {
                initial_state_1 = pids::charge_conjugate_pdg_pid(initial_state_1);
            }
            if self.cc2 {
                initial_state_2 = pids::charge_conjugate_pdg_pid(initial_state_2);
            }

            grid.set_key_value("initial_state_1", &initial_state_1.to_string());
            grid.set_key_value("initial_state_2", &initial_state_2.to_string());
            grid.set_lumis(lumis);
        }

        helpers::write_grid(&self.output, &grid)
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;
    use assert_fs::NamedTempFile;

    const HELP_STR: &str = "pineappl-ops 
A collection of various modifying operations on grids

USAGE:
    pineappl ops [OPTIONS] <INPUT> <OUTPUT>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path of the modified PineAPPL file

OPTIONS:
        --cc1     Charge conjugate the first initial state
        --cc2     Charge conjugate the second initial state
    -h, --help    Print help information
";

    const DEFAULT_STR: &str = "b   etal    disg/detal  scale uncertainty
-+----+----+-----------+--------+--------
0    2 2.25 3.7527620e2   -3.77%    2.71%
1 2.25  2.5 3.4521553e2   -3.79%    2.80%
2  2.5 2.75 3.0001406e2   -3.78%    2.86%
3 2.75    3 2.4257663e2   -3.77%    2.92%
4    3 3.25 1.8093343e2   -3.74%    2.95%
5 3.25  3.5 1.2291115e2   -3.71%    2.98%
6  3.5    4 5.7851018e1   -3.63%    2.97%
7    4  4.5 1.3772029e1   -3.46%    2.85%
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["ops", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn cc1() {
        let output = NamedTempFile::new("cc1.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "ops",
                "--cc1",
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
            .stdout(DEFAULT_STR);
    }

    #[test]
    fn cc2() {
        let output = NamedTempFile::new("cc2.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "ops",
                "--cc2",
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
            .stdout(DEFAULT_STR);
    }
}
