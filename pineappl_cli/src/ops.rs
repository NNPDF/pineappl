use super::helpers::{self, GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use pineappl::lumi::LumiEntry;
use pineappl::pids;
use std::ops::Deref;
use std::path::PathBuf;
use std::process::ExitCode;

/// A collection of various modifying operations on grids.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path of the modified PineAPPL file.
    #[arg(value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// Charge conjugate the first initial state.
    #[arg(long)]
    cc1: bool,
    /// Charge conjugate the second initial state.
    #[arg(long)]
    cc2: bool,
    /// Scale each bin with a different factor.
    #[arg(
        long,
        num_args = 1,
        value_delimiter = ',',
        value_name = "BIN1,BIN2,..."
    )]
    scale_by_bin: Vec<f64>,
}

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<ExitCode> {
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

        if !self.scale_by_bin.is_empty() {
            grid.scale_by_bin(&self.scale_by_bin);
        }

        helpers::write_grid(&self.output, &grid)
    }
}
