use super::helpers::{self, GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::{
    value_parser, Arg, ArgAction, ArgMatches, Args, Command, Error, FromArgMatches, Parser,
    ValueHint,
};
use pineappl::lumi::LumiEntry;
use pineappl::pids;
use std::any::Any;
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
    #[command(flatten)]
    more_args: MoreArgs,
}

struct MoreArgs {
    stack: Vec<(String, Box<dyn Any>)>,
}

impl FromArgMatches for MoreArgs {
    fn from_arg_matches(_: &ArgMatches) -> Result<Self, Error> {
        unreachable!()
    }

    fn from_arg_matches_mut(matches: &mut ArgMatches) -> Result<Self, Error> {
        let mut stack = Vec::new();
        let ids: Vec<_> = matches.ids().map(|id| id.as_str().to_string()).collect();

        for id in ids {
            let value = match id.as_str() {
                "cc1" | "cc2" => Box::new(matches.remove_one::<bool>(&id).unwrap()) as Box<dyn Any>,
                "scale_by_bin" => {
                    Box::new(matches.remove_many::<f64>(&id).unwrap().collect::<Vec<_>>())
                        as Box<dyn Any>
                }
                _ => unreachable!(),
            };

            stack.push((id, value));
        }

        Ok(Self { stack })
    }

    fn update_from_arg_matches(&mut self, _: &ArgMatches) -> Result<(), Error> {
        unreachable!()
    }

    fn update_from_arg_matches_mut(&mut self, _: &mut ArgMatches) -> Result<(), Error> {
        unreachable!()
    }
}

impl Args for MoreArgs {
    fn augment_args(cmd: Command) -> Command {
        cmd.arg(
            Arg::new("cc1")
                .action(ArgAction::SetTrue)
                .help("Charge conjugate the first initial state")
                .long("cc1"),
        )
        .arg(
            Arg::new("cc2")
                .action(ArgAction::SetTrue)
                .help("Charge conjugate the second initial state")
                .long("cc2"),
        )
        .arg(
            Arg::new("scale_by_bin")
                .action(ArgAction::Append)
                .help("Scale each bin with a different factor")
                .long("scale-by-bin")
                .num_args(1)
                .value_delimiter(',')
                .value_name("BIN1,BIN2,...")
                .value_parser(value_parser!(f64)),
        )
    }

    fn augment_args_for_update(_: Command) -> Command {
        unreachable!()
    }
}

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<ExitCode> {
        let mut grid = helpers::read_grid(&self.input)?;

        for (id, value) in &self.more_args.stack {
            match id.as_str() {
                "cc1" | "cc2" => {
                    // default value doesn't have an effect
                    if !value.downcast_ref::<bool>().copied().unwrap() {
                        continue;
                    }

                    let cc1 = id == "cc1";
                    let cc2 = id == "cc2";

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
                                        let (ap, f1) = if cc1 {
                                            pids::charge_conjugate(lumi_id_types, a)
                                        } else {
                                            (a, 1.0)
                                        };
                                        let (bp, f2) = if cc2 {
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

                    if cc1 {
                        initial_state_1 = pids::charge_conjugate_pdg_pid(initial_state_1);
                    }
                    if cc2 {
                        initial_state_2 = pids::charge_conjugate_pdg_pid(initial_state_2);
                    }

                    grid.set_key_value("initial_state_1", &initial_state_1.to_string());
                    grid.set_key_value("initial_state_2", &initial_state_2.to_string());
                    grid.set_lumis(lumis);
                }
                "scale_by_bin" => grid.scale_by_bin(value.downcast_ref::<Vec<_>>().unwrap()),
                _ => unreachable!(),
            }
        }

        helpers::write_grid(&self.output, &grid)
    }
}
