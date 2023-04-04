use super::helpers::{self, GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::builder::{PossibleValuesParser, TypedValueParser};
use clap::{
    value_parser, Arg, ArgAction, ArgMatches, Args, Command, Error, FromArgMatches, Parser,
    ValueHint,
};
use pineappl::fk_table::{FkAssumptions, FkTable};
use pineappl::lumi::LumiEntry;
use pineappl::pids;
use std::any::Any;
use std::fs;
use std::ops::{Deref, RangeInclusive};
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
                "cc1" | "cc2" | "optimize" | "upgrade" => {
                    Box::new(matches.remove_one::<bool>(&id).unwrap()) as Box<dyn Any>
                }
                "delete_bins" => Box::new(
                    matches
                        .remove_many::<RangeInclusive<usize>>(&id)
                        .unwrap()
                        .flatten()
                        .collect::<Vec<_>>(),
                ) as Box<dyn Any>,
                "delete_key" => {
                    Box::new(matches.remove_one::<String>(&id).unwrap()) as Box<dyn Any>
                }
                "merge_bins" => Box::new(matches.remove_one::<RangeInclusive<usize>>(&id).unwrap())
                    as Box<dyn Any>,
                "optimize_fk_table" => {
                    Box::new(matches.remove_one::<FkAssumptions>(&id).unwrap()) as Box<dyn Any>
                }
                "scale" => Box::new(matches.remove_one::<f64>(&id).unwrap()) as Box<dyn Any>,
                "scale_by_bin" | "scale_by_order" => {
                    Box::new(matches.remove_many::<f64>(&id).unwrap().collect::<Vec<_>>())
                        as Box<dyn Any>
                }
                "set_key_value" | "set_key_file" => Box::new(
                    matches
                        .remove_many::<String>(&id)
                        .unwrap()
                        .collect::<Vec<_>>(),
                ) as Box<dyn Any>,
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
            Arg::new("delete_bins")
                .action(ArgAction::Append)
                .help("Delete bins with the specified indices")
                .long("delete-bins")
                .num_args(1)
                .value_delimiter(',')
                .value_name("BIN1-BIN2,...")
                .value_parser(helpers::parse_integer_range),
        )
        .arg(
            Arg::new("delete_key")
                .action(ArgAction::Set)
                .help("Delete an internal key-value pair")
                .long("delete-key")
                .value_name("KEY"),
        )
        .arg(
            Arg::new("merge_bins")
                .action(ArgAction::Set)
                .help("Merge specific bins together")
                .long("merge-bins")
                .value_name("BIN1-BIN2")
                .value_parser(helpers::parse_integer_range),
        )
        .arg(
            Arg::new("optimize")
                .action(ArgAction::SetTrue)
                .help("Optimize internal data structure to minimize memory and disk usage")
                .long("optimize"),
        )
        .arg(
            Arg::new("optimize_fk_table")
                .action(ArgAction::Append)
                .help("Optimize internal data structure of an FkTable to minimize memory and disk usage")
                .long("optimize-fk-table")
                .num_args(1)
                .value_name("OPTIMI")
                .value_parser(
                    PossibleValuesParser::new([
                        "Nf6Ind", "Nf6Sym", "Nf5Ind", "Nf5Sym", "Nf4Ind", "Nf4Sym", "Nf3Ind",
                        "Nf3Sym",
                    ])
                    .try_map(|s| s.parse::<FkAssumptions>()),
                ),
        )
        .arg(
            Arg::new("scale")
                .action(ArgAction::Set)
                .help("Scales all grids with the given factor")
                .long("scale")
                .short('s')
                .value_name("SCALE")
                .value_parser(value_parser!(f64)),
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
        .arg(
            Arg::new("scale_by_order")
                .action(ArgAction::Append)
                .help("Scales all grids with order-dependent factors")
                .long("scale-by-order")
                .num_args(1)
                .value_delimiter(',')
                .value_name("AS,AL,LR,LF")
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("set_key_value")
                .action(ArgAction::Append)
                .allow_hyphen_values(true)
                .help("Set an internal key-value pair")
                .long("set-key-value")
                .num_args(2)
                .value_names(["KEY", "VALUE"]),
        )
        .arg(
            Arg::new("set_key_file")
                .action(ArgAction::Append)
                .allow_hyphen_values(true)
                .help("Set an internal key-value pair, with value being read from a file")
                .long("set-key-file")
                .num_args(2)
                .value_names(["KEY", "FILE"]),
        )
        .arg(
            Arg::new("upgrade")
                .action(ArgAction::SetTrue)
                .help("Convert the file format to the most recent version")
                .long("upgrade"),
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
                "delete_bins" => grid.delete_bins(&value.downcast_ref::<Vec<_>>().unwrap()),
                "delete_key" => {
                    grid.key_values_mut()
                        .remove(value.downcast_ref::<String>().unwrap());
                }
                "merge_bins" => {
                    let range = value.downcast_ref::<RangeInclusive<usize>>().unwrap();
                    grid.merge_bins(*range.start()..(range.end() + 1))?;
                }
                "scale" => grid.scale(value.downcast_ref().copied().unwrap()),
                "optimize" => {
                    if !value.downcast_ref::<bool>().copied().unwrap() {
                        continue;
                    }

                    grid.optimize();
                }
                "optimize_fk_table" => {
                    let assumptions = value.downcast_ref::<FkAssumptions>().copied().unwrap();
                    let mut fk_table = FkTable::try_from(grid)?;
                    fk_table.optimize(assumptions);
                    grid = fk_table.into_grid();
                }
                "scale_by_bin" => grid.scale_by_bin(value.downcast_ref::<Vec<_>>().unwrap()),
                "scale_by_order" => {
                    let scale_by_order = value.downcast_ref::<Vec<_>>().unwrap();
                    grid.scale_by_order(
                        scale_by_order[0],
                        scale_by_order[1],
                        scale_by_order[2],
                        scale_by_order[3],
                        1.0,
                    );
                }
                "set_key_value" => {
                    let key_value = value.downcast_ref::<Vec<String>>().unwrap();
                    grid.set_key_value(&key_value[0], &key_value[1]);
                }
                "set_key_file" => {
                    let key_file = value.downcast_ref::<Vec<String>>().unwrap();
                    grid.set_key_value(&key_file[0], &fs::read_to_string(&key_file[1])?);
                }
                "upgrade" => {
                    if !value.downcast_ref::<bool>().copied().unwrap() {
                        continue;
                    }

                    grid.upgrade();
                }
                _ => unreachable!(),
            }
        }

        helpers::write_grid(&self.output, &grid)
    }
}
