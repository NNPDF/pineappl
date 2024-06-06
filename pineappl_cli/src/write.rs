use super::helpers;
use super::{GlobalConfiguration, Subcommand};
use anyhow::{anyhow, Result};
use clap::builder::{PossibleValuesParser, TypedValueParser};
use clap::{
    value_parser, Arg, ArgAction, ArgMatches, Args, Command, Error, FromArgMatches, Parser,
    ValueHint,
};
use pineappl::bin::BinRemapper;
use pineappl::boc::{Channel, Order};
use pineappl::fk_table::{FkAssumptions, FkTable};
use pineappl::pids::PidBasis;
use std::fs;
use std::ops::RangeInclusive;
use std::path::PathBuf;
use std::process::ExitCode;

/// Write a grid modified by various operations.
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

#[derive(Clone)]
enum OpsArg {
    Cc1(bool),
    Cc2(bool),
    DedupChannels(i64),
    DeleteBins(Vec<RangeInclusive<usize>>),
    DeleteChannels(Vec<RangeInclusive<usize>>),
    DeleteKey(String),
    MergeBins(Vec<RangeInclusive<usize>>),
    Optimize(bool),
    OptimizeFkTable(FkAssumptions),
    Remap(String),
    RemapNorm(f64),
    RemapNormIgnore(Vec<usize>),
    RewriteChannel((usize, Channel)),
    RewriteOrder((usize, Order)),
    RotatePidBasis(PidBasis),
    Scale(f64),
    ScaleByBin(Vec<f64>),
    ScaleByOrder(Vec<f64>),
    SetKeyFile(Vec<String>),
    SetKeyValue(Vec<String>),
    SplitChannels(bool),
    Upgrade(bool),
}

struct MoreArgs {
    args: Vec<OpsArg>,
}

impl FromArgMatches for MoreArgs {
    fn from_arg_matches(_: &ArgMatches) -> Result<Self, Error> {
        unreachable!()
    }

    fn from_arg_matches_mut(matches: &mut ArgMatches) -> Result<Self, Error> {
        let mut args = Vec::new();
        let ids: Vec<_> = matches.ids().map(|id| id.as_str().to_owned()).collect();

        for id in ids {
            let indices: Vec<_> = matches.indices_of(&id).unwrap().collect();
            args.resize(indices.iter().max().unwrap() + 1, None);

            match id.as_str() {
                "cc1" | "cc2" | "optimize" | "split_channels" | "upgrade" => {
                    let arguments: Vec<Vec<_>> = matches
                        .remove_occurrences(&id)
                        .unwrap()
                        .map(Iterator::collect)
                        .collect();
                    assert_eq!(arguments.len(), indices.len());

                    for (index, arg) in indices.into_iter().zip(arguments.into_iter()) {
                        assert_eq!(arg.len(), 1);
                        args[index] = Some(match id.as_str() {
                            "cc1" => OpsArg::Cc1(arg[0]),
                            "cc2" => OpsArg::Cc2(arg[0]),
                            "optimize" => OpsArg::Optimize(arg[0]),
                            "split_channels" => OpsArg::SplitChannels(arg[0]),
                            "upgrade" => OpsArg::Upgrade(arg[0]),
                            _ => unreachable!(),
                        });
                    }
                }
                "remap_norm_ignore" => {
                    for (index, arg) in indices.into_iter().zip(
                        matches
                            .remove_occurrences(&id)
                            .unwrap()
                            .map(Iterator::collect::<Vec<_>>),
                    ) {
                        args[index] = Some(match id.as_str() {
                            "remap_norm_ignore" => OpsArg::RemapNormIgnore(arg),
                            _ => unreachable!(),
                        });
                    }
                }
                "dedup_channels" => {
                    for (index, arg) in indices.into_iter().zip(
                        matches
                            .remove_occurrences(&id)
                            .unwrap()
                            .map(Iterator::collect::<Vec<_>>),
                    ) {
                        assert_eq!(arg.len(), 1);
                        args[index] = Some(match id.as_str() {
                            "dedup_channels" => OpsArg::DedupChannels(arg[0]),
                            _ => unreachable!(),
                        });
                    }
                }
                "delete_key" | "remap" => {
                    for (index, mut arg) in indices.into_iter().zip(
                        matches
                            .remove_occurrences(&id)
                            .unwrap()
                            .map(Iterator::collect::<Vec<_>>),
                    ) {
                        assert_eq!(arg.len(), 1);
                        args[index] = Some(match id.as_str() {
                            "delete_key" => OpsArg::DeleteKey(arg.pop().unwrap()),
                            "remap" => OpsArg::Remap(arg.pop().unwrap()),
                            _ => unreachable!(),
                        });
                    }
                }
                "delete_bins" | "delete_channels" | "merge_bins" => {
                    for (index, arg) in indices.into_iter().zip(
                        matches
                            .remove_occurrences(&id)
                            .unwrap()
                            .map(Iterator::collect),
                    ) {
                        args[index] = Some(match id.as_str() {
                            "delete_bins" => OpsArg::DeleteBins(arg),
                            "delete_channels" => OpsArg::DeleteChannels(arg),
                            "merge_bins" => OpsArg::MergeBins(arg),
                            _ => unreachable!(),
                        });
                    }
                }
                "optimize_fk_table" => {
                    for (index, arg) in indices.into_iter().zip(
                        matches
                            .remove_occurrences(&id)
                            .unwrap()
                            .map(Iterator::collect::<Vec<_>>),
                    ) {
                        assert_eq!(arg.len(), 1);
                        args[index] = Some(match id.as_str() {
                            "optimize_fk_table" => OpsArg::OptimizeFkTable(arg[0]),
                            _ => unreachable!(),
                        });
                    }
                }
                "remap_norm" | "scale" => {
                    for (index, arg) in indices.into_iter().zip(
                        matches
                            .remove_occurrences(&id)
                            .unwrap()
                            .map(Iterator::collect::<Vec<_>>),
                    ) {
                        assert_eq!(arg.len(), 1);
                        args[index] = Some(match id.as_str() {
                            "remap_norm" => OpsArg::RemapNorm(arg[0]),
                            "scale" => OpsArg::Scale(arg[0]),
                            _ => unreachable!(),
                        });
                    }
                }
                "rewrite_channel" | "rewrite_order" => {
                    for (index, arg) in indices.into_iter().zip(
                        matches
                            .remove_occurrences(&id)
                            .unwrap()
                            .map(Iterator::collect::<Vec<String>>),
                    ) {
                        assert_eq!(arg.len(), 2);

                        args[index] = Some(match id.as_str() {
                            "rewrite_channel" => OpsArg::RewriteChannel((
                                str::parse(&arg[0]).unwrap(),
                                str::parse(&arg[1]).unwrap(),
                            )),
                            "rewrite_order" => OpsArg::RewriteOrder((
                                str::parse(&arg[0]).unwrap(),
                                str::parse(&arg[1]).unwrap(),
                            )),
                            _ => unreachable!(),
                        });
                    }
                }
                "rotate_pid_basis" => {
                    for (index, arg) in indices.into_iter().zip(
                        matches
                            .remove_occurrences(&id)
                            .unwrap()
                            .map(Iterator::collect::<Vec<_>>),
                    ) {
                        assert_eq!(arg.len(), 1);
                        args[index] = Some(match id.as_str() {
                            "rotate_pid_basis" => OpsArg::RotatePidBasis(arg[0]),
                            _ => unreachable!(),
                        });
                    }
                }
                "scale_by_bin" | "scale_by_order" => {
                    for (index, arg) in indices.into_iter().zip(
                        matches
                            .remove_occurrences(&id)
                            .unwrap()
                            .map(Iterator::collect::<Vec<_>>),
                    ) {
                        args[index] = Some(match id.as_str() {
                            "scale_by_bin" => OpsArg::ScaleByBin(arg),
                            "scale_by_order" => OpsArg::ScaleByOrder(arg),
                            _ => unreachable!(),
                        });
                    }
                }
                "set_key_file" | "set_key_value" => {
                    for (index, arg) in indices.into_iter().zip(
                        matches
                            .remove_occurrences(&id)
                            .unwrap()
                            .map(Iterator::collect::<Vec<_>>),
                    ) {
                        args[index] = Some(match id.as_str() {
                            "set_key_file" => OpsArg::SetKeyFile(arg),
                            "set_key_value" => OpsArg::SetKeyValue(arg),
                            _ => unreachable!(),
                        });
                    }
                }
                _ => unreachable!(),
            }
        }

        Ok(Self {
            args: args.into_iter().flatten().collect(),
        })
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
                .action(ArgAction::Append)
                .default_missing_value("true")
                .help("Charge conjugate the first initial state")
                .long("cc1")
                .num_args(0..=1)
                .require_equals(true)
                .value_name("ENABLE")
                .value_parser(clap::value_parser!(bool)),
        )
        .arg(
            Arg::new("cc2")
                .action(ArgAction::Append)
                .default_missing_value("true")
                .help("Charge conjugate the second initial state")
                .long("cc2")
                .num_args(0..=1)
                .require_equals(true)
                .value_name("ENABLE")
                .value_parser(clap::value_parser!(bool)),
        )
        .arg(
            Arg::new("dedup_channels")
                .action(ArgAction::Append)
                .default_missing_value("64")
                .help("Deduplicate channels assuming numbers differing by ULPS are the same")
                .long("dedup-channels")
                .num_args(0..=1)
                .require_equals(true)
                .value_name("ULPS")
                .value_parser(value_parser!(i64)),
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
            Arg::new("delete_channels")
                .action(ArgAction::Append)
                .help("Delete channels with the specified indices")
                .long("delete-channels")
                .num_args(1)
                .value_delimiter(',')
                .value_name("CH1-CH2,...")
                .value_parser(helpers::parse_integer_range),
        )
        .arg(
            Arg::new("delete_key")
                .action(ArgAction::Append)
                .help("Delete an internal key-value pair")
                .long("delete-key")
                .value_name("KEY"),
        )
        .arg(
            Arg::new("merge_bins")
                .action(ArgAction::Append)
                .help("Merge specific bins together")
                .long("merge-bins")
                .num_args(1)
                .value_delimiter(',')
                .value_name("BIN1-BIN2,...")
                .value_parser(helpers::parse_integer_range),
        )
        .arg(
            Arg::new("optimize")
                .action(ArgAction::Append)
                .default_missing_value("true")
                .help("Optimize internal data structure to minimize memory and disk usage")
                .long("optimize")
                .num_args(0..=1)
                .require_equals(true)
                .value_name("ENABLE")
                .value_parser(clap::value_parser!(bool)),
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
            Arg::new("remap")
                .action(ArgAction::Append)
                .help("Modify the bin dimensions and widths")
                .long("remap")
                .value_name("REMAPPING"),
        )
        .arg(
            Arg::new("remap_norm")
                .action(ArgAction::Append)
                .help("Modify the bin normalizations with a common factor")
                .long("remap-norm")
                .value_delimiter(',')
                .value_name("NORM")
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("remap_norm_ignore")
                .action(ArgAction::Append)
                .help("Modify the bin normalizations by multiplying with the bin lengths for the given dimensions")
                .long("remap-norm-ignore")
                .num_args(1)
                .value_delimiter(',')
                .value_name("DIM1,...")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("rewrite_channel")
                .action(ArgAction::Append)
                .help("Rewrite the definition of the channel with index IDX")
                .long("rewrite-channel")
                .num_args(2)
                .value_names(["IDX", "CHAN"])
        )
        .arg(
            Arg::new("rewrite_order")
                .action(ArgAction::Append)
                .help("Rewrite the definition of the order with index IDX")
                .long("rewrite-order")
                .num_args(2)
                .value_names(["IDX", "ORDER"])
        )
        .arg(
            Arg::new("rotate_pid_basis")
                .action(ArgAction::Append)
                .help("Rotate the PID basis for this grid")
                .long("rotate-pid-basis")
                .value_name("BASIS")
                .value_parser(
                    PossibleValuesParser::new(["PDG", "EVOL"])
                    .try_map(|s| s.parse::<PidBasis>()),
                ),
        )
        .arg(
            Arg::new("scale")
                .action(ArgAction::Append)
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
            Arg::new("split_channels")
                .action(ArgAction::Append)
                .default_missing_value("true")
                .help("Split the grid such that each channel contains only a single PID combination")
                .long("split-channels")
                .alias("split-lumi")
                .num_args(0..=1)
                .require_equals(true)
                .value_name("ENABLE")
                .value_parser(clap::value_parser!(bool)),
        )
        .arg(
            Arg::new("upgrade")
                .action(ArgAction::Append)
                .default_missing_value("true")
                .help("Convert the file format to the most recent version")
                .long("upgrade")
                .num_args(0..=1)
                .require_equals(true)
                .value_name("ENABLE")
                .value_parser(clap::value_parser!(bool)),
        )
    }

    fn augment_args_for_update(_: Command) -> Command {
        unreachable!()
    }
}

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<ExitCode> {
        let mut grid = helpers::read_grid(&self.input)?;

        for arg in &self.more_args.args {
            match arg {
                OpsArg::Cc1(true) | OpsArg::Cc2(true) => {
                    let cc1 = matches!(arg, OpsArg::Cc1(true));
                    let cc2 = matches!(arg, OpsArg::Cc2(true));

                    let pid_basis = grid.pid_basis();

                    for channel in grid.channels_mut() {
                        *channel = Channel::new(
                            channel
                                .entry()
                                .iter()
                                .map(|&(a, b, f)| {
                                    let (ap, f1) = if cc1 {
                                        pid_basis.charge_conjugate(a)
                                    } else {
                                        (a, 1.0)
                                    };
                                    let (bp, f2) = if cc2 {
                                        pid_basis.charge_conjugate(b)
                                    } else {
                                        (b, 1.0)
                                    };
                                    (ap, bp, f * f1 * f2)
                                })
                                .collect(),
                        );
                    }

                    if cc1 {
                        grid.set_convolution(0, grid.convolutions()[0].charge_conjugate());
                    }
                    if cc2 {
                        grid.set_convolution(1, grid.convolutions()[1].charge_conjugate());
                    }
                }
                OpsArg::DedupChannels(ulps) => {
                    grid.dedup_channels(*ulps);
                }
                OpsArg::DeleteBins(ranges) => {
                    grid.delete_bins(&ranges.iter().flat_map(Clone::clone).collect::<Vec<_>>());
                }
                OpsArg::DeleteChannels(ranges) => {
                    grid.delete_channels(&ranges.iter().flat_map(Clone::clone).collect::<Vec<_>>());
                }
                OpsArg::DeleteKey(key) => {
                    grid.key_values_mut().remove(key);
                }
                OpsArg::MergeBins(ranges) => {
                    // TODO: sort after increasing start indices
                    for range in ranges.iter().rev() {
                        grid.merge_bins(*range.start()..(range.end() + 1))?;
                    }
                }
                OpsArg::Remap(remapping) => grid.set_remapper(str::parse(remapping)?)?,
                OpsArg::RemapNorm(factor) => {
                    let remapper = grid
                        .remapper()
                        .ok_or_else(|| anyhow!("grid does not have a remapper"))?;
                    let normalizations = remapper
                        .normalizations()
                        .iter()
                        .copied()
                        .map(|value| factor * value)
                        .collect();

                    grid.set_remapper(
                        BinRemapper::new(normalizations, remapper.limits().to_vec()).unwrap(),
                    )?;
                }
                OpsArg::RemapNormIgnore(dimensions) => {
                    let remapper = grid
                        .remapper()
                        .ok_or_else(|| anyhow!("grid does not have a remapper"))?;
                    let normalizations = remapper
                        .limits()
                        .chunks_exact(remapper.dimensions())
                        .zip(remapper.normalizations())
                        .map(|(limits, normalization)| {
                            normalization
                                / dimensions
                                    .iter()
                                    .map(|&index| limits[index].1 - limits[index].0)
                                    .product::<f64>()
                        })
                        .collect();

                    grid.set_remapper(
                        BinRemapper::new(normalizations, remapper.limits().to_vec()).unwrap(),
                    )?;
                }
                OpsArg::RewriteChannel((index, new_channel)) => {
                    // TODO: check that `index` is valid
                    grid.channels_mut()[*index] = new_channel.clone();
                }
                OpsArg::RewriteOrder((index, order)) => {
                    grid.orders_mut()[*index] = order.clone();
                }
                OpsArg::RotatePidBasis(pid_basis) => {
                    grid.rotate_pid_basis(*pid_basis);
                }
                OpsArg::Scale(factor) => grid.scale(*factor),
                OpsArg::Optimize(true) => grid.optimize(),
                OpsArg::OptimizeFkTable(assumptions) => {
                    let mut fk_table = FkTable::try_from(grid)?;
                    fk_table.optimize(*assumptions);
                    grid = fk_table.into_grid();
                }
                OpsArg::ScaleByBin(factors) => grid.scale_by_bin(factors),
                OpsArg::ScaleByOrder(factors) => {
                    grid.scale_by_order(factors[0], factors[1], factors[2], factors[3], 1.0);
                }
                OpsArg::SetKeyValue(key_value) => {
                    grid.set_key_value(&key_value[0], &key_value[1]);
                }
                OpsArg::SetKeyFile(key_file) => {
                    grid.set_key_value(&key_file[0], &fs::read_to_string(&key_file[1])?);
                }
                OpsArg::SplitChannels(true) => grid.split_channels(),
                OpsArg::Upgrade(true) => grid.upgrade(),
                OpsArg::Cc1(false)
                | OpsArg::Cc2(false)
                | OpsArg::Optimize(false)
                | OpsArg::SplitChannels(false)
                | OpsArg::Upgrade(false) => {}
            }
        }

        helpers::write_grid(&self.output, &grid)
    }
}
