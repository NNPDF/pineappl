#![warn(clippy::all, clippy::cargo, clippy::nursery, clippy::pedantic)]

mod analyze;
mod channels;
mod convolute;
mod delete;
mod diff;
mod evolve;
mod help;
mod helpers;
mod import;
mod info;
mod merge;
mod obl;
mod ops;
mod optimize;
mod orders;
mod pdfunc;
mod plot;
mod pull;
mod remap;
mod set;
mod subgrids;
mod sum;
mod upgrade;

use anyhow::Result;
use clap::Parser;
use enum_dispatch::enum_dispatch;
use git_version::git_version;
use helpers::{GlobalConfiguration, Subcommand};
use std::process::{ExitCode, Termination};

#[derive(Parser)]
#[command(
    arg_required_else_help = true,
    author,
    about,
    disable_help_subcommand = true,
    name = "pineappl",
    version = git_version!(
        args = ["--always", "--dirty", "--long", "--tags"],
        cargo_prefix = "",
        fallback = "unknown"
    )
)]
struct Opts {
    #[command(flatten)]
    configuration: GlobalConfiguration,
    #[command(subcommand)]
    subcommand: SubcommandEnum,
}

#[enum_dispatch(Subcommand)]
#[derive(Parser)]
enum SubcommandEnum {
    Analyze(analyze::Opts),
    Channels(channels::Opts),
    Convolute(convolute::Opts),
    Delete(delete::Opts),
    Diff(diff::Opts),
    Evolve(evolve::Opts),
    Help(help::Opts),
    Import(import::Opts),
    Info(info::Opts),
    Merge(merge::Opts),
    Obl(obl::Opts),
    Ops(ops::Opts),
    Optimize(optimize::Opts),
    Orders(orders::Opts),
    Pdfunc(pdfunc::Opts),
    Plot(plot::Opts),
    Pull(pull::Opts),
    Remap(remap::Opts),
    Set(set::Opts),
    Subgrids(subgrids::Opts),
    Sum(sum::Opts),
    Upgrade(upgrade::Opts),
}

fn main() -> ExitCode {
    let opts = Opts::parse();

    if opts.configuration.silence_lhapdf {
        lhapdf::set_verbosity(0);
    }

    match opts.subcommand.run(&opts.configuration) {
        Ok(code) => code,
        result @ Err(_) => result.report(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verify_command() {
        use clap::CommandFactory;
        Opts::command().debug_assert();
    }
}
