#![warn(clippy::all, clippy::cargo, clippy::nursery, clippy::pedantic)]

mod channels;
mod convolute;
mod delete;
mod diff;
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
use helpers::Subcommand;

#[derive(Parser)]
#[clap(
    arg_required_else_help = true,
    author,
    about,
    disable_help_subcommand = true,
    name = "pineappl",
    replace("lumis", &["obl", "--lumis"]), // TODO: this is for backwards compatibility, remove it
    version = git_version!(
        args = ["--always", "--dirty", "--long", "--tags"],
        cargo_prefix = "",
        fallback = "unknown"
    )
)]
struct Opts {
    /// Prevents LHAPDF from printing banners.
    #[clap(alias = "silence_lhapdf", long = "silence-lhapdf")]
    silence_lhapdf: bool,
    #[clap(subcommand)]
    subcommand: SubcommandEnum,
}

#[enum_dispatch(Subcommand)]
#[derive(Parser)]
enum SubcommandEnum {
    Channels(channels::Opts),
    Convolute(convolute::Opts),
    Delete(delete::Opts),
    Diff(diff::Opts),
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

fn main() -> Result<()> {
    let opts = Opts::parse();

    if opts.silence_lhapdf {
        lhapdf::set_verbosity(0);
    }

    opts.subcommand.run()
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl v0.5.0-beta.4-46-ge703521-dirty
Christopher Schwan <handgranaten-herbert@posteo.de>
Read, write, and query PineAPPL grids

USAGE:
    pineappl [OPTIONS] <SUBCOMMAND>

OPTIONS:
    -h, --help              Print help information
        --silence-lhapdf    Prevents LHAPDF from printing banners
    -V, --version           Print version information

SUBCOMMANDS:
    channels     Shows the contribution for each partonic channel
    convolute    Convolutes a PineAPPL grid with a PDF set
    diff         Compares the contents of two grids with each other
    info         Shows information about the grid
    lumi         Shows the luminosity function
    merge        Merges one or more PineAPPL grids together
    optimize     Optimizes the internal data structure to minimize memory usage
    orders       Shows the predictions for all bin for each order separately
    pdfunc       Calculates PDF uncertainties
    plot         Creates a matplotlib script plotting the contents of the grid
    pull         Calculates the pull between two different PDF sets
    remap        Modifies the bin dimensions, widths and normalizations
    set          Modifies the internal key-value storage
    subgrids     Print information about the internal subgrid types
    sum          Sums two or more bins of a grid together
    upgrade      Converts the file format to the most recent version
";

    #[test]
    #[ignore] // TODO: fix version string comparison
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .arg("--help")
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn verify_command() {
        use clap::CommandFactory;
        Opts::command().debug_assert();
    }
}
