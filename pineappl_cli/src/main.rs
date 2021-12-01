#![warn(clippy::all, clippy::cargo, clippy::nursery, clippy::pedantic)]

mod channels;
mod convolute;
mod diff;
mod helpers;
mod info;
mod luminosity;
mod merge;
mod optimize;
mod orders;
mod pdf_uncertainty;
mod plot;
mod pull;
mod remap;
mod set;
mod subgrids;
mod sum;
mod upgrade;

use anyhow::Result;
use clap::{AppSettings, Parser};
use git_version::git_version;

#[derive(Parser)]
#[clap(
    author,
    about,
    name = "pineappl",
    setting(AppSettings::DisableHelpSubcommand),
    setting(AppSettings::SubcommandRequiredElseHelp),
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
    subcommand: Subcommands,
}

#[derive(Parser)]
enum Subcommands {
    Channels(channels::Opts),
    Convolute(convolute::Opts),
    Diff(diff::Opts),
    Info(info::Opts),
    Luminosity(luminosity::Opts),
    Merge(merge::Opts),
    Optimize(optimize::Opts),
    Orders(orders::Opts),
    PdfUncertainty(pdf_uncertainty::Opts),
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

    match opts.subcommand {
        Subcommands::Channels(opts) => opts.subcommand()?.printstd(),
        Subcommands::Convolute(opts) => opts.subcommand()?.printstd(),
        Subcommands::Diff(opts) => opts.subcommand()?.printstd(),
        Subcommands::Info(opts) => {
            opts.subcommand()?;
            0
        }
        Subcommands::Luminosity(opts) => opts.subcommand()?.printstd(),
        Subcommands::Merge(opts) => {
            opts.subcommand()?;
            0
        }
        Subcommands::Optimize(opts) => {
            opts.subcommand()?;
            0
        }
        Subcommands::Orders(opts) => opts.subcommand()?.printstd(),
        Subcommands::PdfUncertainty(opts) => opts.subcommand()?.printstd(),
        Subcommands::Plot(opts) => {
            opts.subcommand()?;
            0
        }
        Subcommands::Pull(opts) => opts.subcommand()?.printstd(),
        Subcommands::Remap(opts) => {
            opts.subcommand()?;
            0
        }
        Subcommands::Set(opts) => {
            opts.subcommand()?;
            0
        }
        Subcommands::Subgrids(opts) => opts.subcommand()?.printstd(),
        Subcommands::Sum(opts) => {
            opts.subcommand()?;
            0
        }
        Subcommands::Upgrade(opts) => {
            opts.subcommand()?;
            0
        }
    };

    Ok(())
}
