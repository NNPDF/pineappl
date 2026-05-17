#![allow(clippy::struct_excessive_bools)]
#![allow(clippy::too_many_arguments)]

//! TODO.

mod analyze;
mod channels;
mod convolve;
mod diff;
mod evolve;
mod export;
mod help;
mod helpers;
mod import;
mod merge;
mod orders;
mod plot;
mod pull;
mod read;
mod subgrids;
mod uncert;
mod write;

use anyhow::Result;
use clap::Parser;
use enum_dispatch::enum_dispatch;
use git_version::git_version;
use std::process::ExitCode;

/// TODO.
#[derive(Parser)]
pub struct GlobalConfiguration {
    /// Allow LHAPDF to print banners.
    #[arg(long)]
    pub lhapdf_banner: bool,
    /// Forces negative PDF values to zero.
    #[arg(long)]
    pub force_positive: bool,
    /// Allow extrapolation of PDFs outside their region of validity.
    #[arg(long)]
    pub allow_extrapolation: bool,
    /// Choose the PDF/FF set for the strong coupling.
    #[arg(default_value = "0", long, value_name = "IDX")]
    pub use_alphas_from: usize,
}

/// TODO.
#[enum_dispatch]
pub trait Subcommand {
    /// TODO.
    ///
    /// # Errors
    ///
    /// TODO.
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode>;
}

/// TODO.
#[enum_dispatch(Subcommand)]
#[derive(Parser)]
pub enum SubcommandEnum {
    /// TODO.
    Analyze(analyze::Opts),
    /// TODO.
    Channels(channels::Opts),
    /// TODO.
    Convolve(convolve::Opts),
    /// TODO.
    Diff(diff::Opts),
    /// TODO.
    Evolve(evolve::Opts),
    /// TODO.
    Export(export::Opts),
    /// TODO.
    Help(help::Opts),
    /// TODO.
    Import(import::Opts),
    /// TODO.
    Merge(merge::Opts),
    /// TODO.
    Orders(orders::Opts),
    /// TODO.
    Plot(plot::Opts),
    /// TODO.
    Pull(pull::Opts),
    /// TODO.
    Read(read::Opts),
    /// TODO.
    Subgrids(subgrids::Opts),
    /// TODO.
    Uncert(uncert::Opts),
    /// TODO.
    Write(write::Opts),
}

/// TODO.
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
pub struct Opts {
    /// TODO.
    #[command(flatten)]
    pub configuration: GlobalConfiguration,
    /// TODO.
    #[command(subcommand)]
    pub subcommand: SubcommandEnum,
}
