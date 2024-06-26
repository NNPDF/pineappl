#![allow(clippy::missing_errors_doc)]
#![allow(clippy::struct_excessive_bools)]
#![allow(clippy::too_many_arguments)]
#![allow(missing_docs)]

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

#[enum_dispatch]
pub trait Subcommand {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode>;
}

#[enum_dispatch(Subcommand)]
#[derive(Parser)]
pub enum SubcommandEnum {
    Analyze(analyze::Opts),
    Channels(channels::Opts),
    Convolve(convolve::Opts),
    Diff(diff::Opts),
    Evolve(evolve::Opts),
    Export(export::Opts),
    Help(help::Opts),
    Import(import::Opts),
    Merge(merge::Opts),
    Orders(orders::Opts),
    Plot(plot::Opts),
    Pull(pull::Opts),
    Read(read::Opts),
    Subgrids(subgrids::Opts),
    Uncert(uncert::Opts),
    Write(write::Opts),
}

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
    #[command(flatten)]
    pub configuration: GlobalConfiguration,
    #[command(subcommand)]
    pub subcommand: SubcommandEnum,
}
