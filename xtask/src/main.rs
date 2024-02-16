use anyhow::Result;
use clap::Parser;
use enum_dispatch::enum_dispatch;

mod install_manpages;

#[enum_dispatch]
pub trait Subcommand {
    fn run(&self) -> Result<()>;
}

#[enum_dispatch(Subcommand)]
#[derive(Parser)]
pub enum SubcommandEnum {
    InstallManpages(install_manpages::Opts),
}

#[derive(Parser)]
#[command(about, arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Opts {
    #[command(subcommand)]
    pub subcommand: SubcommandEnum,
}

fn main() -> Result<()> {
    let opts = Opts::parse();

    opts.subcommand.run()
}
