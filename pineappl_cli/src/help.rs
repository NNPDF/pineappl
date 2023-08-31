use super::{GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::Parser;
use itertools::Itertools;
use std::iter;
use std::process::{Command, ExitCode};

/// Display a manpage for selected subcommands.
#[derive(Parser)]
pub struct Opts {
    /// Name of the (chain of) subcommand(s) to read the manpage of.
    subcommand: Vec<String>,
}

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<ExitCode> {
        // TODO: if `man` can't be found display the usual `--help` message
        Command::new("man")
            .arg(
                iter::once("pineappl")
                    .chain(self.subcommand.iter().map(AsRef::as_ref))
                    .join("-"),
            )
            .spawn()?
            .wait()?;

        // TODO: if this fails, print out how to install the manpages

        Ok(ExitCode::SUCCESS)
    }
}
