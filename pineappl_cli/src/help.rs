use super::helpers::{GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::{CommandFactory, Parser};
use clap_mangen::Man;
use is_terminal::IsTerminal;
use std::io;
use std::process::{Command, ExitCode, Stdio};

fn find_subcommand(mut cmd: clap::Command, args: &[String]) -> Option<clap::Command> {
    if !args.is_empty() {
        for s in args.iter() {
            cmd = cmd.find_subcommand(s)?.clone();
        }
    }
    Some(cmd)
}

/// Display a manpage for selected subcommands.
#[derive(Parser)]
pub struct Opts {
    /// Name of the (chain of) subcommand(s) to read the manpage of.
    subcommand: Vec<String>,
}

// TODO: installing manpages requires <https://github.com/rust-lang/cargo/issues/2729>, but we
// could go ahead and install them ourselves to /path/of/pineappl/../share/man/man1/pineappl.1.gz

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<ExitCode> {
        // TODO: SYNOPSIS doesn't show top-level subcommands: instead of `pineappl-convolute` only
        // `convolute` is shown

        // TODO: VERSION shows the version with two `v` in front of it

        let mut cmd = crate::Opts::command();
        cmd.build();

        // TODO: if `unwrap` fails the arguments were wrong
        let man = Man::new(find_subcommand(cmd, &self.subcommand).unwrap());

        if io::stdout().is_terminal() {
            // TODO: if `man` can't be found display the usual `--help` message
            let mut child = Command::new("man")
                .args(["/dev/stdin"])
                .stdin(Stdio::piped())
                .stdout(Stdio::inherit())
                .spawn()?;

            man.render(&mut child.stdin.take().unwrap())?;
            child.wait()?;
        } else {
            man.render(&mut io::stdout())?;
        }

        Ok(ExitCode::SUCCESS)
    }
}
