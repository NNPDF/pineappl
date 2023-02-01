use super::helpers::{GlobalConfiguration, Subcommand};
use anyhow::Result;
use clap::{CommandFactory, Parser};
use clap_mangen::Man;
use std::io::Write;
use std::process::{Command, ExitCode, Stdio};

/// Display a manpage for selected subcommands.
#[derive(Parser)]
pub struct Opts {
    /// Name of the (chain of) subcommand(s) to read the manpage of.
    subcommand: Vec<String>,
}

impl Subcommand for Opts {
    fn run(&self, _: &GlobalConfiguration) -> Result<ExitCode> {
        let cmd = crate::Opts::command();
        let subcmd = if self.subcommand.is_empty() {
            cmd
        } else {
            let mut iter = self.subcommand.iter();
            let mut cmd = cmd;
            while let Some(s) = iter.next() {
                // TODO: if `unwrap` fails the arguments to `help` are wrong
                cmd = cmd.find_subcommand(s).unwrap().clone();
            }
            cmd
        };
        let man = Man::new(subcmd);
        let mut buffer = Vec::new();

        man.render(&mut buffer)?;

        // TODO: check remaining `unwrap`s

        // TODO: is this the best way? What do we do when `man` isn't available (Windows, ...)?
        let mut child = Command::new("man")
            .args(["-l", "-"])
            .stdin(Stdio::piped())
            .stdout(Stdio::inherit())
            .spawn()
            .unwrap();

        child.stdin.take().unwrap().write_all(&buffer).unwrap();
        child.wait().unwrap();

        Ok(ExitCode::SUCCESS)
    }
}
