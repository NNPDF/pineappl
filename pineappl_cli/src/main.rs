#![allow(missing_docs)]

use clap::Parser;
use pineappl_cli::{Opts, Subcommand};
use std::process::{ExitCode, Termination};

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
