//! TODO.

use clap::Parser as _;
use pineappl_cli::{Opts, Subcommand as _};
use std::process::{ExitCode, Termination as _};

fn main() -> ExitCode {
    let opts = Opts::parse();

    if !opts.configuration.lhapdf_banner {
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
        use clap::CommandFactory as _;
        Opts::command().debug_assert();
    }
}
