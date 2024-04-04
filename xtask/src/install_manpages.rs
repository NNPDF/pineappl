use super::Subcommand;
use anyhow::Result;
use clap::CommandFactory;
use clap::{Parser, ValueHint};
use clap_mangen::Man;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

/// Install PineAPPL's manpages.
#[derive(Parser)]
pub struct Opts {
    /// Path to where the directories should installed.
    #[arg(value_hint = ValueHint::DirPath)]
    directory: PathBuf,
}

fn render_manpages(path: &Path, cmd: &clap::Command, version: &str) -> Result<()> {
    let name = cmd
        .get_bin_name()
        .unwrap_or(cmd.get_name())
        .replace(" ", "-");

    Man::new(cmd.clone())
        // pass space, otherwise the ordering of the remaining arguments is incorrect
        .date(" ")
        .manual("PineAPPL CLI Manual")
        .source(format!("PineAPPL {version}"))
        .title(name.to_ascii_uppercase())
        .render(&mut BufWriter::new(File::create(
            path.join(format!("{name}.1")),
        )?))?;

    for subcmd in cmd.get_subcommands() {
        render_manpages(path, subcmd, version)?;
    }

    Ok(())
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let cmd = pineappl_cli::Opts::command();
        let version = cmd
            .get_version()
            // UNWRAP: the command must have a version
            .unwrap();

        // TODO: why does the version string not start with a 'v' on GitHub?
        let version = version.strip_prefix('v').unwrap_or(version).to_string();
        let mut cmd = cmd.version(version.clone());

        // this is needed so subcommands return the correct `bin_name`
        cmd.build();

        render_manpages(&self.directory, &cmd, &version)?;

        Ok(())
    }
}
