use anyhow::Result;
use std::env;
use std::path::Path;

fn main() -> Result<()> {
    let args: Vec<_> = env::args().skip(1).collect();

    match args.as_slice() {
        [cmd, path] if cmd == "install-manpages" => install_manpages(Path::new(path)),
        _ => todo!(),
    }
}

fn install_manpages(path: &Path) -> Result<()> {
    use clap::CommandFactory;
    use clap_mangen::Man;
    use std::fs::File;
    use std::io::BufWriter;

    fn render_manpages(path: &Path, cmd: &clap::Command) -> Result<()> {
        let name = cmd
            .get_bin_name()
            .unwrap_or(cmd.get_name())
            .replace(" ", "-");
        let version = cmd
            .get_version()
            .map(|version| version.strip_prefix('v').unwrap().to_string());

        // set the correct name of the command
        let mut man_cmd = cmd.clone().name(name.clone());

        // if there is a version, remove the 'v'
        if let Some(version) = version {
            man_cmd = man_cmd.version(version);
        }

        Man::new(man_cmd).render(&mut BufWriter::new(File::create(
            path.join(format!("{name}.1")),
        )?))?;

        for subcmd in cmd.get_subcommands() {
            render_manpages(path, subcmd)?;
        }

        Ok(())
    }

    let mut cmd = pineappl_cli::Opts::command();
    // this is needed so subcommands return the correct `bin_name`
    cmd.build();

    render_manpages(path, &cmd)?;

    Ok(())
}
