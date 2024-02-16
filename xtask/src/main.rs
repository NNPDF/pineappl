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

    let cmd = pineappl_cli::Opts::command();
    let version: String = cmd
        .get_version()
        // UNWRAP: the command must have a version
        .unwrap()
        .strip_prefix('v')
        // UNWRAP: the version string must start with a 'v'
        .unwrap()
        .to_string();
    let mut cmd = cmd.version(version.clone());

    // this is needed so subcommands return the correct `bin_name`
    cmd.build();

    render_manpages(path, &cmd, &version)?;

    Ok(())
}
