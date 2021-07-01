use super::helpers;
use anyhow::Result;
use std::fs;

pub fn subcommand(
    input: &str,
    output: &str,
    entries: &[&str],
    entries_from_file: &[&str],
    deletes: &[&str],
) -> Result<()> {
    let mut grid = helpers::read_grid(input)?;

    for key_value in entries.chunks(2) {
        grid.set_key_value(key_value[0], key_value[1]);
    }

    for key_file in entries_from_file.chunks(2) {
        grid.set_key_value(key_file[0], &fs::read_to_string(key_file[1])?);
    }

    for delete in deletes {
        grid.key_values_mut().remove(*delete);
    }

    helpers::write_grid(output, &grid)
}
