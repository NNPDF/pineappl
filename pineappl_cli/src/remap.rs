use super::helpers::{self, Subcommand};
use anyhow::{bail, ensure, Context, Result};
use clap::{Parser, ValueHint};
use itertools::izip;
use itertools::Itertools;
use pineappl::bin::BinRemapper;
use std::path::PathBuf;

/// Modifies the bin dimensions, widths and normalizations.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path of the modified PineAPPL file.
    #[clap(value_parser, value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// Remapping string. See <https://n3pdf.github.io/pineappl/docs/cli-reference.html> for full
    /// reference.
    remapping: String,
    /// Ignore the given observables for differential normalization.
    #[clap(
        alias = "ignore_obs_norm",
        long = "ignore-obs-norm",
        use_value_delimiter = true
    )]
    ignore_obs_norm: Vec<usize>,
    /// Normalization factor in addition to the given bin widths.
    #[clap(
        default_value = "1.0",
        long,
        validator = helpers::validate_pos_non_zero::<f64>
    )]
    norm: f64,
}

fn parse_remapping_string(
    remapping: &str,
    ignore_obs_norm: &[usize],
    norm: f64,
) -> Result<(Vec<f64>, Vec<(f64, f64)>)> {
    let remaps: Result<Vec<Vec<Vec<_>>>> = remapping
        .split(';')
        .map(|string| {
            string
                .split('|')
                .map(|string| {
                    string
                        .find(':')
                        .map_or(string, |index| {
                            let (lhs, rhs) = string.split_at(index);
                            let rhs = &rhs[1..]; // remove ':' which is contained with `split_at`

                            // extract the part that doesn't belong to the ':' specification
                            match (lhs.trim().parse::<usize>(), rhs.trim().parse::<usize>()) {
                                (Err(_), Ok(_)) => lhs,
                                (Ok(_), Err(_)) => rhs,
                                _ => "",
                            }
                        })
                        .split(',')
                        .filter_map(|string| {
                            let string = string.trim();
                            if string.is_empty() {
                                None
                            } else {
                                Some(
                                    string
                                        .parse::<f64>()
                                        .context(format!("unable to parse limit '{string}'")),
                                )
                            }
                        })
                        .collect()
                })
                .collect()
        })
        .collect();
    let mut remaps = remaps?;

    ensure!(
        remaps[0].len() == 1,
        "'|' syntax not meaningful for first dimension"
    );

    // go over `remaps` again, and repeat previous entries as requested with the `|` syntax
    for vec in &mut remaps {
        for i in 1..vec.len() {
            if vec[i].is_empty() {
                ensure!(!vec[i - 1].is_empty(), "empty repetition with '|'");
                vec[i] = vec[i - 1].clone();
            }
        }
    }

    // go over `remaps` again, this time remove bin as requested with the `:N` or `N:` syntax
    for (vec, string) in remaps.iter_mut().zip(remapping.split(';')) {
        for (vec, string) in vec.iter_mut().zip(string.split('|')) {
            let (lhs, rhs) = {
                let split: Vec<_> = string.split(':').collect();

                if split.len() == 1 {
                    // there's no colon
                    continue;
                }

                ensure!(split.len() == 2, "too many ':' found: '{}'", string);

                (split[0].parse::<usize>(), split[1].parse::<usize>())
            };

            if let Ok(num) = rhs {
                vec.truncate(vec.len() - num);
            }

            if let Ok(num) = lhs {
                vec.drain(0..num);
            }

            if lhs.is_err() && rhs.is_err() {
                bail!("unable to parse ':' syntax from: '{}'", string);
            }

            if vec.len() <= 1 {
                bail!("no limits due to ':' syntax");
            }
        }
    }

    let dimensions = remaps.len();
    let mut normalizations = Vec::new();
    let mut limits = Vec::new();
    let mut buffer = Vec::with_capacity(dimensions);
    let mut pipe_indices = vec![0; dimensions];
    let mut last_indices = vec![0; dimensions];

    'looop: for indices in remaps
        .iter()
        .map(|vec| 0..vec.iter().map(|vec| vec.len() - 1).max().unwrap())
        .multi_cartesian_product()
    {
        // calculate `pipe_indices`, which stores the indices for the second dimension of `remaps`
        for d in 0..dimensions - 1 {
            if indices[d] > last_indices[d] {
                for dp in d + 1..dimensions {
                    if remaps[dp].len() != 1 {
                        pipe_indices[dp] += 1;
                    }
                }
            }
        }

        last_indices = indices.clone();

        let mut normalization = 1.0;

        for (d, (remap, &pipe_index, &i)) in izip!(&remaps, &pipe_indices, &indices).enumerate() {
            if let Some(r) = remap.get(pipe_index) {
                if r.len() <= (i + 1) {
                    buffer.clear();

                    // this index doesn't exist
                    continue 'looop;
                }

                let left = r[i];
                let right = r[i + 1];

                buffer.push((left, right));

                if !ignore_obs_norm.iter().any(|dim| *dim == (d + 1)) {
                    normalization *= right - left;
                }
            } else {
                bail!("missing '|' specification: number of variants too small");
            }
        }

        limits.append(&mut buffer);
        normalizations.push(norm * normalization);
    }

    Ok((normalizations, limits))
}

impl Subcommand for Opts {
    fn run(&self) -> Result<u8> {
        let mut grid = helpers::read_grid(&self.input)?;
        let (normalizations, limits) =
            parse_remapping_string(&self.remapping, &self.ignore_obs_norm, self.norm)?;
        grid.set_remapper(BinRemapper::new(normalizations, limits).unwrap())?;
        helpers::write_grid(&self.output, &grid)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "'|' syntax not meaningful for first dimension")]
    fn pipe_syntax_first_dimension() {
        parse_remapping_string("|0,1,2", &[], 1.0).unwrap();
    }

    #[test]
    #[should_panic(expected = "empty repetition with '|'")]
    fn pipe_syntax_first_empty() {
        parse_remapping_string("0,1,2;0,2,4;||", &[], 1.0).unwrap();
    }

    #[test]
    #[should_panic(expected = "too many ':' found: '::'")]
    fn colon_syntax_too_many_colons() {
        parse_remapping_string("0,1,2;0,2,4;1,2,3,4,5|::", &[], 1.0).unwrap();
    }

    #[test]
    #[should_panic(expected = "unable to parse ':' syntax from: '2.5:'")]
    fn colon_syntax_bad_lhs() {
        parse_remapping_string("0,1,2;0,2,4;1,2,3,4,5|2.5:|:3|:3", &[], 1.0).unwrap();
    }

    #[test]
    #[should_panic(expected = "unable to parse ':' syntax from: ':2.5'")]
    fn colon_syntax_bad_rhs() {
        parse_remapping_string("0,1,2;0,2,4;1,2,3,4,5|:2.5|:3|:3", &[], 1.0).unwrap();
    }

    #[test]
    #[should_panic(expected = "no limits due to ':' syntax")]
    fn colon_syntax_no_limits() {
        parse_remapping_string("0,1,2;0,2,4;1,2,3,4,5|:4|:3|:3", &[], 1.0).unwrap();
    }

    #[test]
    #[should_panic(expected = "missing '|' specification: number of variants too small")]
    fn pipe_syntax_too_few_pipes() {
        parse_remapping_string("0,1,2;0,2,4;1,2,3|4,5,6|7,8,9", &[], 1.0).unwrap();
    }
}
