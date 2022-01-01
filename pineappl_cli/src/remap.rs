use super::helpers::{self, Subcommand};
use anyhow::{bail, ensure, Context, Result};
use clap::{Parser, ValueHint};
use itertools::Itertools;
use pineappl::bin::BinRemapper;
use std::path::PathBuf;

/// Modifies the bin dimensions, widths and normalizations.
#[derive(Parser)]
#[clap(
    after_help = "For performance/simplicity reasons Monte Carlo programs (and the PineAPPL
`Grid::fill` method) typically do not support 1) multi-dimensional distributions or 2)
distributions whose bin sizes are not equally sized *during generation*. To work around this
problem a grid with a one-dimensional distribution can be generated instead, and afterwards the
bins can be 'remapped' to an N-dimensional distribution using the limits specified with the
REMAPPING string.

The remapping string uses the following special characters to achieve this (note that some of these
characters must be escaped in certain shells):

- ',': In the easiest case, when mapping bins to 1-dimensional bins, the bin limits can be given by
separating them using ',' as in the following example: '0,0.2,0.4,0.6,1'. This remapping string
will expect the grid to have 4 bins which cover the range from 0 to 0.6 with 3 bins of equal size
and the range 0.6 to 1 with another (single) bin
- ';': If higher-dimensional bins are needed, the bin limits are constructed using a cartesian
product where the bin limits for each dimension is separated using a semicolon, for example:
'0,0.5,1;0,1,2' expects the grid to have 4 bins, whose 2-dimensional bin limits are: (0-0.5;0-1),
(0-0.5;1-2), (00.5-1;0-1), (0.5-1;1-2)
- '|': When constructing bin limits for higher-dimensional distributions, it is sometimes necessary
to have bin limits that are dependent on other dimensions; the pipe '|' can be used to achieve
this. For example the remapping string '0,1,2;0,1,2|0,2,4' expects a grid with a total of 4 bins
and produces the following bin limits: (0-1;0-1), (0-1;1-2), (1-2;0-2), (1-2;2-4). Leaving out the
bin limits means that the previous limits are repeated; for example, the following remapping string
for six bins '0,1,2;0,1,2|0,1,2|0,2,4' can be more succinctly written as '0,1,2;0,1,2||0,2,4'
- ':': Finally, the colon can be used to 'cut' out bins from the left and right for
dimension-dependent bins. For example, the remapping string '0,1,2;0,1,2,3:2|:1||:1|2:' is a more
succinct way of writing the following remapping string: '0,1,2;0,1|0,1,2|0,1,2,3|0,1,2|2,3'

Finally note that the differential cross sections are calculated using the bin sizes (the product
of bin widths of each dimension) given by the remapping string. The option `--ignore-obs-norm` can
be used to remove certain dimensions from the bin size determination, for example `'0,10,20;0,2,4'
--ignore-obs-norm 0` will normalize the bins with a size of `2` because the first dimension (with
index `0` will be ignored)
"
)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path of the modified PineAPPL file.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// Remapping string.
    remapping: String,
    /// Ignore the given observables for differential normalization.
    #[clap(
        alias = "ignore_obs_norm",
        long = "ignore-obs-norm",
        use_delimiter = true
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

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let mut grid = helpers::read_grid(&self.input)?;
        let remaps: Result<Vec<Vec<Vec<_>>>> = self
            .remapping
            .split(';')
            .map(|string| {
                string
                    .split('|')
                    .map(|string| {
                        string
                            .find(':')
                            .map_or(string, |index| {
                                let (lhs, rhs) = string.split_at(index);
                                if lhs.trim().parse::<usize>().is_err() {
                                    lhs
                                } else {
                                    rhs
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
                                            .context(format!("unable to parse limit '{}'", string)),
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

        // go over `remaps` again, and repeat previous entries as requested with the `||` syntax
        for vec in &mut remaps {
            for i in 1..vec.len() {
                if vec[i].is_empty() {
                    ensure!(!vec[i - 1].is_empty(), "empty repetition with '|'");
                    vec[i] = vec[i - 1].clone();
                }
            }
        }

        // go over `remaps` again, this time remove bin as requested with the `:N` or `N:` syntax
        for (vec, string) in remaps.iter_mut().zip(self.remapping.split(';')) {
            for (vec, string) in vec.iter_mut().zip(string.split('|')) {
                let (lhs, rhs) = {
                    let split: Vec<_> = string.split(':').collect();

                    if split.len() == 1 {
                        // there's no colon
                        continue;
                    }

                    ensure!(split.len() == 2, "too many ':' found: '{}'", string);

                    (split[0], split[1])
                };

                let lhs = lhs.parse::<usize>();
                let rhs = rhs.parse::<usize>();

                if let Ok(remove_from_left) = lhs {
                    ensure!(
                        rhs.is_err(),
                        "ambiguity in parsing ':' syntax from: '{}'",
                        string
                    );
                    vec.drain(0..remove_from_left);
                } else if let Ok(remove_from_right) = rhs {
                    ensure!(
                        lhs.is_err(),
                        "ambiguity in parsing ':' syntax from: '{}'",
                        string
                    );
                    vec.truncate(vec.len() - remove_from_right);
                } else {
                    bail!("unable to parse ':' syntax from: '{}'", string);
                }
            }
        }

        let dimensions = remaps.len();
        let max_bins = remaps.iter().fold(1, |bins, vec| {
            bins * vec.iter().map(|vec| vec.len() - 1).max().unwrap()
        });

        let mut normalizations = Vec::with_capacity(max_bins);
        let mut limits = Vec::with_capacity(max_bins * dimensions);
        let mut buffer = Vec::with_capacity(dimensions);

        let mut indices1 = vec![0; dimensions];
        let mut last_indices = vec![0; dimensions];

        'looop: for indices in remaps
            .iter()
            .map(|vec| 0..vec.iter().map(|vec| vec.len() - 1).max().unwrap())
            .multi_cartesian_product()
        {
            for d in 0..dimensions - 1 {
                if indices[d] > last_indices[d] {
                    for dp in d + 1..dimensions {
                        if remaps[dp].len() != 1 {
                            indices1[dp] += 1;
                        }
                    }
                }
            }

            last_indices = indices.clone();

            let mut normalization = 1.0;
            for d in 0..dimensions {
                let index = indices1[d];

                if remaps[d][index].len() <= (indices[d] + 1) {
                    buffer.clear();

                    // this index doesn't exist
                    continue 'looop;
                }

                let left = remaps[d][index][indices[d]];
                let right = remaps[d][index][indices[d] + 1];

                buffer.push((left, right));

                if !self.ignore_obs_norm.iter().any(|dim| *dim == (d + 1)) {
                    normalization *= right - left;
                }
            }

            limits.append(&mut buffer);
            normalizations.push(self.norm * normalization);
        }

        normalizations.shrink_to_fit();
        limits.shrink_to_fit();

        grid.set_remapper(BinRemapper::new(normalizations, limits).unwrap())?;
        helpers::write_grid(&self.output, &grid)
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;
    use assert_fs::NamedTempFile;

    const HELP_STR: &str = "pineappl-remap 
Modifies the bin dimensions, widths and normalizations

USAGE:
    pineappl remap [OPTIONS] <INPUT> <OUTPUT> <REMAPPING>

ARGS:
    <INPUT>        Path to the input grid
    <OUTPUT>       Path of the modified PineAPPL file
    <REMAPPING>    Remapping string

OPTIONS:
    -h, --help
            Print help information

        --ignore-obs-norm <IGNORE_OBS_NORM>
            Ignore the given observables for differential normalization

        --norm <NORM>
            Normalization factor in addition to the given bin widths [default: 1.0]

For performance/simplicity reasons Monte Carlo programs (and the PineAPPL
`Grid::fill` method) typically do not support 1) multi-dimensional distributions or 2)
distributions whose bin sizes are not equally sized *during generation*. To work around this
problem a grid with a one-dimensional distribution can be generated instead, and afterwards the
bins can be 'remapped' to an N-dimensional distribution using the limits specified with the
REMAPPING string.

The remapping string uses the following special characters to achieve this (note that some of these
characters must be escaped in certain shells):

- ',': In the easiest case, when mapping bins to 1-dimensional bins, the bin limits can be given by
separating them using ',' as in the following example: '0,0.2,0.4,0.6,1'. This remapping string
will expect the grid to have 4 bins which cover the range from 0 to 0.6 with 3 bins of equal size
and the range 0.6 to 1 with another (single) bin
- ';': If higher-dimensional bins are needed, the bin limits are constructed using a cartesian
product where the bin limits for each dimension is separated using a semicolon, for example:
'0,0.5,1;0,1,2' expects the grid to have 4 bins, whose 2-dimensional bin limits are: (0-0.5;0-1),
(0-0.5;1-2), (00.5-1;0-1), (0.5-1;1-2)
- '|': When constructing bin limits for higher-dimensional distributions, it is sometimes necessary
to have bin limits that are dependent on other dimensions; the pipe '|' can be used to achieve
this. For example the remapping string '0,1,2;0,1,2|0,2,4' expects a grid with a total of 4 bins
and produces the following bin limits: (0-1;0-1), (0-1;1-2), (1-2;0-2), (1-2;2-4). Leaving out the
bin limits means that the previous limits are repeated; for example, the following remapping string
for six bins '0,1,2;0,1,2|0,1,2|0,2,4' can be more succinctly written as '0,1,2;0,1,2||0,2,4'
- ':': Finally, the colon can be used to 'cut' out bins from the left and right for
dimension-dependent bins. For example, the remapping string '0,1,2;0,1,2,3:2|:1||:1|2:' is a more
succinct way of writing the following remapping string: '0,1,2;0,1|0,1,2|0,1,2,3|0,1,2|2,3'

Finally note that the differential cross sections are calculated using the bin sizes (the product
of bin widths of each dimension) given by the remapping string. The option `--ignore-obs-norm` can
be used to remove certain dimensions from the bin size determination, for example `'0,10,20;0,2,4'
--ignore-obs-norm 0` will normalize the bins with a size of `2` because the first dimension (with
index `0` will be ignored)
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["remap", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    #[ignore]
    fn default() {
        let output = NamedTempFile::new("optimized.pineappl.lz4").unwrap();

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "remap",
                "data/LHCB_WP_7TEV.pineappl.lz4",
                output.path().to_str().unwrap(),
                "0,1,2;0,2,4,6,8",
            ])
            .assert()
            .success()
            .stdout("");

        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "--silence-lhapdf",
                "convolute",
                output.path().to_str().unwrap(),
                "NNPDF31_nlo_as_0118_luxqed",
            ])
            .assert()
            .success()
            .stdout("");
    }
}
