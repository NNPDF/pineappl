use itertools::Itertools;
use pineappl::bin::BinRemapper;
use pineappl::grid::Grid;
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter};

pub fn subcommand(
    input: &str,
    output: &str,
    remapping: &str,
    norm: f64,
    ignore_obs_norm: &[usize],
) -> Result<(), Box<dyn Error>> {
    let mut grid = Grid::read(BufReader::new(File::open(input)?))?;
    let mut remaps = remapping
        .split(';')
        .map(|string| {
            string
                .split('|')
                .map(|string| {
                    string
                        .split(',')
                        .filter_map(|string| {
                            if string.is_empty() {
                                None
                            } else {
                                Some(string.parse::<f64>())
                            }
                        })
                        .collect::<Result<Vec<_>, _>>()
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;

    for vec in &mut remaps {
        for i in 1..vec.len() {
            if vec[i].is_empty() {
                if vec[i - 1].is_empty() {
                    // TODO: return an error
                    todo!();
                }

                vec[i] = vec[i - 1].clone();
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

    'outer: for indices in remaps
        .iter()
        .map(|vec| 0..vec.iter().map(|vec| vec.len() - 1).max().unwrap())
        .multi_cartesian_product()
    {
        let mut normalization = 1.0;
        for d in 0..dimensions {
            let index = if remaps[d].len() == 1 {
                0
            } else {
                indices[d - 1]
            };

            if remaps[d][index].len() <= (indices[d] + 1) {
                buffer.clear();
                continue 'outer;
            }

            let left = remaps[d][index][indices[d]];
            let right = remaps[d][index][indices[d] + 1];

            buffer.push((left, right));

            if !ignore_obs_norm.iter().any(|dim| *dim == (d + 1)) {
                normalization *= right - left;
            }
        }
        limits.append(&mut buffer);
        normalizations.push(norm * normalization);
    }

    grid.set_remapper(BinRemapper::new(normalizations, limits).unwrap())?;
    grid.write(BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)?,
    ))?;

    Ok(())
}
