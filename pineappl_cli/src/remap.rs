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
                                Some(string.parse::<f64>())
                            }
                        })
                        .collect::<Result<Vec<_>, _>>()
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;

    if remaps[0].len() != 1 {
        todo!();
    }

    // go over `remaps` again, and repeat previous entries as requested with the `||` syntax
    for vec in &mut remaps {
        for i in 1..vec.len() {
            if vec[i].is_empty() {
                if vec[i - 1].is_empty() {
                    todo!();
                }

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
                } else if split.len() != 2 {
                    todo!();
                } else {
                    (split[0], split[1])
                }
            };

            if let Ok(remove_from_left) = lhs.parse::<usize>() {
                vec.drain(0..remove_from_left);
            } else if let Ok(remove_from_right) = rhs.parse::<usize>() {
                vec.truncate(vec.len() - remove_from_right);
            } else {
                todo!();
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

    for indices in remaps
        .iter()
        .map(|vec| 0..vec.iter().map(|vec| vec.len() - 1).max().unwrap())
        .multi_cartesian_product()
    {
        let mut normalization = 1.0;
        for d in 0..dimensions {
            let index = indices1[d];

            if remaps[d][index].len() <= (indices[d] + 1) {
                buffer.clear();

                // this index doesn't exist
                break;
            }

            let left = remaps[d][index][indices[d]];
            let right = remaps[d][index][indices[d] + 1];

            buffer.push((left, right));

            if !ignore_obs_norm.iter().any(|dim| *dim == (d + 1)) {
                normalization *= right - left;
            }
        }

        if !buffer.is_empty() {
            limits.append(&mut buffer);
            normalizations.push(norm * normalization);

            println!("{:?} , {:?} -> {:?}", last_indices, indices, indices1);

            for d in 0..dimensions - 1 {
                if indices[d] > last_indices[d] {
                    for dp in d + 1..dimensions {
                        if remaps[dp].len() != 1 {
                            indices1[dp] += 1;
                        }
                    }
                }
            }

            last_indices = indices;
        }
    }

    normalizations.shrink_to_fit();
    limits.shrink_to_fit();

    grid.set_remapper(BinRemapper::new(normalizations, limits).unwrap())?;
    grid.write(BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)?,
    ))?;

    Ok(())
}
