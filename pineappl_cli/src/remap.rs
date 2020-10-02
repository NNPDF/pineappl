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
) -> Result<(), Box<dyn Error>> {
    let mut grid = Grid::read(BufReader::new(File::open(input)?))?;
    let remaps = remapping
        .split(';')
        .map(|string| {
            string
                .split(',')
                .map(str::parse::<f64>)
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;

    let dimensions = remaps.len();
    let bins = remaps.iter().fold(1, |len, vec| len * (vec.len() - 1));

    let mut normalizations = Vec::with_capacity(bins);
    let mut limits = Vec::with_capacity(bins * dimensions);

    for indices in remaps
        .iter()
        .map(|vec| 0..vec.len() - 1)
        .multi_cartesian_product()
    {
        let mut normalization = 1.0;
        for d in 0..dimensions {
            let left = remaps[d][indices[d]];
            let right = remaps[d][indices[d] + 1];

            limits.push((left, right));
            normalization *= right - left;
        }
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
