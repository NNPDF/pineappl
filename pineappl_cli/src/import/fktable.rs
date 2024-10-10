use anyhow::{anyhow, Context, Result};
use flate2::read::GzDecoder;
use ndarray::s;
use pineappl::boc::{Channel, Kinematics, Order, ScaleFuncForm, Scales};
use pineappl::channel;
use pineappl::convolutions::Convolution;
use pineappl::grid::Grid;
use pineappl::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use pineappl::packed_array::PackedArray;
use pineappl::packed_subgrid::PackedQ1X2SubgridV1;
use pineappl::pids::PidBasis;
use pineappl::subgrid::NodeValues;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::iter;
use std::path::Path;
use tar::Archive;

#[derive(Debug, PartialEq)]
enum FkTableSection {
    Sof,
    GridDesc,
    VersionInfo,
    GridInfo,
    FlavourMap,
    TheoryInfo,
    Xgrid,
    FastKernel,
}

fn read_fktable(reader: impl BufRead) -> Result<Grid> {
    let mut section = FkTableSection::Sof;
    let mut flavor_mask = Vec::<bool>::new();
    let mut x_grid = Vec::new();
    let mut grid = None;
    let mut arrays = Vec::new();
    let mut last_bin = 0;

    let mut hadronic = false;
    let mut ndata: u16 = 0;
    let mut nx1 = 0;
    let mut nx2 = 0;
    let mut q0 = 0.0;
    //let mut setname = String::new();

    for line in reader.lines() {
        let line = line?;

        match line.as_str() {
            "{GridDesc___________________________________________________" => {
                assert_eq!(section, FkTableSection::Sof);
                section = FkTableSection::GridDesc;
            }
            "_VersionInfo________________________________________________" => {
                assert_eq!(section, FkTableSection::GridDesc);
                section = FkTableSection::VersionInfo;
            }
            "_GridInfo___________________________________________________" => {
                assert_eq!(section, FkTableSection::VersionInfo);
                section = FkTableSection::GridInfo;
            }
            "{FlavourMap_________________________________________________" => {
                assert_eq!(section, FkTableSection::GridInfo);
                section = FkTableSection::FlavourMap;
            }
            "_TheoryInfo_________________________________________________" => {
                assert_eq!(section, FkTableSection::FlavourMap);
                section = FkTableSection::TheoryInfo;
            }
            "{xGrid______________________________________________________" => {
                assert_eq!(section, FkTableSection::TheoryInfo);
                section = FkTableSection::Xgrid;
            }
            "{FastKernel_________________________________________________" => {
                assert_eq!(section, FkTableSection::Xgrid);
                assert_eq!(nx1, x_grid.len());
                section = FkTableSection::FastKernel;

                nx2 = if hadronic { nx1 } else { 1 };

                // TODO: are FK tables always in the evolution basis?
                let basis = [
                    22, 100, 21, 200, 203, 208, 215, 224, 235, 103, 108, 115, 124, 135,
                ];
                let lumis = if hadronic {
                    flavor_mask
                        .iter()
                        .enumerate()
                        .filter(|&(_, &value)| value)
                        .map(|(index, _)| channel![basis[index / 14], basis[index % 14], 1.0])
                        .collect()
                } else {
                    flavor_mask
                        .iter()
                        .enumerate()
                        .filter(|&(_, &value)| value)
                        .map(|(index, _)| Channel::new(vec![(vec![basis[index]], 1.0)]))
                        .collect()
                };

                let convolutions = if hadronic {
                    vec![Convolution::UnpolPDF(2212); 2]
                } else {
                    vec![Convolution::UnpolPDF(2212)]
                };

                // construct `Grid`
                let fktable = Grid::new(
                    PidBasis::Evol,
                    lumis,
                    vec![Order::new(0, 0, 0, 0, 0)],
                    (0..=ndata).map(Into::into).collect(),
                    // legacy FK-tables only support unpolarized proton PDFs
                    convolutions.clone(),
                    // TODO: what are sensible parameters for FK-tables?
                    if hadronic {
                        vec![
                            Interp::new(
                                1e2,
                                1e8,
                                40,
                                3,
                                ReweightMeth::NoReweight,
                                Map::ApplGridH0,
                                InterpMeth::Lagrange,
                            ),
                            Interp::new(
                                2e-7,
                                1.0,
                                50,
                                3,
                                ReweightMeth::ApplGridX,
                                Map::ApplGridF2,
                                InterpMeth::Lagrange,
                            ),
                            Interp::new(
                                2e-7,
                                1.0,
                                50,
                                3,
                                ReweightMeth::ApplGridX,
                                Map::ApplGridF2,
                                InterpMeth::Lagrange,
                            ),
                        ]
                    } else {
                        vec![
                            Interp::new(
                                1e2,
                                1e8,
                                40,
                                3,
                                ReweightMeth::NoReweight,
                                Map::ApplGridH0,
                                InterpMeth::Lagrange,
                            ),
                            Interp::new(
                                2e-7,
                                1.0,
                                50,
                                3,
                                ReweightMeth::ApplGridX,
                                Map::ApplGridF2,
                                InterpMeth::Lagrange,
                            ),
                        ]
                    },
                    if hadronic {
                        vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2]
                    } else {
                        vec![Kinematics::Scale(0), Kinematics::X1]
                    },
                    // TODO: is this correct?
                    Scales {
                        ren: ScaleFuncForm::NoScale,
                        fac: ScaleFuncForm::Scale(0),
                        frg: ScaleFuncForm::NoScale,
                    },
                );

                grid = Some(fktable);

                arrays = iter::repeat(PackedArray::new(if hadronic {
                    vec![1, nx1, nx2]
                } else {
                    vec![1, nx1]
                }))
                .take(flavor_mask.iter().filter(|&&value| value).count())
                .collect();
            }
            _ => match section {
                FkTableSection::GridInfo => {
                    if let Some((key, value)) = line.split_once(' ') {
                        match key {
                            "*HADRONIC:" => {
                                hadronic = match value {
                                    "0" => false,
                                    "1" => true,
                                    _ => {
                                        unimplemented!("hadronic value: '{value}' is not supported")
                                    }
                                }
                            }
                            "*NDATA:" => ndata = value.parse()?,
                            "*NX:" => nx1 = value.parse()?,
                            "*SETNAME:" => { /*setname = value.to_string()*/ }
                            _ => unimplemented!("grid info key: '{key}' is not supported"),
                        }
                    }
                }
                FkTableSection::FlavourMap => {
                    flavor_mask.extend(line.split_whitespace().map(|token| match token {
                        "1" => true,
                        "0" => false,
                        _ => unimplemented!("flavor map entry: '{token}' is not supported"),
                    }));
                }
                FkTableSection::TheoryInfo => {
                    if let Some(("*Q0:", value)) = line.split_once(' ') {
                        q0 = value.parse()?;
                    }
                }
                FkTableSection::Xgrid => {
                    x_grid.push(line.parse()?);
                }
                FkTableSection::FastKernel => {
                    let tokens: Vec<_> = line.split_whitespace().collect();

                    let (bin, x1, x2) = (
                        tokens[0].parse::<usize>()?,
                        tokens[1].parse::<usize>()?,
                        if hadronic {
                            tokens[2].parse::<usize>()?
                        } else {
                            0
                        },
                    );

                    // if `bin` has changed, we assume that the subgrids in `array` are finished
                    if bin > last_bin {
                        let grid = grid.as_mut().unwrap();

                        for (subgrid, array) in grid
                            .subgrids_mut()
                            .slice_mut(s![0, last_bin, ..])
                            .iter_mut()
                            .zip(arrays.into_iter())
                        {
                            *subgrid = PackedQ1X2SubgridV1::new(
                                array,
                                if hadronic {
                                    vec![
                                        NodeValues::UseThese(vec![q0 * q0]),
                                        NodeValues::UseThese(x_grid.clone()),
                                        NodeValues::UseThese(x_grid.clone()),
                                    ]
                                } else {
                                    vec![
                                        NodeValues::UseThese(vec![q0 * q0]),
                                        NodeValues::UseThese(x_grid.clone()),
                                    ]
                                },
                            )
                            .into();
                        }

                        arrays = iter::repeat(PackedArray::new(if hadronic {
                            vec![1, nx1, nx2]
                        } else {
                            vec![1, nx1]
                        }))
                        .take(flavor_mask.iter().filter(|&&value| value).count())
                        .collect();
                        last_bin = bin;
                    }

                    // we can't handle `last_bin > bin`
                    assert_eq!(last_bin, bin);

                    let grid_values: Vec<f64> = tokens
                        .iter()
                        .skip(if hadronic { 3 } else { 2 })
                        .zip(flavor_mask.iter())
                        .filter(|&(_, &mask)| mask)
                        .map(|(string, _)| {
                            string.parse().with_context(|| {
                                format!("failed to parse floating point number from '{string}'")
                            })
                        })
                        .collect::<Result<_>>()?;

                    assert_eq!(grid_values.len(), arrays.len());

                    for (array, value) in arrays
                        .iter_mut()
                        .zip(grid_values.iter())
                        .filter(|(_, value)| **value != 0.0)
                    {
                        if hadronic {
                            array[[0, x1, x2]] = x_grid[x1] * x_grid[x2] * value;
                        } else {
                            array[[0, x1]] = x_grid[x1] * value;
                        }
                    }
                }
                _ => {}
            },
        }
    }

    assert_eq!(section, FkTableSection::FastKernel);

    let mut grid = grid.unwrap();

    for (subgrid, array) in grid
        .subgrids_mut()
        .slice_mut(s![0, last_bin, ..])
        .iter_mut()
        .zip(arrays.into_iter())
    {
        *subgrid = PackedQ1X2SubgridV1::new(
            array,
            if hadronic {
                vec![
                    NodeValues::UseThese(vec![q0 * q0]),
                    NodeValues::UseThese(x_grid.clone()),
                    NodeValues::UseThese(x_grid.clone()),
                ]
            } else {
                vec![
                    NodeValues::UseThese(vec![q0 * q0]),
                    NodeValues::UseThese(x_grid.clone()),
                ]
            },
        )
        .into();
    }

    Ok(grid)
}

pub fn convert_fktable(input: &Path) -> Result<Grid> {
    let reader = GzDecoder::new(File::open(input)?);

    let mut archive = Archive::new(reader);

    for entry in archive.entries()? {
        let file = entry.unwrap();
        let path = file.header().path().unwrap();

        if let Some(extension) = path.extension() {
            if extension == "dat" {
                return read_fktable(BufReader::new(file));
            }
        }
    }

    Err(anyhow!("no FK table entry found"))
}
