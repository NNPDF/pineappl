use super::bin::{BinLimits, BinRemapper};
use super::boc::{Channel, Kinematics, Order, ScaleFuncForm, Scales};
use super::convolutions::Convolution;
use super::empty_subgrid::EmptySubgridV1;
use super::grid::{Grid, GridError, Mmv4, MoreMembers};
use super::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use super::packed_array::PackedArray;
use super::packed_subgrid::PackedQ1X2SubgridV1;
use super::pids::PidBasis;
use super::subgrid::{Mu2, NodeValues};
use ndarray::Array3;
use pineappl_v0::grid::Grid as GridV0;
use std::io::BufRead;
use std::iter;

pub fn default_interps(convolutions: usize) -> Vec<Interp> {
    let mut interps = vec![Interp::new(
        1e2,
        1e8,
        40,
        3,
        ReweightMeth::NoReweight,
        Map::ApplGridH0,
        InterpMeth::Lagrange,
    )];

    for _ in 0..convolutions {
        interps.push(Interp::new(
            2e-7,
            1.0,
            50,
            3,
            ReweightMeth::ApplGridX,
            Map::ApplGridF2,
            InterpMeth::Lagrange,
        ));
    }

    interps
}

pub fn read_uncompressed_v0(mut reader: impl BufRead) -> Result<Grid, GridError> {
    use pineappl_v0::subgrid::Subgrid as _;

    let grid = GridV0::read(&mut reader).map_err(|err| GridError::Other(err.into()))?;
    let convolutions = read_convolutions_from_metadata(&grid);

    // TODO: read in flexible-scale grids properly
    let mut kinematics = vec![Kinematics::Scale(0), Kinematics::X(0)];
    if convolutions[0].is_some() && convolutions[1].is_some() {
        kinematics.push(Kinematics::X(1));
    }

    assert_eq!(convolutions.len(), 2);

    let result = Grid {
        subgrids: Array3::from_shape_vec(
            grid.subgrids().dim(),
            grid.subgrids()
                .into_iter()
                .map(|subgrid| {
                    if subgrid.is_empty() {
                        EmptySubgridV1.into()
                    } else {
                        let mu2_grid: Vec<_> = subgrid
                            .mu2_grid()
                            .iter()
                            .map(|mu2v0| Mu2 {
                                ren: mu2v0.ren,
                                fac: mu2v0.fac,
                                frg: -1.0,
                            })
                            .collect();

                        let mut dim = vec![mu2_grid.len()];
                        if convolutions[0].is_some() {
                            dim.push(subgrid.x1_grid().len());
                        }
                        if convolutions[1].is_some() {
                            dim.push(subgrid.x2_grid().len());
                        }
                        let mut array = PackedArray::new(dim);

                        if convolutions[0].is_none() {
                            for (index, v) in subgrid.indexed_iter() {
                                array[[index.0, index.2]] = v;
                            }
                        } else if convolutions[1].is_none() {
                            for (index, v) in subgrid.indexed_iter() {
                                array[[index.0, index.1]] = v;
                            }
                        } else {
                            for (index, v) in subgrid.indexed_iter() {
                                array[<[usize; 3]>::from(index)] = v;
                            }
                        }

                        let mut node_values = vec![NodeValues::UseThese(
                            mu2_grid.iter().map(|&Mu2 { ren, .. }| ren).collect(),
                        )];
                        if convolutions[0].is_some() {
                            node_values.push(NodeValues::UseThese(subgrid.x1_grid().into_owned()));
                        }
                        if convolutions[1].is_some() {
                            node_values.push(NodeValues::UseThese(subgrid.x2_grid().into_owned()));
                        }
                        PackedQ1X2SubgridV1::new(array, node_values).into()
                    }
                })
                .collect(),
        )
        // UNWRAP: the dimensions must be the same as in the v0 grid
        .unwrap(),
        channels: grid
            .channels()
            .iter()
            .map(|c| {
                Channel::new(
                    c.entry()
                        .iter()
                        .map(|&(a, b, f)| {
                            let mut pids = Vec::new();
                            if convolutions[0].is_some() {
                                pids.push(a);
                            }
                            if convolutions[1].is_some() {
                                pids.push(b);
                            };
                            (pids, f)
                        })
                        .collect(),
                )
            })
            .collect(),
        // TODO: change this member to something much easier to handle
        bin_limits: BinLimits::new(if grid.remapper().is_none() {
            let limits = &grid.bin_info().limits();
            iter::once(limits[0][0].0)
                .chain(limits.iter().map(|v| v[0].1))
                .collect()
        } else {
            // if there is a BinRemapper this member will likely have no impact
            (0..=grid.bin_info().bins())
                .map(|i| f64::from(u16::try_from(i).unwrap()))
                .collect()
        }),
        orders: grid
            .orders()
            .iter()
            .map(|o| Order {
                alphas: o.alphas,
                alpha: o.alpha,
                logxir: o.logxir,
                logxif: o.logxif,
                logxia: 0,
            })
            .collect(),
        metadata: grid
            .key_values()
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .collect(),
        interps: default_interps(convolutions.len()),
        convolutions: convolutions
            .into_iter()
            //.map(|conv| conv.unwrap_or(Convolution::None))
            .filter_map(|conv| conv)
            .collect(),
        pid_basis: grid
            .key_values()
            .and_then(|kv| kv.get("lumi_id_types"))
            .map_or(PidBasis::Pdg, |lumi_id_types| {
                match lumi_id_types.as_str() {
                    "pdg_mc_ids" => PidBasis::Pdg,
                    "evol" => PidBasis::Evol,
                    _ => panic!("unknown PID basis '{lumi_id_types}'"),
                }
            }),
        more_members: MoreMembers::V4(Mmv4),
        remapper: grid.remapper().map(|r| {
            // UNWRAP: if the old grid could be constructed with the given normalizations
            // and limits we should be able to do the same without error
            BinRemapper::new(r.normalizations().to_vec(), r.limits().to_vec()).unwrap()
        }),
        kinematics,
        scales: Scales {
            // TODO: read in flexible-scale grids properly
            ren: ScaleFuncForm::Scale(0),
            fac: ScaleFuncForm::Scale(0),
            frg: ScaleFuncForm::NoScale,
        },
    };

    assert_eq!(result.bin_info().bins(), grid.bin_info().bins());

    Ok(result)
}

fn read_convolutions_from_metadata(grid: &GridV0) -> Vec<Option<Convolution>> {
    grid.key_values().map_or_else(
        // if there isn't any metadata, we assume two unpolarized proton-PDFs are used
        || vec![Some(Convolution::UnpolPDF(2212)); 2],
        |kv| {
            // file format v0 only supports exactly two convolutions
            (1..=2)
                .map(|index| {
                    // if there are key-value pairs `convolution_particle_1` and
                    // `convolution_type_1` and the same with a higher index, we convert this
                    // metadata into `Convolution`
                    match (
                        kv.get(&format!("convolution_particle_{index}"))
                            .map(|s| s.parse::<i32>()),
                        kv.get(&format!("convolution_type_{index}"))
                            .map(String::as_str),
                    ) {
                        (_, Some("None")) => None,
                        (Some(Ok(pid)), Some("UnpolPDF")) => Some(Convolution::UnpolPDF(pid)),
                        (Some(Ok(pid)), Some("PolPDF")) => Some(Convolution::PolPDF(pid)),
                        (Some(Ok(pid)), Some("UnpolFF")) => Some(Convolution::UnpolFF(pid)),
                        (Some(Ok(pid)), Some("PolFF")) => Some(Convolution::PolFF(pid)),
                        (None, None) => {
                            // if these key-value pairs are missing use the old metadata
                            match kv
                                .get(&format!("initial_state_{index}"))
                                .map(|s| s.parse::<i32>())
                            {
                                Some(Ok(pid)) => {
                                    let condition = !grid.channels().iter().all(|entry| {
                                        entry.entry().iter().all(|&(a, b, _)|
                                            match index {
                                                1 => a,
                                                2 => b,
                                                _ => unreachable!()
                                            } == pid
                                        )
                                    });

                                    condition.then_some(Convolution::UnpolPDF(pid))
                                }
                                None => Some(Convolution::UnpolPDF(2212)),
                                Some(Err(err)) => panic!(
                                    "metadata 'initial_state_{index}' could not be parsed: {err}"
                                ),
                            }
                        }
                        (None, Some(_)) => {
                            panic!("metadata 'convolution_type_{index}' is missing")
                        }
                        (Some(_), None) => {
                            panic!("metadata 'convolution_particle_{index}' is missing")
                        }
                        (Some(Ok(_)), Some(type_)) => {
                            panic!("metadata 'convolution_type_{index} = {type_}' is unknown")
                        }
                        (Some(Err(err)), Some(_)) => {
                            panic!("metadata 'convolution_particle_{index}' could not be parsed: {err}")
                        }
                    }
                })
                .collect()
        },
    )
}
