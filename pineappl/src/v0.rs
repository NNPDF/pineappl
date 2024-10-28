use super::bin::BinRemapper;
use super::boc::{Channel, Kinematics, Order, ScaleFuncForm, Scales};
use super::convolutions::{Conv, ConvType};
use super::grid::{Grid, GridError};
use super::import_subgrid::ImportSubgridV1;
use super::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use super::packed_array::PackedArray;
use super::pids::PidBasis;
use super::subgrid::Mu2;
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

    let mut result = Grid::new(
        grid.key_values()
            .and_then(|kv| kv.get("lumi_id_types"))
            // TODO: use PidBasis::from_str
            .map_or(PidBasis::Pdg, |lumi_id_types| {
                match lumi_id_types.as_str() {
                    "pdg_mc_ids" => PidBasis::Pdg,
                    "evol" => PidBasis::Evol,
                    _ => panic!("unknown PID basis '{lumi_id_types}'"),
                }
            }),
        grid.channels()
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
        grid.orders()
            .iter()
            .map(|o| Order {
                // UNWRAP: there shouldn't be orders with exponents larger than 255
                alphas: o.alphas.try_into().unwrap(),
                alpha: o.alpha.try_into().unwrap(),
                logxir: o.logxir.try_into().unwrap(),
                logxif: o.logxif.try_into().unwrap(),
                logxia: 0,
            })
            .collect(),
        if grid.remapper().is_none() {
            let limits = &grid.bin_info().limits();
            iter::once(limits[0][0].0)
                .chain(limits.iter().map(|v| v[0].1))
                .collect()
        } else {
            // if there is a BinRemapper this member will likely have no impact
            (0..=grid.bin_info().bins())
                .map(|i| f64::from(u16::try_from(i).unwrap()))
                .collect()
        },
        convolutions.clone().into_iter().flatten().collect(),
        default_interps(convolutions.clone().into_iter().flatten().count()),
        kinematics,
        Scales {
            // TODO: read in flexible-scale grids properly
            ren: ScaleFuncForm::Scale(0),
            fac: ScaleFuncForm::Scale(0),
            frg: ScaleFuncForm::NoScale,
        },
    );

    for (new_subgrid, old_subgrid) in result.subgrids_mut().iter_mut().zip(grid.subgrids().iter()) {
        if !old_subgrid.is_empty() {
            let mu2_grid: Vec<_> = old_subgrid
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
                dim.push(old_subgrid.x1_grid().len());
            }
            if convolutions[1].is_some() {
                dim.push(old_subgrid.x2_grid().len());
            }
            let mut array = PackedArray::new(dim);

            if convolutions[0].is_none() {
                for (index, v) in old_subgrid.indexed_iter() {
                    array[[index.0, index.2]] = v;
                }
            } else if convolutions[1].is_none() {
                for (index, v) in old_subgrid.indexed_iter() {
                    array[[index.0, index.1]] = v;
                }
            } else {
                for (index, v) in old_subgrid.indexed_iter() {
                    array[<[usize; 3]>::from(index)] = v;
                }
            }

            let mut node_values = vec![mu2_grid.iter().map(|&Mu2 { ren, .. }| ren).collect()];
            if convolutions[0].is_some() {
                node_values.push(old_subgrid.x1_grid().into_owned());
            }
            if convolutions[1].is_some() {
                node_values.push(old_subgrid.x2_grid().into_owned());
            }
            *new_subgrid = ImportSubgridV1::new(array, node_values).into();
        }
    }

    *result.metadata_mut() = grid
        .key_values()
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .collect();

    if let Some(r) = grid.remapper() {
        result
            .set_remapper(
                BinRemapper::new(r.normalizations().to_vec(), r.limits().to_vec())
                    // UNWRAP: if the old grid could be constructed with the given normalizations
                    // and limits we should be able to do the same without error
                    .unwrap(),
            )
            // UNWRAP: there's a bug if this fails
            .unwrap();
    }

    assert_eq!(result.bin_info().bins(), grid.bin_info().bins());

    Ok(result)
}

fn read_convolutions_from_metadata(grid: &GridV0) -> Vec<Option<Conv>> {
    grid.key_values().map_or_else(
        // if there isn't any metadata, we assume two unpolarized proton-PDFs are used
        || vec![Some(Conv::new(ConvType::UnpolPDF, 2212)); 2],
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
                        (Some(Ok(pid)), Some("UnpolPDF")) => Some(Conv::new(ConvType::UnpolPDF, pid)),
                        (Some(Ok(pid)), Some("PolPDF")) => Some(Conv::new(ConvType::PolPDF, pid)),
                        (Some(Ok(pid)), Some("UnpolFF")) => Some(Conv::new(ConvType::UnpolFF, pid)),
                        (Some(Ok(pid)), Some("PolFF")) => Some(Conv::new(ConvType::PolFF, pid)),
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

                                    condition.then_some(Conv::new(ConvType::UnpolPDF, pid))
                                }
                                None => Some(Conv::new(ConvType::UnpolPDF, 2212)),
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
