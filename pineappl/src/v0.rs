use super::boc::{Bin, BinsWithFillLimits, Channel, Kinematics, Order, ScaleFuncForm, Scales};
use super::convert;
use super::convolutions::{Conv, ConvType};
use super::error::{Error, Result};
use super::grid::Grid;
use super::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use super::packed_array::PackedArray;
use super::pids::PidBasis;
use super::subgrid::{self, ImportSubgridV1};
use pineappl_v0::grid::Grid as GridV0;
use std::io::BufRead;

pub fn default_interps(flexible_scale: bool, convolutions: usize) -> Vec<Interp> {
    let scales = if flexible_scale { 2 } else { 1 };
    let mut interps = Vec::with_capacity(convolutions + scales);
    for _ in 0..scales {
        interps.push(Interp::new(
            1e2,
            1e8,
            40,
            3,
            ReweightMeth::NoReweight,
            Map::ApplGridH0,
            InterpMeth::Lagrange,
        ));
    }

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

pub fn read_uncompressed_v0(mut reader: impl BufRead) -> Result<Grid> {
    use pineappl_v0::subgrid::Subgrid as _;

    let grid = GridV0::read(&mut reader).map_err(|err| Error::Other(err.into()))?;
    let convolutions = read_convolutions_from_metadata(&grid);

    let flexible_scale_grid = grid.subgrids().iter().any(|subgrid| {
        subgrid
            .mu2_grid()
            .iter()
            .any(|mu2v0| !subgrid::node_value_eq(mu2v0.ren, mu2v0.fac))
    });
    let mut kinematics = vec![Kinematics::Scale(0)];
    if flexible_scale_grid {
        kinematics.push(Kinematics::Scale(1));
    }
    kinematics.push(Kinematics::X(0));
    if convolutions[0].is_some() && convolutions[1].is_some() {
        kinematics.push(Kinematics::X(1));
    }

    assert_eq!(convolutions.len(), 2);

    let fill_limits = if grid.remapper().is_none() {
        // if there's no BinRemapper, the limits must one dimensional and we simply take the first
        // dimension
        grid.bin_info()
            .limits()
            .into_iter()
            .map(|limits| limits[0].0)
            .chain(
                grid.bin_info()
                    .limits()
                    .into_iter()
                    .map(|limits| limits[0].1)
                    .next_back(),
            )
            .collect()
    } else {
        // if there's a BinRemapper, we use the canonical fill limits 0, 1, 2, ...
        (0..=grid.bin_info().bins())
            .map(convert::f64_from_usize)
            .collect()
    };
    let bins = BinsWithFillLimits::new(
        grid.bin_info()
            .limits()
            .into_iter()
            .zip(grid.bin_info().normalizations())
            .map(|(limits, normalization)| Bin::new(limits, normalization))
            .collect(),
        fill_limits,
    )
    // UNWRAP: if we could build a v0 grid, we should be able to build v1 grid with panicking
    .unwrap();
    let mut result = Grid::new(
        bins,
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
                            }
                            (pids, f)
                        })
                        .collect(),
                )
            })
            .collect(),
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
        convolutions.clone().into_iter().flatten().collect(),
        default_interps(
            flexible_scale_grid,
            convolutions.clone().into_iter().flatten().count(),
        ),
        kinematics,
        Scales {
            ren: ScaleFuncForm::Scale(0),
            fac: if flexible_scale_grid {
                ScaleFuncForm::Scale(1)
            } else {
                ScaleFuncForm::Scale(0)
            },
            frg: ScaleFuncForm::NoScale,
        },
    );

    for (new_subgrid, old_subgrid) in result.subgrids_mut().iter_mut().zip(grid.subgrids().iter()) {
        if !old_subgrid.is_empty() {
            let scale_node_values: Vec<Vec<_>> = if flexible_scale_grid {
                let mut ren: Vec<_> = old_subgrid
                    .mu2_grid()
                    .iter()
                    .map(|mu2v0| mu2v0.ren)
                    .collect();
                ren.sort_unstable_by(f64::total_cmp);
                ren.dedup_by(subgrid::node_value_eq_ref_mut);
                let mut fac: Vec<_> = old_subgrid
                    .mu2_grid()
                    .iter()
                    .map(|mu2v0| mu2v0.fac)
                    .collect();
                fac.sort_unstable_by(f64::total_cmp);
                fac.dedup_by(subgrid::node_value_eq_ref_mut);
                vec![ren, fac]
            } else {
                vec![old_subgrid
                    .mu2_grid()
                    .iter()
                    .map(|mu2v0| {
                        // TODO: implement importing flexible-scale grids
                        assert!(subgrid::node_value_eq(mu2v0.ren, mu2v0.fac));

                        mu2v0.fac
                    })
                    .collect()]
            };

            let mut dim = if flexible_scale_grid {
                vec![scale_node_values[0].len(), scale_node_values[1].len()]
            } else {
                vec![scale_node_values[0].len()]
            };
            if convolutions[0].is_some() {
                dim.push(old_subgrid.x1_grid().len());
            }
            if convolutions[1].is_some() {
                dim.push(old_subgrid.x2_grid().len());
            }
            let mut array = PackedArray::new(dim);

            if flexible_scale_grid {
                for (index, v) in old_subgrid.indexed_iter() {
                    let index_r = scale_node_values[0]
                        .iter()
                        .position(|&ren| {
                            subgrid::node_value_eq(ren, old_subgrid.mu2_grid()[index.0].ren)
                        })
                        // UNWRAP: if no index can be found, `scale_node_values` isn't sorted
                        // correctly
                        .unwrap();
                    let index_f = scale_node_values[1]
                        .iter()
                        .position(|&fac| {
                            subgrid::node_value_eq(fac, old_subgrid.mu2_grid()[index.0].fac)
                        })
                        // UNWRAP: if no index can be found, `scale_node_values` isn't sorted
                        // correctly
                        .unwrap();
                    if convolutions[0].is_none() {
                        array[[index_r, index_f, index.2]] = v;
                    } else if convolutions[1].is_none() {
                        array[[index_r, index_f, index.1]] = v;
                    } else {
                        array[[index_r, index_f, index.1, index.2]] = v;
                    }
                }
            } else if convolutions[0].is_none() {
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

            let mut node_values = scale_node_values;
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

    assert_eq!(result.bwfl().len(), grid.bin_info().bins());

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
                                    // Addresses: https://github.com/NNPDF/pineappl/issues/334
                                    if grid.channels().is_empty() && pid == 2212 {
                                        Some(Conv::new(ConvType::UnpolPDF, 2212))
                                    } else {
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
