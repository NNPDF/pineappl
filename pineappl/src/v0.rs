use super::bin::{BinLimits, BinRemapper};
use super::boc::{Channel, Order};
use super::convolutions::Convolution;
use super::empty_subgrid::EmptySubgridV1;
use super::grid::{Grid, GridError, Mmv3, MoreMembers};
use super::packed_array::PackedArray;
use super::packed_subgrid::PackedQ1X2SubgridV1;
use super::pids::PidBasis;
use super::subgrid::{Mu2, SubgridParams};
use ndarray::Array3;
use pineappl_v0::grid::Grid as GridV0;
use std::io::BufRead;
use std::iter;

pub fn read_uncompressed_v0(mut reader: impl BufRead) -> Result<Grid, GridError> {
    use pineappl_v0::subgrid::Subgrid as _;

    let grid = GridV0::read(&mut reader).map_err(|err| GridError::Other(err.into()))?;

    // TODO: convert differently if grid only has one convolution
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
                        let x1_grid = subgrid.x1_grid().into_owned();
                        let x2_grid = subgrid.x2_grid().into_owned();
                        let mut array =
                            PackedArray::new([mu2_grid.len(), x1_grid.len(), x2_grid.len()]);
                        for ((o, b, c), v) in subgrid.indexed_iter() {
                            array[[o, b, c]] = v;
                        }
                        PackedQ1X2SubgridV1::new(array, mu2_grid, x1_grid, x2_grid).into()
                    }
                })
                .collect(),
        )
        // UNWRAP: the dimensions must be the same as in the v0 grid
        .unwrap(),
        channels: grid
            .channels()
            .iter()
            .map(|c| Channel::new(c.entry().iter().map(|&(a, b, f)| (vec![a, b], f)).collect()))
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
        // TODO: remove this member
        subgrid_params: SubgridParams::default(),
        metadata: grid
            .key_values()
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .collect(),
        convolutions: read_convolutions_from_metadata(&grid),
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
        // TODO: make these proper members
        more_members: MoreMembers::V3(Mmv3 {
            remapper: grid.remapper().map(|r| {
                // UNWRAP: if the old grid could be constructed with the given normalizations
                // and limits we should be able to do the same without error
                BinRemapper::new(r.normalizations().to_vec(), r.limits().to_vec()).unwrap()
            }),
            // TODO: remove this member
            subgrid_template: EmptySubgridV1.into(),
        }),
    };

    assert_eq!(result.bin_info().bins(), grid.bin_info().bins());

    Ok(result)
}

fn read_convolutions_from_metadata(grid: &GridV0) -> Vec<Convolution> {
    grid.key_values().map_or_else(
            // if there isn't any metadata, we assume two unpolarized proton-PDFs are used
            || vec![Convolution::UnpolPDF(2212), Convolution::UnpolPDF(2212)],
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
                    (_, Some("None")) => Convolution::None,
                    (Some(Ok(pid)), Some("UnpolPDF")) => Convolution::UnpolPDF(pid),
                    (Some(Ok(pid)), Some("PolPDF")) => Convolution::PolPDF(pid),
                    (Some(Ok(pid)), Some("UnpolFF")) => Convolution::UnpolFF(pid),
                    (Some(Ok(pid)), Some("PolFF")) => Convolution::PolFF(pid),
                    (None, None) => {
                        // if these key-value pairs are missing use the old metadata
                        match kv
                            .get(&format!("initial_state_{index}"))
                            .map(|s| s.parse::<i32>())
                        {
                            Some(Ok(pid)) => {
                                let condition = !grid.channels().iter().all(|entry| {
                                    entry.entry().iter().all(|&(a, b, _)| match index {1 => a, 2 => b, _ => unreachable!()} == pid)
                                });

                                if condition {
                                    Convolution::UnpolPDF(pid)
                                } else {
                                    Convolution::None
                                }
                            }
                            None => Convolution::UnpolPDF(2212),
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
