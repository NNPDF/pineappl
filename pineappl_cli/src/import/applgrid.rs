use anyhow::Result;
use float_cmp::assert_approx_eq;
use lhapdf::Pdf;
use pineappl::boc::{BinsWithFillLimits, Channel, Kinematics, Order, ScaleFuncForm, Scales};
use pineappl::convolutions::{Conv, ConvType};
use pineappl::grid::Grid;
use pineappl::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use pineappl::packed_array::PackedArray;
use pineappl::pids::PidBasis;
use pineappl::subgrid::ImportSubgridV1;
use pineappl_applgrid::ffi::{self, grid};
use std::f64::consts::TAU;
use std::pin::Pin;
use std::{iter, ptr};

fn convert_to_pdg_id(pid: usize) -> i32 {
    let pid = i32::try_from(pid).unwrap() - 6;

    match pid {
        -6..=-1 | 1..=6 => pid,
        0 => 21,
        7 => 22,
        _ => unimplemented!("pid = {pid} is not supported"),
    }
}

fn reconstruct_channels(grid: &grid, order: i32) -> Vec<Channel> {
    let pdf = unsafe { &*grid.genpdf(order, false) };
    let nproc: usize = pdf.Nproc().try_into().unwrap();

    let mut channels = vec![Vec::new(); nproc];
    let mut xfx1 = [0.0; 14];
    let mut xfx2 = [0.0; 14];
    let mut results = vec![0.0; nproc];

    for a in 0..=13 {
        xfx1[a] = 1.0;

        if grid.isDIS() {
            unsafe {
                (*pdf).evaluate(xfx1.as_ptr(), ptr::null(), results.as_mut_ptr());
            }

            for i in 0..nproc {
                if results[i] != 0.0 {
                    channels[i].push((vec![convert_to_pdg_id(a)], results[i]));
                }
            }
        } else {
            for b in 0..=13 {
                xfx2[b] = 1.0;

                unsafe {
                    (*pdf).evaluate(xfx1.as_ptr(), xfx2.as_ptr(), results.as_mut_ptr());
                }

                for i in 0..nproc {
                    if results[i] != 0.0 {
                        channels[i]
                            .push((vec![convert_to_pdg_id(a), convert_to_pdg_id(b)], results[i]));
                    }
                }

                xfx2[b] = 0.0;
            }
        }

        xfx1[a] = 0.0;
    }

    channels.into_iter().map(Channel::new).collect()
}

pub fn convert_applgrid(grid: Pin<&mut grid>, alpha: u8) -> Result<Grid> {
    let bin_limits: Vec<_> = (0..=grid.Nobs_internal())
        .map(|i| grid.obslow_internal(i))
        .collect();

    let leading_order: u8 = grid
        .leadingOrder()
        .try_into()
        // UNWRAP: exponents of orders shouldn't be larger than 255
        .unwrap();
    let orders;
    let alphas_factor;

    if grid.calculation() == ffi::grid_CALCULATION::AMCATNLO {
        alphas_factor = 2.0 * TAU;
        orders = if grid.nloops() == 0 {
            vec![Order::new(leading_order, alpha, 0, 0, 0)]
        } else if grid.nloops() == 1 {
            vec![
                Order::new(leading_order + 1, alpha, 0, 0, 0), // NLO
                Order::new(leading_order + 1, alpha, 1, 0, 0), // NLO mur
                Order::new(leading_order + 1, alpha, 0, 1, 0), // NLO muf
                Order::new(leading_order, alpha, 0, 0, 0),     // LO
            ]
        } else {
            unimplemented!("nloops = {} is not supported", grid.nloops());
        };
    } else if grid.calculation() == ffi::grid_CALCULATION::STANDARD {
        alphas_factor = 1.0 / TAU;
        orders = (0..=grid.nloops())
            .map(|power| {
                Order::new(
                    leading_order
                        + u8::try_from(power)
                            // UNWRAP: exponents of orders shouldn't be larger than 255
                            .unwrap(),
                    alpha,
                    0,
                    0,
                    0,
                )
            })
            .collect();
    } else {
        unimplemented!("calculation is not supported");
    }

    // this setting isn't supported
    assert!(!grid.getApplyCorrections());

    for bin in 0..grid.Nobs_internal() {
        // this setting isn't supported either
        assert!(!grid.getApplyCorrection(bin.try_into().unwrap()));
    }

    // this setting isn't supported
    assert_approx_eq!(f64, grid.getDynamicScale(), 0.0, ulps = 4);

    let mut grids = Vec::with_capacity(orders.len());
    let dis = grid.isDIS();

    // from APPLgrid alone we don't know what type of convolution we have
    let convolutions = vec![Conv::new(ConvType::UnpolPDF, 2212); if dis { 1 } else { 2 }];
    // TODO: read out interpolation parameters from APPLgrid
    let mut interps = vec![Interp::new(
        1e2,
        1e8,
        40,
        3,
        ReweightMeth::NoReweight,
        Map::ApplGridH0,
        InterpMeth::Lagrange,
    )];
    for _ in 0..convolutions.len() {
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

    for (i, order) in orders.into_iter().enumerate() {
        let channels = reconstruct_channels(&grid, i.try_into().unwrap());
        let lumis_len = channels.len();
        let mut pgrid = Grid::new(
            BinsWithFillLimits::from_fill_limits(bin_limits.clone())
                // UNWRAP: panic here denotes a possible incompatibility between APPLgrid and
                // PineAPPL
                .unwrap(),
            vec![order],
            channels,
            PidBasis::Pdg,
            convolutions.clone(),
            interps.clone(),
            iter::once(Kinematics::Scale(0))
                .chain((0..convolutions.len()).map(Kinematics::X))
                .collect(),
            Scales {
                ren: ScaleFuncForm::Scale(0),
                fac: ScaleFuncForm::Scale(0),
                frg: ScaleFuncForm::NoScale,
            },
        );

        for bin in 0..grid.Nobs_internal() {
            let igrid = grid.weightgrid(i.try_into().unwrap(), bin);
            let igrid = unsafe { &*igrid };
            let reweight = ffi::igrid_m_reweight(igrid);

            let scale_values: Vec<_> = (0..igrid.Ntau()).map(|i| igrid.getQ2(i)).collect();
            let x1_values: Vec<_> = (0..igrid.Ny1())
                .map(|i| igrid.getx1(i).clamp(0.0, 1.0))
                .collect();
            let x1_weights: Vec<_> = x1_values
                .iter()
                .map(|&x1| if reweight { ffi::weightfun(x1) } else { 1.0 })
                .collect();
            let x2_values: Vec<_> = (0..igrid.Ny2())
                .map(|i| igrid.getx2(i).clamp(0.0, 1.0))
                .collect();
            let x2_weights: Vec<_> = x2_values
                .iter()
                .map(|&x2| {
                    if reweight && !grid.isDIS() {
                        ffi::weightfun(x2)
                    } else {
                        1.0
                    }
                })
                .collect();

            for lumi in 0..lumis_len {
                let matrix = igrid.weightgrid(lumi.try_into().unwrap());

                if matrix.is_null() {
                    continue;
                }

                let matrix = unsafe { &*matrix };

                let mut array = PackedArray::new(if dis {
                    vec![scale_values.len(), x1_values.len()]
                } else {
                    vec![scale_values.len(), x1_values.len(), x2_values.len()]
                });

                for itau in 0..scale_values.len() {
                    for ix1 in 0..x1_values.len() {
                        for ix2 in 0..x2_values.len() {
                            let value = ffi::sparse_matrix_get(
                                matrix,
                                itau.try_into().unwrap(),
                                ix1.try_into().unwrap(),
                                ix2.try_into().unwrap(),
                            );

                            if value != 0.0 {
                                if dis {
                                    array[[itau, ix1]] = value * x1_weights[ix1];
                                } else {
                                    array[[itau, ix1, ix2]] =
                                        value * x1_weights[ix1] * x2_weights[ix2];
                                }
                            }
                        }
                    }
                }

                if !array.is_empty() {
                    pgrid.subgrids_mut()[[0, bin.try_into().unwrap(), lumi]] =
                        ImportSubgridV1::new(
                            array,
                            if dis {
                                vec![scale_values.clone(), x1_values.clone()]
                            } else {
                                vec![scale_values.clone(), x1_values.clone(), x2_values.clone()]
                            },
                        )
                        .into();
                }
            }
        }

        grids.push(pgrid);
    }

    let mut grid0 = grids.remove(0);
    grids
        .into_iter()
        .try_for_each(|g| grid0.merge(g).map_err(anyhow::Error::new))?;

    let combine = ffi::grid_combine(&grid);

    if !combine.is_empty() {
        assert_eq!(combine.iter().sum::<i32>(), grid.Nobs_internal());
    }

    let mut ranges = Vec::new();
    let mut index = 0;

    for bins in combine {
        let bins = usize::try_from(bins).unwrap();
        ranges.push(index..(index + bins));
        index += bins;
    }

    for range in ranges.into_iter().rev() {
        grid0.merge_bins(range)?;
    }

    let mut global = 1.0;

    if !grid.getNormalised() {
        let factor = *grid.run();

        if factor != 0.0 {
            global = 1.0 / factor;
        }
    }

    grid0.scale_by_order(alphas_factor, 1.0, 1.0, 1.0, 1.0, global);

    Ok(grid0)
}

pub fn convolve_applgrid(grid: Pin<&mut grid>, conv_funs: &mut [Pdf]) -> Vec<f64> {
    let nloops = grid.nloops();

    // TODO: add support for convolving an APPLgrid with two functions
    assert_eq!(conv_funs.len(), 1);

    pineappl_applgrid::grid_convolve_with_one(grid, &mut conv_funs[0], nloops, 1.0, 1.0, 1.0)
}
