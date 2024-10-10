use anyhow::Result;
// use itertools::Itertools;
// use ndarray::s;
use pineappl::bin::BinRemapper;
use pineappl::boc::{Channel, Kinematics, Order, ScaleFuncForm, Scales};
use pineappl::convolutions::Convolution;
use pineappl::grid::Grid;
use pineappl::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use pineappl::packed_array::PackedArray;
use pineappl::packed_subgrid::PackedQ1X2SubgridV1;
use pineappl::pids::PidBasis;
use pineappl::subgrid::{Mu2, NodeValues};
use pineappl_fastnlo::ffi::{
    self, fastNLOCoeffAddBase, fastNLOCoeffAddFix, fastNLOCoeffAddFlex, fastNLOLHAPDF,
    fastNLOPDFLinearCombinations, EScaleFunctionalForm,
};
use std::f64::consts::TAU;

fn pid_to_pdg_id(pid: i32) -> i32 {
    match pid {
        -6..=-1 | 1..=6 => pid,
        0 => 21,
        _ => unimplemented!("pid = {pid} is not supported"),
    }
}

fn reconstruct_channels(
    table: &fastNLOCoeffAddBase,
    comb: &fastNLOPDFLinearCombinations,
) -> Vec<Channel> {
    let mut channels = Vec::new();
    let npdf = table.GetNPDF();

    // if there's a (non-empty) PDF coefficient vector reconstruct the channels; the advantage is
    // that we preserve the order of the channels in the PineAPPL grid
    for pdf_entry in 0..ffi::GetPDFCoeffSize(table) {
        let mut entries = Vec::new();

        for entry in ffi::GetPDFCoeff(table, pdf_entry) {
            let mut pids = vec![pid_to_pdg_id(entry.first)];

            if npdf == 2 {
                pids.push(pid_to_pdg_id(entry.second));
            }

            entries.push((pids, 1.0));
        }

        channels.push(Channel::new(entries));
    }

    // if the PDF coefficient vector was empty, we must reconstruct the channels in a different way
    if channels.is_empty() {
        assert_eq!(npdf, 2);

        let nsubproc = table.GetNSubproc().try_into().unwrap();

        let mut xfx1 = [0.0; 13];
        let mut xfx2 = [0.0; 13];

        let mut entries = Vec::new();
        entries.resize(nsubproc, Vec::new());

        for a in 0..13 {
            xfx1[a] = 1.0;

            for b in 0..13 {
                xfx2[b] = 1.0;

                let channel = ffi::CalcPDFLinearCombination(comb, table, &xfx1, &xfx2, false);

                assert!(channel.len() == nsubproc);

                for (i, &l) in channel.iter().enumerate().filter(|(_, &l)| l != 0.0) {
                    let ap = pid_to_pdg_id(i32::try_from(a).unwrap() - 6);
                    let bp = pid_to_pdg_id(i32::try_from(b).unwrap() - 6);

                    entries[i].push((vec![ap, bp], l));
                }

                xfx2[b] = 0.0;
            }

            xfx1[a] = 0.0;
        }

        channels = entries.into_iter().map(Channel::new).collect();
    }

    channels
}

fn convert_coeff_add_fix(
    table: &fastNLOCoeffAddFix,
    comb: &fastNLOPDFLinearCombinations,
    bins: usize,
    alpha: u32,
) -> Grid {
    let table_as_add_base = ffi::downcast_coeff_add_fix_to_base(table);

    // UNWRAP: shouldn't be larger than `2`
    let npdf = usize::try_from(table_as_add_base.GetNPDF()).unwrap();
    assert!(npdf <= 2);

    // TODO: extract the proper convolution PIDs
    let convolutions = vec![Convolution::UnpolPDF(2212); npdf];

    let mut grid = Grid::new(
        PidBasis::Pdg,
        reconstruct_channels(table_as_add_base, comb),
        vec![Order {
            alphas: table_as_add_base.GetNpow().try_into().unwrap(),
            alpha,
            logxir: 0,
            logxif: 0,
            logxia: 0,
        }],
        (0..=bins)
            .map(|limit| u16::try_from(limit).unwrap().into())
            .collect(),
        convolutions,
        // TODO: read out interpolation parameters from fastNLO
        if npdf == 2 {
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
        // TODO: change kinematics for DIS
        if npdf == 2 {
            vec![Kinematics::Scale(0), Kinematics::X1, Kinematics::X2]
        } else {
            vec![Kinematics::Scale(0), Kinematics::X1]
        },
        Scales {
            ren: ScaleFuncForm::Scale(0),
            fac: ScaleFuncForm::Scale(0),
            frg: ScaleFuncForm::NoScale,
        },
    );

    let total_scalenodes: usize = table.GetTotalScalenodes().try_into().unwrap();

    for obs in 0..table_as_add_base.GetNObsBin() {
        let x1_values = ffi::GetXNodes1(table_as_add_base, obs);

        // TODO: is this the correct assumption?
        let x2_values = if table.GetNxtot2(0) == -1 {
            x1_values.clone()
        } else {
            ffi::GetXNodes2(table_as_add_base, obs)
        };

        for subproc in 0..table_as_add_base.GetNSubproc() {
            let factor = table_as_add_base.GetNevt(obs, subproc);

            for j in 0..table.GetTotalScalevars() {
                // TODO: for the time being we only extract the central scale result
                if table.GetScaleFactor(j) != 1.0 {
                    continue;
                }

                let q_values = ffi::GetScaleNodes(table, obs, j);
                let mu2_values: Vec<_> = q_values
                    .iter()
                    .map(|q| Mu2 {
                        ren: q * q,
                        fac: q * q,
                        frg: -1.0,
                    })
                    .collect();

                let mut array =
                    PackedArray::new(vec![mu2_values.len(), x1_values.len(), x2_values.len()]);

                // TODO: figure out what the general case is supposed to be
                assert_eq!(j, 0);

                for mu2_slice in 0..total_scalenodes {
                    let mut ix1: usize = 0;
                    let mut ix2: usize = 0;

                    for ix in 0..table.GetNxmax(obs) {
                        assert_eq!(
                            table.GetXIndex(obs, ix1.try_into().unwrap(), ix2.try_into().unwrap()),
                            ix
                        );

                        let value =
                            table.GetSigmaTilde(obs, j, mu2_slice.try_into().unwrap(), ix, subproc);

                        if value != 0.0 {
                            array[[mu2_slice, ix2, ix1]] =
                                value / factor * x1_values[ix1] * x2_values[ix2];
                        }

                        ix1 += 1;

                        match table.GetNPDFDim() {
                            2 => {
                                if ix1 == x1_values.len() {
                                    ix1 = 0;
                                    ix2 += 1;
                                }
                            }
                            1 => {
                                if ix1 > ix2 {
                                    ix1 = 0;
                                    ix2 += 1;
                                }
                            }
                            n => unimplemented!("GetNPDFDim = {n} is not supported"),
                        }
                    }
                }

                if !array.is_empty() {
                    grid.subgrids_mut()
                        [[0, obs.try_into().unwrap(), subproc.try_into().unwrap()]] =
                        PackedQ1X2SubgridV1::new(
                            array,
                            vec![
                                NodeValues::UseThese(
                                    mu2_values.iter().map(|&Mu2 { ren, .. }| ren).collect(),
                                ),
                                NodeValues::UseThese(x1_values.clone()),
                                NodeValues::UseThese(x2_values.clone()),
                            ],
                        )
                        .into();
                }
            }
        }
    }

    grid
}

fn convert_coeff_add_flex(
    _table: &fastNLOCoeffAddFlex,
    _comb: &fastNLOPDFLinearCombinations,
    _mur_ff: EScaleFunctionalForm,
    _muf_ff: EScaleFunctionalForm,
    _bins: usize,
    _alpha: u32,
    _ipub_units: i32,
) -> Grid {
    todo!()

    // let table_as_add_base = ffi::downcast_coeff_add_flex_to_base(table);

    // let alphas = table_as_add_base.GetNpow().try_into().unwrap();
    // let orders: Vec<_> = [
    //     Order::new(alphas, alpha, 0, 0, 0),
    //     Order::new(alphas, alpha, 1, 0, 0),
    //     Order::new(alphas, alpha, 0, 1, 0),
    //     Order::new(alphas, alpha, 2, 0, 0),
    //     Order::new(alphas, alpha, 0, 2, 0),
    //     Order::new(alphas, alpha, 1, 1, 0),
    // ]
    // .into_iter()
    // .take(match table.GetNScaleDep() {
    //     0..=4 => 1,
    //     5 => 3,
    //     6 => 4,
    //     7 => 6,
    //     _ => unimplemented!(),
    // })
    // .collect();
    // let orders_len = orders.len();

    // let npdf = table_as_add_base.GetNPDF();
    // assert!(npdf <= 2);

    // let convolutions = (0..2)
    //     .map(|index| {
    //         if index < npdf {
    //             Convolution::UnpolPDF(table.GetPDFPDG(index))
    //         } else {
    //             Convolution::None
    //         }
    //     })
    //     .collect();

    // let mut grid = Grid::new(
    //     PidBasis::Pdg,
    //     reconstruct_channels(table_as_add_base, comb, dis_pid),
    //     orders,
    //     (0..=bins)
    //         .map(|limit| u16::try_from(limit).unwrap().into())
    //         .collect(),
    //     convolutions,
    //     // TODO: read out interpolation parameters from fastNLO
    //     vec![
    //         Interp::new(
    //             1e2,
    //             1e8,
    //             40,
    //             3,
    //             ReweightMeth::NoReweight,
    //             Map::ApplGridH0,
    //             InterpMeth::Lagrange,
    //         ),
    //         Interp::new(
    //             2e-7,
    //             1.0,
    //             50,
    //             3,
    //             ReweightMeth::ApplGridX,
    //             Map::ApplGridF2,
    //             InterpMeth::Lagrange,
    //         ),
    //         Interp::new(
    //             2e-7,
    //             1.0,
    //             50,
    //             3,
    //             ReweightMeth::ApplGridX,
    //             Map::ApplGridF2,
    //             InterpMeth::Lagrange,
    //         ),
    //     ],
    //     // TODO: change kinematics for DIS
    //     vec![
    //         Kinematics::Scale(0),
    //         Kinematics::Scale(1),
    //         Kinematics::X1,
    //         Kinematics::X2,
    //     ],
    //     Scales {
    //         ren: todo!(),
    //         fac: todo!(),
    //         frg: ScaleFuncForm::NoScale,
    //     },
    // );

    // let rescale = 0.1_f64.powi(table.GetIXsectUnits() - ipub_units);

    // for obs in 0..bins {
    //     let scale_nodes1 = ffi::GetScaleNodes1(table, obs.try_into().unwrap());
    //     let scale_nodes2 = ffi::GetScaleNodes2(table, obs.try_into().unwrap());
    //     let x1_values = ffi::GetXNodes1(table_as_add_base, obs.try_into().unwrap());
    //     let x2_values = if npdf > 1 {
    //         ffi::GetXNodes2(table_as_add_base, obs.try_into().unwrap())
    //     } else {
    //         vec![1.0]
    //     };

    //     let mu2_values: Vec<_> = scale_nodes1
    //         .iter()
    //         .cartesian_product(scale_nodes2.iter())
    //         .map(|(&s1, &s2)| Mu2 {
    //             ren: mur_ff.compute_scale(s1, s2),
    //             fac: muf_ff.compute_scale(s1, s2),
    //             frg: -1.0,
    //         })
    //         .collect();
    //     let nx = ffi::GetNx(table, obs);

    //     for subproc in 0..table_as_add_base.GetNSubproc() {
    //         let factor = rescale / table_as_add_base.GetNevt(obs.try_into().unwrap(), subproc);
    //         let mut arrays =
    //             vec![
    //                 PackedArray::new(vec![mu2_values.len(), x1_values.len(), x2_values.len()]);
    //                 orders_len
    //             ];

    //         for (mu2_slice, (is1, is2)) in (0..scale_nodes1.len())
    //             .cartesian_product(0..scale_nodes2.len())
    //             .enumerate()
    //         {
    //             let logmur2 = mu2_values[mu2_slice].ren.ln();
    //             let logmuf2 = mu2_values[mu2_slice].fac.ln();
    //             let logs00 = [
    //                 logmur2,
    //                 logmuf2,
    //                 logmur2 * logmur2,
    //                 logmuf2 * logmuf2,
    //                 logmur2 * logmuf2,
    //             ];
    //             let logs10 = [2.0 * logmur2, 0.0, logmuf2];
    //             let logs01 = [0.0, 2.0 * logmuf2, logmur2];

    //             for ix in 0..nx {
    //                 // TODO: is this always correct? Isn't there a member function for it?
    //                 let ix1 = ix % x1_values.len();
    //                 let ix2 = ix / x1_values.len();
    //                 let mut values = [0.0; 6];

    //                 for (index, value) in values.iter_mut().enumerate().take(orders_len) {
    //                     *value = ffi::GetSigmaTilde(table, index, obs, ix, is1, is2, subproc);
    //                 }

    //                 values[0] += values[1..]
    //                     .iter()
    //                     .zip(logs00.iter())
    //                     .map(|(value, log)| value * log)
    //                     .sum::<f64>();
    //                 values[1] += values[3..]
    //                     .iter()
    //                     .zip(logs10.iter())
    //                     .map(|(value, log)| value * log)
    //                     .sum::<f64>();
    //                 values[2] += values[3..]
    //                     .iter()
    //                     .zip(logs01.iter())
    //                     .map(|(value, log)| value * log)
    //                     .sum::<f64>();

    //                 for (value, array) in values
    //                     .iter()
    //                     .copied()
    //                     .zip(arrays.iter_mut())
    //                     .filter(|(value, _)| *value != 0.0)
    //                 {
    //                     array[[mu2_slice, ix1, ix2]] =
    //                         value * factor * x1_values[ix1] * x2_values[ix2];
    //                 }
    //             }
    //         }

    //         for (subgrid, array) in grid
    //             .subgrids_mut()
    //             .slice_mut(s![.., obs, usize::try_from(subproc).unwrap()])
    //             .iter_mut()
    //             .zip(arrays.into_iter())
    //         {
    //             if array.is_empty() {
    //                 continue;
    //             }

    //             *subgrid = PackedQ1X2SubgridV1::new(array, todo!()).into();
    //         }
    //     }
    // }

    // grid
}

pub fn convert_fastnlo_table(file: &fastNLOLHAPDF, alpha: u32) -> Result<Grid> {
    let file_as_reader = ffi::downcast_lhapdf_to_reader(file);
    let file_as_table = ffi::downcast_lhapdf_to_table(file);

    let bins: usize = file_as_table.GetNObsBin().try_into().unwrap();
    let mut grids = Vec::new();

    for id in 0.. {
        // TODO: there doesn't seem to be a better way than trying an index and stopping whenever a
        // NULL pointer is returned
        let coeff_base = file_as_table.GetCoeffTable(id);

        let Some(coeff_base) = (unsafe { coeff_base.as_ref() }) else {
            break;
        };

        if !coeff_base.IsEnabled() {
            break;
        }

        let linear_combinations = ffi::downcast_reader_to_pdf_linear_combinations(file_as_reader);

        let converted = unsafe { ffi::dynamic_cast_coeff_add_fix(coeff_base) };

        if converted.is_null() {
            let converted = unsafe { ffi::dynamic_cast_coeff_add_flex(coeff_base) };

            if converted.is_null() {
                unimplemented!("subclass of fastNLOCoeffBase is not supported");
            } else {
                let mur_ff = file_as_reader.GetMuRFunctionalForm();
                let muf_ff = file_as_reader.GetMuFFunctionalForm();

                grids.push(convert_coeff_add_flex(
                    unsafe { &*converted },
                    linear_combinations,
                    mur_ff,
                    muf_ff,
                    bins,
                    alpha,
                    file_as_table.GetIpublunits(),
                ));
            }
        } else {
            grids.push(convert_coeff_add_fix(
                unsafe { &*converted },
                linear_combinations,
                bins,
                alpha,
            ));
        }
    }

    let mut result = grids.remove(0);
    for grid in grids {
        result.merge(grid)?;
    }

    result.scale_by_order(1.0 / TAU, 1.0, 1.0, 1.0, 1.0);

    let dimensions: usize = file_as_table.GetNumDiffBin().try_into().unwrap();
    let mut limits = Vec::new();
    let normalizations = vec![1.0; bins];
    limits.reserve(2 * dimensions * bins);

    for i in 0..bins {
        for j in 0..dimensions {
            let pair = ffi::GetObsBinDimBounds(
                file_as_table,
                i.try_into().unwrap(),
                j.try_into().unwrap(),
            );

            limits.push((pair.first, pair.second));
        }
    }

    result.set_remapper(BinRemapper::new(normalizations, limits).unwrap())?;

    let labels = ffi::GetDimLabels(file_as_table);

    assert_eq!(labels.len(), dimensions);

    for (dimension, label) in labels.into_iter().enumerate() {
        result
            .metadata_mut()
            .insert(format!("x{}_label", dimension + 1), label);
    }

    result
        .metadata_mut()
        .insert("y_label".to_owned(), ffi::GetXSDescr(file_as_table));
    result.metadata_mut().insert(
        "fastnlo_scenario".to_owned(),
        ffi::GetScDescr(file_as_table).join("\n"),
    );

    Ok(result)
}
