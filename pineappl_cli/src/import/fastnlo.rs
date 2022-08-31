use anyhow::{anyhow, Result};
use itertools::Itertools;
use pineappl::bin::BinRemapper;
use pineappl::grid::{Grid, Order};
use pineappl::import_only_subgrid::ImportOnlySubgridV2;
use pineappl::lumi::LumiEntry;
use pineappl::sparse_array3::SparseArray3;
use pineappl::subgrid::{Mu2, SubgridParams};
use pineappl_fastnlo::ffi::{
    self, fastNLOCoeffAddBase, fastNLOCoeffAddFix, fastNLOCoeffAddFlex, fastNLOLHAPDF,
    fastNLOPDFLinearCombinations, ESMCalculation, ESMOrder, EScaleFunctionalForm,
};
use std::convert::{TryFrom, TryInto};
use std::f64::consts::TAU;

fn pid_to_pdg_id(pid: i32) -> i32 {
    match pid {
        -6..=-1 | 1..=6 => pid,
        0 => 21,
        _ => unimplemented!(),
    }
}

fn create_lumi(table: &fastNLOCoeffAddBase, comb: &fastNLOPDFLinearCombinations) -> Vec<LumiEntry> {
    // TODO: set this to the right value if there's only one PDF
    let lepton_id = if table.GetNPDF() == 2 { 0 } else { 11 };
    let mut lumis = Vec::new();

    // if there's a (non-empty) PDF coefficient vector reconstruct the luminosity function; the
    // advantage is that we preserve the order of the lumi entries in the PineAPPL grid
    for pdf_entry in 0..ffi::GetPDFCoeffSize(table) {
        let mut entries = Vec::new();

        for entry in ffi::GetPDFCoeff(table, pdf_entry) {
            let a = pid_to_pdg_id(entry.first);
            let b = if lepton_id == 0 {
                pid_to_pdg_id(entry.second)
            } else {
                11
            };
            let f = 1.0;

            entries.push((a, b, f));
        }

        lumis.push(LumiEntry::new(entries));
    }

    // if the PDF coefficient vector was empty, we must reconstruct the lumi function
    if lumis.is_empty() {
        let nsubproc = table.GetNSubproc().try_into().unwrap();

        let mut xfx1 = [0.0; 13];
        let mut xfx2 = [0.0; 13];

        let mut entries = Vec::new();
        entries.resize(nsubproc, Vec::new());

        for a in 0..13 {
            xfx1[a] = 1.0;

            for b in 0..13 {
                xfx2[b] = 1.0;

                let lumi = ffi::CalcPDFLinearCombination(comb, table, &xfx1, &xfx2, false);

                assert!(lumi.len() == nsubproc);

                for (i, l) in lumi.iter().copied().enumerate().filter(|(_, l)| *l != 0.0) {
                    let ap = pid_to_pdg_id(i32::try_from(a).unwrap() - 6);
                    let bp = pid_to_pdg_id(i32::try_from(b).unwrap() - 6);

                    entries[i].push((ap, bp, l));
                }

                xfx2[b] = 0.0;
            }

            xfx1[a] = 0.0;
        }

        lumis = entries.into_iter().map(LumiEntry::new).collect();
    }

    lumis
}

fn convert_coeff_add_fix(
    table: &fastNLOCoeffAddFix,
    comb: &fastNLOPDFLinearCombinations,
    bins: usize,
    alpha: u32,
) -> Grid {
    let table_as_add_base = ffi::downcast_coeff_add_fix_to_base(table);

    let mut grid = Grid::new(
        create_lumi(table_as_add_base, comb),
        vec![Order {
            alphas: table_as_add_base.GetNpow().try_into().unwrap(),
            alpha,
            logxir: 0,
            logxif: 0,
        }],
        (0..=bins).map(|limit| limit as f64).collect(),
        SubgridParams::default(),
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
                    })
                    .collect();

                let mut array =
                    SparseArray3::new(mu2_values.len(), x1_values.len(), x2_values.len());

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
                            _ => unreachable!(),
                        }
                    }
                }

                if !array.is_empty() {
                    grid.set_subgrid(
                        0,
                        obs.try_into().unwrap(),
                        subproc.try_into().unwrap(),
                        ImportOnlySubgridV2::new(
                            array,
                            mu2_values,
                            x1_values.clone(),
                            x2_values.clone(),
                        )
                        .into(),
                    );
                }
            }
        }
    }

    grid
}

fn convert_coeff_add_flex(
    table: &fastNLOCoeffAddFlex,
    comb: &fastNLOPDFLinearCombinations,
    mur_ff: EScaleFunctionalForm,
    muf_ff: EScaleFunctionalForm,
    bins: usize,
    alpha: u32,
    ipub_units: i32,
) -> Grid {
    let table_as_add_base = ffi::downcast_coeff_add_flex_to_base(table);

    let mut grid = Grid::new(
        create_lumi(table_as_add_base, comb),
        vec![Order {
            alphas: table_as_add_base.GetNpow().try_into().unwrap(),
            alpha,
            logxir: 0,
            logxif: 0,
        }],
        (0..=bins).map(|limit| limit as f64).collect(),
        SubgridParams::default(),
    );

    let npdf = table_as_add_base.GetNPDF();
    let pdf_pdg1 = table.GetPDFPDG(0).to_string();
    let pdf_pdg2 = if npdf > 1 {
        table.GetPDFPDG(1).to_string()
    } else {
        // TODO: don't hardcode this
        "11".to_string()
    };

    grid.set_key_value("initial_state_1", &pdf_pdg1);
    grid.set_key_value("initial_state_2", &pdf_pdg2);

    let rescale = 0.1_f64.powi(table.GetIXsectUnits() - ipub_units);

    for obs in 0..bins {
        let scale_nodes1 = ffi::GetScaleNodes1(table, obs.try_into().unwrap());
        let scale_nodes2 = ffi::GetScaleNodes2(table, obs.try_into().unwrap());
        let x1_values = ffi::GetXNodes1(table_as_add_base, obs.try_into().unwrap());
        let x2_values = if npdf > 1 {
            ffi::GetXNodes2(table_as_add_base, obs.try_into().unwrap())
        } else {
            vec![1.0]
        };

        let mut mur2_values = Vec::new();

        match mur_ff {
            EScaleFunctionalForm::kScale1 => {
                for s1 in &scale_nodes1 {
                    for _ in 0..scale_nodes2.len() {
                        mur2_values.push(s1 * s1);
                    }
                }
            }
            EScaleFunctionalForm::kScale2 => {
                for _ in 0..scale_nodes1.len() {
                    for s2 in &scale_nodes2 {
                        mur2_values.push(s2 * s2);
                    }
                }
            }
            EScaleFunctionalForm::kQuadraticSum => {
                for s1 in &scale_nodes1 {
                    for s2 in &scale_nodes2 {
                        mur2_values.push(s1 * s1 + s2 * s2);
                    }
                }
            }
            EScaleFunctionalForm::kQuadraticMean => {
                for s1 in &scale_nodes1 {
                    for s2 in &scale_nodes2 {
                        mur2_values.push(0.5 * (s1 * s1 + s2 * s2));
                    }
                }
            }
            _ => {
                // TODO: NYI
                unimplemented!()
            }
        }

        let mut muf2_values = Vec::new();

        match muf_ff {
            EScaleFunctionalForm::kScale1 => {
                for s1 in &scale_nodes1 {
                    for _ in 0..scale_nodes2.len() {
                        muf2_values.push(s1 * s1);
                    }
                }
            }
            EScaleFunctionalForm::kScale2 => {
                for _ in 0..scale_nodes1.len() {
                    for s2 in &scale_nodes2 {
                        muf2_values.push(s2 * s2);
                    }
                }
            }
            EScaleFunctionalForm::kQuadraticSum => {
                for s1 in &scale_nodes1 {
                    for s2 in &scale_nodes2 {
                        muf2_values.push(s1 * s1 + s2 * s2);
                    }
                }
            }
            EScaleFunctionalForm::kQuadraticMean => {
                for s1 in &scale_nodes1 {
                    for s2 in &scale_nodes2 {
                        muf2_values.push(0.5 * (s1 * s1 + s2 * s2));
                    }
                }
            }
            _ => {
                // TODO: NYI
                unimplemented!()
            }
        }

        let mut mu2_values = Vec::new();

        for i in 0..(scale_nodes1.len() * scale_nodes2.len()) {
            mu2_values.push(Mu2 {
                ren: mur2_values[i],
                fac: muf2_values[i],
            });
        }

        let nx = ffi::GetNx(table, obs);

        for subproc in 0..table_as_add_base.GetNSubproc() {
            let factor = rescale / table_as_add_base.GetNevt(obs.try_into().unwrap(), subproc);
            let mut array = SparseArray3::new(mu2_values.len(), x1_values.len(), x2_values.len());

            for (mu2_slice, (is1, is2)) in (0..scale_nodes1.len())
                .cartesian_product(0..scale_nodes2.len())
                .enumerate()
            {
                let logmur2 = mu2_values[mu2_slice].ren.ln();
                let logmuf2 = mu2_values[mu2_slice].fac.ln();

                for ix in 0..nx {
                    // TODO: is this always correct? Isn't there a member function for it?
                    let ix1 = ix % x1_values.len();
                    let ix2 = ix / x1_values.len();
                    let mut value = ffi::GetSigmaTilde(table, 0, obs, ix, is1, is2, subproc);

                    if table.GetNScaleDep() >= 5 {
                        // mur
                        value += logmur2 * ffi::GetSigmaTilde(table, 1, obs, ix, is1, is2, subproc);
                        // muf
                        value += logmuf2 * ffi::GetSigmaTilde(table, 2, obs, ix, is1, is2, subproc);

                        if table.GetNScaleDep() >= 6 {
                            // mur mur
                            value += logmur2
                                * logmur2
                                * ffi::GetSigmaTilde(table, 3, obs, ix, is1, is2, subproc);
                        }

                        if table.GetNScaleDep() >= 7 {
                            // muf muf
                            value += logmuf2
                                * logmuf2
                                * ffi::GetSigmaTilde(table, 4, obs, ix, is1, is2, subproc);
                            // mur muf
                            value += logmur2
                                * logmuf2
                                * ffi::GetSigmaTilde(table, 5, obs, ix, is1, is2, subproc);
                        }
                    }

                    if value != 0.0 {
                        array[[mu2_slice, ix1, ix2]] =
                            value * factor * x1_values[ix1] * x2_values[ix2];
                    }
                }
            }

            if !array.is_empty() {
                grid.set_subgrid(
                    0,
                    obs,
                    subproc.try_into().unwrap(),
                    ImportOnlySubgridV2::new(
                        array,
                        mu2_values.clone(),
                        x1_values.clone(),
                        x2_values.clone(),
                    )
                    .into(),
                );
            }
        }
    }

    grid
}

pub fn convert_fastnlo_table(file: &fastNLOLHAPDF, alpha: u32) -> Result<Grid> {
    let file_as_reader = ffi::downcast_lhapdf_to_reader(file);
    let file_as_table = ffi::downcast_lhapdf_to_table(file);

    let ids: Vec<_> = [
        file_as_reader.ContrId(ESMCalculation::kFixedOrder, ESMOrder::kLeading),
        file_as_reader.ContrId(ESMCalculation::kFixedOrder, ESMOrder::kNextToLeading),
        file_as_reader.ContrId(ESMCalculation::kFixedOrder, ESMOrder::kNextToNextToLeading),
    ]
    .iter()
    .copied()
    .filter(|&id| id >= 0)
    .collect();

    let bins: usize = file_as_table.GetNObsBin().try_into().unwrap();
    let mut grids = Vec::new();

    for id in ids {
        let coeff_table = file_as_table.GetCoeffTable(id);

        if coeff_table.is_null() {
            return Err(anyhow!("could not access coefficient table"));
        }

        let linear_combinations = ffi::downcast_reader_to_pdf_linear_combinations(file_as_reader);

        let converted = unsafe { ffi::dynamic_cast_coeff_add_fix(coeff_table) };

        if converted.is_null() {
            let converted = unsafe { ffi::dynamic_cast_coeff_add_flex(coeff_table) };

            if converted.is_null() {
                // TODO: NYI
                unreachable!();
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
    let normalizations = if file_as_table.IsNorm() {
        ffi::GetBinSize(file_as_table)
    } else {
        vec![1.0; bins]
    };
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

    Ok(result)
}
