use super::helpers::{self, Subcommand};
use anyhow::{anyhow, Result};
use clap::{Parser, ValueHint};
use std::path::PathBuf;

#[cfg(feature = "fastnlo")]
mod fastnlo {
    use super::*;
    use lhapdf::Pdf;
    use pineappl::bin::BinRemapper;
    use pineappl::grid::{Grid, Order};
    use pineappl::import_only_subgrid::ImportOnlySubgridV2;
    use pineappl::lumi::LumiEntry;
    use pineappl::sparse_array3::SparseArray3;
    use pineappl::subgrid::{Mu2, SubgridParams};
    use pineappl_fastnlo::ffi::{
        self, fastNLOCoeffAddBase, fastNLOCoeffAddFix, fastNLOCoeffAddFlex,
        fastNLOPDFLinearCombinations, ESMCalculation, ESMOrder, EScaleFunctionalForm,
    };
    use std::convert::{TryFrom, TryInto};
    use std::f64::consts::TAU;
    use std::path::Path;
    use std::ptr;

    fn pid_to_pdg_id(pid: i32) -> i32 {
        if pid >= -6 && pid <= 6 {
            if pid == 0 {
                21
            } else {
                pid
            }
        } else {
            unreachable!()
        }
    }

    fn create_lumi(
        table: &fastNLOCoeffAddBase,
        comb: &fastNLOPDFLinearCombinations,
    ) -> Vec<LumiEntry> {
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

        let n_obs_bin = table_as_add_base.GetNObsBin();
        let n_subproc = table_as_add_base.GetNSubproc();
        let total_scalevars = table.GetTotalScalevars();
        let total_scalenodes: usize = table.GetTotalScalenodes().try_into().unwrap();

        for obs in 0..n_obs_bin {
            let x1_values = ffi::GetXNodes1(table_as_add_base, obs);

            // TODO: is this the correct assumption?
            let x2_values = if table.GetNxtot2(0) == -1 {
                x1_values.clone()
            } else {
                ffi::GetXNodes2(table_as_add_base, obs)
            };

            for subproc in 0..n_subproc {
                let factor = table.GetNevt(obs, subproc);

                for j in 0..total_scalevars {
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

                    let mut non_zero_subgrid = false;

                    for mu2_slice in 0..total_scalenodes {
                        let mut ix1: usize = 0;
                        let mut ix2: usize = 0;

                        for ix in 0..table.GetNxmax(obs) {
                            assert_eq!(
                                table.GetXIndex(
                                    obs,
                                    ix1.try_into().unwrap(),
                                    ix2.try_into().unwrap()
                                ),
                                ix
                            );

                            let value = table.GetSigmaTilde(
                                obs,
                                j,
                                mu2_slice.try_into().unwrap(),
                                ix,
                                subproc,
                            );

                            if value != 0.0 {
                                non_zero_subgrid = true;
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

                    if non_zero_subgrid {
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
        _table: &fastNLOCoeffAddFlex,
        _comb: &fastNLOPDFLinearCombinations,
        _mur_ff: EScaleFunctionalForm,
        _muf_ff: EScaleFunctionalForm,
        _bins: usize,
        _alpha: u32,
        _ipub_units: i32,
    ) -> Grid {
        //std::vector<uint32_t> order_params = { static_cast <uint32_t> (table->GetNpow()), alpha, 0, 0 };
        //
        //auto* lumi = pineappl_lumi_new();
        //create_lumi(table, comb, lumi);
        //
        //std::vector<double> bin_limits(bins + 1);
        //std::iota(bin_limits.begin(), bin_limits.end(), 0.0);
        //auto* key_vals = pineappl_keyval_new();
        //
        //// flexible grids always have a hadron in initial state 1 ...
        //pineappl_keyval_set_string(key_vals, "initial_state_1",
        //    std::to_string(table->GetPDFPDG(0)).c_str());
        //// and something else at 2
        //pineappl_keyval_set_string(key_vals, "initial_state_2", "11");
        //
        //auto* pgrid = pineappl_grid_new(lumi, 1, order_params.data(), bins, bin_limits.data(),
        //    key_vals);
        //pineappl_keyval_delete(key_vals);
        //pineappl_lumi_delete(lumi);
        //
        //std::size_t n_obs_bin = table->GetNObsBin();
        //std::size_t n_subproc = table->GetNSubproc();
        //
        //auto const& sigma_tildes = table->GetSigmaTildes();
        //
        //auto const rescale = std::pow(0.1, table->GetIXsectUnits() - ipub_units);
        //
        //for (std::size_t obs = 0; obs != n_obs_bin; ++obs)
        //{
        //    auto const& scale_nodes1 = table->GetScaleNodes1(obs);
        //    auto const& scale_nodes2 = table->GetScaleNodes2(obs);
        //    auto const& x1_values = table->GetXNodes1(obs);
        //
        //    std::vector<double> mur2_values;
        //
        //    switch (mur_ff)
        //    {
        //    case fastNLO::kScale1:
        //        for (auto const s1 : scale_nodes1)
        //        {
        //            for (std::size_t i = 0; i != scale_nodes2.size(); ++i)
        //            {
        //                mur2_values.push_back(s1 * s1);
        //            }
        //        }
        //        break;
        //
        //    case fastNLO::kScale2:
        //        for (std::size_t i = 0; i != scale_nodes1.size(); ++i)
        //        {
        //            for (auto const s2 : scale_nodes2)
        //            {
        //                mur2_values.push_back(s2 * s2);
        //            }
        //        }
        //        break;
        //
        //    case fastNLO::kQuadraticSum:
        //        for (auto const s1 : scale_nodes1)
        //        {
        //            for (auto const s2 : scale_nodes2)
        //            {
        //                mur2_values.push_back(s1 * s1 + s2 * s2);
        //            }
        //        }
        //        break;
        //
        //    case fastNLO::kQuadraticMean:
        //        for (auto const s1 : scale_nodes1)
        //        {
        //            for (auto const s2 : scale_nodes2)
        //            {
        //                mur2_values.push_back(0.5 * (s1 * s1 + s2 * s2));
        //            }
        //        }
        //        break;
        //
        //    default:
        //        // TODO: NYI
        //        assert( false );
        //    }
        //
        //    std::vector<double> muf2_values;
        //
        //    switch (muf_ff)
        //    {
        //    case fastNLO::kScale1:
        //        for (auto const s1 : scale_nodes1)
        //        {
        //            for (std::size_t i = 0; i != scale_nodes2.size(); ++i)
        //            {
        //                muf2_values.push_back(s1 * s1);
        //            }
        //        }
        //        break;
        //
        //    case fastNLO::kScale2:
        //        for (std::size_t i = 0; i != scale_nodes1.size(); ++i)
        //        {
        //            for (auto const s2 : scale_nodes2)
        //            {
        //                muf2_values.push_back(s2 * s2);
        //            }
        //        }
        //        break;
        //
        //    case fastNLO::kQuadraticSum:
        //        for (auto const s1 : scale_nodes1)
        //        {
        //            for (auto const s2 : scale_nodes2)
        //            {
        //                muf2_values.push_back(s1 * s1 + s2 * s2);
        //            }
        //        }
        //        break;
        //
        //    case fastNLO::kQuadraticMean:
        //        for (auto const s1 : scale_nodes1)
        //        {
        //            for (auto const s2 : scale_nodes2)
        //            {
        //                muf2_values.push_back(0.5 * (s1 * s1 + s2 * s2));
        //            }
        //        }
        //        break;
        //
        //    default:
        //        // TODO: NYI
        //        assert( false );
        //    }
        //
        //    std::vector<double> mu2_values;
        //
        //    for (std::size_t i = 0; i != scale_nodes1.size() * scale_nodes2.size(); ++i)
        //    {
        //        mu2_values.push_back(mur2_values.at(i));
        //        mu2_values.push_back(muf2_values.at(i));
        //    }
        //
        //    std::vector<double> x2_values{1.0};
        //
        //    for (std::size_t subproc = 0; subproc != n_subproc; ++subproc)
        //    {
        //        auto* subgrid = pineappl_subgrid_new2(mu2_values.size() / 2, mu2_values.data(),
        //            x1_values.size(), x1_values.data(), x2_values.size(), x2_values.data());
        //
        //        auto const factor = rescale / table->GetNevt(obs, subproc);
        //        bool non_zero_subgrid = false;
        //
        //        std::size_t mu2_slice = 0;
        //
        //        for (std::size_t is1 = 0; is1 != scale_nodes1.size(); ++is1)
        //        {
        //            for (std::size_t is2 = 0; is2 != scale_nodes2.size(); ++is2)
        //            {
        //                std::vector<double> slice(x1_values.size());
        //                bool non_zero = false;
        //                auto const logmur2 = std::log(mu2_values.at(2 * mu2_slice + 0));
        //                auto const logmuf2 = std::log(mu2_values.at(2 * mu2_slice + 1));
        //
        //                // flexible scale grids only allow one initial-state hadron
        //                for (std::size_t ix = 0; ix != sigma_tildes.at(0)->at(obs).size(); ++ix)
        //                {
        //                    double value = sigma_tildes.at(0)->at(obs).at(ix).at(is1).at(is2)
        //                        .at(subproc);
        //
        //                    if (table->GetNScaleDep() >= 5)
        //                    {
        //                        // mur
        //                        value += logmur2 *
        //                            sigma_tildes.at(1)->at(obs).at(ix).at(is1).at(is2).at(subproc);
        //                        // muf
        //                        value += logmuf2 *
        //                            sigma_tildes.at(2)->at(obs).at(ix).at(is1).at(is2).at(subproc);
        //
        //                        if (table->GetNScaleDep() >= 6)
        //                        {
        //                            // mur mur
        //                            value += logmur2 * logmur2 *
        //                                sigma_tildes.at(3)->at(obs).at(ix).at(is1).at(is2).at(subproc);
        //                        }
        //
        //                        if (table->GetNScaleDep() >= 7)
        //                        {
        //                            // muf muf
        //                            value += logmuf2 * logmuf2 *
        //                                sigma_tildes.at(4)->at(obs).at(ix).at(is1).at(is2).at(subproc);
        //                            // mur muf
        //                            value += logmur2 * logmuf2 *
        //                                sigma_tildes.at(5)->at(obs).at(ix).at(is1).at(is2).at(subproc);
        //                        }
        //                    }
        //
        //                    if (value != 0.0)
        //                    {
        //                        non_zero = true;
        //                        slice.at(ix) = value * factor * x1_values.at(ix);
        //                    }
        //                }
        //
        //                if (non_zero)
        //                {
        //                    non_zero_subgrid = true;
        //                    pineappl_subgrid_import_mu2_slice(subgrid, mu2_slice, slice.data());
        //                }
        //
        //                ++mu2_slice;
        //            }
        //        }
        //
        //        if (non_zero_subgrid)
        //        {
        //            pineappl_grid_replace_and_delete(pgrid, subgrid, 0, obs, subproc);
        //        }
        //        else
        //        {
        //            pineappl_subgrid_delete(subgrid);
        //        }
        //    }
        //}
        //
        //return pgrid;
        todo!()
    }

    pub fn convert_fastnlo_table(input: &Path, pdfset: &str, alpha: u32) -> Result<Grid> {
        let mut file =
            ffi::make_fastnlo_lhapdf_with_name_file_set(input.to_str().unwrap(), pdfset, 0);
        let file_as_reader = ffi::downcast_lhapdf_to_reader(file.as_ref().unwrap());
        let file_as_table = ffi::downcast_lhapdf_to_table(file.as_ref().unwrap());

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
            assert_ne!(coeff_table, ptr::null_mut());
            let linear_combinations =
                ffi::downcast_reader_to_pdf_linear_combinations(file_as_reader);

            let converted = unsafe { ffi::dynamic_cast_coeff_add_fix(coeff_table) };

            if converted != ptr::null_mut() {
                grids.push(convert_coeff_add_fix(
                    unsafe { &*converted },
                    linear_combinations,
                    bins,
                    alpha,
                ));
            } else {
                let converted = unsafe { ffi::dynamic_cast_coeff_add_flex(coeff_table) };

                if converted != ptr::null_mut() {
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
                } else {
                    // TODO: NYI
                    unreachable!();
                }
            }
        }

        let mut result = grids.remove(0);
        for grid in grids {
            result.merge(grid)?;
        }

        result.scale_by_order(1.0 / TAU, 1.0, 1.0, 1.0, 1.0);
        result.optimize();

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

        result.set_remapper(BinRemapper::new(normalizations.clone(), limits).unwrap())?;

        let file_as_reader_mut = ffi::downcast_lhapdf_to_reader_mut(file.as_mut().unwrap());
        let results = ffi::GetCrossSection(file_as_reader_mut, false);
        let pdf = pdfset
            .parse()
            .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);
        let other_results = helpers::convolute(&result, &pdf, &[], &[], &[], 1);

        let mut different = false;

        for (i, (one, mut two)) in results
            .into_iter()
            .zip(other_results.into_iter())
            .enumerate()
        {
            two *= normalizations[i];

            // catches the case where both results are zero
            if one == two {
                println!(">>> Success!");
                continue;
            }

            if (two / one - 1.0).abs() > 1e-10 {
                println!(
                    ">>> fastNLO: {} PineAPPL: {} fN/P: {} P/fN: {}",
                    one,
                    two,
                    one / two,
                    two / one
                );
                different = true;
            } else {
                println!(">>> Success!");
            }
        }

        if different {
            Err(anyhow!("grids are different"))
        } else {
            Ok(result)
        }
    }
}

/// Converts fastNLO tables to PineAPPL grids.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the converted grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// LHAPDF id or name of the PDF set to check the converted grid with.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset: String,
    /// LO coupling power in alpha.
    #[clap(default_value_t = 0)]
    alpha: u32,
}

impl Subcommand for Opts {
    #[cfg(feature = "fastnlo")]
    fn run(&self) -> Result<()> {
        let grid = fastnlo::convert_fastnlo_table(&self.input, &self.pdfset, self.alpha)?;

        helpers::write_grid(&self.output, &grid)
    }

    #[cfg(not(feature = "fastnlo"))]
    fn run(&self) -> Result<()> {
        Err(anyhow!(
            "you need to install `pineappl` with feature `fastnlo`"
        ))
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-import 
Converts fastNLO tables to PineAPPL grids

USAGE:
    pineappl import <INPUT> <OUTPUT> <PDFSET>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path to the converted grid
    <PDFSET>    LHAPDF id or name of the PDF set to check the converted grid with

OPTIONS:
    -h, --help    Print help information
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["import", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }
}
