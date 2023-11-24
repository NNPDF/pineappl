use anyhow::Result;
use float_cmp::assert_approx_eq;
use lhapdf::Pdf;
use num_complex::Complex;
use pineappl::bin::BinRemapper;
use pineappl::grid::{Grid, GridOptFlags, Ntuple, Order};
use pineappl::lumi::LumiCache;
use pineappl::lumi_entry;
use pineappl::subgrid::{ExtraSubgridParams, Subgrid, SubgridEnum, SubgridParams};
use rand::Rng;
use rand_pcg::Pcg64;
use std::f64::consts::PI;
use std::io::Cursor;
use std::mem;

// If equation numbers are given, they are from Stefan Dittmaier and Max Huber's paper:
//   'Radiative corrections to the neutral-current Drell–Yan process in the Standard Model and its
//    minimal supersymmetric extension' (https://arxiv.org/abs/0911.2329)

// Eq. (2.13) - gamma-gamma contribution to DY lepton pair production
fn int_photo(s: f64, t: f64, u: f64) -> f64 {
    let alpha0: f64 = 1.0 / 137.03599911;
    alpha0.powi(2) / 2.0 / s * (t / u + u / t)
}

// Eq. (2.12) - quark-antiquark contribution to DY lepton pair production
fn int_quark(s: f64, t: f64, u: f64, qq: f64, i3_wq: f64) -> f64 {
    let alphagf: f64 = 1.0 / 132.30818655547878;
    let mw = 80.35198454966643;
    let mz = 91.15348061918276;
    let gw = 2.083799397775285;
    let gz = 2.494266378772824;

    // lepton charge
    let ql: f64 = -1.0;
    // lepton weak isospin
    let i3_wl = -0.5;

    // weak mixing angles
    let cw = (Complex::new(mw * mw, -mw * gw) / Complex::new(mz * mz, -mz * gz)).sqrt();
    let sw = (Complex::new(1.0, 0.0) - cw * cw).sqrt();

    // Eq. (2.8)
    let chi_z = Complex::new(s, 0.0) / Complex::new(s - mz * mz, mz * gz);

    // Eq. (2.7)
    let gp_qqz = -sw / cw * qq;
    let gm_qqz = (i3_wq - sw * sw * qq) / (sw * cw);
    let gp_llz = -sw / cw * ql;
    let gm_llz = (i3_wl - sw * sw * ql) / (sw * cw);

    alphagf.powi(2) / 12.0 / s.powi(3)
        * (2.0 * qq.powi(2) * ql.powi(2) * (t * t + u * u)
            + 2.0
                * qq
                * ql
                * (((gp_qqz * gp_llz + gm_qqz * gm_llz) * u * u
                    + (gp_qqz * gm_llz + gm_qqz * gp_llz) * t * t)
                    * chi_z)
                    .re
            + ((gp_qqz.norm_sqr() * gp_llz.norm_sqr() + gm_qqz.norm_sqr() * gm_llz.norm_sqr())
                * u
                * u
                + (gp_qqz.norm_sqr() * gm_llz.norm_sqr() + gm_qqz.norm_sqr() * gp_llz.norm_sqr())
                    * t
                    * t)
                * chi_z.norm_sqr())
}

struct Psp2to2 {
    s: f64,
    t: f64,
    u: f64,
    x1: f64,
    x2: f64,
    jacobian: f64,
}

fn hadronic_pspgen(rng: &mut impl Rng, mmin: f64, mmax: f64) -> Psp2to2 {
    let smin = mmin * mmin;
    let smax = mmax * mmax;

    let mut jacobian = 1.0;

    let r1 = rng.gen::<f64>();
    let r2 = rng.gen::<f64>();
    let tau0 = smin / smax;
    let tau = tau0.powf(r1);
    let y = tau.powf(1.0 - r2);
    let x1 = y;
    let x2 = tau / y;
    let s = tau * smax;
    jacobian *= tau * tau0.ln().powi(2) * r1;

    // theta integration (in the CMS)
    let cos_theta = 2.0 * rng.gen::<f64>() - 1.0;
    jacobian *= 2.0;

    let t = -0.5 * s * (1.0 - cos_theta);
    let u = -0.5 * s * (1.0 + cos_theta);

    // phi integration
    jacobian *= 2.0 * PI;

    Psp2to2 {
        s,
        t,
        u,
        x1,
        x2,
        jacobian,
    }
}

fn fill_drell_yan_lo_grid(
    rng: &mut impl Rng,
    calls: usize,
    subgrid_type: &str,
    dynamic: bool,
    reweight: bool,
) -> Result<Grid> {
    let lumi = vec![
        // photons
        lumi_entry![22, 22, 1.0],
        // up-antiup
        lumi_entry![2, -2, 1.0; 4, -4, 1.0],
        // antiup-up
        lumi_entry![-2, 2, 1.0; -4, 4, 1.0],
        // down-antidown
        lumi_entry![1, -1, 1.0; 3, -3, 1.0; 5, -5, 1.0],
        // antidown-down
        lumi_entry![-1, 1, 1.0; -3, 3, 1.0; -5, 5, 1.0],
    ];

    let orders = vec![
        // LO
        Order {
            alphas: 0,
            alpha: 2,
            logxir: 0,
            logxif: 0,
        },
        // NLO QCD - won't be filled
        Order {
            alphas: 1,
            alpha: 2,
            logxir: 0,
            logxif: 0,
        },
        Order {
            alphas: 1,
            alpha: 2,
            logxir: 0,
            logxif: 1,
        },
    ];

    // we bin in rapidity from 0 to 2.4 in steps of 0.1
    let bin_limits: Vec<_> = (0..=24).map(|x| x as f64 / 10.0).collect();

    let mut subgrid_params = SubgridParams::default();
    let mut extra = ExtraSubgridParams::default();

    subgrid_params.set_q2_bins(30);
    subgrid_params.set_q2_max(1e6);
    subgrid_params.set_q2_min(1e2);
    subgrid_params.set_q2_order(3);
    subgrid_params.set_reweight(reweight);
    subgrid_params.set_x_bins(50);
    subgrid_params.set_x_max(1.0);
    subgrid_params.set_x_min(2e-7);
    subgrid_params.set_x_order(3);
    extra.set_x2_bins(50);
    extra.set_x2_max(1.0);
    extra.set_x2_min(2e-7);
    extra.set_x2_order(3);
    extra.set_reweight2(reweight);

    // create the PineAPPL grid
    let mut grid = Grid::with_subgrid_type(
        lumi,
        orders,
        bin_limits,
        subgrid_params,
        extra,
        subgrid_type,
    )?;

    // in GeV^2 pbarn
    let hbarc2 = 3.893793721e8;

    for _ in 0..calls {
        // generate a phase-space point
        let Psp2to2 {
            s,
            t,
            u,
            x1,
            x2,
            mut jacobian,
        } = hadronic_pspgen(rng, 10.0, 7000.0);

        let ptl = (t * u / s).sqrt();
        let mll = s.sqrt();
        let yll = 0.5 * (x1 / x2).ln();
        let ylp = (yll + (0.5 * mll / ptl).acosh()).abs();
        let ylm = (yll - (0.5 * mll / ptl).acosh()).abs();

        jacobian *= hbarc2 / (calls as f64);

        // cuts for LO for the invariant-mass slice containing the Z-peak from CMSDY2D11
        if (ptl < 14.0)
            || (yll.abs() > 2.4)
            || (ylp > 2.4)
            || (ylm > 2.4)
            || !(60.0..=120.0).contains(&mll)
        {
            continue;
        }

        let q2 = if dynamic { mll * mll } else { 90.0 * 90.0 };

        // LO photon-photon channel
        let weight = jacobian * int_photo(s, t, u);
        let pto = 0;
        let channel = 0;

        grid.fill(pto, yll.abs(), channel, &Ntuple { x1, x2, q2, weight });

        // LO up-antiup-type channel
        let weight = jacobian * int_quark(s, t, u, 2.0 / 3.0, 0.5);
        let pto = 0;
        let channel = 1;

        grid.fill(pto, yll.abs(), channel, &Ntuple { x1, x2, q2, weight });

        // LO antiup-up-type channel - swap (x1 <-> x2) and (t <-> u)
        let weight = jacobian * int_quark(s, u, t, 2.0 / 3.0, 0.5);
        let pto = 0;
        let channel = 2;

        grid.fill(pto, yll.abs(), channel, &Ntuple { x2, x1, q2, weight });

        // LO down-antidown-type channel
        let weight = jacobian * int_quark(s, t, u, -1.0 / 3.0, -0.5);
        let pto = 0;
        let channel = 3;

        grid.fill(pto, yll.abs(), channel, &Ntuple { x1, x2, q2, weight });

        // LO antidown-down-type channel - swap (x1 <-> x2) and (t <-> u)
        let weight = jacobian * int_quark(s, u, t, -1.0 / 3.0, -0.5);
        let pto = 0;
        let channel = 4;

        grid.fill(pto, yll.abs(), channel, &Ntuple { x2, x1, q2, weight });
    }

    Ok(grid)
}

fn perform_grid_tests(
    subgrid_type: &str,
    dynamic: bool,
    reference: &[f64],
    reference_after_ssd: &[f64],
    x_grid: &[f64],
    reweight: bool,
) -> Result<()> {
    let mut rng = Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);
    let mut grid = fill_drell_yan_lo_grid(&mut rng, 500_000, subgrid_type, dynamic, reweight)?;

    // TEST 1: `merge` and `scale`
    grid.merge(fill_drell_yan_lo_grid(
        &mut rng,
        500_000,
        subgrid_type,
        dynamic,
        reweight,
    )?)?;
    grid.scale(0.5);

    // suppress LHAPDF banners
    lhapdf::set_verbosity(0);

    let pdf_set = "NNPDF31_nlo_as_0118_luxqed";

    assert!(lhapdf::available_pdf_sets().iter().any(|x| x == pdf_set));

    let pdf = Pdf::with_setname_and_member(pdf_set, 0)?;
    let mut xfx = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut alphas = |_| 0.0;

    // TEST 2: `read` and `write`
    let mut file = Cursor::new(Vec::new());
    grid.write(&mut file)?;
    file.set_position(0);
    mem::drop(grid);
    let grid = Grid::read(&mut file)?;

    // TEST 3: `write_lz4`
    let mut file = Cursor::new(Vec::new());
    grid.write_lz4(&mut file)?;
    file.set_position(0);
    mem::drop(grid);
    let mut grid = Grid::read(&mut file)?;

    // TEST 4: `scale_by_order`
    grid.scale_by_order(10.0, 0.5, 10.0, 10.0, 1.0);
    grid.scale_by_order(10.0, 1.0, 10.0, 10.0, 4.0);

    // TEST 5: `convolute`
    let mut lumi_cache = LumiCache::with_one(2212, &mut xfx, &mut alphas);
    let bins = grid.convolute(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0)]);

    for (result, reference) in bins.iter().zip(reference.iter()) {
        assert_approx_eq!(f64, *result, *reference, ulps = 16);
    }

    // TEST 5b: `convolute` with `LumiCache::with_two`
    let mut xfx1 = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut xfx2 = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut alphas2 = |_| 0.0;
    let mut lumi_cache2 = LumiCache::with_two(2212, &mut xfx1, 2212, &mut xfx2, &mut alphas2);
    let bins2 = grid.convolute(&mut lumi_cache2, &[], &[], &[], &[(1.0, 1.0)]);

    for (result, reference) in bins2.iter().zip(reference.iter()) {
        assert_approx_eq!(f64, *result, *reference, ulps = 16);
    }

    mem::drop(lumi_cache2);
    mem::drop(bins2);

    // TEST 6: `convolute_subgrid`
    let bins: Vec<_> = (0..grid.bin_info().bins())
        .map(|bin| {
            (0..grid.lumi().len())
                .map(|channel| {
                    grid.convolute_subgrid(&mut lumi_cache, 0, bin, channel, 1.0, 1.0)
                        .sum()
                })
                .sum()
        })
        .collect();

    for (result, reference) in bins.iter().zip(reference.iter()) {
        assert_approx_eq!(f64, *result, *reference, ulps = 24);
    }

    // TEST 7a: `optimize_using` - tests `symmetrize` for each subgrid type
    grid.optimize_using(GridOptFlags::SYMMETRIZE_CHANNELS);

    // TEST 7b: `optimize`
    grid.optimize();

    assert_eq!(grid.subgrid(0, 0, 0).x1_grid().as_ref(), x_grid);
    assert_eq!(grid.subgrid(0, 0, 0).x2_grid().as_ref(), x_grid);

    // TEST 8: `convolute_subgrid` for the optimized subgrids
    let bins: Vec<_> = (0..grid.bin_info().bins())
        .map(|bin| {
            (0..grid.lumi().len())
                .map(|channel| {
                    grid.convolute_subgrid(&mut lumi_cache, 0, bin, channel, 1.0, 1.0)
                        .sum()
                })
                .sum()
        })
        .collect();

    for (result, reference_after_ssd) in bins.iter().zip(reference_after_ssd.iter()) {
        assert_approx_eq!(f64, *result, *reference_after_ssd, ulps = 24);
    }

    let bins = grid.convolute(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0)]);

    for (result, reference_after_ssd) in bins.iter().zip(reference_after_ssd.iter()) {
        assert_approx_eq!(f64, *result, *reference_after_ssd, ulps = 24);
    }

    // TEST 9: `set_remapper`

    // make a two-dimensional distribution out of it
    grid.set_remapper(BinRemapper::new(
        vec![0.1; 24],
        (0..24)
            .flat_map(|index| {
                let index = f64::from(index);
                vec![(60.0, 120.0), (index * 0.1, (index + 1.0) * 0.1)]
            })
            .collect::<Vec<(f64, f64)>>(),
    )?)?;

    // TEST 10: `merge_bins`

    // trivial merge: first bin is merged into first bin
    grid.merge_bins(0..1)?;

    for (result, reference_after_ssd) in bins.iter().zip(reference_after_ssd.iter()) {
        assert_approx_eq!(f64, *result, *reference_after_ssd, ulps = 24);
    }

    // merge two bins with each other
    for bin in 0..12 {
        grid.merge_bins(bin..bin + 2)?;
    }

    let merged2 = grid.convolute(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0)]);

    for (result, reference_after_ssd) in merged2.iter().zip(
        reference_after_ssd
            .chunks_exact(2)
            .map(|chunk| chunk.iter().sum::<f64>() / 2.0),
    ) {
        assert_approx_eq!(f64, *result, reference_after_ssd, ulps = 32);
    }

    // TEST 11: `delete_bins`

    // delete a few bins from the start
    grid.delete_bins(&[0, 1]);

    let deleted = grid.convolute(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0)]);

    assert_eq!(deleted.len(), 10);

    for (result, reference_after_ssd) in deleted.iter().zip(
        reference_after_ssd
            .chunks_exact(2)
            .map(|chunk| chunk.iter().sum::<f64>() / 2.0)
            .skip(2),
    ) {
        assert_approx_eq!(f64, *result, reference_after_ssd, ulps = 32);
    }

    // delete a few bins from the ending
    grid.delete_bins(&[8, 9]);

    let deleted2 = grid.convolute(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0)]);

    assert_eq!(deleted2.len(), 8);

    for (result, reference_after_ssd) in deleted2.iter().zip(
        reference_after_ssd
            .chunks_exact(2)
            .map(|chunk| chunk.iter().sum::<f64>() / 2.0)
            .skip(2)
            .take(6),
    ) {
        assert_approx_eq!(f64, *result, reference_after_ssd, ulps = 16);
    }

    Ok(())
}

fn generate_grid(subgrid_type: &str, dynamic: bool, reweight: bool) -> Result<Grid> {
    let mut rng = Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);
    fill_drell_yan_lo_grid(&mut rng, 500_000, subgrid_type, dynamic, reweight)
}

const STATIC_REFERENCE: [f64; 24] = [
    270.032948077419,
    266.53042891661994,
    290.18415971698766,
    257.18920263434893,
    238.1024311020346,
    297.4043699884724,
    260.05082619542054,
    239.86180454132813,
    242.8044656290245,
    231.94005409799632,
    224.546005767585,
    186.30387872772653,
    216.22804109717634,
    178.11604807980913,
    171.36239656080207,
    145.5738494644204,
    144.15628290390123,
    96.98346155489186,
    87.02664966484767,
    72.24178056556343,
    63.34788350688179,
    45.150098574761095,
    31.786673143908388,
    10.508579490820983,
];

// numbers are slightly different from `STATIC_REFERENCE` because the static scale detection (SSD)
// removes the Q^2 interpolation error
const STATIC_REFERENCE_AFTER_SSD: [f64; 24] = [
    270.033098990359,
    266.5305770813254,
    290.1843197471277,
    257.18934283371163,
    238.10255893783074,
    297.4045259189172,
    260.05095854131497,
    239.8619225312395,
    242.8045800917593,
    231.94015804152133,
    224.54610060526701,
    186.30395209734667,
    216.22811962577862,
    178.11610617693407,
    171.36244639833768,
    145.57388533569656,
    144.1563115005831,
    96.98347582894642,
    87.02665819940937,
    72.24178332992076,
    63.34788220165245,
    45.15009496754392,
    31.786668442718938,
    10.508577299430973,
];

const DYNAMIC_REFERENCE: [f64; 24] = [
    270.10748821350154,
    266.5411133114593,
    290.17831196525697,
    257.2279175882498,
    238.19001761554745,
    297.4774293868129,
    260.05099925772,
    239.8847562882569,
    242.8256504711324,
    231.93841366438,
    224.55846883004486,
    186.25680820137552,
    216.31250247378915,
    178.08658436751017,
    171.38330537962523,
    145.56521694168075,
    144.1389036069047,
    96.95384922315318,
    87.01668160609148,
    72.24037917990188,
    63.351600701502775,
    45.15589233137642,
    31.78895816514695,
    10.507668133626103,
];

const DYNAMIC_REFERENCE_NO_REWEIGHT: [f64; 24] = [
    269.0260807547221,
    265.6228062891378,
    289.20276155574334,
    256.20635007483935,
    237.32408404334382,
    296.4119556601695,
    259.1570057084532,
    239.0354655497739,
    242.00432881806023,
    231.0153726768181,
    223.77970088158807,
    185.74115061906738,
    215.5756627252798,
    177.31781780498287,
    170.85988900722407,
    145.28303184463408,
    143.54049829751125,
    96.55784208406595,
    86.82029114390203,
    72.05416364844432,
    63.125656746209586,
    45.02250707871406,
    31.699275851444902,
    10.491000144711,
];

#[test]
fn drell_yan_lagrange_static() -> Result<()> {
    perform_grid_tests(
        "LagrangeSubgrid",
        false,
        &STATIC_REFERENCE,
        &STATIC_REFERENCE_AFTER_SSD,
        &[
            0.030521584007828916,
            0.02108918668378717,
            0.014375068581090129,
            0.009699159574043398,
            0.006496206194633799,
            0.004328500638820811,
        ],
        true,
    )
}

#[test]
fn drell_yan_lagrange_v1_static() -> Result<()> {
    perform_grid_tests(
        "LagrangeSubgridV1",
        false,
        &STATIC_REFERENCE,
        &STATIC_REFERENCE, // LagrangeSubgridV1 doesn't have static-scale detection
        &[
            0.030521584007828916,
            0.02108918668378717,
            0.014375068581090129,
            0.009699159574043398,
            0.006496206194633799,
            0.004328500638820811,
        ],
        true,
    )
}

#[test]
fn drell_yan_lagrange_v2_static() -> Result<()> {
    perform_grid_tests(
        "LagrangeSubgridV2",
        false,
        &STATIC_REFERENCE,
        &STATIC_REFERENCE_AFTER_SSD,
        &[
            0.030521584007828916,
            0.02108918668378717,
            0.014375068581090129,
            0.009699159574043398,
            0.006496206194633799,
            0.004328500638820811,
        ],
        true,
    )
}

#[test]
fn drell_yan_lagrange_dynamic() -> Result<()> {
    perform_grid_tests(
        "LagrangeSubgrid",
        true,
        &DYNAMIC_REFERENCE,
        &DYNAMIC_REFERENCE,
        &[
            0.030521584007828916,
            0.02108918668378717,
            0.014375068581090129,
            0.009699159574043398,
            0.006496206194633799,
            0.004328500638820811,
        ],
        true,
    )
}

#[test]
fn drell_yan_lagrange_v1_dynamic() -> Result<()> {
    perform_grid_tests(
        "LagrangeSubgridV1",
        true,
        &DYNAMIC_REFERENCE,
        &DYNAMIC_REFERENCE,
        &[
            0.030521584007828916,
            0.02108918668378717,
            0.014375068581090129,
            0.009699159574043398,
            0.006496206194633799,
            0.004328500638820811,
        ],
        true,
    )
}

#[test]
fn drell_yan_lagrange_v1_dynamic_no_reweight() -> Result<()> {
    perform_grid_tests(
        "LagrangeSubgridV1",
        true,
        &DYNAMIC_REFERENCE_NO_REWEIGHT,
        &DYNAMIC_REFERENCE_NO_REWEIGHT,
        &[
            0.030521584007828916,
            0.02108918668378717,
            0.014375068581090129,
            0.009699159574043398,
            0.006496206194633799,
            0.004328500638820811,
        ],
        false,
    )
}

#[test]
fn drell_yan_lagrange_v2_dynamic() -> Result<()> {
    perform_grid_tests(
        "LagrangeSubgridV2",
        true,
        &DYNAMIC_REFERENCE,
        &DYNAMIC_REFERENCE,
        &[
            0.030521584007828916,
            0.02108918668378717,
            0.014375068581090129,
            0.009699159574043398,
            0.006496206194633799,
            0.004328500638820811,
        ],
        true,
    )
}

#[test]
fn drell_yan_lagrange_v2_dynamic_no_reweight() -> Result<()> {
    perform_grid_tests(
        "LagrangeSubgridV2",
        true,
        &DYNAMIC_REFERENCE_NO_REWEIGHT,
        &DYNAMIC_REFERENCE_NO_REWEIGHT,
        &[
            0.030521584007828916,
            0.02108918668378717,
            0.014375068581090129,
            0.009699159574043398,
            0.006496206194633799,
            0.004328500638820811,
        ],
        false,
    )
}

#[test]
fn drell_yan_lagrange_sparse_dynamic() -> Result<()> {
    perform_grid_tests(
        "LagrangeSparseSubgrid",
        true,
        &DYNAMIC_REFERENCE,
        &DYNAMIC_REFERENCE,
        &[
            0.030521584007828916,
            0.02108918668378717,
            0.014375068581090129,
            0.009699159574043398,
            0.006496206194633799,
            0.004328500638820811,
        ],
        true,
    )
}

#[test]
fn grid_optimize() -> Result<()> {
    let mut grid = generate_grid("LagrangeSubgridV2", false, false)?;

    assert_eq!(grid.orders().len(), 3);
    assert_eq!(grid.lumi().len(), 5);
    assert!(matches!(
        grid.subgrids()[[0, 0, 0]],
        SubgridEnum::LagrangeSubgridV2 { .. }
    ));
    assert_eq!(grid.subgrids()[[0, 0, 0]].x1_grid().len(), 50);
    assert_eq!(grid.subgrids()[[0, 0, 0]].x2_grid().len(), 50);
    assert_eq!(grid.subgrids()[[0, 0, 0]].mu2_grid().len(), 30);

    let mut grid2 = grid.clone();
    grid2.optimize_using(GridOptFlags::OPTIMIZE_SUBGRID_TYPE);

    // `OPTIMIZE_SUBGRID_TYPE` changes the subgrid type ...
    assert!(matches!(
        grid2.subgrids()[[0, 0, 0]],
        SubgridEnum::ImportOnlySubgridV2 { .. }
    ));
    // and the dimensions of the subgrid
    assert_eq!(grid2.subgrids()[[0, 0, 0]].x1_grid().len(), 6);
    assert_eq!(grid2.subgrids()[[0, 0, 0]].x2_grid().len(), 6);
    assert_eq!(grid2.subgrids()[[0, 0, 0]].mu2_grid().len(), 4);

    grid.optimize_using(GridOptFlags::OPTIMIZE_SUBGRID_TYPE | GridOptFlags::STATIC_SCALE_DETECTION);

    assert!(matches!(
        grid.subgrids()[[0, 0, 0]],
        SubgridEnum::ImportOnlySubgridV2 { .. }
    ));
    // if `STATIC_SCALE_DETECTION` is present the `mu2_grid` dimension are better optimized
    assert_eq!(grid.subgrids()[[0, 0, 0]].x1_grid().len(), 6);
    assert_eq!(grid.subgrids()[[0, 0, 0]].x2_grid().len(), 6);
    assert_eq!(grid.subgrids()[[0, 0, 0]].mu2_grid().len(), 1);

    // has no effect for this test
    grid.optimize_using(GridOptFlags::SYMMETRIZE_CHANNELS);

    assert_eq!(grid.orders().len(), 3);
    assert_eq!(grid.lumi().len(), 5);

    grid.optimize_using(GridOptFlags::STRIP_EMPTY_ORDERS);

    assert_eq!(grid.orders().len(), 1);
    assert_eq!(grid.lumi().len(), 5);

    // has no effect for this test
    grid.optimize_using(GridOptFlags::MERGE_SAME_CHANNELS);

    assert_eq!(grid.orders().len(), 1);
    assert_eq!(grid.lumi().len(), 5);

    grid.optimize_using(GridOptFlags::STRIP_EMPTY_CHANNELS);

    assert_eq!(grid.orders().len(), 1);
    assert_eq!(grid.lumi().len(), 3);

    Ok(())
}
