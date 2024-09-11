use anyhow::Result;
use float_cmp::assert_approx_eq;
use lhapdf::Pdf;
use num_complex::Complex;
use pineappl::bin::BinRemapper;
use pineappl::boc::{Kinematics, Order};
use pineappl::channel;
use pineappl::convolutions::{Convolution, LumiCache};
use pineappl::grid::{Grid, GridOptFlags};
use pineappl::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use pineappl::subgrid::{Subgrid, SubgridEnum};
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
    calls: u32,
    dynamic: bool,
    reweight: bool,
) -> Result<Grid> {
    let channels = vec![
        // photons
        channel![22, 22, 1.0],
        // up-antiup
        channel![2, -2, 1.0; 4, -4, 1.0],
        // antiup-up
        channel![-2, 2, 1.0; -4, 4, 1.0],
        // down-antidown
        channel![1, -1, 1.0; 3, -3, 1.0; 5, -5, 1.0],
        // antidown-down
        channel![-1, 1, 1.0; -3, 3, 1.0; -5, 5, 1.0],
    ];

    let orders = vec![
        // LO
        Order {
            alphas: 0,
            alpha: 2,
            logxir: 0,
            logxif: 0,
            logxia: 0,
        },
        // NLO QCD - won't be filled
        Order {
            alphas: 1,
            alpha: 2,
            logxir: 0,
            logxif: 0,
            logxia: 0,
        },
        Order {
            alphas: 1,
            alpha: 2,
            logxir: 0,
            logxif: 1,
            logxia: 0,
        },
    ];

    // we bin in rapidity from 0 to 2.4 in steps of 0.1
    let bin_limits: Vec<_> = (0..=24).map(|x: u32| f64::from(x) / 10.0).collect();

    // the grid represents data with two unpolarized proton PDFs
    let convolutions = vec![Convolution::UnpolPDF(2212), Convolution::UnpolPDF(2212)];

    let reweight = if reweight {
        ReweightMeth::ApplGridX
    } else {
        ReweightMeth::NoReweight
    };

    // define how `Grid::fill` interpolates
    let interps = vec![
        // 1st dimension interpolation parameters
        Interp::new(
            1e2,
            1e6,
            30,
            3,
            ReweightMeth::NoReweight,
            Map::ApplGridH0,
            InterpMeth::Lagrange,
        ),
        // 2nd dimension interpolation parameters
        Interp::new(
            2e-7,
            1.0,
            50,
            3,
            reweight,
            Map::ApplGridF2,
            InterpMeth::Lagrange,
        ),
        // 3rd dimension interpolation parameters
        Interp::new(
            2e-7,
            1.0,
            50,
            3,
            reweight,
            Map::ApplGridF2,
            InterpMeth::Lagrange,
        ),
    ];

    let kinematics = vec![
        // 1st dimension is factorization and at the same time also the renormalization scale
        Kinematics::MU2_RF,
        // 2nd dimension is the parton momentum fraction of the first convolution
        Kinematics::X1,
        // 3rd dimension is the parton momentum fraction of the second convolution
        Kinematics::X2,
    ];

    // create the PineAPPL grid
    let mut grid = Grid::new(
        channels,
        orders,
        bin_limits,
        convolutions,
        interps,
        kinematics,
    );

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

        jacobian *= hbarc2 / f64::from(calls);

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

        grid.fill(pto, yll.abs(), channel, &[q2, x1, x2], weight);

        // LO up-antiup-type channel
        let weight = jacobian * int_quark(s, t, u, 2.0 / 3.0, 0.5);
        let pto = 0;
        let channel = 1;

        grid.fill(pto, yll.abs(), channel, &[q2, x1, x2], weight);

        // LO antiup-up-type channel - swap (x1 <-> x2) and (t <-> u)
        let weight = jacobian * int_quark(s, u, t, 2.0 / 3.0, 0.5);
        let pto = 0;
        let channel = 2;

        grid.fill(pto, yll.abs(), channel, &[q2, x2, x1], weight);

        // LO down-antidown-type channel
        let weight = jacobian * int_quark(s, t, u, -1.0 / 3.0, -0.5);
        let pto = 0;
        let channel = 3;

        grid.fill(pto, yll.abs(), channel, &[q2, x1, x2], weight);

        // LO antidown-down-type channel - swap (x1 <-> x2) and (t <-> u)
        let weight = jacobian * int_quark(s, u, t, -1.0 / 3.0, -0.5);
        let pto = 0;
        let channel = 4;

        grid.fill(pto, yll.abs(), channel, &[q2, x2, x1], weight);
    }

    Ok(grid)
}

fn perform_grid_tests(
    dynamic: bool,
    reference: &[f64],
    reference_after_ssd: &[f64],
    x_grid: &[f64],
    reweight: bool,
) -> Result<()> {
    let mut rng = Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);
    let mut grid = fill_drell_yan_lo_grid(&mut rng, INT_STATS, dynamic, reweight)?;

    // TEST 1: `merge` and `scale`
    grid.merge(fill_drell_yan_lo_grid(
        &mut rng, INT_STATS, dynamic, reweight,
    )?)?;
    grid.scale(0.5);

    // suppress LHAPDF banners
    lhapdf::set_verbosity(0);

    let pdf_set = "NNPDF31_nlo_as_0118_luxqed";

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

    // TEST 5: `convolve`
    let mut lumi_cache = LumiCache::with_one(2212, &mut xfx, &mut alphas);
    let bins = grid.convolve(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

    for (result, reference) in bins.iter().zip(reference.iter()) {
        assert_approx_eq!(f64, *result, *reference, ulps = 16);
    }

    // TEST 5b: `convolve` with `LumiCache::with_two`
    let mut xfx1 = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut xfx2 = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut alphas2 = |_| 0.0;
    let mut lumi_cache2 = LumiCache::with_two(2212, &mut xfx1, 2212, &mut xfx2, &mut alphas2);
    let bins2 = grid.convolve(&mut lumi_cache2, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

    for (result, reference) in bins2.iter().zip(reference.iter()) {
        assert_approx_eq!(f64, *result, *reference, ulps = 16);
    }

    mem::drop(lumi_cache2);
    mem::drop(bins2);

    // TEST 6: `convolve_subgrid`
    let bins: Vec<_> = (0..grid.bin_info().bins())
        .map(|bin| {
            (0..grid.channels().len())
                .map(|channel| {
                    grid.convolve_subgrid(&mut lumi_cache, 0, bin, channel, (1.0, 1.0, 1.0))
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

    assert_eq!(grid.subgrids()[[0, 0, 0]].x1_grid().as_ref(), x_grid);
    assert_eq!(grid.subgrids()[[0, 0, 0]].x2_grid().as_ref(), x_grid);

    // TEST 8: `convolve_subgrid` for the optimized subgrids
    let bins: Vec<_> = (0..grid.bin_info().bins())
        .map(|bin| {
            (0..grid.channels().len())
                .map(|channel| {
                    grid.convolve_subgrid(&mut lumi_cache, 0, bin, channel, (1.0, 1.0, 1.0))
                        .sum()
                })
                .sum()
        })
        .collect();

    for (result, reference_after_ssd) in bins.iter().zip(reference_after_ssd.iter()) {
        assert_approx_eq!(f64, *result, *reference_after_ssd, ulps = 24);
    }

    let bins = grid.convolve(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

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

    let merged2 = grid.convolve(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

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

    let deleted = grid.convolve(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

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

    let deleted2 = grid.convolve(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

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

fn generate_grid(dynamic: bool, reweight: bool) -> Result<Grid> {
    let mut rng = Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);
    fill_drell_yan_lo_grid(&mut rng, 500_000, dynamic, reweight)
}

// number is small enough for the tests to be quick and meaningful
const INT_STATS: u32 = 10_000;

const STATIC_REFERENCE: [f64; 24] = [
    168.3895932203146,
    513.0768026978171,
    371.7128739636697,
    235.8098362992926,
    186.2012477587763,
    216.23045149885752,
    232.04909213040406,
    170.06035328940015,
    103.33958880717893,
    191.8706139862758,
    182.48313944455572,
    171.6188420512558,
    419.8243491958518,
    151.1415299189055,
    226.05268712169604,
    199.6459179085148,
    369.89355607075987,
    77.06480847483716,
    24.490861801031464,
    12.19000587591836,
    23.1109141871296,
    87.06841942116225,
    33.02815936198823,
    2.362076425563693,
];

// numbers are slightly different from `STATIC_REFERENCE` because the static scale detection (SSD)
// removes the Q^2 interpolation error
const STATIC_REFERENCE_AFTER_SSD: [f64; 24] = [
    168.38968474125483,
    513.0770882860817,
    371.7130781560301,
    235.8099644394127,
    186.20134644917377,
    216.23056626541657,
    232.0492126966559,
    170.06043577796652,
    103.3396370886257,
    191.87070884467457,
    182.48320743598543,
    171.61890694318237,
    419.824506574201,
    151.14157394137038,
    226.05274545087767,
    199.64596561172144,
    369.89362418149375,
    77.0648225310743,
    24.490864173669205,
    12.190006763665675,
    23.110914180153898,
    87.06840698213031,
    33.02815852834857,
    2.3620759284762887,
];

const DYNAMIC_REFERENCE: [f64; 24] = [
    167.73075899059725,
    513.3569347058141,
    371.33363988284026,
    235.5300462322826,
    185.89154709174824,
    216.75969380255472,
    231.6876677545832,
    169.6799766165541,
    103.56259318132975,
    191.6874838895875,
    182.00660765470136,
    171.62215911516836,
    420.11223950904895,
    150.81647073656592,
    226.01683150243653,
    199.99645663613603,
    370.57403873938387,
    77.1345922714995,
    24.570096511403527,
    12.150599862369834,
    23.045883969060334,
    87.22986784601223,
    33.0557605051653,
    2.306133253577581,
];

const DYNAMIC_REFERENCE_NO_REWEIGHT: [f64; 24] = [
    167.02560507500345,
    511.63076013139545,
    370.16594705502035,
    234.5118769275643,
    185.2371693622032,
    215.91820762195567,
    230.92167039044898,
    169.015739215293,
    103.19948220954701,
    190.91537687434288,
    181.4974927711281,
    171.15900435348723,
    418.69191536302435,
    150.25485756167274,
    225.25882657551898,
    199.5592968215911,
    369.03016602405376,
    76.82053478041806,
    24.477528662138404,
    12.114147164297156,
    22.96720169916392,
    86.9945451584532,
    32.89496808038666,
    2.2999405116169536,
];

#[test]
fn drell_yan_static() -> Result<()> {
    perform_grid_tests(
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
fn drell_yan_dynamic() -> Result<()> {
    perform_grid_tests(
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
fn drell_yan_dynamic_no_reweight() -> Result<()> {
    perform_grid_tests(
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
fn grid_optimize() -> Result<()> {
    let mut grid = generate_grid(false, false)?;

    assert_eq!(grid.orders().len(), 3);
    assert_eq!(grid.channels().len(), 5);
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
        SubgridEnum::PackedQ1X2SubgridV1 { .. }
    ));
    // and the dimensions of the subgrid
    assert_eq!(grid2.subgrids()[[0, 0, 0]].x1_grid().len(), 6);
    assert_eq!(grid2.subgrids()[[0, 0, 0]].x2_grid().len(), 6);
    assert_eq!(grid2.subgrids()[[0, 0, 0]].mu2_grid().len(), 4);

    grid.optimize_using(GridOptFlags::OPTIMIZE_SUBGRID_TYPE | GridOptFlags::STATIC_SCALE_DETECTION);

    assert!(matches!(
        grid.subgrids()[[0, 0, 0]],
        SubgridEnum::PackedQ1X2SubgridV1 { .. }
    ));
    // if `STATIC_SCALE_DETECTION` is present the `mu2_grid` dimension are better optimized
    assert_eq!(grid.subgrids()[[0, 0, 0]].x1_grid().len(), 6);
    assert_eq!(grid.subgrids()[[0, 0, 0]].x2_grid().len(), 6);
    assert_eq!(grid.subgrids()[[0, 0, 0]].mu2_grid().len(), 1);

    // has no effect for this test
    grid.optimize_using(GridOptFlags::SYMMETRIZE_CHANNELS);

    assert_eq!(grid.orders().len(), 3);
    assert_eq!(grid.channels().len(), 5);

    grid.optimize_using(GridOptFlags::STRIP_EMPTY_ORDERS);

    assert_eq!(grid.orders().len(), 1);
    assert_eq!(grid.channels().len(), 5);

    // has no effect for this test
    grid.optimize_using(GridOptFlags::MERGE_SAME_CHANNELS);

    assert_eq!(grid.orders().len(), 1);
    assert_eq!(grid.channels().len(), 5);

    grid.optimize_using(GridOptFlags::STRIP_EMPTY_CHANNELS);

    assert_eq!(grid.orders().len(), 1);
    assert_eq!(grid.channels().len(), 3);

    Ok(())
}
