#![allow(missing_docs)]

use anyhow::Result;
use float_cmp::assert_approx_eq;
use itertools::izip;
use lhapdf::Pdf;
use num_complex::Complex;
use pineappl::bin::BinRemapper;
use pineappl::boc::{Kinematics, Order, ScaleFuncForm, Scales};
use pineappl::channel;
use pineappl::convolutions::{Conv, ConvType, ConvolutionCache};
use pineappl::grid::{Grid, GridOptFlags};
use pineappl::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use pineappl::pids::PidBasis;
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

// ALLOW: in this example we care more about readability than floating-point accuracy
#[allow(clippy::suboptimal_flops)]
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
    let cos_theta = rng.gen::<f64>().mul_add(2.0, -1.0);
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

fn fill_drell_yan_lo_grid(rng: &mut impl Rng, calls: u32, dynamic: bool, reweight: bool) -> Grid {
    let channels = vec![
        // photons
        channel![1.0 * (22, 22)],
        // up-antiup
        channel![1.0 * (2, -2) + 1.0 * (4, -4)],
        // antiup-up
        channel![1.0 * (-2, 2) + 1.0 * (-4, 4)],
        // down-antidown
        channel![1.0 * (1, -1) + 1.0 * (3, -3) + 1.0 * (5, -5)],
        // antidown-down
        channel![1.0 * (-1, 1) + 1.0 * (-3, 3) + 1.0 * (-5, 5)],
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
    let convolutions = vec![
        Conv::new(ConvType::UnpolPDF, 2212),
        Conv::new(ConvType::UnpolPDF, 2212),
    ];

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
        Kinematics::Scale(0),
        // 2nd dimension is the parton momentum fraction of the first convolution
        Kinematics::X(0),
        // 3rd dimension is the parton momentum fraction of the second convolution
        Kinematics::X(1),
    ];

    let scales = Scales {
        ren: ScaleFuncForm::Scale(0),
        fac: ScaleFuncForm::Scale(0),
        frg: ScaleFuncForm::NoScale,
    };

    // create the PineAPPL grid
    let mut grid = Grid::new(
        // the integers in the channel definition are PDG Monte Carlo IDs
        PidBasis::Pdg,
        channels,
        orders,
        bin_limits,
        convolutions,
        interps,
        kinematics,
        scales,
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

    grid
}

fn perform_grid_tests(
    dynamic: bool,
    reference: &[f64],
    reference_after_ssd: &[f64],
    x_grid: &[f64],
    reweight: bool,
) -> Result<()> {
    let mut rng = Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);
    let mut grid = fill_drell_yan_lo_grid(&mut rng, INT_STATS, dynamic, reweight);

    // TEST 1: `merge` and `scale`
    grid.merge(fill_drell_yan_lo_grid(
        &mut rng, INT_STATS, dynamic, reweight,
    ))?;
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
    grid.scale_by_order(10.0, 0.5, 10.0, 10.0, 1.0, 1.0);
    grid.scale_by_order(10.0, 1.0, 10.0, 10.0, 1.0, 4.0);

    // TEST 5: `convolve`
    let mut convolution_cache = ConvolutionCache::new(
        vec![Conv::new(ConvType::UnpolPDF, 2212)],
        vec![&mut xfx],
        &mut alphas,
    );
    let bins = grid.convolve(&mut convolution_cache, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

    for (result, reference) in bins.iter().zip(reference.iter()) {
        assert_approx_eq!(f64, *result, *reference, ulps = 4);
    }

    // TEST 5b: `convolve` with `ConvolutionCache::with_two`
    let mut xfx1 = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut xfx2 = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut alphas2 = |_| 0.0;
    let mut convolution_cache2 = ConvolutionCache::new(
        vec![
            Conv::new(ConvType::UnpolPDF, 2212),
            Conv::new(ConvType::UnpolPDF, 2212),
        ],
        vec![&mut xfx1, &mut xfx2],
        &mut alphas2,
    );
    let bins2 = grid.convolve(&mut convolution_cache2, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

    for (result, reference) in bins2.iter().zip(reference.iter()) {
        assert_approx_eq!(f64, *result, *reference, ulps = 16);
    }

    mem::drop(convolution_cache2);
    mem::drop(bins2);

    // TEST 6: `convolve_subgrid`
    let bins: Vec<_> = (0..grid.bin_info().bins())
        .map(|bin| {
            (0..grid.channels().len())
                .map(|channel| {
                    grid.convolve_subgrid(&mut convolution_cache, 0, bin, channel, (1.0, 1.0, 1.0))
                        .sum()
                })
                .sum()
        })
        .collect();

    for (result, reference) in bins.iter().zip(reference.iter()) {
        assert_approx_eq!(f64, *result, *reference, ulps = 16);
    }

    // TEST 7a: `optimize_using` - tests `symmetrize` for each subgrid type
    grid.optimize_using(GridOptFlags::SYMMETRIZE_CHANNELS);

    // TEST 7b: `optimize`
    grid.optimize();

    let node_values = grid.subgrids()[[0, 0, 0]].node_values();

    for (&node_value1, &node_value2, &ref_value) in izip!(&node_values[1], &node_values[2], x_grid)
    {
        assert_approx_eq!(f64, node_value1, ref_value, ulps = 4);
        assert_approx_eq!(f64, node_value2, ref_value, ulps = 4);
    }

    // TEST 8: `convolve_subgrid` for the optimized subgrids
    let bins: Vec<_> = (0..grid.bin_info().bins())
        .map(|bin| {
            (0..grid.channels().len())
                .map(|channel| {
                    grid.convolve_subgrid(&mut convolution_cache, 0, bin, channel, (1.0, 1.0, 1.0))
                        .sum()
                })
                .sum()
        })
        .collect();

    for (result, reference_after_ssd) in bins.iter().zip(reference_after_ssd.iter()) {
        assert_approx_eq!(f64, *result, *reference_after_ssd, ulps = 32);
    }

    let bins = grid.convolve(&mut convolution_cache, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

    for (result, reference_after_ssd) in bins.iter().zip(reference_after_ssd.iter()) {
        assert_approx_eq!(f64, *result, *reference_after_ssd, ulps = 64);
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
        assert_approx_eq!(f64, *result, *reference_after_ssd, ulps = 64);
    }

    // merge two bins with each other
    for bin in 0..12 {
        grid.merge_bins(bin..bin + 2)?;
    }

    let merged2 = grid.convolve(&mut convolution_cache, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

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

    let deleted = grid.convolve(&mut convolution_cache, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

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

    let deleted2 = grid.convolve(&mut convolution_cache, &[], &[], &[], &[(1.0, 1.0, 1.0)]);

    assert_eq!(deleted2.len(), 8);

    for (result, reference_after_ssd) in deleted2.iter().zip(
        reference_after_ssd
            .chunks_exact(2)
            .map(|chunk| chunk.iter().sum::<f64>() / 2.0)
            .skip(2)
            .take(6),
    ) {
        assert_approx_eq!(f64, *result, reference_after_ssd, ulps = 32);
    }

    Ok(())
}

fn generate_grid(dynamic: bool, reweight: bool) -> Grid {
    let mut rng = Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);
    fill_drell_yan_lo_grid(&mut rng, 500_000, dynamic, reweight)
}

// number is small enough for the tests to be quick and meaningful
const INT_STATS: u32 = 10_000;

const STATIC_REFERENCE: [f64; 24] = [
    168.38959322031457,
    513.0768026978171,
    371.71287396366995,
    235.8098362992929,
    186.20124775877665,
    216.2304514988581,
    232.04909213040415,
    170.06035328939964,
    103.33958880717813,
    191.87061398627344,
    182.48313944454654,
    171.6188420512471,
    419.82434919583875,
    151.14152991889875,
    226.05268712169035,
    199.64591790851492,
    369.8935560707606,
    77.06480847483719,
    24.490861801031517,
    12.19000587591836,
    23.110914187129627,
    87.06841942116057,
    33.0281593619898,
    2.362076425563692,
];

// numbers are slightly different from `STATIC_REFERENCE` because the node optimization detects the
// static scale and removes the Q^2 interpolation error
const STATIC_REFERENCE_AFTER_SSD: [f64; 24] = [
    168.38968474125505,
    513.0770882860825,
    371.7130781560312,
    235.80996443941294,
    186.20134644917403,
    216.2305662654171,
    232.04921269665587,
    170.060435777966,
    103.33963708862501,
    191.87070884467244,
    182.48320743597665,
    171.6189069431743,
    419.82450657418957,
    151.14157394136424,
    226.0527454508727,
    199.64596561172203,
    369.8936241814945,
    77.06482253107433,
    24.49086417366925,
    12.190006763665671,
    23.110914180153948,
    87.0684069821286,
    33.028158528350204,
    2.3620759284762887,
];

const DYNAMIC_REFERENCE: [f64; 24] = [
    167.73075899059722,
    513.3569347058141,
    371.3336398828405,
    235.53004623228287,
    185.89154709174852,
    216.75969380255526,
    231.6876677545833,
    169.67997661655352,
    103.56259318132902,
    191.68748388958508,
    182.00660765469212,
    171.6221591151596,
    420.1122395090359,
    150.81647073655915,
    226.0168315024308,
    199.99645663613614,
    370.5740387393846,
    77.13459227149951,
    24.570096511403573,
    12.15059986236983,
    23.04588396906037,
    87.2298678460105,
    33.05576050516689,
    2.3061332535775794,
];

const DYNAMIC_REFERENCE_NO_REWEIGHT: [f64; 24] = [
    167.02560507500328,
    511.6307601313951,
    370.1659470550191,
    234.51187692755803,
    185.23716936219668,
    215.91820762194604,
    230.92167039044742,
    169.0157392153116,
    103.19948220956648,
    190.9153768744078,
    181.49749277120952,
    171.15900435357685,
    418.69191536323956,
    150.25485756172543,
    225.25882657556969,
    199.5592968215913,
    369.0301660240426,
    76.82053478041635,
    24.477528662138134,
    12.114147164297055,
    22.967201699163894,
    86.99454515845706,
    32.894968080392196,
    2.2999405116169878,
];

const X_GRID: [f64; 6] = [
    0.030521584007828877,
    0.02108918668378717,
    0.014375068581090129,
    0.009699159574043398,
    0.006496206194633799,
    0.004328500638819831,
];

#[test]
fn drell_yan_static() -> Result<()> {
    perform_grid_tests(
        false,
        &STATIC_REFERENCE,
        &STATIC_REFERENCE_AFTER_SSD,
        &X_GRID,
        true,
    )
}

#[test]
fn drell_yan_dynamic() -> Result<()> {
    perform_grid_tests(true, &DYNAMIC_REFERENCE, &DYNAMIC_REFERENCE, &X_GRID, true)
}

#[test]
fn drell_yan_dynamic_no_reweight() -> Result<()> {
    perform_grid_tests(
        true,
        &DYNAMIC_REFERENCE_NO_REWEIGHT,
        &DYNAMIC_REFERENCE_NO_REWEIGHT,
        &X_GRID,
        false,
    )
}

#[test]
fn grid_optimize() {
    let mut grid = generate_grid(false, false);

    assert_eq!(grid.orders().len(), 3);
    assert_eq!(grid.channels().len(), 5);
    assert!(matches!(
        grid.subgrids()[[0, 0, 0]],
        SubgridEnum::InterpSubgridV1 { .. }
    ));

    let node_values = grid.subgrids()[[0, 0, 0]].node_values();
    assert_eq!(node_values[0].len(), 30);
    assert_eq!(node_values[1].len(), 50);
    assert_eq!(node_values[2].len(), 50);

    let mut grid2 = grid.clone();
    grid2.optimize_using(GridOptFlags::OPTIMIZE_SUBGRID_TYPE);

    // `OPTIMIZE_SUBGRID_TYPE` changes the subgrid type ...
    assert!(matches!(
        grid2.subgrids()[[0, 0, 0]],
        SubgridEnum::ImportSubgridV1 { .. }
    ));
    // and the dimensions of the subgrid
    let node_values = grid2.subgrids()[[0, 0, 0]].node_values();
    assert_eq!(node_values[0].len(), 4);
    assert_eq!(node_values[1].len(), 6);
    assert_eq!(node_values[2].len(), 6);

    grid.optimize_using(GridOptFlags::OPTIMIZE_NODES);

    assert!(matches!(
        grid.subgrids()[[0, 0, 0]],
        SubgridEnum::InterpSubgridV1 { .. }
    ));
    let node_values = grid.subgrids()[[0, 0, 0]].node_values();
    assert_eq!(node_values[0].len(), 1);
    assert_eq!(node_values[1].len(), 6);
    assert_eq!(node_values[2].len(), 6);

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
}
