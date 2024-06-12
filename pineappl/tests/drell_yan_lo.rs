use anyhow::Result;
use float_cmp::assert_approx_eq;
use lhapdf::Pdf;
use num_complex::Complex;
use pineappl::bin::BinRemapper;
use pineappl::boc::Order;
use pineappl::channel;
use pineappl::convolutions::LumiCache;
use pineappl::grid::{Grid, GridOptFlags, Ntuple};
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
    calls: u32,
    subgrid_type: &str,
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
    let bin_limits: Vec<_> = (0..=24).map(|x: u32| f64::from(x) / 10.0).collect();

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
        channels,
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

        grid.fill(
            pto,
            yll.abs(),
            channel,
            &Ntuple {
                x1: x2,
                x2: x1,
                q2,
                weight,
            },
        );

        // LO down-antidown-type channel
        let weight = jacobian * int_quark(s, t, u, -1.0 / 3.0, -0.5);
        let pto = 0;
        let channel = 3;

        grid.fill(pto, yll.abs(), channel, &Ntuple { x1, x2, q2, weight });

        // LO antidown-down-type channel - swap (x1 <-> x2) and (t <-> u)
        let weight = jacobian * int_quark(s, u, t, -1.0 / 3.0, -0.5);
        let pto = 0;
        let channel = 4;

        grid.fill(
            pto,
            yll.abs(),
            channel,
            &Ntuple {
                x1: x2,
                x2: x1,
                q2,
                weight,
            },
        );
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
    let bins = grid.convolve(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0)]);

    for (result, reference) in bins.iter().zip(reference.iter()) {
        assert_approx_eq!(f64, *result, *reference, ulps = 16);
    }

    // TEST 5b: `convolve` with `LumiCache::with_two`
    let mut xfx1 = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut xfx2 = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut alphas2 = |_| 0.0;
    let mut lumi_cache2 = LumiCache::with_two(2212, &mut xfx1, 2212, &mut xfx2, &mut alphas2);
    let bins2 = grid.convolve(&mut lumi_cache2, &[], &[], &[], &[(1.0, 1.0)]);

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
                    grid.convolve_subgrid(&mut lumi_cache, 0, bin, channel, 1.0, 1.0)
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
                    grid.convolve_subgrid(&mut lumi_cache, 0, bin, channel, 1.0, 1.0)
                        .sum()
                })
                .sum()
        })
        .collect();

    for (result, reference_after_ssd) in bins.iter().zip(reference_after_ssd.iter()) {
        assert_approx_eq!(f64, *result, *reference_after_ssd, ulps = 24);
    }

    let bins = grid.convolve(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0)]);

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

    let merged2 = grid.convolve(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0)]);

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

    let deleted = grid.convolve(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0)]);

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

    let deleted2 = grid.convolve(&mut lumi_cache, &[], &[], &[], &[(1.0, 1.0)]);

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
    269.89225458312495,
    266.2168804878282,
    290.0467478314624,
    258.0064918266305,
    239.54186548997865,
    300.17541324377703,
    258.8811221515799,
    238.4064950360576,
    242.5494601562957,
    236.34329830221077,
    230.63243720020898,
    190.03118557029666,
    213.22241277258763,
    177.75582251643334,
    168.07022695390958,
    151.59217101220256,
    143.81017491485716,
    97.09707327367487,
    91.38465432190982,
    73.94464862425771,
    63.859689262732104,
    48.595785504299926,
    27.94818010640803,
    9.343737799674852,
];

// numbers are slightly different from `STATIC_REFERENCE` because the static scale detection (SSD)
// removes the Q^2 interpolation error
const STATIC_REFERENCE_AFTER_SSD: [f64; 24] = [
    269.89240546283145,
    266.2170285827742,
    290.04690782935967,
    258.0066322019259,
    239.54199362567599,
    300.17556967146095,
    258.88125430161745,
    238.40661279174125,
    242.54957458220744,
    236.34340283622035,
    230.63253265929194,
    190.03125927151245,
    213.2224910582812,
    177.7558806305883,
    168.07027678254747,
    151.59220685502618,
    143.81020355582885,
    97.09708758263099,
    91.38466242593998,
    73.94465114837278,
    63.859687905917,
    48.595781165174515,
    27.94817639459665,
    9.343735959243446,
];

const DYNAMIC_REFERENCE: [f64; 24] = [
    269.9662650413552,
    266.2274509325408,
    290.039119030095,
    258.04801305108583,
    239.63561020879277,
    300.2475932636636,
    258.88126161648313,
    238.42709542929794,
    242.5724521248901,
    236.3541498865422,
    230.64832146047578,
    189.999243811704,
    213.2896760201295,
    177.7280865940876,
    168.0886178280483,
    151.59285700593935,
    143.80051106343882,
    97.0715765765853,
    91.38479915098559,
    73.94713838892906,
    63.85622547082087,
    48.61296466751912,
    27.948404940991445,
    9.342761664545428,
];

const DYNAMIC_REFERENCE_NO_REWEIGHT: [f64; 24] = [
    268.8874311488598,
    265.3130436782233,
    289.0614714145284,
    257.02578172672656,
    238.76378338813032,
    299.1756333696102,
    257.98748703027104,
    237.58099891213897,
    241.75215319366012,
    235.41757682699438,
    229.8671307486547,
    189.47964517011536,
    212.56055728623704,
    176.9591711445695,
    167.56523215346917,
    151.30532185043768,
    143.20366078799765,
    96.67453775369947,
    91.18334210163036,
    73.75879631942671,
    63.629606742074984,
    48.47126745674977,
    27.86328933386428,
    9.32654010506528,
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
