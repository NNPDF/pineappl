use anyhow::Result;
use float_cmp::assert_approx_eq;
use lhapdf::Pdf;
use pineappl::bin::BinRemapper;
use pineappl::grid::{Grid, Ntuple, Order};
use pineappl::lumi::LumiCache;
use pineappl::lumi_entry;
use pineappl::subgrid::{ExtraSubgridParams, Subgrid, SubgridParams};
use rand::Rng;
use rand_pcg::Pcg64;
use std::f64::consts::PI;
use std::io::Cursor;
use std::mem;

// If equation numbers are given, they are from Max Huber's PhD thesis:
//   'Radiative corrections to the neutral-current Drell-Yan process'

// Eq. (2.9) - gamma-gamma contribution to DY pair production
fn int_photo(s: f64, t: f64, u: f64) -> f64 {
    let alpha0: f64 = 1.0 / 137.03599911;
    alpha0.powi(2) / 2.0 / s * (t / u + u / t)
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
    ];

    // only LO alpha^2
    let orders = vec![Order {
        alphas: 0,
        alpha: 2,
        logxir: 0,
        logxif: 0,
    }];

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

        let weight = jacobian * int_photo(s, u, t);
        let q2 = if dynamic { mll * mll } else { 90.0 * 90.0 };

        grid.fill(0, yll.abs(), 0, &Ntuple { x1, x2, q2, weight });
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
            grid.convolute_subgrid(&mut lumi_cache, 0, bin, 0, 1.0, 1.0)
                .sum()
        })
        .collect();

    for (result, reference) in bins.iter().zip(reference.iter()) {
        assert_approx_eq!(f64, *result, *reference, ulps = 24);
    }

    // TEST 7: `optimize`
    grid.optimize();

    assert_eq!(grid.subgrid(0, 0, 0).x1_grid().as_ref(), x_grid);
    assert_eq!(grid.subgrid(0, 0, 0).x2_grid().as_ref(), x_grid);

    // TEST 8: `convolute_subgrid` for the optimized subgrids
    let bins: Vec<_> = (0..grid.bin_info().bins())
        .map(|bin| {
            grid.convolute_subgrid(&mut lumi_cache, 0, bin, 0, 1.0, 1.0)
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

const STATIC_REFERENCE: [f64; 24] = [
    5.29438499470369e-1,
    5.407794857747981e-1,
    5.696902048421408e-1,
    5.028293125183927e-1,
    4.91700010813869e-1,
    4.946801648085449e-1,
    4.9188982902741324e-1,
    4.486342876707396e-1,
    4.5078095484478575e-1,
    4.106209790738961e-1,
    3.602582665914635e-1,
    3.275794060602094e-1,
    2.7928887723972295e-1,
    2.498545969916396e-1,
    2.108027399225483e-1,
    1.7799404895027734e-1,
    1.5411875388722898e-1,
    1.1957877908479132e-1,
    9.39935306988927e-2,
    6.719034949888109e-2,
    5.136619446064035e-2,
    3.5716156871884834e-2,
    2.067251421406746e-2,
    7.300411258011377e-3,
];

// numbers are slightly different from `STATIC_REFERENCE` because the static scale detection (SSD)
// removes the Q^2 interpolation error
const STATIC_REFERENCE_AFTER_SSD: [f64; 24] = [
    5.2943850298637385e-1,
    5.4077949153675209e-1,
    5.6969021381443985e-1,
    5.0282932188764318e-1,
    4.9170002525140877e-1,
    4.9468018558433052e-1,
    4.918898576307103e-1,
    4.4863432083118493e-1,
    4.5078099648970371e-1,
    4.1062102518688581e-1,
    3.602583146986339e-1,
    3.2757945699703028e-1,
    2.7928892704491243e-1,
    2.4985464964922635e-1,
    2.1080278995465099e-1,
    1.7799409692247298e-1,
    1.5411880069615053e-1,
    1.1957881962307169e-1,
    9.3993565751843353e-2,
    6.7190377322845676e-2,
    5.1366217946639786e-2,
    3.5716174780312381e-2,
    2.0672525560241309e-2,
    7.3004155738883077e-3,
];

const DYNAMIC_REFERENCE: [f64; 24] = [
    5.093090431949207e-1,
    5.191668797562395e-1,
    5.467930909144878e-1,
    4.845000146366871e-1,
    4.7279670245192884e-1,
    4.7481525429115445e-1,
    4.7060007393610065e-1,
    4.286009571276975e-1,
    4.2997427448811987e-1,
    3.911121604214181e-1,
    3.430494924416351e-1,
    3.1204501485640446e-1,
    2.6656633773109056e-1,
    2.3745772514449243e-1,
    2.0055415662593068e-1,
    1.692229208614707e-1,
    1.4635749293124972e-1,
    1.1348804559115269e-1,
    8.934261002200886e-2,
    6.385631160399488e-2,
    4.866968827833745e-2,
    3.379282482821324e-2,
    1.9494720220040448e-2,
    6.880478270711603e-3,
];

const DYNAMIC_REFERENCE_NO_REWEIGHT: [f64; 24] = [
    5.07858012983822e-1,
    5.175284452774326e-1,
    5.450720052517568e-1,
    4.831265401289246e-1,
    4.714770822390454e-1,
    4.7332381643255367e-1,
    4.6911980551007526e-1,
    4.274224757508723e-1,
    4.2874742764343926e-1,
    3.898631579204727e-1,
    3.419807843167341e-1,
    3.112011822819245e-1,
    2.6582805074531396e-1,
    2.3669405399362703e-1,
    1.999186753344747e-1,
    1.6877228654478332e-1,
    1.4593203448521183e-1,
    1.1312893772255903e-1,
    8.906983809157609e-2,
    6.368194887638147e-2,
    4.852988118458943e-2,
    3.369252022300261e-2,
    1.943724748356927e-2,
    6.8603903480790075e-3,
];

#[test]
fn dy_aa_lagrange_static() -> Result<()> {
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
fn dy_aa_lagrange_v1_static() -> Result<()> {
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
fn dy_aa_lagrange_v2_static() -> Result<()> {
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
fn dy_aa_lagrange_dynamic() -> Result<()> {
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
fn dy_aa_lagrange_v1_dynamic() -> Result<()> {
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
fn dy_aa_lagrange_v1_dynamic_no_reweight() -> Result<()> {
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
fn dy_aa_lagrange_v2_dynamic() -> Result<()> {
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
fn dy_aa_lagrange_v2_dynamic_no_reweight() -> Result<()> {
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
fn dy_aa_lagrange_sparse_dynamic() -> Result<()> {
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
