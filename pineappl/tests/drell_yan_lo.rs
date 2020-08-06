use float_cmp::approx_eq;
use lhapdf::Pdf;
use pineappl::grid::{Grid, Ntuple, Order, SubgridParams};
use pineappl::lumi_entry;
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
) -> anyhow::Result<Grid> {
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
    let bin_limits: Vec<f64> = (0..=24).map(|x| x as f64 / 10.0).collect();

    let mut subgrid_params = SubgridParams::default();

    subgrid_params.set_q2_bins(30);
    subgrid_params.set_q2_max(1e6);
    subgrid_params.set_q2_min(1e2);
    subgrid_params.set_q2_order(3);
    subgrid_params.set_reweight(true);
    subgrid_params.set_x_bins(50);
    subgrid_params.set_x_max(1.0);
    subgrid_params.set_x_min(2e-7);
    subgrid_params.set_x_order(3);

    // create the PineAPPL grid
    let mut grid = Grid::with_subgrid_type(lumi, orders, bin_limits, subgrid_params, subgrid_type)?;

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
            || (mll < 60.0)
            || (mll > 120.0)
        {
            continue;
        }

        let weight = jacobian * int_photo(s, u, t);

        grid.fill(
            0,
            yll.abs(),
            0,
            &Ntuple {
                x1,
                x2,
                q2: 90.0_f64.powi(2),
                weight: weight,
            },
        );
    }

    Ok(grid)
}

#[test]
fn dy_aa_lagrange_subgrid() -> anyhow::Result<()> {
    let mut rng = Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);
    let mut grid = fill_drell_yan_lo_grid(&mut rng, 500_000, "LagrangeSubgrid")?;

    grid.merge(fill_drell_yan_lo_grid(
        &mut rng,
        500_000,
        "LagrangeSubgrid",
    )?)?;
    grid.scale(0.5);

    let pdf_set = "NNPDF31_nlo_as_0118_luxqed";

    assert!(lhapdf::available_pdf_sets().iter().any(|x| x == &pdf_set));

    let pdf = Pdf::with_setname_and_member(&pdf_set, 0);
    let xfx = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let alphas = |_| 0.0;

    // check `read` and `write`
    let mut file = Cursor::new(Vec::new());
    grid.write(&mut file)?;
    file.set_position(0);
    mem::drop(grid);
    let mut grid = Grid::read(&mut file)?;

    // some useless scalings
    grid.scale_by_order(10.0, 0.5, 10.0, 10.0, 1.0);
    grid.scale_by_order(10.0, 1.0, 10.0, 10.0, 4.0);

    let bins = grid.convolute(&xfx, &xfx, &alphas, &[], &[], &[], &[(1.0, 1.0)]);
    let reference = vec![
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

    for (result, reference) in bins.iter().zip(reference.iter()) {
        assert!(approx_eq!(f64, *result, *reference, ulps = 16));
    }

    Ok(())
}

#[test]
fn dy_aa_ntuple_subgrid() -> anyhow::Result<()> {
    let mut rng = Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);
    let mut grid = fill_drell_yan_lo_grid(&mut rng, 500_000, "NtupleSubgrid")?;

    grid.merge(fill_drell_yan_lo_grid(&mut rng, 500_000, "NtupleSubgrid")?)?;
    grid.scale(0.5);

    let pdf_set = "NNPDF31_nlo_as_0118_luxqed";

    assert!(lhapdf::available_pdf_sets().iter().any(|x| x == &pdf_set));

    let pdf = Pdf::with_setname_and_member(&pdf_set, 0);
    let xfx = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let alphas = |_| 0.0;

    // check `read` and `write`
    let mut file = Cursor::new(Vec::new());
    grid.write(&mut file)?;
    file.set_position(0);
    mem::drop(grid);
    let mut grid = Grid::read(&mut file)?;

    // some useless scalings
    grid.scale_by_order(10.0, 0.5, 10.0, 10.0, 1.0);
    grid.scale_by_order(10.0, 1.0, 10.0, 10.0, 4.0);

    let bins = grid.convolute(&xfx, &xfx, &alphas, &[], &[], &[], &[(1.0, 1.0)]);
    let reference = vec![
        5.294102629460027e-1,
        5.407462109703145e-1,
        5.696654765115949e-1,
        5.028229825084712e-1,
        4.917145014197304e-1,
        4.94698595872653e-1,
        4.9190141387427233e-1,
        4.486371523100685e-1,
        4.507732384369339e-1,
        4.1060709945028706e-1,
        3.6024761439304853e-1,
        3.2757335386076103e-1,
        2.7928032956476956e-1,
        2.4983016484405093e-1,
        2.1077131839702393e-1,
        1.7797159986885916e-1,
        1.541103341646812e-1,
        1.195781437556192e-1,
        9.39928943589368e-2,
        6.718800239822122e-2,
        5.13639097541644e-2,
        3.571458571056351e-2,
        2.0671050816841553e-2,
        7.299744833167185e-3,
    ];

    for (result, reference) in bins.iter().zip(reference.iter()) {
        assert!(approx_eq!(f64, *result, *reference, ulps = 16));
    }

    Ok(())
}
