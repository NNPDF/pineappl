use float_cmp::approx_eq;
use lhapdf::Pdf;
use pineappl::grid::{Grid, Ntuple, Order, SubgridParams};
use pineappl::lumi_entry;
use rand::Rng;
use rand_pcg::Pcg64;
use std::f64::consts::PI;

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

fn fill_drell_yan_lo_grid(rng: &mut impl Rng, calls: usize) -> Grid {
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

    // create the PineAPPL grid
    let mut grid = Grid::new(lumi, orders, bin_limits, SubgridParams::default());

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

    grid
}

#[test]
fn dy_aa_convolute_with_nnpdf31_luxqed() {
    let mut rng = Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);
    let grid = fill_drell_yan_lo_grid(&mut rng, 1_000_000);
    let pdf_set = "NNPDF31_nlo_as_0118_luxqed";

    assert!(lhapdf::available_pdf_sets().iter().any(|x| x == &pdf_set));
    let pdf = Pdf::with_setname_and_member(&pdf_set, 0);
    let xfx = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let alphas = |_| 0.0;

    let bins = grid.convolute(&xfx, &xfx, &alphas, &[], &[], &[], &[(1.0, 1.0)]);
    let reference = vec![
        5.279171684656772e-1,
        5.390429345674902e-1,
        5.678642509809768e-1,
        5.013919472729697e-1,
        4.90315459341301e-1,
        4.930991418324766e-1,
        4.9031523058660703e-1,
        4.473911555096199e-1,
        4.4948130313360196e-1,
        4.092807117343194e-1,
        3.591135225778676e-1,
        3.266887202413869e-1,
        2.7850659808559064e-1,
        2.4902890147317383e-1,
        2.1012250445479977e-1,
        1.7751910054649375e-1,
        1.5366372491678373e-1,
        1.1919102476022085e-1,
        9.370244098781093e-2,
        6.700596617428176e-2,
        5.1216211477073045e-2,
        3.560765252853136e-2,
        2.0610595970607874e-2,
        7.278729457396478e-3,
    ];

    for (result, reference) in bins.iter().zip(reference.iter()) {
        assert!(approx_eq!(f64, *result, *reference, ulps = 16));
    }
}
