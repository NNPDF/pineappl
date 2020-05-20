use num_complex::Complex64;
use pineappl::grid::{Grid, Ntuple, Order, SubgridParams};
use pineappl::lumi_entry;
use rand::Rng;
use rand_pcg::Pcg64;
use std::f64::consts::{PI, SQRT_2};

// If equation numbers are given, they are from Max Huber's PhD thesis:
//   'Radiative corrections to the neutral-current Drell-Yan process'

// Eq. (2.8)
fn int_qqbar(s: f64, t: f64, u: f64, qq: f64, i3wq: f64) -> f64 {
    // lepton charges
    let ql = -1.0;
    let i3wl = -0.5;

    // input parameters
    let gf = 1.16637e-5;
    let mw = 80.403;
    let mz = 91.1876;
    let gw = 2.141;
    let gz = 2.4952;

    // derived quantities
    let muw2 = Complex64::new(mw * mw, -gw * mw);
    let muz2 = Complex64::new(mz * mz, -gz * mz);
    let cw2 = muw2 / muz2;
    let sw2 = 1.0 - cw2;
    let cw = cw2.sqrt();
    let sw = sw2.sqrt();
    let alphagf = SQRT_2 / PI * gf * mw * mw * (1.0 - (mw / mz).powi(2));

    // couplings and abbreviations
    let gqqzp = -sw / cw * qq;
    let gllzp = -sw / cw * ql;
    let gqqzm = (i3wq - sw * sw * qq) / (sw * cw);
    let gllzm = (i3wl - sw * sw * ql) / (sw * cw);
    let chiz = s / (s - muz2);

    alphagf.powi(2) / 12.0 / s.powi(3)
        * (2.0 * qq * qq * ql * ql * (t * t + u * u)
            + 2.0
                * qq
                * ql
                * (((gqqzp * gllzp + gqqzm * gllzm) * u * u
                    + (gqqzp * gllzm + gqqzm * gllzp) * t * t)
                    * chiz)
                    .re
            + (((gqqzp * gllzp).norm_sqr() + (gqqzm * gllzm).norm_sqr()) * u * u
                + ((gqqzp * gllzm).norm_sqr() + (gqqzm * gllzp).norm_sqr()) * t * t)
                * chiz.norm_sqr())
}

// Eq. (2.9)
fn int_photo(s: f64, t: f64, u: f64) -> f64 {
    let alpha0: f64 = 1.0 / 137.03599911;
    alpha0.powi(2) / 2.0 / s * (t / u + u / t)
}

struct Psp2to2 {
    s: f64,
    t: f64,
    u: f64,
    ptl: f64,
    mll: f64,
    yll: f64,
    x1: f64,
    x2: f64,
    jacobian: f64,
}

fn hadronic_pspgen(rng: &mut impl Rng, mmin: f64, mmax: f64) -> Psp2to2 {
    let smin = mmin * mmin;
    let smax = mmax * mmax;

    let mz = 91.1876;
    let gz = 2.4952;

    let mut jacobian = 1.0;

    // Eq. (D.8)
    let gmin = ((smin - mz * mz) / (gz * mz)).atan();
    let gmax = ((smax - mz * mz) / (gz * mz)).atan();
    // Eq. (D.7)
    let s = mz * mz + gz * mz * (gmin + (gmax - gmin) * rng.gen::<f64>()).tan();
    // Eq. (D.9)
    jacobian *= (gmin - gmax) / (gz * mz) * ((s - mz * mz).powi(2) + gz * gz * mz * mz);

    // transformation of variables from (x1,x2) to (s,x2)
    let x2 = (s / smax).powf(rng.gen::<f64>());
    let x1 = s / (smax * x2);
    jacobian *= x2 * (smax / s).ln();

    // self-consistency checks
    assert!(x1 >= 0.0);
    assert!(x1 < 1.0);
    assert!(x2 >= 0.0);
    assert!(x2 < 1.0);
    assert!((smax * x1 * x2 / s).max(s / (smax * x1 * x2)) <= 1.0 + 8.0 * std::f64::EPSILON);

    // theta integration (in the CMS)
    let cos_theta = 2.0 * rng.gen::<f64>() - 1.0;
    jacobian *= 2.0;

    let t = -0.5 * s * (1.0 - cos_theta);
    let u = -0.5 * s * (1.0 + cos_theta);

    // self-consistency checks
    assert!(t >= -s);
    assert!(u >= -s);

    // phi integration
    jacobian *= 2.0 * PI;

    let ptl = (t * u / s).sqrt();
    let mll = s.sqrt();
    let yll = (x1 / x2).ln();

    // self-consistency check
    assert!(ptl <= 0.5 * s.sqrt());

    Psp2to2 {
        s,
        t,
        u,
        ptl,
        mll,
        yll,
        x1,
        x2,
        jacobian,
    }
}

fn fill_drell_yan_lo_grid(rng: &mut impl Rng, calls: usize) -> Grid {
    let lumi = vec![
        // down-pair flavors
        lumi_entry![1, -1, 1.0; 3, -3, 1.0; 5, -5, 1.0],
        lumi_entry![-1, 1, 1.0; -3, 3, 1.0; -5, 5, 1.0],
        // up-pair flavor
        lumi_entry![2, -2, 1.0; 4, -4, 1.0],
        lumi_entry![-2, 2, 1.0; -4, 4, 1.0],
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

    for _ in 0..calls {
        // generate a phase-space point
        let Psp2to2 {
            s,
            t,
            u,
            ptl,
            mll,
            yll,
            x1,
            x2,
            jacobian,
        } = hadronic_pspgen(rng, 60.0, 120.0);

        // these cuts should always be fulfilled by the choice of `smin` and `smax` above
        assert!(mll >= 60.0);
        assert!(mll <= 120.0);

        // cuts for LO for the invariant-mass slice containing the Z-peak from CMSDY2D11
        if (ptl < 14.0) || (yll > 2.4) {
            continue;
        }

        let weights = [
            // down-pair flavors
            jacobian * int_qqbar(s, t, u, -0.5, -1.0 / 3.0),
            // down-pair flavors with quarks swapped (t <-> u)
            jacobian * int_qqbar(s, u, t, -0.5, -1.0 / 3.0),
            // up-pair flavors
            jacobian * int_qqbar(s, t, u, 0.5, 2.0 / 3.0),
            // up-pair flavors with quarks swapped (t <-> u)
            jacobian * int_qqbar(s, u, t, 0.5, 2.0 / 3.0),
            // photon pair
            jacobian * int_photo(s, u, t),
        ];

        grid.fill_all(
            0,
            yll,
            &Ntuple {
                x1,
                x2,
                q2: 91.1876_f64.powi(2),
                weight: (),
            },
            &weights,
        );
    }

    grid
}

#[test]
fn dy_fill_speed_1000() {
    fill_drell_yan_lo_grid(
        &mut Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96),
        1000,
    );

    // TODO:
    // - check if the output is meaningful
    // - compare against mg5_aMC
    // - add more tests
}
