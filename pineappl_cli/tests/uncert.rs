use assert_cmd::Command;
use std::num::NonZeroUsize;
use std::thread;

const HELP_STR: &str = "Calculates scale and PDF uncertainties

Usage: pineappl uncert [OPTIONS] <--pdf|--scale-abs[=<SCALES>]|--scale-cov[=<SCALES>]|--scale-env[=<SCALES>]> <INPUT> <PDFSET>

Arguments:
  <INPUT>   Path to the input grid
  <PDFSET>  LHAPDF id or name of the PDF set

Options:
      --pdf                   Calculate the PDF uncertainties
      --scale-abs[=<SCALES>]  Show absolute numbers of the scale-varied results [possible values: 3, 7, 9]
      --scale-cov[=<SCALES>]  Calculate scale uncertainties using the covariance method [possible values: 3, 7, 9]
      --scale-env[=<SCALES>]  Calculate the envelope of results where renormalization and factorization scales varied [possible values: 3, 7, 9]
      --cl <CL>               Confidence level in per cent, for PDF uncertainties [default: 68.26894921370858]
  -i, --integrated            Show integrated numbers (without bin widths) instead of differential ones
  -o, --orders <ORDERS>       Select orders manually
      --threads <THREADS>     Number of threads to utilize [default: {}]
      --digits-abs <ABS>      Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>      Set the number of fractional digits shown for relative numbers [default: 2]
  -h, --help                  Print help
";

const DEFAULT_STR: &str = "b   etal    dsig/detal  PDF central    PDF    
     []        [pb]                    [%]    
-+----+----+-----------+-----------+-----+----
0    2 2.25 7.5459110e2 7.5461655e2 -1.14 1.14
1 2.25  2.5 6.9028342e2 6.9027941e2 -1.16 1.16
2  2.5 2.75 6.0025198e2 6.0022595e2 -1.18 1.18
3 2.75    3 4.8552235e2 4.8548211e2 -1.22 1.22
4    3 3.25 3.6195456e2 3.6191001e2 -1.27 1.27
5 3.25  3.5 2.4586691e2 2.4582640e2 -1.35 1.35
6  3.5    4 1.1586851e2 1.1584074e2 -1.51 1.51
7    4  4.5 2.7517266e1 2.7504644e1 -2.77 2.77
";

const CL_90_STR: &str = "b   etal    dsig/detal  PDF central    PDF    
     []        [pb]                    [%]    
-+----+----+-----------+-----------+-----+----
0    2 2.25 7.5459110e2 7.5461655e2 -1.87 1.87
1 2.25  2.5 6.9028342e2 6.9027941e2 -1.90 1.90
2  2.5 2.75 6.0025198e2 6.0022595e2 -1.95 1.95
3 2.75    3 4.8552235e2 4.8548211e2 -2.00 2.00
4    3 3.25 3.6195456e2 3.6191001e2 -2.08 2.08
5 3.25  3.5 2.4586691e2 2.4582640e2 -2.22 2.22
6  3.5    4 1.1586851e2 1.1584074e2 -2.48 2.48
7    4  4.5 2.7517266e1 2.7504644e1 -4.55 4.55
";

const INTEGRATED_STR: &str = "b   etal       integ    PDF central    PDF    
     []         []                     [%]    
-+----+----+-----------+-----------+-----+----
0    2 2.25 1.8864777e2 1.8865414e2 -1.14 1.14
1 2.25  2.5 1.7257086e2 1.7256985e2 -1.16 1.16
2  2.5 2.75 1.5006300e2 1.5005649e2 -1.18 1.18
3 2.75    3 1.2138059e2 1.2137053e2 -1.22 1.22
4    3 3.25 9.0488640e1 9.0477502e1 -1.27 1.27
5 3.25  3.5 6.1466727e1 6.1456599e1 -1.35 1.35
6  3.5    4 5.7934254e1 5.7920368e1 -1.51 1.51
7    4  4.5 1.3758633e1 1.3752322e1 -2.77 2.77
";

const ORDERS_A2_AS1A2_STR: &str = "b   etal    dsig/detal  PDF central    PDF    
     []        [pb]                    [%]    
-+----+----+-----------+-----------+-----+----
0    2 2.25 7.6246034e2 7.6248591e2 -1.14 1.14
1 2.25  2.5 6.9684577e2 6.9684166e2 -1.16 1.16
2  2.5 2.75 6.0548681e2 6.0546059e2 -1.18 1.18
3 2.75    3 4.8928139e2 4.8924093e2 -1.22 1.22
4    3 3.25 3.6454175e2 3.6449702e2 -1.27 1.27
5 3.25  3.5 2.4754316e2 2.4750254e2 -1.35 1.35
6  3.5    4 1.1667878e2 1.1665095e2 -1.50 1.50
7    4  4.5 2.7737493e1 2.7724826e1 -2.77 2.77
";

const SCALE_ABS_STR: &str =
"b   etal    dsig/detal   (r=1,f=1)   (r=2,f=2)  (r=0.5,f=0.5)  (r=2,f=1)   (r=1,f=2)  (r=0.5,f=1) (r=1,f=0.5)
     []        [pb]        [pb]        [pb]         [pb]         [pb]        [pb]        [pb]        [pb]    
-+----+----+-----------+-----------+-----------+-------------+-----------+-----------+-----------+-----------
0    2 2.25 7.5459110e2 7.5459110e2 7.6745431e2   7.4296019e2 7.4384068e2 7.7529764e2 7.6796494e2 7.2595107e2
1 2.25  2.5 6.9028342e2 6.9028342e2 7.0221920e2   6.7923774e2 6.8058382e2 7.0957480e2 7.0235002e2 6.6417441e2
2  2.5 2.75 6.0025198e2 6.0025198e2 6.1056383e2   5.9046454e2 5.9160658e2 6.1750966e2 6.1100712e2 5.7747295e2
3 2.75    3 4.8552235e2 4.8552235e2 4.9366919e2   4.7761552e2 4.7841022e2 4.9966237e2 4.9437007e2 4.6723687e2
4    3 3.25 3.6195456e2 3.6195456e2 3.6783089e2   3.5611822e2 3.5652780e2 3.7261436e2 3.6870561e2 3.4843600e2
5 3.25  3.5 2.4586691e2 2.4586691e2 2.4967698e2   2.4198770e2 2.4207028e2 2.5316566e2 2.5059003e2 2.3677625e2
6  3.5    4 1.1586851e2 1.1586851e2 1.1746280e2   1.1418227e2 1.1396174e2 1.1930157e2 1.1824058e2 1.1166942e2
7    4  4.5 2.7517266e1 2.7517266e1 2.7787333e1   2.7211003e1 2.7002241e1 2.8306905e1 2.8157972e1 2.6562471e1
";

const SCALE_COV_STR: &str = "b   etal    dsig/detal  7pt scale (cov)
     []        [pb]           [%]      
-+----+----+-----------+-------+-------
0    2 2.25 7.5459110e2   -3.29    3.29
1 2.25  2.5 6.9028342e2   -3.30    3.30
2  2.5 2.75 6.0025198e2   -3.34    3.34
3 2.75    3 4.8552235e2   -3.35    3.35
4    3 3.25 3.6195456e2   -3.35    3.35
5 3.25  3.5 2.4586691e2   -3.34    3.34
6  3.5    4 1.1586851e2   -3.31    3.31
7    4  4.5 2.7517266e1   -3.24    3.24
";

const SCALE_COV_9_STR: &str = "b   etal    dsig/detal  9pt scale (cov)
     []        [pb]           [%]      
-+----+----+-----------+-------+-------
0    2 2.25 7.5459110e2   -4.48    4.48
1 2.25  2.5 6.9028342e2   -4.48    4.48
2  2.5 2.75 6.0025198e2   -4.55    4.55
3 2.75    3 4.8552235e2   -4.57    4.57
4    3 3.25 3.6195456e2   -4.59    4.59
5 3.25  3.5 2.4586691e2   -4.61    4.61
6  3.5    4 1.1586851e2   -4.64    4.64
7    4  4.5 2.7517266e1   -4.67    4.67
";

const SCALE_ENV_STR: &str = "b   etal    dsig/detal  7pt-svar (env) 
     []        [pb]           [%]      
-+----+----+-----------+-------+-------
0    2 2.25 7.5459110e2   -3.80    2.74
1 2.25  2.5 6.9028342e2   -3.78    2.79
2  2.5 2.75 6.0025198e2   -3.79    2.88
3 2.75    3 4.8552235e2   -3.77    2.91
4    3 3.25 3.6195456e2   -3.73    2.95
5 3.25  3.5 2.4586691e2   -3.70    2.97
6  3.5    4 1.1586851e2   -3.62    2.96
7    4  4.5 2.7517266e1   -3.47    2.87
";

const SCALE_ENV_9_STR: &str = "b   etal    dsig/detal  9pt-svar (env) 
     []        [pb]           [%]      
-+----+----+-----------+-------+-------
0    2 2.25 7.5459110e2   -5.61    4.04
1 2.25  2.5 6.9028342e2   -5.54    4.12
2  2.5 2.75 6.0025198e2   -5.53    4.31
3 2.75    3 4.8552235e2   -5.48    4.45
4    3 3.25 3.6195456e2   -5.44    4.59
5 3.25  3.5 2.4586691e2   -5.40    4.73
6  3.5    4 1.1586851e2   -5.37    4.94
7    4  4.5 2.7517266e1   -5.36    5.22
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["uncert", "--help"])
        .assert()
        .success()
        .stdout(
            HELP_STR.replace(
                "{}",
                &thread::available_parallelism()
                    .map_or(1, NonZeroUsize::get)
                    .to_string(),
            ),
        );
}

#[test]
fn pdf_default() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "uncert",
            "--pdf",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}

#[test]
fn pdf_cl_90() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "uncert",
            "--pdf",
            "--cl=90",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(CL_90_STR);
}

#[test]
fn pdf_integrated() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "uncert",
            "--pdf",
            "--integrated",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(INTEGRATED_STR);
}

#[test]
fn pdf_orders_a2_as1a2() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "uncert",
            "--pdf",
            "--orders=a2,as1a2",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ORDERS_A2_AS1A2_STR);
}

#[test]
fn scale_abs() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "uncert",
            "--scale-abs",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_ABS_STR);
}

#[test]
fn scale_cov() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "uncert",
            "--scale-cov",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_COV_STR);
}

#[test]
fn scale_cov_9() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "uncert",
            "--scale-cov=9",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_COV_9_STR);
}

#[test]
fn scale_env() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "uncert",
            "--scale-env",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_ENV_STR);
}

#[test]
fn scale_env_9() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "uncert",
            "--scale-env=9",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_ENV_9_STR);
}
