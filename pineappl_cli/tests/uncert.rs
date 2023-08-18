use assert_cmd::Command;
use std::num::NonZeroUsize;
use std::thread;

const HELP_STR: &str = "Calculates scale and PDF uncertainties

Usage: pineappl uncert [OPTIONS] <--pdf|--scale-abs[=<SCALES>]|--scale-env[=<SCALES>]> <INPUT> <PDFSET>

Arguments:
  <INPUT>   Path to the input grid
  <PDFSET>  LHAPDF id or name of the PDF set

Options:
      --pdf                   Calculate the PDF uncertainties
      --scale-abs[=<SCALES>]  Show absolute numbers of the scale-varied results [possible values: 3, 7, 9]
      --scale-env[=<SCALES>]  Calculate the envelope of results where renormalization and factorization scales varied [possible values: 3, 7, 9]
      --cl <CL>               Confidence level in per cent, for PDF uncertainties [default: 68.26894921370858]
  -i, --integrated            Show integrated numbers (without bin widths) instead of differential ones
  -o, --orders <ORDERS>       Select orders manually
      --threads <THREADS>     Number of threads to utilize [default: {}]
      --digits-abs <ABS>      Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>      Set the number of fractional digits shown for relative numbers [default: 2]
  -h, --help                  Print help
";

const DEFAULT_STR: &str = "b   etal    disg/detal  PDF central    PDF    
     []        [pb]                    [%]    
-+----+----+-----------+-----------+-----+----
0    2 2.25 3.7527620e2 3.7528868e2 -1.14 1.14
1 2.25  2.5 3.4521553e2 3.4521365e2 -1.16 1.16
2  2.5 2.75 3.0001406e2 3.0000102e2 -1.18 1.18
3 2.75    3 2.4257663e2 2.4255656e2 -1.22 1.22
4    3 3.25 1.8093343e2 1.8091118e2 -1.27 1.27
5 3.25  3.5 1.2291115e2 1.2289094e2 -1.35 1.35
6  3.5    4 5.7851018e1 5.7837137e1 -1.50 1.50
7    4  4.5 1.3772029e1 1.3765722e1 -2.76 2.76
";

const CL_90_STR: &str = "b   etal    disg/detal  PDF central    PDF    
     []        [pb]                    [%]    
-+----+----+-----------+-----------+-----+----
0    2 2.25 3.7527620e2 3.7528868e2 -1.88 1.88
1 2.25  2.5 3.4521553e2 3.4521365e2 -1.90 1.90
2  2.5 2.75 3.0001406e2 3.0000102e2 -1.95 1.95
3 2.75    3 2.4257663e2 2.4255656e2 -2.00 2.00
4    3 3.25 1.8093343e2 1.8091118e2 -2.08 2.08
5 3.25  3.5 1.2291115e2 1.2289094e2 -2.22 2.22
6  3.5    4 5.7851018e1 5.7837137e1 -2.48 2.48
7    4  4.5 1.3772029e1 1.3765722e1 -4.54 4.54
";

const INTEGRATED_STR: &str = "b   etal       integ    PDF central    PDF    
     []         []                     [%]    
-+----+----+-----------+-----------+-----+----
0    2 2.25 9.3819050e1 9.3822169e1 -1.14 1.14
1 2.25  2.5 8.6303882e1 8.6303411e1 -1.16 1.16
2  2.5 2.75 7.5003515e1 7.5000256e1 -1.18 1.18
3 2.75    3 6.0644157e1 6.0639140e1 -1.22 1.22
4    3 3.25 4.5233358e1 4.5227794e1 -1.27 1.27
5 3.25  3.5 3.0727788e1 3.0722735e1 -1.35 1.35
6  3.5    4 2.8925509e1 2.8918568e1 -1.50 1.50
7    4  4.5 6.8860146e0 6.8828610e0 -2.76 2.76
";

const ORDERS_A2_AS1A2_STR: &str = "b   etal    disg/detal  PDF central    PDF    
     []        [pb]                    [%]    
-+----+----+-----------+-----------+-----+----
0    2 2.25 3.7918224e2 3.7919477e2 -1.14 1.14
1 2.25  2.5 3.4849530e2 3.4849336e2 -1.16 1.16
2  2.5 2.75 3.0262287e2 3.0260975e2 -1.18 1.18
3 2.75    3 2.4443923e2 2.4441905e2 -1.22 1.22
4    3 3.25 1.8224461e2 1.8222226e2 -1.26 1.26
5 3.25  3.5 1.2371575e2 1.2369548e2 -1.35 1.35
6  3.5    4 5.8295643e1 5.8281739e1 -1.50 1.50
7    4  4.5 1.3881516e1 1.3875186e1 -2.77 2.77
";

const SCALE_ABS_STR: &str =
"b   etal    disg/detal   (r=1,f=1)   (r=2,f=2)  (r=0.5,f=0.5)  (r=2,f=1)   (r=1,f=2)  (r=0.5,f=1) (r=1,f=0.5)
     []        [pb]        [pb]        [pb]         [pb]         [pb]        [pb]        [pb]        [pb]    
-+----+----+-----------+-----------+-----------+-------------+-----------+-----------+-----------+-----------
0    2 2.25 3.7527620e2 3.7527620e2 3.8169721e2   3.6948620e2 3.7004750e2 3.8546011e2 3.8178087e2 3.6114679e2
1 2.25  2.5 3.4521553e2 3.4521553e2 3.5114200e2   3.3973093e2 3.4031501e2 3.5487551e2 3.5131193e2 3.3214336e2
2  2.5 2.75 3.0001406e2 3.0001406e2 3.0511561e2   2.9517901e2 2.9567460e2 3.0860377e2 3.0541248e2 2.8866119e2
3 2.75    3 2.4257663e2 2.4257663e2 2.4665730e2   2.3861256e2 2.3902145e2 2.4965490e2 2.4699938e2 2.3342681e2
4    3 3.25 1.8093343e2 1.8093343e2 1.8387616e2   1.7800964e2 1.7821415e2 1.8627534e2 1.8431630e2 1.7416314e2
5 3.25  3.5 1.2291115e2 1.2291115e2 1.2481060e2   1.2097578e2 1.2099928e2 1.2657016e2 1.2528958e2 1.1835555e2
6  3.5    4 5.7851018e1 5.7851018e1 5.8647563e1   5.7008512e1 5.6897537e1 5.9567473e1 5.9037178e1 5.5752518e1
7    4  4.5 1.3772029e1 1.3772029e1 1.3903642e1   1.3622752e1 1.3512675e1 1.4165115e1 1.4094674e1 1.3296051e1
";

const SCALE_ENV_STR: &str = "b   etal    disg/detal  7pt-svar (env) 
     []        [pb]           [%]      
-+----+----+-----------+-------+-------
0    2 2.25 3.7527620e2   -3.77    2.71
1 2.25  2.5 3.4521553e2   -3.79    2.80
2  2.5 2.75 3.0001406e2   -3.78    2.86
3 2.75    3 2.4257663e2   -3.77    2.92
4    3 3.25 1.8093343e2   -3.74    2.95
5 3.25  3.5 1.2291115e2   -3.71    2.98
6  3.5    4 5.7851018e1   -3.63    2.97
7    4  4.5 1.3772029e1   -3.46    2.85
";

const SCALE_ENV_9_STR: &str = "b   etal    disg/detal  9pt-svar (env) 
     []        [pb]           [%]      
-+----+----+-----------+-------+-------
0    2 2.25 3.7527620e2   -5.55    3.96
1 2.25  2.5 3.4521553e2   -5.55    4.14
2  2.5 2.75 3.0001406e2   -5.53    4.31
3 2.75    3 2.4257663e2   -5.49    4.46
4    3 3.25 1.8093343e2   -5.45    4.60
5 3.25  3.5 1.2291115e2   -5.42    4.76
6  3.5    4 5.7851018e1   -5.37    4.95
7    4  4.5 1.3772029e1   -5.36    5.22
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
            "--silence-lhapdf",
            "uncert",
            "--pdf",
            "--threads=1",
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "--silence-lhapdf",
            "uncert",
            "--pdf",
            "--cl=90",
            "--threads=1",
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "--silence-lhapdf",
            "uncert",
            "--pdf",
            "--integrated",
            "--threads=1",
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "--silence-lhapdf",
            "uncert",
            "--pdf",
            "--orders=a2,as1a2",
            "--threads=1",
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "--silence-lhapdf",
            "uncert",
            "--scale-abs",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_ABS_STR);
}

#[test]
fn scale_env() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--silence-lhapdf",
            "uncert",
            "--scale-env",
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "--silence-lhapdf",
            "uncert",
            "--scale-env=9",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE_ENV_9_STR);
}
