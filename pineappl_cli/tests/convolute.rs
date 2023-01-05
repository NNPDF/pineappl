use assert_cmd::Command;

const HELP_STR: &str = "Convolutes a PineAPPL grid with a PDF set

Usage: pineappl convolute [OPTIONS] <INPUT> <PDFSETS>...

Arguments:
  <INPUT>       Path of the input grid
  <PDFSETS>...  LHAPDF id(s) or name of the PDF set(s)

Options:
  -a, --absolute            Show absolute numbers of the scale variation
  -b, --bins <BINS>...      Selects a subset of bins
  -i, --integrated          Show integrated numbers (without bin widths) instead of differential ones
  -o, --orders <ORDERS>...  Select orders manually
  -s, --scales <SCALES>     Set the number of scale variations [default: 7] [possible values: 1, 3, 7, 9]
      --digits-abs <ABS>    Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>    Set the number of fractional digits shown for relative numbers [default: 2]
  -h, --help                Print help information
";

const DEFAULT_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 3.7527620e2    -3.77     2.71
1 2.25  2.5 3.4521553e2    -3.79     2.80
2  2.5 2.75 3.0001406e2    -3.78     2.86
3 2.75    3 2.4257663e2    -3.77     2.92
4    3 3.25 1.8093343e2    -3.74     2.95
5 3.25  3.5 1.2291115e2    -3.71     2.98
6  3.5    4 5.7851018e1    -3.63     2.97
7    4  4.5 1.3772029e1    -3.46     2.85
";

const FORCE_POSITIVE_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 3.7528881e2    -3.76     2.71
1 2.25  2.5 3.4523428e2    -3.79     2.80
2  2.5 2.75 3.0004398e2    -3.78     2.86
3 2.75    3 2.4262216e2    -3.77     2.91
4    3 3.25 1.8101238e2    -3.74     2.94
5 3.25  3.5 1.2304887e2    -3.69     2.96
6  3.5    4 5.8141170e1    -3.60     2.92
7    4  4.5 1.4224044e1    -3.39     2.74
";

const DEFAULT_MULTIPLE_PDFS_STR: &str =
    "b   etal    disg/detal  scale uncertainty NNPDF31_nlo_as_0118_luxqed 
     []        [pb]            [%]                 [pb] [%]          
-+----+----+-----------+--------+--------+-------------+-------------
0    2 2.25 3.7527620e2    -3.77     2.71   3.7527620e2          0.00
1 2.25  2.5 3.4521553e2    -3.79     2.80   3.4521553e2          0.00
2  2.5 2.75 3.0001406e2    -3.78     2.86   3.0001406e2          0.00
3 2.75    3 2.4257663e2    -3.77     2.92   2.4257663e2          0.00
4    3 3.25 1.8093343e2    -3.74     2.95   1.8093343e2          0.00
5 3.25  3.5 1.2291115e2    -3.71     2.98   1.2291115e2          0.00
6  3.5    4 5.7851018e1    -3.63     2.97   5.7851018e1          0.00
7    4  4.5 1.3772029e1    -3.46     2.85   1.3772029e1          0.00
";

const MULTIPLE_PDFS_WITH_NEW_CONSTRUCTION_STR: &str =
    "b   etal    disg/detal  scale uncertainty NNPDF31_nlo_as_0118_luxqed/1 
     []        [pb]            [%]                  [pb] [%]           
-+----+----+-----------+--------+--------+--------------+--------------
0    2 2.25 3.7527620e2    -3.77     2.71    3.7379477e2          -0.39
1 2.25  2.5 3.4521553e2    -3.79     2.80    3.4316002e2          -0.60
2  2.5 2.75 3.0001406e2    -3.78     2.86    2.9780437e2          -0.74
3 2.75    3 2.4257663e2    -3.77     2.92    2.4059099e2          -0.82
4    3 3.25 1.8093343e2    -3.74     2.95    1.7941535e2          -0.84
5 3.25  3.5 1.2291115e2    -3.71     2.98    1.2195463e2          -0.78
6  3.5    4 5.7851018e1    -3.63     2.97    5.7551676e1          -0.52
7    4  4.5 1.3772029e1    -3.46     2.85    1.3640796e1          -0.95
";

const MULTIPLE_PDFS_WITH_RELABELING_STR: &str =
    "b   etal    disg/detal  scale uncertainty       other      
     []        [pb]            [%]            [pb] [%]     
-+----+----+-----------+--------+--------+-----------+-----
0    2 2.25 3.7527620e2    -3.77     2.71 3.7379477e2 -0.39
1 2.25  2.5 3.4521553e2    -3.79     2.80 3.4316002e2 -0.60
2  2.5 2.75 3.0001406e2    -3.78     2.86 2.9780437e2 -0.74
3 2.75    3 2.4257663e2    -3.77     2.92 2.4059099e2 -0.82
4    3 3.25 1.8093343e2    -3.74     2.95 1.7941535e2 -0.84
5 3.25  3.5 1.2291115e2    -3.71     2.98 1.2195463e2 -0.78
6  3.5    4 5.7851018e1    -3.63     2.97 5.7551676e1 -0.52
7    4  4.5 1.3772029e1    -3.46     2.85 1.3640796e1 -0.95
";

const WRONG_LHAID_STR: &str =
    "error: Invalid value '0' for '<PDFSETS>...': The PDF set for the LHAPDF ID `0` was not found

For more information try '--help'
";

const WRONG_PDFSET_STR: &str =
    "error: Invalid value 'IDONTEXIST' for '<PDFSETS>...': The PDF set `IDONTEXIST` was not found

For more information try '--help'
";

const ABSOLUTE_STR: &str =
"b   etal    disg/detal     (1,1)       (2,2)     (0.5,0.5)     (2,1)       (1,2)      (0.5,1)     (1,0.5)  
     []        [pb]        [pb]        [pb]        [pb]        [pb]        [pb]        [pb]        [pb]    
-+----+----+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------
0    2 2.25 3.7527620e2 3.7527620e2 3.8169721e2 3.6948620e2 3.7004750e2 3.8546011e2 3.8178087e2 3.6114679e2
1 2.25  2.5 3.4521553e2 3.4521553e2 3.5114200e2 3.3973093e2 3.4031501e2 3.5487551e2 3.5131193e2 3.3214336e2
2  2.5 2.75 3.0001406e2 3.0001406e2 3.0511561e2 2.9517901e2 2.9567460e2 3.0860377e2 3.0541248e2 2.8866119e2
3 2.75    3 2.4257663e2 2.4257663e2 2.4665730e2 2.3861256e2 2.3902145e2 2.4965490e2 2.4699938e2 2.3342681e2
4    3 3.25 1.8093343e2 1.8093343e2 1.8387616e2 1.7800964e2 1.7821415e2 1.8627534e2 1.8431630e2 1.7416314e2
5 3.25  3.5 1.2291115e2 1.2291115e2 1.2481060e2 1.2097578e2 1.2099928e2 1.2657016e2 1.2528958e2 1.1835555e2
6  3.5    4 5.7851018e1 5.7851018e1 5.8647563e1 5.7008512e1 5.6897537e1 5.9567473e1 5.9037178e1 5.5752518e1
7    4  4.5 1.3772029e1 1.3772029e1 1.3903642e1 1.3622752e1 1.3512675e1 1.4165115e1 1.4094674e1 1.3296051e1
";

const ABSOLUTE_INTEGRATED_STR: &str =
"b   etal       integ       (1,1)       (2,2)     (0.5,0.5)     (2,1)       (1,2)      (0.5,1)     (1,0.5)  
     []         []          []          []          []          []          []          []          []     
-+----+----+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------
0    2 2.25 9.3819050e1 9.3819050e1 9.5424302e1 9.2371549e1 9.2511874e1 9.6365028e1 9.5445218e1 9.0286698e1
1 2.25  2.5 8.6303882e1 8.6303882e1 8.7785500e1 8.4932734e1 8.5078752e1 8.8718877e1 8.7827983e1 8.3035839e1
2  2.5 2.75 7.5003515e1 7.5003515e1 7.6278903e1 7.3794753e1 7.3918649e1 7.7150944e1 7.6353121e1 7.2165298e1
3 2.75    3 6.0644157e1 6.0644157e1 6.1664326e1 5.9653141e1 5.9755362e1 6.2413725e1 6.1749845e1 5.8356702e1
4    3 3.25 4.5233358e1 4.5233358e1 4.5969040e1 4.4502410e1 4.4553538e1 4.6568835e1 4.6079075e1 4.3540785e1
5 3.25  3.5 3.0727788e1 3.0727788e1 3.1202651e1 3.0243945e1 3.0249820e1 3.1642540e1 3.1322395e1 2.9588888e1
6  3.5    4 2.8925509e1 2.8925509e1 2.9323782e1 2.8504256e1 2.8448768e1 2.9783737e1 2.9518589e1 2.7876259e1
7    4  4.5 6.8860146e0 6.8860146e0 6.9518210e0 6.8113762e0 6.7563375e0 7.0825577e0 7.0473370e0 6.6480254e0
";

const BINS_13567_STR: &str = "b   etal   disg/detal  scale uncertainty
     []       [pb]            [%]       
-+----+---+-----------+--------+--------
1 2.25 2.5 3.4521553e2    -3.79     2.80
3 2.75   3 2.4257663e2    -3.77     2.92
5 3.25 3.5 1.2291115e2    -3.71     2.98
6  3.5   4 5.7851018e1    -3.63     2.97
7    4 4.5 1.3772029e1    -3.46     2.85
";

const INTEGRATED_STR: &str = "b   etal       integ    scale uncertainty
     []         []             [%]       
-+----+----+-----------+--------+--------
0    2 2.25 9.3819050e1    -3.77     2.71
1 2.25  2.5 8.6303882e1    -3.79     2.80
2  2.5 2.75 7.5003515e1    -3.78     2.86
3 2.75    3 6.0644157e1    -3.77     2.92
4    3 3.25 4.5233358e1    -3.74     2.95
5 3.25  3.5 3.0727788e1    -3.71     2.98
6  3.5    4 2.8925509e1    -3.63     2.97
7    4  4.5 6.8860146e0    -3.46     2.85
";

const INTEGRATED_MULTIPLE_PDFS_STR: &str =
    "b   etal       integ    scale uncertainty NNPDF31_nlo_as_0118_luxqed 
     []         []             [%]                  [] [%]           
-+----+----+-----------+--------+--------+-------------+-------------
0    2 2.25 9.3819050e1    -3.77     2.71   9.3819050e1          0.00
1 2.25  2.5 8.6303882e1    -3.79     2.80   8.6303882e1          0.00
2  2.5 2.75 7.5003515e1    -3.78     2.86   7.5003515e1          0.00
3 2.75    3 6.0644157e1    -3.77     2.92   6.0644157e1          0.00
4    3 3.25 4.5233358e1    -3.74     2.95   4.5233358e1          0.00
5 3.25  3.5 3.0727788e1    -3.71     2.98   3.0727788e1          0.00
6  3.5    4 2.8925509e1    -3.63     2.97   2.8925509e1          0.00
7    4  4.5 6.8860146e0    -3.46     2.85   6.8860146e0          0.00
";

const ORDERS_A2_A3_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 3.2092052e2    -9.18     7.92
1 2.25  2.5 2.9427151e2    -8.68     7.41
2  2.5 2.75 2.5490261e2    -8.12     6.84
3 2.75    3 2.0561831e2    -7.55     6.26
4    3 3.25 1.5266481e2    -6.97     5.68
5 3.25  3.5 1.0303603e2    -6.38     5.09
6  3.5    4 4.7938981e1    -5.59     4.31
7    4  4.5 1.1075878e1    -4.60     3.35
";

const WRONG_ORDERS_STR: &str = "error: Invalid value 'a2a2as2' for '--orders <ORDERS>...': unable to parse order; too many couplings in 'a2a2as2'

For more information try '--help'
";

const SCALES_9_STR: &str = "b   etal    disg/detal  scale uncertainty
     []        [pb]            [%]       
-+----+----+-----------+--------+--------
0    2 2.25 3.7527620e2    -5.55     3.96
1 2.25  2.5 3.4521553e2    -5.55     4.14
2  2.5 2.75 3.0001406e2    -5.53     4.31
3 2.75    3 2.4257663e2    -5.49     4.46
4    3 3.25 1.8093343e2    -5.45     4.60
5 3.25  3.5 1.2291115e2    -5.42     4.76
6  3.5    4 5.7851018e1    -5.37     4.95
7    4  4.5 1.3772029e1    -5.36     5.22
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["convolute", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn default() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}

#[test]
fn force_positive() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--force-positive",
            "--silence-lhapdf",
            "convolute",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(FORCE_POSITIVE_STR);
}

#[test]
fn default_multiple_pdfs() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
            "324900=NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_MULTIPLE_PDFS_STR);
}

#[test]
fn multiple_pdfs_with_new_construction() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed/0",
            "NNPDF31_nlo_as_0118_luxqed/1",
        ])
        .assert()
        .success()
        .stdout(MULTIPLE_PDFS_WITH_NEW_CONSTRUCTION_STR);
}

#[test]
fn multiple_pdfs_with_relabeling() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
            "NNPDF31_nlo_as_0118_luxqed/1=other",
        ])
        .assert()
        .success()
        .stdout(MULTIPLE_PDFS_WITH_RELABELING_STR);
}

#[test]
fn wrong_lhaid() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "0",
        ])
        .assert()
        .failure()
        .stderr(WRONG_LHAID_STR);
}

#[test]
fn wrong_pdfset() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "IDONTEXIST",
        ])
        .assert()
        .failure()
        .stderr(WRONG_PDFSET_STR);
}

#[test]
fn absolute() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "--absolute",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ABSOLUTE_STR);
}

#[test]
fn absolute_integrated() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "--absolute",
            "--integrated",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ABSOLUTE_INTEGRATED_STR);
}

#[test]
fn bins_13567() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "--bins=1,3,5-7",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(BINS_13567_STR);
}

#[test]
fn integrated() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "--integrated",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(INTEGRATED_STR);
}

#[test]
fn integrated_multiple_pdfs() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "--integrated",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(INTEGRATED_MULTIPLE_PDFS_STR);
}

#[test]
fn orders_a2_a3() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "--orders=a2,a3",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ORDERS_A2_A3_STR);
}

#[test]
fn wrong_orders() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "--orders=a2a2as2",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(WRONG_ORDERS_STR);
}

#[test]
fn scales_9() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "convolute",
            "--scales=9",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALES_9_STR);
}
