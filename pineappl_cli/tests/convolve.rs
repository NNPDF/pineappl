#![allow(missing_docs)]

use assert_cmd::Command;
use predicates::str;

const HELP_STR: &str = "Convolutes a PineAPPL grid with a PDF set

Usage: pineappl convolve [OPTIONS] <INPUT> <CONV_FUNS>...

Arguments:
  <INPUT>         Path of the input grid
  <CONV_FUNS>...  LHAPDF id(s) or name of the PDF set(s)

Options:
  -b, --bins <BINS>       Selects a subset of bins
  -i, --integrated        Show integrated numbers (without bin widths) instead of differential ones
  -o, --orders <ORDERS>   Select orders manually
      --xir <XIR>         Set the variation of the renormalization scale [default: 1.0]
      --xif <XIF>         Set the variation of the factorization scale [default: 1.0]
      --xia <XIA>         Set the variation of the fragmentation scale [default: 1.0]
      --digits-abs <ABS>  Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>  Set the number of fractional digits shown for relative numbers [default: 2]
  -h, --help              Print help
";

const DEFAULT_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5459110e2
1 2.25  2.5 6.9028342e2
2  2.5 2.75 6.0025198e2
3 2.75    3 4.8552235e2
4    3 3.25 3.6195456e2
5 3.25  3.5 2.4586691e2
6  3.5    4 1.1586851e2
7    4  4.5 2.7517266e1
";

const USE_ALPHAS_FROM_ERROR_STR: &str =
    "expected `use_alphas_from` to be an integer within `[0, 2)`, but got `2`
";

const FORCE_POSITIVE_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.5461571e2
1 2.25  2.5 6.9032107e2
2  2.5 2.75 6.0031056e2
3 2.75    3 4.8561541e2
4    3 3.25 3.6211174e2
5 3.25  3.5 2.4614249e2
6  3.5    4 1.1644853e2
7    4  4.5 2.8422453e1
";

const DEFAULT_MULTIPLE_PDFS_STR: &str = "b   etal    dsig/detal  NNPDF31_nlo_as_0118_luxqed 
     []        [pb]              [pb] [%]          
-+----+----+-----------+-------------+-------------
0    2 2.25 7.5459110e2   7.5459110e2          0.00
1 2.25  2.5 6.9028342e2   6.9028342e2          0.00
2  2.5 2.75 6.0025198e2   6.0025198e2          0.00
3 2.75    3 4.8552235e2   4.8552235e2          0.00
4    3 3.25 3.6195456e2   3.6195456e2          0.00
5 3.25  3.5 2.4586691e2   2.4586691e2          0.00
6  3.5    4 1.1586851e2   1.1586851e2          0.00
7    4  4.5 2.7517266e1   2.7517266e1          0.00
";

const MULTIPLE_PDFS_WITH_NEW_CONSTRUCTION_STR: &str =
    "b   etal    dsig/detal  NNPDF31_nlo_as_0118_luxqed/1 
     []        [pb]               [pb] [%]           
-+----+----+-----------+--------------+--------------
0    2 2.25 7.5459110e2    7.5169067e2          -0.38
1 2.25  2.5 6.9028342e2    6.8612437e2          -0.60
2  2.5 2.75 6.0025198e2    5.9582118e2          -0.74
3 2.75    3 4.8552235e2    4.8155744e2          -0.82
4    3 3.25 3.6195456e2    3.5891650e2          -0.84
5 3.25  3.5 2.4586691e2    2.4395886e2          -0.78
6  3.5    4 1.1586851e2    1.1526800e2          -0.52
7    4  4.5 2.7517266e1    2.7259743e1          -0.94
";

const MULTIPLE_PDFS_WITH_RELABELING_STR: &str = "b   etal    dsig/detal     other mc=1.4   
     []        [pb]          [pb] [%]     
-+----+----+-----------+-----------+------
0    2 2.25 7.5459110e2 7.5169067e2  -0.38
1 2.25  2.5 6.9028342e2 6.8612437e2  -0.60
2  2.5 2.75 6.0025198e2 5.9582118e2  -0.74
3 2.75    3 4.8552235e2 4.8155744e2  -0.82
4    3 3.25 3.6195456e2 3.5891650e2  -0.84
5 3.25  3.5 2.4586691e2 2.4395886e2  -0.78
6  3.5    4 1.1586851e2 1.1526800e2  -0.52
7    4  4.5 2.7517266e1 2.7259743e1  -0.94
";

const TWO_PDFS_WITH_ORDER_SUBSET_STR: &str = "b   etal    dsig/detal  NNPDF31_nlo_as_0118_luxqed/1 
     []        [pb]               [pb] [%]           
-+----+----+-----------+--------------+--------------
0    2 2.25 6.5070305e2    6.4903968e2          -0.26
1 2.25  2.5 5.9601236e2    5.9329838e2          -0.46
2  2.5 2.75 5.1561247e2    5.1256355e2          -0.59
3 2.75    3 4.1534629e2    4.1255233e2          -0.67
4    3 3.25 3.0812719e2    3.0597140e2          -0.70
5 3.25  3.5 2.0807482e2    2.0674262e2          -0.64
6  3.5    4 9.6856769e1    9.6469932e1          -0.40
7    4  4.5 2.2383492e1    2.2182749e1          -0.90
";

const THREE_PDFS_STR: &str =
    "b   etal    dsig/detal  NNPDF31_nlo_as_0118_luxqed/1  NNPDF31_nlo_as_0118_luxqed/2 
     []        [pb]               [pb] [%]                      [pb] [%]           
-+----+----+-----------+--------------+--------------+--------------+--------------
0    2 2.25 7.5459110e2    7.5169067e2          -0.38    7.5748758e2           0.38
1 2.25  2.5 6.9028342e2    6.8612437e2          -0.60    6.9177921e2           0.22
2  2.5 2.75 6.0025198e2    5.9582118e2          -0.74    6.0069567e2           0.07
3 2.75    3 4.8552235e2    4.8155744e2          -0.82    4.8555961e2           0.01
4    3 3.25 3.6195456e2    3.5891650e2          -0.84    3.6199676e2           0.01
5 3.25  3.5 2.4586691e2    2.4395886e2          -0.78    2.4622502e2           0.15
6  3.5    4 1.1586851e2    1.1526800e2          -0.52    1.1661992e2           0.65
7    4  4.5 2.7517266e1    2.7259743e1          -0.94    2.8446007e1           3.38
";

const BINS_13567_STR: &str = "b   etal   dsig/detal 
     []       [pb]    
-+----+---+-----------
1 2.25 2.5 6.9028342e2
3 2.75   3 4.8552235e2
5 3.25 3.5 2.4586691e2
6  3.5   4 1.1586851e2
7    4 4.5 2.7517266e1
";

const INTEGRATED_STR: &str = "b   etal       integ   
     []         []     
-+----+----+-----------
0    2 2.25 1.8864777e2
1 2.25  2.5 1.7257086e2
2  2.5 2.75 1.5006300e2
3 2.75    3 1.2138059e2
4    3 3.25 9.0488640e1
5 3.25  3.5 6.1466727e1
6  3.5    4 5.7934254e1
7    4  4.5 1.3758633e1
";

const INTEGRATED_MULTIPLE_PDFS_STR: &str = "b   etal       integ    NNPDF31_nlo_as_0118_luxqed 
     []         []                [] [%]           
-+----+----+-----------+-------------+-------------
0    2 2.25 1.8864777e2   1.8864777e2          0.00
1 2.25  2.5 1.7257086e2   1.7257086e2          0.00
2  2.5 2.75 1.5006300e2   1.5006300e2          0.00
3 2.75    3 1.2138059e2   1.2138059e2          0.00
4    3 3.25 9.0488640e1   9.0488640e1          0.00
5 3.25  3.5 6.1466727e1   6.1466727e1          0.00
6  3.5    4 5.7934254e1   5.7934254e1          0.00
7    4  4.5 1.3758633e1   1.3758633e1          0.00
";

const ORDERS_A2_A3_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 6.4283381e2
1 2.25  2.5 5.8945001e2
2  2.5 2.75 5.1037764e2
3 2.75    3 4.1158725e2
4    3 3.25 3.0554001e2
5 3.25  3.5 2.0639857e2
6  3.5    4 9.6046494e1
7    4  4.5 2.2163265e1
";

const WRONG_ORDERS_STR: &str = "error: invalid value 'a2a3as2' for '--orders <ORDERS>': exponent of 'a' has already been set to '2'

For more information, try '--help'.
";

const XIR_XIF_STR: &str = "b   etal    dsig/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 7.6241231e2
1 2.25  2.5 6.9755130e2
2  2.5 2.75 6.0636076e2
3 2.75    3 4.9019741e2
4    3 3.25 3.6518490e2
5 3.25  3.5 2.4783934e2
6  3.5    4 1.1656958e2
7    4  4.5 2.7565811e1
";

const THREE_CONVOLUTIONS_STR: &str = "b                  pT                  dsig/dpT (pol)
                 [GeV]                    [pb/GeV]   
-+-----------------+------------------+--------------
0 5.108395099639893  6.045444965362549    2.2605116e3
1 6.045444965362549  6.982494831085205    1.0361301e3
2 6.982494831085205  7.992245197296143    4.8947508e2
3 7.992245197296143  8.960753917694092    2.4023939e2
4 8.960753917694092  9.929026126861572    1.2464463e2
5 9.929026126861572 11.660773754119873    5.2680349e1
";

const NO_CHANNELS_GRID_STR: &str =
    "b    Q2                            x                             F2d    
   [GeV^2]                        []                             []     
--+---+---+------------------------+------------------------+-----------
 0   5   5                0.0000005                0.0000005 0.0000000e0
 1   5   5 0.0000019407667236782136 0.0000019407667236782136 0.0000000e0
 2   5   5  0.000007533150951473337  0.000007533150951473337 0.0000000e0
 3   5   5  0.000029240177382128657  0.000029240177382128657 0.0000000e0
 4   5   5   0.00011349672651536727   0.00011349672651536727 0.0000000e0
 5   5   5   0.00044054134013486355   0.00044054134013486355 0.0000000e0
 6   5   5    0.0017099759466766963    0.0017099759466766963 0.0000000e0
 7   5   5     0.006637328831200572     0.006637328831200572 0.0000000e0
 8   5   5      0.02576301385940815      0.02576301385940815 0.0000000e0
 9   5   5                      0.1                      0.1 0.0000000e0
10   5   5                     0.18                     0.18 0.0000000e0
11   5   5                     0.26                     0.26 0.0000000e0
12   5   5      0.33999999999999997      0.33999999999999997 0.0000000e0
13   5   5      0.42000000000000004      0.42000000000000004 0.0000000e0
14   5   5                      0.5                      0.5 0.0000000e0
15   5   5                     0.58                     0.58 0.0000000e0
16   5   5                     0.66                     0.66 0.0000000e0
17   5   5                     0.74                     0.74 0.0000000e0
18   5   5                     0.82                     0.82 0.0000000e0
19   5   5                      0.9                      0.9 0.0000000e0
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["convolve", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn default() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--use-alphas-from=1",
            "convolve",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed,NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}

#[test]
fn use_alphas_from_error() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--use-alphas-from=2",
            "convolve",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed,NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(str::contains(USE_ALPHAS_FROM_ERROR_STR));
}

#[test]
fn force_positive() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--force-positive",
            "convolve",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
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
        .args([
            "convolve",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
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
        .args([
            "convolve",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
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
        .args([
            "convolve",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
            "NNPDF31_nlo_as_0118_luxqed/1=other mc=1.4",
        ])
        .assert()
        .success()
        .stdout(MULTIPLE_PDFS_WITH_RELABELING_STR);
}

#[test]
fn two_pdfs_with_order_subset() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            "--orders=a2",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed/0",
            "NNPDF31_nlo_as_0118_luxqed/1",
        ])
        .assert()
        .success()
        .stdout(TWO_PDFS_WITH_ORDER_SUBSET_STR);
}

#[test]
fn three_pdfs() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed/0",
            "NNPDF31_nlo_as_0118_luxqed/1",
            "NNPDF31_nlo_as_0118_luxqed/2",
        ])
        .assert()
        .success()
        .stdout(THREE_PDFS_STR);
}

#[test]
fn bins_13567() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            "--bins=1,3,5-7",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
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
        .args([
            "convolve",
            "--integrated",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
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
        .args([
            "convolve",
            "--integrated",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
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
        .args([
            "convolve",
            "--orders=a2,a3",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
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
        .args([
            "convolve",
            "--orders=a2a3as2",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(WRONG_ORDERS_STR);
}

#[test]
fn xir_xif() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            "--xir=2.34",
            "--xif=1.79",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(XIR_XIF_STR);
}

#[test]
fn three_convolutions() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            "../test-data/SIHP-PP-POLARIZED-STAR-NLO.pineappl.lz4",
            "NNPDFpol11_100+p,MAPFF10NLOPIsum+f",
        ])
        .assert()
        .success()
        .stdout(THREE_CONVOLUTIONS_STR);
}

#[test]
fn no_channels_grid() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "convolve",
            "../test-data/NNPDF_POS_F2D_40.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(NO_CHANNELS_GRID_STR);
}
