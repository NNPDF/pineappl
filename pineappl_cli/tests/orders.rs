use assert_cmd::Command;

const HELP_STR: &str = "Shows the predictions for all bin for each order separately

Usage: pineappl orders [OPTIONS] <INPUT> <CONV_FUNS>

Arguments:
  <INPUT>      Path to the input grid
  <CONV_FUNS>  LHAPDF ID(s) or name(s) of the PDF(s)/FF(s)

Options:
  -a, --absolute               Show absolute numbers of each perturbative order
  -i, --integrated             Show integrated numbers (without bin widths) instead of differential ones
  -n, --normalize <NORMALIZE>  Normalize contributions to the specified orders
      --digits-abs <ABS>       Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>       Set the number of fractional digits shown for relative numbers [default: 2]
  -h, --help                   Print help
";

const DEFAULT_STR: &str = "b   etal    dsig/detal  O(as^0 a^2) O(as^1 a^2) O(as^0 a^3)
     []        [pb]         [%]         [%]         [%]    
-+----+----+-----------+-----------+-----------+-----------
0    2 2.25 7.5459110e2      100.00       17.17       -1.21
1 2.25  2.5 6.9028342e2      100.00       16.92       -1.10
2  2.5 2.75 6.0025198e2      100.00       17.43       -1.02
3 2.75    3 4.8552235e2      100.00       17.80       -0.91
4    3 3.25 3.6195456e2      100.00       18.31       -0.84
5 3.25  3.5 2.4586691e2      100.00       18.97       -0.81
6  3.5    4 1.1586851e2      100.00       20.47       -0.84
7    4  4.5 2.7517266e1      100.00       23.92       -0.98
";

const ABSOLUTE_STR: &str = "b   etal    dsig/detal  O(as^0 a^2) O(as^1 a^2)  O(as^0 a^3) 
     []        [pb]        [pb]        [pb]         [pb]     
-+----+----+-----------+-----------+-----------+-------------
0    2 2.25 7.5459110e2 6.5070305e2 1.1175729e2  -7.8692484e0
1 2.25  2.5 6.9028342e2 5.9601236e2 1.0083341e2  -6.5623495e0
2  2.5 2.75 6.0025198e2 5.1561247e2 8.9874343e1  -5.2348261e0
3 2.75    3 4.8552235e2 4.1534629e2 7.3935106e1  -3.7590420e0
4    3 3.25 3.6195456e2 3.0812719e2 5.6414554e1  -2.5871885e0
5 3.25  3.5 2.4586691e2 2.0807482e2 3.9468336e1  -1.6762487e0
6  3.5    4 1.1586851e2 9.6856769e1 1.9822014e1 -8.1027456e-1
7    4  4.5 2.7517266e1 2.2383492e1 5.3540011e0 -2.2022770e-1
";

const ABSOLUTE_INTEGRATED_STR: &str =
    "b   etal       integ    O(as^0 a^2) O(as^1 a^2)  O(as^0 a^3) 
     []         []          []          []           []      
-+----+----+-----------+-----------+-----------+-------------
0    2 2.25 1.8864777e2 1.6267576e2 2.7939322e1  -1.9673121e0
1 2.25  2.5 1.7257086e2 1.4900309e2 2.5208352e1  -1.6405874e0
2  2.5 2.75 1.5006300e2 1.2890312e2 2.2468586e1  -1.3087065e0
3 2.75    3 1.2138059e2 1.0383657e2 1.8483776e1 -9.3976050e-1
4    3 3.25 9.0488640e1 7.7031799e1 1.4103638e1 -6.4679713e-1
5 3.25  3.5 6.1466727e1 5.2018706e1 9.8670841e0 -4.1906217e-1
6  3.5    4 5.7934254e1 4.8428384e1 9.9110071e0 -4.0513728e-1
7    4  4.5 1.3758633e1 1.1191746e1 2.6770006e0 -1.1011385e-1
";

const INTEGRATED_STR: &str = "b   etal       integ    O(as^0 a^2) O(as^1 a^2) O(as^0 a^3)
     []         []          [%]         [%]         [%]    
-+----+----+-----------+-----------+-----------+-----------
0    2 2.25 1.8864777e2      100.00       17.17       -1.21
1 2.25  2.5 1.7257086e2      100.00       16.92       -1.10
2  2.5 2.75 1.5006300e2      100.00       17.43       -1.02
3 2.75    3 1.2138059e2      100.00       17.80       -0.91
4    3 3.25 9.0488640e1      100.00       18.31       -0.84
5 3.25  3.5 6.1466727e1      100.00       18.97       -0.81
6  3.5    4 5.7934254e1      100.00       20.47       -0.84
7    4  4.5 1.3758633e1      100.00       23.92       -0.98
";

const NORMALIZE_A2_AS1A2_STR: &str = "b   etal    dsig/detal  O(as^0 a^2) O(as^1 a^2) O(as^0 a^3)
     []        [pb]         [%]         [%]         [%]    
-+----+----+-----------+-----------+-----------+-----------
0    2 2.25 7.5459110e2       85.34       14.66       -1.03
1 2.25  2.5 6.9028342e2       85.53       14.47       -0.94
2  2.5 2.75 6.0025198e2       85.16       14.84       -0.86
3 2.75    3 4.8552235e2       84.89       15.11       -0.77
4    3 3.25 3.6195456e2       84.52       15.48       -0.71
5 3.25  3.5 2.4586691e2       84.06       15.94       -0.68
6  3.5    4 1.1586851e2       83.01       16.99       -0.69
7    4  4.5 2.7517266e1       80.70       19.30       -0.79
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["orders", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn default() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "orders",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}

#[test]
fn absolute() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "orders",
            "--absolute",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
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
        .args([
            "orders",
            "--absolute",
            "--integrated",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ABSOLUTE_INTEGRATED_STR);
}

#[test]
fn integrated() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "orders",
            "--integrated",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(INTEGRATED_STR);
}

#[test]
fn normalize_a2_as1a2() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "orders",
            "--normalize=a2,as1a2",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(NORMALIZE_A2_AS1A2_STR);
}
