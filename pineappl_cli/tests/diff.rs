use assert_cmd::Command;
use assert_fs::NamedTempFile;

const HELP_STR: &str = "Compares the numerical content of two grids with each other

Usage: pineappl diff [OPTIONS] <INPUT1> <INPUT2> <PDFSET>

Arguments:
  <INPUT1>  Path to the first grid
  <INPUT2>  Path to the second grid
  <PDFSET>  LHAPDF id or name of the PDF set

Options:
      --ignore-orders         Ignore differences in the orders and sum them
      --ignore-bin-limits     Ignore bin limits (but not number of bins)
      --ignore-lumis          Ignore differences in the luminosity functions
      --orders1 <ORDERS1>...  Select orders of the first grid
      --orders2 <ORDERS2>...  Select orders of the second grid
      --scale1 <SCALE1>       Scale all results of the first grid [default: 1.0]
      --scale2 <SCALE2>       Scale all results of the second grid [default: 1.0]
      --digits-abs <ABS>      Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>      Set the number of fractional digits shown for relative numbers [default: 3]
      --force-positive        Forces negative PDF values to zero
  -h, --help                  Print help information
";

const ORDERS1_A2_ORDERS2_A2_STR: &str = "b    x1               O(as^0 a^2)          
-+----+----+-----------+-----------+-------
0    2 2.25 3.2482657e2 3.2482657e2 0.000e0
1 2.25  2.5 2.9755128e2 2.9755128e2 0.000e0
2  2.5 2.75 2.5751142e2 2.5751142e2 0.000e0
3 2.75    3 2.0748091e2 2.0748091e2 0.000e0
4    3 3.25 1.5397599e2 1.5397599e2 0.000e0
5 3.25  3.5 1.0384063e2 1.0384063e2 0.000e0
6  3.5    4 4.8383606e1 4.8383606e1 0.000e0
7    4  4.5 1.1185365e1 1.1185365e1 0.000e0
";

const ORDERS1_A2_A2AS1_ORDERS2_A2_A2AS1_STR: &str =
    "b    x1               O(as^0 a^2)                     O(as^1 a^2)          
-+----+----+-----------+-----------+-------+-----------+-----------+-------
0    2 2.25 3.2482657e2 3.2482657e2 0.000e0 5.4355679e1 5.4355679e1 0.000e0
1 2.25  2.5 2.9755128e2 2.9755128e2 0.000e0 5.0944018e1 5.0944018e1 0.000e0
2  2.5 2.75 2.5751142e2 2.5751142e2 0.000e0 4.5111446e1 4.5111446e1 0.000e0
3 2.75    3 2.0748091e2 2.0748091e2 0.000e0 3.6958317e1 3.6958317e1 0.000e0
4    3 3.25 1.5397599e2 1.5397599e2 0.000e0 2.8268620e1 2.8268620e1 0.000e0
5 3.25  3.5 1.0384063e2 1.0384063e2 0.000e0 1.9875123e1 1.9875123e1 0.000e0
6  3.5    4 4.8383606e1 4.8383606e1 0.000e0 9.9120372e0 9.9120372e0 0.000e0
7    4  4.5 1.1185365e1 1.1185365e1 0.000e0 2.6961509e0 2.6961509e0 0.000e0
";

const ORDERS1_A2_A2AS1_IGNORE_ORDERS_STR: &str = "b    x1                   diff              
-+----+----+-----------+-----------+--------
0    2 2.25 3.7918224e2 3.7527620e2 1.041e-2
1 2.25  2.5 3.4849530e2 3.4521553e2 9.501e-3
2  2.5 2.75 3.0262287e2 3.0001406e2 8.696e-3
3 2.75    3 2.4443923e2 2.4257663e2 7.678e-3
4    3 3.25 1.8224461e2 1.8093343e2 7.247e-3
5 3.25  3.5 1.2371575e2 1.2291115e2 6.546e-3
6  3.5    4 5.8295643e1 5.7851018e1 7.686e-3
7    4  4.5 1.3881516e1 1.3772029e1 7.950e-3
";

const SCALE2_2_STR: &str = "b    x1                O(as^0 a^2)                         O(as^0 a^3)                         O(as^1 a^2)           
-+----+----+-----------+-----------+---------+-------------+-------------+---------+-----------+-----------+---------
0    2 2.25 3.2482657e2 6.4965313e2 -5.000e-1  -3.9060418e0  -7.8120836e0 -5.000e-1 5.4355679e1 1.0871136e2 -5.000e-1
1 2.25  2.5 2.9755128e2 5.9510256e2 -5.000e-1  -3.2797697e0  -6.5595394e0 -5.000e-1 5.0944018e1 1.0188804e2 -5.000e-1
2  2.5 2.75 2.5751142e2 5.1502284e2 -5.000e-1  -2.6088069e0  -5.2176138e0 -5.000e-1 4.5111446e1 9.0222892e1 -5.000e-1
3 2.75    3 2.0748091e2 4.1496182e2 -5.000e-1  -1.8626008e0  -3.7252015e0 -5.000e-1 3.6958317e1 7.3916633e1 -5.000e-1
4    3 3.25 1.5397599e2 3.0795198e2 -5.000e-1  -1.3111794e0  -2.6223588e0 -5.000e-1 2.8268620e1 5.6537240e1 -5.000e-1
5 3.25  3.5 1.0384063e2 2.0768125e2 -5.000e-1 -8.0459807e-1  -1.6091961e0 -5.000e-1 1.9875123e1 3.9750247e1 -5.000e-1
6  3.5    4 4.8383606e1 9.6767212e1 -5.000e-1 -4.4462513e-1 -8.8925027e-1 -5.000e-1 9.9120372e0 1.9824074e1 -5.000e-1
7    4  4.5 1.1185365e1 2.2370731e1 -5.000e-1 -1.0948700e-1 -2.1897400e-1 -5.000e-1 2.6961509e0 5.3923018e0 -5.000e-1
";

const ORDERS_DIFFER_STR: &str = "Error: selected orders differ
";

const BIN_LIMITS_DIFFER_STR: &str = "Error: bins limits differ
";

const BIN_NUMBER_DIFFERS_STR: &str = "Error: number of bins differ
";

const LUMIS_DIFFER_STR: &str = "Error: luminosities differ
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["diff", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn orders1_a2_orders2_a2() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "diff",
            "--orders1=a2",
            "--orders2=a2",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ORDERS1_A2_ORDERS2_A2_STR);
}

#[test]
fn orders1_a2_a2as1_orders2_a2_a2as1() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "diff",
            "--orders1=a2,a2as1",
            "--orders2=a2,a2as1",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ORDERS1_A2_A2AS1_ORDERS2_A2_A2AS1_STR);
}

#[test]
fn orders1_a2_a2as1_ignore_orders() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "diff",
            "--orders1=a2,a2as1",
            "--ignore-orders",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ORDERS1_A2_A2AS1_IGNORE_ORDERS_STR);
}

#[test]
fn scale2_2() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "diff",
            "--scale2=2",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(SCALE2_2_STR);
}

#[test]
fn orders_differ() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "diff",
            "--orders1=a2",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(ORDERS_DIFFER_STR)
        .stdout("");
}

#[test]
fn bin_limits_differ() {
    let output = NamedTempFile::new("remapped.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "remap",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
            "0,1,2,3,4,5,6,7,8",
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "diff",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(BIN_LIMITS_DIFFER_STR)
        .stdout("");
}

#[test]
fn bin_number_differs() {
    let output = NamedTempFile::new("remapped.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "delete",
            "--bins=0,1",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "diff",
            "--ignore-bin-limits",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(BIN_NUMBER_DIFFERS_STR)
        .stdout("");
}

#[test]
fn lumis_differ() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&[
            "--silence-lhapdf",
            "diff",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "data/LHCB_WP_7TEV_new.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(LUMIS_DIFFER_STR)
        .stdout("");
}
