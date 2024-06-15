use assert_cmd::Command;
use assert_fs::NamedTempFile;

const HELP_STR: &str = "Compares the numerical content of two grids with each other

Usage: pineappl diff [OPTIONS] <INPUT1> <INPUT2> <CONV_FUNS>...

Arguments:
  <INPUT1>        Path to the first grid
  <INPUT2>        Path to the second grid
  <CONV_FUNS>...  LHAPDF ID(s) or name(s) of the PDF(s)/FF(s)

Options:
      --ignore-orders      Ignore differences in the orders and sum them
      --ignore-bin-limits  Ignore bin limits (but not number of bins)
      --ignore-channels    Ignore differences in the channel definition
      --orders1 <ORDERS1>  Select orders of the first grid
      --orders2 <ORDERS2>  Select orders of the second grid
      --scale1 <SCALE1>    Scale all results of the first grid [default: 1.0]
      --scale2 <SCALE2>    Scale all results of the second grid [default: 1.0]
      --digits-abs <ABS>   Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>   Set the number of fractional digits shown for relative numbers [default: 3]
  -h, --help               Print help
";

const ORDERS1_A2_ORDERS2_A2_STR: &str = "b    x1               O(as^0 a^2)          
-+----+----+-----------+-----------+-------
0    2 2.25 6.5070305e2 6.5070305e2 0.000e0
1 2.25  2.5 5.9601236e2 5.9601236e2 0.000e0
2  2.5 2.75 5.1561247e2 5.1561247e2 0.000e0
3 2.75    3 4.1534629e2 4.1534629e2 0.000e0
4    3 3.25 3.0812719e2 3.0812719e2 0.000e0
5 3.25  3.5 2.0807482e2 2.0807482e2 0.000e0
6  3.5    4 9.6856769e1 9.6856769e1 0.000e0
7    4  4.5 2.2383492e1 2.2383492e1 0.000e0
";

const ORDERS1_A2_A2AS1_ORDERS2_A2_A2AS1_STR: &str =
    "b    x1               O(as^0 a^2)                     O(as^1 a^2)          
-+----+----+-----------+-----------+-------+-----------+-----------+-------
0    2 2.25 6.5070305e2 6.5070305e2 0.000e0 1.1175729e2 1.1175729e2 0.000e0
1 2.25  2.5 5.9601236e2 5.9601236e2 0.000e0 1.0083341e2 1.0083341e2 0.000e0
2  2.5 2.75 5.1561247e2 5.1561247e2 0.000e0 8.9874343e1 8.9874343e1 0.000e0
3 2.75    3 4.1534629e2 4.1534629e2 0.000e0 7.3935106e1 7.3935106e1 0.000e0
4    3 3.25 3.0812719e2 3.0812719e2 0.000e0 5.6414554e1 5.6414554e1 0.000e0
5 3.25  3.5 2.0807482e2 2.0807482e2 0.000e0 3.9468336e1 3.9468336e1 0.000e0
6  3.5    4 9.6856769e1 9.6856769e1 0.000e0 1.9822014e1 1.9822014e1 0.000e0
7    4  4.5 2.2383492e1 2.2383492e1 0.000e0 5.3540011e0 5.3540011e0 0.000e0
";

const ORDERS1_A2_A2AS1_IGNORE_ORDERS_STR: &str = "b    x1                   diff               
-+----+----+-----------+-----------+---------
0    2 2.25 7.6246034e2 7.5459110e2 -1.032e-2
1 2.25  2.5 6.9684577e2 6.9028342e2 -9.417e-3
2  2.5 2.75 6.0548681e2 6.0025198e2 -8.646e-3
3 2.75    3 4.8928139e2 4.8552235e2 -7.683e-3
4    3 3.25 3.6454175e2 3.6195456e2 -7.097e-3
5 3.25  3.5 2.4754316e2 2.4586691e2 -6.772e-3
6  3.5    4 1.1667878e2 1.1586851e2 -6.944e-3
7    4  4.5 2.7737493e1 2.7517266e1 -7.940e-3
";

const SCALE2_2_STR: &str = "b    x1               O(as^0 a^2)                       O(as^0 a^3)                       O(as^1 a^2)          
-+----+----+-----------+-----------+-------+-------------+-------------+-------+-----------+-----------+-------
0    2 2.25 6.5070305e2 1.3014061e3 1.000e0  -7.8692484e0  -1.5738497e1 1.000e0 1.1175729e2 2.2351458e2 1.000e0
1 2.25  2.5 5.9601236e2 1.1920247e3 1.000e0  -6.5623495e0  -1.3124699e1 1.000e0 1.0083341e2 2.0166682e2 1.000e0
2  2.5 2.75 5.1561247e2 1.0312249e3 1.000e0  -5.2348261e0  -1.0469652e1 1.000e0 8.9874343e1 1.7974869e2 1.000e0
3 2.75    3 4.1534629e2 8.3069258e2 1.000e0  -3.7590420e0  -7.5180840e0 1.000e0 7.3935106e1 1.4787021e2 1.000e0
4    3 3.25 3.0812719e2 6.1625439e2 1.000e0  -2.5871885e0  -5.1743770e0 1.000e0 5.6414554e1 1.1282911e2 1.000e0
5 3.25  3.5 2.0807482e2 4.1614964e2 1.000e0  -1.6762487e0  -3.3524974e0 1.000e0 3.9468336e1 7.8936673e1 1.000e0
6  3.5    4 9.6856769e1 1.9371354e2 1.000e0 -8.1027456e-1  -1.6205491e0 1.000e0 1.9822014e1 3.9644028e1 1.000e0
7    4  4.5 2.2383492e1 4.4766985e1 1.000e0 -2.2022770e-1 -4.4045540e-1 1.000e0 5.3540011e0 1.0708002e1 1.000e0
";

const ORDERS_DIFFER_STR: &str = "Error: selected orders differ
";

const BIN_LIMITS_DIFFER_STR: &str = "Error: bins limits differ
";

const BIN_NUMBER_DIFFERS_STR: &str = "Error: number of bins differ
";

const CHANNELS_DIFFER_STR: &str = "Error: channels differ
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["diff", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn orders1_a2_orders2_a2() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "diff",
            "--orders1=a2",
            "--orders2=a2",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .args([
            "diff",
            "--orders1=a2,a2as1",
            "--orders2=a2,a2as1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .args([
            "diff",
            "--orders1=a2,a2as1",
            "--ignore-orders",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .args([
            "diff",
            "--scale2=2",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .args([
            "diff",
            "--orders1=a2",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .args([
            "write",
            "--remap=0,1,2,3,4,5,6,7,8",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "diff",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
        .args([
            "write",
            "--delete-bins=0,1",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "diff",
            "--ignore-bin-limits",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(BIN_NUMBER_DIFFERS_STR)
        .stdout("");
}

#[test]
fn channels_differ() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "diff",
            "../test-data/LHCB_WP_7TEV_old.pineappl.lz4",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(CHANNELS_DIFFER_STR)
        .stdout("");
}
