use assert_cmd::Command;

const HELP_STR: &str = "Convolutes a PineAPPL grid with a PDF set

Usage: pineappl convolute [OPTIONS] <INPUT> <PDFSETS>...

Arguments:
  <INPUT>       Path of the input grid
  <PDFSETS>...  LHAPDF id(s) or name of the PDF set(s)

Options:
  -b, --bins <BINS>       Selects a subset of bins
  -i, --integrated        Show integrated numbers (without bin widths) instead of differential ones
  -o, --orders <ORDERS>   Select orders manually
      --digits-abs <ABS>  Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>  Set the number of fractional digits shown for relative numbers [default: 2]
  -h, --help              Print help
";

const DEFAULT_STR: &str = "b   etal    disg/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 3.7527620e2
1 2.25  2.5 3.4521553e2
2  2.5 2.75 3.0001406e2
3 2.75    3 2.4257663e2
4    3 3.25 1.8093343e2
5 3.25  3.5 1.2291115e2
6  3.5    4 5.7851018e1
7    4  4.5 1.3772029e1
";

const FORCE_POSITIVE_STR: &str = "b   etal    disg/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 3.7528881e2
1 2.25  2.5 3.4523428e2
2  2.5 2.75 3.0004398e2
3 2.75    3 2.4262216e2
4    3 3.25 1.8101238e2
5 3.25  3.5 1.2304887e2
6  3.5    4 5.8141170e1
7    4  4.5 1.4224044e1
";

const DEFAULT_MULTIPLE_PDFS_STR: &str = "b   etal    disg/detal  NNPDF31_nlo_as_0118_luxqed 
     []        [pb]              [pb] [%]          
-+----+----+-----------+-------------+-------------
0    2 2.25 3.7527620e2   3.7527620e2          0.00
1 2.25  2.5 3.4521553e2   3.4521553e2          0.00
2  2.5 2.75 3.0001406e2   3.0001406e2          0.00
3 2.75    3 2.4257663e2   2.4257663e2          0.00
4    3 3.25 1.8093343e2   1.8093343e2          0.00
5 3.25  3.5 1.2291115e2   1.2291115e2          0.00
6  3.5    4 5.7851018e1   5.7851018e1          0.00
7    4  4.5 1.3772029e1   1.3772029e1          0.00
";

const MULTIPLE_PDFS_WITH_NEW_CONSTRUCTION_STR: &str =
    "b   etal    disg/detal  NNPDF31_nlo_as_0118_luxqed/1 
     []        [pb]               [pb] [%]           
-+----+----+-----------+--------------+--------------
0    2 2.25 3.7527620e2    3.7379477e2          -0.39
1 2.25  2.5 3.4521553e2    3.4316002e2          -0.60
2  2.5 2.75 3.0001406e2    2.9780437e2          -0.74
3 2.75    3 2.4257663e2    2.4059099e2          -0.82
4    3 3.25 1.8093343e2    1.7941535e2          -0.84
5 3.25  3.5 1.2291115e2    1.2195463e2          -0.78
6  3.5    4 5.7851018e1    5.7551676e1          -0.52
7    4  4.5 1.3772029e1    1.3640796e1          -0.95
";

const MULTIPLE_PDFS_WITH_RELABELING_STR: &str = "b   etal    disg/detal     other mc=1.4   
     []        [pb]          [pb] [%]     
-+----+----+-----------+-----------+------
0    2 2.25 3.7527620e2 3.7379477e2  -0.39
1 2.25  2.5 3.4521553e2 3.4316002e2  -0.60
2  2.5 2.75 3.0001406e2 2.9780437e2  -0.74
3 2.75    3 2.4257663e2 2.4059099e2  -0.82
4    3 3.25 1.8093343e2 1.7941535e2  -0.84
5 3.25  3.5 1.2291115e2 1.2195463e2  -0.78
6  3.5    4 5.7851018e1 5.7551676e1  -0.52
7    4  4.5 1.3772029e1 1.3640796e1  -0.95
";

const TWO_PDFS_WITH_ORDER_SUBSET_STR: &str = "b   etal    disg/detal  NNPDF31_nlo_as_0118_luxqed/1 
     []        [pb]               [pb] [%]           
-+----+----+-----------+--------------+--------------
0    2 2.25 3.2482657e2    3.2399084e2          -0.26
1 2.25  2.5 2.9755128e2    2.9619757e2          -0.45
2  2.5 2.75 2.5751142e2    2.5598573e2          -0.59
3 2.75    3 2.0748091e2    2.0608487e2          -0.67
4    3 3.25 1.5397599e2    1.5289865e2          -0.70
5 3.25  3.5 1.0384063e2    1.0317243e2          -0.64
6  3.5    4 4.8383606e1    4.8189863e1          -0.40
7    4  4.5 1.1185365e1    1.1083933e1          -0.91
";

const THREE_PDFS_STR: &str =
    "b   etal    disg/detal  NNPDF31_nlo_as_0118_luxqed/1  NNPDF31_nlo_as_0118_luxqed/2 
     []        [pb]               [pb] [%]                      [pb] [%]           
-+----+----+-----------+--------------+--------------+--------------+--------------
0    2 2.25 3.7527620e2    3.7379477e2          -0.39    3.7670804e2           0.38
1 2.25  2.5 3.4521553e2    3.4316002e2          -0.60    3.4596819e2           0.22
2  2.5 2.75 3.0001406e2    2.9780437e2          -0.74    3.0023944e2           0.08
3 2.75    3 2.4257663e2    2.4059099e2          -0.82    2.4259549e2           0.01
4    3 3.25 1.8093343e2    1.7941535e2          -0.84    1.8095662e2           0.01
5 3.25  3.5 1.2291115e2    1.2195463e2          -0.78    1.2308878e2           0.14
6  3.5    4 5.7851018e1    5.7551676e1          -0.52    5.8229565e1           0.65
7    4  4.5 1.3772029e1    1.3640796e1          -0.95    1.4234097e1           3.36
";

const WRONG_LHAID_STR: &str =
    "error: invalid value '0' for '<PDFSETS>...': The PDF set for the LHAPDF ID `0` was not found

For more information, try '--help'.
";

const WRONG_PDFSET_STR: &str =
    "error: invalid value 'IDONTEXIST' for '<PDFSETS>...': The PDF set `IDONTEXIST` was not found

For more information, try '--help'.
";

const BINS_13567_STR: &str = "b   etal   disg/detal 
     []       [pb]    
-+----+---+-----------
1 2.25 2.5 3.4521553e2
3 2.75   3 2.4257663e2
5 3.25 3.5 1.2291115e2
6  3.5   4 5.7851018e1
7    4 4.5 1.3772029e1
";

const INTEGRATED_STR: &str = "b   etal       integ   
     []         []     
-+----+----+-----------
0    2 2.25 9.3819050e1
1 2.25  2.5 8.6303882e1
2  2.5 2.75 7.5003515e1
3 2.75    3 6.0644157e1
4    3 3.25 4.5233358e1
5 3.25  3.5 3.0727788e1
6  3.5    4 2.8925509e1
7    4  4.5 6.8860146e0
";

const INTEGRATED_MULTIPLE_PDFS_STR: &str = "b   etal       integ    NNPDF31_nlo_as_0118_luxqed 
     []         []                [] [%]           
-+----+----+-----------+-------------+-------------
0    2 2.25 9.3819050e1   9.3819050e1          0.00
1 2.25  2.5 8.6303882e1   8.6303882e1          0.00
2  2.5 2.75 7.5003515e1   7.5003515e1          0.00
3 2.75    3 6.0644157e1   6.0644157e1          0.00
4    3 3.25 4.5233358e1   4.5233358e1          0.00
5 3.25  3.5 3.0727788e1   3.0727788e1          0.00
6  3.5    4 2.8925509e1   2.8925509e1          0.00
7    4  4.5 6.8860146e0   6.8860146e0          0.00
";

const ORDERS_A2_A3_STR: &str = "b   etal    disg/detal 
     []        [pb]    
-+----+----+-----------
0    2 2.25 3.2092052e2
1 2.25  2.5 2.9427151e2
2  2.5 2.75 2.5490261e2
3 2.75    3 2.0561831e2
4    3 3.25 1.5266481e2
5 3.25  3.5 1.0303603e2
6  3.5    4 4.7938981e1
7    4  4.5 1.1075878e1
";

const WRONG_ORDERS_STR: &str = "error: invalid value 'a2a2as2' for '--orders <ORDERS>': unable to parse order; too many couplings in 'a2a2as2'

For more information, try '--help'.
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["convolute", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn default() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
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
        .args([
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
        .args([
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
        .args([
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
        .args([
            "--silence-lhapdf",
            "convolute",
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "--silence-lhapdf",
            "convolute",
            "--orders=a2",
            "data/LHCB_WP_7TEV.pineappl.lz4",
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
            "--silence-lhapdf",
            "convolute",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed/0",
            "NNPDF31_nlo_as_0118_luxqed/1",
            "NNPDF31_nlo_as_0118_luxqed/2",
        ])
        .assert()
        .success()
        .stdout(THREE_PDFS_STR);
}

#[test]
fn wrong_lhaid() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
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
        .args([
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
fn bins_13567() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
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
        .args([
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
        .args([
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
        .args([
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
        .args([
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
