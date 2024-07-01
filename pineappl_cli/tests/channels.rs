use assert_cmd::Command;

const HELP_STR: &str = "Shows the contribution for each partonic channel

Usage: pineappl channels [OPTIONS] <INPUT> <CONV_FUNS>

Arguments:
  <INPUT>      Path to the input grid
  <CONV_FUNS>  LHAPDF ID(s) or name(s) of the PDF(s)/FF(s)

Options:
  -a, --absolute             Show absolute numbers of each contribution
  -l, --limit <LIMIT>        The maximum number of channels displayed [default: 10]
  -i, --integrated           Show integrated numbers (without bin widths) instead of differential ones
      --channels <CHANNELS>  Show only the listed channels
  -o, --orders <ORDERS>      Select orders manually
      --dont-sort            Do not sort the channels according to their size
      --digits-abs <ABS>     Set the number of fractional digits shown for absolute numbers [default: 7]
      --digits-rel <REL>     Set the number of fractional digits shown for relative numbers [default: 2]
  -h, --help                 Print help
";

const DEFAULT_STR: &str = "b   etal    c  size  c  size  c size  c size c size
     []        [%]      [%]      [%]    [%]    [%] 
-+----+----+-+------+-+------+-+-----+-+----+-+----
0    2 2.25 0 111.32 3  -8.05 1 -3.31 4 0.02 2 0.01
1 2.25  2.5 0 112.20 3  -8.85 1 -3.38 4 0.02 2 0.01
2  2.5 2.75 0 113.10 3  -9.56 1 -3.57 4 0.02 2 0.01
3 2.75    3 0 113.98 3 -10.17 1 -3.84 4 0.02 2 0.01
4    3 3.25 0 114.83 3 -10.58 1 -4.27 4 0.01 2 0.01
5 3.25  3.5 0 115.62 3 -10.81 1 -4.84 4 0.02 2 0.01
6  3.5    4 0 116.26 3 -10.48 1 -5.80 2 0.01 4 0.01
7    4  4.5 0 115.88 3  -8.61 1 -7.29 4 0.01 2 0.01
";

const ABSOLUTE_STR: &str =
    "b   etal    c dsig/detal  c  dsig/detal  c  dsig/detal  c  dsig/detal  c  dsig/detal 
     []          [pb]           [pb]           [pb]           [pb]           [pb]    
-+----+----+-+-----------+-+------------+-+------------+-+------------+-+------------
0    2 2.25 0 8.4002759e2 3 -6.0727462e1 1 -2.4969360e1 4 1.7176328e-1 2 8.8565923e-2
1 2.25  2.5 0 7.7448295e2 3 -6.1109036e1 1 -2.3319483e1 4 1.4518685e-1 2 8.3802762e-2
2  2.5 2.75 0 6.7891182e2 3 -5.7385834e1 1 -2.1436419e1 4 1.1534278e-1 2 4.7074109e-2
3 2.75    3 0 5.5341626e2 3 -4.9385114e1 1 -1.8639887e1 4 7.2943823e-2 2 5.8147927e-2
4    3 3.25 0 4.1562095e2 3 -3.8287410e1 1 -1.5462782e1 4 4.9352954e-2 2 3.4452663e-2
5 3.25  3.5 0 2.8427837e2 3 -2.6578788e1 1 -1.1889878e1 4 3.8564621e-2 2 1.8643688e-2
6  3.5    4 0 1.3470473e2 3 -1.2142190e1 1 -6.7199873e0 2 1.3223117e-2 4 1.2734974e-2
7    4  4.5 0 3.1886258e1 3 -2.3686722e0 1 -2.0056686e0 4 3.4154203e-3 2 1.9334685e-3
";

const ABSOLUTE_INTEGRATED_STR: &str =
    "b   etal    c    integ    c    integ     c    integ     c    integ     c    integ    
     []           []             []             []             []             []     
-+----+----+-+-----------+-+------------+-+------------+-+------------+-+------------
0    2 2.25 0 2.1000690e2 3 -1.5181865e1 1 -6.2423401e0 4 4.2940819e-2 2 2.2141481e-2
1 2.25  2.5 0 1.9362074e2 3 -1.5277259e1 1 -5.8298709e0 4 3.6296712e-2 2 2.0950691e-2
2  2.5 2.75 0 1.6972795e2 3 -1.4346459e1 1 -5.3591048e0 4 2.8835695e-2 2 1.1768527e-2
3 2.75    3 0 1.3835407e2 3 -1.2346279e1 1 -4.6599717e0 4 1.8235956e-2 2 1.4536982e-2
4    3 3.25 0 1.0390524e2 3 -9.5718524e0 1 -3.8656954e0 4 1.2338239e-2 2 8.6131657e-3
5 3.25  3.5 0 7.1069592e1 3 -6.6446971e0 1 -2.9724695e0 4 9.6411553e-3 2 4.6609220e-3
6  3.5    4 0 6.7352364e1 3 -6.0710951e0 1 -3.3599937e0 2 6.6115583e-3 4 6.3674872e-3
7    4  4.5 0 1.5943129e1 3 -1.1843361e0 1 -1.0028343e0 4 1.7077102e-3 2 9.6673424e-4
";

const LIMIT_3_STR: &str = "b   etal    c  size  c  size  c size 
     []        [%]      [%]      [%] 
-+----+----+-+------+-+------+-+-----
0    2 2.25 0 111.32 3  -8.05 1 -3.31
1 2.25  2.5 0 112.20 3  -8.85 1 -3.38
2  2.5 2.75 0 113.10 3  -9.56 1 -3.57
3 2.75    3 0 113.98 3 -10.17 1 -3.84
4    3 3.25 0 114.83 3 -10.58 1 -4.27
5 3.25  3.5 0 115.62 3 -10.81 1 -4.84
6  3.5    4 0 116.26 3 -10.48 1 -5.80
7    4  4.5 0 115.88 3  -8.61 1 -7.29
";

const BAD_LIMIT_STR: &str = "error: invalid value '0' for '--limit <LIMIT>': 0 is not in 1..=65535

For more information, try '--help'.
";

const LUMIS_0123_STR: &str = "b   etal    c  size  c  size  c size  c size
     []        [%]      [%]      [%]    [%] 
-+----+----+-+------+-+------+-+-----+-+----
0    2 2.25 0 111.32 3  -8.05 1 -3.31 2 0.01
1 2.25  2.5 0 112.20 3  -8.85 1 -3.38 2 0.01
2  2.5 2.75 0 113.10 3  -9.56 1 -3.57 2 0.01
3 2.75    3 0 113.98 3 -10.17 1 -3.84 2 0.01
4    3 3.25 0 114.83 3 -10.58 1 -4.27 2 0.01
5 3.25  3.5 0 115.62 3 -10.81 1 -4.84 2 0.01
6  3.5    4 0 116.26 3 -10.48 1 -5.80 2 0.01
7    4  4.5 0 115.88 3  -8.61 1 -7.29 2 0.01
";

const ORDERS_A2_AS1A2_STR: &str = "b   etal    c  size  c  size  c size  c size c size
     []        [%]      [%]      [%]    [%]    [%] 
-+----+----+-+------+-+------+-+-----+-+----+-+----
0    2 2.25 0 111.24 3  -7.96 1 -3.27 2 0.00 4 0.00
1 2.25  2.5 0 112.12 3  -8.77 1 -3.35 2 0.00 4 0.00
2  2.5 2.75 0 113.02 3  -9.48 1 -3.54 2 0.00 4 0.00
3 2.75    3 0 113.90 3 -10.09 1 -3.81 2 0.00 4 0.00
4    3 3.25 0 114.74 3 -10.50 1 -4.24 2 0.00 4 0.00
5 3.25  3.5 0 115.54 3 -10.74 1 -4.80 2 0.00 4 0.00
6  3.5    4 0 116.17 3 -10.41 1 -5.76 2 0.00 4 0.00
7    4  4.5 0 115.77 3  -8.54 1 -7.23 2 0.00 4 0.00
";

const DONT_SORT_ABSOLUTE_STR: &str =
    "b   etal    c dsig/detal  c  dsig/detal  c  dsig/detal  c  dsig/detal  c  dsig/detal 
     []          [pb]           [pb]           [pb]           [pb]           [pb]    
-+----+----+-+-----------+-+------------+-+------------+-+------------+-+------------
0    2 2.25 0 8.4002759e2 1 -2.4969360e1 2 8.8565923e-2 3 -6.0727462e1 4 1.7176328e-1
1 2.25  2.5 0 7.7448295e2 1 -2.3319483e1 2 8.3802762e-2 3 -6.1109036e1 4 1.4518685e-1
2  2.5 2.75 0 6.7891182e2 1 -2.1436419e1 2 4.7074109e-2 3 -5.7385834e1 4 1.1534278e-1
3 2.75    3 0 5.5341626e2 1 -1.8639887e1 2 5.8147927e-2 3 -4.9385114e1 4 7.2943823e-2
4    3 3.25 0 4.1562095e2 1 -1.5462782e1 2 3.4452663e-2 3 -3.8287410e1 4 4.9352954e-2
5 3.25  3.5 0 2.8427837e2 1 -1.1889878e1 2 1.8643688e-2 3 -2.6578788e1 4 3.8564621e-2
6  3.5    4 0 1.3470473e2 1 -6.7199873e0 2 1.3223117e-2 3 -1.2142190e1 4 1.2734974e-2
7    4  4.5 0 3.1886258e1 1 -2.0056686e0 2 1.9334685e-3 3 -2.3686722e0 4 3.4154203e-3
";

const DONT_SORT_STR: &str = "b   etal    c  size  c size  c size c  size  c size
     []        [%]      [%]    [%]     [%]     [%] 
-+----+----+-+------+-+-----+-+----+-+------+-+----
0    2 2.25 0 111.32 1 -3.31 2 0.01 3  -8.05 4 0.02
1 2.25  2.5 0 112.20 1 -3.38 2 0.01 3  -8.85 4 0.02
2  2.5 2.75 0 113.10 1 -3.57 2 0.01 3  -9.56 4 0.02
3 2.75    3 0 113.98 1 -3.84 2 0.01 3 -10.17 4 0.02
4    3 3.25 0 114.83 1 -4.27 2 0.01 3 -10.58 4 0.01
5 3.25  3.5 0 115.62 1 -4.84 2 0.01 3 -10.81 4 0.02
6  3.5    4 0 116.26 1 -5.80 2 0.01 3 -10.48 4 0.01
7    4  4.5 0 115.88 1 -7.29 2 0.01 3  -8.61 4 0.01
";

const MISSING_CONV_FUN_STR: &str = "error: the following required arguments were not provided:
  <CONV_FUNS>

Usage: pineappl channels <INPUT> <CONV_FUNS>

For more information, try '--help'.
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["channels", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn default() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "channels",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}

#[test]
fn missing_conv_fun() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["channels", "../test-data/LHCB_WP_7TEV.pineappl.lz4"])
        .assert()
        .failure()
        .stderr(MISSING_CONV_FUN_STR)
        .stdout("");
}

#[test]
fn absolute() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "channels",
            "--absolute",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
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
            "channels",
            "--absolute",
            "--integrated",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ABSOLUTE_INTEGRATED_STR);
}

#[test]
fn limit_3() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "channels",
            "--limit=3",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(LIMIT_3_STR);
}

#[test]
fn bad_limit() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "channels",
            "--limit=0",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .failure()
        .stderr(BAD_LIMIT_STR);
}

#[test]
fn channels_0123() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "channels",
            "--channels=0-3",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(LUMIS_0123_STR);
}

#[test]
fn orders_a2_as1a2() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "channels",
            "--orders=a2,as1a2",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(ORDERS_A2_AS1A2_STR);
}

#[test]
fn dont_sort_absolute() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "channels",
            "--absolute",
            "--dont-sort",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DONT_SORT_ABSOLUTE_STR);
}

#[test]
fn dont_sort() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "channels",
            "--dont-sort",
            "../test-data/LHCB_WP_7TEV.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DONT_SORT_STR);
}
