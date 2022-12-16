use assert_cmd::Command;

const HELP_STR: &str = "pineappl-obl 
Shows information about orders (o), bins (b), or luminosities (l) of a grid

USAGE:
    pineappl obl <--orders|--orders-spaces|--orders-long|--bins|--lumis|--fktable> <INPUT>

ARGS:
    <INPUT>    Path to the input grid

OPTIONS:
    -b, --bins             Show the bins of a grid
        --fktable          Check if input is an FK table
    -h, --help             Print help information
    -l, --lumis            Show the luminsities a grid
    -o, --orders           Show the orders of a grid, stripping zero powers
        --orders-long      Show the orders of a grid, including zero powers
        --orders-spaces    Show the orders of a grid, replacing zero powers with spaces
";

const BINS_STR: &str = "b   etal    norm
-+----+----+----
0    2 2.25 0.25
1 2.25  2.5 0.25
2  2.5 2.75 0.25
3 2.75    3 0.25
4    3 3.25 0.25
5 3.25  3.5 0.25
6  3.5    4  0.5
7    4  4.5  0.5
";

const LUMIS_STR: &str = "l    entry        entry
-+------------+------------
0 1 × ( 2, -1) 1 × ( 4, -3)
1 1 × ( 0, -3) 1 × ( 0, -1)
2 1 × (22, -3) 1 × (22, -1)
3 1 × ( 2,  0) 1 × ( 4,  0)
4 1 × ( 2, 22) 1 × ( 4, 22)
";

const ORDERS_STR: &str = "o      order
-+----------------
0 O(a^2)
1 O(as^1 a^2)
2 O(as^1 a^2 lr^1)
3 O(as^1 a^2 lf^1)
4 O(a^3)
5 O(a^3 lr^1)
6 O(a^3 lf^1)
";

const ORDERS_LONG_STR: &str = "o         order
-+---------------------
0 O(as^0 a^2 lr^0 lf^0)
1 O(as^1 a^2 lr^0 lf^0)
2 O(as^1 a^2 lr^1 lf^0)
3 O(as^1 a^2 lr^0 lf^1)
4 O(as^0 a^3 lr^0 lf^0)
5 O(as^0 a^3 lr^1 lf^0)
6 O(as^0 a^3 lr^0 lf^1)
";

const ORDERS_SPACES_STR: &str = "o         order
-+---------------------
0 O(     a^2          )
1 O(as^1 a^2          )
2 O(as^1 a^2 lr^1     )
3 O(as^1 a^2      lf^1)
4 O(     a^3          )
5 O(     a^3 lr^1     )
6 O(     a^3      lf^1)
";

const FKTABLE_STR: &str = "no
multiple orders detected
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["obl", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn bins() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["obl", "--bins", "data/LHCB_WP_7TEV.pineappl.lz4"])
        .assert()
        .success()
        .stdout(BINS_STR);
}

#[test]
fn lumis() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["obl", "--lumis", "data/LHCB_WP_7TEV.pineappl.lz4"])
        .assert()
        .success()
        .stdout(LUMIS_STR);
}

#[test]
fn orders() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["obl", "--orders", "data/LHCB_WP_7TEV.pineappl.lz4"])
        .assert()
        .success()
        .stdout(ORDERS_STR);
}

#[test]
fn orders_long() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["obl", "--orders-long", "data/LHCB_WP_7TEV.pineappl.lz4"])
        .assert()
        .success()
        .stdout(ORDERS_LONG_STR);
}

#[test]
fn orders_spaces() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["obl", "--orders-spaces", "data/LHCB_WP_7TEV.pineappl.lz4"])
        .assert()
        .success()
        .stdout(ORDERS_SPACES_STR);
}

#[test]
fn fktable() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(&["obl", "--fktable", "data/LHCB_WP_7TEV.pineappl.lz4"])
        .assert()
        .success()
        .stdout(FKTABLE_STR);
}
