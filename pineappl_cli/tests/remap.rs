use assert_cmd::Command;
use assert_fs::NamedTempFile;

const HELP_STR: &str = "Modifies the bin dimensions, widths and normalizations

Usage: pineappl remap [OPTIONS] <INPUT> <OUTPUT> <REMAPPING>

Arguments:
  <INPUT>      Path to the input grid
  <OUTPUT>     Path of the modified PineAPPL file
  <REMAPPING>  Remapping string. See <https://nnpdf.github.io/pineappl/docs/cli-reference.html> for full reference

Options:
      --ignore-obs-norm <OBS1,OBS2,...>
          Ignore the given observables for differential normalization
      --norm <NORM>
          Normalization factor in addition to the given bin widths [default: 1.0]
  -h, --help
          Print help
";

const DEFAULT_STR: &str = "b etal  x2  x3  disg/detal  scale uncertainty
   []   []  []     [pb]            [%]       
-+--+--+-+-+-+-+-----------+--------+--------
0  0  1 0 2 1 2 1.8763810e1    -3.77     2.71
1  0  1 0 2 2 3 1.7260776e1    -3.79     2.80
2  0  1 0 2 3 4 1.5000703e1    -3.78     2.86
3  0  1 0 2 4 5 1.2128831e1    -3.77     2.92
4  0  1 2 4 1 2 9.0466716e0    -3.74     2.95
5  1  2 0 2 8 9 6.1455576e0    -3.71     2.98
6  1  2 2 4 3 4 5.7851018e0    -3.63     2.97
7  1  2 2 4 4 5 1.3772029e0    -3.46     2.85
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["remap", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn default() {
    let output = NamedTempFile::new("optimized.pineappl.lz4").unwrap();

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "remap",
            "--ignore-obs-norm=2",
            "--norm=5",
            "data/LHCB_WP_7TEV.pineappl.lz4",
            output.path().to_str().unwrap(),
            "0,1,2;0,2,4;1,2,3,4,5|:3|5:1,2,3,4,5,8,9|2:2",
        ])
        .assert()
        .success()
        .stdout("");

    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "--silence-lhapdf",
            "convolute",
            output.path().to_str().unwrap(),
            "NNPDF31_nlo_as_0118_luxqed",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}
