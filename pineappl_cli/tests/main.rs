#![allow(missing_docs)]

use assert_cmd::Command;

const HELP_STR: &str = "Read, write, and query PineAPPL grids

Usage: pineappl [OPTIONS] <COMMAND>

Commands:
  analyze   Perform various analyses with grids
  channels  Shows the contribution for each partonic channel
  convolve  Convolutes a PineAPPL grid with a PDF set
  diff      Compares the numerical content of two grids with each other
  evolve    Evolve a grid with an evolution kernel operator to an FK table
  export    Converts PineAPPL grids to APPLgrid files
  help      Display a manpage for selected subcommands
  import    Converts APPLgrid/fastNLO/FastKernel files to PineAPPL grids
  merge     Merges one or more PineAPPL grids together
  orders    Shows the predictions for all bin for each order separately
  plot      Creates a matplotlib script plotting the contents of the grid
  pull      Calculates the pull between two different PDF sets
  read      Read out information of a grid
  subgrids  Print information about the internal subgrid types
  uncert    Calculate scale and convolution function uncertainties
  write     Write a grid modified by various operations

Options:
      --lhapdf-banner          Allow LHAPDF to print banners
      --force-positive         Forces negative PDF values to zero
      --allow-extrapolation    Allow extrapolation of PDFs outside their region of validity
      --use-alphas-from <IDX>  Choose the PDF/FF set for the strong coupling [default: 0]
  -h, --help                   Print help
  -V, --version                Print version
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .arg("--help")
        .assert()
        .success()
        .stdout(HELP_STR);
}
