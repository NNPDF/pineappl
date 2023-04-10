use assert_cmd::Command;

const HELP_STR: &str = "Read, write, and query PineAPPL grids

Usage: pineappl [OPTIONS] <COMMAND>

Commands:
  analyze    Perform various analyses with grids
  channels   Shows the contribution for each partonic channel
  convolute  Convolutes a PineAPPL grid with a PDF set
  diff       Compares the numerical content of two grids with each other
  evolve     Evolve a grid with an evolution kernel operator to an FK table
  help       Display a manpage for selected subcommands
  import     Converts APPLgrid/fastNLO/FastKernel files to PineAPPL grids
  merge      Merges one or more PineAPPL grids together
  orders     Shows the predictions for all bin for each order separately
  pdfunc     Calculates PDF uncertainties
  plot       Creates a matplotlib script plotting the contents of the grid
  pull       Calculates the pull between two different PDF sets
  read       Read out information of a grid
  subgrids   Print information about the internal subgrid types
  write      Write a grid modified by various operations

Options:
      --silence-lhapdf  Prevents LHAPDF from printing banners
      --force-positive  Forces negative PDF values to zero
  -h, --help            Print help
  -V, --version         Print version
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
