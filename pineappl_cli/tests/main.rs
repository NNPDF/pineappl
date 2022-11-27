use assert_cmd::Command;
use git_version::git_version;

const HELP_STR: &str = "
Christopher Schwan <handgranaten-herbert@posteo.de>
Read, write, and query PineAPPL grids

USAGE:
    pineappl [OPTIONS] <SUBCOMMAND>

OPTIONS:
    -h, --help              Print help information
        --silence-lhapdf    Prevents LHAPDF from printing banners
    -V, --version           Print version information

SUBCOMMANDS:
    analyze      Perform various analyses with grids
    channels     Shows the contribution for each partonic channel
    convolute    Convolutes a PineAPPL grid with a PDF set
    delete       Deletes parts from a PineAPPL grid
    diff         Compares the numerical content of two grids with each other
    evolve       Evolve a grid with an evolution kernel operator to an FK table
    import       Converts APPLgrid/fastNLO/FastKernel files to PineAPPL grids
    info         Shows information about the grid
    merge        Merges one or more PineAPPL grids together
    obl          Shows information about orders (o), bins (b), or luminosities (l) of a grid
    ops          A collection of various modifying operations on grids
    optimize     Optimizes the internal data structure to minimize memory usage
    orders       Shows the predictions for all bin for each order separately
    pdfunc       Calculates PDF uncertainties
    plot         Creates a matplotlib script plotting the contents of the grid
    pull         Calculates the pull between two different PDF sets
    remap        Modifies the bin dimensions, widths and normalizations
    set          Modifies the internal key-value storage
    subgrids     Print information about the internal subgrid types
    sum          Sums two or more bins of a grid together
    upgrade      Converts the file format to the most recent version
";

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .arg("--help")
        .assert()
        .success()
        .stdout(format!(
            "pineappl {}{}",
            git_version!(
                args = ["--always", "--dirty", "--long", "--tags"],
                cargo_prefix = "cargo:",
                fallback = "unknown"
            ),
            HELP_STR
        ));
}
