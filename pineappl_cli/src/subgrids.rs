use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{ArgGroup, Parser, ValueHint};
use pineappl::subgrid::Mu2;
use pineappl::subgrid::{Subgrid, SubgridEnum};
use prettytable::{cell, row};
use std::path::PathBuf;

/// Print information about the internal subgrid types.
#[derive(Parser)]
#[clap(group = ArgGroup::new("show").multiple(true).required(true))]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Show empty subgrids.
    #[clap(long = "show-empty")]
    show_empty: bool,
    /// Show the subgrid type.
    #[clap(group = "show", long = "type")]
    type_: bool,
    /// Show the renormalization grid values.
    #[clap(group = "show", long)]
    mur: bool,
    /// Show the factorization grid values.
    #[clap(group = "show", long)]
    muf: bool,
    /// Show the x1 grid values.
    #[clap(group = "show", long)]
    x1: bool,
    /// Show the x2 grid values.
    #[clap(group = "show", long)]
    x2: bool,
    /// Show grid statistics (figures are the number of entries).
    #[clap(group = "show", long)]
    stats: bool,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let grid = helpers::read_grid(&self.input)?;
        let mut table = helpers::create_table();
        let mut titles = row![c => "o", "b", "l"];

        if self.type_ {
            titles.add_cell(cell!(c->"type"));
        }
        if self.mur {
            titles.add_cell(cell!(c->"mur"));
        }
        if self.muf {
            titles.add_cell(cell!(c->"muf"));
        }
        if self.x1 {
            titles.add_cell(cell!(c->"x1"));
        }
        if self.x2 {
            titles.add_cell(cell!(c->"x2"));
        }
        if self.stats {
            titles.add_cell(cell!(c->"total"));
            titles.add_cell(cell!(c->"allocated"));
            titles.add_cell(cell!(c->"zeros"));
            titles.add_cell(cell!(c->"overhead"));
        }
        table.set_titles(titles);

        for ((order, bin, lumi), subgrid) in grid.subgrids().indexed_iter() {
            if !self.show_empty && subgrid.is_empty() {
                continue;
            }

            let row = table.add_empty_row();

            row.add_cell(cell!(l->&format!("{}", order)));
            row.add_cell(cell!(l->&format!("{}", bin)));
            row.add_cell(cell!(l->&format!("{}", lumi)));

            if self.type_ {
                row.add_cell(cell!(l->
                    match subgrid {
                        SubgridEnum::LagrangeSubgridV1(_) => "LagrangeSubgridV1",
                        SubgridEnum::NtupleSubgridV1(_) => "NtupleSubgridV1",
                        SubgridEnum::LagrangeSparseSubgridV1(_) => "LagrangeSparseSubgridV1",
                        SubgridEnum::LagrangeSubgridV2(_) => "LagrangeSubgridV2",
                        SubgridEnum::ImportOnlySubgridV1(_) => "ImportOnlySubgridV1",
                        SubgridEnum::ImportOnlySubgridV2(_) => "ImportOnlySubgridV2",
                        SubgridEnum::EmptySubgridV1(_) => "EmptySubgridV1",
                    }
                ));
            }
            if self.mur {
                let values: Vec<_> = subgrid
                    .mu2_grid()
                    .iter()
                    .map(|Mu2 { ren, fac: _ }| format!("{:.3}", ren.sqrt()))
                    .collect();

                row.add_cell(cell!(l->&values.join(", ")));
            }
            if self.muf {
                let values: Vec<_> = subgrid
                    .mu2_grid()
                    .iter()
                    .map(|Mu2 { ren: _, fac }| format!("{:.3}", fac.sqrt()))
                    .collect();

                row.add_cell(cell!(l->&values.join(", ")));
            }
            if self.x1 {
                let values: Vec<_> = subgrid
                    .x1_grid()
                    .iter()
                    .map(|x| format!("{:.3e}", x))
                    .collect();

                row.add_cell(cell!(l->&values.join(", ")));
            }
            if self.x2 {
                let values: Vec<_> = subgrid
                    .x2_grid()
                    .iter()
                    .map(|x| format!("{:.3e}", x))
                    .collect();

                row.add_cell(cell!(l->&values.join(", ")));
            }
            if self.stats {
                let stats = subgrid.stats();
                row.add_cell(cell!(r->&stats.total.to_string()));
                row.add_cell(cell!(r->&stats.allocated.to_string()));
                row.add_cell(cell!(r->&stats.zeros.to_string()));
                row.add_cell(cell!(r->&stats.overhead.to_string()));
            }
        }

        table.printstd();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-subgrids 
Print information about the internal subgrid types

USAGE:
    pineappl subgrids [OPTIONS] <--type|--mur|--muf|--x1|--x2|--stats> <INPUT>

ARGS:
    <INPUT>    Path to the input grid

OPTIONS:
    -h, --help          Print help information
        --muf           Show the factorization grid values
        --mur           Show the renormalization grid values
        --show-empty    Show empty subgrids
        --stats         Show grid statistics (figures are the number of entries)
        --type          Show the subgrid type
        --x1            Show the x1 grid values
        --x2            Show the x2 grid values
";

    const MUF_STR: &str = "o b l  muf
-+-+-+------
0 0 0 80.352
0 1 0 80.352
0 2 0 80.352
0 3 0 80.352
0 4 0 80.352
0 5 0 80.352
0 6 0 80.352
0 7 0 80.352
1 0 0 80.352
1 0 1 80.352
1 0 3 80.352
1 1 0 80.352
1 1 1 80.352
1 1 3 80.352
1 2 0 80.352
1 2 1 80.352
1 2 3 80.352
1 3 0 80.352
1 3 1 80.352
1 3 3 80.352
1 4 0 80.352
1 4 1 80.352
1 4 3 80.352
1 5 0 80.352
1 5 1 80.352
1 5 3 80.352
1 6 0 80.352
1 6 1 80.352
1 6 3 80.352
1 7 0 80.352
1 7 1 80.352
1 7 3 80.352
3 0 0 80.352
3 0 1 80.352
3 0 3 80.352
3 1 0 80.352
3 1 1 80.352
3 1 3 80.352
3 2 0 80.352
3 2 1 80.352
3 2 3 80.352
3 3 0 80.352
3 3 1 80.352
3 3 3 80.352
3 4 0 80.352
3 4 1 80.352
3 4 3 80.352
3 5 0 80.352
3 5 1 80.352
3 5 3 80.352
3 6 0 80.352
3 6 1 80.352
3 6 3 80.352
3 7 0 80.352
3 7 1 80.352
3 7 3 80.352
4 0 0 80.352
4 0 2 80.352
4 0 4 80.352
4 1 0 80.352
4 1 2 80.352
4 1 4 80.352
4 2 0 80.352
4 2 2 80.352
4 2 4 80.352
4 3 0 80.352
4 3 2 80.352
4 3 4 80.352
4 4 0 80.352
4 4 2 80.352
4 4 4 80.352
4 5 0 80.352
4 5 2 80.352
4 5 4 80.352
4 6 0 80.352
4 6 2 80.352
4 6 4 80.352
4 7 0 80.352
4 7 2 80.352
4 7 4 80.352
6 0 0 80.352
6 0 2 80.352
6 0 4 80.352
6 1 0 80.352
6 1 2 80.352
6 1 4 80.352
6 2 0 80.352
6 2 2 80.352
6 2 4 80.352
6 3 0 80.352
6 3 2 80.352
6 3 4 80.352
6 4 0 80.352
6 4 2 80.352
6 4 4 80.352
6 5 0 80.352
6 5 2 80.352
6 5 4 80.352
6 6 0 80.352
6 6 2 80.352
6 6 4 80.352
6 7 0 80.352
6 7 2 80.352
6 7 4 80.352
";

    const MUR_STR: &str = "o b l  mur
-+-+-+------
0 0 0 80.352
0 1 0 80.352
0 2 0 80.352
0 3 0 80.352
0 4 0 80.352
0 5 0 80.352
0 6 0 80.352
0 7 0 80.352
1 0 0 80.352
1 0 1 80.352
1 0 3 80.352
1 1 0 80.352
1 1 1 80.352
1 1 3 80.352
1 2 0 80.352
1 2 1 80.352
1 2 3 80.352
1 3 0 80.352
1 3 1 80.352
1 3 3 80.352
1 4 0 80.352
1 4 1 80.352
1 4 3 80.352
1 5 0 80.352
1 5 1 80.352
1 5 3 80.352
1 6 0 80.352
1 6 1 80.352
1 6 3 80.352
1 7 0 80.352
1 7 1 80.352
1 7 3 80.352
3 0 0 80.352
3 0 1 80.352
3 0 3 80.352
3 1 0 80.352
3 1 1 80.352
3 1 3 80.352
3 2 0 80.352
3 2 1 80.352
3 2 3 80.352
3 3 0 80.352
3 3 1 80.352
3 3 3 80.352
3 4 0 80.352
3 4 1 80.352
3 4 3 80.352
3 5 0 80.352
3 5 1 80.352
3 5 3 80.352
3 6 0 80.352
3 6 1 80.352
3 6 3 80.352
3 7 0 80.352
3 7 1 80.352
3 7 3 80.352
4 0 0 80.352
4 0 2 80.352
4 0 4 80.352
4 1 0 80.352
4 1 2 80.352
4 1 4 80.352
4 2 0 80.352
4 2 2 80.352
4 2 4 80.352
4 3 0 80.352
4 3 2 80.352
4 3 4 80.352
4 4 0 80.352
4 4 2 80.352
4 4 4 80.352
4 5 0 80.352
4 5 2 80.352
4 5 4 80.352
4 6 0 80.352
4 6 2 80.352
4 6 4 80.352
4 7 0 80.352
4 7 2 80.352
4 7 4 80.352
6 0 0 80.352
6 0 2 80.352
6 0 4 80.352
6 1 0 80.352
6 1 2 80.352
6 1 4 80.352
6 2 0 80.352
6 2 2 80.352
6 2 4 80.352
6 3 0 80.352
6 3 2 80.352
6 3 4 80.352
6 4 0 80.352
6 4 2 80.352
6 4 4 80.352
6 5 0 80.352
6 5 2 80.352
6 5 4 80.352
6 6 0 80.352
6 6 2 80.352
6 6 4 80.352
6 7 0 80.352
6 7 2 80.352
6 7 4 80.352
";

    const STATS_STR: &str = "o b l total allocated zeros overhead
-+-+-+-----+---------+-----+--------
0 0 0  2500       987    24      102
0 1 0  2500       974    62      102
0 2 0  2500      1002    48      102
0 3 0  2500       910    36      102
0 4 0  2500       932    56      102
0 5 0  2500       845    31      102
0 6 0  2500       907    62      102
0 7 0  2500       696    32      102
1 0 0  2500      1012     0      102
1 0 1  2500      1002     0      102
1 0 3  2500      1011     0      102
1 1 0  2500      1014     0      102
1 1 1  2500      1007     0      102
1 1 3  2500      1011     0      102
1 2 0  2500      1030     0      102
1 2 1  2500      1020     4      102
1 2 3  2500      1027     0      102
1 3 0  2500      1028     0      102
1 3 1  2500      1010     0      102
1 3 3  2500      1027     0      102
1 4 0  2500      1031     0      102
1 4 1  2500       996     5      102
1 4 3  2500      1027     0      102
1 5 0  2500      1037     0      102
1 5 1  2500       980     9      102
1 5 3  2500      1030     0      102
1 6 0  2500      1031     0      102
1 6 1  2500      1009     4      102
1 6 3  2500      1028     0      102
1 7 0  2500       982     0      102
1 7 1  2500       972    14      102
1 7 3  2500       971     0      102
3 0 0  2500      1012     0      102
3 0 1  2500       954     5      102
3 0 3  2500       996     1      102
3 1 0  2500      1011    13      102
3 1 1  2500       898     0      102
3 1 3  2500       997    14      102
3 2 0  2500      1029     2      102
3 2 1  2500       963     7      102
3 2 3  2500      1015     0      102
3 3 0  2500      1028     0      102
3 3 1  2500       947     2      102
3 3 3  2500      1025     3      102
3 4 0  2500      1031     0      102
3 4 1  2500       975    21      102
3 4 3  2500      1017     0      102
3 5 0  2500      1037     0      102
3 5 1  2500       970    32      102
3 5 3  2500      1025     0      102
3 6 0  2500      1031     0      102
3 6 1  2500      1005    29      102
3 6 3  2500      1019     0      102
3 7 0  2500       982     0      102
3 7 1  2500       970    72      102
3 7 3  2500       925     3      102
4 0 0  2500      1012     6      102
4 0 2  2500       839     6      102
4 0 4  2500       837     4      102
4 1 0  2500      1003    14      102
4 1 2  2500       822    18      102
4 1 4  2500       905    11      102
4 2 0  2500      1025     1      102
4 2 2  2500       818     7      102
4 2 4  2500       896     3      102
4 3 0  2500      1008     2      102
4 3 2  2500       830    14      102
4 3 4  2500       883    10      102
4 4 0  2500      1021    12      102
4 4 2  2500       786    19      102
4 4 4  2500       885    18      102
4 5 0  2500      1021     7      102
4 5 2  2500       778    21      102
4 5 4  2500       814    12      102
4 6 0  2500      1031     3      102
4 6 2  2500       790    64      102
4 6 4  2500       855    27      102
4 7 0  2500       977    21      102
4 7 2  2500       806   174      102
4 7 4  2500       621     1      102
6 0 0  2500      1000    15      102
6 0 2  2500       619     7      102
6 0 4  2500       736     8      102
6 1 0  2500       985    27      102
6 1 2  2500       633    10      102
6 1 4  2500       793     5      102
6 2 0  2500      1020     5      102
6 2 2  2500       690    58      102
6 2 4  2500       765    17      102
6 3 0  2500       999     9      102
6 3 2  2500       719    72      102
6 3 4  2500       786    12      102
6 4 0  2500      1015    12      102
6 4 2  2500       723    97      102
6 4 4  2500       773    19      102
6 5 0  2500      1015     7      102
6 5 2  2500       759   115      102
6 5 4  2500       723    13      102
6 6 0  2500      1031     3      102
6 6 2  2500       727   106      102
6 6 4  2500       788    33      102
6 7 0  2500       977    22      102
6 7 2  2500       774   205      102
6 7 4  2500       537     0      102
";

    const TYPE_STR: &str = "o b l        type
-+-+-+-------------------
0 0 0 ImportOnlySubgridV1
0 1 0 ImportOnlySubgridV1
0 2 0 ImportOnlySubgridV1
0 3 0 ImportOnlySubgridV1
0 4 0 ImportOnlySubgridV1
0 5 0 ImportOnlySubgridV1
0 6 0 ImportOnlySubgridV1
0 7 0 ImportOnlySubgridV1
1 0 0 ImportOnlySubgridV1
1 0 1 ImportOnlySubgridV1
1 0 3 ImportOnlySubgridV1
1 1 0 ImportOnlySubgridV1
1 1 1 ImportOnlySubgridV1
1 1 3 ImportOnlySubgridV1
1 2 0 ImportOnlySubgridV1
1 2 1 ImportOnlySubgridV1
1 2 3 ImportOnlySubgridV1
1 3 0 ImportOnlySubgridV1
1 3 1 ImportOnlySubgridV1
1 3 3 ImportOnlySubgridV1
1 4 0 ImportOnlySubgridV1
1 4 1 ImportOnlySubgridV1
1 4 3 ImportOnlySubgridV1
1 5 0 ImportOnlySubgridV1
1 5 1 ImportOnlySubgridV1
1 5 3 ImportOnlySubgridV1
1 6 0 ImportOnlySubgridV1
1 6 1 ImportOnlySubgridV1
1 6 3 ImportOnlySubgridV1
1 7 0 ImportOnlySubgridV1
1 7 1 ImportOnlySubgridV1
1 7 3 ImportOnlySubgridV1
3 0 0 ImportOnlySubgridV1
3 0 1 ImportOnlySubgridV1
3 0 3 ImportOnlySubgridV1
3 1 0 ImportOnlySubgridV1
3 1 1 ImportOnlySubgridV1
3 1 3 ImportOnlySubgridV1
3 2 0 ImportOnlySubgridV1
3 2 1 ImportOnlySubgridV1
3 2 3 ImportOnlySubgridV1
3 3 0 ImportOnlySubgridV1
3 3 1 ImportOnlySubgridV1
3 3 3 ImportOnlySubgridV1
3 4 0 ImportOnlySubgridV1
3 4 1 ImportOnlySubgridV1
3 4 3 ImportOnlySubgridV1
3 5 0 ImportOnlySubgridV1
3 5 1 ImportOnlySubgridV1
3 5 3 ImportOnlySubgridV1
3 6 0 ImportOnlySubgridV1
3 6 1 ImportOnlySubgridV1
3 6 3 ImportOnlySubgridV1
3 7 0 ImportOnlySubgridV1
3 7 1 ImportOnlySubgridV1
3 7 3 ImportOnlySubgridV1
4 0 0 ImportOnlySubgridV1
4 0 2 ImportOnlySubgridV1
4 0 4 ImportOnlySubgridV1
4 1 0 ImportOnlySubgridV1
4 1 2 ImportOnlySubgridV1
4 1 4 ImportOnlySubgridV1
4 2 0 ImportOnlySubgridV1
4 2 2 ImportOnlySubgridV1
4 2 4 ImportOnlySubgridV1
4 3 0 ImportOnlySubgridV1
4 3 2 ImportOnlySubgridV1
4 3 4 ImportOnlySubgridV1
4 4 0 ImportOnlySubgridV1
4 4 2 ImportOnlySubgridV1
4 4 4 ImportOnlySubgridV1
4 5 0 ImportOnlySubgridV1
4 5 2 ImportOnlySubgridV1
4 5 4 ImportOnlySubgridV1
4 6 0 ImportOnlySubgridV1
4 6 2 ImportOnlySubgridV1
4 6 4 ImportOnlySubgridV1
4 7 0 ImportOnlySubgridV1
4 7 2 ImportOnlySubgridV1
4 7 4 ImportOnlySubgridV1
6 0 0 ImportOnlySubgridV1
6 0 2 ImportOnlySubgridV1
6 0 4 ImportOnlySubgridV1
6 1 0 ImportOnlySubgridV1
6 1 2 ImportOnlySubgridV1
6 1 4 ImportOnlySubgridV1
6 2 0 ImportOnlySubgridV1
6 2 2 ImportOnlySubgridV1
6 2 4 ImportOnlySubgridV1
6 3 0 ImportOnlySubgridV1
6 3 2 ImportOnlySubgridV1
6 3 4 ImportOnlySubgridV1
6 4 0 ImportOnlySubgridV1
6 4 2 ImportOnlySubgridV1
6 4 4 ImportOnlySubgridV1
6 5 0 ImportOnlySubgridV1
6 5 2 ImportOnlySubgridV1
6 5 4 ImportOnlySubgridV1
6 6 0 ImportOnlySubgridV1
6 6 2 ImportOnlySubgridV1
6 6 4 ImportOnlySubgridV1
6 7 0 ImportOnlySubgridV1
6 7 2 ImportOnlySubgridV1
6 7 4 ImportOnlySubgridV1
";

    const TYPE_SHOW_EMPTY_STR: &str = "o b l        type
-+-+-+-------------------
0 0 0 ImportOnlySubgridV1
0 0 1 EmptySubgridV1
0 0 2 EmptySubgridV1
0 0 3 EmptySubgridV1
0 0 4 EmptySubgridV1
0 1 0 ImportOnlySubgridV1
0 1 1 EmptySubgridV1
0 1 2 EmptySubgridV1
0 1 3 EmptySubgridV1
0 1 4 EmptySubgridV1
0 2 0 ImportOnlySubgridV1
0 2 1 EmptySubgridV1
0 2 2 EmptySubgridV1
0 2 3 EmptySubgridV1
0 2 4 EmptySubgridV1
0 3 0 ImportOnlySubgridV1
0 3 1 EmptySubgridV1
0 3 2 EmptySubgridV1
0 3 3 EmptySubgridV1
0 3 4 EmptySubgridV1
0 4 0 ImportOnlySubgridV1
0 4 1 EmptySubgridV1
0 4 2 EmptySubgridV1
0 4 3 EmptySubgridV1
0 4 4 EmptySubgridV1
0 5 0 ImportOnlySubgridV1
0 5 1 EmptySubgridV1
0 5 2 EmptySubgridV1
0 5 3 EmptySubgridV1
0 5 4 EmptySubgridV1
0 6 0 ImportOnlySubgridV1
0 6 1 EmptySubgridV1
0 6 2 EmptySubgridV1
0 6 3 EmptySubgridV1
0 6 4 EmptySubgridV1
0 7 0 ImportOnlySubgridV1
0 7 1 EmptySubgridV1
0 7 2 EmptySubgridV1
0 7 3 EmptySubgridV1
0 7 4 EmptySubgridV1
1 0 0 ImportOnlySubgridV1
1 0 1 ImportOnlySubgridV1
1 0 2 EmptySubgridV1
1 0 3 ImportOnlySubgridV1
1 0 4 EmptySubgridV1
1 1 0 ImportOnlySubgridV1
1 1 1 ImportOnlySubgridV1
1 1 2 EmptySubgridV1
1 1 3 ImportOnlySubgridV1
1 1 4 EmptySubgridV1
1 2 0 ImportOnlySubgridV1
1 2 1 ImportOnlySubgridV1
1 2 2 EmptySubgridV1
1 2 3 ImportOnlySubgridV1
1 2 4 EmptySubgridV1
1 3 0 ImportOnlySubgridV1
1 3 1 ImportOnlySubgridV1
1 3 2 EmptySubgridV1
1 3 3 ImportOnlySubgridV1
1 3 4 EmptySubgridV1
1 4 0 ImportOnlySubgridV1
1 4 1 ImportOnlySubgridV1
1 4 2 EmptySubgridV1
1 4 3 ImportOnlySubgridV1
1 4 4 EmptySubgridV1
1 5 0 ImportOnlySubgridV1
1 5 1 ImportOnlySubgridV1
1 5 2 EmptySubgridV1
1 5 3 ImportOnlySubgridV1
1 5 4 EmptySubgridV1
1 6 0 ImportOnlySubgridV1
1 6 1 ImportOnlySubgridV1
1 6 2 EmptySubgridV1
1 6 3 ImportOnlySubgridV1
1 6 4 EmptySubgridV1
1 7 0 ImportOnlySubgridV1
1 7 1 ImportOnlySubgridV1
1 7 2 EmptySubgridV1
1 7 3 ImportOnlySubgridV1
1 7 4 EmptySubgridV1
2 0 0 EmptySubgridV1
2 0 1 EmptySubgridV1
2 0 2 EmptySubgridV1
2 0 3 EmptySubgridV1
2 0 4 EmptySubgridV1
2 1 0 EmptySubgridV1
2 1 1 EmptySubgridV1
2 1 2 EmptySubgridV1
2 1 3 EmptySubgridV1
2 1 4 EmptySubgridV1
2 2 0 EmptySubgridV1
2 2 1 EmptySubgridV1
2 2 2 EmptySubgridV1
2 2 3 EmptySubgridV1
2 2 4 EmptySubgridV1
2 3 0 EmptySubgridV1
2 3 1 EmptySubgridV1
2 3 2 EmptySubgridV1
2 3 3 EmptySubgridV1
2 3 4 EmptySubgridV1
2 4 0 EmptySubgridV1
2 4 1 EmptySubgridV1
2 4 2 EmptySubgridV1
2 4 3 EmptySubgridV1
2 4 4 EmptySubgridV1
2 5 0 EmptySubgridV1
2 5 1 EmptySubgridV1
2 5 2 EmptySubgridV1
2 5 3 EmptySubgridV1
2 5 4 EmptySubgridV1
2 6 0 EmptySubgridV1
2 6 1 EmptySubgridV1
2 6 2 EmptySubgridV1
2 6 3 EmptySubgridV1
2 6 4 EmptySubgridV1
2 7 0 EmptySubgridV1
2 7 1 EmptySubgridV1
2 7 2 EmptySubgridV1
2 7 3 EmptySubgridV1
2 7 4 EmptySubgridV1
3 0 0 ImportOnlySubgridV1
3 0 1 ImportOnlySubgridV1
3 0 2 EmptySubgridV1
3 0 3 ImportOnlySubgridV1
3 0 4 EmptySubgridV1
3 1 0 ImportOnlySubgridV1
3 1 1 ImportOnlySubgridV1
3 1 2 EmptySubgridV1
3 1 3 ImportOnlySubgridV1
3 1 4 EmptySubgridV1
3 2 0 ImportOnlySubgridV1
3 2 1 ImportOnlySubgridV1
3 2 2 EmptySubgridV1
3 2 3 ImportOnlySubgridV1
3 2 4 EmptySubgridV1
3 3 0 ImportOnlySubgridV1
3 3 1 ImportOnlySubgridV1
3 3 2 EmptySubgridV1
3 3 3 ImportOnlySubgridV1
3 3 4 EmptySubgridV1
3 4 0 ImportOnlySubgridV1
3 4 1 ImportOnlySubgridV1
3 4 2 EmptySubgridV1
3 4 3 ImportOnlySubgridV1
3 4 4 EmptySubgridV1
3 5 0 ImportOnlySubgridV1
3 5 1 ImportOnlySubgridV1
3 5 2 EmptySubgridV1
3 5 3 ImportOnlySubgridV1
3 5 4 EmptySubgridV1
3 6 0 ImportOnlySubgridV1
3 6 1 ImportOnlySubgridV1
3 6 2 EmptySubgridV1
3 6 3 ImportOnlySubgridV1
3 6 4 EmptySubgridV1
3 7 0 ImportOnlySubgridV1
3 7 1 ImportOnlySubgridV1
3 7 2 EmptySubgridV1
3 7 3 ImportOnlySubgridV1
3 7 4 EmptySubgridV1
4 0 0 ImportOnlySubgridV1
4 0 1 EmptySubgridV1
4 0 2 ImportOnlySubgridV1
4 0 3 EmptySubgridV1
4 0 4 ImportOnlySubgridV1
4 1 0 ImportOnlySubgridV1
4 1 1 EmptySubgridV1
4 1 2 ImportOnlySubgridV1
4 1 3 EmptySubgridV1
4 1 4 ImportOnlySubgridV1
4 2 0 ImportOnlySubgridV1
4 2 1 EmptySubgridV1
4 2 2 ImportOnlySubgridV1
4 2 3 EmptySubgridV1
4 2 4 ImportOnlySubgridV1
4 3 0 ImportOnlySubgridV1
4 3 1 EmptySubgridV1
4 3 2 ImportOnlySubgridV1
4 3 3 EmptySubgridV1
4 3 4 ImportOnlySubgridV1
4 4 0 ImportOnlySubgridV1
4 4 1 EmptySubgridV1
4 4 2 ImportOnlySubgridV1
4 4 3 EmptySubgridV1
4 4 4 ImportOnlySubgridV1
4 5 0 ImportOnlySubgridV1
4 5 1 EmptySubgridV1
4 5 2 ImportOnlySubgridV1
4 5 3 EmptySubgridV1
4 5 4 ImportOnlySubgridV1
4 6 0 ImportOnlySubgridV1
4 6 1 EmptySubgridV1
4 6 2 ImportOnlySubgridV1
4 6 3 EmptySubgridV1
4 6 4 ImportOnlySubgridV1
4 7 0 ImportOnlySubgridV1
4 7 1 EmptySubgridV1
4 7 2 ImportOnlySubgridV1
4 7 3 EmptySubgridV1
4 7 4 ImportOnlySubgridV1
5 0 0 EmptySubgridV1
5 0 1 EmptySubgridV1
5 0 2 EmptySubgridV1
5 0 3 EmptySubgridV1
5 0 4 EmptySubgridV1
5 1 0 EmptySubgridV1
5 1 1 EmptySubgridV1
5 1 2 EmptySubgridV1
5 1 3 EmptySubgridV1
5 1 4 EmptySubgridV1
5 2 0 EmptySubgridV1
5 2 1 EmptySubgridV1
5 2 2 EmptySubgridV1
5 2 3 EmptySubgridV1
5 2 4 EmptySubgridV1
5 3 0 EmptySubgridV1
5 3 1 EmptySubgridV1
5 3 2 EmptySubgridV1
5 3 3 EmptySubgridV1
5 3 4 EmptySubgridV1
5 4 0 EmptySubgridV1
5 4 1 EmptySubgridV1
5 4 2 EmptySubgridV1
5 4 3 EmptySubgridV1
5 4 4 EmptySubgridV1
5 5 0 EmptySubgridV1
5 5 1 EmptySubgridV1
5 5 2 EmptySubgridV1
5 5 3 EmptySubgridV1
5 5 4 EmptySubgridV1
5 6 0 EmptySubgridV1
5 6 1 EmptySubgridV1
5 6 2 EmptySubgridV1
5 6 3 EmptySubgridV1
5 6 4 EmptySubgridV1
5 7 0 EmptySubgridV1
5 7 1 EmptySubgridV1
5 7 2 EmptySubgridV1
5 7 3 EmptySubgridV1
5 7 4 EmptySubgridV1
6 0 0 ImportOnlySubgridV1
6 0 1 EmptySubgridV1
6 0 2 ImportOnlySubgridV1
6 0 3 EmptySubgridV1
6 0 4 ImportOnlySubgridV1
6 1 0 ImportOnlySubgridV1
6 1 1 EmptySubgridV1
6 1 2 ImportOnlySubgridV1
6 1 3 EmptySubgridV1
6 1 4 ImportOnlySubgridV1
6 2 0 ImportOnlySubgridV1
6 2 1 EmptySubgridV1
6 2 2 ImportOnlySubgridV1
6 2 3 EmptySubgridV1
6 2 4 ImportOnlySubgridV1
6 3 0 ImportOnlySubgridV1
6 3 1 EmptySubgridV1
6 3 2 ImportOnlySubgridV1
6 3 3 EmptySubgridV1
6 3 4 ImportOnlySubgridV1
6 4 0 ImportOnlySubgridV1
6 4 1 EmptySubgridV1
6 4 2 ImportOnlySubgridV1
6 4 3 EmptySubgridV1
6 4 4 ImportOnlySubgridV1
6 5 0 ImportOnlySubgridV1
6 5 1 EmptySubgridV1
6 5 2 ImportOnlySubgridV1
6 5 3 EmptySubgridV1
6 5 4 ImportOnlySubgridV1
6 6 0 ImportOnlySubgridV1
6 6 1 EmptySubgridV1
6 6 2 ImportOnlySubgridV1
6 6 3 EmptySubgridV1
6 6 4 ImportOnlySubgridV1
6 7 0 ImportOnlySubgridV1
6 7 1 EmptySubgridV1
6 7 2 ImportOnlySubgridV1
6 7 3 EmptySubgridV1
6 7 4 ImportOnlySubgridV1
";

    const X1_STR: &str = "o b l                                                                                                                                                                                                                                                        x1
-+-+-+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
0 0 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 1 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 2 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 3 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 4 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 5 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 6 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 7 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 0 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 0 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 0 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 1 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 1 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 1 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 2 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 2 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 2 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 3 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 3 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 3 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 4 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 4 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 4 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 5 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 5 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 5 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 6 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 6 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 6 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 7 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 7 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 7 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 0 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 0 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 0 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 1 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 1 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 1 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 2 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 2 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 2 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 3 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 3 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 3 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 4 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 4 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 4 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 5 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 5 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 5 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 6 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 6 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 6 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 7 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 7 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 7 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 0 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 0 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 0 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 1 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 1 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 1 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 2 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 2 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 2 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 3 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 3 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 3 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 4 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 4 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 4 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 5 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 5 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 5 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 6 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 6 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 6 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 7 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 7 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 7 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 0 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 0 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 0 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 1 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 1 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 1 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 2 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 2 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 2 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 3 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 3 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 3 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 4 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 4 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 4 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 5 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 5 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 5 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 6 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 6 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 6 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 7 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 7 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 7 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
";

    const X2_STR: &str = "o b l                                                                                                                                                                                                                                                        x2
-+-+-+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
0 0 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 1 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 2 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 3 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 4 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 5 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 6 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
0 7 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 0 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 0 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 0 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 1 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 1 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 1 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 2 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 2 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 2 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 3 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 3 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 3 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 4 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 4 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 4 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 5 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 5 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 5 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 6 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 6 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 6 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 7 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 7 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
1 7 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 0 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 0 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 0 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 1 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 1 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 1 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 2 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 2 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 2 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 3 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 3 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 3 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 4 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 4 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 4 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 5 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 5 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 5 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 6 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 6 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 6 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 7 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 7 1 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
3 7 3 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 0 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 0 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 0 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 1 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 1 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 1 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 2 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 2 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 2 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 3 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 3 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 3 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 4 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 4 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 4 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 5 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 5 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 5 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 6 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 6 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 6 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 7 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 7 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
4 7 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 0 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 0 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 0 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 1 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 1 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 1 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 2 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 2 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 2 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 3 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 3 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 3 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 4 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 4 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 4 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 5 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 5 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 5 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 6 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 6 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 6 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 7 0 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 7 2 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
6 7 4 1.000e0, 9.309e-1, 8.628e-1, 7.956e-1, 7.296e-1, 6.648e-1, 6.015e-1, 5.398e-1, 4.799e-1, 4.222e-1, 3.669e-1, 3.144e-1, 2.651e-1, 2.195e-1, 1.780e-1, 1.411e-1, 1.091e-1, 8.228e-2, 6.048e-2, 4.341e-2, 3.052e-2, 2.109e-2, 1.438e-2, 9.699e-3, 6.496e-3, 4.329e-3, 2.874e-3, 1.903e-3, 1.259e-3, 8.314e-4, 5.488e-4, 3.621e-4, 2.388e-4, 1.575e-4, 1.038e-4, 6.844e-5, 4.511e-5, 2.974e-5, 1.960e-5, 1.292e-5, 8.517e-6, 5.614e-6, 3.700e-6, 2.439e-6, 1.608e-6, 1.060e-6, 6.984e-7, 4.604e-7, 3.034e-7, 2.000e-7
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["subgrids", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }

    #[test]
    fn muf() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["subgrids", "--muf", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(MUF_STR);
    }

    #[test]
    fn mur() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["subgrids", "--mur", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(MUR_STR);
    }

    #[test]
    fn stats() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["subgrids", "--stats", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(STATS_STR);
    }

    #[test]
    fn type_() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["subgrids", "--type", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(TYPE_STR);
    }

    #[test]
    fn type_show_empty() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&[
                "subgrids",
                "--show-empty",
                "--type",
                "data/LHCB_WP_7TEV.pineappl.lz4",
            ])
            .assert()
            .success()
            .stdout(TYPE_SHOW_EMPTY_STR);
    }

    #[test]
    fn x1() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["subgrids", "--x1", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(X1_STR);
    }

    #[test]
    fn x2() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["subgrids", "--x2", "data/LHCB_WP_7TEV.pineappl.lz4"])
            .assert()
            .success()
            .stdout(X2_STR);
    }
}
