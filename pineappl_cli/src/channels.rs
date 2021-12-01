use super::helpers;
use anyhow::Result;
use clap::Parser;
use lhapdf::Pdf;
use prettytable::{cell, Row, Table};
use std::ops::Range;

/// Shows the contribution for each partonic channel.
#[derive(Parser)]
#[clap(name = "channels")]
pub struct Opts {
    /// Path to the input grid.
    input: String,
    /// LHAPDF id or name of the PDF set.
    #[clap(validator = helpers::validate_pdfset)]
    pdfset: String,
    /// Show absolute numbers of each contribution.
    #[clap(long, short)]
    absolute: bool,
    /// The maximum number of channels displayed.
    #[clap(
        default_value = "10",
        long,
        short,
        validator = helpers::validate_pos_non_zero::<usize>
    )]
    limit: usize,
    /// Show integrated numbers (without bin widths) instead of differential ones.
    #[clap(long, requires = "absolute", short)]
    integrated: bool,
    /// Show only the listed channels.
    #[clap(
        conflicts_with = "limit",
        long,
        multiple_values = true,
        parse(try_from_str = helpers::try_parse_integer_range),
        use_delimiter = true
    )]
    lumis: Vec<Range<usize>>,
    /// Select orders manually.
    #[clap(
        long,
        multiple_values = true,
        short,
        parse(try_from_str = helpers::parse_order),
    )]
    orders: Vec<(u32, u32)>,
}

impl Opts {
    pub fn subcommand(&self) -> Result<Table> {
        let grid = helpers::read_grid(&self.input)?;
        let pdf = self.pdfset.parse().map_or_else(
            |_| Pdf::with_setname_and_member(&self.pdfset, 0),
            Pdf::with_lhaid,
        );
        let lumis: Vec<_> = self
            .lumis
            .iter()
            .cloned()
            .flat_map(|range| range.collect::<Vec<_>>())
            .collect();
        let limit = if lumis.is_empty() {
            grid.lumi().len().min(self.limit)
        } else {
            lumis.iter().filter(|lumi| **lumi < lumis.len()).count()
        };

        let results: Vec<_> = (0..grid.lumi().len())
            .map(|lumi| {
                let mut lumi_mask = vec![false; grid.lumi().len()];
                lumi_mask[lumi] = true;
                helpers::convolute(&grid, &pdf, &self.orders, &[], &lumi_mask, 1)
            })
            .collect();

        let bin_info = grid.bin_info();
        let left_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.left(i))
            .collect();
        let right_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.right(i))
            .collect();
        let normalizations = bin_info.normalizations();

        let labels = helpers::labels(&grid);
        let (y_label, x_labels) = labels.split_last().unwrap();
        let mut title = Row::empty();
        title.add_cell(cell!(c->"bin"));
        for x_label in x_labels {
            let mut cell = cell!(c->&x_label);
            cell.set_hspan(2);
            title.add_cell(cell);
        }
        for _ in 0..limit {
            title.add_cell(cell!(c->"lumi"));
            title.add_cell(
            cell!(c->if self.absolute { if self.integrated { "integ" } else { y_label } } else { "size" }),
        );
        }

        let mut table = helpers::create_table();
        table.set_titles(title);

        for bin in 0..bin_info.bins() {
            let row = table.add_empty_row();

            row.add_cell(cell!(r->&format!("{}", bin)));

            for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                row.add_cell(cell!(r->&format!("{}", left[bin])));
                row.add_cell(cell!(r->&format!("{}", right[bin])));
            }

            if self.absolute {
                let mut values: Vec<_> = results
                    .iter()
                    .enumerate()
                    .map(|(lumi, vec)| {
                        (
                            lumi,
                            if self.integrated {
                                normalizations[bin] * vec[bin]
                            } else {
                                vec[bin]
                            },
                        )
                    })
                    .collect();

                // sort using the absolute value in descending order
                values.sort_unstable_by(|(_, left), (_, right)| {
                    right.abs().partial_cmp(&left.abs()).unwrap()
                });

                for (lumi, value) in values
                    .iter()
                    .filter(|(lumi, _)| lumis.is_empty() || lumis.iter().any(|l| l == lumi))
                    .take(limit)
                {
                    row.add_cell(cell!(r->&format!("#{}", lumi)));
                    row.add_cell(cell!(r->&format!("{:.7e}", value)));
                }
            } else {
                let sum: f64 = results.iter().map(|vec| vec[bin]).sum();
                let mut percentages: Vec<_> = results
                    .iter()
                    .enumerate()
                    .map(|(lumi, vec)| (lumi, vec[bin] / sum * 100.0))
                    .collect();

                // sort using the absolute value in descending order
                percentages.sort_unstable_by(|(_, left), (_, right)| {
                    right.abs().partial_cmp(&left.abs()).unwrap()
                });

                for (lumi, percentage) in percentages
                    .iter()
                    .filter(|(lumi, _)| lumis.is_empty() || lumis.iter().any(|l| l == lumi))
                    .take(limit)
                {
                    row.add_cell(cell!(r->&format!("#{}", lumi)));
                    row.add_cell(cell!(r->&format!("{:.2}%", percentage)));
                }
            }
        }

        Ok(table)
    }
}
