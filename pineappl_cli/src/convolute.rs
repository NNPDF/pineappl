use super::helpers::{self, Subcommand};
use anyhow::Result;
use clap::{Parser, ValueHint};
use lhapdf::Pdf;
use prettytable::{cell, Row};
use std::ops::Range;
use std::path::PathBuf;

/// Convolutes a PineAPPL grid with a PDF set.
#[derive(Parser)]
#[clap(name = "convolute")]
pub struct Opts {
    /// Path of the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// LHAPDF id(s) or name of the PDF set(s).
    #[clap(min_values = 1, validator = helpers::validate_pdfset)]
    pdfsets: Vec<String>,
    /// Show absolute numbers of the scale variation.
    #[clap(long, short)]
    absolute: bool,
    /// Selects a subset of bins.
    #[clap(
        long,
        short,
        multiple_values = true,
        parse(try_from_str = helpers::try_parse_integer_range),
        use_delimiter = true
    )]
    bins: Vec<Range<usize>>,
    /// Show integrated numbers (without bin widths) instead of differential ones.
    #[clap(long, short)]
    integrated: bool,
    /// Select orders manually.
    #[clap(
        long,
        multiple_values = true,
        parse(try_from_str = helpers::parse_order),
        short,
        use_delimiter = true
    )]
    orders: Vec<(u32, u32)>,
    /// Set the number of scale variations.
    #[clap(default_value = "7", long, possible_values = ["1", "3", "7", "9"], short)]
    scales: usize,
}

impl Subcommand for Opts {
    fn run(&self) -> Result<()> {
        let grid = helpers::read_grid(&self.input)?;
        let pdf = self.pdfsets[0].parse().map_or_else(
            |_| Pdf::with_setname_and_member(&self.pdfsets[0], 0),
            Pdf::with_lhaid,
        );
        let bins: Vec<_> = self
            .bins
            .iter()
            .cloned()
            .flat_map(|range| range.collect::<Vec<_>>())
            .collect();

        let results = helpers::convolute(&grid, &pdf, &self.orders, &bins, &[], self.scales);

        let other_results: Vec<f64> = self.pdfsets[1..]
            .iter()
            .flat_map(|pdfset| {
                let pdf = pdfset
                    .parse()
                    .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);
                helpers::convolute(&grid, &pdf, &[], &bins, &[], 1)
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
        title.add_cell(cell!(c->if self.integrated { "integ" } else { y_label }));

        if self.absolute {
            for scale in &helpers::SCALES_VECTOR[0..self.scales] {
                title.add_cell(cell!(c->&format!("({},{})", scale.0, scale.1)));
            }
        } else {
            title.add_cell(cell!(c->"neg unc"));
            title.add_cell(cell!(c->"pos unc"));
        }

        for other in self.pdfsets[1..].iter() {
            let mut cell = cell!(c->other);
            cell.set_hspan(2);
            title.add_cell(cell);
        }

        let mut table = helpers::create_table();
        table.set_titles(title);

        for (index, values) in results.chunks_exact(self.scales).enumerate() {
            let min_value = values
                .iter()
                .min_by(|left, right| left.partial_cmp(right).unwrap())
                .unwrap();
            let max_value = values
                .iter()
                .max_by(|left, right| left.partial_cmp(right).unwrap())
                .unwrap();
            let bin = if bins.is_empty() { index } else { bins[index] };

            let row = table.add_empty_row();

            row.add_cell(cell!(r->&format!("{}", bin)));
            for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                row.add_cell(cell!(r->&format!("{}", left[bin])));
                row.add_cell(cell!(r->&format!("{}", right[bin])));
            }
            row.add_cell(cell!(r->&format!("{:.7e}", if self.integrated { values[0] * normalizations[bin] } else { values[0] })));

            if self.absolute {
                for &value in values.iter() {
                    row.add_cell(cell!(r->&format!("{:.7e}", if self.integrated { value * normalizations[bin] } else { value })));
                }
            } else {
                row.add_cell(cell!(r->&format!("{:.2}%", (min_value / values[0] - 1.0) * 100.0)));
                row.add_cell(cell!(r->&format!("{:.2}%", (max_value / values[0] - 1.0) * 100.0)));
            }

            let bins = if bins.is_empty() {
                bin_info.bins()
            } else {
                self.bins.len()
            };

            for &other in other_results.iter().skip(index).step_by(bins) {
                row.add_cell(cell!(r->&format!("{:.7e}", if self.integrated { other * normalizations[bin] } else { other })));
                row.add_cell(cell!(r->&format!("{:.2}%", (other / values[0] - 1.0) * 100.0)));
            }
        }

        table.printstd();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-convolute 

Convolutes a PineAPPL grid with a PDF set

USAGE:
    pineappl convolute [OPTIONS] <INPUT> [--] [PDFSETS]...

ARGS:
    <INPUT>         Path of the input grid
    <PDFSETS>...    LHAPDF id(s) or name of the PDF set(s)

OPTIONS:
    -a, --absolute              Show absolute numbers of the scale variation
    -b, --bins <BINS>...        Selects a subset of bins
    -h, --help                  Print help information
    -i, --integrated            Show integrated numbers (without bin widths) instead of differential
                                ones
    -o, --orders <ORDERS>...    Select orders manually
    -s, --scales <SCALES>       Set the number of scale variations [default: 7] [possible values: 1,
                                3, 7, 9]
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["convolute", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }
}
