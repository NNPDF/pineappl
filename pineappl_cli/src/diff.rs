use lhapdf::Pdf;
use pineappl::grid::Grid;
use prettytable::{cell, Row, Table};
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;

use super::helpers;

pub fn subcommand(
    input1: &str,
    input2: &str,
    pdfset: &str,
    ignore_orders: bool,
) -> Result<Table, Box<dyn Error>> {
    let grid1 = Grid::read(BufReader::new(File::open(input1)?))?;
    let grid2 = Grid::read(BufReader::new(File::open(input2)?))?;
    let pdf = pdfset
        .parse()
        .map_or_else(|_| Pdf::with_setname_and_member(pdfset, 0), Pdf::with_lhaid);

    let mut table = helpers::create_table();

    if grid1.bin_info() == grid2.bin_info() {
        if ignore_orders {
            let bin_info = grid1.bin_info();
            let left_limits: Vec<_> = (0..bin_info.dimensions())
                .map(|i| bin_info.left(i))
                .collect();
            let right_limits: Vec<_> = (0..bin_info.dimensions())
                .map(|i| bin_info.right(i))
                .collect();

            let mut title = Row::empty();
            title.add_cell(cell!(c->"bin"));
            for i in 0..bin_info.dimensions() {
                let mut cell = cell!(c->&format!("x{}", i + 1));
                cell.set_hspan(2);
                title.add_cell(cell);
            }

            table.set_titles(title);

            let results1 = helpers::convolute(&grid1, &pdf, &[], &[], &[], 1);
            let results2 = helpers::convolute(&grid2, &pdf, &[], &[], &[], 1);

            for (bin, (result1, result2)) in results1.iter().zip(results2.iter()).enumerate() {
                let row = table.add_empty_row();

                row.add_cell(cell!(r->bin));
                for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                    row.add_cell(cell!(r->&format!("{}", left[bin])));
                    row.add_cell(cell!(r->&format!("{}", right[bin])));
                }

                row.add_cell(cell!(r->&format!("{:.7e}", result1)));
                row.add_cell(cell!(r->&format!("{:.7e}", result2)));
                row.add_cell(cell!(r->&format!("{:.3e}",
                    if result1 == result2 { 0.0 } else { result1 / result2 - 1.0 })));
            }
        } else {
            let orders1: HashSet<_> = grid1
                .orders()
                .iter()
                .filter(|order| (order.logxir == 0) && (order.logxif == 0))
                .collect();
            let orders2: HashSet<_> = grid2
                .orders()
                .iter()
                .filter(|order| (order.logxir == 0) && (order.logxif == 0))
                .collect();

            let mut diff1: Vec<_> = orders1.difference(&orders2).collect();
            diff1.sort();
            let diff1 = diff1;
            let mut diff2: Vec<_> = orders2.difference(&orders1).collect();
            diff2.sort();
            let diff2 = diff2;

            if !diff1.is_empty() || !diff2.is_empty() {
                print!("--- Orders: ");
                for order in &diff1 {
                    if order.alphas == 0 {
                        print!("O(a^{}) ", order.alpha);
                    } else if order.alpha == 0 {
                        print!("O(as^{}) ", order.alphas);
                    } else {
                        print!("O(as^{} a^{}) ", order.alphas, order.alpha);
                    }
                }
                println!();
                print!("+++ Orders: ");
                for order in &diff2 {
                    if order.alphas == 0 {
                        print!("O(a^{}) ", order.alpha);
                    } else if order.alpha == 0 {
                        print!("O(as^{}) ", order.alphas);
                    } else {
                        print!("O(as^{} a^{}) ", order.alphas, order.alpha);
                    }
                }
                println!();
                println!();
            }

            let bin_info = grid1.bin_info();
            let left_limits: Vec<_> = (0..bin_info.dimensions())
                .map(|i| bin_info.left(i))
                .collect();
            let right_limits: Vec<_> = (0..bin_info.dimensions())
                .map(|i| bin_info.right(i))
                .collect();

            let mut title = Row::empty();
            title.add_cell(cell!(c->"bin"));
            for i in 0..bin_info.dimensions() {
                let mut cell = cell!(c->&format!("x{}", i + 1));
                cell.set_hspan(2);
                title.add_cell(cell);
            }

            let mut orders: Vec<_> = orders1.intersection(&orders2).collect();
            orders.sort();
            let orders = orders;

            for order in &orders {
                let mut cell = cell!(c->&format!("O(as^{} a^{})", order.alphas, order.alpha));
                cell.set_hspan(3);
                title.add_cell(cell);
            }

            table.set_titles(title);

            let order_results1: Vec<Vec<f64>> = orders
                .iter()
                .map(|order| {
                    let mut order_mask = vec![false; grid1.orders().len()];
                    order_mask[grid1.orders().iter().position(|o| o == **order).unwrap()] = true;
                    helpers::convolute(&grid1, &pdf, &order_mask, &[], &[], 1)
                })
                .collect();
            let order_results2: Vec<Vec<f64>> = orders
                .iter()
                .map(|order| {
                    let mut order_mask = vec![false; grid2.orders().len()];
                    order_mask[grid2.orders().iter().position(|o| o == **order).unwrap()] = true;
                    helpers::convolute(&grid2, &pdf, &order_mask, &[], &[], 1)
                })
                .collect();

            for bin in 0..bin_info.bins() {
                let row = table.add_empty_row();

                row.add_cell(cell!(r->bin));
                for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                    row.add_cell(cell!(r->&format!("{}", left[bin])));
                    row.add_cell(cell!(r->&format!("{}", right[bin])));
                }

                for (result1, result2) in order_results1.iter().zip(order_results2.iter()) {
                    let result1 = result1[bin];
                    let result2 = result2[bin];
                    row.add_cell(cell!(r->&format!("{:.7e}", result1)));
                    row.add_cell(cell!(r->&format!("{:.7e}", result2)));
                    row.add_cell(cell!(r->&format!("{:.3e}",
                        if result1 == result2 { 0.0 } else { result1 / result2 - 1.0 })));
                }
            }
        }
    } else {
        let bin_info = grid1.bin_info();
        let left_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.left(i))
            .collect();
        let right_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.right(i))
            .collect();

        println!("--- Bin limits:");
        for bin in 0..bin_info.bins() {
            for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                print!("{} ", left[bin]);
                print!("{} ", right[bin]);
            }
            println!();
        }

        let bin_info = grid2.bin_info();
        let left_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.left(i))
            .collect();
        let right_limits: Vec<_> = (0..bin_info.dimensions())
            .map(|i| bin_info.right(i))
            .collect();

        println!("+++ Bin limits:");
        for bin in 0..bin_info.bins() {
            for (left, right) in left_limits.iter().zip(right_limits.iter()) {
                print!("{} ", left[bin]);
                print!("{} ", right[bin]);
            }
            println!();
        }
    }

    Ok(table)
}
