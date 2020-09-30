use itertools::Itertools;
use pineappl::grid::Grid;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;

pub fn subcommand(input: &str, mode: &str) -> Result<(), Box<dyn Error>> {
    let grid = Grid::read(BufReader::new(File::open(input)?))?;

    let mut sorted_grid_orders: Vec<_> = grid
        .orders()
        .iter()
        .filter(|order| (order.logxir == 0) && (order.logxif == 0))
        .collect();
    sorted_grid_orders.sort();

    let orders = sorted_grid_orders
        .into_iter()
        .group_by(|order| order.alphas + order.alpha)
        .into_iter()
        .map(|mut iter| match mode {
            "qcd" => iter.1.next().unwrap(),
            "ew" => iter.1.last().unwrap(),
            _ => unreachable!(),
        })
        .map(|order| {
            if order.alphas == 0 {
                format!("a{}", order.alpha)
            } else if order.alpha == 0 {
                format!("as{}", order.alphas)
            } else {
                format!("as{}a{}", order.alphas, order.alpha)
            }
        })
        .collect::<Vec<_>>()
        .join(",");

    println!("{}", orders);

    Ok(())
}
