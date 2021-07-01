use super::helpers;
use anyhow::Result;
use itertools::Itertools;

pub fn subcommand_qcd_ew(input: &str, mode: &str) -> Result<()> {
    let grid = helpers::read_grid(input)?;

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

pub fn subcommand_get(input: &str, key: &str) -> Result<()> {
    let mut grid = helpers::read_grid(input)?;

    grid.upgrade();

    if let Some(key_values) = grid.key_values() {
        if let Some(value) = key_values.get(key) {
            println!("{}", value);
        }
    } else {
        unreachable!();
    }

    Ok(())
}

pub fn subcommand_keys(input: &str) -> Result<()> {
    let mut grid = helpers::read_grid(input)?;

    grid.upgrade();

    if let Some(key_values) = grid.key_values() {
        let mut vector = key_values.iter().collect::<Vec<_>>();
        vector.sort();

        for (key, _) in &vector {
            println!("{}", key);
        }
    } else {
        unreachable!();
    }

    Ok(())
}

pub fn subcommand_show(input: &str) -> Result<()> {
    let mut grid = helpers::read_grid(input)?;

    grid.upgrade();

    if let Some(key_values) = grid.key_values() {
        let mut vector = key_values.iter().collect::<Vec<_>>();
        vector.sort();

        for (key, value) in &vector {
            println!("{}: {}", key, value);
        }
    } else {
        unreachable!();
    }

    Ok(())
}
