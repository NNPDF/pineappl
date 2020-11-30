use pineappl::grid::{Grid, Order};
use pineappl::lumi_entry;
use pineappl::subgrid::SubgridParams;
use std::convert::TryFrom;
use std::env;
use std::fs;
use yaml_rust::YamlLoader;

fn main() {
    let args: Vec<_> = env::args().collect();

    if args.len() != 3 {
        println!("usage: {} <input-yaml> <output-pineappl>", args[0]);
        return;
    }

    let input = &args[1];
    //let output = &args[2];

    let yaml = YamlLoader::load_from_str(&fs::read_to_string(input).unwrap()).unwrap();

    assert_eq!(yaml.len(), 1);

    let doc = &yaml[0];

    let interpolation_is_log = doc["interpolation_is_log"].as_bool().unwrap();
    let interpolation_polynomial_degree =
        usize::try_from(doc["interpolation_polynomial_degree"].as_i64().unwrap()).unwrap();
    let interpolation_xgrid = doc["interpolation_xgrid"]
        .as_vec()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect::<Vec<_>>();
    let pids = doc["pids"].as_vec().unwrap();
    //let xif = doc["xiF"].as_f64().unwrap();

    assert_eq!(interpolation_is_log, false);
    //assert_eq!(xif, 1.0);

    let lepton_pid = 11;
    let lumi: Vec<_> = pids
        .iter()
        .map(|pid| {
            lumi_entry![
                i32::try_from(pid.as_i64().unwrap()).unwrap(),
                lepton_pid,
                1.0
            ]
        })
        .collect();
    // TODO: right now there's only a leading order
    let orders = vec![Order::new(0, 2, 0, 0)];
    let bins = doc["F2total"].as_vec().unwrap().iter().count();
    let bin_limits: Vec<_> = (0..=bins).map(|limit| limit as f64).collect();
    let mut params = SubgridParams::default();

    //params.set_q2_bins();
    //params.set_q2_max();
    //params.set_q2_min();
    //params.set_q2_order();
    params.set_reweight(false);
    params.set_x_bins(interpolation_xgrid.len());
    params.set_x_max(*interpolation_xgrid.last().unwrap());
    params.set_x_min(*interpolation_xgrid.first().unwrap());
    params.set_x_order(interpolation_polynomial_degree);

    // TODO: check that the x-grid points are the really the same generated from the f2 function

    let _ = Grid::new(lumi, orders, bin_limits, params);

    // TODO: loop over observables and create subgrids

    // TODO: write a BinRemapper to set the correct values of the observables

    // TODO: write the grid to `output`

    //grid.write();
}
