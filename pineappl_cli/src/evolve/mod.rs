#[cfg(feature = "evolve")]
mod eko;

use super::helpers::{self, ConvFuns, ConvoluteMode};
use super::pdf_backend::PdfBackend;
use super::{GlobalConfiguration, Subcommand};
use anyhow::{Result, anyhow};
use clap::{Parser, ValueHint};
use pineappl::fk_table::FkTable;
use pineappl::grid::Grid;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

/// Evolve a grid with an evolution kernel operator to an FK table.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[arg(value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the evolution kernel operator(s).
    #[arg(value_name = "EKO1,...")]
    ekos: String,
    /// Path to the converted grid.
    #[arg(value_hint = ValueHint::FilePath)]
    output: PathBuf,
    /// LHAPDF ID(s) or name of the PDF(s)/FF(s).
    conv_funs: ConvFuns,
    /// Relative threshold between the table and the converted grid when comparison fails.
    #[arg(default_value = "1e-3", long)]
    accuracy: f64,
    /// Set the number of fractional digits shown for absolute numbers.
    #[arg(default_value_t = 7, long, value_name = "ABS")]
    digits_abs: usize,
    /// Set the number of fractional digits shown for relative numbers.
    #[arg(default_value_t = 7, long, value_name = "REL")]
    digits_rel: usize,
    /// Select which orders to evolve.
    #[arg(
        long,
        num_args = 1,
        short,
        value_delimiter = ',',
        value_parser = helpers::parse_order
    )]
    orders: Vec<(u8, u8)>,
    /// Rescale the renormalization scale with this factor.
    #[arg(default_value_t = 1.0, long)]
    xir: f64,
    /// Rescale the factorization scale with this factor.
    #[arg(default_value_t = 1.0, long)]
    xif: f64,
    /// Rescale the fragmentation scale with this factor.
    #[arg(default_value_t = 1.0, long)]
    xia: f64,
}

impl Subcommand for Opts {
    fn run(&self, cfg: &GlobalConfiguration) -> Result<ExitCode> {
        use prettytable::row;

        let grid = helpers::read_grid(&self.input)?;
        let mut conv_funs = helpers::create_conv_funs(&self.conv_funs, cfg.pdf_backend)?;
        let results = helpers::convolve_scales(
            &grid,
            &mut conv_funs,
            &self.conv_funs.conv_types,
            &self.orders,
            &[],
            &[],
            &[(self.xir, self.xif, self.xia)],
            ConvoluteMode::Normal,
            cfg,
        );

        let eko_paths: Vec<_> = self.ekos.split(',').map(Path::new).collect();
        let fk_table = evolve_grid(
            &grid,
            &eko_paths,
            &*conv_funs[cfg.use_alphas_from],
            &self.orders,
            self.xir,
            self.xif,
            self.xia,
        )?;

        let evolved_results = helpers::convolve_scales(
            fk_table.grid(),
            &mut conv_funs,
            &self.conv_funs.conv_types,
            &[],
            &[],
            &[],
            &[(1.0, 1.0, 1.0)],
            ConvoluteMode::Normal,
            cfg,
        );

        // if both grids don't have the same number of bins there's a bug in the program
        assert_eq!(results.len(), evolved_results.len());

        let mut table = helpers::create_table();
        table.set_titles(row![c => "b", "Grid", "FkTable", "rel. diff"]);

        let mut different = false;

        for (bin, (one, two)) in results
            .into_iter()
            .zip(evolved_results.into_iter())
            .enumerate()
        {
            #[expect(clippy::float_cmp, reason = "here we really need an exact comparison")]
            let rel_diff = if one == two { 0.0 } else { two / one - 1.0 };

            if rel_diff.abs() > self.accuracy {
                different = true;
            }

            table.add_row(row![
                bin.to_string(),
                r->format!("{:.*e}", self.digits_abs, one),
                r->format!("{:.*e}", self.digits_abs, two),
                r->format!("{:.*e}", self.digits_rel, rel_diff)
            ]);
        }

        table.printstd();

        if different {
            Err(anyhow!("grids are different"))
        } else {
            helpers::write_grid(&self.output, fk_table.grid())
        }
    }
}

#[cfg(feature = "evolve")]
fn evolve_grid(
    grid: &Grid,
    ekos: &[&Path],
    use_alphas_from: &dyn PdfBackend,
    orders: &[(u8, u8)],
    xir: f64,
    xif: f64,
    xia: f64,
) -> Result<FkTable> {
    use eko::EkoSlices;
    use pineappl::evolution::AlphasTable;

    let order_mask: Vec<_> = grid
        .orders()
        .iter()
        .map(|order| {
            orders.is_empty()
                || orders
                    .iter()
                    .any(|other| (order.alphas == other.0) && (order.alpha == other.1))
        })
        .collect();

    let mut eko_slices: Vec<_> = ekos
        .iter()
        .map(|eko| EkoSlices::new(eko))
        .collect::<Result<_, _>>()?;
    let eko_slices: Vec<_> = eko_slices.iter_mut().collect();
    let alphas_table = AlphasTable::from_grid(grid, xir, &|q2| use_alphas_from.alphas_q2(q2));

    Ok(grid.evolve(eko_slices, &order_mask, (xir, xif, xia), &alphas_table)?)
}

#[cfg(not(feature = "evolve"))]
fn evolve_grid(
    _: &Grid,
    _: &[&Path],
    _: &dyn PdfBackend,
    _: &[(u8, u8)],
    _: f64,
    _: f64,
    _: f64,
) -> Result<FkTable> {
    Err(anyhow!(
        "you need to install `pineappl` with feature `evolve`"
    ))
}
