//! Module containing all traits and supporting structures for grids.

/// Structure for the metadata of a subgrid.
pub struct Subgrid {
    /// Exponent of the strong coupling.
    pub alphas: u32,
    /// Exponent of the electromagnetic coupling.
    pub alpha: u32,
    /// Exponent of the logarithm of the scale factor of the renomalization scale.
    pub logxir: u32,
    /// Exponent of the logarithm of the scale factor of the factorization scale.
    pub logxif: u32,
}

/// Trait each grid must implement.
pub trait Grid {
    /// Fills the grid with events for the parton momentum fractions `x1` and `x2`, the scale `q2`,
    /// the observable `obs`. The events are stored in `weights` and must be ordered as the
    /// corresponding luminosity function was created. The weights are filled into the subgrid with
    /// the given index.
    fn fill(&mut self, x1: f64, x2: f64, q2: f64, obs: f64, weights: &[f64], subgrid: usize);
}

#[cfg(test)]
mod test {
    //use super::*;

    #[test]
    fn test() {}
}
