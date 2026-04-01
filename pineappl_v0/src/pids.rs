//! TODO

/// Particle ID bases. In `PineAPPL` every particle is identified using a particle identifier
/// (PID), which is represented as an `i32`. The values of this `enum` specify how this value is
/// interpreted.
pub enum PidBasis {
    /// This basis uses the [particle data group](https://pdg.lbl.gov/) (PDG) PIDs. For a complete
    /// definition see the section 'Monte Carlo Particle Numbering Scheme' of the PDG Review, for
    /// instance the [2023 review](https://pdg.lbl.gov/2023/mcdata/mc_particle_id_contents.html).
    Pdg,
    /// This basis specifies the evolution basis, which is the same as [`PidBasis::Pdg`], except
    /// the following values have a special meaning: `100`, `103`, `108`, `115`, `124`, `135`,
    /// `200`, `203`, `208`, `215`, `224`, `235`.
    Evol,
}
