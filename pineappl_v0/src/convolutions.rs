//! Module for everything related to luminosity functions.

use super::pids;

/// Data type that indentifies different types of convolutions.
#[derive(Debug, Eq, PartialEq)]
pub enum Convolution {
    // TODO: eventually get rid of this value
    /// No convolution.
    None,
    /// Unpolarized parton distribution function. The integer denotes the type of hadron with a PDG
    /// MC ID.
    UnpolPDF(i32),
    /// Polarized parton distribution function. The integer denotes the type of hadron with a PDG
    /// MC ID.
    PolPDF(i32),
    /// Unpolarized fragmentation function. The integer denotes the type of hadron with a PDG MC
    /// ID.
    UnpolFF(i32),
    /// Polarized fragmentation function. The integer denotes the type of hadron with a PDG MC ID.
    PolFF(i32),
}

impl Convolution {
    /// Return the convolution if the PID is charged conjugated.
    #[must_use]
    pub const fn charge_conjugate(&self) -> Self {
        match *self {
            Self::None => Self::None,
            Self::UnpolPDF(pid) => Self::UnpolPDF(pids::charge_conjugate_pdg_pid(pid)),
            Self::PolPDF(pid) => Self::PolPDF(pids::charge_conjugate_pdg_pid(pid)),
            Self::UnpolFF(pid) => Self::UnpolFF(pids::charge_conjugate_pdg_pid(pid)),
            Self::PolFF(pid) => Self::PolFF(pids::charge_conjugate_pdg_pid(pid)),
        }
    }

    /// Return the PID of the convolution if it has any.
    #[must_use]
    pub const fn pid(&self) -> Option<i32> {
        match *self {
            Self::None => None,
            Self::UnpolPDF(pid) | Self::PolPDF(pid) | Self::UnpolFF(pid) | Self::PolFF(pid) => {
                Some(pid)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn convolution_charge_conjugate() {
        assert_eq!(Convolution::None.charge_conjugate(), Convolution::None);
        assert_eq!(
            Convolution::UnpolPDF(2212).charge_conjugate(),
            Convolution::UnpolPDF(-2212)
        );
        assert_eq!(
            Convolution::PolPDF(2212).charge_conjugate(),
            Convolution::PolPDF(-2212)
        );
        assert_eq!(
            Convolution::UnpolFF(2212).charge_conjugate(),
            Convolution::UnpolFF(-2212)
        );
        assert_eq!(
            Convolution::PolFF(2212).charge_conjugate(),
            Convolution::PolFF(-2212)
        );
    }

    #[test]
    fn convolution_pid() {
        assert_eq!(Convolution::None.pid(), None);
        assert_eq!(Convolution::UnpolPDF(2212).pid(), Some(2212));
        assert_eq!(Convolution::PolPDF(2212).pid(), Some(2212));
        assert_eq!(Convolution::UnpolFF(2212).pid(), Some(2212));
        assert_eq!(Convolution::PolFF(2212).pid(), Some(2212));
    }
}
