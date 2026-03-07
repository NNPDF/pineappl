//! PDF backend abstraction layer.
//!
//! This module provides a unified interface for different PDF interpolation backends,
//! currently supporting `LHAPDF` and `NeoPDF`. It allows runtime selection of the backend
//! and provides type-safe access to PDF metadata.

use anyhow::{anyhow, Context, Result};
use std::fmt;
use std::str::FromStr;

/// Constant for 1-sigma confidence level (for compatibility with `lhapdf::CL_1_SIGMA`).
pub const CL_1_SIGMA: f64 = 68.268_949_213_708_58;

/// The type of PDF set (space-like for PDFs, time-like for FFs).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SetType {
    /// Space-like PDF (standard parton distribution function).
    #[default]
    SpaceLike,
    /// Time-like PDF (fragmentation function).
    TimeLike,
}

impl fmt::Display for SetType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SpaceLike => write!(f, "SpaceLike"),
            Self::TimeLike => write!(f, "TimeLike"),
        }
    }
}

/// Method for handling negative PDF values.
#[derive(Clone, Copy, Debug, Default)]
pub enum ForcePositive {
    /// No clipping - return values as-is.
    #[default]
    None,
    /// Clip negative values to zero.
    ClipNegative,
}

/// Available PDF backends.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Backend {
    /// `LHAPDF` backend (C++ library with Rust bindings).
    #[default]
    Lhapdf,
    /// `NeoPDF` backend (pure Rust implementation).
    Neopdf,
}

impl FromStr for Backend {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "lhapdf" => Ok(Self::Lhapdf),
            "neopdf" => Ok(Self::Neopdf),
            _ => Err(anyhow!(
                "unknown PDF backend '{}'; must be 'lhapdf' or 'neopdf'",
                s
            )),
        }
    }
}

impl fmt::Display for Backend {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Lhapdf => write!(f, "lhapdf"),
            Self::Neopdf => write!(f, "neopdf"),
        }
    }
}

/// Unified PDF backend interface.
///
/// This trait provides a common interface for PDF interpolation backends,
/// allowing the CLI to work with different implementations transparently.
pub trait PdfBackend: Send {
    /// Evaluates xf(x, Q^2) for a given flavor.
    ///
    /// # Arguments
    /// * `id` - The Monte Carlo PDG flavor ID.
    /// * `x` - The momentum fraction.
    /// * `q2` - The squared energy scale.
    ///
    /// # Returns
    /// The PDF value xf(x, Q^2).
    ///
    /// # TODO
    /// Extend to support `NeoPDF` multi-parameters interpolation.
    fn xfx_q2(&self, id: i32, x: f64, q2: f64) -> f64;

    /// Evaluates the strong coupling constant `alpha_s(Q^2)`.
    fn alphas_q2(&self, q2: f64) -> f64;

    /// Returns the minimum valid x value.
    fn x_min(&mut self) -> f64;

    /// Returns the maximum valid x value.
    fn x_max(&mut self) -> f64;

    /// Returns the hadron particle ID (PDG code).
    fn particle_id(&self) -> i32;

    /// Returns whether this is a polarized PDF.
    fn is_polarized(&self) -> bool;

    /// Returns the PDF set type (space-like or time-like).
    fn set_type(&self) -> SetType;

    /// Sets the method for handling negative PDF values.
    fn set_force_positive(&mut self, method: ForcePositive);

    /// Returns the error type of the PDF set (e.g., "replicas", "hessian").
    fn error_type(&self) -> String;
}

/// Unified PDF set interface for uncertainty calculations.
pub trait PdfSetBackend {
    /// Returns the number of members in this set.
    fn num_members(&self) -> usize;

    /// Creates all PDF members in the set.
    fn mk_pdfs(&self) -> Result<Vec<Box<dyn PdfBackend>>>;

    /// Calculates the uncertainty for a set of values.
    fn uncertainty(&self, values: &[f64], cl: f64, alternative: bool) -> Result<Uncertainty>;
}

/// Uncertainty information from a PDF set.
#[derive(Clone, Debug)]
pub struct Uncertainty {
    /// Central value.
    pub central: f64,
    /// Negative error (absolute value).
    pub errminus: f64,
    /// Positive error (absolute value).
    pub errplus: f64,
}

// ============================================================================
// LHAPDF Backend Implementation
// ============================================================================

/// LHAPDF backend wrapper.
pub struct LhapdfPdf {
    pdf: lhapdf::Pdf,
}

impl LhapdfPdf {
    /// Creates a new LHAPDF PDF by set name and member index.
    pub fn with_setname_and_member(setname: &str, member: i32) -> Result<Self> {
        Ok(Self {
            pdf: lhapdf::Pdf::with_setname_and_member(setname, member).with_context(|| {
                format!("failed to load LHAPDF set '{setname}' member {member}")
            })?,
        })
    }

    /// Creates a new LHAPDF PDF by LHAID.
    pub fn with_lhaid(lhaid: i32) -> Result<Self> {
        Ok(Self {
            pdf: lhapdf::Pdf::with_lhaid(lhaid)
                .with_context(|| format!("failed to load LHAPDF with LHAID {lhaid}"))?,
        })
    }
}

impl PdfBackend for LhapdfPdf {
    fn xfx_q2(&self, id: i32, x: f64, q2: f64) -> f64 {
        self.pdf.xfx_q2(id, x, q2)
    }

    fn alphas_q2(&self, q2: f64) -> f64 {
        self.pdf.alphas_q2(q2)
    }

    fn x_min(&mut self) -> f64 {
        self.pdf.x_min()
    }

    fn x_max(&mut self) -> f64 {
        self.pdf.x_max()
    }

    fn particle_id(&self) -> i32 {
        self.pdf
            .set()
            .entry("Particle")
            .map_or(Ok(2212), |s| s.parse())
            .unwrap_or(2212)
    }

    fn is_polarized(&self) -> bool {
        self.pdf
            .set()
            .entry("Polarized")
            .is_some_and(|s| s.eq_ignore_ascii_case("true"))
    }

    fn set_type(&self) -> SetType {
        self.pdf
            .set()
            .entry("SetType")
            .map_or(SetType::SpaceLike, |s| {
                if s.eq_ignore_ascii_case("timelike") {
                    SetType::TimeLike
                } else {
                    SetType::SpaceLike
                }
            })
    }

    fn set_force_positive(&mut self, method: ForcePositive) {
        match method {
            ForcePositive::None => self.pdf.set_force_positive(0),
            ForcePositive::ClipNegative => self.pdf.set_force_positive(1),
        }
    }

    fn error_type(&self) -> String {
        self.pdf
            .set()
            .entry("ErrorType")
            .map(|s| s.to_owned())
            .unwrap_or_default()
    }
}

/// LHAPDF set wrapper.
pub struct LhapdfSet {
    set: lhapdf::PdfSet,
}

impl LhapdfSet {
    /// Creates a new LHAPDF set by name.
    pub fn new(setname: &str) -> Result<Self> {
        Ok(Self {
            set: lhapdf::PdfSet::new(setname)
                .with_context(|| format!("failed to load LHAPDF set '{setname}'"))?,
        })
    }
}

impl PdfSetBackend for LhapdfSet {
    fn num_members(&self) -> usize {
        self.set
            .entry("NumMembers")
            .and_then(|s| s.parse().ok())
            .unwrap_or(1)
    }

    fn mk_pdfs(&self) -> Result<Vec<Box<dyn PdfBackend>>> {
        let pdfs = self.set.mk_pdfs()?;
        Ok(pdfs
            .into_iter()
            .map(|pdf| Box::new(LhapdfPdf { pdf }) as Box<dyn PdfBackend>)
            .collect())
    }

    fn uncertainty(&self, values: &[f64], cl: f64, alternative: bool) -> Result<Uncertainty> {
        let unc = self.set.uncertainty(values, cl, alternative)?;
        Ok(Uncertainty {
            central: unc.central,
            errminus: unc.errminus,
            errplus: unc.errplus,
        })
    }
}

// ============================================================================
// NeoPDF Backend Implementation
// ============================================================================

/// `NeoPDF` backend wrapper.
pub struct NeopdfPdf {
    pdf: neopdf::pdf::PDF,
}

impl NeopdfPdf {
    /// Loads a NeoPDF member by set name and member index.
    pub fn load(pdf_name: &str, member: usize) -> Self {
        Self {
            pdf: neopdf::pdf::PDF::load(pdf_name, member),
        }
    }
}

impl PdfBackend for NeopdfPdf {
    fn xfx_q2(&self, id: i32, x: f64, q2: f64) -> f64 {
        self.pdf.xfxq2(id, &[x, q2])
    }

    fn alphas_q2(&self, q2: f64) -> f64 {
        self.pdf.alphas_q2(q2)
    }

    fn x_min(&mut self) -> f64 {
        self.pdf.metadata().x_min
    }

    fn x_max(&mut self) -> f64 {
        self.pdf.metadata().x_max
    }

    fn particle_id(&self) -> i32 {
        self.pdf.metadata().hadron_pid
    }

    fn is_polarized(&self) -> bool {
        self.pdf.metadata().polarised
    }

    fn set_type(&self) -> SetType {
        match self.pdf.metadata().set_type {
            neopdf::metadata::SetType::SpaceLike => SetType::SpaceLike,
            neopdf::metadata::SetType::TimeLike => SetType::TimeLike,
        }
    }

    fn set_force_positive(&mut self, method: ForcePositive) {
        match method {
            ForcePositive::None => {
                self.pdf
                    .set_force_positive(neopdf::gridpdf::ForcePositive::NoClipping);
            }
            ForcePositive::ClipNegative => {
                self.pdf
                    .set_force_positive(neopdf::gridpdf::ForcePositive::ClipNegative);
            }
        }
    }

    fn error_type(&self) -> String {
        self.pdf.metadata().error_type.clone()
    }
}

/// `NeoPDF` set wrapper.
pub struct NeopdfSet {
    pdf_name: String,
    num_members: usize,
    error_type: String,
}

impl NeopdfSet {
    /// Creates a new NeoPDF set by name.
    pub fn new(pdf_name: &str) -> Self {
        // Load member 0 to get metadata
        let pdf0 = neopdf::pdf::PDF::load(pdf_name, 0);
        let metadata = pdf0.metadata();

        Self {
            pdf_name: pdf_name.to_owned(),
            num_members: metadata.num_members as usize,
            error_type: metadata.error_type.clone(),
        }
    }
}

impl PdfSetBackend for NeopdfSet {
    fn num_members(&self) -> usize {
        self.num_members
    }

    fn mk_pdfs(&self) -> Result<Vec<Box<dyn PdfBackend>>> {
        let pdfs = neopdf::pdf::PDF::load_pdfs(&self.pdf_name);
        Ok(pdfs
            .into_iter()
            .map(|pdf| Box::new(NeopdfPdf { pdf }) as Box<dyn PdfBackend>)
            .collect())
    }

    fn uncertainty(&self, values: &[f64], cl: f64, _alternative: bool) -> Result<Uncertainty> {
        let central = values[0];

        let (errminus, errplus) = if self.error_type.to_lowercase().contains("replicas")
            || self.error_type.to_lowercase().contains("monte carlo")
        {
            let n = values.len() as f64;
            let mean: f64 = values.iter().skip(1).sum::<f64>() / (n - 1.0);
            let variance: f64 = values
                .iter()
                .skip(1)
                .map(|&v| (v - mean).powi(2))
                .sum::<f64>()
                / (n - 1.0);
            let std_dev = variance.sqrt();

            let scale = if (cl - 68.27).abs() < 0.1 {
                1.0
            } else if (cl - 95.45).abs() < 0.1 {
                2.0
            } else {
                cl / 100.0 * 1.645 // approximate
            };

            (std_dev * scale, std_dev * scale)
        } else {
            let mut sum_sq = 0.0;
            for pair in values[1..].chunks(2) {
                if pair.len() == 2 {
                    let diff = (pair[0] - pair[1]).abs() / 2.0;
                    sum_sq += diff.powi(2);
                }
            }
            let err = sum_sq.sqrt();
            (err, err)
        };

        Ok(Uncertainty {
            central,
            errminus,
            errplus,
        })
    }
}

// ============================================================================
// Factory Functions
// ============================================================================

/// Creates a single PDF from a set name and member index.
pub fn create_pdf(name: &str, member: usize, backend: Backend) -> Result<Box<dyn PdfBackend>> {
    match backend {
        Backend::Lhapdf => {
            if let Ok(lhaid) = name.parse::<i32>() {
                Ok(Box::new(LhapdfPdf::with_lhaid(lhaid)?))
            } else {
                Ok(Box::new(LhapdfPdf::with_setname_and_member(
                    name,
                    member.try_into().unwrap(),
                )?))
            }
        }
        // NOTE: `NeoPDF` doesn't support reading LHAID yet.
        Backend::Neopdf => Ok(Box::new(NeopdfPdf::load(name, member))),
    }
}

/// Creates a PDF set for uncertainty calculations.
pub fn create_pdf_set(name: &str, backend: Backend) -> Result<Box<dyn PdfSetBackend>> {
    match backend {
        Backend::Lhapdf => {
            // Try parsing as LHAID first
            let setname = if let Ok(lhaid) = name.parse::<i32>() {
                lhapdf::lookup_pdf(lhaid)
                    .map(|(set, _)| set)
                    .ok_or_else(|| anyhow!("no PDF set found for LHAID = {lhaid}"))?
            } else {
                name.to_owned()
            };
            Ok(Box::new(LhapdfSet::new(&setname)?))
        }
        Backend::Neopdf => Ok(Box::new(NeopdfSet::new(name))),
    }
}
