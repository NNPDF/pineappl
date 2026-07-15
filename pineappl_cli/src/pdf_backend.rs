//! PDF backend abstraction layer.
//!
//! This module provides a unified interface for different PDF interpolation backends,
//! currently supporting `LHAPDF` and `NeoPDF`. It allows runtime selection of the backend
//! and provides type-safe access to PDF metadata.

use anyhow::{Context as _, Result, anyhow};
use enum_dispatch::enum_dispatch;
// use std::fmt;
use lhapdf::PdfSet;
use neopdf::gridpdf::ForcePositive as NeoPDFForcePositive;
use neopdf::pdf::PDF;
use std::str::FromStr;

pub use neopdf::uncertainty::CL_1_SIGMA;

/// Method for handling negative PDF values.
#[derive(Clone, Copy, Debug)]
pub enum ForcePositive {
    /// No clipping - return values as-is.
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
                "unknown PDF backend '{s}'; must be 'lhapdf' or 'neopdf'"
            )),
        }
    }
}

/// Unified PDF backend interface.
///
/// This trait provides a common interface for PDF interpolation backends,
/// allowing the CLI to work with different implementations transparently.
#[enum_dispatch]
pub trait PdfBackend: Send {
    /// Evaluates the strong coupling constant `alpha_s(Q^2)`.
    fn alphas_q2(&self, q2: f64) -> f64;

    // /// Returns the error type of the PDF set (e.g., "replicas", "hessian").
    // fn error_type(&self) -> String;

    // /// Returns whether this is a polarized PDF.
    // fn is_polarized(&self) -> bool;

    /// Returns the hadron particle ID (PDG code).
    fn particle_id(&self) -> i32;

    /// Sets the method for handling negative PDF values.
    fn set_force_positive(&mut self, method: ForcePositive);

    // /// Returns the PDF set type (space-like or time-like).
    // fn set_type(&self) -> SetType;

    /// Returns the maximum valid x value.
    fn x_max(&mut self) -> f64;

    /// Returns the minimum valid x value.
    fn x_min(&mut self) -> f64;

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
}

/// TODO.
#[expect(
    clippy::large_enum_variant,
    reason = "LhapdfPdf is artificially small due to indirection"
)]
#[enum_dispatch(PdfBackend)]
pub enum ConvFunBackend {
    /// TODO.
    LhapdfPdf,
    /// TODO.
    NeopdfPdf,
}

/// Unified PDF set interface for uncertainty calculations.
#[enum_dispatch]
pub trait PdfSetBackend {
    /// Creates all PDF members in the set.
    ///
    /// # Errors
    ///
    /// TODO.
    fn mk_pdfs(&self) -> Result<Vec<ConvFunBackend>>;

    /// Calculates the uncertainty for a set of values.
    ///
    /// # Errors
    ///
    /// TODO.
    fn uncertainty(&self, values: &[f64], cl: f64, alternative: bool) -> Result<Uncertainty>;
}

/// TODO.
#[enum_dispatch(PdfSetBackend)]
pub enum ConvFunSetBackend {
    /// TODO.
    LhapdfSet,
    /// TODO.
    NeopdfSet,
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

/// LHAPDF backend wrapper.
pub struct LhapdfPdf {
    pdf: lhapdf::Pdf,
}

impl LhapdfPdf {
    /// Creates a new LHAPDF PDF by LHAID.
    ///
    /// # Errors
    ///
    /// TODO.
    pub fn with_lhaid(lhaid: i32) -> Result<Self> {
        Ok(Self {
            pdf: lhapdf::Pdf::with_lhaid(lhaid)
                .with_context(|| format!("failed to load LHAPDF with LHAID {lhaid}"))?,
        })
    }

    /// Creates a new LHAPDF PDF by set name and member index.
    ///
    /// # Errors
    ///
    /// TODO.
    pub fn with_setname_and_member(setname: &str, member: i32) -> Result<Self> {
        Ok(Self {
            pdf: lhapdf::Pdf::with_setname_and_member(setname, member).with_context(|| {
                format!("failed to load LHAPDF set '{setname}' member {member}")
            })?,
        })
    }
}

impl PdfBackend for LhapdfPdf {
    fn alphas_q2(&self, q2: f64) -> f64 {
        self.pdf.alphas_q2(q2)
    }

    fn particle_id(&self) -> i32 {
        self.pdf
            .set()
            .entry("Particle")
            .map_or(Ok(2212), |s| s.parse())
            .unwrap_or(2212)
    }

    fn set_force_positive(&mut self, method: ForcePositive) {
        match method {
            ForcePositive::None => self.pdf.set_force_positive(0),
            ForcePositive::ClipNegative => self.pdf.set_force_positive(1),
        }
    }

    fn x_max(&mut self) -> f64 {
        self.pdf.x_max()
    }

    fn x_min(&mut self) -> f64 {
        self.pdf.x_min()
    }

    fn xfx_q2(&self, id: i32, x: f64, q2: f64) -> f64 {
        self.pdf.xfx_q2(id, x, q2)
    }
}

/// LHAPDF set wrapper.
pub struct LhapdfSet {
    set: lhapdf::PdfSet,
}

impl LhapdfSet {
    /// Creates a new LHAPDF set by name.
    ///
    /// # Errors
    ///
    /// TODO.
    pub fn new(setname: &str) -> Result<Self> {
        Ok(Self {
            set: PdfSet::new(setname)
                .with_context(|| format!("failed to load LHAPDF set '{setname}'"))?,
        })
    }
}

impl PdfSetBackend for LhapdfSet {
    fn mk_pdfs(&self) -> Result<Vec<ConvFunBackend>> {
        let pdfs = self.set.mk_pdfs()?;
        Ok(pdfs
            .into_iter()
            .map(|pdf| LhapdfPdf { pdf }.into())
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

/// `NeoPDF` backend wrapper.
pub struct NeopdfPdf {
    pdf: PDF,
}

impl NeopdfPdf {
    /// Loads a NeoPDF member by set name and member index.
    #[must_use]
    pub fn load(pdf_name: &str, member: usize) -> Self {
        Self {
            pdf: PDF::load(pdf_name, member),
        }
    }

    /// Loads a NeoPDF member by its LHAPDF ID (LHAID).
    #[must_use]
    pub fn load_by_lhaid(lhaid: u32) -> Self {
        Self {
            pdf: PDF::load_by_lhaid(lhaid),
        }
    }
}

impl PdfBackend for NeopdfPdf {
    fn alphas_q2(&self, q2: f64) -> f64 {
        self.pdf.alphas_q2(q2)
    }

    fn particle_id(&self) -> i32 {
        self.pdf.metadata().hadron_pid
    }

    fn set_force_positive(&mut self, method: ForcePositive) {
        match method {
            ForcePositive::None => {
                self.pdf.set_force_positive(NeoPDFForcePositive::NoClipping);
            }
            ForcePositive::ClipNegative => {
                self.pdf
                    .set_force_positive(NeoPDFForcePositive::ClipNegative);
            }
        }
    }

    fn x_max(&mut self) -> f64 {
        self.pdf.metadata().x_max
    }

    fn x_min(&mut self) -> f64 {
        self.pdf.metadata().x_min
    }

    fn xfx_q2(&self, id: i32, x: f64, q2: f64) -> f64 {
        self.pdf.xfxq2(id, &[x, q2])
    }
}

/// `NeoPDF` set wrapper.
pub struct NeopdfSet {
    pdf_name: String,
    // num_members: usize,
    error_type: String,
}

impl NeopdfSet {
    /// Creates a new NeoPDF set by name.
    #[must_use]
    pub fn new(pdf_name: &str) -> Self {
        // Load member 0 to get metadata
        let pdf0 = PDF::load(pdf_name, 0);
        let metadata = pdf0.metadata();

        Self {
            pdf_name: pdf_name.to_owned(),
            // num_members: metadata.num_members as usize,
            error_type: metadata.error_type.clone(),
        }
    }
}

impl PdfSetBackend for NeopdfSet {
    // fn num_members(&self) -> usize {
    //     self.num_members
    // }

    fn mk_pdfs(&self) -> Result<Vec<ConvFunBackend>> {
        let pdfs = PDF::load_pdfs(&self.pdf_name);
        Ok(pdfs
            .into_iter()
            .map(|pdf| NeopdfPdf { pdf }.into())
            .collect())
    }

    fn uncertainty(&self, values: &[f64], cl: f64, alternative: bool) -> Result<Uncertainty> {
        let unc =
            neopdf::uncertainty::uncertainty(values, &self.error_type, CL_1_SIGMA, cl, alternative)
                .map_err(|e| anyhow::anyhow!(e))?;

        Ok(Uncertainty {
            central: unc.central,
            errminus: unc.errminus,
            errplus: unc.errplus,
        })
    }
}

/// Creates a single PDF from a set name and member index.
///
/// # Errors
///
/// TODO.
pub fn create_pdf(name: &str, member: usize, backend: Backend) -> Result<ConvFunBackend> {
    match backend {
        Backend::Lhapdf => {
            if let Ok(lhaid) = name.parse::<i32>() {
                Ok(LhapdfPdf::with_lhaid(lhaid)?.into())
            } else {
                Ok(LhapdfPdf::with_setname_and_member(name, member.try_into().unwrap())?.into())
            }
        }
        Backend::Neopdf => {
            if let Ok(lhaid) = name.parse::<u32>() {
                Ok(NeopdfPdf::load_by_lhaid(lhaid).into())
            } else {
                Ok(NeopdfPdf::load(name, member).into())
            }
        }
    }
}

/// Creates a PDF set for uncertainty calculations.
///
/// # Errors
///
/// TODO.
pub fn create_pdf_set(name: &str, backend: Backend) -> Result<ConvFunSetBackend> {
    match backend {
        Backend::Lhapdf => {
            // Try parsing as LHAID first
            let setname = if let Ok(lhaid) = name.parse::<i32>() {
                lhapdf::lookup_pdf(lhaid)
                    .map(|(set, _)| set)
                    .ok_or_else(|| anyhow!("no convolution function for LHAID = `{lhaid}` found"))?
            } else {
                name.to_owned()
            };
            Ok(LhapdfSet::new(&setname)?.into())
        }
        Backend::Neopdf => Ok(NeopdfSet::new(name).into()),
    }
}
