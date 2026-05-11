//! Error types for the `pineappl` library.

use thiserror::Error;

/// Catch-all error for this crate.
#[derive(Debug, Error)]
pub enum Error {
    /// An error that originates in this crate.
    #[error("{0}")]
    General(String),
    /// Error that does not originate from this crate.
    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

/// Convenient result type using [`Error`].
pub type Result<T> = std::result::Result<T, Error>;
