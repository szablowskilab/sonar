use super::design;
use derive_more::From;

/// Alias for sonar CLI results.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors that can occur during CLI execution.
#[derive(Debug, From)]
pub enum Error {
    #[from]
    ParseRange(design::ParseRangeError),
    #[from]
    ParseFasta(needletail::errors::ParseError),
    #[from]
    Utf8(std::str::Utf8Error),
    #[from]
    Io(std::io::Error),
    #[from]
    Sonar(sonar::prelude::Error),
    NoReferencePath,
    #[from]
    AdapterNotFound(design::AdapterNotFoundError),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::ParseRange(e) => write!(f, "{}", e),
            Error::ParseFasta(e) => write!(f, "{}", e),
            Error::Utf8(e) => write!(f, "{}", e),
            Error::Io(e) => write!(f, "{}", e),
            Error::Sonar(e) => write!(f, "{}", e),
            Error::NoReferencePath => write!(f, "no reference path provided"),
            Error::AdapterNotFound(e) => write!(f, "{}", e),
        }
    }
}
impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::ParseRange(e) => Some(e),
            Error::ParseFasta(e) => Some(e),
            Error::Utf8(e) => Some(e),
            Error::Io(e) => Some(e),
            Error::Sonar(e) => Some(e),
            Error::NoReferencePath => None,
            Error::AdapterNotFound(e) => Some(e),
        }
    }
}
