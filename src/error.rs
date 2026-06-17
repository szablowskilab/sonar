/// Alias for sonar results.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors that can occur during library generation.
#[derive(Debug)]
pub enum Error {
    InvalidStopCount,
    InvalidFrame,
    InvalidSesLengthRange,
    InvalidStopWindowRange,
    InvalidRegion {
        start: usize,
        end: usize,
    },
    SesLengthOutsideRange {
        mean_length: usize,
        start: usize,
        end: usize,
    },
    SesLongerThanTarget {
        ses_length: usize,
        target_length: usize,
    },
}

impl std::error::Error for Error {}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidStopCount => write!(f, "stop-count must be >= 1"),
            Error::InvalidFrame => write!(f, "frame must be 0, 1, or 2"),
            Error::InvalidSesLengthRange => write!(f, "invalid ses-length range"),
            Error::InvalidStopWindowRange => write!(f, "invalid stop-window range"),
            Error::InvalidRegion { start, end } => {
                write!(f, "Invalid region: start={} end={}", start, end)
            }
            Error::SesLengthOutsideRange {
                mean_length,
                start,
                end,
            } => {
                write!(
                    f,
                    "SES length {} is outside the range [{}, {}]",
                    mean_length, start, end
                )
            }
            Error::SesLongerThanTarget {
                ses_length,
                target_length,
            } => {
                write!(
                    f,
                    "SES length {} is longer than target length {}",
                    ses_length, target_length
                )
            }
        }
    }
}
