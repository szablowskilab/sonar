use std::{fs, ops::RangeInclusive, path::PathBuf};

use clap::Parser;
use needletail::{parse_fastx_file, sequence::normalize};

use crate::{error::Result, table};
use sonar::{Candidate, DesignParams, generate_candidates};

pub fn generate(mut args: Args) -> Result<()> {
    let mut fasta_reader = parse_fastx_file(&args.target)?;
    let table_path = args.table.take();
    let fasta_path = args.fasta.take();
    let design_rules: DesignParams = args.into();
    let mut candidates: Vec<Candidate> = Vec::new();

    while let Some(Ok(record)) = fasta_reader.next() {
        let id = str::from_utf8(record.id())?;
        let seq =
            normalize(record.raw_seq(), design_rules.allow_iupac).unwrap_or(record.seq().to_vec());

        let candidate_sublibrary = generate_candidates(id, &seq, &design_rules)?;
        candidates.extend(candidate_sublibrary);
    }

    let table = table::render_tsv(&candidates);

    if let Some(path) = table_path {
        fs::write(path, table)?;
    } else {
        print!("{table}");
    }

    if let Some(path) = fasta_path {
        fs::write(path, render_fasta(&candidates))?;
    }

    Ok(())
}

fn render_fasta(candidates: &[Candidate]) -> String {
    let mut out = String::new();
    for candidate in candidates {
        out.push('>');
        out.push_str(&candidate.id);
        out.push('\n');
        out.push_str(&candidate.final_sesrna);
        out.push('\n');
    }
    out
}

#[derive(Debug, Parser)]
pub struct Args {
    /// Input target FASTA file
    pub target: PathBuf,

    /// 1-based inclusive target subregion, e.g. 10:350
    #[clap(long, short = 'r', value_parser = parse_range)]
    pub region: Option<RangeInclusive<usize>>,

    /// sesRNA length range in nucleotides as min:max
    #[clap(long, short = 'l', value_parser = parse_range, default_value = "200:300")]
    pub ses_length: RangeInclusive<usize>,

    /// Number of designed TAG stop codons to include
    #[clap(long, short = 's', default_value_t = 1)]
    pub stop_count: usize,

    /// allowed design stop position window as min:max from the sesRNA 5' end
    #[clap(long, short = 'w', value_parser = parse_range, default_value = "80:220")]
    pub stop_window: RangeInclusive<usize>,

    /// Minimum distance between designed TAG and any converted stop codons
    #[clap(long, short = 'd', default_value_t = 10)]
    pub min_stop_distance: usize,

    /// Translation frame for stop and ATG checks
    #[clap(long, short = 'o', default_value_t = 0)]
    pub frame: usize,

    /// Max number of candidates to output
    #[clap(long, short = 'm', default_value_t = 100)]
    pub max_candidates: usize,

    /// Do not spread out candidate designs across the target sequence
    #[clap(long, short = 'x', default_value_t = false)]
    pub no_spread: bool,

    /// Allow IUPAC nucleotide codes in the target sequence
    #[clap(long, short = 'n', default_value_t = false)]
    pub allow_iupac: bool,

    /// Write table to this path instead of stdout
    #[clap(long, short = 't')]
    pub table: Option<PathBuf>,

    /// Write FASTA output to this path instead of stdout
    #[clap(long, short = 'f')]
    pub fasta: Option<PathBuf>,
}

impl From<Args> for DesignParams {
    fn from(args: Args) -> Self {
        Self {
            region: args.region,
            ses_length: args.ses_length,
            stop_count: args.stop_count,
            stop_window: args.stop_window,
            min_stop_distance: args.min_stop_distance,
            frame: args.frame,
            max_candidates: args.max_candidates,
            no_spread: args.no_spread,
            allow_iupac: args.allow_iupac,
        }
    }
}

/// Parse a range string in the format `x:y` where `x` and `y` are 1-based inclusive range endpoints.
/// This value parser will return an error if the format is not `x:y`, if `x` is greater than `y`,
/// or if `x` or `y` are not valid positive integers.
pub fn parse_range(value: &str) -> std::result::Result<RangeInclusive<usize>, String> {
    let parts: Vec<_> = value.split(':').collect();
    if parts.len() != 2 {
        return Err(ParseRangeError::Format {
            recieved: value.to_string(),
        }
        .to_string());
    }

    let min: usize = parts[0]
        .parse()
        .map_err(|e: std::num::ParseIntError| e.to_string())?;
    let max: usize = parts[1]
        .parse()
        .map_err(|e: std::num::ParseIntError| e.to_string())?;
    if min == 0 {
        return Err(ParseRangeError::Min {
            value: min.to_string(),
        }
        .to_string());
    }
    if max < min {
        return Err(ParseRangeError::Max {
            value: max.to_string(),
        }
        .to_string());
    }

    Ok(min..=max)
}

/// Additional errors that can occur when parsing a range string.
#[derive(Debug)]
pub enum ParseRangeError {
    Format { recieved: String },
    Min { value: String },
    Max { value: String },
}

impl std::error::Error for ParseRangeError {}

impl std::fmt::Display for ParseRangeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseRangeError::Format { recieved } => {
                write!(f, "invalid range format: Expected: x:y, got: {}", recieved)
            }
            ParseRangeError::Min { value } => write!(f, "invalid range min: {}", value),
            ParseRangeError::Max { value } => write!(f, "invalid range max: {}", value),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_valid_range() {
        let result = parse_range("1:10");
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 1..=10);
    }

    #[test]
    fn reject_invalid_range() {
        let result = parse_range("1:x");
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "invalid digit found in string".to_string()
        );

        let result = parse_range("invalid");
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            ParseRangeError::Format {
                recieved: "invalid".to_string()
            }
            .to_string()
        );
    }

    #[test]
    fn args_map_no_spread_to_design_params() {
        let args = Args {
            target: PathBuf::from("target.fa"),
            region: None,
            ses_length: 200..=300,
            stop_count: 1,
            stop_window: 80..=220,
            min_stop_distance: 10,
            frame: 0,
            max_candidates: 100,
            no_spread: true,
            allow_iupac: false,
            table: None,
            fasta: None,
        };

        let design_params = DesignParams::from(args);

        assert!(design_params.no_spread);
    }
}
