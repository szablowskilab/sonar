use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};

use clap::{Parser, builder};

use crate::error::Result;

pub fn map_to_ref(args: Args) -> Result<()> {
    let mut mapping_writer: Box<dyn Write> = if let Some(output) = &args.output {
        let mapping_file = File::create(output)?;
        Box::new(BufWriter::new(mapping_file))
    } else {
        Box::new(BufWriter::new(std::io::stdout()))
    };

    let _sensor_mappings = sonar::map_with_writer(
        args.sensor_path,
        args.reference_path,
        args.save_index,
        args.best_n,
        args.allow_iupac,
        &mut mapping_writer,
        None,
    )?;

    mapping_writer.flush()?;

    Ok(())
}

#[derive(Debug, Parser)]
pub struct Args {
    /// Input sensor candidates FASTA file.
    pub sensor_path: PathBuf,

    /// Reference transcriptome FASTA or index file.
    pub reference_path: PathBuf,

    /// Number of high scoring secondary alignments to consider.
    #[clap(long, short = 'm', default_value_t = 100, value_parser = builder::RangedU64ValueParser::<i32>::new().range(1..))]
    pub best_n: i32,

    /// Allow IUPAC characters in sensor sequences
    #[clap(long, short = 'n', default_value_t = false)]
    pub allow_iupac: bool,

    /// Optionally save reference index to the provided path
    #[clap(long, short = 'i')]
    pub save_index: Option<PathBuf>,

    /// Output mapping file
    #[clap(long, short = 'o')]
    pub output: Option<PathBuf>,
}
