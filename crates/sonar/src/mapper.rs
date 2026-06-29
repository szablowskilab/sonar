//! Alignment mapping for sensor specificity checks.

use std::{io::Write, num::NonZeroI32, path::Path};

use log::warn;
use minimap2::{Aligner, Mapping};
use needletail::{Sequence as _, parse_fastx_file};

use crate::error::{Error, Result};

/// Formatter for accepted mappings.
pub type MappingFormatter<'a> = dyn Fn(&Mapping, &mut dyn Write) -> std::io::Result<()> + 'a;

/// Map sensor candidates to a reference.
///
/// ## Errors
///
/// This function will return an error if the aligner cannot be initialized.
///
/// If the alignment fails, a warning is logged and the sensor candidate is not
/// included in the output.
pub fn map<P: AsRef<Path>>(
    sensor_path: P,
    ref_path: P,
    save_index: Option<P>,
    best_n: i32,
    allow_iupac: bool,
) -> Result<Vec<Mapping>> {
    map_inner(sensor_path, ref_path, save_index, best_n, allow_iupac, None)
}

/// Map sensor candidates to a reference and write accepted mappings.
///
/// The writer receives one formatted row per mapping after the `max_hits`
/// filter has accepted the sensor candidate. If `format` is `None`, mappings
/// are written with [`mapping_to_paf`].
pub fn map_with_writer<P>(
    sensor_path: P,
    ref_path: P,
    save_index: Option<P>,
    best_n: i32,
    allow_iupac: bool,
    writer: &mut dyn Write,
    format: Option<&MappingFormatter<'_>>,
) -> Result<Vec<Mapping>>
where
    P: AsRef<Path>,
{
    let format = format.unwrap_or(&mapping_to_paf);

    map_inner(
        sensor_path,
        ref_path,
        save_index,
        best_n,
        allow_iupac,
        Some((writer, format)),
    )
}

fn map_inner<P: AsRef<Path>>(
    sensor_path: P,
    ref_path: P,
    save_index: Option<P>,
    best_n: i32,
    allow_iupac: bool,
    mut output: Option<(&mut dyn Write, &MappingFormatter<'_>)>,
) -> Result<Vec<Mapping>> {
    let mut sensors = parse_fastx_file(sensor_path)?;
    let mut sensor_mappings: Vec<Mapping> = Vec::new();

    let mut aligner = Aligner::builder().sr(); // short read alignment without splicing
    aligner.mapopt.unset_no_print_2nd(); // allow for secondary alignments
    aligner.mapopt.best_n = best_n;
    aligner.mapopt.pri_ratio = 0.0;
    let aligner = aligner
        .with_index(ref_path, save_index)
        .map_err(|e| Error::FailedAlignerInit { source: e })?;

    while let Some(sensor) = sensors.next() {
        let record = match sensor {
            Ok(record) => record,
            Err(err) => {
                warn!("Failed to parse sensor: {}", err);
                continue;
            }
        };

        let seq = record.normalize(allow_iupac);

        let mappings = match aligner.map(&seq, false, false, None, None, Some(record.id())) {
            Ok(mappings) => mappings,
            Err(err) => {
                let id = String::from_utf8_lossy(record.id());
                warn!("Failed to map sensor {}: {}", id, err);
                continue;
            }
        };

        if let Some((writer, format)) = &mut output {
            for mapping in &mappings {
                format(mapping, *writer)?;
            }
        }

        sensor_mappings.extend(mappings);
    }

    Ok(sensor_mappings)
}

/// Write a tab-delimited row for an accepted mapping.
pub fn mapping_to_paf(mapping: &Mapping, writer: &mut dyn Write) -> std::io::Result<()> {
    let fallback_query_len = NonZeroI32::new(mapping.query_end - mapping.query_start + 1)
        .unwrap_or(NonZeroI32::new(1).unwrap());
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        mapping
            .query_name
            .as_deref()
            .map_or("", std::string::String::as_str),
        mapping.query_len.unwrap_or(fallback_query_len),
        mapping.query_start,
        mapping.query_end,
        mapping.strand,
        mapping
            .target_name
            .as_deref()
            .map_or("", std::string::String::as_str),
        mapping.target_len,
        mapping.target_start,
        mapping.target_end,
        mapping.match_len,
        mapping.block_len,
        mapping.mapq,
    )
}

#[cfg(test)]
mod tests {
    use std::{num::NonZeroI32, sync::Arc};

    use minimap2::Strand;

    use super::*;

    #[test]
    fn default_formatter_writes_tsv_row() {
        let mapping = Mapping {
            query_name: Some(Arc::new("sensor1".to_string())),
            query_len: NonZeroI32::new(20),
            query_start: 1,
            query_end: 19,
            strand: Strand::Forward,
            target_name: Some(Arc::new("target1".to_string())),
            target_len: 100,
            target_start: 10,
            target_end: 28,
            target_id: 0,
            match_len: 18,
            block_len: 18,
            mapq: 60,
            is_primary: true,
            is_supplementary: false,
            is_spliced: false,
            trans_strand: None,
            alignment: None,
            segment_id: 0,
        };

        let mut output = Vec::new();

        mapping_to_paf(&mapping, &mut output).unwrap();

        assert_eq!(
            String::from_utf8(output).unwrap(),
            "sensor1\t20\t1\t19\t+\ttarget1\t100\t10\t28\t18\t18\t60\n"
        );
    }

    #[test]
    fn accepts_custom_formatter() {
        let prefix = "sensor";
        let format = |mapping: &Mapping, writer: &mut dyn Write| {
            writeln!(writer, "{}\t{}", prefix, mapping.query_start)
        };
        let mapping = Mapping {
            query_name: None,
            query_len: NonZeroI32::new(20),
            query_start: 4,
            query_end: 19,
            strand: Strand::Forward,
            target_name: None,
            target_len: 100,
            target_start: 10,
            target_end: 28,
            target_id: 0,
            match_len: 18,
            block_len: 18,
            mapq: 60,
            is_primary: true,
            is_supplementary: false,
            is_spliced: false,
            trans_strand: None,
            alignment: None,
            segment_id: 0,
        };
        let mut output = Vec::new();

        let format: &MappingFormatter<'_> = &format;
        format(&mapping, &mut output).unwrap();

        assert_eq!(String::from_utf8(output).unwrap(), "sensor\t4\n");
    }
}
