use std::{cmp::Ordering, collections::HashSet, ops::RangeInclusive};

use crate::error::{Error, Result};

pub struct DesignParams {
    /// 1-based inclusive target subregion, e.g. 10:350.
    pub region: Option<RangeInclusive<usize>>,

    /// sesRNA length range in nucleotides as min:max.
    pub ses_length: RangeInclusive<usize>,

    /// Number of designed TAG stop codons to include.
    pub stop_count: usize,

    /// Allowed design stop position window as min:max from the sesRNA 5' end.
    pub stop_window: RangeInclusive<usize>,

    /// Minimum distance between designed TAG and any converted stop codons.
    pub min_stop_distance: usize,

    /// Translation frame for stop and ATG checks.
    pub frame: usize,

    /// Max number of candidates to output.
    pub max_candidates: usize,

    /// Do not spread out candidate designs across the target sequence.
    pub no_spread: bool,

    /// Allow IUPAC nucleotide codes in the target sequence.
    pub allow_iupac: bool,
}

/// Records conversion of a stop codon to a non-stop.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StopEdit {
    pub position: usize,
    pub from: String,
    pub to: String,
}

/// A candidate RNA sensor.
#[derive(Debug, Clone, PartialEq)]
pub struct Candidate {
    pub id: String,
    pub target_id: String,
    pub target_start: usize,
    pub target_end: usize,
    pub ses_length: usize,
    pub seed_target_pos: Vec<usize>,
    pub seed_sequence: String,
    pub designed_stop_pos: Vec<usize>,
    pub final_sesrna: String,
    pub edited_stops: Vec<StopEdit>,
    pub downstream_atg: bool,
    pub gc_content: f64,
    pub score: f64,
    pub fail_reason: String,
    pub designed_stop_count: usize,
    pub window_index: usize,
    pub seed_index: usize,
}

/// Generate a library of candidate RNA sensors from the given target information and design parameters.
pub fn generate_candidates(
    id: &str,
    target_seq: &[u8],
    design_params: &DesignParams,
) -> Result<Vec<Candidate>> {
    validate_rules(design_params)?;

    let (region_start0, region_end1) = get_ses_region(target_seq, design_params.region.as_ref())?;
    let ses_length = get_mean_ses_length(&design_params.ses_length, target_seq.len())?;
    let seeds = find_seeds(target_seq, region_start0, region_end1);
    let mut candidates = Vec::new();
    let mut candidate_index = 0usize;
    let mut seen_candidates = HashSet::new();

    for (seed_index, &seed_start0) in seeds.iter().enumerate() {
        let placements = valid_window_starts(
            seed_start0,
            region_start0,
            region_end1,
            ses_length,
            design_params.stop_window.clone(),
        );

        for (window_index, &window_start0) in placements.iter().enumerate() {
            let window_end0 = window_start0 + ses_length;
            if window_end0 > target_seq.len() {
                continue;
            }

            let window_seq = &target_seq[window_start0..window_end0];
            let mut ses = reverse_complement_bytes(window_seq);
            let mut eligible_stops = eligible_designed_stops(
                target_seq,
                &ses,
                window_start0,
                window_end0,
                &design_params.stop_window,
            );
            if eligible_stops.is_empty() {
                continue;
            }

            sort_stops_by_center(&mut eligible_stops, &design_params.stop_window);
            eligible_stops.truncate(design_params.stop_count);
            eligible_stops.sort_by_key(|stop| stop.stop0);

            let designed_stops = eligible_stops
                .iter()
                .map(|stop| stop.stop0)
                .collect::<Vec<_>>();
            if !seen_candidates.insert((window_start0, designed_stops.clone())) {
                continue;
            }

            for &stop_pos in &designed_stops {
                ses[stop_pos..stop_pos + 3].copy_from_slice(b"TAG");
            }

            let edits = convert_in_frame_stops(&mut ses, design_params.frame, &designed_stops);
            if has_close_edit(&edits, &designed_stops, design_params.min_stop_distance) {
                continue;
            }
            if has_downstream_atg(&ses, design_params.frame, &designed_stops) {
                continue;
            }

            candidate_index += 1;
            let gc = gc_content(&ses);
            let edit_count = edits.len();
            candidates.push(Candidate {
                id: format!("cand_{candidate_index:04}"),
                target_id: id.to_string(),
                target_start: window_start0 + 1,
                target_end: window_end0,
                ses_length: ses.len(),
                seed_target_pos: eligible_stops
                    .iter()
                    .map(|stop| stop.seed_start0 + 1)
                    .collect(),
                seed_sequence: "CCA".to_string(),
                designed_stop_pos: designed_stops.iter().map(|stop0| stop0 + 1).collect(),
                final_sesrna: String::from_utf8(ses).expect("sequence should remain ASCII"),
                edited_stops: edits,
                downstream_atg: false,
                gc_content: gc,
                score: score_candidate(
                    gc,
                    edit_count,
                    &designed_stops,
                    design_params.stop_window.clone(),
                ),
                fail_reason: String::new(),
                designed_stop_count: designed_stops.len(),
                window_index,
                seed_index,
            });
        }
    }

    candidates.sort_by(compare_candidates);
    limit_candidates(
        &mut candidates,
        design_params.max_candidates,
        design_params.no_spread,
        region_start0,
        region_end1,
    );

    Ok(candidates)
}

fn validate_rules(rules: &DesignParams) -> Result<()> {
    if rules.stop_count < 1 {
        return Err(Error::InvalidStopCount);
    }
    if rules.frame > 2 {
        return Err(Error::InvalidFrame);
    }
    if *rules.ses_length.start() == 0 || rules.ses_length.end() < rules.ses_length.start() {
        return Err(Error::InvalidSesLengthRange);
    }
    if *rules.stop_window.start() == 0 || rules.stop_window.end() < rules.stop_window.start() {
        return Err(Error::InvalidStopWindowRange);
    }

    Ok(())
}

/// Get the SES region as a 0-based start and exclusive end from a 1-based
/// inclusive region range supplied by the user.
fn get_ses_region(
    target_seq: &[u8],
    region: Option<&RangeInclusive<usize>>,
) -> Result<(usize, usize)> {
    let mut start = 0usize;
    let mut end = target_seq.len();
    if let Some(region) = region {
        if *region.start() < 1 || *region.end() > end || region.start() > region.end() {
            return Err(Error::InvalidRegion {
                start: *region.start(),
                end: *region.end(),
            });
        }

        start = *region.start() - 1;
        end = *region.end();
    }

    Ok((start, end))
}

/// Get the average SES length within the given range, ensuring it does not exceed the target length.
fn get_mean_ses_length(
    ses_length_range: &RangeInclusive<usize>,
    target_length: usize,
) -> Result<usize> {
    let mean_length = (*ses_length_range.start() + *ses_length_range.end()) / 2;
    if mean_length < *ses_length_range.start() || mean_length > *ses_length_range.end() {
        Err(Error::SesLengthOutsideRange {
            mean_length,
            start: *ses_length_range.start(),
            end: *ses_length_range.end(),
        })
    } else if mean_length > target_length {
        Err(Error::SesLongerThanTarget {
            ses_length: mean_length,
            target_length,
        })
    } else {
        Ok(mean_length)
    }
}

/// For candidate comparison and sorting.
fn compare_candidates(a: &Candidate, b: &Candidate) -> Ordering {
    match b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal) {
        Ordering::Equal => match a.edited_stops.len().cmp(&b.edited_stops.len()) {
            Ordering::Equal => {
                let a_gc = (a.gc_content - 50.0).abs();
                let b_gc = (b.gc_content - 50.0).abs();
                match a_gc.partial_cmp(&b_gc).unwrap_or(Ordering::Equal) {
                    Ordering::Equal => a.id.cmp(&b.id),
                    other => other,
                }
            }
            other => other,
        },
        other => other,
    }
}

fn limit_candidates(
    candidates: &mut Vec<Candidate>,
    max_candidates: usize,
    no_spread: bool,
    region_start0: usize,
    region_end1: usize,
) {
    if candidates.len() <= max_candidates {
        return;
    }

    if no_spread {
        candidates.truncate(max_candidates);
    } else {
        *candidates =
            select_spread_candidates(candidates, max_candidates, region_start0, region_end1);
    }
}

fn select_spread_candidates(
    candidates: &[Candidate],
    max_candidates: usize,
    region_start0: usize,
    region_end1: usize,
) -> Vec<Candidate> {
    if max_candidates == 0 || candidates.is_empty() {
        return Vec::new();
    }

    let region_len = region_end1.saturating_sub(region_start0);
    let mut selected_indices = Vec::new();
    let mut selected_index_set = HashSet::new();

    if region_len > 0 {
        for bin_index in 0..max_candidates {
            let bin_start = region_start0 + (region_len * bin_index) / max_candidates;
            let bin_end = region_start0 + (region_len * (bin_index + 1)) / max_candidates;

            if let Some(candidate_index) =
                candidates
                    .iter()
                    .enumerate()
                    .find_map(|(candidate_index, candidate)| {
                        (!selected_index_set.contains(&candidate_index)
                            && candidate_center0(candidate) >= bin_start
                            && candidate_center0(candidate) < bin_end)
                            .then_some(candidate_index)
                    })
            {
                selected_indices.push(candidate_index);
                selected_index_set.insert(candidate_index);
            }
        }
    }

    for candidate_index in 0..candidates.len() {
        if selected_indices.len() == max_candidates {
            break;
        }
        if selected_index_set.insert(candidate_index) {
            selected_indices.push(candidate_index);
        }
    }

    let mut selected = selected_indices
        .into_iter()
        .map(|index| candidates[index].clone())
        .collect::<Vec<_>>();
    selected.sort_by(|a, b| {
        a.target_start
            .cmp(&b.target_start)
            .then_with(|| a.target_end.cmp(&b.target_end))
            .then_with(|| a.id.cmp(&b.id))
    });
    selected
}

fn candidate_center0(candidate: &Candidate) -> usize {
    ((candidate.target_start - 1) + candidate.target_end) / 2
}

/// Find CCA seeds in the target sequence.
pub fn find_seeds(target: &[u8], start0: usize, end1: usize) -> Vec<usize> {
    let mut seeds = Vec::new();
    let end1 = end1.min(target.len());
    for i in start0..end1.saturating_sub(2) {
        if &target[i..i + 3] == b"CCA" {
            seeds.push(i);
        }
    }
    seeds
}

/// Returns a sorted list of valid window starting positions for the given seed and stop window.
pub fn valid_window_starts(
    seed_start0: usize,
    allowed_start0: usize,
    allowed_end1: usize,
    length: usize,
    stop_window: RangeInclusive<usize>,
) -> Vec<usize> {
    let mut starts = Vec::new();
    for stop_pos1 in stop_window {
        let stop_pos0 = stop_pos1 as isize - 1;
        let window_start0 = seed_start0 as isize + 3 - length as isize + stop_pos0;
        if window_start0 < 0 {
            continue;
        }

        let window_start0 = window_start0 as usize;
        if window_start0 < allowed_start0 || window_start0 + length > allowed_end1 {
            continue;
        }

        starts.push(window_start0);
    }

    starts.sort_unstable();
    starts.dedup();
    starts
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct DesignedStop {
    seed_start0: usize,
    stop0: usize,
}

fn eligible_designed_stops(
    target: &[u8],
    ses: &[u8],
    window_start0: usize,
    window_end0: usize,
    stop_window: &RangeInclusive<usize>,
) -> Vec<DesignedStop> {
    find_seeds(target, window_start0, window_end0)
        .into_iter()
        .filter_map(|seed_start0| {
            let stop0 = seed_to_stop_pos0(seed_start0, window_start0, ses.len())?;
            if stop0 + 3 > ses.len() || !stop_window.contains(&(stop0 + 1)) {
                return None;
            }
            if &ses[stop0..stop0 + 3] != b"TGG" {
                return None;
            }
            Some(DesignedStop { seed_start0, stop0 })
        })
        .collect()
}

fn seed_to_stop_pos0(seed_start0: usize, window_start0: usize, ses_length: usize) -> Option<usize> {
    let stop0 = ses_length as isize - ((seed_start0 as isize - window_start0 as isize) + 3);
    if stop0 < 0 {
        None
    } else {
        Some(stop0 as usize)
    }
}

fn sort_stops_by_center(stops: &mut [DesignedStop], stop_window: &RangeInclusive<usize>) {
    let center = (*stop_window.start() + *stop_window.end()) as f64 / 2.0;
    stops.sort_by(|a, b| {
        let a_distance = ((a.stop0 + 1) as f64 - center).abs();
        let b_distance = ((b.stop0 + 1) as f64 - center).abs();
        a_distance
            .partial_cmp(&b_distance)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.stop0.cmp(&b.stop0))
            .then_with(|| a.seed_start0.cmp(&b.seed_start0))
    });
}

/// Converts in-frame stop codons (that are not designed from CCA seeds) to non-stop codons.
pub fn convert_in_frame_stops(
    sequence: &mut [u8],
    frame: usize,
    designed_stop_starts0: &[usize],
) -> Vec<StopEdit> {
    let mut edits = Vec::new();
    for i in (frame..sequence.len().saturating_sub(2)).step_by(3) {
        if designed_stop_starts0.contains(&i) {
            continue;
        }

        let codon = &sequence[i..i + 3];
        let replacement = match codon {
            b"TAG" | b"TAA" => Some(b"TAC".as_slice()),
            b"TGA" => Some(b"TCA".as_slice()),
            _ => None,
        };
        let Some(replacement) = replacement else {
            continue;
        };

        let from = String::from_utf8(codon.to_vec()).expect("ASCII codon");
        sequence[i..i + 3].copy_from_slice(replacement);
        edits.push(StopEdit {
            position: i + 1,
            from,
            to: String::from_utf8(replacement.to_vec()).expect("ASCII codon"),
        });
    }
    edits
}

/// Returns whether there is a close stop codon edit within the given distance of any designed stop.
pub fn has_close_edit(
    edits: &[StopEdit],
    designed_stop_starts0: &[usize],
    min_distance: usize,
) -> bool {
    edits.iter().any(|edit| {
        designed_stop_starts0
            .iter()
            .any(|&stop0| edit.position.abs_diff(stop0 + 1) < min_distance)
    })
}

/// Returns whether there is a downstream ATG codon after the first designed stop.
pub fn has_downstream_atg(sequence: &[u8], frame: usize, designed_stop_starts0: &[usize]) -> bool {
    let Some(&first_stop) = designed_stop_starts0.iter().min() else {
        return false;
    };
    for i in (frame..sequence.len().saturating_sub(2)).step_by(3) {
        if i <= first_stop {
            continue;
        }
        if &sequence[i..i + 3] == b"ATG" {
            return true;
        }
    }
    false
}

/// Calculate GC content of the given sequence.
pub fn gc_content(sequence: &[u8]) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }
    let gc = sequence
        .iter()
        .filter(|&&base| matches!(base, b'G' | b'C'))
        .count();
    100.0 * (gc as f64) / (sequence.len() as f64)
}

/// Scores the candidate sensor based on GC content, edit count, and number of designed stop codons.
pub fn score_candidate(
    gc: f64,
    edit_count: usize,
    designed_stops0: &[usize],
    stop_window: RangeInclusive<usize>,
) -> f64 {
    let center = (*stop_window.start() + *stop_window.end()) as f64 / 2.0;
    let stop_penalty = if designed_stops0.is_empty() {
        0.0
    } else {
        designed_stops0
            .iter()
            .map(|stop0| ((*stop0 + 1) as f64 - center).abs())
            .sum::<f64>()
            / designed_stops0.len() as f64
    };
    let additional_stop_bonus = designed_stops0.len().saturating_sub(1) as f64 * 10.0;
    let gc_penalty = (gc - 50.0).abs();
    1000.0 + additional_stop_bonus
        - (gc_penalty * 10.0)
        - ((edit_count as f64) * 25.0)
        - stop_penalty
}

/// Takes the reverse complement of the given sequence as a byte slice.
pub fn reverse_complement_bytes(sequence: &[u8]) -> Vec<u8> {
    let mut rev = vec![0_u8; sequence.len()];

    for (i, &base) in sequence.iter().enumerate() {
        rev[sequence.len() - 1 - i] = complement(base);
    }
    rev
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => base,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_rules() -> DesignParams {
        DesignParams {
            region: None,
            ses_length: 120..=120,
            stop_count: 1,
            stop_window: 1..=1,
            min_stop_distance: 10,
            frame: 0,
            max_candidates: 10,
            no_spread: false,
            allow_iupac: false,
        }
    }

    fn target_with_cca_seeds(length: usize, seed_offsets0: &[usize]) -> String {
        let mut target = vec![b'A'; length];
        for &offset in seed_offsets0 {
            target[offset..offset + 3].copy_from_slice(b"CCA");
        }
        String::from_utf8(target).unwrap()
    }

    fn test_candidate(id: &str, target_start: usize, target_end: usize, score: f64) -> Candidate {
        Candidate {
            id: id.to_string(),
            target_id: "target1".to_string(),
            target_start,
            target_end,
            ses_length: target_end - target_start + 1,
            seed_target_pos: Vec::new(),
            seed_sequence: "CCA".to_string(),
            designed_stop_pos: Vec::new(),
            final_sesrna: String::new(),
            edited_stops: Vec::new(),
            downstream_atg: false,
            gc_content: 50.0,
            score,
            fail_reason: String::new(),
            designed_stop_count: 1,
            window_index: 0,
            seed_index: 0,
        }
    }

    #[test]
    fn region_range_is_zero_indexed() {
        let target = b"ACTCCATAGAGTCCA";
        assert_eq!(get_ses_region(target, Some(&(4..=6))).unwrap(), (3, 6));
    }

    #[test]
    fn mean_ses_length_uses_range_average() {
        assert_eq!(get_mean_ses_length(&(200..=300), 500).unwrap(), 250);
    }

    #[test]
    fn get_all_valid_seeds() {
        let target = b"ACTCCATAGAGTCCA";
        let all_seeds = find_seeds(target, 0, target.len());
        let single_seed = find_seeds(target, 2, 6);
        let no_seeds = find_seeds(b"ACGTCCG", 0, target.len());

        assert_eq!(all_seeds, vec![3, 12]);
        assert_eq!(single_seed, vec![3]);
        assert_eq!(no_seeds, vec![]);
    }

    #[test]
    fn reverse_complement_sanity_check() {
        let got = String::from_utf8(reverse_complement_bytes(b"ACGTCCA")).unwrap();
        assert_eq!(got, "TGGACGT");
    }

    #[test]
    fn no_spread_limits_to_top_scoring_candidates() {
        let mut candidates = vec![
            test_candidate("cand_0003", 201, 220, 30.0),
            test_candidate("cand_0001", 1, 20, 100.0),
            test_candidate("cand_0002", 101, 120, 50.0),
        ];

        candidates.sort_by(compare_candidates);
        limit_candidates(&mut candidates, 2, true, 0, 300);

        let ids = candidates
            .iter()
            .map(|candidate| candidate.id.as_str())
            .collect::<Vec<_>>();
        assert_eq!(ids, vec!["cand_0001", "cand_0002"]);
    }

    #[test]
    fn spread_limits_to_best_candidate_per_region_bin() {
        let mut candidates = vec![
            test_candidate("early_second", 31, 50, 90.0),
            test_candidate("late", 281, 300, 40.0),
            test_candidate("middle", 141, 160, 50.0),
            test_candidate("early_best", 1, 20, 100.0),
        ];

        candidates.sort_by(compare_candidates);
        limit_candidates(&mut candidates, 3, false, 0, 300);

        let ids = candidates
            .iter()
            .map(|candidate| candidate.id.as_str())
            .collect::<Vec<_>>();
        assert_eq!(ids, vec!["early_best", "middle", "late"]);
    }

    #[test]
    fn spread_fills_empty_bins_with_best_unselected_candidates() {
        let mut candidates = vec![
            test_candidate("early_second", 31, 50, 90.0),
            test_candidate("late", 281, 300, 80.0),
            test_candidate("early_best", 1, 20, 100.0),
            test_candidate("early_third", 61, 80, 70.0),
        ];

        candidates.sort_by(compare_candidates);
        limit_candidates(&mut candidates, 3, false, 0, 300);

        let ids = candidates
            .iter()
            .map(|candidate| candidate.id.as_str())
            .collect::<Vec<_>>();
        assert_eq!(ids, vec!["early_best", "early_second", "late"]);
    }

    #[test]
    fn convert_inframe_stops_preserve_designed_stops() {
        let mut sequence = b"TAGTAATGA".to_vec();
        let edits = convert_in_frame_stops(&mut sequence, 0, &[0]);
        assert_eq!(String::from_utf8(sequence).unwrap(), "TAGTACTCA");
        assert_eq!(edits.len(), 2);
        assert_eq!(edits[0].position, 4);
        assert_eq!(edits[0].from, "TAA");
        assert_eq!(edits[0].to, "TAC");
        assert_eq!(edits[1].position, 7);
        assert_eq!(edits[1].from, "TGA");
        assert_eq!(edits[1].to, "TCA");
    }

    #[test]
    fn convert_stop_in_selected_frame() {
        let mut sequence = b"ATAGAA".to_vec();
        let edits = convert_in_frame_stops(&mut sequence, 1, &[]);
        assert_eq!(String::from_utf8(sequence).unwrap(), "ATACAA");
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].position, 2);
    }

    #[test]
    fn has_downstream_start() {
        assert!(has_downstream_atg(b"TAGCAAATG", 0, &[0]));
        assert!(!has_downstream_atg(b"ATGTAGCAA", 0, &[3]));
    }

    #[test]
    fn has_downstream_start_after_first_designed_stop() {
        assert!(has_downstream_atg(b"TAGAAAATGAAATAG", 0, &[0, 12]));
    }

    #[test]
    fn has_too_close_edit() {
        let edits = vec![StopEdit {
            position: 8,
            from: "TAA".to_string(),
            to: "TAC".to_string(),
        }];
        assert!(has_close_edit(&edits, &[9], 3));
        assert!(!has_close_edit(&edits, &[20], 3));
    }

    #[test]
    fn gen_candidate_lib() {
        let target_seq = format!("{}CCA", "A".repeat(117));
        let candidates =
            generate_candidates("target1", target_seq.as_bytes(), &default_rules()).unwrap();
        assert_eq!(candidates.len(), 1);

        let cand = &candidates[0];
        assert_eq!(cand.target_start, 1);
        assert_eq!(cand.target_end, 120);
        assert_eq!(cand.seed_target_pos, vec![118]);
        assert_eq!(cand.designed_stop_pos, vec![1]);
        assert!(cand.final_sesrna.starts_with("TAG"));
        assert!(!cand.downstream_atg);

        assert_eq!(cand.id, "cand_0001");
        assert_eq!(cand.designed_stop_count, 1);
    }

    #[test]
    fn region_limits_candidate_windows() {
        let target_seq = format!("{}CCA{}CCA", "A".repeat(117), "A".repeat(117));
        let mut rules = default_rules();
        rules.region = Some(121..=240);
        let candidates = generate_candidates("target1", target_seq.as_bytes(), &rules).unwrap();
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].target_start, 121);
        assert_eq!(candidates[0].seed_target_pos, vec![238]);
    }

    #[test]
    fn generate_candidate_with_multiple_designed_stops() {
        let target_seq = target_with_cca_seeds(20, &[6, 10]);
        let mut rules = default_rules();
        rules.ses_length = 20..=20;
        rules.stop_count = 2;
        rules.stop_window = 8..=12;

        let candidates = generate_candidates("target1", target_seq.as_bytes(), &rules).unwrap();
        assert_eq!(candidates.len(), 1);

        let cand = &candidates[0];
        assert_eq!(cand.designed_stop_pos, vec![8, 12]);
        assert_eq!(cand.seed_target_pos, vec![11, 7]);
        assert_eq!(cand.designed_stop_count, 2);
        assert_eq!(&cand.final_sesrna[7..10], "TAG");
        assert_eq!(&cand.final_sesrna[11..14], "TAG");
    }

    #[test]
    fn select_most_centered_stops_up_to_stop_count() {
        let target_seq = target_with_cca_seeds(24, &[7, 10, 13]);
        let mut rules = default_rules();
        rules.ses_length = 24..=24;
        rules.stop_count = 2;
        rules.stop_window = 9..=15;

        let candidates = generate_candidates("target1", target_seq.as_bytes(), &rules).unwrap();
        assert_eq!(candidates.len(), 1);

        let cand = &candidates[0];
        assert_eq!(cand.designed_stop_pos, vec![9, 12]);
        assert_eq!(cand.seed_target_pos, vec![14, 11]);
        assert_eq!(cand.designed_stop_count, 2);
    }

    #[test]
    fn accept_fewer_designed_stops_than_stop_count() {
        let target_seq = target_with_cca_seeds(24, &[7, 10, 13]);
        let mut rules = default_rules();
        rules.ses_length = 24..=24;
        rules.stop_count = 4;
        rules.stop_window = 9..=15;

        let candidates = generate_candidates("target1", target_seq.as_bytes(), &rules).unwrap();
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].designed_stop_pos, vec![9, 12, 15]);
        assert_eq!(candidates[0].designed_stop_count, 3);
    }

    #[test]
    fn reject_candidate_with_downstream_start() {
        let target_seq = format!("CAT{}CCA", "A".repeat(114));
        let mut rules = default_rules();
        rules.stop_window = 4..=4;
        let candidates = generate_candidates("target2", target_seq.as_bytes(), &rules).unwrap();
        assert!(candidates.is_empty());
    }

    #[test]
    fn reject_candidate_with_downstream_start_between_designed_stops() {
        let target_seq =
            String::from_utf8(reverse_complement_bytes(b"TGGAAAATGAAATGGAAA")).unwrap();
        let mut rules = default_rules();
        rules.ses_length = 18..=18;
        rules.stop_count = 2;
        rules.stop_window = 1..=13;

        let candidates = generate_candidates("target2", target_seq.as_bytes(), &rules).unwrap();
        assert!(candidates.is_empty());
    }

    #[test]
    fn reject_candidate_with_nearby_edited_stop() {
        let target_seq = format!("{}TTA{}CCA", "A".repeat(4), "A".repeat(110));
        let mut rules = default_rules();
        rules.stop_window = 4..=4;
        let candidates = generate_candidates("target3", target_seq.as_bytes(), &rules).unwrap();
        assert!(candidates.is_empty());
    }
}
