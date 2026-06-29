use sonar::Candidate;

pub fn render_tsv(candidates: &[Candidate]) -> String {
    let mut out = String::new();
    out.push_str(tsv_header());
    out.push('\n');
    for candidate in candidates {
        out.push_str(&candidate_tsv(candidate));
        out.push('\n');
    }
    out
}

fn candidate_tsv(candidate: &Candidate) -> String {
    let designed_stops = candidate
        .designed_stop_pos
        .iter()
        .map(|pos| pos.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let seed_positions = candidate
        .seed_target_pos
        .iter()
        .map(|pos| pos.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let edits = candidate
        .edited_stops
        .iter()
        .map(|edit| format!("{}:{}>{}", edit.position, edit.from, edit.to))
        .collect::<Vec<_>>()
        .join(";");

    [
        candidate.id.clone(),
        candidate.target_id.clone(),
        candidate.target_start.to_string(),
        candidate.target_end.to_string(),
        candidate.ses_length.to_string(),
        seed_positions,
        designed_stops,
        edits,
        format!("{:.2}", candidate.gc_content),
        format!("{:.2}", candidate.score),
    ]
    .join("\t")
}

fn tsv_header() -> &'static str {
    "candidate_id\ttarget_id\ttarget_start\ttarget_end\tses_length\tseed_target_pos\tdesigned_stop_pos\tedited_stops\tgc_content\tscore"
}

#[cfg(test)]
mod tests {
    use super::*;
    use sonar::{Candidate, StopEdit};

    #[test]
    fn render_candidate_table() {
        let candidate = Candidate {
            id: "cand_0001".to_string(),
            target_id: "target1".to_string(),
            target_start: 1,
            target_end: 120,
            ses_length: 120,
            seed_target_pos: vec![118, 100],
            designed_stop_pos: vec![1, 19],
            final_sesrna: "TAG".to_string(),
            edited_stops: vec![StopEdit {
                position: 10,
                from: "TAA".to_string(),
                to: "TAC".to_string(),
            }],
            gc_content: 50.0,
            score: 1000.0,
            designed_stop_count: 2,
            window_index: 0,
            seed_index: 0,
        };

        let table = render_tsv(&[candidate]);
        assert!(table.starts_with("candidate_id\ttarget_id"));
        assert!(table.contains("cand_0001\ttarget1"));
        assert!(table.contains("\t118,100\t1,19\t"));
        assert!(table.contains("10:TAA>TAC"));
    }
}
