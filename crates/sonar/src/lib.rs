//! RNA sensor design
//!
//! CellREADR and RADAR approaches allow for translation of an effector
//! gene in response to a target mRNA. Designed sense-edit-switch (ses) RNA
//! sensors are made for a specific target mRNA sequence. This crate provides
//! utilities for generating sesRNA sensor designs.
//!
//! ## Tools
//!
//! - [`generate_candidates`] - Initial candidate sensor library design
//! - [`map`] - Map candidate sensors to reference transcriptome

mod design;

pub mod error;
pub mod mapper;
pub mod prelude;

pub use design::{Candidate, DesignParams, StopEdit, generate_candidates, sanitize_target_id};
pub use mapper::{map, map_with_writer, mapping_to_paf};
