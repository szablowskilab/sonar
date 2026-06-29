//! RNA sensor design
//!
//! CellREADR and RADAR approaches allow for translation of an effector
//! gene in response to a target mRNA. Designed sense-edit-switch (ses) RNA
//! sensors are made for a specific target mRNA sequence. This crate provides
//! utilities for generating sesRNA sensor designs.
//!
//! ## Tools
//!
//! - [`generate_ses_lib`] - Initial candidate sensor library design
//! - [`map`] - Map candidate sensors to reference transcriptome

mod design;

pub mod error;
pub mod prelude;
pub mod specificity;

pub use design::{Candidate, DesignParams, StopEdit, generate_candidates};
