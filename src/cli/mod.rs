pub mod error;
pub mod sensors;
pub mod table;

use clap::{Parser, Subcommand};

use sensors::Args as SensorArgs;

#[derive(Debug, Parser)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Generate an RNA sensor library for a given target
    #[clap(arg_required_else_help = true)]
    Gen(SensorArgs),
    // /// Filter RNA sensor candidates
    // #[clap(arg_required_else_help = true)]
    // Filter(SensorArgs),
}
