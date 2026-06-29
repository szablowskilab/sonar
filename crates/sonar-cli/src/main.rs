mod design;
mod error;
mod mapper;
mod table;

use clap::{Parser, Subcommand};

pub use crate::error::Result;

#[derive(Debug, Parser)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Generate an RNA sensor library for a given target
    #[clap(arg_required_else_help = true)]
    Gen(design::Args),

    /// Map RNA sensor candidates to a reference transcriptome for specificity filtering
    #[clap(arg_required_else_help = true)]
    Map(mapper::Args),

    /// Append adapters to an RNA sensor library
    #[clap(arg_required_else_help = true)]
    Append(design::AppendArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Gen(args) => design::generate(args),
        Commands::Map(args) => mapper::map_to_ref(args),
        Commands::Append(args) => design::append_adapters(args),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory as _;

    #[test]
    fn verify_cli() {
        Cli::command().debug_assert();
    }
}
