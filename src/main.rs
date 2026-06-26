mod cli;

use clap::Parser as _;

pub use cli::error::Result;

use cli::{Cli, Commands};

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Gen(args) => cli::sensors::generate_candidates(args),
    }
}
