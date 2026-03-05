use crate::consts::*;
use clap::{Args, Parser, Subcommand};
use clap::builder::styling::{AnsiColor, Effects, Styles};

fn custom_styles() -> Styles {
    Styles::styled()
        .header(AnsiColor::Blue.on_default() | Effects::BOLD)
        .usage(AnsiColor::White.on_default() | Effects::BOLD)
        .literal(AnsiColor::White.on_default() | Effects::BOLD)
        .placeholder(AnsiColor::White.on_default())
}

#[derive(Parser)]
#[command(author, version, about, long_about = None, styles=custom_styles())]
#[command(propagate_version = true)]
pub struct Cli {
    #[clap(subcommand)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    Complete(CompleteArgs), // Run all and compare
    Sampling(SampleArgs), //Sampling --> parses genomes and calculates conservation
    Mums(MUMArgs), //Call --> takes the genomes and sequencing data and does the viral variation analysis
}

#[derive(Args, Default)]
#[clap(about="Run all seeding methods against a set of sequences and evaluate differences", arg_required_else_help = true)]
pub struct CompleteArgs {

    //sequence inputs
    #[clap(num_args=1.., short='g', long="genomes", help_heading = "SEQUENCE INPUT", help="Genome files (fasta format) to be tested for seeding performance")]
    pub genomes: Vec<String>

}

#[derive(Args, Default)]
#[clap(about="Run sampling-based methods against a set of sequences and get statistics", arg_required_else_help = true)]
pub struct SampleArgs {

    //sequence inputs
    #[clap(num_args=1.., short='g', long="genomes", help_heading = "SEQUENCE INPUT", help="Genome files to be tested for seeding performance")]
    pub genomes: Vec<String>

}

#[derive(Args, Default)]
#[clap(about="Run MUM-based methods against a set of sequences and get statistics", arg_required_else_help = true)]
pub struct MUMArgs {

    //sequence inputs
    #[clap(num_args=1.., short='g', long="genomes", help_heading = "SEQUENCE INPUT", help="Genome files to be tested for seeding performance")]
    pub genomes: Vec<String>

}

pub fn parse_args() -> Cli {
    Cli::parse()
}
