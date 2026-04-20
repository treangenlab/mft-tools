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
    DefineMapping(DefineMappingArgs), //Define a mapping table given kmer size and alphabet size and/or reduced  
    Optimize(OptimizeOrderArgs), // Optimize an ordering given a mapping table and desired transition nucleotides balance
    MFT(MinFrameTransitionArgs), //Call --> takes the genomes and sequencing data and does the viral variation analysis
}

#[derive(Args, Default)]
#[clap(about="Run hill-climbing optimization on ordering for MFT performance", arg_required_else_help = true)]
pub struct OptimizeOrderArgs {

    //sequence inputs
    #[clap(num_args=1.., short='g', long="genomes", help_heading = "SEQUENCE INPUT", help="Genome files (fasta format) to be tested for seeding performance")]
    pub genomes: Vec<String>

}

#[derive(Args, Default)]
#[clap(about="Function to build mapping table from kmers of size k to reduced alphabet", arg_required_else_help = true)]
pub struct DefineMappingArgs {

    //sequence inputs
    #[clap(num_args=1.., short='g', long="genomes", help_heading = "SEQUENCE INPUT", help="Genome files to be tested for seeding performance")]
    pub genomes: Vec<String>

}

#[derive(Args, Default)]
#[clap(about="Get min frame transformation of a sequence from fasta format", arg_required_else_help = true)]
pub struct MinFrameTransitionArgs {

    //sequence inputs
    #[clap(num_args=1.., short='g', long="genomes", help_heading = "SEQUENCE INPUT", help="Genome files to be tested for seeding performance")]
    pub genomes: Vec<String>

}

pub fn parse_args() -> Cli {
    Cli::parse()
}
