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
    pub genomes: Vec<String>,

        //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,

}

#[derive(Args, Default)]
#[clap(about="Function to build mapping table from kmers of size k to reduced alphabet", arg_required_else_help = true)]
pub struct DefineMappingArgs {

    // K-mer size
    #[clap(short='k', long="kmer-size", help_heading="ALGORITHM", help="Length of the k-mer to be mapped")]
    pub k: usize,

    // Reduced alphabet size
    #[clap(short='a', long="alphabet-size", help_heading="ALGORITHM", help="Size of the reduced alphabet (e.g., 20 for same size as amino acids)")]
    pub alphabet_size: Option<usize>,

    // Spaced seed string
    #[clap(short='s', long="spaced-seed", help_heading="ALGORITHM", help="Spaced seed pattern with 1s being match, 0 mismatch (e.g., 11011). Length must match k.")]
    pub spaced_seed: Option<String>,

    // Output location
    #[clap(short='o', long="output", help_heading="Output", help="Path where the resulting mapping table will be stored")]
    pub output: String,

    //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,

}

#[derive(Args, Default)]
#[clap(about="Get min frame transformation of a sequence from fasta format", arg_required_else_help = true)]
pub struct MinFrameTransitionArgs {

    //sequence inputs
    #[clap(num_args=1.., short='g', long="genomes", help_heading = "SEQUENCE INPUT", help="Genome files to be tested for seeding performance")]
    pub genomes: Vec<String>,

        //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,

}

pub fn parse_args() -> Cli {
    Cli::parse()
}
