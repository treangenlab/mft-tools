use crate::consts::*;
use clap::builder::styling::{AnsiColor, Effects, Styles};
use clap::{Args, Parser, Subcommand};

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
    DefineOrder(DefineOrderArgs), // Define an ordering to be used for optimization or transformation
    Evaluate(EvaluateArgs), // Evaluate the theoretical SNP masking rate of an ordering/mapping table pairing
    Optimize(OptimizeOrderArgs), // Optimize an ordering given a mapping table and desired transition nucleotides balance
    Transform(MinFrameTransitionArgs), //Call --> takes the genomes and sequencing data and does the viral variation analysis
}

#[derive(Args, Default)]
#[clap(
    about = "Function to evaluate the theoretical SNP masking rate of an ordering/mapping table pairing",
    arg_required_else_help = true
)]
pub struct EvaluateArgs {
    // K-mer size
    #[clap(
        short = 'k',
        long = "kmer-size",
        help_heading = "MFT Algorithm Params",
        default_value_t = 3,
        help = "Length of the k-mer (must match mapping table)"
    )]
    pub k: usize,

    // window size 
    #[clap(
        short = 'w',
        long = "window-size",
        help_heading = "MFT Algorithm Params",
        default_value_t = 5,
        help = "Length of the sequence window (w >= k)"
    )]
    pub w: usize,   

    // Mapping Table path
    #[clap(
        short = 'm',
        long = "mapping-table",
        help_heading = "MFT Mapping Table",
        help = "Path to the tab-delimited k-mer mapping file"
    )]
    pub mapping_table: String,

    // Ordering path
    #[clap(
        short = 'b',
        long = "ordering",
        help_heading = "MFT Ordering",
        help = "Path to the k-mer ordering file"
    )]
    pub ordering: String,

    // Number of samples for masking rate estimation
    #[clap(
        short = 'n',
        long = "n-samples",
        help_heading = "Evaluation Parameters",
        default_value_t = 10000,
        help = "Number of windows to sample for masking rate calculation"
    )]
    pub n_samples: usize,

    // Output path containing resulting information
    #[clap(
        short = 'o',
        long = "output",
        help_heading = "Output",
        default_value = "evaluation.txt",
        help = "Path to save the resulting evaluation information"
    )]
    pub output: String,

    //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,
}

#[derive(Args, Default)]
#[clap(
    about = "Function to define a k-mer ordering to be used for MFT transformation or baseline for further optimization",
    arg_required_else_help = true
)]
pub struct DefineOrderArgs {
    // k-mer size
    #[clap(
        short = 'k',
        long = "kmer-size",
        help_heading = "MFT Algorithm Params",
        default_value_t = 3,
        help = "Length of the k-mer (must match mapping table)"
    )]
    pub k: usize,

    // Random ordering seed
    #[clap(
        short = 's',
        long = "seed",
        help_heading = "Ordering Parameters",
        help = "Seed for random ordering"
    )]
    pub seed: Option<u64>,

    // Alphabetical ordering
    #[clap(
        short = 'a',
        long = "alphabetical",
        help_heading = "Ordering Parameters",
        help = "Flag to use alphabetical ordering instead of random"
    )]
    pub alphabetical: bool,

    // Output location
    #[clap(
        short = 'o',
        long = "output",
        help_heading = "Output",
        default_value = "order.txt",
        help = "Path where the resulting k-mer ordering will be stored"
    )]
    pub output: String,

    //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,
}

#[derive(Args, Default)]
#[clap(
    about = "Run hill-climbing optimization on ordering for MFT performance",
    arg_required_else_help = true
)]
pub struct OptimizeOrderArgs {
    // K-mer size
    #[clap(
        short = 'k',
        long = "kmer-size",
        help_heading = "MFT Algorithm Params",
        default_value_t = 3,
        help = "Length of the k-mer (must match mapping table)"
    )]
    pub k: usize,

    // Window size
    #[clap(
        short = 'w',
        long = "window-size",
        help_heading = "MFT Algorithm Params",
        default_value_t = 5,
        help = "Length of the sequence window (w >= k)"
    )]
    pub w: usize,

    // Mapping Table path
    #[clap(
        short = 'm',
        long = "mapping-table",
        help_heading = "MFT Algorithm Params",
        help = "Path to the tab-delimited k-mer mapping file"
    )]
    pub mapping_table: String,

    // Optional: Base Order path
    #[clap(
        short = 'b',
        long = "base-order",
        help_heading = "MFT Algorithm Params",
        help = "Initial k-mer ordering file to start optimization from"
    )]
    pub base_order: Option<String>,

    // Optional: Number of iterations
    #[clap(
        short = 'i',
        long = "iterations",
        help_heading = "Optimization Parameters",
        default_value_t = 4000,
        help = "Number of hill-climbing iterations"
    )]
    pub num_iterations: usize,

    // Optional: Number of iterations
    #[clap(
        short = 's',
        long = "samples",
        help_heading = "Optimization Parameters",
        default_value_t = 10000,
        help = "Number of samples used to calculate the masking rate at each iteration"
    )]
    pub num_samples: usize,

    // Output path
    #[clap(
        short = 'o',
        long = "output",
        help_heading = "Output",
        default_value = "optimized_order.txt",
        help = "Path to save the optimized k-mer ordering"
    )]
    pub output: String,

    //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,
}

#[derive(Args, Default)]
#[clap(
    about = "Function to build mapping table from kmers of size k to reduced alphabet",
    arg_required_else_help = true
)]
pub struct DefineMappingArgs {
    // K-mer size
    #[clap(
        short = 'k',
        long = "kmer-size",
        help_heading = "MFT Algorithm Params",
        help = "Length of the k-mer to be mapped"
    )]
    pub k: usize,

    // Reduced alphabet size
    #[clap(
        short = 'a',
        long = "alphabet-size",
        help_heading = "MFT Algorithm Params",
        help = "Size of the reduced alphabet (e.g., 20 for same size as amino acids)"
    )]
    pub alphabet_size: Option<usize>,

    // Spaced seed string
    #[clap(
        short = 's',
        long = "spaced-seed",
        help_heading = "MFT Algorithm Params",
        help = "Spaced seed pattern with 1s being match, 0 mismatch (e.g., 11011). Length must match k."
    )]
    pub spaced_seed: Option<String>,

    // Output location
    #[clap(
        short = 'o',
        long = "output",
        help_heading = "Output",
        help = "Path where the resulting mapping table will be stored"
    )]
    pub output: String,

    //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,
}

#[derive(Args, Default)]
#[clap(
    about = "Get min frame transformation of a sequence from fasta format",
    arg_required_else_help = true
)]
pub struct MinFrameTransitionArgs {
    //sequence inputs
    #[clap(num_args=1.., short='g', long="genomes", help_heading = "Sequence Input", help="Sequence files to be converted to MFT")]
    pub genomes: Vec<String>,

    // K-mer size
    #[clap(
        short = 'k',
        long = "kmer-size",
        help_heading = "MFT Algorithm Params",
        default_value_t = 3,
        help = "Length of the k-mer (must match mapping table)"
    )]
    pub k: usize,

    // Window size
    #[clap(
        short = 'w',
        long = "window-size",
        help_heading = "MFT Algorithm Params",
        default_value_t = 5,
        help = "Length of the sequence window (w >= k)"
    )]
    pub w: usize,

    // Mapping Table path
    #[clap(
        short = 'm',
        long = "mapping-table",
        help_heading = "MFT Mapping Table",
        help = "Path to the tab-delimited k-mer mapping file"
    )]
    pub mapping_table: String,

    // Optional: Base Order path
    #[clap(
        short = 'b',
        long = "order",
        help_heading = "MFT Ordering",
        help = "Path to file containing pre-existing k-mer ordering. If not defined a random ordering will be used"
    )]
    pub order: Option<String>,

    #[clap(
        short = 's',
        long = "seed",
        help_heading = "MFT Ordering",
        help = "Seed for random ordering (will not be used if -b is provided)"
    )]
    pub seed: Option<u64>,

    // Output location
    #[clap(
        short = 'o',
        long = "output",
        help_heading = "Output",
        default_value = "mft.fa",
        help = "Path where the resulting fasta files will be outputted too"
    )]
    pub output: String,

    //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,
}

pub fn parse_args() -> Cli {
    Cli::parse()
}
