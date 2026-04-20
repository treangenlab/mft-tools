use crate::cli::*;
use log::*;
use std::path::Path;

pub fn check_args(args: &DefineMappingArgs){
    let output_level;
    if args.verbose {
        output_level = log::LevelFilter::Trace;
    } else {
        output_level = log::LevelFilter::Info;
    }

    simple_logger::SimpleLogger::new()
        .with_level(output_level)
        .init()
        .unwrap();

    // limits for kmer sizes
    if args.k < 3 || args.k > 15 {
        error!("Invalid k-mer size: {}. Must be between 3 and 15.", args.k);
        std::process::exit(1);
    }

    // Mutual Exclusion -- alphabet size and Spaced Seed
    if args.alphabet_size.is_some() && args.spaced_seed.is_some() {
        error!("Provide either --alphabet-size or --spaced-seed, not both.");
        std::process::exit(1);
    }
    if args.alphabet_size.is_none() && args.spaced_seed.is_none() {
        error!("You must provide either --alphabet-size or --spaced-seed.");
        std::process::exit(1);
    }

    // Alphabet size check: 4-64 (could update later)
    if let Some(size) = args.alphabet_size {
        if size < 4 || size > 64 {
            error!("Alphabet size {} is out of bounds. Must be 4-64.", size);
            std::process::exit(1);
        }
    }

    // Spaced seed checks: length == k and weight (number of 1s/Xs) between 3-5
    if let Some(ref seed) = args.spaced_seed {
        if seed.len() != args.k {
            error!("Spaced seed length ({}) must match k ({}).", seed.len(), args.k);
            std::process::exit(1);
        }

        // Assuming '1' or 'X' represents a match (weight)
        let weight = seed.chars().filter(|&c| c == '1' || c == 'X').count();
        if weight < 3 || weight > 5 {
            error!("Spaced seed weight is {}. Must be between 3 and 5.", weight);
            std::process::exit(1);
        }
    }

    // Output path check: Ensure the parent directory exists
    let path = Path::new(&args.output);
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() && !parent.exists() {
            error!("Output directory does not exist: {:?}", parent);
            std::process::exit(1);
        }
    }
}

pub fn define_mapping(args: DefineMappingArgs) {

    check_args(&args);

    info!("Mapping!");
    return;
}