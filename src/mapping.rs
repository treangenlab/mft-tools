use crate::cli::*;
use crate::consts::*;

use log::*;
use std::fs::File;
use std::path::Path;
use std::io::{BufWriter, Write};
use std::collections::HashMap;
use itertools::Itertools;

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
        if weight < 2 || weight > 3 {
            error!("Spaced seed weight is {}. Must be between 2 and 3.", weight);
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

    info!("Creating mapping for k={}", args.k);


    let file = File::create(&args.output).expect("Unable to create output file");
    let mut writer = BufWriter::new(file);

    if let Some(ref seed) = args.spaced_seed {
        generate_spaced_mapping(&args, seed, &mut writer);
    } else {
        // Placeholder for fixed alphabet size logic which is stil tbd
        println!("Alphabet size logic not yet implemented.");
    }
    return;
}

fn generate_spaced_mapping(args: &DefineMappingArgs, seed: &str, writer: &mut BufWriter<File>) {
    let bases = ['A', 'C', 'G', 'T'];
    let mut pattern_to_char = HashMap::new();
    let mut next_char_idx = 0;
    
    // Collect into a Vec once for faster indexing
    let alphabet: Vec<char> = MAP_ALPHABET.chars().collect();

    // Enumerate every possible k-mer
    for kmer_vec in std::iter::repeat(bases.iter()).take(args.k).multi_cartesian_product() {
        let kmer: String = kmer_vec.iter().cloned().collect();

        // Extract key based on spaced seed
        let key: String = kmer.chars()
            .enumerate()
            .filter(|(i, _)| {
                // Get seed char safely; check_args already verified length
                seed.chars().nth(*i).map_or(false, |c| c == '1' || c == 'X')
            })
            .map(|(_, c)| c)
            .collect();

        // Map key to char
        let mapped_char = *pattern_to_char.entry(key).or_insert_with(|| {
            if next_char_idx >= alphabet.len() {
                error!("Exhausted alphabet! Spaced seed creates more than {} patterns.", alphabet.len());
                std::process::exit(1);
            }
            let c = alphabet[next_char_idx];
            next_char_idx += 1;
            c
        });

        if let Err(e) = writeln!(writer, "{}\t{}", kmer, mapped_char) {
            error!("Failed to write to file: {}", e);
            std::process::exit(1);
        }

    }

    info!(
        "Wrote mapping of kmers (k={}) with pattern {} to reduced alphabet with {} characters", 
        args.k, seed, next_char_idx
    );
}