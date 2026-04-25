use crate::cli::*;
use crate::utils::*;
use log::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand::seq::SliceRandom;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn check_args(args: &DefineOrderArgs) {
    let output_level = if args.verbose {
        log::LevelFilter::Trace
    } else {
        log::LevelFilter::Info
    };

    let _ = simple_logger::SimpleLogger::new()
        .with_level(output_level)
        .init();

    if args.k < 3 || args.k > 15 {
        error!("Invalid k-mer size: {}. Must be between 3 and 15.", args.k);
        std::process::exit(1);
    }

    if args.seed.is_some() && args.alphabetical {
        error!("Cannot specify both --seed (-s) and --alphabetical (-a). Choose one.");
        std::process::exit(1);
    }

    let out_path = Path::new(&args.output);
    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() && !parent.exists() {
            error!("Output directory does not exist: {:?}", parent);
            std::process::exit(1);
        }
    }
}

pub fn define_order(args: DefineOrderArgs) {
    check_args(&args);

    let kmers = if args.alphabetical {
        info!("Using alphabetical k-mer ordering");
        generate_all_kmers(args.k)
    } else {
        let seed = match args.seed {
            Some(s) => s,
            None => {
                let s: u64 = rand::random();
                info!("No seed provided. Using randomly generated seed: {}", s);
                s
            }
        };
        let mut rng = StdRng::seed_from_u64(seed);
        let mut ordered = generate_all_kmers(args.k);
        ordered.shuffle(&mut rng);
        ordered
    };

    let out_file = File::create(&args.output).unwrap_or_else(|_| {
        error!("Failed to create output file: {}", args.output);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(out_file);
    for kmer in &kmers {
        writeln!(writer, "{}", kmer).expect("Failed to write k-mer to output");
    }
    writer.flush().ok();

    info!("Wrote {} k-mers to {}", kmers.len(), args.output);
}
