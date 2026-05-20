use crate::cli::*;
use crate::utils::*;
use log::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rand::RngExt;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;

pub fn load_transition_weights(path: &str) -> FxHashMap<String, f64> {
    let file = File::open(path).unwrap_or_else(|_| {
        error!("Could not open transition weights file: {}", path);
        std::process::exit(1);
    });

    let reader = BufReader::new(file);
    let mut weights = FxHashMap::default();

    for (index, line) in reader.lines().enumerate() {
        let line = line.unwrap_or_else(|e| {
            error!("Error reading line {} in weights file: {}", index + 1, e);
            std::process::exit(1);
        });
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 2 {
            error!(
                "Malformed line {} in weights file (expected transition<tab>weight): {}",
                index + 1,
                line
            );
            std::process::exit(1);
        }
        let transition = parts[0].trim().to_string();
        let weight: f64 = parts[1].trim().parse().unwrap_or_else(|_| {
            error!(
                "Invalid weight on line {} of weights file: {}",
                index + 1,
                parts[1]
            );
            std::process::exit(1);
        });
        if weight < 0.0 {
            error!("Negative weight on line {}: {}", index + 1, weight);
            std::process::exit(1);
        }
        weights.insert(transition, weight);
    }

    if weights.is_empty() {
        error!("Transition weights file is empty or has no valid entries: {}", path);
        std::process::exit(1);
    }

    let total: f64 = weights.values().sum();
    if (total - 1.0).abs() > 1e-4 {
        error!(
            "Transition weights must sum to 1.0 (current sum: {:.6}): {}",
            total, path
        );
        std::process::exit(1);
    }

    info!("Loaded {} transition weights from {}", weights.len(), path);
    weights
}

pub fn compute_weighted_score(
    transition_rates: &FxHashMap<String, f64>,
    weights: &FxHashMap<String, f64>,
) -> f64 {
    weights
        .iter()
        .filter_map(|(label, &weight)| transition_rates.get(label).map(|&rate| weight * rate))
        .sum()
}

pub fn check_args(args: &OptimizeOrderArgs) {
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

    if args.w <= args.k {
        error!(
            "Invalid window size: {}. Window size (w) must be strictly greater than k-mer size (k={}).",
            args.w, args.k
        );
        std::process::exit(1);
    }

    if args.num_iterations < 1 || args.num_iterations > 10000 {
        error!(
            "Number of iterations {} is out of bounds. Must be between 500 and 10000.",
            args.num_iterations
        );
        std::process::exit(1);
    }

    if !Path::new(&args.mapping_table).exists() {
        error!("Mapping table file does not exist: {}", args.mapping_table);
        std::process::exit(1);
    }

    if let Some(ref base_path) = args.base_order {
        if !Path::new(base_path).exists() {
            error!("Base order file provided does not exist: {}", base_path);
            std::process::exit(1);
        }
    }

    if let Some(ref weights_path) = args.transition_weights {
        if !Path::new(weights_path).exists() {
            error!("Transition weights file does not exist: {}", weights_path);
            std::process::exit(1);
        }
    }

    let out_path = Path::new(&args.output);
    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() && !parent.exists() {
            error!("Output directory does not exist: {:?}", parent);
            std::process::exit(1);
        }
    }
}

pub fn optimize(args: OptimizeOrderArgs) {
    check_args(&args);
    let mut args = args;
    args.output = ensure_extension(&args.output, "txt");

    let weights = args.transition_weights.as_deref().map(load_transition_weights);

    info!("Starting hill-climbing optimization");
    info!(
        "Parameters: k={}, w={}, iterations={}",
        args.k, args.w, args.num_iterations
    );
    if weights.is_some() {
        info!("Scoring mode: weighted transition masking rate");
    } else {
        info!("Scoring mode: global masking rate");
    }

    let mapping = load_mapping(&args.mapping_table);
    let (mut current_order, _seed) = initialize_order(args.k, &args.base_order);

    let mut current_map = get_order_map(&current_order);

    let (initial_global_rate, mut curr_entropy, initial_transition_rates) = calculate_masking_rate(
        &current_map,
        &mapping,
        args.k,
        args.w,
        Some(args.num_samples),
    );

    let mut best_score = match &weights {
        Some(w) => compute_weighted_score(&initial_transition_rates, w),
        None => initial_global_rate,
    };
    let mut best_transition_rates = initial_transition_rates;

    info!("Initial Score: {:.4}%", best_score * 100.0);

    let mut rng = rand::rng();

    for i in 1..=args.num_iterations {
        let len = current_order.len();
        let idx1 = rng.random_range(0..len);
        let idx2 = rng.random_range(0..len);

        current_order.swap(idx1, idx2);
        let test_map = get_order_map(&current_order);

        let (test_global_rate, test_entropy, test_transition_rates) =
            calculate_masking_rate(&test_map, &mapping, args.k, args.w, Some(args.num_samples));

        let test_score = match &weights {
            Some(w) => compute_weighted_score(&test_transition_rates, w),
            None => test_global_rate,
        };

        if test_score > best_score {
            let diff = test_score - best_score;
            best_score = test_score;
            best_transition_rates = test_transition_rates;
            curr_entropy = test_entropy;
            info!(
                " Trial {:4}: Improvement -- New Score: {:.4}% (+{:.4}%) -- New Entropy: {:.2}",
                i,
                best_score * 100.0,
                diff * 100.0,
                curr_entropy
            );
        } else {
            current_order.swap(idx1, idx2);
        }
    }

    save_order(&args.output, &current_order);
}

pub fn initialize_order(k: usize, base_order_path: &Option<String>) -> (Vec<String>, u64) {
    if let Some(path) = base_order_path {
        info!("Loading initial ordering from: {}", path);
        let order = load_order(path);
        (order, 0)
    } else {
        let mut kmers = generate_all_kmers(k);
        let seed = rand::random::<u64>();
        info!(
            "No base order provided. Initializing random order with seed: {}",
            seed
        );
        let mut rng = StdRng::seed_from_u64(seed);
        kmers.shuffle(&mut rng);
        (kmers, seed)
    }
}
