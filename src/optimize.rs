use crate::cli::*;
use crate::utils::*;
use log::*;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use rustc_hash::{FxHashMap, FxHashSet};

use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::RngExt;
use itertools::Itertools;

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
        error!("Invalid window size: {}. Window size (w) must be strictly greater than k-mer size (k={}).", args.w, args.k);
        std::process::exit(1);
    }

    if args.num_iterations < 1 || args.num_iterations > 10000 {
        error!("Number of iterations {} is out of bounds. Must be between 500 and 10000.", args.num_iterations);
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

    info!("Starting hill-climbing optimization");
    info!("Parameters: k={}, w={}, iterations={}", args.k, args.w, args.num_iterations);

    let mapping = load_mapping(&args.mapping_table);
    let (mut current_order, _seed) = initialize_order(args.k, &args.base_order);

    let mut current_map = get_order_map(&current_order);

    let (mut best_rate, mut curr_entropy, mut best_transition_rates) = calculate_masking_rate(
        &current_map, &mapping, args.k, args.w, Some(args.num_samples)
    );

    info!("Initial Masking Rate: {:.4}%", best_rate * 100.0);

    let mut rng = rand::rng();

    for i in 1..=args.num_iterations {
        let len = current_order.len();
        let idx1 = rng.random_range(0..len);
        let idx2 = rng.random_range(0..len);
        
        current_order.swap(idx1, idx2);
        let test_map = get_order_map(&current_order);

        let (test_rate, test_entropy, test_transition_rates) = calculate_masking_rate(
            &test_map, &mapping, args.k, args.w, Some(args.num_samples)
        );

        if test_rate > best_rate {
            let diff = test_rate - best_rate;
            best_rate = test_rate;
            best_transition_rates = test_transition_rates;
            curr_entropy = test_entropy;
            info!(" Trial {:4}: Improvement -- New Masking Rate: {:.4}% (+{:.4}%) -- New Entropy: {:.2}", i, best_rate * 100.0, diff * 100.0, curr_entropy);
        } else {
            current_order.swap(idx1, idx2);
        }
    }

    save_order(&args.output, &current_order);
}

pub fn initialize_order(
    k: usize, 
    base_order_path: &Option<String>
) -> (Vec<String>, u64) {
    if let Some(path) = base_order_path {
        info!("Loading initial ordering from: {}", path);
        let content = std::fs::read_to_string(path).expect("Failed to read base order file");
        let order: Vec<String> = content.lines().map(|s| s.trim().to_string()).collect();
        (order, 0) 
    } else {
        let mut kmers = generate_all_kmers(k);
        let seed = rand::random::<u64>();
        info!("No base order provided. Initializing random order with seed: {}", seed);
        let mut rng = StdRng::seed_from_u64(seed);
        kmers.shuffle(&mut rng);
        (kmers, seed)
    }
}

fn generate_all_kmers(k: usize) -> Vec<String> {
    let bases = ['A', 'C', 'G', 'T'];
    std::iter::repeat(bases.iter())
        .take(k)
        .multi_cartesian_product()
        .map(|v| v.into_iter().collect())
        .collect()
}

pub fn calculate_masking_rate(
    order_map: &FxHashMap<String, usize>,
    mapping: &FxHashMap<String, char>,
    k: usize,
    w: usize,
    n_samples: Option<usize>,
) -> (f64, f64, FxHashMap<String, f64>) {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut total_trials = 0;
    let mut flicker_events = 0;

    let mut transformation_counts: FxHashMap<char, usize> = FxHashMap::default();

    // Calculate the number of windows for a given k size
    let total_possible_windows = 4_usize.pow(w as u32);
    let num_to_eval = match n_samples {
        Some(n) => n.min(total_possible_windows),
        None => total_possible_windows,
    };

    // Initialize tracking for the 12 transition types
    // Maps (from_base, to_base) -> (total_trials, flicker_events)
    let mut mutation_stats: FxHashMap<(u8, u8), (usize, usize)> = FxHashMap::default();
    for &b1 in &bases {
        for &b2 in &bases {
            if b1 != b2 {
                mutation_stats.insert((b1, b2), (0, 0));
            }
        }
    }

    let mut rng = rand::rng();

    // generate the windows to evaluate using hash-set approach
    let windows_to_eval: Vec<Vec<u8>> = if num_to_eval == total_possible_windows && w <= 10 {
        // Safe to generate all exactly
        std::iter::repeat(bases.iter().copied())
            .take(w)
            .multi_cartesian_product()
            .collect()
    } else {
        // Subsample using a HashSet to guarantee uniqueness without generating the whole universe
        let mut seen = FxHashSet::with_capacity_and_hasher(num_to_eval, Default::default());
        while seen.len() < num_to_eval {
            let window: Vec<u8> = (0..w).map(|_| bases[rng.random_range(0..4)]).collect();
            seen.insert(window);
        }
        seen.into_iter().collect()
    };

    // Evaluate masking rates over those windows
    for window in windows_to_eval {
        let original_res = minimizer(&window, order_map, mapping, k);
        *transformation_counts.entry(original_res).or_insert(0) += 1;

        for i in 0..w {
            let orig_base = window[i];
            for &mutation_base in &bases {
                if mutation_base == orig_base { continue; }

                let mut mut_window = window.clone();
                mut_window[i] = mutation_base;

                let mut_res = minimizer(&mut_window, order_map, mapping, k);

                // info!("{} {} {} {}", String::from_utf8_lossy(&window), original_res, String::from_utf8_lossy(&mut_window), mut_res);
                total_trials += 1;
                
                // Track the specific mutation type
                if let Some(stats) = mutation_stats.get_mut(&(orig_base, mutation_base)) {
                    stats.0 += 1; // Increment total trials for this specific transition
                    
                    if mut_res != original_res {
                        flicker_events += 1;
                        stats.1 += 1; // Increment flicker events for this specific transition
                    }
                }
            }
        }
    }

    // 5. Final Calculations
    let global_rate = 1.0 - (flicker_events as f64 / total_trials as f64);

    let mut entropy = 0.0;
    for &count in transformation_counts.values() {
        let p = count as f64 / num_to_eval as f64;
        if p > 0.0 {
            entropy -= p * p.log2();
        }
    }

    // Compile the 12 transition rates into a friendly format
    let mut mutation_rates = FxHashMap::default();
    for ((from, to), (trials, flickers)) in mutation_stats {
        let label = format!("{}->{}", from as char, to as char);
        let rate = if trials > 0 {
            1.0 - (flickers as f64 / trials as f64)
        } else {
            0.0
        };
        mutation_rates.insert(label, rate);
    }

    (global_rate, entropy, mutation_rates)
}

fn minimizer(
    window: &[u8],
    order_map: &FxHashMap<String, usize>,
    mapping: &FxHashMap<String, char>,
    k: usize
) -> char {
    let mut best_rank = usize::MAX;
    let mut best_kmer: Option<String> = None;

    for i in 0..=(window.len() - k) {
        let kmer = String::from_utf8_lossy(&window[i..i+k]).to_string();
        
        if let Some(&rank) = order_map.get(&kmer) {
            if rank < best_rank {
                best_rank = rank;
                best_kmer = Some(kmer);
            }
        }
    }

    match best_kmer {
        Some(k) => *mapping.get(&k).unwrap_or(&'X'),
        None => 'X',
    }
}

pub fn save_order(path: &str, order: &[String]) {
    let file = File::create(path).expect("Unable to create output file");
    let mut writer = BufWriter::new(file);

    for kmer in order {
        writeln!(writer, "{}", kmer).expect("Failed to write k-mer to ordering file");
    }
    info!("Optimized ordering saved to {}", path);
}