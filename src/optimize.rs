use crate::cli::*;
use crate::utils::*;
use log::*;
use std::path::Path;

use rand::RngExt;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;

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
    info!(
        "Parameters: k={}, w={}, iterations={}",
        args.k, args.w, args.num_iterations
    );

    let mapping = load_mapping(&args.mapping_table);
    let (mut current_order, _seed) = initialize_order(args.k, &args.base_order);

    let mut current_map = get_order_map(&current_order);

    let (mut best_rate, mut curr_entropy, mut best_transition_rates) = calculate_masking_rate(
        &current_map,
        &mapping,
        args.k,
        args.w,
        Some(args.num_samples),
    );

    info!("Initial Masking Rate: {:.4}%", best_rate * 100.0);

    let mut rng = rand::rng();

    for i in 1..=args.num_iterations {
        let len = current_order.len();
        let idx1 = rng.random_range(0..len);
        let idx2 = rng.random_range(0..len);

        current_order.swap(idx1, idx2);
        let test_map = get_order_map(&current_order);

        let (test_rate, test_entropy, test_transition_rates) =
            calculate_masking_rate(&test_map, &mapping, args.k, args.w, Some(args.num_samples));

        if test_rate > best_rate {
            let diff = test_rate - best_rate;
            best_rate = test_rate;
            best_transition_rates = test_transition_rates;
            curr_entropy = test_entropy;
            info!(
                " Trial {:4}: Improvement -- New Masking Rate: {:.4}% (+{:.4}%) -- New Entropy: {:.2}",
                i,
                best_rate * 100.0,
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
