use crate::cli::*;
use crate::utils::*;
use log::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

fn check_args(args: &EvaluateArgs) {
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

    if !Path::new(&args.mapping_table).exists() {
        error!("Mapping table file does not exist: {}", args.mapping_table);
        std::process::exit(1);
    }

    if !Path::new(&args.ordering).exists() {
        error!("Order file does not exist: {}", args.ordering);
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

pub fn evaluate(args: EvaluateArgs) {
    check_args(&args);
    let mut args = args;
    args.output = ensure_extension(&args.output, "txt");

    info!(
        "Evaluating masking rate: k={}, w={}, samples={}",
        args.k, args.w, args.n_samples
    );

    let mapping = load_mapping(&args.mapping_table);
    let order = load_order(&args.ordering);
    let order_map = get_order_map(&order);

    let (global_rate, entropy, mutation_rates) =
        calculate_masking_rate(&order_map, &mapping, args.k, args.w, Some(args.n_samples));

    info!("Global masking rate: {:.4}%", global_rate * 100.0);
    info!("Entropy: {:.4}", entropy);

    let out_file = File::create(&args.output).unwrap_or_else(|_| {
        error!("Failed to create output file: {}", args.output);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(out_file);

    writeln!(writer, "global_masking_rate\t{:.6}", global_rate).unwrap();
    writeln!(writer, "entropy\t{:.6}", entropy).unwrap();

    let mut transitions: Vec<(String, f64)> = mutation_rates.into_iter().collect();
    transitions.sort_by(|a, b| a.0.cmp(&b.0));
    for (label, rate) in transitions {
        writeln!(writer, "{}\t{:.6}", label, rate).unwrap();
    }

    writer.flush().ok();
    info!("Evaluation results saved to {}", args.output);
}
