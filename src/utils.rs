use itertools::Itertools;
use log::*;
use rand::RngExt;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

pub fn load_mapping(path: &str) -> FxHashMap<String, char> {
    let file = File::open(path).unwrap_or_else(|_| {
        error!("Could not open mapping file: {}", path);
        std::process::exit(1);
    });

    let reader = BufReader::new(file);
    let mut mapping = FxHashMap::default();

    for (index, line) in reader.lines().enumerate() {
        let line = line.unwrap_or_else(|e| {
            error!("Error reading line {} in mapping file: {}", index + 1, e);
            std::process::exit(1);
        });

        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() == 2 {
            let kmer = parts[0].to_string();
            let mapped_char = parts[1].chars().next().unwrap_or('X');
            mapping.insert(kmer, mapped_char);
        } else {
            warn!("Skipping malformed line {}: {}", index + 1, line);
        }
    }

    info!("Loaded mapping table with {} entries", mapping.len());
    mapping
}

pub fn get_order_map(order: &[String]) -> FxHashMap<String, usize> {
    order
        .iter()
        .enumerate()
        .map(|(i, k)| (k.clone(), i))
        .collect()
}

pub fn generate_all_kmers(k: usize) -> Vec<String> {
    let bases = [b'A', b'C', b'G', b'T'];
    let n = 1 << (2 * k);
    (0..n)
        .map(|i| {
            let mut kmer = vec![0u8; k];
            for j in 0..k {
                kmer[k - 1 - j] = bases[(i >> (2 * j)) & 3];
            }
            String::from_utf8(kmer).unwrap()
        })
        .collect()
}

pub fn load_order(path: &str) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|_| {
            error!("Could not open order file: {}", path);
            std::process::exit(1);
        })
        .lines()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

pub fn save_order(path: &str, order: &[String]) {
    let file = File::create(path).unwrap_or_else(|_| {
        error!("Unable to create output file: {}", path);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(file);
    for kmer in order {
        writeln!(writer, "{}", kmer).expect("Failed to write k-mer to ordering file");
    }
    info!("Ordering saved to {}", path);
}

fn minimizer(
    window: &[u8],
    order_map: &FxHashMap<String, usize>,
    mapping: &FxHashMap<String, char>,
    k: usize,
) -> char {
    let mut best_rank = usize::MAX;
    let mut best_kmer: Option<String> = None;

    for i in 0..=(window.len() - k) {
        let kmer = String::from_utf8_lossy(&window[i..i + k]).to_string();
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

    let total_possible_windows = 4_usize.pow(w as u32);
    let num_to_eval = match n_samples {
        Some(n) => n.min(total_possible_windows),
        None => total_possible_windows,
    };

    let mut mutation_stats: FxHashMap<(u8, u8), (usize, usize)> = FxHashMap::default();
    for &b1 in &bases {
        for &b2 in &bases {
            if b1 != b2 {
                mutation_stats.insert((b1, b2), (0, 0));
            }
        }
    }

    let mut rng = rand::rng();

    let windows_to_eval: Vec<Vec<u8>> = if num_to_eval == total_possible_windows && w <= 10 {
        std::iter::repeat(bases.iter().copied())
            .take(w)
            .multi_cartesian_product()
            .collect()
    } else {
        let mut seen = FxHashSet::with_capacity_and_hasher(num_to_eval, Default::default());
        while seen.len() < num_to_eval {
            let window: Vec<u8> = (0..w).map(|_| bases[rng.random_range(0..4)]).collect();
            seen.insert(window);
        }
        seen.into_iter().collect()
    };

    for window in windows_to_eval {
        let original_res = minimizer(&window, order_map, mapping, k);
        *transformation_counts.entry(original_res).or_insert(0) += 1;

        for i in 0..w {
            let orig_base = window[i];
            for &mutation_base in &bases {
                if mutation_base == orig_base {
                    continue;
                }

                let mut mut_window = window.clone();
                mut_window[i] = mutation_base;
                let mut_res = minimizer(&mut_window, order_map, mapping, k);

                total_trials += 1;

                if let Some(stats) = mutation_stats.get_mut(&(orig_base, mutation_base)) {
                    stats.0 += 1;
                    if mut_res != original_res {
                        flicker_events += 1;
                        stats.1 += 1;
                    }
                }
            }
        }
    }

    let global_rate = 1.0 - (flicker_events as f64 / total_trials as f64);

    let mut entropy = 0.0;
    for &count in transformation_counts.values() {
        let p = count as f64 / num_to_eval as f64;
        if p > 0.0 {
            entropy -= p * p.log2();
        }
    }

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

pub fn check_fasta(file: &str) -> bool {
    if file.ends_with(".fa")
        || file.ends_with(".fasta")
        || file.ends_with(".fa.gz")
        || file.ends_with("fasta.gz")
        || file.ends_with("fna")
        || file.ends_with("fna.gz")
    {
        return true;
    }
    return false;
}
