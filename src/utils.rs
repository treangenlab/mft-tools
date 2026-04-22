use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use log::*;


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
    order.iter().enumerate().map(|(i, k)| (k.clone(), i)).collect()
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