use crate::cli::*;
use crate::utils::*;
use log::*;
use needletail::parse_fastx_file;
use rand::seq::SliceRandom;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rustc_hash::FxHashMap;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn check_args(args: &MinFrameTransitionArgs) {
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

    if let Some(ref order_path) = args.order {
        if !Path::new(order_path).exists() {
            error!("Order file provided does not exist: {}", order_path);
            std::process::exit(1);
        }
    }

    for fasta_file in &args.genomes {
        if !check_fasta(&fasta_file) {
            error!(
                "{} does not appear to be a fasta file (must be .fa(.gz)/.fasta(.gz)/.fna(.gz))",
                &fasta_file
            );
            std::process::exit(1)
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

/// Converts a DNA string to a 2-bit encoded integer.
fn kmer_to_idx(kmer: &str) -> usize {
    kmer.as_bytes().iter().fold(0, |acc, &b| {
        let val = match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0, // Ns and invalids map to A (TODO: handle this case)
        };
        (acc << 2) | val
    })
}

/// Computes the reverse complement of a 2-bit encoded k-mer integer.
fn rev_comp_kmer(mut kmer: usize, k: usize) -> usize {
    let mut rev = 0;
    for _ in 0..k {
        // Extract last 2 bits, XOR with 3 to complement (A<->T, C<->G), 
        // then push into the new integer.
        rev = (rev << 2) | ((kmer & 3) ^ 3);
        kmer >>= 2;
    }
    rev
}


/// Builds the flat arrays required for execution of simd-minimizers
pub fn prepare_luts(
    order_map: &FxHashMap<String, usize>,
    mapping: &FxHashMap<String, char>,
    k: usize,
) -> (Vec<u32>, Vec<char>) {
    let lut_size = 1 << (2 * k);
    let mut order_lut = vec![u32::MAX; lut_size];
    let mut rank_to_char_lut = vec!['?'; order_map.len()];

    for (kmer, &rank) in order_map {
        let idx = kmer_to_idx(kmer);

        if idx < lut_size {
            order_lut[idx] = rank as u32;
        }

        if let Some(&mapped_char) = mapping.get(kmer) {
            rank_to_char_lut[rank] = mapped_char;
        }
    }

    (order_lut, rank_to_char_lut)
}

fn generate_all_kmers(k: usize) -> Vec<String> {
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

pub fn mft(args: MinFrameTransitionArgs) {
    check_args(&args);

    let raw_mapping = load_mapping(&args.mapping_table);

    let raw_order_map = if let Some(ref order_path) = args.order {
        let order_lines = std::fs::read_to_string(order_path)
            .expect("Failed to read order file")
            .lines()
            .map(String::from)
            .collect::<Vec<String>>();
        get_order_map(&order_lines)
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
        let mut kmers = generate_all_kmers(args.k);
        kmers.shuffle(&mut rng);
        get_order_map(&kmers)
    };

    let (order_lut, rank_to_char_lut) = prepare_luts(&raw_order_map, &raw_mapping, args.k);

    let out_file = File::create(&args.output).expect("Failed to create output file");
    let mut buf_writer = BufWriter::new(out_file);

    let window_kmers = args.w - args.k + 1;
    let k_mask = (1 << (2 * args.k)) - 1;

    for genome_path in &args.genomes {
        let mut reader = parse_fastx_file(genome_path).unwrap_or_else(|_| {
            panic!("Failed to open or parse FASTA/FASTQ file: {}", genome_path);
        });

        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid record");
            let seq = seqrec.seq(); // This is the Cow<[u8]>

            if seq.len() < args.k { continue; }

            let mut ranks = Vec::with_capacity(seq.len() - args.k + 1);
            let mut current_kmer = 0usize;

            // Single pass: Calculate canonical rank (min of fwd and rev k-mer)
            for (i, &b) in seq.iter().enumerate() {
                let val = match b {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    b'T' | b't' => 3,
                    _ => 0,
                };
                current_kmer = ((current_kmer << 2) | val) & k_mask;

                if i >= args.k - 1 {
                    // Forward rank
                    let fwd_rank = *order_lut.get(current_kmer).unwrap_or(&u32::MAX);
                    
                    let mut rev_kmer = 0usize;
                    let mut tmp = current_kmer;
                    for _ in 0..args.k {
                        rev_kmer = (rev_kmer << 2) | ((tmp & 3) ^ 3);
                        tmp >>= 2;
                    }
                    let rev_rank = *order_lut.get(rev_kmer).unwrap_or(&u32::MAX);

                    // Canonical Rank is the minimum of both strands at this position
                    ranks.push(fwd_rank.min(rev_rank));
                }
            }

            // --- Sliding Window Logic ---
            let mut transformed_sequence = String::with_capacity(ranks.len());
            let mut window = VecDeque::with_capacity(window_kmers);
            let n = ranks.len();
            let mut r = 0;

            for l in 0..n {
                while r < n && r < l + window_kmers {
                    while let Some(&idx) = window.back() {
                        if ranks[idx] > ranks[r] { window.pop_back(); } 
                        else { break; }
                    }
                    window.push_back(r);
                    r += 1;
                }
                while let Some(&idx) = window.front() {
                    if idx < l { window.pop_front(); } 
                    else { break; }
                }

                let min_idx = *window.front().unwrap();
                let rank = ranks[min_idx];
                let c = rank_to_char_lut.get(rank as usize).copied().unwrap_or('?');
                transformed_sequence.push(c);
            }

            writeln!(buf_writer, ">{}", String::from_utf8_lossy(seqrec.id())).unwrap();
            for chunk in transformed_sequence.as_bytes().chunks(60) {
                buf_writer.write_all(chunk).expect("Failed to write sequence chunk");
                buf_writer.write_all(b"\n").expect("Failed to write newline");
            }
        }
    }
    buf_writer.flush().ok();
}
