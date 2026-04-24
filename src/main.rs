use std::time::Instant;

use crate::consts::MFT_VERSION;

pub mod cli;
pub mod consts;
pub mod mapping;
pub mod mft;
pub mod optimize;
pub mod utils;

use cli::*;

fn main() {
    println!("Min-Frame Transformation v{}", MFT_VERSION);
    println!("Developed by Ryan Doughty (Rice University)");
    println!("Correspondence: rdd4@rice.edu, treangen@rice.edu\n");

    let start = Instant::now();

    let args = cli::parse_args();
    match args.mode {
        Mode::DefineMapping(define_mapping_args) => mapping::define_mapping(define_mapping_args),
        Mode::Optimize(optimize_args) => optimize::optimize(optimize_args),
        Mode::Transform(mft_args) => mft::mft(mft_args),
    }

    let end = Instant::now();
    eprintln!(
        "\nMFT v{} finished in {}s",
        MFT_VERSION,
        end.duration_since(start).as_secs_f32()
    );
}
