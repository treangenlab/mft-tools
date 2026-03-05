use std::time::Instant;

use crate::consts::SEEDING_VERSION;

pub mod consts;
pub mod cli;


fn main() {
    println!("MGA Seeding Analysis v{}", SEEDING_VERSION);
    println!("Developed by Ryan Doughty (Rice University)");
    println!("Correspondence: rdd4@rice.edu, treangen@rice.edu\n");

    let start = Instant::now();

    let args = cli::parse_args();

    let end = Instant::now();
    eprintln!("\nSeeding Analysis v{} finished in {}s", SEEDING_VERSION, end.duration_since(start).as_secs_f32());
}

