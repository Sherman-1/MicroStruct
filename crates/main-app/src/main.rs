mod args;

use pdb_io::parse_pdb;
use metrics::{radius_of_gyration, bounding_box_volume, contact_order};
use rayon::prelude::*;
use std::path::Path;
use std::fs;
use rayon::ThreadPoolBuilder;
use args::parse_arguments;


fn process_pdb_file(file_path: &str) -> (String, f64, f64, f64) {
    let pdb = parse_pdb(file_path).expect("Failed to parse PDB file");
    let rg = radius_of_gyration(&pdb) as f64;
    let vol = bounding_box_volume(&pdb) as f64;
    let co = contact_order(&pdb) as f64;

    // Extract the stem (filename without extension)
    let file_stem = Path::new(file_path)
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("unknown")
        .to_string();

    (file_stem, rg, vol, co)
}


fn main() {

    let config = parse_arguments();

    println!(
        "Using {} CPUs\nReading from: {}\nSaving results to: {}",
        config.num_cpus, config.pdb_dir, config.output_file
    );

    fs::write(&config.output_file, "ID;Radius;Volume;Contact\n").expect("Failed to write CSV header");

    let pool = ThreadPoolBuilder::new()
        .num_threads(config.num_cpus)
        .build()
        .expect("Failed to create thread pool");

    let pdb_files: Vec<String> = fs::read_dir(&config.pdb_dir)
        .expect("Failed to read PDB directory")
        .filter_map(|entry| {
            let path = entry.ok()?.path();
            if path.extension()?.to_str()? == "pdb" {
                Some(path.to_str()?.to_string())
            } else {
                None
            }
        })
        .collect();

    if pdb_files.is_empty() {
        eprintln!("No PDB files found in directory: {}", config.pdb_dir);
        return;
    };
        
    let results: Vec<(String, f64, f64, f64)> = pool.install(|| {
        pdb_files
            .par_iter()
            .map(|file_path| process_pdb_file(file_path))
            .collect()
    });
        
    let output: String = results.iter()
        .map(|(id, rg, vol, co)| format!("{};{:.4};{:.4};{:.4}\n", id, rg, vol, co))
        .collect();

    fs::write(&config.output_file, output).expect("Failed to write output file");

    println!("Processing complete. Results saved to {}", &config.output_file);
}