mod args;

use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::sync::mpsc;
use std::thread;
use rand::prelude::IndexedRandom;
use rand::rng; 

use metrics::{get_ca_atoms, radius_of_gyration, bounding_box_volume, contact_order, plddt_statistics};
use pdb_io::parse_pdb;
use args::parse_arguments;


pub fn process_pdb_file(file_path: &str) -> (String, f32, f32, f32, f32, f32, f32, f32, usize) {
    // Parse the PDB file.
    let pdb = parse_pdb(file_path).expect("Failed to parse PDB file");

    let ca_atoms = get_ca_atoms(&pdb);

    let rg = radius_of_gyration(&ca_atoms);
    let vol = bounding_box_volume(&ca_atoms);
    let co = contact_order(&ca_atoms);
    let (mean_plddt, plddt_50, plddt_70, plddt_90) = plddt_statistics(&ca_atoms);

    let file_stem = Path::new(file_path)
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("unknown")
        .to_string();

    let length = ca_atoms.len();

    (file_stem, rg, vol, co, mean_plddt, plddt_50, plddt_70, plddt_90, length)
}

fn append_to_file(file_path: &str, output: &str) {
    let mut file = OpenOptions::new()
        .append(true)   // Enable appending
        .create(true)   // Create file if it doesn't exist
        .open(file_path)
        .expect("Failed to open file for appending");

    writeln!(file, "{}", output).expect("Failed to write to file");
}

fn main() {

    let config = parse_arguments();

    println!(
        "Using {} CPUs\nReading from: {}\nSaving results to: {}",
        config.num_cpus, config.pdb_dir, config.output_file
    );

    fs::write(&config.output_file, "ID;Gyration_Radius;Box_Volume;Contact_Order;mean_pLDDT;pLDDT_50;pLDDT_70;pLDDT_90;seq_len\n")
        .expect("Failed to write CSV header");

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

    let files_to_process: Vec<String> = if let Some(n) = config.subset {
        let mut _rng = rng(); // Ensure mutable RNG
        pdb_files
            .choose_multiple(&mut _rng, n.min(pdb_files.len())) // Trait method explicitly used
            .cloned() // Convert references to owned Strings
            .collect() // Collect into a Vec<String>
    } else {
        pdb_files.clone() // Ensure ownership consistency
    };


    // Create a channel for inter-process communication.
    let (tx, rx) = mpsc::channel();

    // Spawn a process for each file.
    let files_clone = files_to_process.clone();
    for file in files_clone {
        let tx = tx.clone();
        let output_file = config.output_file.clone();

        thread::spawn(move || {
            let result = process_pdb_file(&file);
            let output = format!(
                "{};{:.4};{:.4};{:.4};{:.4};{:.4};{:.4};{:.4};{}",
                result.0, result.1, result.2, result.3, result.4, result.5, result.6, result.7, result.8
            );

            append_to_file(&output_file, &output);

            tx.send(()).expect("Failed to send completion signal");
        });
    }

    // Wait for all processes to finish.
    for _ in 0..files_to_process.len() {
        rx.recv().expect("Failed to receive completion signal");
    }

    println!("Processing complete. Results saved to {}", &config.output_file);
}