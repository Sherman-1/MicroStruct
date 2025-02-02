mod args;

use metrics::{radius_of_gyration, bounding_box_volume, contact_order, plddt_statistics, get_ca_atoms, calc_sasa_from_parsed_pdb};
use rayon::prelude::*;
use std::path::Path;
use std::fs;
use rand::prelude::*; 
use rayon::ThreadPoolBuilder;
use args::parse_arguments;
use std::fs::OpenOptions;
use std::io::Write;
use pdb_io::{parse_pdb};


pub fn process_pdb_file(file_path: &str) -> (String, f32, f32, f32, f32, f32, f32, f32, usize) {
    // Parse the PDB file.
    let pdb = parse_pdb(file_path).expect("Failed to parse PDB file");

    let ca_atoms = get_ca_atoms(&pdb);

    let rg = radius_of_gyration(&ca_atoms);
    let vol = bounding_box_volume(&ca_atoms);
    let co = contact_order(&ca_atoms);
    let (mean_plddt, plddt_50, plddt_70, plddt_90) = plddt_statistics(&ca_atoms);

    // Extract the file stem for naming purposes.
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

    fs::write(&config.output_file, "ID;Gyration_Radius;Box_Volume;Contact_Order;mean_pLDDT;pLDDT_50;pLDDT_70;pLDDT_90;seq_len\n").expect("Failed to write CSV header");

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

    let files_to_process = if let Some(n) = config.subset {
        let mut _rng = rand::rng();
        pdb_files
            .choose_multiple(&mut _rng, n.min(pdb_files.len()))
            .cloned()
            .collect()
    } else {
        pdb_files
    };
        
    let results: Vec<_> = pool.install(|| {
        files_to_process
            .par_iter()
            .map(|file| process_pdb_file(file))
            .collect()
    });
        
    let output: String = results
        .iter()
        .map(|(id, rg, vol, co, pmean, p50, p70, p90, len)| {
            format!("{};{:.4};{:.4};{:.4};{:.4};{:.4};{:.4};{:.4};{}", id, rg, vol, co, pmean, p50, p70, p90, len)
        })
        .collect::<Vec<String>>() 
        .join("\n"); 

    
    append_to_file(&config.output_file, &output);

    println!("Processing complete. Results saved to {}", &config.output_file);
}