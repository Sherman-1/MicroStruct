use std::fs;
use std::sync::mpsc;
use std::thread;

mod main {
    include!("../src/main.rs");
}
use main::process_pdb_file;

#[test]
fn test() {
    // Read all PDB files from the data directory
    let pdb_files: Vec<String> = fs::read_dir("../../data/")
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

    let (tx, rx) = mpsc::channel();
    let num_files = pdb_files.len();

    for file in pdb_files {
        let tx = tx.clone();
        
        thread::spawn(move || {
            let result = process_pdb_file(&file);
            let output = format!(
                "{};{:.4};{:.4};{:.4};{:.4};{:.4};{:.4};{:.4};{}",
                result.0, result.1, result.2, result.3, result.4, result.5, result.6, result.7, result.8
            );

            tx.send(output).expect("Failed to send result");
        });
    }

    let mut results = Vec::new();
    for _ in 0..num_files {
        let output = rx.recv().expect("Failed to receive result");
        results.push(output);
    }
    
    for output in results {
        println!("{}", output);
    }
}

