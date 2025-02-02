use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug)]
pub struct CaCoordinate {
    pub residue_seq: i32,
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[derive(Debug)]
pub struct ParsedPDB {
    pub ca_coords: Vec<CaCoordinate>,
}

/// Parses a PDB file, extracting only (x, y, z) coordinates for Cα (alpha carbon) atoms.
/// Skips lines with missing or invalid coordinate fields.
pub fn parse_pdb<P: AsRef<Path>>(pdb_path: P) -> Result<ParsedPDB, Box<dyn Error>> {
    let file = File::open(pdb_path)?;
    let reader = BufReader::new(file);

    let mut ca_coords = Vec::with_capacity(10_000);

    for line_result in reader.lines() {
        let line = line_result?;

        if line.starts_with("ATOM ") {
            if let Some("CA") = line.get(12..16).map(str::trim) {
                let residue_seq = line.get(22..26)
                    .map(str::trim)
                    .ok_or_else(|| "Missing residue sequence".to_string())?
                    .parse::<i32>()
                    .map_err(|_| "Failed to parse residue sequence".to_string())?;

                let x = line.get(30..38)
                    .map(str::trim)
                    .ok_or_else(|| "Missing x coordinate".to_string())?
                    .parse::<f32>()
                    .map_err(|_| "Failed to parse x coordinate".to_string())?;

                let y = line.get(38..46)
                    .map(str::trim)
                    .ok_or_else(|| "Missing y coordinate".to_string())?
                    .parse::<f32>()
                    .map_err(|_| "Failed to parse y coordinate".to_string())?;

                let z = line.get(46..54)
                    .map(str::trim)
                    .ok_or_else(|| "Missing z coordinate".to_string())?
                    .parse::<f32>()
                    .map_err(|_| "Failed to parse z coordinate".to_string())?;

                ca_coords.push(CaCoordinate {
                    residue_seq,
                    x,
                    y,
                    z,
                });
            }
        }
    }

    Ok(ParsedPDB { ca_coords })
}
