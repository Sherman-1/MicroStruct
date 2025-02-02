use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug)]
pub struct AtomCoordinate {
    pub atom_name: String,
    pub residue_seq: i32,
    pub residue_name: String,
    pub bfactor: f32,
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[derive(Debug)]
pub struct ParsedPDB {
    pub atoms: Vec<AtomCoordinate>,
}


pub fn parse_pdb<P: AsRef<Path>>(pdb_path: P) -> Result<ParsedPDB, Box<dyn Error>> {

    /*
    1-4   "ATOM"                          left   character
    7-11  Atom serial number               right  integer
    13-16 Atom name                        left   * character
    17    Alternate location indicator     character
    18-20 Residue name                     right  character
    22    Chain identifier                 character
    23-26 Residue sequence number          right  integer
    27    Code for insertions of residues  character
    31-38 X orthogonal Angstrom coordinate right  floating
    39-46 Y orthogonal Angstrom coordinate right  floating
    47-54 Z orthogonal Angstrom coordinate right  floating
    55-60 Occupancy                        right  floating
    61-66 Temperature factor               right  floating
    73-76 Segment identifier (optional)    left   character
    77-78 Element symbol                   right  character
    79-80 Charge (optional)                character

    Taken from https://www.biostat.jhsph.edu/~iruczins/teaching/260.655/links/pdbformat.pdf
    */


    let file = File::open(pdb_path)?;
    let reader = BufReader::new(file);

    // Pre-allocate a vector for performance.
    let mut atoms = Vec::with_capacity(10_000);

    for line_result in reader.lines() {
        let line = line_result?;

        // Only process lines that start with "ATOM ".
        if line.starts_with("ATOM ") {
            // Parse the atom name (columns 13-16)
            let atom_name = line.get(12..16)
                .map(str::trim)
                .ok_or_else(|| "Missing atom name".to_string())?
                .to_string();

            // Parse the residue name (columns 18-20)
            let residue_name = line.get(17..20)
                .map(str::trim)
                .ok_or_else(|| "Missing residue name".to_string())?
                .to_string();

            // Parse the residue sequence number (columns 23-26)
            let residue_seq = line.get(22..26)
                .map(str::trim)
                .ok_or_else(|| "Missing residue sequence".to_string())?
                .parse::<i32>()
                .map_err(|_| "Failed to parse residue sequence".to_string())?;

            // Parse the X coordinate (columns 31-38)
            let x = line.get(30..38)
                .map(str::trim)
                .ok_or_else(|| "Missing x coordinate".to_string())?
                .parse::<f32>()
                .map_err(|_| "Failed to parse x coordinate".to_string())?;

            // Parse the Y coordinate (columns 39-46)
            let y = line.get(38..46)
                .map(str::trim)
                .ok_or_else(|| "Missing y coordinate".to_string())?
                .parse::<f32>()
                .map_err(|_| "Failed to parse y coordinate".to_string())?;

            // Parse the Z coordinate (columns 47-54)
            let z = line.get(46..54)
                .map(str::trim)
                .ok_or_else(|| "Missing z coordinate".to_string())?
                .parse::<f32>()
                .map_err(|_| "Failed to parse z coordinate".to_string())?;

            // Parse the B-factor (columns 61-66)
            let bfactor = line.get(60..66)
                .map(str::trim)
                .ok_or_else(|| "Missing B-factor".to_string())?
                .parse::<f32>()
                .map_err(|_| "Failed to parse B-factor".to_string())?;

            atoms.push(AtomCoordinate {
                atom_name,
                residue_seq,
                residue_name,
                bfactor,
                x,
                y,
                z,
            });
        }
    }

    Ok(ParsedPDB { atoms })
}  