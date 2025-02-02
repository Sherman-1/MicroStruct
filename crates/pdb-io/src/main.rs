use std::error::Error;
use pdb_io::parse_pdb; 

fn main() -> Result<(), Box<dyn Error>> {
   
    // Very very dégueulasse
    let parsed_pdb = parse_pdb("/store/EQUIPES/BIM/MEMBERS/simon.herman/MicroStruct/test.pdb")?;
    
    println!("Parsed {} alpha-carbons", parsed_pdb.ca_coords.len());
    Ok(())
}
