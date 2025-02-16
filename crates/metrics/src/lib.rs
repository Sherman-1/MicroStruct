//pub mod sasa;
use pdbtbx::*;

/*
pub fn radius_of_gyration(pdb: &PDB) -> f32 {
    // 1. Collect Cα coordinates from the PDB.
    let mut ca_coords = Vec::new();

    for model in pdb.models() {
        for chain in model.chains() {
            for residue in chain.residues() {
                for atom in residue.atoms() {
  
                    if atom.record_name() == "ATOM" && atom.atom_name().trim() == "CA" {
                        let xyz = atom.xyz();
                        
                        ca_coords.push((xyz[0] as f32, xyz[1] as f32, xyz[2] as f32));
                    }
                }
            }
        }
    }

    let n = ca_coords.len();
    if n == 0 {
        
        return 0.0;
    }

    let (sum_x, sum_y, sum_z) = ca_coords.iter().fold((0.0, 0.0, 0.0), |(sx, sy, sz), &(x, y, z)| {
        (sx + x, sy + y, sz + z)
    });
    let cx = sum_x / n as f32;
    let cy = sum_y / n as f32;
    let cz = sum_z / n as f32;

    let sum_dist_sq = ca_coords.iter().fold(0.0, |acc, &(x, y, z)| {
        let dx = x - cx;
        let dy = y - cy;
        let dz = z - cz;
        acc + (dx * dx + dy * dy + dz * dz)
    });

    
    (sum_dist_sq / n as f32).sqrt()
}


pub fn plddt_statistics(pdb: &PDB) -> (f32, f32, f32, f32) {

    let all_ca_bfactors: Vec<f32> = pdb
        .models()
        .flat_map(|model| model.chains())
        .flat_map(|chain| chain.residues())
        .flat_map(|res| res.atoms())
        .filter(|atom| atom.record_name() == "ATOM" && atom.atom_name().trim() == "CA")
        .map(|atom| atom.bfactor() as f32)
        .collect();

    let n = all_ca_bfactors.len();

    if n == 0 {
        return (0.0, 0.0, 0.0, 0.0);
    }

    let sum_plddt: f32 = all_ca_bfactors.iter().sum();

    let (count_50, count_70, count_90) = all_ca_bfactors.iter().fold((0usize, 0usize, 0usize), |(c50, c70, c90), &val| {
        (
            c50 + (val > 50.0) as usize,
            c70 + (val > 70.0) as usize,
            c90 + (val > 90.0) as usize,
        )
    });

    let mean_plddt = sum_plddt / n as f32;

    let fraction_50 = count_50 as f32 / n as f32 * 100.0;
    let fraction_70 = count_70 as f32 / n as f32 * 100.0;
    let fraction_90 = count_90 as f32 / n as f32 * 100.0;

    (mean_plddt, fraction_50, fraction_70, fraction_90)

}

pub fn bounding_box_volume(pdb: &PDB) -> f32 {
    
    let ca_coords: Vec<(f32, f32, f32)> = pdb
        .models()
        .flat_map(|model| model.chains())
        .flat_map(|chain| chain.residues())
        .flat_map(|residue| residue.atoms())
        
        .filter(|atom| atom.record_name() == "ATOM" && atom.atom_name().trim() == "CA")
        .map(|atom| {
            let (x, y, z) = atom.pos();
            (x as f32, y as f32, z as f32)
        })
        .collect();

    let n = ca_coords.len();
    if n == 0 {
        return 0.0;
    }

    
    let mut x_min = f32::MAX;
    let mut x_max = f32::MIN;
    let mut y_min = f32::MAX;
    let mut y_max = f32::MIN;
    let mut z_min = f32::MAX;
    let mut z_max = f32::MIN;

    
    for &(x, y, z) in &ca_coords {
        if x < x_min { x_min = x; }
        if x > x_max { x_max = x; }
        if y < y_min { y_min = y; }
        if y > y_max { y_max = y; }
        if z < z_min { z_min = z; }
        if z > z_max { z_max = z; }
    }

    let dx = x_max - x_min;
    let dy = y_max - y_min;
    let dz = z_max - z_min;
    dx * dy * dz
}

pub fn contact_order(pdb: &PDB) -> f32 {
    
    let ca_data: Vec<(f32, f32, f32, i32)> = pdb
        .models()
        .flat_map(|model| model.chains())
        .flat_map(|chain| chain.residues())
        .flat_map(|residue| {
            let rid = residue.serial_number(); // or residue.sequence_number()
            residue.atoms().filter_map(move |atom| {
                if atom.record_name() == "ATOM" && atom.atom_name().trim() == "CA" {
                    let (x, y, z) = atom.pos();
                    
                    Some((x as f32, y as f32, z as f32, rid))
                } else {
                    None
                }
            })
        })
        .collect();

    let l = ca_data.len();
    if l < 20 {
        
        return 0.0;
    }

    
    let mut sum_seqsep = 0.0;
    let mut n_contacts = 0.0;

    for i in 0..l {
        let (x_i, y_i, z_i, res_i) = ca_data[i];
        for j in (i + 1)..l {
            let (x_j, y_j, z_j, res_j) = ca_data[j];
 
            let dx = x_i - x_j;
            let dy = y_i - y_j;
            let dz = z_i - z_j;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt(); 

            if dist <= 8.0 {
                let seq_sep = (res_j - res_i).abs() as f32;
                sum_seqsep += seq_sep;
                n_contacts += 1.0;
            }
        }
    }

    if n_contacts == 0.0 {
        0.0
    } else {
       
        sum_seqsep / ((l as f32) * n_contacts) * 100.0
    }
}


pub fn return_len(pdb: &PDB) -> usize {

    pdb
        .models()
        .flat_map(|model| model.chains())
        .flat_map(|chain| chain.residues())
        .flat_map(|residue| residue.atoms())
        .filter(|atom| atom.record_name() == "ATOM" && atom.atom_name().trim() == "CA")
        .count()
}

*/