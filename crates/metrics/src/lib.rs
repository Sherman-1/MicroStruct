use pdb_io::ParsedPDB;
use pdb_io::AtomCoordinate;
//pub mod sasa;


pub fn get_ca_atoms<'a>(pdb: &'a ParsedPDB) -> Vec<&'a AtomCoordinate> {
    pdb.atoms
        .iter()
        .filter(|atom| atom.atom_name.trim() == "CA")
        .collect()
}
pub fn radius_of_gyration(ca_atoms: &[&AtomCoordinate]) -> f32 {
    let n = ca_atoms.len();
    if n == 0 {
        return 0.0;
    }

    // Compute the centroid.
    let (sum_x, sum_y, sum_z) = ca_atoms.iter().fold((0.0, 0.0, 0.0), |(sx, sy, sz), ca| {
        (sx + ca.x, sy + ca.y, sz + ca.z)
    });
    let cx = sum_x / n as f32;
    let cy = sum_y / n as f32;
    let cz = sum_z / n as f32;

    // Compute the sum of squared distances.
    let sum_dist_sq = ca_atoms.iter().fold(0.0, |acc, ca| {
        let dx = ca.x - cx;
        let dy = ca.y - cy;
        let dz = ca.z - cz;
        acc + dx * dx + dy * dy + dz * dz
    });

    (sum_dist_sq / n as f32).sqrt()
}

pub fn plddt_statistics(ca_atoms: &[&AtomCoordinate]) -> (f32, f32, f32, f32) {
    /* Returns:
       - mean pLDDT
       - fraction of residues over 50 pLDDT
       - fraction of residues over 70 pLDDT
       - fraction of residues over 90 pLDDT    
    */
    let n = ca_atoms.len();
    if n == 0 {
        return (0.0, 0.0, 0.0, 0.0);
    }

    let (sum_plddt, count_50, count_70, count_90) = ca_atoms
        .iter()
        .fold((0.0, 0, 0, 0), |(sum, c50, c70, c90), ca| {
            (
                sum + ca.bfactor,
                c50 + (ca.bfactor > 50.0) as usize,
                c70 + (ca.bfactor > 70.0) as usize,
                c90 + (ca.bfactor > 90.0) as usize,
            )
        });

    let mean_plddt = sum_plddt / n as f32;
    let fraction_50 = count_50 as f32 / n as f32 * 100.0;
    let fraction_70 = count_70 as f32 / n as f32 * 100.0;
    let fraction_90 = count_90 as f32 / n as f32 * 100.0;

    (mean_plddt, fraction_50, fraction_70, fraction_90)
}

pub fn bounding_box_volume(ca_atoms: &[&AtomCoordinate]) -> f32 {
    let n = ca_atoms.len();
    if n == 0 {
        return 0.0;
    }

    let mut x_min = f32::MAX;
    let mut x_max = f32::MIN;
    let mut y_min = f32::MAX;
    let mut y_max = f32::MIN;
    let mut z_min = f32::MAX;
    let mut z_max = f32::MIN;

    for ca in ca_atoms {
        if ca.x < x_min { x_min = ca.x; }
        if ca.x > x_max { x_max = ca.x; }
        if ca.y < y_min { y_min = ca.y; }
        if ca.y > y_max { y_max = ca.y; }
        if ca.z < z_min { z_min = ca.z; }
        if ca.z > z_max { z_max = ca.z; }
    }

    let dx = x_max - x_min;
    let dy = y_max - y_min;
    let dz = z_max - z_min;

    dx * dy * dz
}

pub fn contact_order(ca_atoms: &[&AtomCoordinate]) -> f32 {
    let l = ca_atoms.len();
    if l < 20 {
        return 0.0;
    }

    let mut sum_seqsep = 0.0;
    let mut n = 0.0; // Number of contacts detected

    for i in 0..l {
        for j in (i + 1)..l {
            let dx = ca_atoms[i].x - ca_atoms[j].x;
            let dy = ca_atoms[i].y - ca_atoms[j].y;
            let dz = ca_atoms[i].z - ca_atoms[j].z;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            if dist <= 8.0 {
                let seq_sep = (ca_atoms[j].residue_seq - ca_atoms[i].residue_seq).abs() as f32;
                sum_seqsep += seq_sep;
                n += 1.0;
            }
        }
    }

    if n == 0.0 {
        0.0
    } else {
        sum_seqsep / ((l as f32) * n) * 100.0
    }
}

pub fn return_len(ca_atoms: &[&AtomCoordinate]) -> usize {
    ca_atoms.len()
}
