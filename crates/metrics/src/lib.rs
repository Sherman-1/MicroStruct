use pdb_io::ParsedPDB;

pub fn radius_of_gyration(pdb: &ParsedPDB) -> f32 {
    let n = pdb.ca_coords.len();
    if n == 0 {
        return 0.0;
    }

    // Sum all coordinates
    let (sum_x, sum_y, sum_z) = pdb
        .ca_coords
        .iter()
        .fold((0.0, 0.0, 0.0), |(sx, sy, sz), ca| {
            (sx + ca.x, sy + ca.y, sz + ca.z)
        });

    let cx = sum_x / n as f32;
    let cy = sum_y / n as f32;
    let cz = sum_z / n as f32;

    // Compute sum of square distances
    let mut sum_dist_sq = 0.0;
    for ca in &pdb.ca_coords {
        let dx = ca.x - cx;
        let dy = ca.y - cy;
        let dz = ca.z - cz;
        sum_dist_sq += dx * dx + dy * dy + dz * dz;
    }

    (sum_dist_sq / n as f32).sqrt()
}

pub fn plddt_statistics(pdb: &ParsedPDB) -> (f32, f32, f32, f32) {
    /* Returns:
        - mean pLDDT
        - fraction of residues over 50 pLDDT
        - fraction of residues over 70 pLDDT
        - fraction of residues over 90 pLDDT    
    */

    let n = pdb.ca_coords.len();
    if n == 0 {
        return (0.0, 0.0, 0.0, 0.0);
    }

    let (sum_plddt, count_50, count_70, count_90) = pdb
        .ca_coords
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
    let fraction_50 = count_50 as f32 / n as f32;
    let fraction_70 = count_70 as f32 / n as f32;
    let fraction_90 = count_90 as f32 / n as f32;

    (mean_plddt, fraction_50, fraction_70, fraction_90)

}


pub fn bounding_box_volume(pdb: &ParsedPDB) -> f32 {
    let n = pdb.ca_coords.len();
    if n == 0 {
        return 0.0;
    }

    let mut x_min = f32::MAX;
    let mut x_max = f32::MIN;
    let mut y_min = f32::MAX;
    let mut y_max = f32::MIN;
    let mut z_min = f32::MAX;
    let mut z_max = f32::MIN;

    for ca in &pdb.ca_coords {
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

pub fn contact_order(pdb: &ParsedPDB) -> f32 {

    // Based on https://en.wikipedia.org/wiki/Contact_order definition 

    let coords = &pdb.ca_coords;
    let l = coords.len();
    if l < 20 {
        return 0.0;
    }

    let mut sum_seqsep = 0.0;
    let mut n = 0.0; // Number of contacts detected 

    for i in 0..l {
        for j in (i + 1)..l {
            let dx = coords[i].x - coords[j].x;
            let dy = coords[i].y - coords[j].y;
            let dz = coords[i].z - coords[j].z;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            if dist <= 8.0 {
                // If residues are in contact, compute the distance within the 1D sequence 
                let seq_sep = (coords[j].residue_seq - coords[i].residue_seq).abs() as f32;
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
