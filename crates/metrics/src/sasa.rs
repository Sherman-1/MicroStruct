use ndarray::{array, Array1, Array2};
use std::f64::consts::PI;
use pdb_io::{ParsedPDB, AtomCoordinate};

/// Mode for SASA computation: either use all atoms for both target and occluders,
/// or compute SASA only at CA positions (while still using all atoms as occluders).
pub enum SasaMode {
    Full,
    CAOnly,
}

fn golden_spiral(num_pts: usize, _radius: f64) -> Array2<f64> {
    let mut points = Array2::<f64>::zeros((num_pts, 3));
    for i in 0..num_pts {
        let index = i as f64 + 0.5;
        let phi = (1.0 - 2.0 * index / (num_pts as f64)).acos();
        let theta = PI * (1.0 + 5.0_f64.sqrt()) * index;
        points[[i, 0]] = theta.cos() * phi.sin();
        points[[i, 1]] = theta.sin() * phi.sin();
        points[[i, 2]] = phi.cos();
    }
    points
}

/// Generalized SASA calculation.
/// 
/// * `target_xyz` and `target_radii` are the coordinates and radii for the atoms
///   for which SASA will be computed.
/// * `occluder_xyz` and `occluder_radii` are used for the accessibility check.
/// * `n` is the number of sampling points per atom.
fn calc_sasa_general(
        target_xyz: &Array2<f64>,
        target_radii: &Array1<f64>,
        occluder_xyz: &Array2<f64>,
        occluder_radii: &Array1<f64>,
        n: usize,
    ) -> Array1<f64> {
    let n_target = target_xyz.shape()[0];
    let pts = golden_spiral(n, 1.0);
    let mut sasa = Array1::<f64>::zeros(n_target);

    for (i, atom_center) in target_xyz.outer_iter().enumerate() {
        let inflated_radius = target_radii[i] + 1.4;
        let r_i = inflated_radius + 1e-5;
        let mut count_accessible = 0;

        // For each sample point on the atom's inflated sphere…
        for pt in pts.outer_iter() {
            let sp = array![
                atom_center[0] + r_i * pt[0],
                atom_center[1] + r_i * pt[1],
                atom_center[2] + r_i * pt[2]
            ];
            let mut accessible = true;
            // Check against every occluder.
            for (j, occ_center) in occluder_xyz.outer_iter().enumerate() {
                let threshold = occluder_radii[j] + 1.4;
                let dx = sp[0] - occ_center[0];
                let dy = sp[1] - occ_center[1];
                let dz = sp[2] - occ_center[2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if d <= threshold {
                    accessible = false;
                    break;
                }
            }
            if accessible {
                count_accessible += 1;
            }
        }
        let fraction_outside = count_accessible as f64 / n as f64;
        sasa[i] = fraction_outside * (4.0 * PI * inflated_radius * inflated_radius);
    }
    sasa
}

/// A simple heuristic to assign a van der Waals radius based on the atom's name.
fn get_vdw_radius(atom: &AtomCoordinate) -> f64 {
    let name = atom.atom_name.trim();
    match name.chars().next().unwrap_or(' ') {
        'H' => 1.2,
        'C' => 1.7,
        'N' => 1.55,
        'O' => 1.52,
        'S' => 1.8,
        _   => 1.5,
    }
}

/// Computes SASA from a ParsedPDB. The user can choose to compute SASA for:
/// - All atoms (mode = SasaMode::Full), or
/// - Only for CA atoms (mode = SasaMode::CAOnly),
/// while using all atoms as occluders.
pub fn calc_sasa_from_parsed_pdb(pdb: &ParsedPDB, n: usize, mode: SasaMode) -> Array1<f64> {
    // Build occluder arrays from all atoms.
    let num_full = pdb.atoms.len();
    let mut occ_coords = Array2::<f64>::zeros((num_full, 3));
    let mut occ_radii = Array1::<f64>::zeros(num_full);
    for (i, atom) in pdb.atoms.iter().enumerate() {
        occ_coords[[i, 0]] = atom.x as f64;
        occ_coords[[i, 1]] = atom.y as f64;
        occ_coords[[i, 2]] = atom.z as f64;
        occ_radii[i] = get_vdw_radius(atom);
    }

    // Depending on mode, choose the target atoms.
    let (target_coords, target_radii) = match mode {
        SasaMode::Full => (occ_coords.clone(), occ_radii.clone()),
        SasaMode::CAOnly => {
            // Filter only CA atoms for SASA computation.
            let ca_atoms: Vec<&AtomCoordinate> = pdb.atoms
                .iter()
                .filter(|atom| atom.atom_name.trim() == "CA")
                .collect();
            let n_ca = ca_atoms.len();
            let mut coords = Array2::<f64>::zeros((n_ca, 3));
            let mut radii = Array1::<f64>::zeros(n_ca);
            for (i, atom) in ca_atoms.iter().enumerate() {
                coords[[i, 0]] = atom.x as f64;
                coords[[i, 1]] = atom.y as f64;
                coords[[i, 2]] = atom.z as f64;
                radii[i] = get_vdw_radius(atom);
            }
            (coords, radii)
        }
    };

    calc_sasa_general(&target_coords, &target_radii, &occ_coords, &occ_radii, n)
}
