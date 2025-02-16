
use pdbtbx::*;

#[test]
fn test_metrics_basic() {

    let Ok((pdb, _error)) = pdbtbx::open("/store/EQUIPES/BIM/MEMBERS/simon.herman/MicroStruct/data/example.cif");


    let all_ca_bfactors: Vec<f32> = pdb
        .models()
        .flat_map(|model| model.chains())
        .flat_map(|chain| chain.residues())
        .flat_map(|res| res.atoms())
        .filter(|atom| atom.record_name() == "ATOM" && atom.atom_name().trim() == "CA")
        .map(|atom| atom.bfactor() as f32)
        .collect();
    
}
