use pdb_io::{parse_pdb};
use metrics::{radius_of_gyration, bounding_box_volume, contact_order, plddt_statistics};

#[test]
fn test_metrics_basic() {
    // Vraiment vraiment crade
    let pdb = parse_pdb("/store/EQUIPES/BIM/MEMBERS/simon.herman/MicroStruct/test.pdb").expect("File sadly not found ... ");
    let rg = radius_of_gyration(&pdb);
    let vol = bounding_box_volume(&pdb);
    let co = contact_order(&pdb);
    let plddt = plddt_statistics(&pdb);

    assert!(rg >= 0.0);
    assert!(vol >= 0.0);
    assert!(co >= 0.0);

    println!("Radius : {rg}, Volume : {vol}, Contact : {co}, pLDDT : {plddt:?}");
}
