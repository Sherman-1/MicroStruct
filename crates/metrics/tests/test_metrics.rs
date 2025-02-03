use pdb_io::parse_pdb;
use metrics::{radius_of_gyration, 
    bounding_box_volume, 
    contact_order, 
    plddt_statistics, 
    get_ca_atoms};
//use metrics::sasa::{calc_sasa_from_parsed_pdb,SasaMode};

#[test]
fn test_metrics_basic() {

    let pdb = parse_pdb("/store/EQUIPES/BIM/MEMBERS/simon.herman/MicroStruct/test.pdb").expect("File sadly not found ... ");
    
    let ca_atoms = get_ca_atoms(&pdb);
    let rg = radius_of_gyration(&ca_atoms);
    let vol = bounding_box_volume(&ca_atoms);
    let co = contact_order(&ca_atoms);
    let plddt = plddt_statistics(&ca_atoms);

    //let sasa = calc_sasa_from_parsed_pdb(&pdb, 10, SasaMode::Full);


    assert!(rg >= 0.0);
    assert!(vol >= 0.0);
    assert!(co >= 0.0);

    println!("Radius : {rg}, Volume : {vol}, Contact : {co}, pLDDT : {plddt:?}");
    

}
