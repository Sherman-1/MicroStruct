use clap::{Arg, Command};

pub struct Config {
    pub num_cpus: usize,
    pub pdb_dir: String,
    pub output_file: String,
    pub subset: Option<usize>
}

pub fn parse_arguments() -> Config {
    let matches = Command::new("PDB Processing App")
        .version("1.0")
        .author("Simon Herman")
        .about("Processes PDB files and computes metrics")
        .arg(
            Arg::new("cpus")
                .short('c')
                .long("cpus")
                .value_name("CPUS")
                .help("Number of CPUs to use")
                .default_value("40") 
        )
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .value_name("INPUT_DIR")
                .help("Path to the PDB directory")
                .default_value("/datas/SIMON/AFDB_20_to_100/pdbs/") 
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("OUTPUT_FILE")
                .help("Path to the output CSV file")
                .default_value("/store/EQUIPES/BIM/MEMBERS/simon.herman/MicroStruct/results.csv")
        )
        .arg(
            Arg::new("subset")
                .short('s')
                .long("subset")
                .value_name("SUBSET")
                .help("Number of random PDB files to process (optional)")
                .required(false) 
        )
        .get_matches();

            
    let num_cpus: usize = matches.get_one::<String>("cpus")
        .unwrap()
        .parse::<usize>()
        .expect("Invalid CPU count");
    let pdb_dir = matches.get_one::<String>("input")
        .unwrap()
        .to_string();
    let output_file = matches.get_one::<String>("output")
        .unwrap()
        .to_string();
    let subset: Option<usize> = matches.get_one::<String>("subset")
        .and_then(|s| s.parse::<usize>().ok()); 


    Config {
        num_cpus,
        pdb_dir,
        output_file,
        subset
    }
}
