# MicroStruct

my-protein-analysis/
├── Cargo.toml
├── README.md
├── data/        <-- Potential directory for sample test data or reference data
└── crates/
    ├── pdb-io/  <-- Crate for parsing PDB files
    ├── metrics/ <-- Crate for computing geometry or contact metrics
    └── main-app/ <-- The main binary crate
        ├── src/
        │   └── main.rs
        └── Cargo.toml