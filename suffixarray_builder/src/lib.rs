pub mod binary;

use std::error::Error;
use clap::{Parser, ValueEnum};

#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be concatenated using a `#`.
    #[arg(short, long)]
    pub database_file: String,
    #[arg(short, long)]
    /// The taxonomy to be used as a tsv file. This is a preprocessed version of the NCBI taxonomy.
    pub taxonomy: String,
    /// Output file to store the built index.
    #[arg(short, long)]
    pub output: String,
    /// The sample rate used on the suffix array (default value 1, which means every value in the SA is used)
    #[arg(long, default_value_t = 1)]
    pub sample_rate: u8,
    #[arg(short, long, value_enum, default_value_t = SAConstructionAlgorithm::LibSais)]
    pub construction_algorithm: SAConstructionAlgorithm,
}

#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SAConstructionAlgorithm {
    LibDivSufSort,
    LibSais,
}

pub fn build_sa(data: &Vec<u8>, construction_algorithm: &SAConstructionAlgorithm, sample_rate: u8) -> Result<Vec<i64>, Box<dyn Error>> {

    let mut sa = match construction_algorithm {
        SAConstructionAlgorithm::LibSais => libsais64_rs::sais64(data),
        SAConstructionAlgorithm::LibDivSufSort => {
            libdivsufsort_rs::divsufsort64(data)
        }
    }.ok_or("Building suffix array failed")?;

    // make the SA sparse and decrease the vector size if we have sampling (== sampling_rate > 1)
    if sample_rate > 1 {
        let mut current_sampled_index = 0;
        for i in 0..sa.len() {
            let current_sa_val = sa[i];
            if current_sa_val % sample_rate as i64 == 0 {
                sa[current_sampled_index] = current_sa_val;
                current_sampled_index += 1;
            }
        }
        // make shorter
        sa.resize(current_sampled_index, 0);
    }

    Ok(sa)
}