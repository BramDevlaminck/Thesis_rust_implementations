use clap::Parser;
use suffixarray_builder::{Arguments, build_sa};
use suffixarray_builder::binary::write_binary;
use tsv_utils::get_text_from_database_file;
use tsv_utils::taxon_id_calculator::{AggregationMethod, TaxonIdCalculator};

fn main() {
    let args = Arguments::parse();
    let Arguments { database_file, taxonomy, output, sample_rate, construction_algorithm } = args;
    let taxon_id_calculator = TaxonIdCalculator::new(&taxonomy, AggregationMethod::LcaStar);
    
    // read input
    let data = get_text_from_database_file(&database_file, &*taxon_id_calculator);
    if let Err(err) = data {
        eprintln!("{}", err);
        std::process::exit(1);
    }
    let data = data.unwrap();
    
    // calculate sa
    let sa = build_sa(&data.into_bytes(), &construction_algorithm, sample_rate);
    if let Err(err) = sa {
        eprintln!("{}", err);
        std::process::exit(1);
    }
    let sa = sa.unwrap();
    
    // output the build SA
    if let Err(err) = write_binary(sample_rate, &sa, &output) {
        eprintln!("{}", err);
        std::process::exit(1);
    };
}