pub mod taxon_id_calculator;

use std::fs::File;
use std::io;
use std::io::BufRead;
use std::ops::Add;
use std::path::Path;
use umgap::taxon::TaxonId;

// END_CHARACTER should ALWAYS be lexicographically than SEPARATION_CHARACTER
// otherwise the building of the suffix array will not happen correctly
pub static SEPARATION_CHARACTER: u8 = b'-';
pub static END_CHARACTER: u8 = b'$';


// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

/// The useful information about a protein for our use case
pub struct Protein {
    pub sequence: String,
    pub id: TaxonId,
}

/// Parse the given database tsv file into a Vector of Proteins with the data from the tsv file
pub fn get_proteins_from_database_file(database_file: &str) -> Vec<Protein> {
    let mut proteins: Vec<Protein> = vec![];
    if let Ok(lines) = read_lines(database_file) {
        for line in lines.into_iter().flatten() {
            let [_, _, protein_id_str, _, _, protein_sequence]: [&str; 6] = line.splitn(6, '\t').collect::<Vec<&str>>().try_into().unwrap();
            let protein_id_as_taxon_id = protein_id_str.parse::<TaxonId>().expect("Could not parse id of protein to usize!");
            proteins.push(
                Protein {
                    sequence: protein_sequence.to_uppercase(),
                    id: protein_id_as_taxon_id,
                }
            )
        }
    } else {
        eprintln!("Database file {} could not be opened!", database_file);
        std::process::exit(1);
    }
    proteins
}

/// Joins the sequences of all the proteins together with `#` as delimiter and adds a `$` at the end
pub fn proteins_to_concatenated_string(proteins: &[Protein]) -> String {
    proteins
        .iter()
        .map(|prot| prot.sequence.clone())
        .collect::<Vec<String>>()
        .join(&(SEPARATION_CHARACTER as char).to_string())
        .add(&(END_CHARACTER as char).to_string())
}