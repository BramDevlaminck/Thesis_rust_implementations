pub mod taxon_id_calculator;

use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::from_utf8;
use bytelines::ByteLines;
use umgap::taxon::TaxonId;
use sa_mappings::proteins::TERMINATION_CHARACTER;
use sa_mappings::taxonomy::TaxonAggregator;
use crate::taxon_id_calculator::{TaxonIdVerifier};

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

pub struct Proteins {
    pub input_string: Vec<u8>,
    pub proteins: Vec<Protein>
}

impl Proteins {

    pub fn get_sequence(&self, protein: &Protein) ->&str {
        let (begin, size) = protein.sequence;
        std::str::from_utf8( &self.input_string[begin..begin + size as usize]).unwrap() // should never fail since the input string will always be utf8
    }
}

#[derive(Debug)]
/// The useful information about a protein for our use case
pub struct Protein {
    pub uniprot_id: String,
    pub sequence: (usize, u32),
    pub id: TaxonId,
}

/// Parse the given database tsv file into a Vector of Proteins with the data from the tsv file
pub fn get_proteins_from_database_file(database_file: &str, taxon_id_calculator: &dyn TaxonIdVerifier) -> Result<Proteins, Box<dyn Error>> {
    let mut input_string: String = "".to_string();
    let mut proteins: Vec<Protein> = vec![];
    let mut begin_index: usize = 0;
    let lines = read_lines(database_file)?;
    for line in lines.into_iter().map_while(Result::ok) {
        let parts: Vec<String> = line.split('\t').map(str::to_string).collect();
        let [uniprot_id, protein_id_str, protein_sequence]: [String; 3] = parts.try_into().map_err(|e| DatabaseFormatError{ error: e})?;
        let protein_id_as_taxon_id = protein_id_str.parse::<TaxonId>()?;
        // if the taxon ID is not a valid ID in our NCBI taxonomy, skip this protein
        if !taxon_id_calculator.taxon_id_exists(protein_id_as_taxon_id) {
            // eprintln!("Skipped protein with taxon id {}!", protein_id_as_taxon_id);
            continue;
        }

        if begin_index != 0 {
            input_string.push(SEPARATION_CHARACTER as char);
        }
        input_string.push_str(&protein_sequence.to_uppercase());
        proteins.push(
            Protein {
                uniprot_id,
                sequence: (begin_index, protein_sequence.len() as u32),
                id: protein_id_as_taxon_id,
            }
        );
        begin_index += protein_sequence.len() + 1;
    }
    input_string.push(END_CHARACTER as char);
    Ok(Proteins {
        input_string: input_string.into_bytes(),
        proteins
    })
}

/// Parse the given database tsv file into a String that has all the proteins concatenated as 1 large text
pub fn get_text_from_database_file(database_file: &str, taxon_aggregator: &TaxonAggregator) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut input_string: String = String::new();

    let file = File::open(database_file)?;
    
    // Read the lines as bytes, since the input string is not guaranteed to be utf8
    // because of the encoded functional annotations
    let mut lines = ByteLines::new(BufReader::new(file));

    while let Some(Ok(line)) = lines.next() {
        let mut fields = line.split(|b| *b == b'\t');

        // only get the taxon id and sequence from each line, we don't need the other parts
        fields.next();
        let taxon_id = from_utf8(fields.next().unwrap())?.parse::<TaxonId>()?;
        let sequence = from_utf8(fields.next().unwrap())?;
        fields.next();

        if !taxon_aggregator.taxon_exists(taxon_id) {
            continue;
        }

        input_string.push_str(&sequence.to_uppercase());
        input_string.push(SEPARATION_CHARACTER.into());
    }

    input_string.pop();
    input_string.push(TERMINATION_CHARACTER.into());
    
    input_string.shrink_to_fit();
    Ok(input_string.into_bytes())
}

#[derive(Debug)]
struct DatabaseFormatError {
    error: Vec<String>
}

impl std::fmt::Display for DatabaseFormatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Expected the protein database file to have the following fields separated by a tab: <Uniprot_accession> <protein id> <sequence>\nBut tried to unpack following vector in 3 variables: {:?}", self.error)
    }
}

impl Error for DatabaseFormatError {}