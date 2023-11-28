mod searcher;

use std::error::Error;
use std::io;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};
use clap::{arg, Parser, ValueEnum};
use get_size::GetSize;

use tsv_utils::{get_proteins_from_database_file, proteins_to_concatenated_string, read_lines};
use tsv_utils::taxon_id_calculator::TaxonIdCalculator;
use crate::searcher::Searcher;

/// Enum that represents the 2 kinds of search that we support
/// - Search until match and return boolean that indicates if there is a match
/// - Search until match, if there is a match return the min and max index in the SA that matches
/// - Search until match, if there is a match search the whole subtree to find all matching proteins
/// - Search until match, there we can immediately retrieve the taxonId that represents all the children
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SearchMode {
    Match,
    MinMaxBound,
    AllOccurrences,
    TaxonId,
}

#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be concatenated using a `#`.
    #[arg(short, long)]
    database_file: String,
    #[arg(short, long)]
    search_file: Option<String>,
    /// `match` will only look if there is match.
    /// `all-occurrences` will search for the match and look for all the different matches in the subtree.
    /// `min-max-bound` will search for the match and retrieve the minimum and maximum index in the SA that contains a suffix that matches.
    /// `Taxon-id` will search for the matching taxon id using lca*
    #[arg(short, long, value_enum)]
    mode: Option<SearchMode>,
    #[arg(short, long)]
    /// The taxonomy to be used as a tsv file. This is a preprocessed version of the NCBI taxonomy.
    taxonomy: String,
}

pub fn run(args: Arguments) -> Result<(), Box<dyn Error>> {
    let proteins = get_proteins_from_database_file(&args.database_file);
    // construct the sequence that will be used to build the tree
    let data = proteins_to_concatenated_string(&proteins);
    let u8_text = data.as_bytes();

    let sa = libdivsufsort_rs::divsufsort64(&u8_text.to_vec()).ok_or("Building suffix array failed")?;

    println!("{}", sa.len());
    println!("{}", sa[0]);

    let mut current_protein_index: u32 = 0;
    let mut index_to_protein: Vec<Option<u32>> = vec![];
    for &char in u8_text.iter() {
        if char == b'-' || char == b'$' {
            current_protein_index += 1;
            index_to_protein.push(None);
        } else {
            index_to_protein.push(Some(current_protein_index));
        }
    }
    println!("{}", sa.get_size() + u8_text.get_size() + index_to_protein.get_size()); // print mem size of structures

    let taxon_id_calculator = TaxonIdCalculator::new(&args.taxonomy);

    let searcher = Searcher::new(u8_text, &sa, &index_to_protein, &proteins, &taxon_id_calculator);
    execute_search(searcher, &args);
    Ok(())
}

/// Perform the search as set with the commandline arguments
fn execute_search(mut searcher: Searcher, args: &Arguments) {
    let mode = args.mode.as_ref().unwrap();
    // let verbose = args.verbose;
    let mut verbose_output: Vec<String> = vec![];
    if let Some(search_file) = &args.search_file {
        // File `search_file` must exist in the current path
        if let Ok(lines) = read_lines(search_file) {
            for line in lines.into_iter().flatten() {
                handle_search_word(&mut searcher, line, mode);
            }
        } else {
            eprintln!("File {} could not be opened!", search_file);
            std::process::exit(1);
        }
    } else {
        loop {
            print!("Input your search string: ");
            io::stdout().flush().unwrap();
            let mut word = String::new();

            if io::stdin().read_line(&mut word).is_err() {
                continue;
            }
            handle_search_word(&mut searcher, word, mode);
        }
    }
    verbose_output.iter().for_each(|val| println!("{}", val));
}

fn time_execution(searcher: &mut Searcher, f: &dyn Fn(&mut Searcher) -> bool) -> (bool, f64) {
    let start_ms = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
    let found = f(searcher);
    let end_ms = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
    (found, end_ms - start_ms, )
}


/// Executes the kind of search indicated by the commandline arguments
fn handle_search_word(searcher: &mut Searcher, word: String, search_mode: &SearchMode) {
    let word = match word.strip_suffix('\n') {
        None => word,
        Some(stripped) => String::from(stripped)
    }.to_uppercase();
    match *search_mode {
        SearchMode::Match => println!("{}", searcher.search_if_match(word.as_bytes())),
        SearchMode::MinMaxBound => {
            let (found, min_bound, max_bound) = searcher.search_bounds(word.as_bytes());
            println!("{found};{min_bound};{max_bound}");
        }
        SearchMode::AllOccurrences => {
            let results = searcher.search_protein(word.as_bytes());
            println!("found {} matches", results.len());
            results.iter()
                .for_each(|res| println!("* {}", res.sequence));
        }
        SearchMode::TaxonId => {
            match searcher.search_taxon_id(word.as_bytes()) {
                Some(taxon_id) => println!("{}", taxon_id),
                None => println!("/"),
            }
        }
    }
}
