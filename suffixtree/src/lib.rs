use std::io;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{arg, Parser, ValueEnum};

use tsv_utils::{get_proteins_from_database_file, Protein, proteins_to_concatenated_string, read_lines};

use crate::searcher::Searcher;
use crate::tree_taxon_id_calculator::TreeTaxonIdCalculator;
use crate::tree::Tree;
use crate::tree_builder::{TreeBuilder, UkkonenBuilder};

mod tree_builder;
mod tree;
mod cursor;
mod search_cursor;
mod searcher;
mod tree_taxon_id_calculator;


/// Enum that represents the 2 kinds of search that we support
/// - Search until match and return boolean that indicates if there is a match
/// - Search until match, if there is a match search the whole subtree to find all matching proteins
/// - Search until match, there we can immediately retrieve the taxonId that represents all the children
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SearchMode {
    Match,
    AllOccurrences,
    TaxonId,
}

#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be concatenated using a `#`.
    #[arg(short, long)]
    database_file: String,
    /// A file that contains sequences that we want to search in the tree. Every line contains a new sequence.
    #[arg(short, long)]
    search_file: Option<String>,
    /// This will only build the tree and stop after that is completed. Used during benchmarking.
    #[arg(long)]
    build_only: bool,
    /// `match` will only look if there is match.
    /// While `all-occurrences` will search for the match and look for all the different matches in the subtree.
    /// `Taxon-id` will search for the matching taxon id using lca*
    #[arg(short, long, value_enum)]
    mode: Option<SearchMode>,
    /// This will change the output to <found (0 or 1)>;<protein length>;<search time in ms>
    /// The given num will be used to run the search x times and the average of these x runs will be given as search time
    #[arg(short, long)]
    verbose: Option<u8>,
    #[arg(short, long)]
    /// The taxonomy to be used as a tsv file. This is a preprocessed version of the NCBI taxonomy.
    taxonomy: String,
    #[arg(long)]
    /// Prints all the taxon ids in the tree in pre-order traversal
    print_tree_taxon_ids: bool,
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
fn handle_search_word(searcher: &mut Searcher, word: String, search_mode: &SearchMode, verbose: Option<u8>, verbose_output: &mut Vec<String>) {
    let word = match word.strip_suffix('\n') {
        None => word,
        Some(stripped) => String::from(stripped)
    }.to_uppercase();
    if let Some(num_iter) = verbose {
        let mut found_total: bool = false;
        let mut total_time: f64 = 0.0;
        for _ in 0..num_iter {
            let (found, execution_time) = match *search_mode {
                SearchMode::Match => time_execution(searcher, &|searcher| searcher.search_if_match(word.as_bytes())),
                SearchMode::AllOccurrences => time_execution(searcher, &|searcher| !searcher.find_all_suffix_indices(word.as_bytes()).is_empty()),
                SearchMode::TaxonId => time_execution(searcher, &|searcher| searcher.search_taxon_id(word.as_bytes()).is_some()),
            };
            total_time += execution_time;
            found_total = found;
        }
        let avg = total_time / (num_iter as f64);

        verbose_output.push(format!("{};{};{}", found_total as u8, word.len(), avg));
    } else {
        match *search_mode {
            SearchMode::Match => println!("{}", searcher.search_if_match(word.as_bytes())),
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
}


/// Main run function that executes all the logic with the received arguments
pub fn run(args: Arguments) {
    let proteins = get_proteins_from_database_file(&args.database_file);
    // construct the sequence that will be used to build the tree
    let data = proteins_to_concatenated_string(&proteins);

    // build the tree
    let mut tree = Tree::new(&data, UkkonenBuilder::new());
    // fill in the Taxon Ids in the tree using the LCA implementations from UMGAP
    let tree_taxon_id_calculator = TreeTaxonIdCalculator::new(&args.taxonomy);
    tree_taxon_id_calculator.calculate_taxon_ids(&mut tree, &proteins);

    // print the taxon ids of the tree if flag set
    if args.print_tree_taxon_ids {
        TreeTaxonIdCalculator::output_all_taxon_ids(&tree);
    }

    // option that only builds the tree, but does not allow for querying (easy for benchmark purposes)
    if args.build_only {
        return;
    } else if args.mode.is_none() {
        eprintln!("search mode expected!");
        std::process::exit(1);
    }

    let searcher = Searcher::new(&tree, data.as_bytes(), &proteins, &tree_taxon_id_calculator);
    execute_search(searcher, &args)
}

/// Perform the search as set with the commandline arguments
fn execute_search(mut searcher: Searcher, args: &Arguments) {
    let mode = args.mode.as_ref().unwrap();
    let verbose = args.verbose;
    let mut verbose_output: Vec<String> = vec![];
    if let Some(search_file) = &args.search_file {
        // File `search_file` must exist in the current path
        if let Ok(lines) = read_lines(search_file) {
            for line in lines.into_iter().flatten() {
                handle_search_word(&mut searcher, line, mode, verbose, &mut verbose_output);
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
            handle_search_word(&mut searcher, word, mode, verbose, &mut verbose_output);
        }
    }
    verbose_output.iter().for_each(|val| println!("{}", val));
}