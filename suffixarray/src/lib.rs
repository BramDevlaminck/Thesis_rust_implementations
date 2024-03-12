use std::error::Error;
use std::num::NonZeroUsize;

use clap::{arg, Parser, ValueEnum};
use rayon::prelude::*;
use suffixarray_builder::{build_sa, SAConstructionAlgorithm};
use suffixarray_builder::binary::{load_binary, write_binary};

use tsv_utils::taxon_id_calculator::{AggregationMethod, TaxonIdCalculator};
use tsv_utils::{get_proteins_from_database_file, read_lines};

use crate::searcher::Searcher;
use crate::suffix_to_protein_index::{
    DenseSuffixToProtein, SparseSuffixToProtein, SuffixToProteinIndex, SuffixToProteinMappingStyle,
};
use crate::util::get_time_ms;

pub mod searcher;
pub mod suffix_to_protein_index;
pub mod util;
/// Enum that represents the 5 kinds of search that we support
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SearchMode {
    Match,
    MinMaxBound,
    AllOccurrences,
    TaxonId,
    Analyses
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
    /// `Analyses` will return all the Unipept analyses results
    #[arg(short, long, value_enum)]
    mode: Option<SearchMode>,
    #[arg(short, long)]
    /// The taxonomy to be used as a tsv file. This is a preprocessed version of the NCBI taxonomy.
    taxonomy: String,
    /// This will only build the tree and stop after that is completed. Used during benchmarking.
    #[arg(long)]
    build_only: bool,
    /// Output file to store the built index.
    #[arg(short, long)]
    output: Option<String>,
    /// The sample rate used on the suffix array (default value 1, which means every value in the SA is used)
    #[arg(long, default_value_t = 1)]
    sample_rate: u8,
    /// Set the style used to map back from the suffix to the protein. 2 options <sparse> or <dense>. Dense is default
    /// Dense uses O(n) memory with n the size of the input text, and takes O(1) time to find the mapping
    /// Sparse uses O(m) memory with m the number of proteins, and takes O(log m) to find the mapping
    #[arg(long, value_enum, default_value_t = SuffixToProteinMappingStyle::Sparse)]
    suffix_to_protein_mapping: SuffixToProteinMappingStyle,
    #[arg(long)]
    load_index: Option<String>,
    #[arg(short, long, value_enum, default_value_t = SAConstructionAlgorithm::LibSais)]
    construction_algorithm: SAConstructionAlgorithm,
    /// Assume the resulting taxon ID is root (1) whenever a peptide matches >= cutoff proteins
    #[arg(long, default_value_t = 10000)]
    cutoff: usize,
    #[arg(long)]
    threads: Option<NonZeroUsize>,
}

pub fn run(mut args: Arguments) -> Result<(), Box<dyn Error>> {
    let start_reading_proteins_ms = get_time_ms()?;
    let taxon_id_calculator = TaxonIdCalculator::new(&args.taxonomy, AggregationMethod::LcaStar);
    // println!("taxonomy calculator built");

    let proteins = get_proteins_from_database_file(&args.database_file, &*taxon_id_calculator)?;

    // construct the sequence that will be used to build the tree
    // println!("read all proteins");
    let current = get_time_ms()?;
    // println!("Time spent for reading: {}", current - start_reading_proteins_ms);

    let sa = match &args.load_index {
        // load SA from file
        Some(index_file_name) => {
            let start_loading_ms = get_time_ms()?;
            let (sample_rate, sa) = load_binary(index_file_name)?;
            args.sample_rate = sample_rate;
            let end_loading_ms = get_time_ms()?;
            // println!("Loading the SA took {} ms and loading the proteins + SA took {} ms", end_loading_ms - start_loading_ms, end_loading_ms - start_reading_proteins_ms);
            // TODO: some kind of security check that the loaded database file and SA match
            sa
        }
        // build the SA
        None => {
            build_sa(&proteins.input_string, &args.construction_algorithm, args.sample_rate)?
        }
    };

    if let Some(output) = &args.output {
        // println!("storing index to file {}", output);
        write_binary(args.sample_rate, &sa, output)?;
        // println!("Index written away");
    }

    // option that only builds the tree, but does not allow for querying (easy for benchmark purposes)
    if args.build_only {
        return Ok(());
    } else if args.mode.is_none() {
        eprintln!("search mode expected!");
        std::process::exit(1);
    }

    // build the right mapping index, use box to be able to store both types in this variable
    let suffix_index_to_protein: Box<dyn SuffixToProteinIndex> =
        match args.suffix_to_protein_mapping {
            SuffixToProteinMappingStyle::Dense => {
                Box::new(DenseSuffixToProtein::new(&proteins.input_string))
            }
            SuffixToProteinMappingStyle::Sparse => {
                Box::new(SparseSuffixToProtein::new(&proteins.input_string))
            }
        };
    // println!("mapping built");

    let searcher = Searcher::new(
        sa,
        args.sample_rate,
        suffix_index_to_protein, 
        proteins,
        *taxon_id_calculator,
    );

    execute_search(&searcher, &args)?;
    Ok(())
}



/// Perform the search as set with the commandline argumentsc
fn execute_search(searcher: &Searcher, args: &Arguments) -> Result<(), Box<dyn Error>> {
    let mode = args.mode.as_ref().ok_or("No search mode provided")?;
    let cutoff = args.cutoff;
    let search_file = args
        .search_file
        .as_ref()
        .ok_or("No peptide file provided to search in the database")?;

    let start_time = get_time_ms()?;
    let lines = read_lines(search_file)?;
    let all_peptides: Vec<String> = lines.map_while(Result::ok).collect();

    // Explicitly set the number of threads to use if the commandline argument was set
    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads.get())
            .build_global()?;
    }

    all_peptides
        .par_iter()
        // calculate the results
        .map(|peptide| search_peptide(searcher, peptide, mode, cutoff))
        // output the results, collect is needed to store order so the output is in the right sequential order
        .collect::<Vec<String>>()// TODO: this collect that makes the output again sequential is possibly unneeded since we also output the corresponding peptide (but make sure this still makes the right peptide;taxon-id mapping)
        .iter()
        .enumerate()
        .for_each(|(index, res)| println!("{};{}", all_peptides[index], res));

    let end_time = get_time_ms()?;

    // output to other channel to prevent integrating it into the actual output
    eprintln!(
        "Spend {} ms to search the whole file",
        end_time - start_time
    );

    Ok(())
}

/// Executes the kind of search indicated by the commandline arguments
pub fn search_peptide(
    searcher: &Searcher,
    word: &str,
    search_mode: &SearchMode,
    cutoff: usize,
) -> String {
    let word = word.strip_suffix('\n').unwrap_or(word).to_uppercase();

    // words that are shorter than the sample rate are not searchable
    if word.len() < searcher.sample_rate as usize {
        println!("/ (word too short short for SA sample size)");
        return String::new();
    }

    match *search_mode {
        SearchMode::Match => format!("{}", searcher.search_if_match(word.as_bytes())),
        SearchMode::MinMaxBound => {
            let (found, min_bound, max_bound) = searcher.search_bounds(word.as_bytes());
            format!("{found};{min_bound};{max_bound};")
        }
        SearchMode::AllOccurrences => {
            let results = searcher.search_protein(word.as_bytes());
            let number_of_proteins = results.len();
            let peptide_length = word.len();
            format!("{peptide_length};{number_of_proteins};") // TODO: return all the matching protein strings perhaps?
        }
        SearchMode::TaxonId => {
            let (cutoff_used, suffixes) = searcher.search_matching_suffixes(word.as_bytes(), cutoff);
            let result = if cutoff_used {
                Some(1)
            } else {
                let proteins = searcher.retrieve_proteins(&suffixes);
                searcher.retrieve_lca(&proteins)
            };

            if let Some(id) = result {
                format!("{id}")
            } else {
                "/".to_string()
            }
        }
        SearchMode::Analyses => {
            let (cutoff_used, suffixes) = searcher.search_matching_suffixes(word.as_bytes(), cutoff);
            let proteins = searcher.retrieve_proteins(&suffixes);
            let result = if cutoff_used {
                Some(1)
            } else {
                searcher.retrieve_lca(&proteins)
            };

            if let Some(id) = result {
                let annotations = Searcher::get_uniprot_and_taxa_ids(&proteins);
                format!("{id};{:?}", annotations)
            } else {
                "/".to_string()
            }
        }
    }
}

/// Custom trait implemented by types that have a value that represents NULL
pub trait Nullable<T> {
    const NULL: T;

    fn is_null(&self) -> bool;
}

impl Nullable<u32> for u32 {
    const NULL: u32 = u32::MAX;

    fn is_null(&self) -> bool {
        *self == Self::NULL
    }
}
