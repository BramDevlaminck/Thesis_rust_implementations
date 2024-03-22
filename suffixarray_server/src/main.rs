use std::error::Error;
use std::sync::Arc;

use axum::{http::StatusCode, Json, Router};
use axum::extract::{DefaultBodyLimit, State};
use axum::routing::{get, post};
use clap::Parser;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use suffixarray::searcher::Searcher;
use suffixarray::suffix_to_protein_index::SparseSuffixToProtein;
use suffixarray_builder::binary::load_binary;
use tsv_utils::get_proteins_from_database_file;
use tsv_utils::taxon_id_calculator::{AggregationMethod, TaxonIdCalculator};

#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be concatenated using a `#`.
    #[arg(short, long)]
    database_file: String,
    #[arg(short, long)]
    index_file: String,
    #[arg(short, long)]
    /// The taxonomy to be used as a tsv file. This is a preprocessed version of the NCBI taxonomy.
    taxonomy: String,
}

/// Function used by serde to place a default value in the cutoff field of the input
fn default_cutoff() -> usize {
    10000
}

#[derive(Debug, Deserialize, Serialize)]
#[allow(non_snake_case)]
struct InputData {
    peptides: Vec<String>,
    #[serde(default = "default_cutoff")] // default value is 10000
    cutoff: usize,
    #[serde(default = "bool::default")] // default value is false
    equalize_I_and_L: bool,
}

#[derive(Debug, Serialize)]
struct OutputData {
    result: Vec<SearchResult>,
}

#[derive(Debug, Serialize)]
struct SearchResult {
    sequence: String,
    lca: usize,
    taxa: Vec<usize>,
    uniprot_accessions: Vec<String>,
    cutoff_used: bool
}

// basic handler that responds with a static string
async fn root() -> &'static str {
    "Server is online"
}

/// Executes the search of 1 peptide
pub fn search_peptide(
    searcher: &Searcher,
    peptide: &str,
    cutoff: usize,
    equalize_i_and_l: bool
) -> Option<(usize, bool, Vec<String>, Vec<usize>)> {
    let mut peptide = peptide.to_uppercase();
    if equalize_i_and_l { // translate L to an I if we equalize them
        peptide = peptide.chars()
            .map(|character| match character {
                'L' => 'I',
                _ => character,
            })
            .collect()
    }

    // words that are shorter than the sample rate are not searchable
    if peptide.len() < searcher.sample_rate as usize {
        return None;
    }

    let (cutoff_used, suffixes) = searcher.search_matching_suffixes(peptide.as_bytes(), cutoff, equalize_i_and_l);
    let proteins = searcher.retrieve_proteins(&suffixes);
    let (uniprot_acc, taxa) = Searcher::get_uniprot_and_taxa_ids(&proteins);
    if cutoff_used {
        Some((1, cutoff_used, uniprot_acc, taxa))
    } else {
        let id = searcher.retrieve_lca(&proteins);
        id.map(|id_unwrapped| (id_unwrapped, cutoff_used, uniprot_acc, taxa))
    }
}

/// Search all the peptides in the given data json using the searcher.
async fn search(
    State(searcher): State<Arc<Searcher>>,
    data: Json<InputData>,
) -> Result<Json<OutputData>, StatusCode> {
    
    
    let res: Vec<SearchResult> = data
        .peptides
        .par_iter()
        // calculate the results
        .map(|peptide| search_peptide(&searcher, peptide, data.cutoff, data.equalize_I_and_L))
        .enumerate()
        // remove the peptides that did not match any proteins
        .filter(|(_, search_result)| search_result.is_some())
        // transform search results into output
        .map(|(index, search_results)| {
            let (taxon_id, cutoff_used, uniprot_acc, taxa) = search_results.unwrap();
            SearchResult {
                sequence: data.peptides[index].clone(),
                lca: taxon_id,
                taxa,
                uniprot_accessions: uniprot_acc,
                cutoff_used
            }
        })
        .collect();

    // Create output JSON
    let output_data = OutputData { result: res };

    Ok(Json(output_data))
}

#[tokio::main]
async fn main() {
    let args = Arguments::parse();
    if let Err(err) = start_server(args).await {
        eprintln!("{}", err);
        std::process::exit(1);
    }
}

async fn start_server(args: Arguments) -> Result<(), Box<dyn Error>> {
    let Arguments {
        database_file,
        index_file,
        taxonomy,
    } = args;
    
    let (sample_rate, sa) = load_binary(&index_file)?;

    let taxon_id_calculator = TaxonIdCalculator::new(&taxonomy, AggregationMethod::LcaStar);
    let proteins = get_proteins_from_database_file(&database_file, &*taxon_id_calculator)?;
    let suffix_index_to_protein = Box::new(SparseSuffixToProtein::new(&proteins.input_string));

    let searcher = Arc::new(Searcher::new(
        sa,
        sample_rate,
        suffix_index_to_protein,
        proteins,
        *taxon_id_calculator,
    ));

    // build our application with a route
    let app = Router::new()
        // `GET /` goes to `root`
        .route("/", get(root))
        // `POST /search` goes to `calculate` and set max payload size to 5 MB
        .route("/search", post(search)).layer(DefaultBodyLimit::max(5 * 10_usize.pow(6)))
        .with_state(searcher);

    let listener = tokio::net::TcpListener::bind("0.0.0.0:3000").await?;
    println!("server is ready...");
    axum::serve(listener, app).await?;

    Ok(())
}
