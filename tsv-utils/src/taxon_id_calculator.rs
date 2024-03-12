use umgap::{agg, rmq, taxon};
use umgap::agg::MultiThreadSafeAggregator;
use umgap::taxon::{TaxonId, TaxonList, TaxonTree};
pub struct TaxonIdCalculator {
    snapping: Vec<Option<TaxonId>>,
    aggregator: Box<dyn MultiThreadSafeAggregator>,
    taxon_list: TaxonList
}

pub trait TaxonIdVerifier: Send {

    fn taxon_id_exists(&self, id: TaxonId) -> bool;

}

pub enum AggregationMethod {
    Lca,
    LcaStar
}

impl TaxonIdCalculator {
    pub fn new(ncbi_taxonomy_fasta_file: &str, aggregation_method: AggregationMethod) -> Box<Self> {
        let taxons = taxon::read_taxa_file(ncbi_taxonomy_fasta_file).unwrap();
        let taxon_tree = TaxonTree::new(&taxons);
        let by_id = TaxonList::new(taxons);
        let snapping = taxon_tree.snapping(&by_id, true);
        
        let aggregator: Box<dyn MultiThreadSafeAggregator> = match aggregation_method { 
            AggregationMethod::Lca => Box::new(rmq::mix::MixCalculator::new(taxon_tree, 1.0)),
            AggregationMethod::LcaStar => Box::new(rmq::lca::LCACalculator::new(taxon_tree)),
        };

        Box::new(Self {
            snapping,
            aggregator,
            taxon_list: by_id
        })
    }

    /// Snaps the given taxon ID using the ncbi taxonomy that was provided to the TaxonIdCalculator.
    pub fn snap_taxon_id(&self, id: TaxonId) -> TaxonId {
        self.snapping[id].unwrap_or_else(|| panic!("Could not snap taxon with id {id}"))
    }

    /// Calculates the aggregate of a vector containing TaxonIds
    pub fn get_aggregate(&self, ids: Vec<TaxonId>) -> TaxonId {
        let count = agg::count(ids.into_iter().map(|it| (it, 1.0)));
        self.aggregator.aggregate(&count).unwrap_or_else(|_| panic!("Could not aggregate following taxon ids: {:?}", &count))
    }

}

impl TaxonIdVerifier for TaxonIdCalculator {
    fn taxon_id_exists(&self, id: TaxonId) -> bool {
        self.taxon_list.get(id).is_some()
    }
}