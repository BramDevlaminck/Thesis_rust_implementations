use umgap::taxon::TaxonId;

use tsv_utils::taxon_id_calculator::{AggregationMethod, TaxonIdCalculator, TaxonIdVerifier};

use crate::Protein;
use crate::tree::{NodeIndex, Nullable, Tree};

pub struct TreeTaxonIdCalculator {
    taxon_id_calculator: Box<TaxonIdCalculator>,
}

impl TreeTaxonIdCalculator {
    
    /// Always aggregate using LCA and not LCA*
    /// since we precalculate the aggregate efficiently in the tree
    /// and LCA* gives incorrect results the way we precalculate the taxa
    pub fn new(ncbi_taxonomy_fasta_file: &str) -> Self {
        Self {
            taxon_id_calculator: TaxonIdCalculator::new(ncbi_taxonomy_fasta_file, AggregationMethod::Lca)
        }
    }

    /// Calculates the taxon ids by only using the leaves in the tree
    pub fn calculate_taxon_ids_leaf(&self, tree: &mut Tree, proteins: &Vec<Protein>) {
        self.calculate_taxon_ids_leaf_recursive(tree, proteins, 0);
    }

    /// Calculate the taxon id for current_node by only using the leaves in the tree
    fn calculate_taxon_ids_leaf_recursive(&self, tree: &mut Tree, proteins: &Vec<Protein>, current_node_index: NodeIndex) {
        let current_node = &mut tree.arena[current_node_index];
        // we are in a leaf
        if !current_node.suffix_index.is_null() {
            current_node.taxon_id = proteins[current_node.suffix_index].id;
            return;
        }

        let taxon_id = self.get_aggregate(Self::get_taxon_id_leaves_under_node(tree, proteins, current_node_index));

        let current_node = &mut tree.arena[current_node_index];
        current_node.taxon_id = taxon_id;

        for child in current_node.children {
            if !child.is_null() {
                self.calculate_taxon_ids_leaf_recursive(tree, proteins, child);
            }
        }
    }
    
    /// Fill in all the taxon ids for the complete tree
    pub fn calculate_taxon_ids(&self, tree: &mut Tree, proteins: &Vec<Protein>) {
        self.recursive_calculate_taxon_ids(tree, proteins, 0);
    }

    /// The recursive function called by `calculate_taxon_ids`
    /// This function traverses the tree post order and fills in the taxon ids
    fn recursive_calculate_taxon_ids(&self, tree: &mut Tree, proteins: &Vec<Protein>, current_node_index: NodeIndex) -> TaxonId {
        let current_node = &mut tree.arena[current_node_index];
        // we are in a leaf
        if !current_node.suffix_index.is_null() {
            current_node.taxon_id = proteins[current_node.suffix_index].id;
            return current_node.taxon_id;
        }

        let mut taxon_ids = vec![];
        for child in current_node.children {
            if !child.is_null() {
                taxon_ids.push(self.recursive_calculate_taxon_ids(tree, proteins, child));
            }
        }

        let taxon_id = self.get_aggregate(taxon_ids);

        let current_node = &mut tree.arena[current_node_index];
        current_node.taxon_id = taxon_id;
        taxon_id
    }

    /// returns the taxon ids of all the leaves that are under current_node
    fn get_taxon_id_leaves_under_node(tree: &Tree, proteins: &Vec<Protein>, current_node_index: NodeIndex) -> Vec<TaxonId> {
        let mut ids: Vec<TaxonId> = vec![];
        Self::get_taxon_id_leaves_recursive(tree, proteins, &mut ids, current_node_index);
        ids
    }

    /// recursive function that adds al the taxon ids of the leaves under current_node to taxon_ids
    fn get_taxon_id_leaves_recursive(tree: &Tree, proteins: &Vec<Protein>, taxon_ids: &mut Vec<TaxonId>, current_node_index: NodeIndex) {
        let current_node = &tree.arena[current_node_index];
        // we are in a leaf
        if !current_node.suffix_index.is_null() {
            taxon_ids.push(proteins[current_node.suffix_index].id);
            return;
        }

        for child in current_node.children {
            if !child.is_null() {
                Self::get_taxon_id_leaves_recursive(tree, proteins, taxon_ids, child);
            }
        }
    }

    /// Helper function to print out the taxon ids of the tree in pre-order traversal
    pub fn output_all_taxon_ids(tree: &Tree) {
        Self::output_all_taxon_ids_recursive(tree, 0)
    }

    fn output_all_taxon_ids_recursive(tree: &Tree, current_node_index: TaxonId) {
        let current_node = &tree.arena[current_node_index];

        println!("{}", current_node.taxon_id);
        print!("children: ");

        // we are in a leaf
        if !current_node.suffix_index.is_null() {
            println!("/"); // indicate that there are no children since this is a leaf
            return;
        }

        for child in current_node.children {
            if !child.is_null() {
                print!("{};", tree.arena[child].taxon_id);
            }
        }
        println!();

        for child in current_node.children {
            if !child.is_null() {
                Self::output_all_taxon_ids_recursive(tree, child);
            }
        }
    }

    /// Snaps the given taxon ID using the ncbi taxonomy that was provided to the TaxonIdCalculator.
    pub fn snap_taxon_id(&self, id: TaxonId) -> TaxonId {
        self.taxon_id_calculator.snap_taxon_id(id)
    }

    /// Calculates the aggregate of a vector containing TaxonIds
    fn get_aggregate(&self, ids: Vec<TaxonId>) -> TaxonId {
        self.taxon_id_calculator.get_aggregate(ids)
    }
}

impl TaxonIdVerifier for TreeTaxonIdCalculator {
    fn taxon_id_exists(&self, id: TaxonId) -> bool {
        self.taxon_id_calculator.taxon_id_exists(id)
    }
}


#[cfg(test)]
mod test {
    use crate::Protein;
    use crate::tree_taxon_id_calculator::TreeTaxonIdCalculator;
    use crate::tree::{MAX_CHILDREN, Node, NodeIndex, Nullable, Range, Tree};

    #[test]
    fn test_calculate_taxon_ids() {
        // the tree structure we are building in this test
        // with the expected taxon id between parentheses under the node id
        //              0
        //             (1)
        //             /|\
        //            / | \
        //           /  |  \
        //          /   |   \
        //         1    2    3
        //        (6)  (20) (2)
        //        / \
        //       /   \
        //      4     5
        //     (10)  (9)

        let test_taxonomy_file = "../small_taxonomy.tsv";
        let mut tree = Tree {
            arena: vec![
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
            ]
        };

        // set some child structure
        tree.arena[0].children[0] = 1;
        tree.arena[0].children[3] = 2;
        tree.arena[0].children[5] = 3;
        tree.arena[1].children[1] = 4;
        tree.arena[1].children[5] = 5;

        // set the parents to match what we set in the children
        tree.arena[1].parent = 0;
        tree.arena[2].parent = 0;
        tree.arena[3].parent = 0;
        tree.arena[4].parent = 1;
        tree.arena[5].parent = 1;

        // set suffix ids in the leaves
        tree.arena[4].suffix_index = 0;
        tree.arena[5].suffix_index = 1;
        tree.arena[2].suffix_index = 2;
        tree.arena[3].suffix_index = 3;


        let proteins = vec![
            Protein {
                uniprot_id: "Q0".to_string(),
                sequence: (0, 3),
                id: 10,
            },
            Protein {
                uniprot_id: "Q1".to_string(),
                sequence: (4, 7),
                id: 9,
            },
            Protein {
                uniprot_id: "Q2".to_string(),
                sequence: (8, 11),
                id: 20,
            },
            Protein {
                uniprot_id: "Q3".to_string(),
                sequence: (12, 13),
                id: 2,
            },
        ];

        TreeTaxonIdCalculator::new(test_taxonomy_file).calculate_taxon_ids(&mut tree, &proteins);

        assert_eq!(tree.arena[0].taxon_id, 1);
        assert_eq!(tree.arena[1].taxon_id, 6);
        assert_eq!(tree.arena[2].taxon_id, 20);
        assert_eq!(tree.arena[3].taxon_id, 2);
        assert_eq!(tree.arena[4].taxon_id, 10);
        assert_eq!(tree.arena[5].taxon_id, 9);
    }
}