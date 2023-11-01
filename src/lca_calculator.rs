use umgap::{agg, rmq::lca::LCACalculator, taxon};
use umgap::agg::Aggregator;
use umgap::taxon::TaxonId;

use crate::Protein;
use crate::tree::{NodeIndex, Nullable, Tree};

pub fn calculate_lcas(tree: &mut Tree, ncbi_taxonomy_fasta_file: &str, proteins: &[Protein]) {
    let taxons = taxon::read_taxa_file(ncbi_taxonomy_fasta_file).unwrap();
    let taxon_tree = taxon::TaxonTree::new(&taxons);
    let by_id = taxon::TaxonList::new(taxons);
    let snapping = taxon_tree.snapping(&by_id, false);

    let calculator = LCACalculator::new(taxon_tree);

    stack_calculate_lcas(tree, proteins, &calculator, &snapping);
}

fn stack_calculate_lcas(tree: &mut Tree, proteins: &[Protein], calculator: &dyn Aggregator, snapping: &[Option<TaxonId>]) {
    let mut stack: Vec<(NodeIndex, bool)> = vec![(0, false)];
    let mut stack_calculated_children: Vec<Vec<TaxonId>> = vec![vec![]];
    while let Some((node_index, visited)) = stack.pop() {
        // children already visited => calculate lca* for this node and "return from the recursion"
        if visited {
            let taxon_ids = stack_calculated_children.pop().unwrap();
            let new_taxon_id = get_lca(calculator, snapping, taxon_ids);
            let current_node = &mut tree.arena[node_index];
            current_node.taxon_id = new_taxon_id;
            stack_calculated_children.last_mut().unwrap().push(new_taxon_id);
            continue;
        }

        // base case for leaves
        let current_node = &mut tree.arena[node_index];
        if !current_node.suffix_index.is_null() {
            let taxon_id = snap_taxon_id(snapping, proteins[current_node.suffix_index].id);
            current_node.taxon_id = taxon_id;
            stack_calculated_children.last_mut().unwrap().push(taxon_id);
            continue;
        }

        // visit the children
        stack_calculated_children.push(vec![]);
        stack.push((node_index, true));
        for child in tree.arena[node_index].children {
            if !child.is_null() {
                stack.push((child, false));
            }
        }
    }
}

fn snap_taxon_id(snapping: &[Option<TaxonId>], id: TaxonId) -> TaxonId {
    snapping[id].unwrap_or_else(|| panic!("Could not snap taxon id {id}"))
}


fn get_lca(calculator: &dyn Aggregator, snapping: &[Option<TaxonId>], ids: Vec<TaxonId>) -> TaxonId {
    let count = agg::count(ids.into_iter().filter(|&id| id != 0).map(|it| (it, 1.0)));
    // let counts = agg::filter(counts, args.lower_bound); TODO: used in umgap, but probably not needed here?
    let aggregate = calculator.aggregate(&count).unwrap_or_else(|_| panic!("Could not aggregate following taxon ids: {:?}", &count));
    snap_taxon_id(snapping, aggregate)
}

#[cfg(test)]
mod test {
    use crate::lca_calculator::calculate_lcas;
    use crate::Protein;
    use crate::tree::{MAX_CHILDREN, Node, NodeIndex, Nullable, Range, Tree};

    #[test]
    fn test_calculate_lcas() {
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

        let test_taxonomy_file = "small_taxonomy.tsv";
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
                sequence: "ABC".to_string(),
                id: 10,
            },
            Protein {
                sequence: "XYZ".to_string(),
                id: 9,
            },
            Protein {
                sequence: "XAZ".to_string(),
                id: 20,
            },
            Protein {
                sequence: "W".to_string(),
                id: 2,
            },
        ];

        calculate_lcas(&mut tree, test_taxonomy_file, &proteins);

        assert_eq!(tree.arena[0].taxon_id, 1);
        assert_eq!(tree.arena[1].taxon_id, 6);
        assert_eq!(tree.arena[2].taxon_id, 20);
        assert_eq!(tree.arena[3].taxon_id, 2);
        assert_eq!(tree.arena[4].taxon_id, 10);
        assert_eq!(tree.arena[5].taxon_id, 9);
    }
}