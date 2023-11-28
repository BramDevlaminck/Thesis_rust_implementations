use umgap::taxon::TaxonId;
use crate::Protein;
use crate::search_cursor::SearchCursor;
use crate::taxon_id_calculator::TaxonIdCalculator;
use crate::tree::{Node, Nullable, Tree};

pub struct Searcher<'a> {
    cursor: SearchCursor<'a>,
    original_input_string: &'a [u8],
    proteins: &'a Vec<Protein>,
    taxon_id_calculator: &'a TaxonIdCalculator, // used to snap taxon_id
}


impl<'a> Searcher<'a> {
    pub fn new(tree: &'a Tree, original_input_string: &'a [u8], proteins: &'a Vec<Protein>, taxon_id_calculator: &'a TaxonIdCalculator) -> Self {
        Self {
            cursor: SearchCursor::new(tree),
            original_input_string,
            proteins,
            taxon_id_calculator,
        }
    }

    /// Return true as first value of the tuple if we have a valid match until the end
    /// the second value of the tuple is the index of the last current node in the arena during search
    fn find_end_node(&mut self, search_string: &[u8]) -> (bool, &'a Node) {
        if search_string.is_empty() {
            return (true, &self.cursor.tree.arena[0]);
        }
        let string_length = search_string.len();
        let mut index_in_string: usize = 0;

        while self.cursor.next(search_string[index_in_string], self.original_input_string).is_some() {
            index_in_string += 1;
            if index_in_string == string_length {
                let end_node = self.cursor.current_node;
                self.cursor.reset(); // prepare cursor for next search
                return (true, end_node);
            }
        }

        let end_node = self.cursor.current_node;
        self.cursor.reset(); // prepare cursor for next search

        (false, end_node)
    }


    pub fn search_protein(&mut self, search_string: &[u8]) -> Vec<&Protein> {
        let suffix_indices_list = self.find_all_suffix_indices(search_string);

        let mut solutions_list: Vec<&Protein> = vec![];
        suffix_indices_list.iter().for_each(|index| {
            solutions_list.push(&self.proteins[*index]);
        });
        solutions_list
    }

    pub fn find_all_suffix_indices(&mut self, search_string: &[u8]) -> Vec<usize> {
        let (match_found, end_node) = self.find_end_node(search_string);
        if !match_found {
            return vec![];
        }
        let mut suffix_indices_list: Vec<usize> = vec![];
        let mut stack = vec![end_node];
        while let Some(current_node) = stack.pop() {
            if !current_node.suffix_index.is_null() {
                suffix_indices_list.push(current_node.suffix_index);
            } else {
                current_node.children.iter().for_each(|&child| {
                    if !child.is_null() {
                        stack.push(&self.cursor.tree.arena[child]);
                    }
                });
            }
        }
        suffix_indices_list
    }

    pub fn search_if_match(&mut self, search_string: &[u8]) -> bool {
        self.find_end_node(search_string).0
    }

    pub fn search_taxon_id(&mut self, search_string: &[u8]) -> Option<TaxonId> {
        let (match_found, end_node) = self.find_end_node(search_string);
        if match_found {
            Some(self.taxon_id_calculator.snap_taxon_id(end_node.taxon_id))
        } else {
            None
        }
    }
}