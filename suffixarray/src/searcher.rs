use std::cmp::{min};
use umgap::taxon::TaxonId;
use tsv_utils::Protein;
use tsv_utils::taxon_id_calculator::TaxonIdCalculator;

pub struct Searcher<'a> {
    original_input_string: &'a [u8],
    sa: &'a Vec<i64>,
    suffix_index_to_protein: &'a Vec<Option<u32>>,
    proteins: &'a Vec<Protein>,
    taxon_id_calculator: &'a TaxonIdCalculator
}

impl <'a> Searcher<'a> {

    pub fn new(original_input_string: &'a [u8], sa: &'a Vec<i64>, suffix_index_to_protein: &'a Vec<Option<u32>>, proteins: &'a Vec<Protein>, taxon_id_calculator: &'a TaxonIdCalculator) -> Self {
        Self {
            original_input_string,
            sa,
            suffix_index_to_protein,
            proteins,
            taxon_id_calculator
        }
    }

    fn compare(&self, search_string: &[u8], suffix: i64, skip: usize, compare_fn: fn(usize, usize) -> bool) -> (bool, usize) {
        let mut index = (suffix as usize) + skip;
        let mut index_in_search_string = skip;
        let mut is_cond_or_equal = false;
        // match as long as possible
        while index_in_search_string < search_string.len()
            && search_string[index_in_search_string] == self.original_input_string[index] {
            index += 1;
            index_in_search_string += 1;
        }
        // check if match found OR current search string is smaller lexicographically
        if index_in_search_string == search_string.len() || compare_fn(search_string[index_in_search_string] as usize, self.original_input_string[index] as usize) {
            is_cond_or_equal = true;
        }
        (is_cond_or_equal, index_in_search_string)
    }

    // note: $ cannot be found since this is the first value in the array and the way the bounds are changed with center, index 0 is never checked
    fn binary_search_min_match(&self, search_string: &[u8]) -> (bool, usize) {
        let mut left: usize = 0;
        let mut right: usize = self.original_input_string.len();
        let mut lcp_left: usize = 0;
        let mut lcp_right: usize = 0;
        let mut found = false;

        // repeat until search window is minimum size OR we matched the whole search string last iteration
        while right - left > 1 {
            let center = (left + right) / 2;
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center) = self.compare(search_string, self.sa[center], skip, |a, b| a < b);
            found |= lcp_center == search_string.len();
            if retval {
                right = center;
                lcp_right = lcp_center;
            } else {
                left = center;
                lcp_left = lcp_center;
            }
        }
        (found, right)
    }

    // note: $ cannot be found since this is the first value in the array and the way the bounds are changed with center, index 0 is never checked
    fn binary_search_max_match(&self, search_string: &[u8]) -> (bool, usize) {
        let mut left: usize = 0;
        let mut right: usize = self.original_input_string.len();
        let mut lcp_left: usize = 0;
        let mut lcp_right: usize = 0;
        let mut found = false;

        // repeat until search window is minimum size OR we matched the whole search string last iteration
        while right - left > 1 {
            let center = (left + right) / 2;
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center) = self.compare(search_string, self.sa[center], skip, |a, b| a > b);
            found |= lcp_center == search_string.len();
            if retval {
                left = center;
                lcp_left = lcp_center;
            } else {
                right = center;
                lcp_right = lcp_center;
            }
        }
        (found, left)
    }

    // note: $ cannot be found since this is the first value in the array and the way the bounds are changed with center, index 0 is never checked
    fn binary_search_match(&self, search_string: &[u8]) -> bool {
        let mut left: usize = 0;
        let mut right: usize = self.original_input_string.len();
        let mut lcp_left: usize = 0;
        let mut lcp_right: usize = 0;
        let mut found = false;

        // repeat until search window is minimum size OR we matched the whole search string last iteration
        while !found && right - left > 1 {
            let center = (left + right) / 2;
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center) = self.compare(search_string, self.sa[center], skip, |a, b| a < b);
            found |= lcp_center == search_string.len();
            if retval {
                right = center;
                lcp_right = lcp_center;
            } else {
                left = center;
                lcp_left = lcp_center;
            }
        }
        found
    }

    pub fn search_bounds(&self, search_string: &[u8]) -> (bool, usize, usize) {
        let (found, min_bound) = self.binary_search_min_match(search_string);
        if !found {
            return (false, 0, self.original_input_string.len());
        }
        let (_, max_bound) = self.binary_search_max_match(search_string);

        (true, min_bound, max_bound + 1) // add 1 to max bound to have exclusive end
    }

    pub fn search_if_match(&self, search_string: &[u8]) -> bool {
        self.binary_search_match(search_string)
    }

    pub fn search_protein(&self, search_string: &[u8]) -> Vec<&Protein> {
        let (found, min_bound, max_bound) = self.search_bounds(search_string);
        let mut res = vec![];
        if !found {
            return res;
        }
        for i in min_bound..max_bound {
            let suffix_index = self.sa[i] as usize;
            if let Some(protein_index) = self.suffix_index_to_protein[suffix_index] {
                res.push(&self.proteins[protein_index as usize]);
            }
        }
        res
    }

    pub fn search_taxon_id(&self, search_string: &[u8]) -> Option<TaxonId> {
        let taxon_ids: Vec<TaxonId> = self.search_protein(search_string).into_iter().map(|prot| prot.id).collect();
        match taxon_ids.is_empty() {
            true => None,
            false => Some(self.taxon_id_calculator.snap_taxon_id(self.taxon_id_calculator.get_aggregate(taxon_ids)))
        }

    }

}