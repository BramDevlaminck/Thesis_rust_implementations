use std::cmp::{min};
use umgap::taxon::TaxonId;
use tsv_utils::{Protein, Proteins};
use tsv_utils::taxon_id_calculator::TaxonIdCalculator;
use crate::Nullable;
use crate::suffix_to_protein_index::SuffixToProteinIndex;

pub struct Searcher<'a> {
    sa: &'a Vec<i64>,
    pub sample_rate: u8,
    suffix_index_to_protein: &'a dyn SuffixToProteinIndex,
    proteins: &'a Proteins,
    taxon_id_calculator: &'a TaxonIdCalculator
}

impl <'a> Searcher<'a> {

    pub fn new(sa: &'a Vec<i64>, sample_rate: u8, suffix_index_to_protein: &'a dyn SuffixToProteinIndex, proteins: &'a Proteins, taxon_id_calculator: &'a TaxonIdCalculator) -> Self {
        Self {
            sa,
            sample_rate,
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
            && index < self.proteins.input_string.len()
            && search_string[index_in_search_string] == self.proteins.input_string[index] {
            index += 1;
            index_in_search_string += 1;
        }
        // check if match found OR current search string is smaller lexicographically (and the empty search string should not be found)
        if !search_string.is_empty()
            && (
            index_in_search_string == search_string.len()
            || (index < self.proteins.input_string.len() && compare_fn(search_string[index_in_search_string] as usize, self.proteins.input_string[index] as usize))
        ) {
            is_cond_or_equal = true;
        }
        (is_cond_or_equal, index_in_search_string)
    }

    // note: $ cannot be found since this is the first value in the array and the way the bounds are changed with center, index 0 is never checked
    fn binary_search_min_match(&self, search_string: &[u8]) -> (bool, usize) {
        let mut left: usize = 0;
        let mut right: usize = self.sa.len();
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

        // handle edge case to search at index 0
        if right == 1 && left == 0 {
            let (retval, lcp_center) = self.compare(search_string, self.sa[0], min(lcp_left, lcp_right), |a, b| a < b);
            found |= lcp_center == search_string.len();
            if retval {
                right = 0;
            }
        }
        (found, right)
    }

    fn binary_search_max_match(&self, search_string: &[u8]) -> (bool, usize) {
        let mut left: usize = 0;
        let mut right: usize = self.sa.len();
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

    fn binary_search_match(&self, search_string: &[u8]) -> bool {
        let mut left: usize = 0;
        let mut right: usize = self.sa.len();
        let mut lcp_left: usize = 0;
        let mut lcp_right: usize = 0;
        let mut found = false;

        // repeat until search window is minimum size OR we matched the whole search string last iteration
        while !found && right - left > 1 {
            let center = (left + right) / 2;
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center) = self.compare(search_string, self.sa[center], skip, |a, b| a < b);
            found = lcp_center == search_string.len();
            if retval {
                right = center;
                lcp_right = lcp_center;
            } else {
                left = center;
                lcp_left = lcp_center;
            }
        }

        // handle edge case to search at index 0
        if !found && right == 1 && left == 0 {
            let (_, lcp_center) = self.compare(search_string, self.sa[0], min(lcp_left, lcp_right), |a, b| a < b);
            found = lcp_center == search_string.len();
        }

        found
    }

    pub fn search_bounds(&self, search_string: &[u8]) -> (bool, usize, usize) {
        let (found, min_bound) = self.binary_search_min_match(search_string);
        if !found {
            return (false, 0, self.sa.len());
        }
        let (_, max_bound) = self.binary_search_max_match(search_string);

        (true, min_bound, max_bound + 1) // add 1 to max bound to have exclusive end
    }

    pub fn search_if_match(&self, search_string: &[u8]) -> bool {
        // case where SA is not sparse
        if self.sample_rate == 1 {
            return self.binary_search_match(search_string)
        }

        for skip in 0..self.sample_rate as usize {
            let (found, min_bound, max_bound) = self.search_bounds(&search_string[skip..]);
            // if the shorter part is matched, see if what goes before the matched suffix matches the unmatched part of the prefix
            if found {
                let unmatched_prefix = &search_string[..skip];
                // try all the partially matched suffixes
                for sa_index in min_bound..max_bound {
                    let suffix = self.sa[sa_index] as usize;
                    if suffix >= skip && unmatched_prefix == &self.proteins.input_string[suffix - skip..suffix] {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// Search all the suffixes that search string matches with
    #[inline]
    pub fn search_matching_suffixes(&self, search_string: &[u8], max_matches: usize) -> Vec<i64> {
        let mut matching_suffixes: Vec<i64> = vec![];
        let mut skip: usize = 0;
        while matching_suffixes.len() < max_matches && skip < self.sample_rate as usize {
            let (found, min_bound, max_bound) = self.search_bounds(&search_string[skip..]);
            // if the shorter part is matched, see if what goes before the matched suffix matches the unmatched part of the prefix
            if found {
                let unmatched_prefix = &search_string[..skip];
                // try all the partially matched suffixes and store the matching suffixes in an array (stop when our max number of matches is reached)
                let mut sa_index = min_bound;
                while sa_index < max_bound && matching_suffixes.len() < max_matches {
                    let suffix = self.sa[sa_index] as usize;
                    // if skip is 0, then we already checked the complete match during bound search, otherwise check if the skipped part also matches
                    if skip == 0 || (suffix >= skip && unmatched_prefix == &self.proteins.input_string[suffix - skip..suffix]) {
                        matching_suffixes.push((suffix - skip) as i64);
                    }
                    sa_index += 1;
                }
            }
            skip += 1;
        }
        matching_suffixes
    }

    /// get all the proteins matching with the given suffixes
    #[inline]
    pub fn retrieve_proteins(&self, suffixes: &Vec<i64>) -> Vec<&Protein> {
        let mut res = vec![];
        for &suffix in suffixes {
            let protein_index = self.suffix_index_to_protein.suffix_to_protein(suffix);
            if !protein_index.is_null() {
                res.push(&self.proteins.proteins[protein_index as usize]);
            }
        }
        res
    }

    /// Search all the Proteins that a given search_string matches with
    pub fn search_protein(&self, search_string: &[u8]) -> Vec<&Protein> {
        let matching_suffixes = self.search_matching_suffixes(search_string, usize::MAX);
        self.retrieve_proteins(&matching_suffixes)
    }

    #[inline]
    pub fn retrieve_taxon_id(&self, proteins: &[&Protein]) -> Option<TaxonId> {
        let taxon_ids: Vec<TaxonId> = proteins.iter().map(|prot| prot.id).collect();
        match taxon_ids.is_empty() {
            true => None,
            false => Some(self.taxon_id_calculator.snap_taxon_id(self.taxon_id_calculator.get_aggregate(taxon_ids)))
        }
    }

    pub fn search_taxon_id(&self, search_string: &[u8]) -> Option<TaxonId> {
        self.retrieve_taxon_id(&self.search_protein(search_string))
    }

}