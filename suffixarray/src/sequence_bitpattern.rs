
/// bitvector indicating which combinations of I and L we have already visited
/// we use this as if we were counting binary and I = 0 and L = 1.
/// e.g. the string IACLUIL has the pattern "0101" (we ignore every character that is not I or L),
/// in that representation, the least significant bit is to the left, but in a computer it is to the right, so we have "1010"
/// we know that this as an integer is "10", this means that the 10th bit should be set to a 1 in the visited_pattern vector
pub struct SequenceBitPattern {
    pattern: Vec<u64>
}

impl SequenceBitPattern {
    
    pub fn new(il_locations: &Vec<usize>) -> Self {
        Self {
            pattern: vec![0; (2_usize.pow(il_locations.len() as u32) / 64) + 1]
        }
    }
    
    /// Translates a sequence (containing I's and L's) to a bit pattern, with the limit that there are not allowed to be more than 64 il_locations
    /// (this will already explode too much and make the pattern unusable)
    fn il_sequence_to_bit_number(sequence: &Vec<u8>, il_locations: &Vec<usize>) -> u64 {
        let mut pattern: u64 = 0;

        for (index, &location) in il_locations.iter().enumerate() {
            if sequence[location] == b'L' {
                pattern |= 1 << index
            }
        }

        pattern
    }
    
    pub fn check_if_contains_and_add(&mut self, sequence: &Vec<u8>, il_locations: &Vec<usize>) -> bool {
        let bit_nr = Self::il_sequence_to_bit_number(sequence, il_locations);
        
        let array_index = (bit_nr / 64) as usize;
        let bit_index = bit_nr % 64;
        
        let exists = (self.pattern[array_index] & (1 << bit_index)) != 0;

        self.pattern[array_index] |= 1 << bit_index;
        exists
    }
}


#[cfg(test)]
mod tests {
    use crate::sequence_bitpattern::SequenceBitPattern;

    #[test]
    fn sequence_to_bit_pattern_test() {
        let sequence = "ILACILAAICIL".to_string().into_bytes();
        let il_locations = vec![0, 1, 4, 5, 8, 10, 11];
        let expected_output = 74; // 0101001 when we see the sequence, we need to reverse it to get the integer value 1001010
        assert_eq!(SequenceBitPattern::il_sequence_to_bit_number(&sequence, &il_locations), expected_output);

    } 
    
    #[test]
    fn check_if_contains_and_add_test() {
        let mut pattern = SequenceBitPattern::new(&vec![1, 3, 5, 6]);
        let sequence = "IILAJCL".to_string().into_bytes();
        let il_locations = &vec![0, 1, 2, 6];
        assert_eq!(pattern.check_if_contains_and_add(&sequence, &il_locations), false);
        // second time we add should return that the pattern already exists
        assert_eq!(pattern.check_if_contains_and_add(&sequence, &il_locations), true);
        // check internal state
        assert_eq!(pattern.pattern.len(), 1);
        assert_eq!(pattern.pattern[0], 4096)
    }
    
}
