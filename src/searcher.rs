use crate::cursor::CursorIterator;
use crate::read_only_cursor::ReadOnlyCursor;
use crate::tree::Tree;

pub struct Searcher<'a> {
    cursor: ReadOnlyCursor<'a>,
    original_input_string: &'a [u8],
}


impl<'a> Searcher<'a> {
    pub fn new(tree: &'a Tree, original_input_string: &'a [u8]) -> Self {
        Self {
            cursor: ReadOnlyCursor::new(tree),
            original_input_string,
        }
    }

    /// Return true as first value of the tuple if we have a valid match until the end
    /// the second value of the tuple is the index of the last current node in the arena during search
    fn find_end_node(&mut self, search_string: &[u8]) -> (bool, usize) {
        if search_string.is_empty() {
            return (true, 0);
        }
        let string_length = search_string.len();
        let mut index_in_string: usize = 0;
        let mut ret_value = self.cursor.next(search_string[0] as char, search_string);

        while ret_value == CursorIterator::Ok && index_in_string + 1 < string_length {
            index_in_string += 1;
            ret_value = self.cursor.next(search_string[index_in_string] as char, search_string);
        }

        if index_in_string == string_length - 1 && ret_value == CursorIterator::Ok {
            return (true, self.cursor.current_node_index_in_arena);
        }

        (false, self.cursor.current_node_index_in_arena)
    }


    pub fn search_protein(&mut self, search_string: &[u8]) -> Vec<String> {
        let (match_found, end_node) = self.find_end_node(search_string);
        if !match_found {
            return vec![];
        }
        let mut suffix_indices_list: Vec<usize> = vec![];
        let mut stack = vec![end_node];
        while let Some(current_node_index) = stack.pop() {
            let current_node = &self.cursor.tree.arena[current_node_index];
            if let Some(suffix_index) = current_node.suffix_index {
                suffix_indices_list.push(suffix_index);
            } else {
                current_node.children.iter().for_each(|child| {
                    if let Some(child_index) = child {
                        stack.push(*child_index);
                    }
                });
            }
        }

        self.cursor.reset();

        let dollar_u8 = b'$';
        let hashtag_u8 = b'#';
        let mut solutions_list: Vec<String> = vec![];
        suffix_indices_list.iter().for_each(|index| {
            let mut begin = *index;
            let mut end = *index;
            while begin > 0 && self.original_input_string[begin - 1] != hashtag_u8 {
                begin -= 1;
            }

            while self.original_input_string[end] != dollar_u8 && self.original_input_string[end] != hashtag_u8 {
                end += 1;
            }
            let substring = &self.original_input_string[begin..end];
            solutions_list.push(String::from_utf8_lossy(substring).into_owned())
        });
        solutions_list
    }
}