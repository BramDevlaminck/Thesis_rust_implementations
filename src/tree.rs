use crate::tree_builder::TreeBuilder;

pub const MAX_CHILDREN: usize = 28;

/// Custom trait implemented by types that have a value that represents NULL
pub trait Nullable<T> {
    const NULL: T;

    fn is_null(&self) -> bool;
}

/// Type that represents the index of a node in the arena part of the tree
pub type NodeIndex = usize;

impl Nullable<NodeIndex> for NodeIndex {
    /// Use usize::MAX as NULL value since this will in practice never be reached.
    /// It is not possible to create 2^64-1 nodes (on a 64-bit machine). 
    /// This would simply never fit in memory
    const NULL: NodeIndex = usize::MAX;

    fn is_null(&self) -> bool {
        *self == NodeIndex::NULL
    }
}


#[derive(Debug)]
pub struct Tree {
    pub arena: Vec<Node>,
}

impl Tree {
    pub fn new(data: &str, builder: impl TreeBuilder) -> Self {
        builder.build(
            data,
            Tree {
                arena: vec![Node::create_root()],
            },
        )
    }
}

#[derive(Debug)]
pub struct Node {
    pub range: Range,
    pub children: [NodeIndex; MAX_CHILDREN],
    pub parent: NodeIndex,
    pub link: NodeIndex,
    pub suffix_index: NodeIndex,
}

impl Node {
    pub fn create_root() -> Self {
        Node {
            range: Range::new(0, 0),
            children: [NodeIndex::NULL; MAX_CHILDREN],
            parent: NodeIndex::NULL,
            link: NodeIndex::NULL,
            suffix_index: NodeIndex::NULL,
        }
    }

    /// Returns a tuple that contains the index of the new node in the arena and a reference to that node
    pub fn new(range: Range, parent: NodeIndex, children: [NodeIndex; MAX_CHILDREN], link: NodeIndex, suffix_index: usize) -> Node {
        Node {
            range,
            children,
            parent,
            link,
            suffix_index,
        }
    }

    pub fn new_with_child_tuples(range: Range, parent: NodeIndex, children_tuples: Vec<(u8, NodeIndex)>, link: NodeIndex, suffix_index: usize) -> Node {
        let mut node = Node {
            range,
            children: [NodeIndex::NULL; MAX_CHILDREN],
            parent,
            link,
            suffix_index,
        };
        children_tuples.iter().for_each(|(char, child)| node.add_child(*char, *child));
        node
    }

    fn char_to_child_index(character: u8) -> usize {
        if character == b'#' {
            26
        } else if character == b'$' {
            27
        } else {
            character as usize - 65 // 65 is 'A' in ascii
        }
    }

    pub fn add_child(&mut self, character: u8, child: NodeIndex) {
        self.children[Self::char_to_child_index(character)] = child;
    }

    // TODO: mogelijks rekening houden met feit dat er tekens ingeput kunnen worden die niet in de array zitten?
    //  (bv ?*^), dit gaat er nu voor zorgen dat alles crasht als we zo'n teken mee geven
    pub fn get_child(&self, character: u8) -> NodeIndex {
        self.children[Self::char_to_child_index(character)]
    }

    pub fn set_new_children(&mut self, new_children: Vec<(u8, NodeIndex)>) {
        self.children = [NodeIndex::NULL; MAX_CHILDREN];
        new_children.iter().for_each(|(character, child)| self.add_child(*character, *child));
    }
}

#[derive(Debug)]
pub struct Range {
    pub start: usize,
    pub end: usize,
}

impl Range {
    pub fn new(start: usize, end: usize) -> Self {
        Range { start, end }
    }
    pub fn length(&self) -> usize {
        self.end - self.start
    }
}