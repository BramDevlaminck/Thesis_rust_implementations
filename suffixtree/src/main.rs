use clap::Parser;

use suffixtree::{Arguments, run};

fn main() {
    let args = Arguments::parse();
    run(args);
}
