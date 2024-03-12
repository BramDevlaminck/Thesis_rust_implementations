use clap::Parser;
use suffixarray::{Arguments, run};

fn main() {
    let args = Arguments::parse();
    if let Err(error) = run(args) {
        eprintln!("{}", error);
        std::process::exit(1);
    };
}
