use std::error::Error;
use clap::Parser;
use suffixarray::{Arguments, run};

fn main() -> Result<(), Box<dyn Error>>{
    let args = Arguments::parse();
    run(args)?;
    Ok(())
}
