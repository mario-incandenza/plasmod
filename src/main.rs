use anyhow::Result;
use clap::Parser;
use plasmod::*;
use std::path::PathBuf;

#[derive(Debug, Parser)]
struct Cli {
    #[arg(
        short,
        long,
        default_value_t = false,
        help = "use delete operator instead of refskip to fill in gaps"
    )]
    use_del: bool,
    ref_len: usize,
    bam_path: PathBuf,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();
    halve(args.ref_len, &args.bam_path, args.use_del)?;
    Ok(())
}
