use std::error::Error;
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::str::from_utf8;
use anyhow::Result;
use noodles::fasta;
use flate2::read::MultiGzDecoder;
use poasta::api::PoastaConsensus;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    let seqs: Vec<Vec<u8>> = load_fasta(&args[1])?;
    let consensus = seqs.iter()
        .map(|seq| (seq.as_ref(), None))
        .build_consensus_with_defaults()?;
    for seq in seqs {
        println!("{}", from_utf8(&seq)?);
    }
    println!();
    println!("{}", from_utf8(&consensus)?);
    Ok(())
}

fn load_fasta(sequences_fname: &str) -> Result<Vec<Vec<u8>>, Box<dyn Error>> {
    let sequences_fname = Path::new(sequences_fname);

    // Let's read the sequences from the given FASTA
    let is_gzipped = sequences_fname
        .file_name()
        .map(|v| v.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false);

    // Check if we have a gzipped file
    let reader_inner: Box<dyn std::io::BufRead> = if is_gzipped {
        Box::new(
            File::open(sequences_fname)
                .map(MultiGzDecoder::new)
                .map(BufReader::new)?,
        )
    } else {
        Box::new(File::open(sequences_fname).map(BufReader::new)?)
    };
    let mut reader = fasta::io::Reader::new(reader_inner);

    let mut seqs: Vec<_> = Vec::new();
    for result in reader.records() {
        let record = result?;
        seqs.push(record.sequence().as_ref().to_vec());
    }
    Ok(seqs)
}
