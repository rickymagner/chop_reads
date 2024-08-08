use std::path::PathBuf;
use std::sync::Mutex;
use rayon::prelude::*;
use std::sync::mpsc::sync_channel;
use rust_htslib::bam as hts_bam;
use rust_htslib::bam::{Format, Read};
use std::time::Instant;
use clap::Parser;
use rayon::iter::{IntoParallelIterator, ParallelBridge, ParallelIterator};
use rust_htslib::bam::header::HeaderRecord;
use chop_reads::alignment_chopper::AlignmentChopper;


#[derive(Parser, Debug)]
struct Cli {
    /// Input file to chop records from
    #[arg(short, long)]
    input: PathBuf,

    /// Path to reference file to use with crams
    #[arg(short, long)]
    reference: Option<PathBuf>,

    /// Path to write output to
    #[arg(short, long)]
    output: PathBuf,

    /// Length of chunks to split records into
    #[arg(short='s', long)]
    chunk_size: u32,

    /// Min record length to include in chopped outputs when handling final chunk
    #[arg(long, default_value_t=0)]
    min_length: u32,

    /// Toggle whether to skip softclipped bases at edges of record
    #[arg(long, default_value_t=true)]
    skip_clipped_bases: bool,

    /// Read group value to use for new split records
    #[arg(short='g', long)]
    read_group: Option<String>,

    /// Sample name to use for new read group
    #[arg(short='n', long, requires("read_group"))]
    sample_name: Option<String>,

    /// Number of threads to use
    // #[arg(short, long, default_value_t=1)]
    // threads: u32,
}

fn main() {
    let now = Instant::now();

    let args = Cli::parse();

    let mut hts_reader = hts_bam::Reader::from_path(args.input).unwrap();
    let mut header = hts_bam::header::Header::from_template(hts_reader.header());

    if let Some(rg) = &args.read_group {
        let mut header_record = HeaderRecord::new(b"RG");
        header_record.push_tag(b"ID", rg);
        if let Some(sn) = args.sample_name {
            header_record.push_tag(b"SM", &sn);
        }

        header.push_record(&header_record);
    }

    let mut hts_writer = hts_bam::Writer::from_path(args.output, &header, Format::Bam).unwrap();

    let mut alignment_chopper = AlignmentChopper::new(args.chunk_size, args.min_length, args.skip_clipped_bases, args.read_group.clone());

    let mut record = hts_bam::Record::new();
    while let Some(r) = hts_reader.read(&mut record) {
        r.expect("Failed to parse record");
        for cr in alignment_chopper.chop_read(&record) {
            hts_writer.write(&cr).expect("Cannot write record.");
        }
    }

    println!("Runtime: {}s", now.elapsed().as_secs());
}
