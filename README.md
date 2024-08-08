# chop_reads

This is a command line tool for chopping alignments in a BAM into smaller pieces, while preserving the CIGAR logic.

## Documentation

See the output from `--help` for supported usage.

```
Usage: chop-reads [OPTIONS] --input <INPUT> --output <OUTPUT> --chunk-size <CHUNK_SIZE>

Options:
  -i, --input <INPUT>              Input file to chop records from
  -r, --reference <REFERENCE>      Path to reference file to use with crams
  -o, --output <OUTPUT>            Path to write output to
  -s, --chunk-size <CHUNK_SIZE>    Length of chunks to split records into
      --min-length <MIN_LENGTH>    Min record length to include in chopped outputs when handling final chunk [default: 0]
  -g, --read-group <READ_GROUP>    Read group value to use for new split records
  -n, --sample-name <SAMPLE_NAME>  Sample name to use for new read group
  -h, --help                       Print help
```