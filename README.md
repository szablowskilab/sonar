# Sonar

**Experimental** RNA sensor design tool.

## Usage
The `gen` command will generate a candidate sensor library from a given FASTA file containing
the target sequence(s).

```bash
sonar gen target.fa 
```

This will produce a columnar output (default `tsv`) containing the candidate sensor information
like the target ID, start position, end position, GC content, etc. The sensor sequence information
is output in FASTA format. The output table and FASTA file paths can be set with the `--table` and
`--fasta` options. For a full list of options, see the [Commands](#commands)
section or use `sonar gen --help`.

## Installation

### Source

To build from source, make sure you have `cargo` installed. The provided Nix flake
can be used to create a new environment with all the necessary dependencies. The
binary can then be built using `cargo build --release`. For example:

```bash
nix develop
cargo build --release
```

## Commands

### Generating sensor candidates

```bash
Usage: sonar gen [OPTIONS] <TARGET>

Arguments:
  <TARGET>  Input target FASTA file

Options:
      --region <REGION>
          1-based inclusive target subregion, e.g. 10:350
      --ses-length <SES_LENGTH>
          sesRNA length range in nucleotides as min:max [default: 200:300]
      --stop-count <STOP_COUNT>
          Number of designed TAG stop codons to include [default: 1]
      --stop-window <STOP_WINDOW>
          allowed design stop position window as min:max from the sesRNA 5' end [default: 80:220]
      --min-stop-distance <MIN_STOP_DISTANCE>
          Minimum distance between designed TAG and any converted stop codons [default: 10]
      --frame <FRAME>
          Translation frame for stop and ATG checks [default: 0]
      --max-candidates <MAX_CANDIDATES>
          Max number of candidates to output [default: 100]
      --no-spread
          Do not spread out candidate designs across the target sequence
      --allow-iupac
          Allow IUPAC nucleotide codes in the target sequence
      --table <TABLE>
          Write table to this path instead of stdout
      --fasta <FASTA>
          Write FASTA output to this path instead of stdout
  -h, --help
          Print help
```
