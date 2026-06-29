# Sonar

**Experimental** RNA sensor design tool.

## Usage

The workflow is generally
1. Generate candidate sensor library (`sonar gen`)
2. RNA secondary structure and accessibility prediction (`RNAfold`, `RNAplfold`)
    - Generally, RNAfold can be used to compute the whole-sensor self folding / ensemble diversity.
    - RNAplfold is used to determine the free-sensor accessability around the edited stop codons.
3. Interaction metrics (`IntaRNA`)
    - sensor-target interaction energy and whether the interaction covers the edited stop regions.
4. Non-specificity checks / transcriptome mapping (`sonar spec`)


The `gen` command will generate a candidate sensor library from a given FASTA file containing
the target sequence(s). For example,

```bash
sonar gen target.fa 
```

This will produce a columnar output (default `tsv`) containing the candidate sensor information
like the target ID, start position, end position, GC content, etc. The sensor sequence information
is output in FASTA format. The output table and FASTA file paths can be set with the `--table` and
`--fasta` options. For a full list of options, see the [Commands](#commands)
section or use `sonar gen --help`.

## Installation

```bash
cargo install sonar-rna

# or with binstall if you have it available
cargo binstall sonar-rna
```

### GitHub

The CLI is available in the [GitHub releases](https://github.com/szablowskilab/sonar/releases). Download the binary
for your platform and add it to your PATH.

### Source

To build from source, make sure you have `cargo` installed. The provided Nix flake
can be used to create a new environment with all the necessary dependencies. The
binary can then be built using `cargo install`. For example:

```bash
nix develop
cargo install --path crates/sonar-cli
```

## Commands

### gen

```
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

### spec

```txt
Map RNA sensor candidates to a reference transcriptome for specificity filtering

Usage: sonar spec [OPTIONS] <SENSOR_PATH> <REFERENCE_PATH>

Arguments:
  <SENSOR_PATH>     Input sensor candidates FASTA file
  <REFERENCE_PATH>  Reference transcriptome FASTA or index file

Options:
  -m, --best-n <BEST_N>          Number of high scoring secondary alignments to consider [default: 100]
  -n, --allow-iupac              Allow IUPAC characters in sensor sequences
  -i, --save-index <SAVE_INDEX>  Optionally save reference index to the provided path
  -o, --output <OUTPUT>          Output mapping file
  -h, --help                     Print help
```

## Library

Sonar is also available as a library that you can use in your Rust projects. To add it to your project,
use `cargo add sonar-rna` or add the following to your `Cargo.toml`:

```toml
[dependencies]
sonar = { package = "sonar-rna", version = "0.1.0-alpha.1" }
```

Please refer to the [API documentation](https://docs.rs/sonar-rna) for more information.

## TODO

In the future, I would like to integrate folding, interaction prediction, and non-specificity checks
directly to the library. However, this is low priority since existing tools already exist. A middle
ground might be to add a `filter` command that allows for re-ranking and filtering based on RNAfold,
RNAplfold, IntaRNA, and mm2 output for example.
