# Pyxis
A tool for detecting enriched motifs in ChIP peaks from known motif position weight matrices (PWMs). Implemented as part of Advanced Bioinformatics Laboratory (CSE 185) final project at UC San Diego.

## Dependencies
- pyfaidx
- numpy
- pandas
- scipy
- seqlogo

`pip install [--user] pyfaidx numpy pandas scipy seqlogo`

(Note: Please specify `--user` if you do not have root access.)

## Installation

```
git clone https://github.com/evelynsq/pyxis.git
cd pyxis
python setup.py install [--user]
```

## Basic Usage

The basic usage of `pyxis` is:
`pyxis peaks.bed ref.fa pwms.motifs [-o output.tsv] [other options]`

To run `pyxis` on a small test example:
`pyxis ...`

## Pyxis options
- -b FILE, --background FILE: Use specified BED file for background peaks in motif-finding. Default: Background sequences are randomly generated from the reference genome.
- -p PVAL, --pval PVAL: Threshold p-value for enrichment significance. Default: 1e-5.
- -o FILE, --output FILE: Write output to specified file. Default: output is written to stdout.

## File format

### Known PWMs MOTIFS file

### Peaks BED file

### Reference Genome FASTA file

### Output TSV file


## Contributors
This repository was created by Evelyn Quan with inspiration from [mypileup](https://github.com/gymreklab/cse185-demo-project).
Feel free to submit a pull request with any proposed corrections or suggestions.
