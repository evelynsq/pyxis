# Pyxis
A tool for detecting enriched motifs from ChIP-seq peaks using known motif position weight matrices (PWMs). Implemented as part of Advanced Bioinformatics Laboratory (CSE 185) final project at UC San Diego.

## Dependencies
`pyxis` requires the following libraries to be installed:

- pyfaidx
- numpy
- pandas
- scipy
- seqlogo

You can install these with `pip`:
```
pip install [--user] pyfaidx numpy pandas scipy seqlogo
```
(Note: Please specify `--user` option to install locally if you do not have root access.)

## Installation
Once the required libraries have been installed, `pyxis` can now be installed with the following:
```
git clone https://github.com/evelynsq/pyxis.git
cd pyxis
python setup.py install [--user]
```
If the installation was successful, you should see a helpful message when typing `pyxis --help`.

## Basic Usage
The basic usage of `pyxis` is:
```
pyxis peaks.bed ref.fa pwms.motifs [-o output.tsv] [other options]
```

To run `pyxis` on a small test example using the provided files:
```
pyxis example-files/peaks.bed example-files/ref.fa example-files/test.pwms
```
The output below should be produced:
```
motif_name  pval  log_pval  num_peaks_motif  pct_peak_motif  num_bg_motif  pct_bg_motif  enriched_status
(fillers - Under Construction)
```
To compare to the output of the analogous command `findMotifsGenome.pl` from [HOMER](http://homer.ucsd.edu/homer/):
```
findMotifsGenome.pl \
 example-files/peaks.bed \
 example-files/ref.fa \
 example-files/test.pwms
```

## pyxis Options
- `-h`, `--help`: Show help message and exit.
- `-o FILE`, --`output FILE`: Write output to specified file. Default: output is written to stdout.
- `-b FILE`, `--background FILE`: Use specified BED file for background peaks in motif-finding. Default: Background sequences are randomly generated from the reference genome.
- `-p PVAL`, --`pval PVAL`: Threshold p-value for enrichment significance. Default: 1e-5.
- `-s`, `--seqlogo`: Generate sequence logo for enriched motif. Default: True.
- `--version`: Print version and exit.

## File format

### Known Motif PWMs file
This tab-delimited input file should contain the PWMs of known motifs, but transposed with columns representing A, C, G and T nucleotides and the rows being the motif sequence. The header for each motif consists of the motif name, as well as an alphabet length `alength` (the size of 4 for DNA's bases) and `w` number of positions in the motif specifying the overall dimensions for its PWM.

```
motif_name    alength=4      w=18
PWM[i,A]      PWM[i,C]      PWM[i,G]    PWM[i,T]
PWM[i+1,A]    PWM[i+1,C]    PWM[i+1,G]  PWM[i+1,T]
```

### Peaks BED file
The input peaks file follows a standard [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html) with the following columns:
1. `chr`: Chromosome number.
2. `start`: Starting position in sequence.
3. `end`: Ending position in sequence. Note: Position numbers are 0-based with [start,end).
4. `peakID`: Identifier name of peak.
5. `strand`: + for forward, - for negative strand.

The file is tab-delimited and should look something like this:
```
chr    start    end    peakID    score    strand
9       30      56   MUSC.0001   0-1000    +/-
```

### Reference Genome FASTA file
The reference genome file is in [FASTA format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/). Rows starting with a carat (">") 
specify the chromosome number, followed by the chromosome's nucleotide sequence on the next line.

```
>chr6
ACTGTTTACTACTATCGACCCGTAAATC....
>chr7
TCGTTGATACCGTACATGCGATCGGATCGA....
```

### Output TSV file
The output file `pyxis_enrichments.tsv` contains several resulting statistics from our motif enrichment, including:
1. `motif_name`: identifier for motif.
2. `pval`: p-value for the motif being significantly enriched in the input peaks.
3. `log_pval`: log-10 value of the p-value.
4. `num_peak_motif`: number of peaks that matched this motif
5. `pct_peak_motif`: percentage of all peaks that matched motif
6. `num_bg_motif`: number of background sequences that matched this motif
7. `pct_bg_motif`: percentage of all peaks that matched motif
8. `enriched_status`: significance of enrichment in this motif, given a user-specified p-value threshold or otherwise a default of 1e-5. 
```
motif_name                pval    log_pval    num_peak_motif   pct_peak_motif    num_bg_motif    pct_bg_motif    enriched_status
ATF1_HUMAN.H11MO.0.B    1.23E-4   -3.91009          5               0.2                3              0.1            Non-sig             
```

### Motif Sequence Logo
If the option `-s` or `--seqlogo` is specified, a motif sequence logo for each motif will be printed. You can read more about `seqlogo` [here](https://pypi.org/project/seqlogo/).

## Contributors
This repository was created by Evelyn Quan with inspiration from [mypileup](https://github.com/gymreklab/cse185-demo-project).
Feel free to submit a pull request with any proposed corrections or suggestions.

## Testing
Under construction.
