## File Information

This directory contains example test files to be run with the following command for testing `pyxis` usage:
```
pyxis example-files/peaks.bed example-files/ref.fa example-files/test.pwms [-b example-files/background.bed] [-s]
```

Something similar to the output below should be produced:
```
motif_name      	pval   log_pval  num_peak_motif  pct_peak_motif  num_bg_motif    pct_bg_motif    enriched_status
BACH2_HUMAN.H11MO.0.A   1.0     0.0     	3       	0.3     	2		0.2		Non-sig
ALX3_HUMAN.H11MO.0.D    1.0     0.0     	2       	0.2     	3       	0.3     	Non-sig
ELK1_HUMAN.H11MO.0.B    1.0     0.0     	0       	0.0     	1      		0.1     	Non-sig
KAISO_HUMAN.H11MO.1.A   1.0     0.0     	2       	0.2     	1       	0.1    		Non-sig
MLX_HUMAN.H11MO.0.D     1.0     0.0     	1       	0.1     	1       	0.1   		Non-sig
```

To compare to the output of the analogous command `findMotifsGenome.pl` from [HOMER](http://homer.ucsd.edu/homer/):
```
findMotifsGenome.pl \
 example-files/peaks.bed \
 example-files/ref.fa \
 homer-results/example-data
```

## File formats

### Known Motif PWMS file
This tab-delimited input file should contain the PWMs of known motifs, but transposed with columns representing A, C, G and T nucleotides and the rows being the motif sequence. The header for each motif consists of the motif name, as well as an alphabet length `alength` (the size of 4 for DNA's bases) and `w` number of positions in the motif specifying the overall dimensions for its PWM.

```
motif_name	alength=4	w=18
PWM[i,A]	PWM[i,C]	PWM[i,G]	PWM[i,T]
PWM[i+1,A]	PWM[i+1,C]	PWM[i+1,G]	PWM[i+1,T]
```

### Peaks BED file
The input peaks file follows a standard [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html) with the following columns:
1. `chrom`: Chromosome number.
2. `start`: Starting position in sequence.
3. `end`: Ending position in sequence. Note: Position numbers are 0-based with [start,end).
4. `peakID`: Identifier name of peak.
5. `.`: Column 5 not used.
6. `strand`: + for forward, - for negative strand.
7. `score`: Score value for peak (a number between 0-1000).

The file is tab-delimited and should look something like this:
```
chrom	start	end	peakID		.	strand	score
chr9	30	56	MUSC.0001	.	+	900
```

### Reference Genome FASTA file
The reference genome file is in [FASTA format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/). Rows starting with a carat (">") 
specify the chromosome number, followed by the chromosome's nucleotide sequence on the next line.

```
>chr6
ACTGTTTACTACTATCGACCCGTAAATCAT....
>chr7
TCGTTGATACCGTACATGCGATCGGATCGA....
```

### Output TSV file
The output file `pyxis_enrichments.tsv` contains several resulting statistics from our motif enrichment, including:
1. `motif_name`: identifier for motif.
2. `pval`: p-value for the motif being significantly enriched in the input peaks.
3. `log_pval`: log-10 of the p-value.
4. `num_peak_motif`: number of peaks that matched this motif.
5. `pct_peak_motif`: fraction of all peaks that matched motif.
6. `num_bg_motif`: number of background sequences that matched this motif.
7. `pct_bg_motif`: fraction of all peaks that matched motif.
8. `enriched_status`: significance of enrichment in this motif, given a user-specified p-value threshold or otherwise a default of 1e-5. 
```
motif_name		pval	   log_pval	num_peak_motif	pct_peak_motif	num_bg_motif	pct_bg_motif	enriched_status
ATF1_HUMAN.H11MO.0.B	1.23E-4	   -3.91009	5		0.2		3		0.1		Non-sig             
```

### Motif Sequence Logo
If the option `-s` or `--seqlogo` is specified, a motif sequence logo for each motif will be printed. You can read more about `seqlogo` [here](https://pypi.org/project/seqlogo/).

### Background BED file
The user-specified background file also follows the same format as the above [input peaks BED file](#peaks-bed-file).

**Note**: The background has the same number of sequences specified as the number of peaks in `peaks.bed`, and the background sequences are also the same lengths as the peak sequences. This is the standard assumption for any commands run with `pyxis` including a user-specified background peaks file.

