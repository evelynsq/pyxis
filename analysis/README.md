## Pyxis Analysis

The following notebook performs analyses of `pyxis` results to HOMER's `findMotifsGenome.pl`, generating visualizations comparing results, runtimes, and sequence logos.

Note: As mentioned in `BenchmarkPyxis.ipynb`, the reference genome `hg19.fa` required for running commands with the ISL1 dataset is not included in the `pyxis` repository due to its large size. 

If you are interested in running these commands with the ISL1 dataset (as well as several of the commands in `BenchmarkPyxis.ipynb`), please install `hg19.fa.gz` [here](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and move it into the `benchmark` directory of `pyxis`. Be sure to also decompress the file before running `pyxis` on it.

You will also not see the images for the sequence logos of `pyxis` for the ISL1 Dataset on the notebook unless you run the command yourself to generate these logos.
