## Benchmarking pyxis against HOMER findMotifsGenome.pl

The following notebook contains tests for benchmarking `pyxis` functions/commands for runtime and memory usage, both against example test files (several provided in the `example-files` directory) and also against a real [ChIP peaks ISL1 transcription factor dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5838045) (provided in the `analysis` directory). 

There are also comparisons drawn to [HOMER](http://homer.ucsd.edu/homer/ngs/peakMotifs.html)'s analogous command, `findMotifsGenome.pl`.

### Runtime

For testing `pyxis` functions' average runtimes, use this general format command as seen in `BenchmarkPyxis`:

```
%timeit readbed = pyxis.myutils.ReadBED(f_peaks_dir, test_ref)
```

We can also use this general command in Bash to run `%timeit` with `pyxis` and `HOMER` commands, though be aware that these runtimes may be large and take an hour or more depending on the size of the dataset and number of loop you are running on the command.

```
python -m timeit -n 10 -s 'import os' 'os.system("pyxis example-files/peaks.bed example-files/ref.fa example-files/test.pwms > /dev/null")'
```

I would suggest specifying `-r 1` for `findMotifsGenome.pl` with the ISL1 dataset.

### Memory

To test for memory consumption, you can clone the [memusg repo](https://github.com/jhclark/memusg.git), move the script into the outer `pyxis` directory, and then run your command with `pyxis` as you typically would, but through memusg. For example:

```
./memusg pyxis example-files/peaks.bed example-files/ref.fa example-files/test.pwms > /dev/null
```

For visualizations generated relating to these benchmarking results, please check `PyxisAnalysis.ipynb` in the `analysis` directory.

