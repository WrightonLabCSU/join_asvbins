# Comparing 16s(ASV) to Bins

We want to connect our 16S rRNA sequences to our bins to the best ability of
the data. There are many ways to go about this but we will focus on the uses
of barrnap and BLAST to extract 16s information from the bins themselves.
This method is the most accessible and well documented as of yet, however
it also has known drawback. To understand this method and its us here it is
best to look at the "Comparing ASV to Bins Protocol, By: Katherine Kokkinias
and Kai Leleiwi" for more detail.

## Install

The instalation should be no more complex than:

```
wget https://raw.githubusercontent.com/rmFlynn/16s_to_bins_project/main/environment.yaml
conda env create -f environment.yaml -n join_asvbins
conda activate join_asvbins
```

## Use Example
```
join_asvbins \
        -b path/to/bins/folder/or/file.fa \
        -a /path/to/asv/file.fa \
        -o /path/to/output \
        -t 32 # Threads
```

The most important comand line options are:

```
  -b BINS, --bins BINS  The bin that you would like to match asvs to. This can be an fna file that has all the bins
                        combided or a directory of bins in seperate fa files, but you must run the rename script
                        before you use this tool

  -a ASV_SEQS, --asv_seqs ASV_SEQS
                        The asvs you would like to atach to your bins.

  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        The folder where you would like the temporary files and the final output to be stored

  -t THREADS, --threads THREADS
                        The number of threads that will be used by the program and subprocess.

  --no_clean            Specifies that the directory should NOT be cleaned of results of pass runs. If your run is
                        interrupted this will allow you to to pickup. where you left off. Use at your own risk.

```
But there are many more, use `-h` to see all of them.




## To do
 * Split mmseqs into multi steps
 * More general tests
 * Add a test to run the full snake pipline.
 * Add a beter file path dose not exist mesage, snakemakes dose not handle this with this setup.
 * Check that verbos works like quiet
