# Comparing 16s(ASV) to Bins

This tool provides a method for connecting 16S rRNA sequences to a set of bins.
The code is in essence a wrapper around a Snakemake pipeline that uses barnap and MMseqs2, with the option of substituting BLAST for MMseqs2.
This is not a perfect system, and is also a work in progress.

## Install

The tool can run on as few as 1 core but it will utilize as many cores as specified by the -t argument.  MMseqs2 is much faster than BLAST but will require more memory, often in the range of 40 - 60 gigabytes depending on the size of the target data set. Using blast will decrease the memory requirements but will also require a long run time.

## Install

The installation should be no more complex than:

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
        -t 20 # Threads
```

The most important command line options are:

```
  -b BINS, --bins BINS  The bin that you would like to match asvs to. This can be an fna file that has all the bins
                        combined or a directory of bins in separate fa files, but you must run the rename script
                        before you use this tool

  -a ASV_SEQS, --asv_seqs ASV_SEQS
                        The asvs you would like to attach to your bins.

  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        The folder where you would like the temporary files and the final output to be stored

  -t THREADS, --threads THREADS
                        The number of threads that will be used by the program and subprocess.

  --no_clean            Specifies that the directory should NOT be cleaned of results of pass runs. If your run is
                        interrupted this will allow you to to pick up. where you left off. Use at your own risk.

```
But there are many more, use `join_asvbins -h` to see all of them. A options section will also be added to the wiki soon

## How it works

![example_dag](./images/salmonella_dag.pdf)


## To do
 * Split mmseqs into multi steps
 * More general tests
 * Add a test to run the full snake pipeline.
 * Add a better file path does not exist message, snakemakes does not handle this with this setup.
 * Check that verbos works like quiet



# Current Testers in the Wrighton Lab
Please use the `--generic_16s` argument to run the program, pulling data from git LFS is currently problematic in python. This will be fixed in the future but for now you can contact me and I will provide a clustered SILVA dataset for reproducibility.



# Current Testers in the Wrighton Lab
Please use the `--generic_16s` argument to run the program, pulling data from git LFS is currently problematic in python. This will be fixed in the future but for now you can contact me and I will provide a clustered SILVA dataset for reproducibility.
