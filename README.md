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
wget https://github.com/rmFlynn/16s_to_bins_project/blob/main/environment.yaml
conda env create -f environment.yaml -n join_asvbins
conda activate join_asvbins
```

## Use Example
```
/bin/time --verbose  --output=/home/projects/16s_to_bins_project/results/salmonella_small_time_post_filter.txt \
        join_asvbins \
        -b data/salmonella_small/combined_bins.fna \
        -a data/salmonella/ASVs_r1-r4_dna-sequences.fasta \
        -o ./results/salmonella_small_t1 \
        -t 2 #
```
The most important comand line options are:

`--no_clean`

But there are many more



## To do
 * Split mmseqs into multi steps
 * More general tests
 * Add a test to run the full snake pipline.
