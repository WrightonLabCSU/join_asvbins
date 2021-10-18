# Comparing 16s(ASV) to Bins

We want to connect our 16S rRNA sequences to our bins to the best ability of
the data. There are many ways to go about this but we will focus on the uses
of barrnap and BLAST to extract 16s information from the bins themselves.
This method is the most accessible and well documented as of yet, however
it also has known drawback. To understand this method and its us here it is
best to look at the "Comparing ASV to Bins Protocol, By: Katherine Kokkinias
and Kai Leleiwi" for more detail.

## Install


## Use Example
 '''
join_asvbins \
        -b data/salmonella_small/combined_bins.fna \
        -a data/salmonella/ASVs_r1-r4_dna-sequences.fasta \
        -w results/salmonella_small_mmseqs4 \
        -o results/salmonella_small_mmseqs4 \
        -t 32 \
 '''

## To do
 * TODO LOTS.
      conda install -c conda-forge gcc_linux-64
