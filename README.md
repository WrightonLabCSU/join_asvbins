# Comparing 16s(ASV) to Bins

We want to connect our 16S rRNA sequences to our bins to the best ability of
the data. There are many ways to go about this but we will focus on the uses
of barrnap and BLAST to extract 16s information from the bins themselves.
This method is the most accessible and well documented as of yet, however
it also has known drawback. To understand this method and its us here it is
best to look at the "Comparing ASV to Bins Protocol, By: Katherine Kokkinias
and Kai Leleiwi" for more detail.

## Test Data Sets

The fallowing data sets will be used for testing the program and will reside
in the data folder organized by name

 * M3C4D3_v2_megahit from /home/projects/Wetlands/All_genomes/miller_lab_JGI_CSP_bins
 This files was chosen because I wanted to use wetlands data as a test case,
 and it was good to have some data that was assembled with MEGAHIT. using
 this file is mostly random and will be probably be replaced.

 * traing_c8a_idba from /home/projects/Training/Rory_working/all_hw/bining_data/scaffold.fa.metabat-bins/
 This is the only set of bins that I have myself created, and I will probably recreate them.

 * Salmonella
 /home/projects-wrighton/NIH_Salmonella/Salmonella/Metagenomes/MAG_database/16S_From_All_Bins

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
