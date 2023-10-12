# Sacer-RADseq
Code used in McCullough et al. in review "Rapidly Radiating Clade of ‘Great Speciators’ Harbors Rampant Incomplete Lineage Sorting and Low Levels of Gene Flow (Aves: Todiramphus kingfishers)" 

This repository has code for data processing and analyses. 
See Dryad repository (link to include here when generated) for supplementary information and datasets. 

We created a dataset will all samples (n=63) and seven different subsets for PCAs/ADMIXTURE etc. 
Some analyses were run on the UNM Center for Advanced Research computing, which requires a slurm submission script to the cluster and used a conda environment that had stacks, bwa, and samtools installed. Other code was run in R. Color editing for PCAs/admixture plots was done by hand in Adobe Illustrator.  

1. 01_align: folder with a slurm that aligns demultiplexed reads to Mariana Kingfisher reference genome and generation of bam files using SAMtools
3. 02_stacks: Stacks to build catalog of RAD loci and creation of VCF files for all samples and 7 subsets. Includes popmaps
4. 03_SNPfiltr: R code for SNPFiltr for the full dataset and subsets
5. 
