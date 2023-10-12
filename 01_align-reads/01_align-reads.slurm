#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=normal
#SBATCH --time=48:00:00
#SBATCH --job-name sacerGBS
#SBATCH --output=sacerGBS_o%j
#SBATCH --error=sacerGBS_e%j
#SBATCH --mail-type=FAIL,END

# Loading conda environment (set up elsewhere, includes stacks, bwa, and samtools)
module load miniconda3/4.10.3-an4v
eval "$(conda shell.bash hook)"
conda activate stacks-env

# Variable to represent working directory
src=$SLURM_SUBMIT_DIR

# Path to genome for Todiramphus chloris collaris genome from Eliason et al. 2023 
bwa_db=$src/reference-genome/todChl.scaffolds.full_mask

### Alignment section ###

## Index fasta
bwa index -p $src/reference-genome/todChl.scaffolds.full_mask $src/reference-genome/todChl.scaffolds.full_mask.fa

## Align
## takes in a sample list of all samples, "sample_list.txt". This file is in the github folder for this step
## Raw read data is available on genbank (see sampling table S1)
while read indiv
do
        echo ${indiv}
        bwa mem -t 8 $bwa_db $src/raw_reads/${indiv}.fq.gz > $src/sam_files/${indiv}.sam &&
        samtools view -bS $src/sam_files/${indiv}.sam > $src/bam_files/${indiv}_unsort.bam &&
        echo "${indiv} converted to bam" &&
        samtools sort  $src/bam_files/${indiv}_unsort.bam -o $src/bam_files/${indiv}.bam &&
        echo "${indiv} sorted" &&
        rm $src/sam_files/${indiv}.sam &&
        rm $src/bam_files/${indiv}_unsort.bam
done < $src/sample_list.txt


