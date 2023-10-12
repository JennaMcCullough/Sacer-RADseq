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

# Sacer63: the full dataset 
gstacks -I $src/bam_files/ -M $src/popmaps/sacer63.txt -O $src/populations_out/sacer63/ -t 8
populations -P $src/populations_out/sacer63 -M $src/popmaps/sacer63.txt -O $src/populations_out/pop_sacer63_50per --vcf -R .5 --fstats --smooth -t 8
mv $src/populations_out/pop_sacer63_50per/populations.snps.vcf $src/populations_out/outputs/sacer63_50per.vcf

# the rest of the subsets follow the same idea (pulling from the same bam_files) but 
# with different population maps. I've included the different population maps in this 
# github folder. name the output files appropriately. 
# the output files (i.e., sacer63_50per.vcf) are what goes into R for SNPfiltr  
# 
# subsets are as follows: 
# sacer63 = all samples
# sacer45 = only sacer 
# PolySacer26 = all polynesian clade sacer (w/ pealei singleton)
# PolySacer25 = same above without pealei 
# marinus13 = only marinus
# MelanSacer14 = melanesian sacer (w/ ornatus singleton)
# MelanSacer13 = same as above without ornatus
# Micro10 = Micronesian and Eastern Pacific non-sacer taxa 