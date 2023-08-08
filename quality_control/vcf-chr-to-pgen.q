#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=vcf-chr-to-pgen-%a
#SBATCH --output=vcf-chr-to-pgen-%a.out
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=alejandro.ochoa@duke.edu
#SBATCH --mail-type=END,FAIL

module load Plink/2.00a2LM

# parallelize by chr
chr=$SLURM_ARRAY_TASK_ID

# minor intermediate
name=chr$chr

# remove SNPs that don't "PASS", otherwise just reformats them for merging at a later step
time plink2 --vcf raw/$name.dose.vcf.gz --var-filter --make-pgen vzs --out $name

module purge
