#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=30

module load R/4.1.1-rhel8

# $base: selects base data: NOTE: must be defined outside!!!
time Rscript prs-new-07-ld-matched-snps.R $base

module unload R/4.1.1-rhel8
