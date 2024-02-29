#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=30
##SBATCH --mail-user=tiffany.tu@duke.edu
##SBATCH --mail-type=END,FAIL

module load R/4.1.1-rhel8

# $base $train: selects base and training data types: NOTE: must be defined outside!!!
time Rscript prs-new-07-ld-matched-snps.R $base $train

module unload R/4.1.1-rhel8
