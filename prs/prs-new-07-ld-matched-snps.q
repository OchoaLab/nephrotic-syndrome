#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=30
##SBATCH --mail-user=tiffany.tu@duke.edu
##SBATCH --mail-type=END,FAIL

module load R/4.1.1-rhel8

# $type: selects base data subtype (for old cases only): NOTE: must be defined outside!!!
time Rscript prs-new-07-ld-matched-snps.R $type

module unload R/4.1.1-rhel8
