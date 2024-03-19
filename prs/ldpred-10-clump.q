#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=10

module load R/4.1.1-rhel8

# $base: selects base data: NOTE: must be defined outside!!!
time Rscript ldpred-10-clump.R $base

module unload R/4.1.1-rhel8
