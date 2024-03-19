#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=1

module load R/4.1.1-rhel8

# $base $train: selects base and training data types: NOTE: must be defined outside!!!
time Rscript ldpred-11-ct-fit.R $base $train

module unload R/4.1.1-rhel8
