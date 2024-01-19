#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=ldpred-06-lassosum
#SBATCH --output=ldpred-06-lassosum.out
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=10
##SBATCH --mail-user=alejandro.ochoa@duke.edu
##SBATCH --mail-type=END,FAIL

module load R/4.1.1-rhel8

time Rscript ldpred-06-lassosum.R 

module unload R/4.1.1-rhel8
