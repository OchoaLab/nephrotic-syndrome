#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=ldpred-07-lassosum-fit
#SBATCH --output=ldpred-07-lassosum-fit.out
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=1
##SBATCH --mail-user=alejandro.ochoa@duke.edu
##SBATCH --mail-type=END,FAIL

module load R/4.1.1-rhel8

time Rscript ldpred-07-lassosum-fit.R 

module unload R/4.1.1-rhel8
