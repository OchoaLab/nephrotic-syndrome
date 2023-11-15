#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=ldpred
#SBATCH --output=ldpred.out
## --mem: 256G was submittable but didn't run (not enough resources), 512G didn't submit at all; 64G was not enough
## with one thread, 64G was not enough (OOM), 200 didn't run (not enough resources)
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=alejandro.ochoa@duke.edu
#SBATCH --mail-type=END,FAIL

module load R/4.1.1-rhel8

time Rscript ldpred.R 

module unload R/4.1.1-rhel8
