#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=1
##SBATCH --mail-user=alejandro.ochoa@duke.edu
##SBATCH --mail-type=END,FAIL

module load R/4.1.1-rhel8

# $base $train $test: selects base, training, testing data types: NOTE: must be defined outside!!!
time Rscript ldpred-02-score.R $base $train $test

module unload R/4.1.1-rhel8
