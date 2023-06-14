#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=srns_ssns_or_batch-%a
#SBATCH --output=srns_ssns_or_batch-%a.out
#SBATCH --mem=2G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=alejandro.ochoa@duke.edu
#SBATCH --mail-type=END,FAIL

module load R/4.0.0

batch=$SLURM_ARRAY_TASK_ID
time Rscript gmmat_OR_srns_ssns_batch.R $batch

module unload R/4.0.0
