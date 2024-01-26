#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=gmmat-base
#SBATCH --output=gmmat-base.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=20
##SBATCH --mail-user=tiffany.tu@duke.edu
##SBATCH --mail-type=END,FAIL

module load R/4.1.1-rhel8

time Rscript prs-new-03-gmmat.R 

module unload R/4.1.1-rhel8
