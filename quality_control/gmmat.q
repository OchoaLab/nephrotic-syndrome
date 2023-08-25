#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=gmmat
#SBATCH --output=gmmat.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=92
##SBATCH --mail-user=tiffany.tu@duke.edu
##SBATCH --mail-type=END,FAIL

module load R/4.1.1-rhel8

time Rscript gmmat.R 

module unload R/4.1.1-rhel8
