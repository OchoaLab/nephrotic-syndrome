#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --mem=180G
#SBATCH --ntasks-per-node=24

export LD_LIBRARY_PATH=/opt/apps/rhel8/lapack/lib64:$LD_LIBRARY_PATH
module load R/4.1.1-rhel8

# saige with binary traits
echo "step1"
time Rscript saige_step1_alex.R --bfile $name
echo "step2"
time Rscript saige_step2_alex.R --bfile $name
module unload R/4.1.1-rhel8
