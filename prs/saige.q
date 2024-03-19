#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=prs_ssns_ctrl
#SBATCH --output=prs_ssns_ctrl_saige.out
#SBATCH --mem=180G
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

export LD_LIBRARY_PATH=/opt/apps/rhel8/lapack/lib64:$LD_LIBRARY_PATH
module load R/4.1.1-rhel8

# saige with binary traits
time Rscript saige_step1.R -d ssns_ctrl
time Rscript saige_step2.R -d ssns_ctrl
module unload R/4.1.1-rhel8
