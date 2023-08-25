#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=grm-ns-ctrl
#SBATCH --output=grm-ns-ctrl.out
##SBATCH --job-name=grm-ssns-ctrl
##SBATCH --output=grm-ssns-ctrl.out
##SBATCH --job-name=grm-srns-ctrl
##SBATCH --output=grm-srns-ctrl.out
##SBATCH --job-name=grm-ssns-srns
##SBATCH --output=grm-ssns-srns.out
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=6
##SBATCH --mail-user=tiffany.tu@duke.edu
##SBATCH --mail-type=END,FAIL

# select subtype
dir=. # == ns_ctrl
# dir=ssns_ctrl
# dir=srns_ctrl
# dir=ssns_srns

# run!
time gcta64 --bfile $dir/mac20 --make-grm --out $dir/mac20
time gcta64 --grm $dir/mac20 --pca 10 --out $dir/mac20
