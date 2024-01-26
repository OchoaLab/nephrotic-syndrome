#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=6
##SBATCH --mail-user=tiffany.tu@duke.edu
##SBATCH --mail-type=END,FAIL

# $dir: selects dataset/subtype: NOTE: must be defined outside!!!

# run!
time gcta64 --bfile $dir/mac20 --make-grm --out $dir/mac20
time gcta64 --grm $dir/mac20 --pca 10 --out $dir/mac20
