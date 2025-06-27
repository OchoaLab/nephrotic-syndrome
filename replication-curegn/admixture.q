#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=admixture_maf10
#SBATCH --output=admixture_maf10.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=10

# name of dataset
name=curegn_tgp_merge_maf10
# number of ancestries
k=5

# run admixture! (must be on path!)
# use 10 threads (same as specified to sbatch above!)
time admixture -j10 $name.bed $k > curegn_tgp_merge_maf10_.$k.log

# # compress big outputs
# gzip $name.$k.{P,Q}
# 
# # results go in subdirectory
# mkdir results/
# cd results/
# mv $name.$k.{P.gz,Q.gz,log} results/
