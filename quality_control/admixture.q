#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=admixture
#SBATCH --output=admixture.out
#SBATCH --mem=4G
#SBATCH --ntasks-per-node=10

# name of dataset
name=array-clean
# number of ancestries
k=5

# run admixture! (must be on path!)
# use 10 threads (same as specified to sbatch above!)
time admixture -j10 $name.bed $k > array-clean.$k.log

# compress big outputs
gzip $name.$k.{P,Q}

# results go in subdirectory
mkdir admixture/
mv $name.$k.{P.gz,Q.gz,log} admixture/
