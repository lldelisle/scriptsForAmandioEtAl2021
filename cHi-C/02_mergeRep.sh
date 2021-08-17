#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 1
#SBATCH --time 6:00:00
#SBATCH --job-name mergeRep
#SBATCH --chdir /scratch/ldelisle/AmandioEtAl2021/cHi-C/toGEO

gitHubDirectory=$1

path="$PWD/"

nbOfThreads=1

module purge

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda activate amandio2021

# Merge replicates
for f in *rep1*cool; do
  if [ ! -e ${f/rep1/merge} ]; then
    tmpdir=$(mktemp -d)
    hicConvertFormat --matrices $f --inputFormat cool --load_raw_values --outputFormat cool -o ${tmpdir}/$f
    hicConvertFormat --matrices ${f/rep1/rep2} --inputFormat cool --load_raw_values --outputFormat cool -o ${tmpdir}/${f/rep1/rep2}
    hicSumMatrices --matrices ${tmpdir}/$f ${tmpdir}/${f/rep1/rep2} -o ${f/rep1/merge}
    cooler balance --mad-max 5 --min-nnz 10 --min-count 0 --ignore-diags 2 --tol 1e-05 --max-iters 200 --cis-only -f ${f/rep1/merge}
  fi
done
