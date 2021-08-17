#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 1
#SBATCH --time 1:00:00
#SBATCH --job-name getDep
#SBATCH --chdir /scratch/ldelisle/AmandioEtAl2021/RNAseq

pathWithDependencies=$1

module purge
cd $pathWithDependencies
if [ ! -e rnaseq_rscripts ]; then
  git clone https://github.com/lldelisle/rnaseq_rscripts.git
fi
cd rnaseq_rscripts
git checkout "8683b0ce6e70789a5cf0300d1ab3fe71d583b81e"
cd ..
if [ ! -e toolBoxForMutantAndWTGenomes ]; then
  git clone https://github.com/lldelisle/toolBoxForMutantAndWTGenomes.git
fi
cd toolBoxForMutantAndWTGenomes
git checkout "62ab52414d81dcf796498509b6ca3445553317e8"
