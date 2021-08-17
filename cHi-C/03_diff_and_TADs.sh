#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 1
#SBATCH --time 6:00:00
#SBATCH --job-name findTADs
#SBATCH --chdir /scratch/ldelisle/AmandioEtAl2021/cHi-C/

gitHubDirectory=$1

path="$PWD/"

nbOfThreads=1

module purge

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda activate amandio2021

mkdir -p TADs

cd TADs

# hicFindTADs
for cool in ${path}/toGEO/E9.5*merge*; do
  hicFindTADs -m ${cool} --minBoundaryDistance 200000 \
    --correctForMultipleTesting fdr --outPrefix $(basename $cool .cool).240kb \
    --minDepth 240000 --maxDepth 480000 --step 480000
done

cd ..

mkdir -p diff

cd diff

# Compare to corresponding WT
for f in $(ls ${path}/toGEO/*.cool | grep -v rep | grep -v WT); do
  wt=$(echo $f | awk -F "_" -v OFS="_" '{$3="WT"; print}')
  output=WT_minus_${f}
  if [ ! -e $output ]; then
    hicCompareMatrices --matrices ${wt} ${f} \
      --operation diff --outFileName ${output}
  fi
done

cd ..

mkdir -p correct
cd correct
hicPlotDistVsCounts --matrices ${path}/toGEO/E9.5*merge* --plotFile ${gitHubDirectory}/cHi-C/distvscount.png --maxdepth 80000 --plotsize 8 4
for f in ${path}/toGEO/E9.5*merge*; do
  hicTransform --matrix ${f} --outFileName $(basename ${f} .cool)_hicTransform.cool --method obs_exp_non_zero --perChromosome
done

# Compare to corresponding WT
for f in $(ls *.cool | grep -v WT); do
  wt=$(echo $f | awk -F "_" -v OFS="_" '{$3="WT"; print}')
  output=WT_minus_${f}
  if [ ! -e $output ]; then
    hicCompareMatrices --matrices ${wt} ${f} \
      --operation diff --outFileName ${output}
  fi
done
