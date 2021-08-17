#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 1
#SBATCH --time 6:00:00
#SBATCH --job-name vC
#SBATCH --chdir /scratch/ldelisle/AmandioEtAl2021/cHi-C/

gitHubDirectory=$1

path="$PWD/"

module purge

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda activate amandio2021

mkdir -p vC

cd vC

# The Digester_File_chr2.txt was obtained by:
# First on galaxy:
# hicup_digester --re1 '^GATC' --genome 'mm10' mm10_UCSC.fa 
# Then:
# grep "chr2" Digester_File.txt > Digester_File_chr2.txt

while read l; do
  sample=$(echo $l | awk '{print $1}')
  coo=$(echo $l | awk '{print $2}')
  output=$(echo $l | awk '{print $3}')
  if [ -e $output ]; then
    echo "$output already exists"
  else
    if [ ! -e ${sample}.validPair ]; then
      if [ -e ${path}/toGEO/${sample}_validPairs_filtered.csort.gz ]; then
        zcat ${path}/toGEO/${sample}_validPairs_filtered.csort.gz > ${sample}.validPair
      else
        cat ${path}/toGEO/${sample/merge/rep}*_validPairs_filtered.csort.gz | gunzip -c > ${sample}.validPair
      fi
    fi
    python  ${gitHubDirectory}/scripts/fromFragFileAndValidPairsToVirtualCaptureC.py \
      --validPair ${sample}.validPair \
      --colChr1 3 --colChr2 7 --colFrag1 5 --colFrag2 9 --lineToSkipInValidPair 0 \
      --fragmentFile ${gitHubDirectory}/cHi-C/Digester_File_chr2.txt \
      --colForChr 1 --colForStart 2 --colForEnd 3 --colForID 4 \
      --lineToSkipInFragmentFile 0 --viewpointCoo ${coo} --output ${output}
  fi
done < ${gitHubDirectory}/cHi-C/vC_table.txt
