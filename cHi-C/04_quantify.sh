#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 1
#SBATCH --time 20:00
#SBATCH --job-name quantify
#SBATCH --chdir /scratch/ldelisle/AmandioEtAl2021/cHi-C/

gitHubDirectory=$1

path="$PWD/"

module purge

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda activate amandio2021

mkdir -p quantif

cd quantif
my_regions=${gitHubDirectory}/cHi-C/trunk_S7.bedpe
c_cool=${path}/correct/E9.5_trunk_WT_cHiC_merge_10kb_hicTransform.cool
allFiles=""
for m_cool in ${path}/correct/E9.5*Del*merge*.cool; do
  python ${gitHubDirectory}/scripts/quantifyCooler.py \
    --bedpe $my_regions --cool1 $m_cool --cool2 $c_cool --raw \
    | awk -v name1=$(basename $m_cool .cool) -v name2=$(basename $c_cool .cool) 'BEGIN{print "Region\t"name1"\t"name2"\t"name1"relativeToWT\tpval_Mann_Whitney_U_test_"name1"\tpval_Wilcoxon_signed_rank_test"}NR>1{print}' > $(basename $m_cool .cool)VS$(basename $c_cool .cool)_$(basename ${my_regions} .bedpe).txt
  allFiles="$allFiles $(basename $m_cool .cool)VS$(basename $c_cool .cool)_$(basename ${my_regions} .bedpe).txt"
done
paste $allFiles | cut -f 1,4,5,10,11,16,17 > ${my_regions/.bedpe/_quantif.txt}


my_regions=${gitHubDirectory}/cHi-C/trunk_S8.bedpe
c_cool=${path}/toGEO/E9.5_trunk_WT_cHiC_merge_10kb.cool
allFiles=""
for m_cool in ${path}/toGEO/E9.5*Del*merge*.cool; do
  python ${gitHubDirectory}/scripts/quantifyCooler.py \
    --bedpe $my_regions --cool1 $m_cool --cool2 $c_cool \
    | awk -v name1=$(basename $m_cool .cool) -v name2=$(basename $c_cool .cool) 'BEGIN{print "Region\t"name1"\t"name2"\t"name1"relativeToWT\tpval_Mann_Whitney_U_test_"name1"\tpval_Wilcoxon_signed_rank_test"}NR>1{print}' > $(basename $m_cool .cool)VS$(basename $c_cool .cool)_$(basename ${my_regions} .bedpe).txt
  allFiles="$allFiles $(basename $m_cool .cool)VS$(basename $c_cool .cool)_$(basename ${my_regions} .bedpe).txt"
done
paste $allFiles | cut -f 1,4,5,10,11,16,17 > ${my_regions/.bedpe/_quantif.txt}

for my_regions in ${gitHubDirectory}/cHi-C/PFL*S*.bedpe; do
  c_cool=${path}/toGEO/E12.5_PFL_WT_cHiC_10kb.cool
  allFiles=""
  for m_cool in ${path}/toGEO/E12.5_PFL_*Del*.cool; do
    python ${gitHubDirectory}/scripts/quantifyCooler.py \
      --bedpe $my_regions --cool1 $m_cool --cool2 $c_cool \
      | awk -v name1=$(basename $m_cool .cool) -v name2=$(basename $c_cool .cool) 'BEGIN{print "Region\t"name1"\t"name2"\t"name1"relativeToWT\tpval_Mann_Whitney_U_test_"name1"\tpval_Wilcoxon_signed_rank_test"}NR>1{print}' > $(basename $m_cool .cool)VS$(basename $c_cool .cool)_$(basename ${my_regions} .bedpe).txt
    allFiles="$allFiles $(basename $m_cool .cool)VS$(basename $c_cool .cool)_$(basename ${my_regions} .bedpe).txt"
  done
  paste $allFiles | cut -f 1,4,5,10,11,16,17 > ${my_regions/.bedpe/_quantif.txt}
done
