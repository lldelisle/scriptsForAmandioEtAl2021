#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 32
#SBATCH --time 1:00:00
#SBATCH --job-name quantif
#SBATCH --chdir /scratch/ldelisle/AmandioEtAl2021/ChIPM/

gitHubDirectory=$1

path="$PWD/"
pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10

module purge
module load gcc/7.4.0 # required for bowtie2, samtools, star and bedtools
module load samtools/1.9
module load bedtools2/2.27.1

inputBED="CTCF_HoxD_500bp.bed"

if [ ! -e $inputBED ]; then
  echo -e "chr2\t73644685\t75724063" > regionPlot.bed
  bedtools intersect -a ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3.bed -b regionPlot.bed -wa -u > CTCF_regionPlot.bed
  bedtools slop -i CTCF_regionPlot.bed -g ${pathForFasta}${genome}.fa.fai -b 250 > ${inputBED}
fi

bwFiles="$(ls ${path}/toGEO/*RAD21*MAnorm.bw)"

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda activate amandio2021

if [ ! -e results.txt ]; then
  multiBigwigSummary BED-file --BED $inputBED -b $bwFiles -out results.npz -p 32 -v --outRawCounts results.txt
fi

# R version 4.1.0
Rscript ${gitHubDirectory}/ChIPmentation/plotQuantifMAnorm.R ${gitHubDirectory}
