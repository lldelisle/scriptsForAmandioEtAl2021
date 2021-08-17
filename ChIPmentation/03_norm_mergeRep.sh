#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 20G
#SBATCH --cpus-per-task 1
#SBATCH --time 2:00:00
#SBATCH --array=1-11
#SBATCH --job-name mergeRep
#SBATCH --chdir /scratch/ldelisle/AmandioEtAl2021/ChIPM/

gitHubDirectory=$1

# This script perform average between replicates

path="$PWD/"
pathForTable="${gitHubDirectory}/ChIPmentation/replicates.txt"

pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10

module purge
module load gcc/7.4.0 # required for bowtie2, samtools, star and bedtools
module load bedtools2/2.27.1


sampleName=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
samples=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{split($2,a,",");for(j in a){print a[j]}}')
n=$(echo $samples | wc -w)
sample="${sampleName}_neq$n"

mkdir -p ${path}/${sample}

pathResults=${path}/${sample}/
echo $sample
cd $pathResults
allBDG=""
for s in $samples; do
  zcat ${path}bedGraphs/${s}_macs_SR200_norm.bedGraph.gz > normalized_${s}.bdg
  allBDG="$allBDG normalized_${s}.bdg"
done
if [ $n -ne 1 ]; then
  bedtools unionbedg -i $allBDG | awk -v n=$sample -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"mean of SR200 norm macs2 of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}{n=NF-3;sum=0;for(i=4;i<=NF;i++){sum+=$i};if(sum!=0){print $1,$2,$3,sum/n}}' > ${sample}.bedgraph
else
  cat $allBDG | awk -v n=$sample -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"mean of SR200 norm macs2 of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}$4!=0{print}' > ${sample}.bedgraph
fi
bedGraphToBigWig ${sample}.bedgraph ${pathForFasta}${genome}.fa.fai ${sample}.bw
gzip ${sample}.bedgraph

mkdir -p ${path}/toGEO
cp ${sample}.bw ${path}/toGEO
