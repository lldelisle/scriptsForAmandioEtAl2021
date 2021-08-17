#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 32
#SBATCH --time 4:00:00
#SBATCH --array=1-4
#SBATCH --job-name ChIP
#SBATCH --chdir /scratch/ldelisle/AmandioEtAl2021/ChIP/

gitHubDirectory=$1

# This script remove adapters 
# make alignment with Bowtie2
# select MAPQ30 alignments
# Compute coverage SR200
# Normalize by nb of tags

path="$PWD/"
pathForTable="${gitHubDirectory}/ChIP/sra_table.txt"
pathForFastq="${path}/fastq/"
pathForIndex="/home/ldelisle/genomes/bowtie2/"
pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10

nbOfThreads=32

module purge
module load gcc/7.4.0 # required for bowtie2, samtools, star and bedtools
module load bowtie2/2.3.5
module load samtools/1.9
module load bedtools2/2.27.1
# cutadapt is version 1.16 working with python 3.6.1 built with intel 17.0.2

sample=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
sra=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
fastqFile=${sample}.fastq.gz
adapterSeq="TruSeq"

indexPath=${pathForIndex}${genome}

pathResults=${path}/${sample}/
echo $sample
mkdir -p $pathResults
cd $pathResults

# Cutadapt
if [ ! -e ${path}allFinalFiles/reports/${sample}_report-cutadapt.txt ]; then
  fastq=${pathForFastq}/$fastqFile
  if [ ! -e $fastq ]; then
    mkdir -p $pathForFastq
    cd $pathForFastq
    module load sra-toolkit/2.9.6
    fasterq-dump -o ${sample}.fastq ${sra}
    gzip ${sample}.fastq
    cd $pathResults
  fi
  if [ ! -s $fastq ]; then
    echo "FASTQ IS EMPTY"
    exit 1
  fi
  if [ $adapterSeq = "TruSeq" ]; then
    cutadapt -j $nbOfThreads -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 15 -o ${pathResults}${sample}-cutadapt.fastq.gz $fastq > ${pathResults}${sample}_report-cutadapt.txt
  else
    if [ $adapterSeq = "Nextera" ]; then
      cutadapt -j $nbOfThreads -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -q 30 -m 15 -o ${pathResults}${sample}-cutadapt.fastq.gz $fastq > ${pathResults}${sample}_report-cutadapt.txt
    else
      echo "YOU NEED TO WRITE THE CODE"
      exit 1
    fi
  fi
  mkdir -p ${path}allFinalFiles/reports/
  cp ${pathResults}${sample}_report-cutadapt.txt ${path}allFinalFiles/reports/
fi

# Mapping
if [ ! -e ${path}allFinalFiles/reports/${sample}_mapping_stats.txt ];then
  bowtie2 -p $nbOfThreads -x $indexPath -U ${pathResults}${sample}-cutadapt.fastq.gz 2> ${pathResults}${sample}_mapping_stats.txt  | samtools view --threads $nbOfThreads -Su - | samtools sort --threads $nbOfThreads -o ${pathResults}${sample}_mapped_sorted.bam
  mkdir -p ${path}allFinalFiles/reports/
  cp ${pathResults}${sample}_mapping_stats.txt ${path}allFinalFiles/reports/
fi

# MAPQ30
if [ ! -e ${pathResults}${sample}_mapped_sorted_q30.bam ]; then
  samtools view --threads $nbOfThreads -b ${pathResults}${sample}_mapped_sorted.bam -q 30 > ${pathResults}${sample}_mapped_sorted_q30.bam
fi

mkdir -p ${path}bedGraphs
mkdir -p ${path}bw
if [ ! -e ${pathForFasta}${genome}.fa.fai ]; then
  samtools faidx ${pathForFasta}${genome}.fa
fi

# Use macs2 with fixed fragment size of 200bp
if [ ! -e ${path}bedGraphs/${sample}_macs_SR200.bedGraph.gz ]; then
  macs2 callpeak -t ${sample}_mapped_sorted_q30.bam -n ${sample}_macs_SR200 --call-summits --nomodel --extsize 200 -B 2> ${pathResults}${sample}_macs_SR200.log
  cp ${pathResults}${sample}_macs_SR200.log ${path}allFinalFiles/reports/
  bash ${gitHubDirectory}/scripts/fromMacs2BdgToSimplifiedBdgAndBw.sh ${pathResults}${sample}_macs_SR200_treat_pileup.bdg ${path}bedGraphs/${sample}_macs_SR200 "macs2 SR200 of ${sample}" ${pathForFasta}${genome}.fa.fai &
fi
wait

if [ ! -e ${path}bedGraphs/${sample}_macs_SR200_norm.bedGraph.gz ]; then
  nbtags=$(grep "total tags in treatment" ${pathResults}${sample}_macs_SR200.log | awk '{print $NF}')
  zcat ${path}bedGraphs/${sample}_macs_SR200.bedGraph.gz | awk -v s=$sample -v n=$nbtags -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\""s" SR200 normalized by million tags\" visibility=full autoScale=on windowingFunction=mean"}NR>1{$4=$4/n*1e6; print}' > ${path}bedGraphs/${sample}_macs_SR200_norm.bedGraph
  bedGraphToBigWig ${path}bedGraphs/${sample}_macs_SR200_norm.bedGraph ${pathForFasta}${genome}.fa.fai ${path}bedGraphs/${sample}_macs_SR200_norm.bw
  gzip ${path}bedGraphs/${sample}_macs_SR200_norm.bedGraph &
fi

wait

mkdir -p ${path}/toGEO
cp ${path}bedGraphs/${sample}_macs_SR200_norm.bw ${path}/toGEO/${sample}.bw
cp ${pathResults}${sample}_macs_SR200_peaks.narrowPeak ${path}/toGEO/${sample}.narrowPeak
echo "Everything is done"
find . -size 0 -delete


