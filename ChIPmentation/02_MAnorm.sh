#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 1
#SBATCH --time 2:00:00
#SBATCH --array=1-26
#SBATCH --job-name MAnorm
#SBATCH --chdir /scratch/ldelisle/AmandioEtAl2021/ChIPM/

gitHubDirectory=$1

path="$PWD/"
pathForTable="${gitHubDirectory}/ChIPmentation/MAnorm.txt"
pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10

module purge
module load gcc/7.4.0 #required for bowtie2, samtools, star and bedtools
module load bedtools2/2.27.1

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda activate amandio2021

nbOfThreads=1

myExpeName=$(cat $pathForTable | awk -F "\t" -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
name1=$(cat $pathForTable | awk -F "\t" -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
name2=$(cat $pathForTable | awk -F "\t" -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}')
w=$(cat $pathForTable | awk -F "\t" -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $4}')
p1=$(cat $pathForTable | awk -F "\t" -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $5}')
p2=$(cat $pathForTable | awk -F "\t" -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $6}')
bg1gz=$(cat $pathForTable | awk -F "\t" -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $7}')
bg2gz=$(cat $pathForTable | awk -F "\t" -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $8}')

if [[ "$name1" = *"HH"* ]]; then
  genome=galGal6
fi

if [ -z $w ]; then
  w=100
fi

if [ -z "$p1" ]; then
  if [ -e ${path}${name1}/${name1}_macs_SR200_peaks.narrowPeak ]; then
    p1=${path}${name1}/${name1}_macs_SR200_peaks.narrowPeak
  else
    echo "No p1"
    exit 1
  fi
fi

if [ -z "$p2" ]; then
  if [ -e ${path}${name2}/${name2}_macs_SR200_peaks.narrowPeak ]; then
    p2=${path}${name2}/${name2}_macs_SR200_peaks.narrowPeak
  else
    echo "No p2"
    exit 1
  fi
fi


if [ -e ${path}${name1}/${name1}_mapped_sorted_q30_rmdup.bed ]; then
  r1=${path}${name1}/${name1}_mapped_sorted_q30_rmdup.bed
else
  if [ -e ${path}/input/${name1}.bed ]; then
    r1=${path}/input/${name1}.bed
  else
    echo "No r1"
    exit 1
  fi
fi

if [ -e ${path}${name2}/${name2}_mapped_sorted_q30_rmdup.bed ]; then
  r2=${path}${name2}/${name2}_mapped_sorted_q30_rmdup.bed
else
  if [ -e ${path}/input/${name2}.bed ]; then
    r2=${path}/input/${name2}.bed
  else
    echo "No r2"
    exit 1
  fi
fi

if [ -z $bg1gz ]; then
  if [ -e ${path}bedGraphs/${name1}_macs_SR200.bedGraph.gz ]; then
    bg1gz=${path}bedGraphs/${name1}_macs_SR200.bedGraph.gz
  else
    echo "No bg1gz"
    exit 1
  fi
fi
namebg1=$(basename $bg1gz .bedGraph.gz)

if [ -z $bg2gz ]; then
  if [ -e ${path}bedGraphs/${name2}_macs_SR200.bedGraph.gz ]; then
    bg2gz=${path}bedGraphs/${name2}_macs_SR200.bedGraph.gz
  else
    echo "No bg2gz"
    exit 1
  fi
fi

namebg2=$(basename $bg2gz .bedGraph.gz)
    
if [ -e ${r1}_splitReads.bed ]; then
    r1=${r1}_splitReads.bed
fi
if [ -e ${r2}_splitReads.bed ]; then
    r2=${r2}_splitReads.bed
fi

# When BAMPE is specified bam to bed provide 1 line per pair and 1 need 1 line per read.
for file in $r1 $r2; do
    nbField=$(head -n1 $file | awk '{print NF}')
    if [ $nbField -eq 10 ]; then
        cat $file | awk -v OFS="\t" '{
            print $1,$2,$3,$7"/1",$8,$9
            print $4,$5,$6,$7"/2",$8,$9
        }' > ${file}_splitReads.bed
    fi
done
if [ -e ${r1}_splitReads.bed ]; then
    r1=${r1}_splitReads.bed
fi
if [ -e ${r2}_splitReads.bed ]; then
    r2=${r2}_splitReads.bed
fi

# There are duplicated entries for peak with multiple summits
# The expected input of MAnorm is chr start end summit_abs_pos (not relative)
if [ ! -e ${p1}_filteredconvertedForMAnorm.bed ]; then
  cat $p1 | awk -v OFS="\t" 'BEGIN{best=""}!($4~/[a-z]$/){if(best!=""){print best;best=""};print $1, $2, $3, $2+$10}$4~/[a-z]$/{if($4~/a$/){if(best!=""){print best};best=$1"\t"$2"\t"$3"\t"$2+$10;bestScore=$5}else{if($5>bestScore){bestScore=$5;best=$1"\t"$2"\t"$3"\t"$2+$10}}}' > ${p1}_filteredconvertedForMAnorm.bed
fi
p1="${p1}_filteredconvertedForMAnorm.bed"
if [ ! -e ${p2}_filteredconvertedForMAnorm.bed ]; then
  cat $p2 | awk -v OFS="\t" 'BEGIN{best=""}!($4~/[a-z]$/){if(best!=""){print best;best=""};print $1, $2, $3, $2+$10}$4~/[a-z]$/{if($4~/a$/){if(best!=""){print best};best=$1"\t"$2"\t"$3"\t"$2+$10;bestScore=$5}else{if($5>bestScore){bestScore=$5;best=$1"\t"$2"\t"$3"\t"$2+$10}}}' > ${p2}_filteredconvertedForMAnorm.bed
fi
p2="${p2}_filteredconvertedForMAnorm.bed"

if [ ! -e manorm_${myExpeName}.log ]; then
    manorm --p1 "$p1" --p2 "$p2" --r1 "$r1" --r2 "$r2" -s -w $w --name1 $name1 --name2 $name2 -o $myExpeName 2> manorm_${myExpeName}.log
fi

if [ ! -e ${pathForFasta}${genome}.fa.fai ]; then
  samtools faidx ${pathForFasta}${genome}.fa
fi
fileWithSizes=${pathForFasta}${genome}.fa.fai

if [ ! -e ${namebg1}_withPseudoCount_normedby${name2}_in${myExpeName}.bw ]; then
    intersect=$(cat manorm_${myExpeName}.log | awk '/M-A model/{match($0, /([\-0-9.]+) \* A ([+-]) ([0-9.]+)/);split(substr($0, RSTART, RLENGTH),a," ");print a[4]a[5]}')
    slope=$(cat manorm_${myExpeName}.log | awk '/M-A model/{match($0, /([\-0-9.]+) \* A ([+-]) ([0-9.]+)/);split(substr($0, RSTART, RLENGTH),a," ");print a[1]}')
    pseudoCount=$(cat  manorm_${myExpeName}.log | awk '/Peak width to extend from summit/{print 1000/(2*$NF)}')
    zcat "$bg2gz" | awk -F "\t" -v OFS="\t" -v ps=$pseudoCount -v name=$name2 -v me=$myExpeName 'BEGIN{
        print "track type=bedGraph name=\"macs2 profile of "name" with pseudocount in "me"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"
    }
    NF==4{
        $4=ps*(1+$4);print
    }' > ${namebg2}_withPseudoCount_in${myExpeName}.bedGraph
    bedGraphToBigWig ${namebg2}_withPseudoCount_in${myExpeName}.bedGraph $fileWithSizes ${namebg2}_withPseudoCount_in${myExpeName}.bw
    gzip ${namebg2}_withPseudoCount_in${myExpeName}.bedGraph &

    bedtools unionbedg -i "$bg1gz" "$bg2gz" | awk -v i=$intersect -v s=$slope -v OFS="\t" -v name=$name1 -v name2=$name2 -v ps=$pseudoCount -v me=$myExpeName 'BEGIN{
        print "track type=bedGraph name=\"macs2 profile of "name" with pseudocount normalized with "name2" in "me"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"
    }
    $4>0{
        $4=ps*(1+$4)
        $5=ps*(1+$4)
        log2NormValue=log($4)/log(2) - (i + s*log($4*$5)/(2*log(2)))
        print $1,$2,$3,exp(log2NormValue*log(2))
    }' > ${namebg1}_withPseudoCount_normedby${name2}_in${myExpeName}.bedGraph

    bedGraphToBigWig ${namebg1}_withPseudoCount_normedby${name2}_in${myExpeName}.bedGraph $fileWithSizes ${namebg1}_withPseudoCount_normedby${name2}_in${myExpeName}.bw
    gzip ${namebg1}_withPseudoCount_normedby${name2}_in${myExpeName}.bedGraph &
fi
mkdir -p ${path}/toGEO
cp ${namebg1}_withPseudoCount_normedby${name2}_in${myExpeName}.bw ${path}/toGEO/${name1}_MAnorm.bw
cp ${namebg2}_withPseudoCount_in${myExpeName}.bw ${path}/toGEO/${name2}_MAnorm.bw

wait
