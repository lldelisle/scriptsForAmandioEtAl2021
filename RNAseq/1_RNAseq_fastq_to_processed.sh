#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 36
#SBATCH --time 10:00:00
#SBATCH --array=1-36
#SBATCH --job-name RNAseqMerged102
#SBATCH --chdir /scratch/ldelisle/AmandioEtAl2021/RNAseq

gitHubDirectory=$1

# This script make alignment with STAR ENCODE parameters
# Evaluate FPKM with cufflinks
# Coverage with bedtools

path="$PWD/"
pathForFastq="$path/fastq/"
pathForTable="${gitHubDirectory}/RNAseq/sraTable.txt"
pathForFasta="/home/ldelisle/genomes/fasta/"
pathForIndex="/scratch/ldelisle/genomes/STARIndex_2.7.0e/"
genome=mm10
ensemblVersion=102
versionOfCufflinks="2.2.1"
# I tryed 36 and it failed:
# BAMoutput.cpp:27:BAMoutput: exiting because of *OUTPUT FILE* error: could not create output file ./_STARtmp//BAMsort/19/48
# SOLUTION: check that the path exists and you have write permission for this file. Also check ulimit -n and increase it to allow more open files.
nbOfThreads=16

module purge
module load gcc/7.4.0 #required for bowtie2, samtools, star and bedtools
module load openblas/0.3.6-openmp
module load r/3.6.0
module load star/2.7.0e
module load samtools/1.9
module load bedtools2/2.27.1
#cutadapt is working with python 3.6.1 built with intel 17.0.2 version 1.16

sample=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
sra=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
fastqFile=${sample}.fastq.gz
adapterSeq="TruSeq"

if [ ! -z $newGenome ]; then
  genome=$newGenome
fi

if [ $genome = "mm10" ]; then
  versionOfGtf="Mus_musculus.GRCm38.$ensemblVersion"
  gtfFile=${path}mergeOverlapGenesOfFilteredTranscriptsOf${versionOfGtf}_ExonsOnly_UCSC.gtf
else
  echo "unknown genome"
  exit 1
fi
indexPath=${pathForIndex}${genome}

#For the first one, The gtf and the file for DEXseq are prepared.
if [ $SLURM_ARRAY_TASK_ID == 1 ] && [ ! -e ${pathForFasta}${genome}.fa ];then
  cat ${pathForFasta}$genome/*.fa.gz > ${pathForFasta}${genome}.fa.gz
  gunzip ${pathForFasta}${genome}.fa.gz
fi
if [ $SLURM_ARRAY_TASK_ID == 1 ] && [ ! -e $gtfFile ];then
  wget "https://zenodo.org/record/4596490/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsOnly_UCSC.gtf.gz?download=1" -O ${gtfFile}.gz
  gunzip ${gtfFile}.gz
fi
if [ $SLURM_ARRAY_TASK_ID == 1 ] && [ ! -e ${path}/MTmouse.gtf ];then
  echo -e "chrM\tchrM_gene\texon\t0\t16299\t.\t+\t.\tgene_id \"chrM_gene_plus\"; transcript_id \"chrM_tx_plus\"; exon_id \"chrM_ex_plus\";">MTmouse.gtf
  echo -e "chrM\tchrM_gene\texon\t0\t16299\t.\t-\t.\tgene_id \"chrM_gene_minus\"; transcript_id \"chrM_tx_minus\"; exon_id \"chrM_ex_minus\";" >>MTmouse.gtf
fi

mkdir -p ${path}STAR/$sample

pathResults=${path}STAR/${sample}/
echo $sample
cd $pathResults

if [ ! -e ${path}allFinalFiles/reports/${sample}_report-cutadapt.txt ]; then
  fastq=${pathForFastq}/$fastqFile
  if [ ! -e $fastqR1 ]; then
    cd $pathForFastq
    fasterq-dump -o ${sample}.fastq ${sra}
    gzip ${sample}.fastq
    cd $pathResults
  fi
  if [ ! -s $fastq ]; then
    echo "FASTQ IS EMPTY"
    exit 1
  fi
  if [ $adapterSeq = "TruSeq" ]; then
    cutadapt -j $nbOfThreads -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 15 -o ${pathResults}${sample}-cutadapt.fastq.gz ${fastq} > ${pathResults}${sample}_report-cutadapt.txt
  else
    echo "YOU NEED TO WRITE THE CODE"
    exit 1
  fi
  mkdir -p ${path}allFinalFiles/reports/
  cp ${pathResults}${sample}_report-cutadapt.txt ${path}allFinalFiles/reports/
fi

if [ ! -e ${path}allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam ];then
  if [[ $sample = "ENCODE"* ]]; then 
    #I need to add --outSAMstrandField intronMotif because it is not stranded library
    STAR --runThreadN $nbOfThreads --genomeDir ${indexPath} --readFilesIn ${pathResults}${sample}-cutadapt.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate  --outSAMstrandField intronMotif  --sjdbOverhang '99' --sjdbGTFfile $gtfFile  --quantMode GeneCounts  --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1
  else
    STAR --runThreadN $nbOfThreads --genomeDir ${indexPath} --readFilesIn ${pathResults}${sample}-cutadapt.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate  --sjdbOverhang '99' --sjdbGTFfile $gtfFile  --quantMode GeneCounts  --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1
  fi
  mkdir -p ${path}allFinalFiles/reports/
  cp Log.final.out ${path}allFinalFiles/reports/${sample}_STAR_logFinal.txt
  mkdir -p ${path}allFinalFiles/bam/
  cp ${pathResults}Aligned.sortedByCoord.out.bam ${path}allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam
  samtools index ${path}allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam
fi

if [ -e ${pathForFasta}${genome}.fa ] && [ -e ${path}MTmouse.gtf ] && [ -e $gtfFile ] && [ -s ${pathResults}Aligned.sortedByCoord.out.bam ];then
 if [ ! -e ${path}allFinalFiles/FPKM_${sample}_isoforms.txt ];then
   echo "export PATH=$PATH:/home/ldelisle/softwares/cufflinks-${versionOfCufflinks}.Linux_x86_64" >cufflinks_${sample}.sh
   echo "mkdir -p ${pathResults}cufflinksWOMT" >>cufflinks_${sample}.sh
   echo "cufflinks -p 10 -o ${pathResults}cufflinksWOMT --max-bundle-length 10000000 --multi-read-correct --library-type \"fr-firststrand\" -b ${pathForFasta}${genome}.fa  --no-effective-length-correction -M ${path}MTmouse.gtf -G $gtfFile ${pathResults}Aligned.sortedByCoord.out.bam" >>cufflinks_${sample}.sh
   echo "" >>cufflinks_${sample}.sh
   echo "mkdir -p ${path}allFinalFiles" >>cufflinks_${sample}.sh
   echo "cp ${pathResults}cufflinksWOMT/genes.fpkm_tracking ${path}allFinalFiles/FPKM_${sample}.txt" >>cufflinks_${sample}.sh
   echo "cp ${pathResults}cufflinksWOMT/isoforms.fpkm_tracking ${path}allFinalFiles/FPKM_${sample}_isoforms.txt" >>cufflinks_${sample}.sh
   
   if [[ $sample = "ENCODE"* ]]; then
     sed -i 's/fr-firststrand/fr-unstranded/g' cufflinks_${sample}.sh
   fi
   echo "Launching cufflinks"
   bash cufflinks_${sample}.sh &
 fi
else
 echo "cufflinks not launch because some files are missing."
fi

if { [ ! -e accepted_hits_unique_${sample}.bam ] && [ -s ${pathResults}Aligned.sortedByCoord.out.bam ] ;} || [ -e tmp.header ] ;then
 echo "Compute uniquely aligned"
 samtools view -H Aligned.sortedByCoord.out.bam >  tmp.header
 samtools view -@ 5 Aligned.sortedByCoord.out.bam | grep  -w "NH:i:1" | cat tmp.header - | samtools view -@ 5 -b > accepted_hits_unique_${sample}.bam
 rm tmp.header
fi

if [ ! -e ${path}allFinalFiles/htseqCount_${sample}.txt ] && [ -e ReadsPerGene.out.tab ];then 
  mkdir -p ${path}allFinalFiles
  echo "write htseqCount"
  if [[ $sample = "ENCODE"* ]]; then
    cat ReadsPerGene.out.tab | awk '{print $1"\t"$2}' > ${path}allFinalFiles/htseqCount_${sample}.txt
  else
    cat ReadsPerGene.out.tab | awk '{print $1"\t"$4}' > ${path}allFinalFiles/htseqCount_${sample}.txt
  fi
fi

#This is to make the bedgraph of coverage
mkdir -p ${path}bedGraphs
if { [ ! -e ${path}bedGraphs/${sample}_Uniq_fwd.bedGraph.gz ] && [ -e accepted_hits_unique_${sample}.bam ];} || [ -e tmp.header.uf ] ;then
  #strand + corresponds to reverse strand due to TruSeq
  echo "Building uniq fwd reads bedGraph"
  echo "track type=bedGraph name=\"${sample} Uniq reads forward\"">tmp.header.uf
  echo "bedtools genomecov -ibam accepted_hits_unique_${sample}.bam -bg -split -strand - | cat tmp.header.uf - > ${path}bedGraphs/${sample}_Uniq_fwd.bedGraph">uf.sh
  echo "cat ${path}bedGraphs/${sample}_Uniq_fwd.bedGraph | grep -v track | LC_COLLATE=C sort -k1,1 -k2,2n > ${path}bedGraphs/${sample}_Uniq_fwd_sorted.bedGraph">>uf.sh
  echo "bedGraphToBigWig ${path}bedGraphs/${sample}_Uniq_fwd_sorted.bedGraph ${pathForFasta}${genome}.fa.fai ${path}bedGraphs/${sample}_Uniq_fwd.bw">>uf.sh
  echo "gzip ${path}bedGraphs/${sample}_Uniq_fwd.bedGraph">>uf.sh
  echo "rm tmp.header.uf">>uf.sh
  echo "touch uf.done">>uf.sh
  bash uf.sh &
fi
if { [ ! -e ${path}bedGraphs/${sample}_Uniq_rev.bedGraph.gz ] && [ -e accepted_hits_unique_${sample}.bam ];} || [ -e tmp.header.ur ];then
  echo "Building uniq rev reads bedGraph"
  echo "track type=bedGraph name=\"${sample} Uniq reads reverse\"">tmp.header.ur
  echo "bedtools genomecov -ibam accepted_hits_unique_${sample}.bam -bg -split -strand + | cat tmp.header.ur - > ${path}bedGraphs/${sample}_Uniq_rev.bedGraph">ur.sh
  echo "cat ${path}bedGraphs/${sample}_Uniq_rev.bedGraph | grep -v track | LC_COLLATE=C sort -k1,1 -k2,2n > ${path}bedGraphs/${sample}_Uniq_rev_sorted.bedGraph">>ur.sh
  echo "bedGraphToBigWig ${path}bedGraphs/${sample}_Uniq_rev_sorted.bedGraph ${pathForFasta}${genome}.fa.fai ${path}bedGraphs/${sample}_Uniq_rev.bw">>ur.sh
  echo "gzip ${path}bedGraphs/${sample}_Uniq_rev.bedGraph">>ur.sh
  echo "rm tmp.header.ur">>ur.sh
  echo "touch ur.done">>ur.sh
  bash ur.sh &
fi


wait
echo "Everything is done"
find . -size 0 -delete

mkdir -p ${path}/toGEO
cp ${path}bedGraphs/${sample}_Uniq_rev.bw ${path}/toGEO
cp ${path}bedGraphs/${sample}_Uniq_fwd.bw ${path}/toGEO
