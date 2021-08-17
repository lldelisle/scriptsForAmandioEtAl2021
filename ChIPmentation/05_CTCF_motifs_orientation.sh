
gitHubDirectory=""
pathForChIPGEO="/scratch/ldelisle/AmandioEtAl2021/ChIPM/toGEO/"
pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10
pathForWtFasta="${pathForFasta}/${genome}.fa"

# Go to a working directory
mkdir -p ${pathForChIPGEO}/../wd
cd ${pathForChIPGEO}/../wd

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda activate amandio2021

sample="E10.5_trunk_WT_ChIPM_CTCF"

for rep in 1 2 3; do
  # Summits around Hox clusters are extended
  cat ${pathForChIPGEO}/${sample}_rep${rep}.narrowPeak | awk -v OFS="\t" '($1=="chr2" && $2 < 75800000 && $3 > 73640001)||($1=="chr6" && $2 < 53186209 && $3 > 51264260)||($1=="chr11" && $2 < 96815432 && $3 > 96044583)||($1=="chr15" && $2 < 103235773 && $3 > 102634907){print $1, $2 + $10 - 250, $2 + $10 + 250, $4}' > ${sample}_rep${rep}_plotted_summitextended.bed
  # The fasta is extracted
  bedtools getfasta -fi ${pathForWtFasta} -bed ${sample}_rep${rep}_plotted_summitextended.bed > ${sample}_rep${rep}_plotted_summitextended.fa
done

# On the website of http://insulatordb.uthsc.edu/ in CTCFBS Prediction Tool the fasta is uploaded.
# The table output is copied into a file with extension _insdb_output.txt

for OUTPUT in *_insdb_output.txt; do
  # The orientation of the CTCF motif is deduced and put to . if score is negative
  # EMBL motifs are excluded as very short.
  cat $OUTPUT | grep -v EMBL | awk -v OFS="\t" '
NR > 1{
  split($3, a, ":|-")
  cur_chr = a[1]
  cur_start = a[2] + $4 + $5/2
  cur_motif = $1
  cur_score = $7
  cur_strand = $6
  if ($3 in scores_per_region){
    if (scores_per_region[$3] < cur_score){
      scores_per_region[$3] = cur_score
      infos_per_region[$3]["chr"] = cur_chr
      infos_per_region[$3]["start"] = cur_start
      infos_per_region[$3]["motif"] = cur_motif
      infos_per_region[$3]["score"] = cur_score
      infos_per_region[$3]["strand"] = cur_strand
    }
  } else {
    scores_per_region[$3] = cur_score
    infos_per_region[$3]["chr"] = cur_chr
    infos_per_region[$3]["start"] = cur_start
    infos_per_region[$3]["motif"] = cur_motif
    infos_per_region[$3]["score"] = cur_score
    infos_per_region[$3]["strand"] = cur_strand
  }
}
END{
  for (region in infos_per_region){
    name = infos_per_region[region]["chr"]"_"infos_per_region[region]["start"]
    if (!(name in written)){
      if (infos_per_region[region]["score"] < 0){
        cur_strand = "."
      } else {
        if (infos_per_region[region]["motif"] == "EMBL_M1"){
          if (infos_per_region[region]["strand"] == "+"){
            cur_strand = "-"
          } else {
            cur_strand = "+"
          }
        } else {
          cur_strand = infos_per_region[region]["strand"]
        }
      }
      printf "%s\t%d\t%d\t%s\t%f\t%s\n",
        infos_per_region[region]["chr"], \
        infos_per_region[region]["start"], \
        infos_per_region[region]["start"] + 1, \
        region"_"infos_per_region[region]["motif"], \
        infos_per_region[region]["score"], cur_strand
      written[name] = 1
    }
  }
}' | sort -k1,1 -k2,3n > ${OUTPUT/_insdb_output.txt/_noEMBL.bed}
done


# merge motifs: only the motifs that are present in all files
mkdir -p ${gitHubDirectory}/ChIPmentation/outputs
cat ${sample}_rep*_summitextended_noEMBL.bed | sort -k1,1  -k2,2n | awk -v nrep=3 '
{
    if(cur_chr == $1 && cur_start >= $2 - 1 && cur_start <= $2 + 1 && cur_strand == $6){
        timesFound=timesFound + 1
        if(cur_score < $5){
            bestLine = $0
        }
    } else {
        if (NR > 1 && timesFound == nrep){
            print bestLine
        }
        bestLine = $0
        cur_chr = $1
        cur_start = $2
        cur_end = $3
        cur_strand = $6
        timesFound = 1
    }
}' > ${gitHubDirectory}/ChIPmentation/outputs/${sample}_rep1and2and3.bed

# Create a bed with rgb field corresponding to motif orientation of CTCF:
f=${gitHubDirectory}/ChIPmentation/outputs/${sample}_rep1and2and3.bed
awk -F "\t" -v OFS="\t" -v colorPos="236,28,36" -v colorNeg="46,49,145" -v colorOther="0,0,0" '
{
  if ($6 == "+"){
    color = colorPos
  } else if ($6 == "-" ){
    color = colorNeg
  } else {
    color = colorOther
  }
  print $1, $2, $3, $4, $5, $6, $2, $2, color
}' ${f} > ${f/.bed/_colored.bed}

# Create a table with the sequence for SugFig4
for rep in 1 2 3; do
  awk '
NR == FNR{
  if ($1 == "chr2" && $3<74767377 && $2>74650810){
    split($4, a, "_")
    pos[a[1]]=$2
    motif[a[1]] = a[2]"_"a[3]; score[a[1]] = $5
  }
}
NR != FNR{
  if($3 in motif && $1 == motif[$3]){
    print pos[$3]"\t"$1"\t"$2"\t"$7
  }
}' ${f} ${sample}_rep${rep}_plotted_summitextended_insdb_output.txt
done | sort > ${gitHubDirectory}/ChIPmentation/outputs/ctcf_scores_seq.txt
