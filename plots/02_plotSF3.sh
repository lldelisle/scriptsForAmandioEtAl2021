cd /scratch/ldelisle/AmandioEtAl2021/plots/
gitHubDirectory=""
mainPath="$PWD/../"
pathChIPM="${mainPath}/ChIPM/"

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda activate amandio2021

# Download data for Chicken:
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3182nnn/GSM3182452/suppl/GSM3182452%5FHH20FL%5FCTCF%2Ebedgraph%2Egz"
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3182nnn/GSM3182452/suppl/GSM3182452%5FHH20FL%5FCTCF%5FnarrowPeaks%2Ebed%2Egz"

# Download data for Human:
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733672/suppl/GSM733672%5Fhg19%5FwgEncodeBroadHistoneH1hescCtcfStdPk%2EbroadPeak%2Egz"
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733672/suppl/GSM733672%5Fhg19%5FwgEncodeBroadHistoneH1hescCtcfStdSig%2EbigWig"

# Go to a working directory
mkdir -p wd
cd wd

# Get corresponding fasta:
wget "https://hgdownload.soe.ucsc.edu/goldenPath/galGal5/bigZips/galGal5.fa.gz"
wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr2.fa.gz"
gunzip galGal5.fa.gz
gunzip chr2.fa.gz

# For the chicken:
sample="GSM3182452_HH20FL_CTCF"
# Summits around HoxD cluster are extended
zcat ../${sample}_narrowPeaks.bed.gz | awk -v OFS="\t" '($1=="chr7" && $2 < 16790000 && $3 > 15696357){print $1, $2 + $10 - 250, $2 + $10 + 250, $4}' > ${sample}_plotted_summitextended.bed
# The fasta is extracted
bedtools getfasta -fi galGal5.fa -bed ${sample}_plotted_summitextended.bed > ${sample}_plotted_summitextended.fa

# For the human:
sample="GSM733672_hg19_wgEncodeBroadHistoneH1hescCtcf"
# Peaks around HoxD cluster are extended or shrinked to 500bp
zcat ../${sample}StdPk.broadPeak.gz | awk -v OFS="\t" '($1=="chr2" && $2 < 178663278 && $3 > 175313779){printf("%s\t%d\t%d\t%s\n", $1, ($2 + $3) / 2 - 250, ($2 + $3) / 2 + 250, $4)}' > ${sample}_plotted_summitextended.bed
# The fasta is extracted
bedtools getfasta -fi chr2.fa -bed ${sample}_plotted_summitextended.bed > ${sample}_plotted_summitextended.fa

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
}' | sort -k1,1 -k2,3n > ${gitHubDirectory}/plots/annotations/${OUTPUT/_insdb_output.txt/.bed}
done

# For the human:
sample="GSM733672_hg19_wgEncodeBroadHistoneH1hescCtcf"
# Create a bed with rgb field corresponding to motif orientation of CTCF:
f=${gitHubDirectory}/plots/annotations/${sample}.bed
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

# For the chicken:
sample="GSM3182452_HH20FL_CTCF"
# Create a bed with rgb field corresponding to motif orientation of CTCF (colors are inverted):
f=${gitHubDirectory}/plots/annotations/${sample}.bed
awk -F "\t" -v OFS="\t" -v colorNeg="236,28,36" -v colorPos="46,49,145" -v colorOther="0,0,0" '
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

cd ..

# Convert bedgraph to bigwig:
gunzip GSM3182452_HH20FL_CTCF.bedgraph.gz
bedGraphToBigWig GSM3182452_HH20FL_CTCF.bedgraph galGal5.fa.fai GSM3182452_HH20FL_CTCF.bw

# Download gtf:
wget "http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"
wget "http://ftp.ensembl.org/pub/release-94/gtf/gallus_gallus/Gallus_gallus.Gallus_gallus-5.0.94.gtf.gz"

# Extract protein coding:
zcat Homo_sapiens.GRCh37.75.gtf.gz | awk '$3 == "exon" {print}' | grep protein_coding | grep -v "HOXD3-002" | grep -v "HOXD3-005" | grep -v "HOXD3-003" | grep -v "HOXD10-002" | grep -v "HOXD10-003" | grep -v "HOXD11-003" > Human_prot_coding.gtf
cat Human_prot_coding.gtf | awk '($1=="2" && $4 < 178663278 && $5 > 175313779){print}' > Human_prot_coding_around_HoxD.gtf
zcat Gallus_gallus.Gallus_gallus-5.0.94.gtf.gz | grep protein_coding > Chicken_prot_coding.gtf
cat Chicken_prot_coding.gtf | awk '($1=="7" && $4 < 16790000 && $5 > 15696357){print}' > Chicken_prot_coding_around_HoxD.gtf
for specie in Human Chicken; do
  python ${gitHubDirectory}/scripts/fromgtfTobed12.py --output ${specie}_prot_coding_around_HoxD.bed --mergeTranscriptsAndOverlappingExons ${specie}_prot_coding_around_HoxD.gtf
  cat ${specie}_prot_coding_around_HoxD.bed | awk -v OFS="\t" '
BEGIN{
  hoxd=0
}
{
  if ($4~/HOXD/){
    if (hoxd==0){
      start=$2
      hoxd=1
    } else{
      end=$3
    }
  } else {
    if (hoxd == 1){
      print $1,start,end,"HOXD","0",".",start,start,"0,0,0"
      hoxd=0
    }
    print $1,$2,$3,$4,"0",".",$2,$2,"166,168,171"
  }
}' > ${gitHubDirectory}/plots/annotations/${specie}_simplified_genes_around_HoxD.bed

  # Split the gtf in Hox and non Hox:
  grep -i "Hox" ${specie}_prot_coding.gtf > ${specie}_Hox.gtf
  grep -i -v "Hox" ${specie}_prot_coding.gtf > ${specie}_nonHox.gtf
done

# Figure S3
bws=("${pathChIPM}/toGEO/E10.5_trunk_WT_ChIPM_CTCF_MAnorm_neq3.bw" 'GSM733672_hg19_wgEncodeBroadHistoneH1hescCtcfStdSig.bigWig' 'GSM3182452_HH20FL_CTCF.bw')
ctcfs=("${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed" "${gitHubDirectory}/plots/annotations/GSM733672_hg19_wgEncodeBroadHistoneH1hescCtcf_colored.bed" "${gitHubDirectory}/plots/annotations/GSM3182452_HH20FL_CTCF_colored.bed")
species=('' 'Human_' 'Chicken_')
regionsZO=('chr2:73779626-75669724' 'chr2:175895200-178092461' 'chr7:15889800-16780000')
maxsZO=('1610' '400' '60')
regionsZI=('chr2:74650810-74767377' 'chr2:176941129-177057814' 'chr7:16323377-16402356')
maxsZI=('1000' '400' '60')
for i in 0 1 2; do
  specie=${species[$i]}
  if [ "${specie}" = "Chicken_" ]; then
    flag="--decreasingXAxis"
  else
    flag=""
  fi
  ini_file=FigS3${specie}_zoomout.ini
  scalebarPos=$(echo ${regionsZO[$i]} | awk '{split($0, a, "-"); print a[2] - 50000}')
  if [ "${specie}" = "Chicken_" ]; then
    scalebarPos=$(echo ${regionsZO[$i]} | awk '{split($0, a, ":|-"); print a[2] + 50000}')
  fi
  echo "[spacer]
height = 3.6

[scalebar]
file_type = scalebar
height = 0.3
where = top
x_center = $scalebarPos
size = 50000
fontsize = 8

[bigwig]
file = ${bws[$i]}
height = 3
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = ${maxsZO[$i]}
summary_method = max
type= line
file_type = bigwig

[spacer]
height = 0.15

[CTCF peaks]
file = ${ctcfs[$i]}
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5

[spacer]
height = 0.05

[genes]
file = ${gitHubDirectory}/plots/annotations/${specie}simplified_genes_around_HoxD.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.7

[genes]
file = ${gitHubDirectory}/plots/annotations/${specie}simplified_genes_around_HoxD.bed
display = collapsed
color = none
border_color = none
labels = true
height = 0.7
fontsize = 5
" > ${ini_file}

  pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.pdf} --region ${regionsZO[$i]} --width 20.5 ${flag} --height 4

  ini_file=FigS3${specie}_zoomin.ini
  scalebarPos=$(echo ${regionsZI[$i]} | awk '{split($0, a, "-"); print a[2] - 5000}')
  if [ "${specie}" = "Chicken_" ]; then
    scalebarPos=$(echo ${regionsZI[$i]} | awk '{split($0, a, ":|-"); print a[2] + 5000}')
  fi
  echo "[spacer]
height = 1.6

[scalebar]
file_type = scalebar
height = 0.3
where = top
x_center = ${scalebarPos}
size = 5000
fontsize = 8

[bigwig]
file = ${bws[$i]}
height = 3
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = ${maxsZI[$i]}
summary_method = max
type= line
file_type = bigwig

[spacer]
height = 0.15

[CTCF peaks]
file = ${ctcfs[$i]}
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5

[spacer]
height = 0.05


[Hox]
file = ${specie}Hox.gtf
merge_transcripts = true
color = black
border_color = none
height = 0.4
display = collapsed
labels = False
style = flybase 
color_utr = black

[nonHox]
file = ${specie}nonHox.gtf
merge_transcripts = true
color = grey
border_color = none
display = collapsed
labels = False
style = flybase 
color_utr = grey
color_backbone = grey
overlay_previous = yes

[spacer]
height = 0.2

[gene names]
file = ${specie}prot_coding.gtf
merge_transcripts = true
prefered_name = gene_name
color = none
border_color = none
color_utr = none
fontsize = 5
height = 0.3
display = collapsed
style = UCSC
line_width = 0
" > ${ini_file}

  pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.pdf} --region ${regionsZI[$i]} --width 20.5 ${flag} --height 3
done


pdflatex -jobname FigS3 &> /dev/null <<-EOF
\\documentclass[preview]{standalone}
\\usepackage{graphicx}
\\begin{document}
\\includegraphics{FigS3_zoomout.pdf}
\\includegraphics{FigS3_zoomin.pdf}
\\includegraphics{FigS3Human__zoomout.pdf}
\\includegraphics{FigS3Human__zoomin.pdf}
\\includegraphics{FigS3Chicken__zoomout.pdf}
\\includegraphics{FigS3Chicken__zoomin.pdf}
\\end{document}
EOF




for i in 0 1 2; do
  specie=${species[$i]}
  if [ "${specie}" = "Chicken_" ]; then
    flag="--decreasingXAxis"
  else
    flag=""
  fi
  ini_file=FigS3${specie}_zoomout_Rita.ini
  scalebarPos=$(echo ${regionsZO[$i]} | awk '{split($0, a, "-"); print a[2] - 50000}')
  if [ "${specie}" = "Chicken_" ]; then
    scalebarPos=$(echo ${regionsZO[$i]} | awk '{split($0, a, ":|-"); print a[2] + 50000}')
  fi
  echo "[spacer]
height = 3.6

[scalebar]
file_type = scalebar
height = 0.3
where = top
x_center = $scalebarPos
size = 50000
fontsize = 8

[bigwig]
file = ${bws[$i]}
height = 3
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = ${maxsZO[$i]}
summary_method = max
type= line
file_type = bigwig

[spacer]
height = 0.15

[CTCF peaks]
file = ${ctcfs[$i]}
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5

[spacer]
height = 0.05

[genes]
file = ${gitHubDirectory}/plots/annotations/${specie}simplified_genes_around_HoxD.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.7

[spacer]
height = 0.7
" > ${ini_file}

  pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.pdf} --region ${regionsZO[$i]} --width 20.5 ${flag} --height 4

  ini_file=FigS3${specie}_zoomin_Rita.ini
  scalebarPos=$(echo ${regionsZI[$i]} | awk '{split($0, a, "-"); print a[2] - 5000}')
  if [ "${specie}" = "Chicken_" ]; then
    scalebarPos=$(echo ${regionsZI[$i]} | awk '{split($0, a, ":|-"); print a[2] + 5000}')
  fi
  echo "[spacer]
height = 1.6

[scalebar]
file_type = scalebar
height = 0.3
where = top
x_center = ${scalebarPos}
size = 5000
fontsize = 8

[bigwig]
file = ${bws[$i]}
height = 3
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = ${maxsZI[$i]}
summary_method = max
type= line
file_type = bigwig

[spacer]
height = 0.15

[CTCF peaks]
file = ${ctcfs[$i]}
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5

[spacer]
height = 0.05


[Hox]
file = ${specie}Hox.gtf
merge_transcripts = true
color = black
border_color = none
height = 0.4
display = collapsed
labels = False
style = flybase 
color_utr = black

[nonHox]
file = ${specie}nonHox.gtf
merge_transcripts = true
color = grey
border_color = none
display = collapsed
labels = False
style = flybase 
color_utr = grey
color_backbone = grey
overlay_previous = yes

[spacer]
height = 0.2

[spacer]
height = 0.3
" > ${ini_file}

  pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.pdf} --region ${regionsZI[$i]} --width 20.5 ${flag} --height 3
done


pdflatex -jobname FigS3_ForRita &> /dev/null <<-EOF
\\documentclass[preview]{standalone}
\\usepackage{graphicx}
\\begin{document}
\\includegraphics{FigS3_zoomout_Rita.pdf}
\\includegraphics{FigS3_zoomin_Rita.pdf}
\\includegraphics{FigS3Human__zoomout_Rita.pdf}
\\includegraphics{FigS3Human__zoomin_Rita.pdf}
\\includegraphics{FigS3Chicken__zoomout_Rita.pdf}
\\includegraphics{FigS3Chicken__zoomin_Rita.pdf}
\\end{document}
EOF
