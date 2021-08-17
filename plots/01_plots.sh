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
#SBATCH --chdir /scratch/ldelisle/AmandioEtAl2021/plots/

gitHubDirectory=$1

path="$PWD/"
mainPath="$PWD/../"
pathRNAseq="${mainPath}/RNAseq/"
pathCHiC="${mainPath}/cHi-C/"
pathChIPM="${mainPath}/ChIPM/"
pathChIP="${mainPath}/ChIP/"

module purge

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda activate amandio2021

# From gtf deduce a simplified gene annotations
gtfFile="${pathRNAseq}/merge*.gtf"
grep protein_coding ${gtfFile} | grep -v "Hoxa3-201" | grep -v "Hoxd3-203" > prot_coding.gtf
cat prot_coding.gtf | awk '$1=="chr2" && $4 < 77000000 && $5 > 73000000{print}' > prot_coding_around_HoxD.gtf
python ${gitHubDirectory}/scripts/fromgtfTobed12.py --output prot_coding_around_HoxD.bed --mergeTranscriptsAndOverlappingExons prot_coding_around_HoxD.gtf
cat prot_coding_around_HoxD.bed | awk -v OFS="\t" '
BEGIN{
  hoxd=0
}
{
  if ($4~/Hoxd/){
    if (hoxd==0){
      start=$2
      hoxd=1
    } else{
      end=$3
    }
  } else {
    if (hoxd == 1){
      print $1,start,end,"HoxD","0",".",start,start,"0,0,0"
      hoxd=0
    }
    print $1,$2,$3,$4,"0",".",$2,$2,"166,168,171"
  }
}' > ${gitHubDirectory}/plots/annotations/simplified_genes_around_HoxD.bed

# Split the gtf in Hox and non Hox:
grep "Hox" prot_coding.gtf > Hox.gtf
grep -v "Hox" prot_coding.gtf > nonHox.gtf

# Get the left pos:
cat prot_coding_around_HoxD.bed | awk -v OFS="\t" '{print $1, $2, $2 +1, $4}' > prot_coding_around_HoxD_start.bed

# More roughtly for FigS2:
cat prot_coding.gtf | grep "exon_number \"1\"" | awk -F "\t" -v OFS="\t" '$7 == "-"{$4=$5; print}$7 == "+"{$5=$4; print}' > tss.gtf


# Get HoxD and CS3840:
grep HoxD ${gitHubDirectory}/plots/annotations/simplified_genes_around_HoxD.bed > HoxD_CS3840.bed
echo -e "chr2\t75137845\t75152903\tCS38-40\t0\t.\t75137845\t75137845\t166,168,171" >> HoxD_CS3840.bed

# For FigS2 we quantify the signal in peaks:
ln -s ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed .
ctcf_motifs=E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
multiBigwigSummary BED-file -b ${pathChIPM}/toGEO/E10.5_trunk_WT_ChIPM_CTCF_MAnorm_neq3.bw -o tmp.npz --BED $ctcf_motifs --outRawCounts tmp.txt
awk -v OFS="\t" '
NR>1{
  if(NR == FNR){
    # This is the output of multiBigwigSummary
    scores[$1"_"$2] = $4
  } else {
    $4 = $4"_"$5
    $5 = scores[$1"_"$2]
    print
  }
}' tmp.txt $ctcf_motifs > ${ctcf_motifs/.bed/_cov.bed}

cat  ${ctcf_motifs/.bed/_cov.bed} | awk '$6 == "+"{print}' > ${ctcf_motifs/.bed/_cov_pos.bed}
cat  ${ctcf_motifs/.bed/_cov.bed} | awk '$6 == "-"{print}' > ${ctcf_motifs/.bed/_cov_neg.bed}

echo -e "chr11\t96189088\t96376579
chr15\t102916367\t103042207
chr2\t74650810\t74767377" > Hox_pos.bed

echo -e "chr6\t52151994\t52278047" > Hox_neg.bed

# Common ini_files:
echo ""

# Figure 1
ini_file="Fig1A.ini"
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 74000000
size = 200000

[CHiC_E9.5_WT_10kb]
file = ${pathCHiC}/toGEO/E9.5_trunk_WT_cHiC_merge_10kb.cool
colormap = Greys
rasterize = false
depth = 1200000
height = 10
min_value = 0
max_value = 0.02
show_masked_bins = false
file_type = hic_matrix

[spacer]
height = 0.25

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.7

[genes]
file = ${gitHubDirectory}/plots/annotations/simplified_genes_around_HoxD.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.7
" > ${ini_file}

pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.png} --region chr2:73779626-75669724 --width 20 --height 6 --dpi 500

ini_file="Fig1B.ini"
echo "[Genes_file]
file = prot_coding.gtf
color = black
border_color = none
height = 1
display = collapsed
labels = False
style = flybase 
color_utr = black
" > ${ini_file}

pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.pdf} --region chr2:74665763-74767377 --width 18 --height 0.5 

ini_file="Fig1C.ini"
echo "[scalebar]
file_type = scalebar
height = 0.3
where = top
x_center = 74657000
size = 10000
fontsize = 8

[CTCF_WT]
file = ${pathChIPM}/toGEO/E10.5_trunk_WT_ChIPM_CTCF_MAnorm_neq3.bw
height = 3
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 1000
summary_method = max

[spacer]
height = 0.15

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.3

[spacer]
height = 0.05

[Hox]
file = Hox.gtf
merge_transcripts = true
color = black
border_color = none
height = 0.3
display = collapsed
labels = False
style = flybase 
color_utr = black

[nonHox]
file = nonHox.gtf
merge_transcripts = true
color = grey
border_color = none
height = 0.3
display = collapsed
labels = False
style = flybase 
color_utr = grey
color_backbone = grey
overlay_previous = yes

[spacer]
height = 0.2

[gene names]
file = prot_coding_around_HoxD_start.bed
color = none
border_color = none
color_utr = none
fontsize = 5
height = 0.3
display = collapsed
style = UCSC
line_width = 0
" > ${ini_file}

pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.pdf} --region chr2:74650810-74767377 --width 20 --height 2.5

## Figure 2
ini_file="Fig2A.ini"
echo "[gene names]
file = prot_coding_around_HoxD_start.bed
color = none
border_color = none
color_utr = none
fontsize = 5
height = 0.3
display = collapsed
style = UCSC
line_width = 0

[spacer]
height = 0.2

[Hox]
file = Hox.gtf
merge_transcripts = true
color = black
border_color = none
height = 0.4
display = collapsed
labels = False
style = flybase 
color_utr = black

[nonHox]
file = nonHox.gtf
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
height = 0.25

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.8


[spacer]
height = 0.25

[scalebar]
file_type = scalebar
height = 0.2
where = top
x_center = 74657000
size = 10000
fontsize = 6

[CTCF_WT]
file = ${pathChIPM}/toGEO/E10.5_trunk_WT_ChIPM_CTCF_MAnorm_neq3.bw
height = 2
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 1000
summary_method = max
type= line
line_width = 0.5

[RAD21_WT]
file = ${pathChIPM}/toGEO/E10.5_trunk_WT_ChIPM_RAD21_MAnorm_neq4.bw
height = 2
color = gray 
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 1500
summary_method = max
type= line
line_width = 0.5
orientation = inverted
" > ${ini_file}

pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.pdf}  --region chr2:74650810-74767377 --width 19 --height 2

for prot in CTCF RAD21; do
  ini_file="Fig2_${prot}.ini"
  echo "[scalebar]
file_type = scalebar
height = 0.3
where = top
x_center = 74657000
size = 10000
fontsize = 6
" > ${ini_file}
  for genotype in "WT" "DelCBS1" "DelCBS1-2" "DelCBS1-3" "DelCBS1-4" "DelCBS1-5"; do
    mergedFile=$(ls ${pathChIPM}/toGEO/E10.5_trunk_${genotype}_ChIPM_${prot}_MAnorm_neq*.bw)
    if [ -z "$mergedFile" ]; then
      mergedFile=${pathChIPM}/toGEO/E10.5_trunk_${genotype}_ChIPM_${prot}_MAnorm.bw
    fi
    echo "[${prot}_${genotype}]
file = ${mergedFile}
height = 2
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 1500
summary_method = max

[spacer]
height = 0.15
" >> ${ini_file}
  done
  echo "[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5

[spacer]
height = 0.15

[Hox]
file = Hox.gtf
merge_transcripts = true
color = black
border_color = none
height = 0.4
display = collapsed
labels = False
style = flybase 
color_utr = black

[nonHox]
file = nonHox.gtf
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
file = prot_coding_around_HoxD_start.bed
color = none
border_color = none
color_utr = none
fontsize = 5
height = 0.3
display = collapsed
style = UCSC
line_width = 0
" >> ${ini_file}

  pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.pdf} --region chr2:74650810-74767377 --width 9 --height 6
done


# Figure 3
ini_file="Fig3.ini"
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 74144685
size = 200000

[spacer]
height = 0.1

[CHiC_E9.5_WT_10kb]
file = ${pathCHiC}/toGEO/E9.5_trunk_WT_cHiC_merge_10kb.cool
colormap = ['white', 'darkblue', (.3333, 0, 0)]
rasterize = false
depth = 1200000
min_value = 0
max_value = 0.04
show_masked_bins = false
file_type = hic_matrix

[spacer]
height = 0.1

[tads]
file = ${pathCHiC}/TADs/E9.5_trunk_WT_cHiC_merge_10kb.240kb_domains.bed
display = interleaved
height = 0.2
color= darkgray
labels = False
line_width=0

[spacer]
height = 0.1

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.4

[spacer]
height = 0.05

[genes]
file = nonHox.gtf
merge_transcripts = true
color = #a6a8ab
color_utr = #a6a8ab
border_color = none
height = 0.4
display = collapsed
labels = False
style = flybase

[hoxD]
file = HoxD_CS3840.bed
color = bed_rgb
color_utr = bed_rgb
border_color = none
overlay_previous = yes
labels = False

[hoxD]
file = HoxD_CS3840.bed
color = none
color_utr = none
border_color = none
height = 0.3
" > ${ini_file}
for geno in "1" "1-3" "1-5"; do
  echo "[spacer]
height = 0.2

[Diff_WT_DelCBS${geno}_10kb]
file = ${pathCHiC}/diff/WT_minus_E9.5_trunk_DelCBS${geno}_cHiC_merge_10kb.cool
title = Diff_WT_DelCBS${geno}_10kb
colormap = ['darkred', 'white', 'darkblue']
rasterize = false
depth = 1200000
min_value = -0.00001
max_value = 0.00001
show_masked_bins = false
file_type = hic_matrix
" >> ${ini_file}
done
echo "[spacer]
height = 0.1

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.4

[spacer]
height = 0.1

[genes]
file = nonHox.gtf
merge_transcripts = true
color = #a6a8ab
color_utr = #a6a8ab
border_color = none
height = 0.4
display = collapsed
labels = False
style = flybase

[hoxD]
file = HoxD_CS3840.bed
color = bed_rgb
color_utr = bed_rgb
border_color = none
overlay_previous = yes
labels = False

[hoxD]
file = HoxD_CS3840.bed
color = none
color_utr = none
border_color = none
height = 0.3
" >> ${ini_file}

pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.png} --region chr2:73779626-75669724 --width 19.5 --height 23 --dpi 500

## Figure 4Abottom
ini_file="Fig4Abottom.ini"
echo "[Genes_file]
file = prot_coding.gtf
color = black
border_color = none
height = 1
display = collapsed
labels = False
style = flybase 
color_utr = black

[spacer]
height = 0.5

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 1
" > ${ini_file}

pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.pdf} --region chr2:74665763-74767377 --width 18 --height 0.5 


## Figure 6 and 7
figs=('6' '7')
tissues=('PFL' 'DFL')
for i in 0 1; do
  tissue=${tissues[$i]}
  ini_file="Fig${figs[$i]}Atop.ini"
  echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 74000000
size = 200000

[CB_CHIC_PL]
file = ${pathCHiC}/toGEO/E12.5_${tissue}_WT_cHiC_10kb.cool
colormap = Greys
rasterize = false
depth = 1200000
height = 10
min_value = 0
max_value = 0.02
show_masked_bins = false
file_type = hic_matrix

[spacer]
height = 0.25

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.7

[genes]
file = ${gitHubDirectory}/plots/annotations/simplified_genes_around_HoxD.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.7
" > ${ini_file}

  pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.png} \
    --region chr2:73779626-75669724 --width 20 --height 6 --dpi 500

  ini_file="Fig${figs[$i]}CE.ini"
  echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 74144685
size = 200000
" > ${ini_file}
  for geno in "1" "1-3" "1-5"; do
    echo "[spacer]
height = 0.1

[Diff_WT_DelCBS${geno}_10kb]
file = ${pathCHiC}/diff/WT_minus_E12.5_${tissue}_DelCBS${geno}_cHiC_10kb.cool
title = Diff_WT_DelCBS${geno}_10kb
colormap = ['darkred', 'white', 'darkblue']
rasterize = false
depth = 1200000
min_value = -0.00001
max_value = 0.00001
show_masked_bins = false
file_type = hic_matrix
" >> ${ini_file}
  done
  echo "[spacer]
height = 0.05

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.4

[spacer]
height = 0.05

[genes]
file = nonHox.gtf
merge_transcripts = true
color = #a6a8ab
color_utr = #a6a8ab
border_color = none
height = 0.4
display = collapsed
labels = False
style = flybase

[hoxD]
file = HoxD_CS3840.bed
color = bed_rgb
color_utr = bed_rgb
border_color = none
overlay_previous = yes
labels = False

[hoxD]
file = HoxD_CS3840.bed
color = none
color_utr = none
border_color = none
height = 0.3
" >> ${ini_file}

  pyGenomeTracks --tracks ${ini_file}  -out ${ini_file/.ini/.png} --region chr2:73779626-75669724 --width 18 --height 15 --dpi 500
done


## Figure S1
ini_file="FigS1.ini"
echo "[scalebar]
file_type = scalebar
height = 0.2
where = top
#x_center = 74550000
x_center = 74657000
#scalebar_start_position = 74466370
#scalebar_end_position = 74767370
size = 10000
fontsize = 6

[spacer]
height = 0.25
" > ${ini_file}

for tissue in E12_brain E10_post_trunk E12_PFL E12_DFL; do
  echo "[CTCF ${tissue}]
file = ${pathChIP}/toGEO/${tissue}_ChIP_CTCF.bw
title = ${tissue}
height = 4
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 3
summary_method = max
type= line

[spacer]
height = 0.25
" >> ${ini_file}
done
echo "[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5

[spacer]
height = 0.05

[Hox]
file = Hox.gtf
merge_transcripts = true
color = black
border_color = none
height = 0.5
display = collapsed
labels = False
style = flybase 
color_utr = black

[nonHox]
file = nonHox.gtf
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
height = 0.1

[gene names]
file = prot_coding_around_HoxD_start.bed
color = none
border_color = none
color_utr = none
fontsize = 5
height = 0.3
display = collapsed
style = UCSC
line_width = 0
" >> ${ini_file}
pyGenomeTracks --tracks ${ini_file} -out ${ini_file/.ini/.pdf} --region chr2:74650810-74767377 --width 18 --height 8


## Figure S1 Rita
ini_file="FigS1_Rita.ini"
echo "[scalebar]
file_type = scalebar
height = 0.2
where = top
#x_center = 74550000
x_center = 74657000
#scalebar_start_position = 74466370
#scalebar_end_position = 74767370
size = 10000
fontsize = 6

[spacer]
height = 0.25
" > ${ini_file}

for tissue in E12_brain E10_post_trunk E12_PFL E12_DFL; do
  echo "[CTCF ${tissue}]
file = ${pathChIP}/toGEO/${tissue}_ChIP_CTCF.bw
height = 4
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 3
summary_method = max
type= line

[spacer]
height = 0.25
" >> ${ini_file}
done
echo "[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5

[spacer]
height = 0.05

[Hox]
file = Hox.gtf
merge_transcripts = true
color = black
border_color = none
height = 0.5
display = collapsed
labels = False
style = flybase 
color_utr = black

[nonHox]
file = nonHox.gtf
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
height = 0.4
" >> ${ini_file}
pyGenomeTracks --tracks ${ini_file} -out ${ini_file/.ini/.pdf} --region chr2:74650810-74767377 --width 18 --height 8

## Figure S2
ini_file="FigS2.ini"
echo "[scalebar]
file_type = scalebar
height = 0.2
where = top
x_center = 74658000
size = 10000
fontsize = 6

[scalebar]
file_type = scalebar
height = 0.2
where = top
x_center = 52271000
size = 10000
fontsize = 6
overlay_previous = yes

[scalebar]
file_type = scalebar
height = 0.2
where = top
x_center = 96196000
size = 10000
fontsize = 6
overlay_previous = yes

[scalebar]
file_type = scalebar
height = 0.2
where = top
x_center = 102923000
size = 10000
fontsize = 6
overlay_previous = yes

[CTCF_WT]
file = ${pathChIPM}/toGEO/E10.5_trunk_WT_ChIPM_CTCF_MAnorm_neq3.bw
title = Wt_CTCF
height = 3
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 1110
summary_method = max
type = line

[spacer]
height = 0.08

[CTCF peaks]
file = E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored_cov_pos.bed
display = collapsed
color = ['white', (0.93, 0.11, 0.14)]
min_value = 0
max_value = 500
border_color = none
labels = false
arrowhead_fraction = 0.01
title = all 3 summits +-250bp

[CTCF peaks]
file = E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored_cov_neg.bed
display = collapsed
color = ['white', (0.18, 0.19, 0.57)]
min_value = 0
max_value = 500
border_color = none
labels = false
arrowhead_fraction = 0.01
overlay_previous = yes

[spacer]
height = 0.08

[Genes_file]
file = Hox.gtf
file_type = gtf
prefered_name = gene_name
merge_transcripts = true
color = black
border_color = black
height = 0.3
display = collapsed
labels = false
style = flybase 
color_utr = black

[Genes_file]
file = nonHox.gtf
file_type = gtf
prefered_name = gene_name
merge_transcripts = true
color = grey
border_color = grey
display = collapsed
labels = false
style = flybase 
color_utr = grey
color_backbone = grey
overlay_previous = yes

[spacer]
height = 0.2

[Genes_file]
file = tss.gtf
file_type = gtf
prefered_name = gene_name
merge_transcripts = true
color = none
border_color = none
color_utr = none
fontsize = 5
height = 0.3
display = collapsed
style = UCSC
line_width = 0
" > ${ini_file}

pyGenomeTracks --tracks ${ini_file} --out ${ini_file/.ini/.pdf} --BED Hox_pos.bed --width 19 --height 3

cat ${ini_file} | awk '{if ($0~/pos/){gsub("pos", "neg", $0)}else{gsub("neg", "pos", $0)}; print}' > ${ini_file/.ini/_rev.ini}

pyGenomeTracks --tracks ${ini_file/.ini/_rev.ini} --out ${ini_file/.ini/.pdf} --BED Hox_neg.bed --width 19 --height 3 --decreasingXAxis

pdflatex -jobname FigS2 &> /dev/null <<-EOF
\\documentclass[preview]{standalone}
\\usepackage{graphicx}
\\begin{document}
\\includegraphics{FigS2_chr6-52151994-52278047.pdf}
\\includegraphics{FigS2_chr11-96189088-96376579.pdf}
\\includegraphics{FigS2_chr15-102916367-103042207.pdf}
\\includegraphics{FigS2_chr2-74650810-74767377.pdf}
\\end{document}
EOF

## Figure S4
ini_file="FigS4A.ini"
echo "[gene names]
file = prot_coding_around_HoxD_start.bed
color = none
border_color = none
color_utr = none
fontsize = 5
height = 0.3
display = collapsed
style = UCSC
line_width = 0

[spacer]
height = 0.1

[Genes_file]
file = Hox.gtf
file_type = gtf
prefered_name = gene_name
merge_transcripts = true
color = black
border_color = black
height = 0.4
display = collapsed
labels = false
style = flybase 
color_utr = black

[Genes_file]
file = nonHox.gtf
file_type = gtf
prefered_name = gene_name
merge_transcripts = true
color = grey
border_color = grey
display = collapsed
labels = false
style = flybase 
color_utr = grey
color_backbone = grey
overlay_previous = yes

[spacer]
height = 0.25


[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.8

[spacer]
height = 0.25

[scalebar]
file_type = scalebar
height = 0.2
where = top
x_center = 74657000
size = 10000
fontsize = 6

[CTCF_WT]
file = ${pathChIPM}/toGEO/E10.5_trunk_WT_ChIPM_CTCF_MAnorm_neq3.bw
height = 2
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 1000
summary_method = max
type= line
line_width = 0.5

[Rad21_WT]
file = ${pathChIPM}/toGEO/E10.5_trunk_WT_ChIPM_RAD21_MAnorm_neq4.bw
height = 2
color = gray 
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
#max_value = 10
max_value = 1500
summary_method = max
type= line
line_width = 0.5
orientation = inverted
" > ${ini_file}

pyGenomeTracks --tracks ${ini_file} --out ${ini_file/.ini/.pdf} --region chr2:74650810-74767377 --width 19 --height 2.2

## Figure S6
ini_file="FigS6A.ini"
echo "" > ${ini_file}
for geno in "1" "1-3" "1-5"; do
  echo "[Diff_WT_DelCBS${geno}_10kb_hicTransform]
file = ${pathCHiC}/correct/WT_minus_E9.5_trunk_DelCBS${geno}_cHiC_merge_10kb_hicTransform.cool
title = Diff_WT_DelCBS${geno}_10kb_hicTransform
colormap = ['darkred', 'white', 'darkblue']
rasterize = false
depth = 1200000
min_value = -0.000002
max_value = 0.000002
show_masked_bins = false
file_type = hic_matrix

[quantification]
file = ${gitHubDirectory}/cHi-C/trunk_S6.bedpe
color = black
line_width = 3
overlay_previous = share-y
file_type = links
links_type = loops

[spacer]
height = 0.1
" >> ${ini_file}
done

echo "[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.6

[spacer]
height = 0.1

[Hox]
file = Hox.gtf
merge_transcripts = true
color = black
border_color = none
height = 0.3
display = collapsed
labels = False
style = flybase 
color_utr = black

[nonHox]
file = nonHox.gtf
merge_transcripts = true
color = grey
border_color = none
height = 0.3
display = collapsed
labels = False
style = flybase 
color_utr = grey
color_backbone = grey
overlay_previous = yes

[spacer]
height = 0.2

[gene names]
file = prot_coding_around_HoxD_start.bed
color = none
border_color = none
color_utr = none
fontsize = 5
height = 0.3
display = collapsed
style = UCSC
line_width = 0
" >> ${ini_file}

pyGenomeTracks --tracks ${ini_file} --out ${ini_file/.ini/.png} --region chr2:74650810-74767377 --width 9.5 --height 9.5  --dpi 500

ini_file="FigS6C.ini"
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
size = 10000
line_width=0.5
fontsize = 5

[spacer]
height = 0.1
" > ${ini_file}
genos=('WT' 'DelCBS1' 'DelCBS1-3' 'DelCBS1-5')
colors=('#08306B' '#2171B5' '#6BAED6' '#C6DBEF')
i=0
echo "[insulation_E9.5_${genos[$i]}_10kb]
file = ${pathCHiC}/TADs/E9.5_trunk_${genos[$i]}_cHiC_merge_10kb.240kb_tad_score.bm
type = line
height = 4
color= ${colors[$i]}
file_type=bedgraph
min_value = -1.8
max_value = -1
" >> ${ini_file}
for i in 1 2 3; do
  echo "[insulation_E9.5_${genos[$i]}_10kb]
file = ${pathCHiC}/TADs/E9.5_trunk_${genos[$i]}_cHiC_merge_10kb.240kb_tad_score.bm
type = line
color= ${colors[$i]}
file_type=bedgraph
overlay_previous = share-y
" >> ${ini_file}
done
for i in 0 1 2 3; do
  echo "[spacer]
height = 0.1

[TAD_E9.5_${genos[$i]}_10kb]
file = ${pathCHiC}/TADs/E9.5_trunk_${genos[$i]}_cHiC_merge_10kb.240kb_domains.bed
display = interleaved
height = 0.3
color= ${colors[$i]}
labels = False
line_width=0
" >> ${ini_file}
done

echo "[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.25

[spacer]
height = 0.1

[Hox]
file = Hox.gtf
merge_transcripts = true
color = black
border_color = none
height = 0.25
display = collapsed
labels = False
style = flybase 
color_utr = black

[nonHox]
file = nonHox.gtf
merge_transcripts = true
color = grey
border_color = none
height = 0.3
display = collapsed
labels = False
style = flybase 
color_utr = grey
color_backbone = grey
overlay_previous = yes

[spacer]
height = 0.2

[gene names]
file = prot_coding_around_HoxD_start.bed
color = none
border_color = none
color_utr = none
fontsize = 5
height = 0.25
display = collapsed
style = UCSC
line_width = 0
" >> ${ini_file}

pyGenomeTracks --tracks ${ini_file} --dpi 500 --width 7 --height 5 --region chr2:74580020-74780798 -o ${ini_file/.ini/.png}

## Figure S7
ini_file="FigS7A.ini"
echo "[scalebar]
file_type = scalebar
x_center = 74570000
height = 0.5
where = top
size = 200000
fontsize = 5

[spacer]
height = 0.1
" > ${ini_file}

for genotype in "WT" "DelCBS1" "DelCBS1-3" "DelCBS1-5"; do
  echo "[CHiC_E9.5_${genotype}_10kb]
file = ${pathCHiC}/toGEO/E9.5_trunk_${genotype}_cHiC_merge_10kb.cool
title = CHiC_E9.5_${genotype}_merge_iced_10kb
colormap = ['white', 'darkblue', (.33, 0, 0)]
rasterize = false
depth = 1200000
min_value = 0
max_value = 0.04
show_masked_bins = false
file_type = hic_matrix

[quantification]
file = ${gitHubDirectory}/cHi-C/trunk_S7.bedpe
color = black
line_width = 0.5
overlay_previous = share-y
file_type = links
links_type = loops

[spacer]
height = 0.1
" >> ${ini_file}
done
echo "[spacer]
height = 0.1

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.4

[spacer]
height = 0.1

[genes]
file = nonHox.gtf
merge_transcripts = true
color = #a6a8ab
color_utr = #a6a8ab
border_color = none
height = 0.4
display = collapsed
labels = False
style = flybase

[hoxD]
file = HoxD_CS3840.bed
color = bed_rgb
color_utr = bed_rgb
border_color = none
overlay_previous = yes
labels = False

[hoxD]
file = HoxD_CS3840.bed
color = none
color_utr = none
border_color = none
height = 0.3
fontsize = 5
" >> ${ini_file}

pyGenomeTracks --tracks ${ini_file} --dpi 500 --width 9 --height 12 \
  --region chr2:74469563-75691603 -o ${ini_file/.ini/.png}

ini_file="FigS7C.ini"
viewpoints=('Hoxd4' 'Hoxd4' 'Hoxd9' 'border' 'Hoxd4' 'Hoxd9' 'border')
mutants=('1' '1-3' '1-3' '1-3' '1-5' '1-5' '1-5')
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
size = 200000

[spacer]
height = 0.1
"
for i in "${!viewpoints[@]}"; do
  mutant=DelCBS${mutants[$i]}
  if [ "$prev_mut" = "$mutant" ]; then
    spacer=0.1
  else
    spacer=0.2
  fi
  prev_mut=$mutant
  echo "[spacer]
height = ${spacer}

[${viewpoints[$i]} WT]
file = ${pathCHiC}/vC/E9.5_trunk_WT_cHiC_merge_vC_${viewpoints[$i]}.bedgraph
color = blue
height = 2
title =${viewpoints[$i]} WT vs ${mutant}
min_value = 0
max_value = 4
alpha = 0.9
summary_method = max
number_of_bins = 400

[${viewpoints[$i]} ${mutant}]
file = ${pathCHiC}/vC/E9.5_trunk_${mutant}_cHiC_merge_vC_${viewpoints[$i]}.bedgraph
color = red
alpha = 0.5
summary_method = max
number_of_bins = 400
overlay_previous = share-y
" >> ${ini_file}
done

echo "[spacer]
height = 0.1

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.4

[spacer]
height = 0.1

[genes]
file = nonHox.gtf
merge_transcripts = true
color = #a6a8ab
color_utr = #a6a8ab
border_color = none
height = 0.4
display = collapsed
labels = False
style = flybase

[hoxD]
file = HoxD_CS3840.bed
color = bed_rgb
color_utr = bed_rgb
border_color = none
overlay_previous = yes
labels = False

[hoxD]
file = HoxD_CS3840.bed
color = none
color_utr = none
border_color = none
height = 0.3
fontsize = 5
" >> ${ini_file}

pyGenomeTracks --tracks ${ini_file} --dpi 500 --width 9 --height 8 \
  --region chr2:74469563-75691603 -o ${ini_file/.ini/.png}



## Figure S9
ini_file="FigS9E.ini"
echo "[gene names]
file = prot_coding_around_HoxD_start.bed
color = none
border_color = none
color_utr = none
fontsize = 5
height = 0.3
display = collapsed
style = UCSC
line_width = 0

[spacer]
height = 0.1

[Genes_file]
file = Hox.gtf
file_type = gtf
prefered_name = gene_name
merge_transcripts = true
color = black
border_color = black
height = 0.4
display = collapsed
labels = false
style = flybase 
color_utr = black

[Genes_file]
file = nonHox.gtf
file_type = gtf
prefered_name = gene_name
merge_transcripts = true
color = grey
border_color = grey
display = collapsed
labels = false
style = flybase 
color_utr = grey
color_backbone = grey
overlay_previous = yes

[scalebar]
file_type = scalebar
height = 0.5
where = top
x_center = 74657000
size = 10000
fontsize = 6

[CTCF_WT]
file = ${pathChIPM}/toGEO/E10.5_trunk_WT_ChIPM_CTCF_neq3.bw
height = 2
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 10
summary_method = max
type= line:0.3

[Rad21_WT]
file = ${pathChIPM}/toGEO/E10.5_trunk_WT_ChIPM_RAD21_neq4.bw
color = gray 
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 10
summary_method = max
type= line:0.3
overlay_previous = share-y

[spacer]
height = 0.1

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5

[spacer]
height = 1

[scalebar]
file_type = scalebar
height = 0.3
where = top
x_center = 74657000
size = 10000
fontsize = 6
" > ${ini_file}
for prot in CTCF RAD21; do
  echo "[spacer]
height = 0.15

[CTCF_deld8d9]
file = ${pathChIPM}/toGEO/E10.5_trunk_DelCBS2_ChIPM_${prot}.bw
title = ${prot}
height = 2
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 10
summary_method = max
" >> ${ini_file}
done
echo "[spacer]
height = 0.1

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5

[spacer]
height = 0.15

[Genes_file]
file = Hox.gtf
file_type = gtf
prefered_name = gene_name
merge_transcripts = true
color = black
border_color = black
height = 0.4
display = collapsed
labels = false
style = flybase 
color_utr = black

[Genes_file]
file = nonHox.gtf
file_type = gtf
prefered_name = gene_name
merge_transcripts = true
color = grey
border_color = grey
display = collapsed
labels = false
style = flybase 
color_utr = grey
color_backbone = grey
overlay_previous = yes

[gene names]
file = prot_coding_around_HoxD_start.bed
color = none
border_color = none
color_utr = none
fontsize = 5
height = 0.5
display = collapsed
style = UCSC
line_width = 0
" >> ${ini_file}

pyGenomeTracks --tracks ${ini_file} -out ${ini_file/.ini/.pdf} \
  --region chr2:74650810-74767377 --width 11.7 --height 4

## Figure S13
ini_file=FigS13A.ini
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
size = 200000
" > ${ini_file}

for genotype in "WT" "DelCBS1" "DelCBS1-3" "DelCBS1-5"; do
  echo "[spacer]
height = 0.1

[CHiC_PFL_${genotype}_10kb]
file = ${pathCHiC}/toGEO/E12.5_PFL_${genotype}_cHiC_10kb.cool
title = CHiC_PFL_${genotype}_iced_10kb
colormap = ['white', 'darkblue', (.33, 0, 0)]
rasterize = false
depth = 1200000
min_value = 0
max_value = 0.035
show_masked_bins = false
file_type = hic_matrix

[quantification]
file = ${gitHubDirectory}/cHi-C/PFL_S13.bedpe
color = black
line_width = 0.5
overlay_previous = share-y
file_type = links
links_type = loops
" >> ${ini_file}
done
echo "[spacer]
height = 0.05

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2

[spacer]
height = 0.05

[genes]
file = nonHox.gtf
merge_transcripts = true
color = #a6a8ab
color_utr = #a6a8ab
border_color = none
height = 0.2
display = collapsed
labels = False
style = flybase

[hoxD]
file = HoxD_CS3840.bed
color = bed_rgb
color_utr = bed_rgb
border_color = none
overlay_previous = yes
labels = False

[hoxD]
file = HoxD_CS3840.bed
color = none
color_utr = none
border_color = none
height = 0.2
fontsize = 5
" >> ${ini_file}

pyGenomeTracks --tracks ${ini_file} --dpi 500 --width 9 --height 12 \
  --region chr2:74469563-75691603 -o ${ini_file/.ini/.png}

ini_file="FigS13C.ini"
viewpoints=('Hoxd4' 'Hoxd8' 'Hoxd8' 'Hoxd8')
mutants=('1' '1' '1-3' '1-5')
echo "[scalebar]
file_type = scalebar
height = 0.5
where = top
size = 200000

[spacer]
height = 0.1
" > ${ini_file}
for i in "${!viewpoints[@]}"; do
  mutant=DelCBS${mutants[$i]}
  if [ ! "$i" = "0" ]; then
    echo "[spacer]
height = 0.2
" >> ${ini_file}
  fi
echo "[${viewpoints[$i]} WT]
file = ${pathCHiC}/vC/E12.5_PFL_WT_cHiC_vC_${viewpoints[$i]}.bedgraph
color = blue
height = 2
title =${viewpoints[$i]} WT vs ${mutant}
min_value = 0
max_value = 4
alpha = 0.9
summary_method = max
number_of_bins = 400

[${viewpoints[$i]} ${mutant}]
file = ${pathCHiC}/vC/E12.5_PFL_${mutant}_cHiC_vC_${viewpoints[$i]}.bedgraph
color = red
alpha = 0.5
summary_method = max
number_of_bins = 400
overlay_previous = share-y
" >> ${ini_file}
done

echo "[spacer]
height = 0.05

[CTCF peaks]
file = ${gitHubDirectory}/ChIPmentation/outputs/E10.5_trunk_WT_ChIPM_CTCF_rep1and2and3_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2

[spacer]
height = 0.05

[genes]
file = nonHox.gtf
merge_transcripts = true
color = #a6a8ab
color_utr = #a6a8ab
border_color = none
height = 0.2
display = collapsed
labels = False
style = flybase

[hoxD]
file = HoxD_CS3840.bed
color = bed_rgb
color_utr = bed_rgb
border_color = none
overlay_previous = yes
labels = False

[hoxD]
file = HoxD_CS3840.bed
color = none
color_utr = none
border_color = none
height = 0.3
fontsize = 5
" >> ${ini_file}

pyGenomeTracks --tracks ${ini_file} --dpi 500 --width 9 --height 6 \
  --region chr2:74469563-75691603 -o ${ini_file/.ini/.png}
