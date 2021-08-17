#!/usr/bash
inputBedGraph=$1 # The output of macs2 for example ${pathResults}${sample}_macs_likeATAC_treat_pileup.bdg
outputBasename=$2 # The outputBasename of both bedgraph and bw, for example ${path}bedGraphs/${sample}_macs_likeATAC
nameOfTheTrack=$3 # For example "macs2 likeATAC of $sample"
pathForChromSizes=$4 # For example $pathForFasta/${genome}.fa.fai

# First put the chromosome sizes in a string
chromSizes=`cat "$pathForChromSizes" | tr "\t" ":" | tr "\n" ","`
# There are 2 issues with converting bedgraph from macs2 to bw:
# 1) The chromosome names are not sorted
# 2) Sometimes the coverage goes over the size of the chromosome
# Quicker than bedtools sort is to output each chromosome in a different file and then do a cat.
myTmp=`mktemp`
cat "$inputBedGraph" | awk -v h="${nameOfTheTrack}" -v tmp=$myTmp -v s=$chromSizes 'BEGIN{
  print "track type=bedGraph name=\""h"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"
  OFS="\t"
  split(s,mychrs,",");
  for (i in mychrs){
    split(mychrs[i],nameAndSize,":");
    chromSizes[nameAndSize[1]] = nameAndSize[2]
  }
}
$4!=0 && $1 != "track"{
  if ($2<chromSizes[$1]){
    if($3>chromSizes[$1]){
      $3=chromSizes[$1]
    }
    print $0 > tmp"_"$1
  }
}' > "${outputBasename}.bedGraph"
cat ${myTmp}_* >> "${outputBasename}.bedGraph"
bedGraphToBigWig "${outputBasename}.bedGraph" "${pathForChromSizes}" "${outputBasename}.bw"
if [ -e "${outputBasename}.bedGraph.gz" ]; then
  echo "${outputBasename}.bedGraph.gz already exists and will be removed."
  rm "${outputBasename}.bedGraph.gz"
fi
gzip "${outputBasename}.bedGraph"
