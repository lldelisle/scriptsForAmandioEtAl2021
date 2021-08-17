options(stringsAsFactors = F)
rm(list = ls())

# Get arguments
gitHubDirectory <- commandArgs(TRUE)[1]

# Load dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("VennDiagram")
safelyLoadAPackageInCRANorBioconductor("gplots")

summary.df <- read.delim(file.path(gitHubDirectory, "RNAseq", "output", "DESeq2_pairwise", "summary_significant.txt"))
genotypes <- list('PFL' = paste0("Del(CBS1", c("", paste0("-", 2:5)), ")"),
                  'DFL' = paste0("Del(CBS1-", 2:4, ")"))
input <- list()
for (tissue in names(genotypes)){
  input[[tissue]] <- lapply(genotypes[[tissue]], function(geno){summary.df$gene_short_name[!is.na(summary.df[, paste0(tissue, "_", gsub("\\(|\\)|-", ".", geno), "_padj")])]})
  names(input[[tissue]]) <- genotypes[[tissue]]
  tmp <- venn(input[[tissue]])
  print(tissue)
  print(attr(tmp, "intersections")[[paste(genotypes[[tissue]], collapse = ":")]])
}
# [1] "PFL"
# [1] "Hoxd3" "Hoxd4"
# [1] "DFL"
# [1] "Hoxd11" "Hoxd10" "Hoxd9" 

# Plot each venn diagram:
# PFL:
venn.diagram(
  x = input[["PFL"]],
  filename = file.path(gitHubDirectory, "RNAseq", "output", "DESeq2_pairwise", 'PFL_All_paper.png'),
  output=FALSE,
  # Circles
  lwd = 3,
  col= c("#31A354", "#74C476", "#A1D99B", "#C7E9C0", "#EDF8E9"),
  # Numbers
  cex = 1.5,
  fontfamily = "Arial",
  # Set names
  cat.cex = 1.5 ,
  cat.default.pos = "outer",
  cat.pos = c(-70, 1, -80,150, 10),
  cat.dist = c(0.3, 0.2, 0.3, 0.21, 0.2),
  cat.fontfamily = "Arial",
  ##Output features
  imagetype="png" ,
  resolution = 500,
  compression = "lzw")

venn.diagram(
  x = input[["DFL"]],
  filename = file.path(gitHubDirectory, "RNAseq", "output", "DESeq2_pairwise", 'DFL_All_paper.png'),
  # Output features
  imagetype="png" ,
  resolution = 500,
  compression = "lzw",
  # Circles
  lwd = 3,
  col= c("#DF65B0", "#C994C7", "#D4B9DA"),
  # Numbers
  cex = 1.5,
  fontfamily = "Arial",
  #Set names
  cat.cex = 1.5,
  cat.default.pos = "outer",
  cat.pos = c(330, 30, 0),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.fontfamily = "Arial",
  rotation = 1
)

