options(stringsAsFactors = F)
rm(list = ls())

# Get arguments
gitHubDirectory <- commandArgs(TRUE)[1]

path <- commandArgs(TRUE)[2]

# Load dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("DESeq2")
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
source(file.path(gitHubDirectory, "RNAseq/functions.R"))

# Fixed variables:
tableWithCounts <- file.path(path, "mergedTables/AllHTSeqCounts_subset.txt")
gtf.file <- list.files(path, pattern = "merge.*UCSC.gtf", full.names = T)
samplesPlan <- file.path(gitHubDirectory, "RNAseq/samplesPlan.txt")
factorToStudy <- "genotype"
pathForDESeq2 <- file.path(gitHubDirectory, "RNAseq/output/DESeq2_pairwise/")
log2FC.threshold <- log2(1.5)

# Prepare inputs
samplesPlanDF<-read.delim(samplesPlan)
rownames(samplesPlanDF) <- samplesPlanDF$sample
factorizedSP <- samplesPlanDF
for (cn in colnames(factorizedSP)){
  uniqVal <- unique(factorizedSP[,cn])
  factorizedSP[,cn] <- factor(factorizedSP[,cn], levels=uniqVal)
}
samplesToPlot <- samplesPlanDF$sample
count.table <- read.delim(tableWithCounts)
colnames(count.table)[1] <- "gene_id"
colnames(count.table) <- gsub("CBS1\\.", "CBS1-", colnames(count.table))


if (!dir.exists(pathForDESeq2)){
  dir.create(pathForDESeq2, recursive = T)
}
# Prepare a big table with the results of all DESeq2
gtf <- readGFF(gtf.file)
big.annot <- unique(gtf[, c("gene_id", "gene_name", "seqid", "gene_biotype")])
colnames(big.annot) <- c("gene_id", "gene_short_name", "chr", "gene_biotype")
count.table <- merge(count.table, big.annot, all.x = T)
count.table <- subset(count.table, gene_biotype == "protein_coding")
rownames(count.table) <- count.table$gene_id
rm(gtf)

big.annot <- subset(big.annot, gene_biotype == "protein_coding")

# For each tissue
for (tis in unique(samplesPlanDF$tissue)){
  big.annot[, tis] <- ""
  for (geno in setdiff(unique(samplesPlanDF$genotype), "WT")){
    # Select the samples
    new.samples.plan <- subset(samplesPlanDF, tissue == tis & genotype %in% c("WT", geno))
    # Run or read DESeq2 results with Wald test threshold of FC at 1.5
    if ( ! file.exists(paste0(pathForDESeq2,"/",tis,"_", geno, "vsWT_DESeq2significant.txt"))){
      signif <- simpleDeseqAna(count.table, factorToStudy, paste0(pathForDESeq2,"/",tis,"_", geno, "vsWT_"),
                               new.samples.plan,
                               LRT = F, lfcT = log2FC.threshold, writeRLOG = F,
                               theta = c(0.15, 0.99))
    } else {
      signif <- read.delim(paste0(pathForDESeq2,"/",tis,"_", geno, "vsWT_DESeq2significant.txt"))
    }
    # Add results to the dataframe
    big.annot[, paste0(tis, "_", geno, "_l2fc")] <- signif$log2FoldChange[match(big.annot$gene_id, signif$gene_id)]
    big.annot[, paste0(tis, "_", geno, "_padj")] <- signif$padj[match(big.annot$gene_id, signif$gene_id)]
    big.annot[big.annot$gene_id %in% signif$gene_id, tis] <- paste0(big.annot[big.annot$gene_id %in% signif$gene_id, tis], "_", geno)
  }
  big.annot[, tis] <- gsub("^_", "", big.annot[, tis])
}
write.table(big.annot, file.path(pathForDESeq2, "summary_significant.txt"), sep = "\t", quote = F, row.names = F)
