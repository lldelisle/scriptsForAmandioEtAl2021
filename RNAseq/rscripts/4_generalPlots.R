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
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("pheatmap")
safelyLoadAPackageInCRANorBioconductor("RColorBrewer")
safelyLoadAPackageInCRANorBioconductor("reshape")
safelyLoadAPackageInCRANorBioconductor("rtracklayer")

# Fixed variables:
tableWithNormalizedExpression <- file.path(path, "mergedTables/AllCufflinks_Simplified.txt")
samplesPlan <- file.path(gitHubDirectory, "RNAseq/samplesPlan.txt")
outputFolder <- file.path(gitHubDirectory, "RNAseq/output/plots_FPKM")
fixedColors<-list(genotype=c('WT'="#A6D854",'Del(CBS1)'="#E78AC3",'Del(CBS1-2)'="#8DA0CB",'Del(CBS1-3)'="#FC8D62",'Del(CBS1-4)'="#66C2A5",'Del(CBS1-5)'="#FFD92F"), 
                  tissue=c('DFL'="#9E1F63",'PFL'="#006838"),
                  replicate=c('rep1'="#C6DBEF",'rep2'="#9ECAE1",'rep3'="#6BAED6"), 
                  litter= c("1"="#F0F0F0","2"="#D9D9D9","3"="#BDBDBD","4"="#969696","5"="#737373","6"="#525252"))
fixedColormaps <- list(PFL="Greens",
                       DFL="PuRd")
chrom.to.exclude <- c("chrM", "chrX", "chrY")
epsilon <- 0.0000001
# GTF file
gtf.file <- list.files(path, pattern = "merge.*UCSC.gtf", full.names = T)

if (!dir.exists(outputFolder)){
  dir.create(outputFolder, recursive = T)
}

# Prepare inputs
samplesPlanDF<-read.delim(samplesPlan)
rownames(samplesPlanDF) <- samplesPlanDF$sample
factorizedSP <- samplesPlanDF
for (cn in colnames(factorizedSP)){
  uniqVal <- unique(factorizedSP[,cn])
  factorizedSP[,cn] <- factor(factorizedSP[,cn], levels=uniqVal)
}
expressionDF <- read.delim(tableWithNormalizedExpression)
colnames(expressionDF) <- gsub("CBS1\\.", "CBS1-", colnames(expressionDF))
chrom.expressionDF <- sapply(strsplit(expressionDF$locus, ":"), "[", 1)
expressionDF <- subset(expressionDF, ! chrom.expressionDF %in% chrom.to.exclude)
metaCols <- which(sapply(colnames(expressionDF),function(cn){class(expressionDF[,cn]) != "numeric"}))
samplesToPlot <- intersect(colnames(expressionDF), samplesPlanDF$sample)
expressionDF[, samplesToPlot] <- list(NULL)
colnames(expressionDF) <- gsub("^FPKM_", "", colnames(expressionDF))
gtf <- readGFF(gtf.file)
prot.cod.genes <- gtf$gene_id[gtf$gene_biotype == "protein_coding"]
expressionDF <- subset(expressionDF, gene_id %in% prot.cod.genes)
# Results from DESeq2
summary.df <- read.delim(file.path(gitHubDirectory, "RNAseq/output/DESeq2_pairwise/summary_significant.txt"))
p.vals <- melt(summary.df[grep("Hoxd", summary.df$gene_short_name), c("gene_short_name", grep("padj", colnames(summary.df), value = T))])
p.vals$gene <- p.vals$gene_short_name

# First PCA and correlation for all samples:
samplesToPlot <- intersect(colnames(expressionDF), samplesPlanDF$sample)
data <- expressionDF[,samplesToPlot]
sumperline <- apply(data,1,sum)
nonZdata <- data[sumperline != 0,]
ldata <- log2(nonZdata + 1)
# Select only the top 500 variant genes
rldata <- ldata[order(apply(ldata,1,var),decreasing = T)[1:min(nrow(ldata), 500 )],]
# PCA
sample.pca <- prcomp(t(rldata),
                     center = TRUE,
                     scale. = FALSE)
new.df <- data.frame(factorizedSP[samplesToPlot, ], sample.pca$x[samplesToPlot,])
var <- round((sample.pca$sdev) ^ 2 / sum(sample.pca$sdev ^ 2) * 100)
new.df$genolit <- factor(paste0(new.df$genotype, " & ", new.df$litter), levels = unique(paste0(new.df$genotype, " & ", new.df$litter)))
temp.df <- unique(new.df[, c("genolit", "genotype", "litter")])
shapes.genolit <- c(21:25,24)[temp.df$litter]
names(shapes.genolit) <- temp.df$genolit
fill.genolit <- fixedColors[["genotype"]][temp.df$genotype]
names(fill.genolit) <- temp.df$genolit
i <- 1
j <- 2
png(paste0(outputFolder,"/PC1-PC2_all.png"),width= 20 ,height= 20, res = 200, units = "cm" )
g <- ggplot(new.df, aes_string(paste0("PC", i), paste0("PC", j))) +
  geom_point( aes(fill=genolit, color=tissue, shape=genolit), stroke=2 ,size=3) +
  theme_grey(base_size = 20) +
  xlab(paste0("PC", i, ": ",var[i],"% variance")) +
  ylab(paste0("PC", j, ": ",var[j],"% variance")) +
  scale_fill_manual(values=fill.genolit)+
  scale_color_manual(values=fixedColors[["tissue"]]) +
  scale_shape_manual(values=shapes.genolit) +
  guides(color=guide_legend(override.aes = list(shape = 21))) +
  labs(fill = "genotype & litter", shape = "genotype & litter")
print(g)
dev.off()

# Clustering
sampleDists <- dist(t(rldata), method="euclidean")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- samplesToPlot
colnames(sampleDistMatrix) <- samplesToPlot
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
cols<-c( 2,3,5 ) # to set the variables to annotate 
annot <- NULL
annot<-factorizedSP[,cols]

png(paste0(outputFolder,"/CorrelationMatrix_EuclWard_all.png"), width= 21 ,height= 20, res = 100, units = "cm", type = "cairo")
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  cellwidth = 10,
  cellheight = 10,
  annotation = annot,
  annotation_colors = fixedColors,
  main="Euclidean distance - ward clustering",
  clustering_method="ward.D2",
  col = colors,
  show_colnames = FALSE,
  show_rownames = FALSE,
  fontsize = 15
)
dev.off()

for (tis in c("DFL", "PFL")){
  samplesToPlot <- intersect(colnames(expressionDF), samplesPlanDF$sample[samplesPlanDF$tissue == tis])
  data <- expressionDF[,samplesToPlot]
  ldata <- log2(data + 1)
  dfGene <- data.frame(gene_short_name = intersect(paste0("Hoxd", 1:13), expressionDF$gene_short_name))
  colOfGeneID <- colnames(dfGene)[1]
  # Make a summary plot with all HoxD genes
  subdf <- ldata[match(dfGene[, colOfGeneID], expressionDF[,colOfGeneID]), ]
  subdf$gene <- factor(dfGene[, colOfGeneID], levels = dfGene[, colOfGeneID])
  subdf.gg <- melt(subdf, variable_name = "sample")
  subdf.gg <- merge(subdf.gg, factorizedSP)
  # Get the p-values
  p.vals.current <- aggregate(list(max.value = subdf.gg$value), by = subdf.gg[, c("genotype", "tissue", "gene")], FUN = max)
  p.vals.current$variable <- paste0(p.vals.current$tissue, "_", gsub("\\(|\\)|-", ".", p.vals.current$genotype), "_padj")
  p.vals.current <- merge(p.vals, p.vals.current)
  p.vals.current$gene <- factor(p.vals.current$gene, levels = dfGene[, colOfGeneID])
  p.vals.current$value[is.na(p.vals.current$value)] <- 1
  p.vals.current$p.sign <- symnum(p.vals.current$value, cutpoints = c(0, 0.05, 1), 
                                  symbols = c("*",""))
  scale.pval <- 1.06
  ggplot(subdf.gg, aes(x = genotype, y = value, color = genotype)) +
    stat_summary(aes(fill = genotype), fun = "mean", geom = "bar") + # Plot a bar for the average value and color
    geom_point(shape = 1, color = "gray70") + # select shape of the points and color
    geom_text(data = p.vals.current,
              aes(x = genotype, y = max.value + 0.5, label = as.character(p.sign)),
              color = "black") + # Put the sign for p-value
    facet_grid(. ~ gene, # Split the panel by gene 
               switch = "x") + # Put the label below
    ylab("log2(1 + FPKM)") + 
    xlab("") +
    scale_fill_manual("",
                      values = rev(brewer.pal(n = 6, name = fixedColormaps[[tis]]))) + # Use the colors defined on top
    scale_color_manual("",
                       values = rev(brewer.pal(n = 6, name = fixedColormaps[[tis]]))) + # Use the colors defined on top
    theme_minimal() + # Use white background
    theme(legend.position = "top",
          legend.text.align = 0,
          legend.title = element_text(size = 14), #set size of the title of the legend 
          legend.text = element_text (size = 11), #set size of the legend 
          axis.text.x = element_blank(), # Remove the genotype name below each bar
          axis.ticks.x = element_blank(), # Remove the tick for each bar
          panel.grid.major.x = element_blank(), # Remove the vertical grey bar
          strip.text = element_text(face = "italic", # Put all names of panels in italic
                                    size = 11 ), # size of the genes 
          plot.title = element_text(hjust = 0.5, # Put the title centered
                                    size = 10)) + # Size of title
    scale_y_continuous(expand = c(0, 0), # Do not put extra around the minimum and maximum value
                       limits = c(0, 1.05 * max(subdf[, 1:(ncol(subdf) - 1)]) + 0.5)) + # Set the limit to 0 and 105% of max + 0.5
    ggtitle(paste(tis, "HoxD")) # Title
  ggsave(paste0(outputFolder,"/HoxD_", tis, ".png"), units = "cm",width= 15 ,height= 10, dpi=200, device = "png")
}
