options(stringsAsFactors = F)
# Dependencies:
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("reshape")

# Get arguments
gitHubDirectory <- commandArgs(TRUE)[1]

# First get the quantification result:
quantif <- read.delim("results.txt", quote = "\'")

# I process the colnames of the data frame:
colnames(quantif) <- gsub("^X\\.", "", colnames(quantif))
colnames(quantif) <- gsub("_MAnorm.bw$", "", colnames(quantif))
# I extract from the samples name the info
meta.data <- as.data.frame(matrix(sapply(strsplit(colnames(quantif)[4:ncol(quantif)], "_"), function(v){v[1:6]}),
                           ncol = 6, byrow = T))
colnames(meta.data) <- c("stage", "tissue", "genotype", "technique", "prot", "rep")
# I add a column with the samples name:
meta.data$name <- colnames(quantif)[4:ncol(quantif)]

# I change the structure of the quantif to be compatible with ggplot (1 row per observation)
quantif.m <- melt(quantif, id.vars = colnames(quantif)[1:3], variable_name = "name")

# I add the metadata
quantif.m <- merge(quantif.m, meta.data)

# Reshape the genotype:
quantif.m$genotype[grepl("Del", quantif.m$genotype)] <- paste0(gsub("Del", "Del(", gsub("\\.", "-", quantif.m$genotype[grepl("Del", quantif.m$genotype)])), ")")

# reorder the samples to plot
quantif.m$genotype <- factor(quantif.m$genotype, levels = c("WT", "Del(CBS1)", "Del(CBS1-2)",
                                                            "Del(CBS1-3)", "Del(CBS1-4)", "Del(CBS1-5)")) 


# I plot only the peaks between Evx2 and Hoxd1
my.labeller <- paste0("CBS", 8:1)
names(my.labeller) <- sort(unique(subset(quantif.m, start >= 74655616 & start < 74765142)$start))
ggplot(subset(quantif.m, start >= 74655616 & start < 74765142), aes(x = genotype, y = value)) +
  geom_point(aes(shape = rep), color = "black") +
  facet_grid( prot ~ start, labeller = labeller(start = my.labeller)) +
  theme_light() +
  expand_limits(y = 0) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, size = 15),
        strip.text = element_text(size = 15, color = "black"), # change size of the names of panels 
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +  # Size of ylab
  scale_shape_discrete("replicates")
ggsave(file.path(gitHubDirectory, "ChIPmentation", "outputs", "RAD21_quantification_HoxD_Cluster.pdf"),
       width = 3 + 1.5 * length(my.labeller)) 

# write table for paper with values plotted on the graph
quantif.m.plotted <- subset(quantif.m, start >= 74655616 & start < 74765142)
quantif.m.plotted$name <- NULL
quantif.m.plotted$technique <- NULL
quantif.m.plotted <- cbind(data.frame(name = my.labeller[as.character(quantif.m.plotted$start)]), quantif.m.plotted)
write.table(quantif.m.plotted, file.path(gitHubDirectory, "ChIPmentation", "outputs", "RAD21_quantification_HoxD_Cluster.txt"), sep = "\t", quote = F, row.names = F)  
