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
safelyLoadAPackageInCRANorBioconductor("reshape")

quantif.files <- list.files(path = file.path(gitHubDirectory, "cHi-C"), pattern = "quantif.txt", full.names = T)
names(quantif.files) <- sapply(strsplit(basename(quantif.files), "_"), "[[", 2)
quantif.tables <- lapply(quantif.files, read.delim)
quantif.tables.rel <- lapply(names(quantif.tables), function(fig){
  df <- quantif.tables[[fig]]
  new.df <- df[, c(1, grep("relativeToWT", colnames(df)))]
  wt.name <- paste(c(strsplit(colnames(new.df)[2], "_")[[1]][1:2], "WT"), collapse = "_")
  new.df[, wt.name] <- 1
  new.df$fig <- fig
  return(new.df)
})
quantif.tables.melt <- do.call(rbind, lapply(quantif.tables.rel, melt))
p.vals <- lapply(names(quantif.tables), function(fig){
  df <- quantif.tables[[fig]]
  new.df <- df[, c(1, grep("pval", colnames(df)))]
  new.df$fig <- fig
  return(new.df)
})
p.vals.melt <- do.call(rbind, lapply(p.vals, melt))
p.vals.melt$variable <- paste0(gsub("pval_Mann_Whitney_U_test_", "", p.vals.melt$variable), "relativeToWT")
colnames(p.vals.melt)[4] <- "pval"
quantif.tables.melt <- merge(quantif.tables.melt, p.vals.melt, all.x = T)
quantif.tables.melt$pval[is.na(quantif.tables.melt$pval)] <- 1

my.samples <- as.character(unique(quantif.tables.melt$variable))
my.labeller <- gsub("DelC", "Del(C", gsub("\\.", "-", sapply(strsplit(my.samples, "_"), "[[", 3)))
my.labeller[my.labeller != "WT"] <- paste0(my.labeller[my.labeller != "WT"], ")")
quantif.tables.melt$geno <- factor(my.labeller[quantif.tables.melt$variable],
                                   levels = c("WT", "Del(CBS1)", "Del(CBS1-3)", "Del(CBS1-5)"))
quantif.tables.melt$variable <- factor(quantif.tables.melt$variable,
                                       levels = unique(quantif.tables.melt$variable[order(quantif.tables.melt$geno)]))
my.colors <- c("#A1D99B", "#EDF8E9", "#31A354", "#006D2C",
               '#6BAED6', '#C6DBEF', '#2171B5', '#08306B',
               '#6BAED6', '#C6DBEF', '#2171B5')
names(my.colors) <- c("E12.5_PFL_DelCBS1.3_cHiC_10kbrelativeToWT", "E12.5_PFL_DelCBS1.5_cHiC_10kbrelativeToWT", "E12.5_PFL_DelCBS1_cHiC_10kbrelativeToWT", "E12.5_PFL_WT",
                      "E9.5_trunk_DelCBS1.3_cHiC_merge_10kb_hicTransformrelativeToWT", "E9.5_trunk_DelCBS1.5_cHiC_merge_10kb_hicTransformrelativeToWT", "E9.5_trunk_DelCBS1_cHiC_merge_10kb_hicTransformrelativeToWT", "E9.5_trunk_WT",
                      "E9.5_trunk_DelCBS1.3_cHiC_merge_10kbrelativeToWT", "E9.5_trunk_DelCBS1.5_cHiC_merge_10kbrelativeToWT", "E9.5_trunk_DelCBS1_cHiC_merge_10kbrelativeToWT")
region.labeller <- c("Hoxd11-d8\nT-DOMa", "Hoxd11-d1\nCS38-40", "Hoxd11-d1\nT-DOMb",
                     "Hoxd11-d1\nborder", "HoxD", "Hoxd1-Hoxd8", "HoxD / CS38-40", "HoxD / T-DOMb")
names(region.labeller) <- c("d11-d8_T-DOMa","d11-d1_CS39extended", "d11-d1_T-DOMb",
                            "d11-d1_borderCTCFextended", "HoxDEvx_HoxDEvx", "d1-d8_d1-d8", "HoxDEvx_CS39extended", "HoxDEvx_T-DOMb")
quantif.tables.melt$Region <- factor(quantif.tables.melt$Region, levels = names(region.labeller))
for (cur.fig in unique(quantif.tables.melt$fig)){
  cur.df <- subset(quantif.tables.melt, fig == cur.fig)
  cur.colors <- my.colors[as.character(unique(cur.df$variable))]
  names(cur.colors) <- unique(cur.df$geno)
  ggplot(cur.df, aes(x = geno, y = value)) +
    geom_col(aes(fill = geno)) +
    geom_text(data = subset(cur.df, pval < 0.05), label = "*") +
    facet_grid(. ~ Region, switch = "x", labeller = labeller(Region = region.labeller)) +
    scale_fill_manual("",
                      values = cur.colors) + # Use the colors defined on top
    xlab("") +
    ylab("Average signal ratio") +
    theme_minimal() + # Use white background
    theme(legend.position = "top",
          legend.text.align = 0,
          axis.text.x = element_blank(), # Remove the genotype name below each bar
          axis.ticks.x = element_blank(), # Remove the tick for each bar
          panel.grid.major.x = element_blank(), # Remove the vertical grey bar
          strip.text = element_text(size = 11 ), # size of the regions
          )
  ggsave(file.path(gitHubDirectory, "cHi-C", paste0(cur.fig, ".png")),
         width = 7, height = 5)
}
