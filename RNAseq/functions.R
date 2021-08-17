# Wrapper for DESeq2
simpleDeseqAna <- function(count.table, factorForAna, pathOutput,
                           samplesPlan, LRT=F, pvalT=0.05,
                           lfcT=1.5, writeRLOG=F, ...){
  # Checking the conditions
  if(!(factorForAna%in%colnames(samplesPlan))){
    stop("This is not part of the column names.")
  }
  if(length(levels(samplesPlan[, factorForAna])) == 1){
    stop("The factor you chose have only 1 value. The analysis is not possible.")
  }
  if(length(levels(samplesPlan[, factorForAna])) > 2 & ! LRT){
    print("The factor you chose have more than 2 values. LRT will be applied.")
    LRT <- T
  }
  # Lauching DESeq2
  cmd<-paste("dds <- DESeqDataSetFromMatrix(countData = count.table[, match(rownames(samplesPlan), colnames(count.table))],
             colData = samplesPlan,
             design = ~ ",factorForAna,")",sep='')
  eval(parse(text = cmd))
  print("Genes that are never expressed are removed")
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  if(LRT){
    dds <- DESeq(dds, minReplicatesForReplace=Inf, test="LRT", reduced= ~ 1)
  } else{
    dds <- DESeq(dds, minReplicatesForReplace=Inf)
  }
  res<-results(dds, ...)
  resOrdered<-res[order(res$padj),]
  # Subsetting the annotation file
  ann <- count.table[, colnames(count.table) %in% c("gene_id", "gene_short_name", "locus")]
  annot.df <- data.frame(ann[match(rownames(resOrdered), ann$gene_id), ])
  resToExport <- data.frame(annot.df, counts(dds, normalized = TRUE)[rownames(resOrdered), ], resOrdered)
  write.table(resToExport, file=paste(pathOutput, "DESeq2Results.txt", sep=''), sep='\t', row.names = F, quote=F)
  dfDiffExp <- subset(resToExport, resToExport$padj < pvalT & abs(resToExport$log2FoldChange) > lfcT)
  write.table(dfDiffExp, file=paste0(pathOutput, "DESeq2significant.txt"), sep='\t', row.names = F, quote=F)
  rld <- rlog(dds)
  rlogdata <- assay(rld)
  if(writeRLOG){
    resToExport2 <- data.frame(annot.df, rlogdata[rownames(resOrdered), ])
    write.table(resToExport2, file=paste0(pathOutput, "rlog.txt"), sep='\t', row.names = F, quote=F)
  }
  return(invisible(dfDiffExp))
}

##From Anouk
# Normalize FPKM values using the stable rank genes
normalizeHKRank<-function(FPKMCuff,nbgenes,chrToRm=c("chrM")){
  expdata<-FPKMCuff[,grep("FPKM_",colnames(FPKMCuff))]
  haveRowNames<-F
  if(all(c("gene_id","locus")%in%colnames(FPKMCuff))){
    rownames(expdata)<-paste(FPKMCuff$gene_id,FPKMCuff$locus,sep='__')
    haveRowNames<-T
  } else {
    rownames(expdata)<-paste(1:nrow(FPKMCuff),rep("NA",nrow(FPKMCuff)),sep='__')
  }
  
  if(!all(is.na(chrToRm))){
    #I check that the rownames have info:
    if(!haveRowNames){
      cat("The cufflinks table did not have the gene_id and locus info. It is impossible to remove the specified chr.\n")
    } else {
      #I remove the chrToRm genes
      linesWithChrToRm<-unlist(sapply(chrToRm,function(c){grep(paste0("__",c,":"),rownames(expdata))}))
      if(length(linesWithChrToRm)>0){
        cat(length(linesWithChrToRm),"lines removed.\n")
        expdata<-expdata[-linesWithChrToRm,]
      } else{
        cat("There is no gene in the chromosomes",chrToRm,"\n")
      }    
    }
  }
  
  ### get genes that are expressed in all samples
  expdata.nonzero=expdata[which(apply(expdata,1,min)>0 ),]
  
  ## transform the expression levels into ranks
  
  expdata.ranks=apply(expdata.nonzero,2,rank)
  
  ## compute the variance of the ranks for each gene
  
  expdata.ranks.var=apply(expdata.ranks,1,var,na.rm=T)
  
  ## rank the genes according to their variance
  
  expdata.nonzero.consrank=rank(expdata.ranks.var)
  
  ## compute the median rank over all samples, for each gene
  
  median.rank=apply(expdata.ranks,1,median,na.rm=T)
  
  ## get the genes which have a median rank in the 25%-75% range
  
  interquartile=median.rank>(0.25*length(median.rank)) & median.rank<(0.75*length(median.rank))
  
  ## get the house-keeping genes: the nbgenes genes with the most conserved ranks, among those that fall in the interquartile range
  
  hkgenes=names(sort(expdata.nonzero.consrank[interquartile])[1:nbgenes])
  
  ## compute the normalization coefficient for each sample =  the median of the RPKM of the hkgenes in that sample
  
  normcoeff=apply(expdata.nonzero[hkgenes,],2,median,na.rm=T)
  
  ## we want to bring all of the medians at the average of the medians (computed on all expressed genes)
  
  ##  normcoeff=normcoeff/mean(apply(expdata.nonzero,2,median))
  
  normcoeff=normcoeff/mean(normcoeff)
  
  ## finally, normalize the data
  
  hkgenesEnsID<-unlist(strsplit(hkgenes,"__"))[seq(1,2*length(hkgenes),2)]
  
  if(haveRowNames && "gene_short_name"%in%colnames(FPKMCuff)){
    hkgenesName<-FPKMCuff$gene_short_name[match(hkgenesEnsID,FPKMCuff$gene_id)]
  } else {
    hkgenesName<-NA
  }
  
  FPKMCuff.norm<-data.frame(FPKMCuff[,grep("FPKM_",colnames(FPKMCuff),invert=T)],t(t(FPKMCuff[,grep("FPKM_",colnames(FPKMCuff))])/normcoeff))
  
  results=list("hkgenesName"=hkgenesName,"hkgenesENSID"=hkgenesEnsID,"normcoeff"=normcoeff,"normData"=FPKMCuff.norm)
  
}


# From alko989
# On https://stackoverflow.com/questions/18509527/first-letter-to-upper-case
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
