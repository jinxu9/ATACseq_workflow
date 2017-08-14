library("RColorBrewer")                                                                                                                                                            
library("gplots")
countsTable <- read.table( "ECs_merge_peaks.bed.Brain.Skin.Liver.lung.readscount", stringsAsFactors=TRUE )
rownames( countsTable ) <- countsTable$V4
#peak<-countsTable[,4]
countsTable <- countsTable[ , 8:20]
colData <- data.frame(condition=factor(c("B","B","B","B","noB","noB","noB","noB","noB","noB","noB","noB","noB")),tissue=factor(c("B","B","B","B","S","S","S","Li","Li","Li","Lu","Lu","Lu")),samples=factor(c("Brain1","Brain2","Brain3","Brain4","Skin1","Skin2","Skin3","Liver1","Liver2","Liver3","Lung1","Lung2","Lung3")))
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~tissue)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
write.table(as.data.frame(res.BS),"DESseq2_differentialPeaks.txt", row.name=T,sep="\t",quote=F)
