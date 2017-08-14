library("DESeq2")                                                                                                                                                                  
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
pdf("BBB_EC_DESeq_BvsS.pdf")

rld <-rlogTransformation(dds)
#write.table(as.data.frame(rld),"Normalized_log2_matrixtxt", row.name=T,sep="\t",quote=F)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
distsRL <- dist(t(assay(rld))) 
mat <- as.matrix(distsRL)
#samples<-as.factor(c("Brain1","Brain2","Brain3","Brain4","Skin1","Skin2","Skin3","Liver1","Liver2","Liver3","Lung1","Lung2","Lung2"))
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, samples, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(15, 15))
print(plotPCA(rld, intgroup=c("condition","samples")))

dev.off() 


