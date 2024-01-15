library(DESeq2) 
library(geneplotter)
library(ggthemes)
library(pheatmap)
library(GO.db) 
library(org.Hs.eg.db) 
library(topGO)
library(GSEABase)
library(clusterProfiler) 
library(biomaRt)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(scatterplot3d)
library(ggrepel)
library(DESeq2)
library(patchwork)
library(cowplotl)
library(ggplotify)
library(tinyarray)
library(stringr)
geneID_symbol=read.csv("C:\\Users\\Administrator\\OneDrive\\课题\\生信\\脚本\\adrenal\\14_bulk\\Macaca_fascicularis_6.0.104.chr.filtered.csv")
data.path="D:/7-Zip/01_adrenal/13_bulk/bulk/MK"
data.list=list.files(data.path,pattern = '*txt',full.names = T)
data.list
setwd("D:/7-Zip/01_adrenal/13_bulk/bulk/MK")

HTseq.handling=function(x,n){
  df=read.table(x,header = F)
  name=substr(x,n,nchar(x)-11)
  colnames(df)=c("ensembl_gene_id",name)
  df=df[-((nrow(df)-4):nrow(df)),]
  return(df)
}
tt=lapply(data.list, HTseq.handling,40)
head(tt[[1]])
tail(tt[[1]])
head(tt[[7]])
data.name=unique(substr(data.list,40,42))
allmerge=function(x){
  merge.data=merge(x[[1]],x[[2]],by='ensembl_gene_id')
  for (i in 3:length(x)) {
    merge.data=merge(merge.data,x[[i]],by='ensembl_gene_id')
  }
  return(merge.data)
}
data.all=allmerge(tt)
write.csv(data.all,paste(data.path='merge_sample.csv',sep = ''))



# DESeq2
dataFrame <-read.csv('merge_sample.csv',row.names = 1)
countMatrix <- as.matrix(dataFrame[2:8])
rownames(countMatrix)<-dataFrame[,1]
condition <- factor(c(rep("OM",4), rep("YM",3)), levels = c("YM","OM"))
sampleNames <- colnames(countMatrix)
colData <- data.frame(Samplename= sampleNames, condition)
rownames(colData) <- sampleNames
dds <- DESeqDataSetFromMatrix(countMatrix, colData, design= ~ condition)
dds <- dds[rowSums(counts(dds)) > 1,]
dds2 <- DESeq(dds)
res <- results(dds2)
rld <- rlog(dds2)
res <- res[order(res$pvalue),]
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
names(resdata)[1]<-"GeneID"
head(resdata)

DEG<-merge(geneID_symbol,resdata,by="GeneID")
rows <- rownames(unique(DEG['gene_name']))
DEG <- DEG[rows,]
DEG = na.omit(DEG)
DEG = subset(DEG,str_sub(gene_name,1,5) != "ENSMF")
head(DEG)
resdata <- DEG
dim(DEG)
# result
setwd('D:/7-Zip/01_adrenal/13_bulk/bulk/MK')

diff_gene <-subset(resdata, padj < 0.05 & (log2FoldChange >= 0.58 | log2FoldChange <= -0.58))
up_gene<-subset(diff_gene,diff_gene$log2FoldChange >= 0.58)
down_gene<-subset(diff_gene,diff_gene$log2FoldChange <= -0.58)
write.csv(resdata, paste("Adrenal_all_gene_anno.csv",sep=""),row.names = F)
write.csv(diff_gene, paste("padj_Adrenal_diff_gene_anno.csv",sep=""),row.names = F)
write.csv(up_gene,paste("padj_Adrenal_up_gene_anno.csv",sep=""),row.names = F)
write.csv(down_gene,paste("padj_Adrenal_down_gene_anno.csv",sep=""),row.names = F)




diff_gene <-subset(resdata, pvalue < 0.05 & (log2FoldChange >= 0.58 | log2FoldChange <= -0.58))
up_gene<-subset(diff_gene,diff_gene$log2FoldChange >= 0.58)
down_gene<-subset(diff_gene,diff_gene$log2FoldChange <= -0.58)
write.csv(resdata, paste("Adrenal_all_gene_anno.csv",sep=""),row.names = F)
write.csv(diff_gene, paste("Adrenal_diff_gene_anno.csv",sep=""),row.names = F)
write.csv(up_gene,paste("Adrenal_up_gene_anno.csv",sep=""),row.names = F)
write.csv(down_gene,paste("Adrenal_down_gene_anno.csv",sep=""),row.names = F)

# figure
# Euclidean distance plot
sampleDists <- dist(t(assay(rld))) #样品距离,欧氏距离,t转置
sampleDistMatrix <- as.matrix(sampleDists)#样品间距离的矩阵
sampleDist_colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sampleDistMatrix
pdf(paste("Adrenal_Euclidean_distances_heatmap.pdf",sep=""),width=8, height=8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=sampleDist_colors)
dev.off()

# sample_PCA
pca_data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar")) # why?
pdf(paste("Adrenal_PCA_plot.pdf",sep=""),width=8, height=8)
ggplot(pca_data, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+theme_bw()
dev.off()


# volcano plot
Volcano_data <-DEG
Volcano_data$threshold <- as.factor(ifelse(Volcano_data$padj < 0.05 & abs(Volcano_data$log2FoldChange) >=0.58,ifelse(Volcano_data$log2FoldChange >= 0.58 ,'Up regulated','Down regulated'),'Unchanged')) 
volcano_sum<-as.data.frame(t(summary(Volcano_data$threshold)))
Volcano_tile <- paste0(paste("Adrenal","_volcano plot",sep=""),'\nCutoff for padj is 0.05','\nThe number of up regulated gene is ',volcano_sum$`Up regulated`,'\nThe number of down regulated gene is ',volcano_sum$`Down regulated`,'\nThe number of unchanged gene is ',volcano_sum$Unchanged )
Volcano_data$lg10<--log10(Volcano_data$padj)
Volcano_data <- subset(Volcano_data,threshold != 'NA')


pdf(paste("Adrenal_volcano plot.pdf",sep=""), width = 6, height = 5.5)
(ggplot( data=Volcano_data,
         aes(x=log2FoldChange, y =lg10,colour=threshold)) +
    scale_color_manual(values=c("deepskyblue4", "grey80","firebrick4"))+
    geom_point(alpha=0.5, size=1.8) +
    # xlim(c(-5,5)) +ylim(c(0,50))+
    ggtitle(Volcano_tile) +
    theme_bw() + 
    geom_vline(xintercept=c(-0.58,0.58),lty=4,col="grey",lwd=0.6)+
    geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
    theme(plot.title = element_text(hjust = 0.5,size=12),legend.title = element_blank())+
    labs(x="log2 (fold change)",y="-log10 (padj)")+
    geom_text_repel(
      data = Volcano_data[-log10(Volcano_data$padj)>=4&abs(Volcano_data$log2FoldChange)>=0.58,],
      aes(x = log2FoldChange, y = lg10, label = gene_name),
      size = 4.5,box.padding = unit(0.5, 'lines'),point.padding = unit(0.8, 'lines'), 
      segment.color = 'grey', show.legend = F,max.overlaps = getOption("ggrepel.max.overlaps", default = 600000)))
dev.off()

