

########Gene ste score######
library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggpubr)
library(ggplot2)
library(ggridges)
library(reshape2)
library(reshape)
library(pheatmap)
library(ggplot2)
library(introdataviz)
adrenal <- readRDS('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_adrenal_celltype_final.rds')
gene_list=read.csv("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/10_gene_score/CREM_target0.35.csv")
DefaultAssay(adrenal) <- "RNA"
adrenal <- NormalizeData(object =adrenal, normalization.method = "LogNormalize")
adrenal$group <- adrenal$age
object=colnames(gene_list)
CT.col <- c('#379F7C','#F3A15E')
celltypes <- c("ZG","Transitional cell","ZF", "ZR","Chromaffin cell","Neuron", "FB", "SMC", "EC","Adi", "TC", "Mac")
setwd("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/10_gene_score/")
for (i in object ){ 
  genes <- select(gene_list,i)
  genes <- genes[!is.na(genes)]
  genes <- list(as.character(genes))
  adrenal <- AddModuleScore(adrenal, features=genes, seed=15, name="score")
  data<- FetchData(adrenal,vars = c("celltype","group","sample","score1"))
  data$celltype <- factor(data$celltype, 
                          levels= celltypes , ordered=TRUE)
  data$group <- factor(data$group, 
                       levels=c("Young","Old"), ordered=TRUE)
  data$sample <- factor(data$sample, 
                        levels=c('YM1','YM2','YM3','OM1','OM2','OM3','OM4'), ordered=TRUE)
  pdf(paste0(i,".pdf"), width = 10, height =5)
#####box plot###
  p <- ggplot(data,aes(celltype,score1,fill=group))+
    geom_boxplot(color="black",outlier.size = 0,outlier.color = NA,
                 notch = T,notchwidth = 0.6,lwd=0.1)+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..),outlier.size = 0,notch = T,notchwidth = 0.6,
                 color="white"
    )+
    scale_fill_manual(values=CT.col)+
    scale_color_manual(values=c('white','white'))+
    #scale_x_discrete(breaks=cell_order,labels=cell_order)+
    #theme_bw()+
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,color="black",size=15,face="plain"),
          axis.text.y=element_text(size=15,color="black",face="plain"),
          axis.title.x = element_blank(),
          axis.title.y=element_text(size=15,color="black",face="plain"),
          #axis.line=element_line(color="black",size = unit(0.25,'inches')),
          panel.grid=element_blank(),
          panel.background=element_blank())+
    stat_compare_means(aes(label = paste0( ..p.format..)),angle=45,label.y=max(data$score1)-0.01 )+
    labs(title=i,x="celltype", y="score", fill="group")
  print(p)
###violin plot###
    p <- ggplot(data,aes(group,score1,colour=group))+
       geom_violin(aes(fill=group),cex=0.5,trim = FALSE)+
	   scale_fill_manual(values=CT.col) +
	   geom_boxplot(width=0.2,cex=0.5,outlier.shape = NA)+
       stat_compare_means(label = "p.format", label.x = 1.5,
                     comparisons=list(c("Young","Old"))) +
       facet_wrap(~celltype, scales = "free",nrow=2)+
       scale_color_manual(values=c('#000000','#000000','#000000')) +
       theme(axis.text.x = element_text(angle = 15, hjust = 1))+
       theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))+
    labs(title=i,x="celltype", y="score", fill="group")
  print(p)
  ####split violin plot####
p1=ggplot(data,aes(x=celltype,y = score1,fill=group))+theme_bw() +theme(panel.grid = element_blank())+
    introdataviz::geom_split_violin(trim = T,colour=NA)+
    geom_point(stat = 'summary',fun=mean,size=1,position = position_dodge(width = 0.3))+
    scale_fill_manual(values =  c('#379F7C','#F3A15E'))+
    stat_summary(fun.min = function(x){quantile(x)[2]},fun.max = function(x){quantile(x)[4]},
                 geom = 'errorbar',color='black',width=0,size=0.7,position = position_dodge(width = 0.3))+ylab("Score")+
  theme(axis.title.x = element_text(size = 0,color = 'black'),axis.title.y = element_text(size = 12,colour = 'black'),axis.text.x = element_text(size=12,angle = 90,hjust = 0.5,color = 'black'),axis.text.y = element_text(size=10,color = 'black'),legend.text=element_text(size=10,color = 'black'),legend.title=element_text(size=0))+
    ggpubr::stat_compare_means(aes(group = group),label = "p.format",method = "wilcox.test",hide.ns = F)
    print(p1)
  gs_means <- data %>%
    group_by(group) %>%
    summarise(heightIn = median(score1))
  ####绘制在两个x轴上的密度图
  p2 <-ggplot(data, aes(score1,group,fill=group)) +
    geom_density_ridges(alpha = 0.9, na.rm=TRUE,color = NA)+
    stat_compare_means(comparisons=c("Young","Old"),method = "t.test",label.y = 1.5)+
    scale_fill_manual(values=CT.col) +
    geom_vline(data = gs_means,aes(xintercept = heightIn,colour = group), linetype = "dashed", size = 1)+
    scale_color_manual(values=CT.col)+
    theme_ridges() +
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,color="black",size=15,face="plain"),
          axis.text.y=element_text(size=15,color="black",face="plain"),
          axis.title.x = element_text(size=15,color="black",face="plain"),
          axis.title.y=element_text(size=15,color="black",face="plain"),
          #axis.line=element_line(color="black",size = unit(0.25,'inches')),
          panel.grid=element_blank(),
          panel.background=element_blank())+
    labs(title=i,x="score", y="density", fill="group")
  print(p2)
  dev.off()
}




###cell markers####
library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggpubr)
library(ggplot2)
library(ggridges)
library(reshape2)
library(reshape)
library(pheatmap)
library(ggplot2)
adrenal <- readRDS('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_adrenal_celltype_final.rds')
gene_list=read.csv("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/top_av0.58_gene.csv")
DefaultAssay(adrenal) <- "RNA"
adrenal <- NormalizeData(object =adrenal, normalization.method = "LogNormalize")
adrenal$group <- adrenal$age

CT.col <- c('#379F7C','#F3A15E')
celltypes <- c("ZG","Transitional cell","ZF", "ZR","Chromaffin cell","Neuron", "FB", "SMC", "EC","Adi", "TC", "Mac")
setwd("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/10_gene_score/")
for (i in celltypes ){ 
  genes <- subset(gene_list, cluster == i)
  genes <- genes$gene
  genes <- list(as.character(genes))
  tmp <- subset(adrenal, celltype == i)
  tmp <- AddModuleScore(tmp, features=genes, seed=15, name="score")
  data<- FetchData(tmp,vars = c("group","celltype","score1","sample"))
  assign(i,data)
}


data <- rbind(get("ZG"),get("Transitional cell"),get("ZF"), get("ZR"),get("Chromaffin cell"), get("FB"), get("SMC"), get("EC"), get("TC"), get("Mac"), get("Neuron"),get("Adi"))
data$celltype <- factor(data$celltype, 
                          levels= celltypes , ordered=TRUE)
  data$group <- factor(data$group, 
                       levels=c("Young","Old"), ordered=TRUE)
  pdf(paste0("cell_identity_celltype_all_markers.pdf"), width = 8, height =8)
  ###box plot
  p <- ggplot(data,aes(celltype,score1,fill=group))+
    geom_boxplot(color="black",outlier.size = 0,outlier.color = NA,
                 notch = T,notchwidth = 0.6,lwd=0.1)+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..),outlier.size = 0,notch = T,notchwidth = 0.6,
                 color="white"
    )+
    scale_fill_manual(values=CT.col)+
    scale_color_manual(values=c('white','white'))+
    #scale_x_discrete(breaks=cell_order,labels=cell_order)+
    #theme_bw()+
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,color="black",size=15,face="plain"),
          axis.text.y=element_text(size=15,color="black",face="plain"),
          axis.title.x = element_blank(),
          axis.title.y=element_text(size=15,color="black",face="plain"),
          #axis.line=element_line(color="black",size = unit(0.25,'inches')),
          panel.grid=element_blank(),
          panel.background=element_blank())+
    stat_compare_means(aes(label = paste0( ..p.format..)),angle=45,label.y=max(data$score1)-0.01 )+
    theme_bw()+
    labs(title="Cell identity",x="celltype", y="score", fill="group")
  print(p)
  ###violin plot
    p <- ggplot(data,aes(group,score1,colour=group))+
       geom_violin(aes(fill=group),cex=0.5,trim = FALSE)+
	   scale_fill_manual(values=CT.col) +
	   geom_boxplot(width=0.2,cex=0.5,outlier.shape = NA)+
       stat_compare_means(label = "p.format", label.x = 1.5,
                     comparisons=list(c("Young","Old"))) +
       facet_wrap(~celltype, scales = "free",nrow=3)+
       scale_color_manual(values=c('#000000','#000000','#000000')) +
       theme(axis.text.x = element_text(angle = 15, hjust = 1))+
       theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'none',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))+
        theme_bw()
print(p)
