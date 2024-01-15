library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(SeuratData)
library(patchwork)
set.seed(2)
adrenal<- readRDS("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_adrenal_celltype_final.rds")
celltypes <- c("ZG","Transitional cell","ZF", "ZR","Chromaffin cell","Neuron", "FB", "SMC", "EC","Adi", "TC", "Mac")
DefaultAssay(adrenal) <- "RNA"

adrenal <- NormalizeData(object = adrenal, normalization.method = "LogNormalize")

adrenal@meta.data$celltype<-adrenal@active.ident
adrenal$group <- substr(adrenal$sample,1,2)
OM_YM <- data.frame()
Idents(adrenal) <- paste(adrenal$celltype, adrenal$group, sep='_')
for (cell in celltypes){
  tmp <- FindMarkers(adrenal, ident.1=paste0(cell, '_OM'), ident.2=paste0(cell, '_YM'))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$compare <- 'OM_YM'
  OM_YM <- rbind(OM_YM, tmp)
  print(paste0(cell, ' is finished'))
}


OM_YM <- subset(OM_YM,str_sub(gene,1,7) != "ENSMFAG" )
###log2FC 0.35
write.csv(OM_YM,file=paste("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/08_DEG/OM_YM.csv",sep=""))
OM_YM.deg <- subset(OM_YM, p_val_adj<0.05 & abs(avg_log2FC)>0.35)
table(OM_YM.deg$celltype)
adrenal.up <- subset(OM_YM.deg,avg_log2FC>0)
adrenal.down <- subset(OM_YM.deg,avg_log2FC<0)
write.csv(adrenal.up,file = '/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/08_DEG/OM_YM_deg_up_0.35.csv')
write.csv(adrenal.down,file = '/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/08_DEG/OM_YM_deg_down_0.35.csv')
write.csv(OM_YM.deg,file=paste("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/08_DEG/","OM_YM_deg_0.35.csv",sep=""))



#########################The number of DEGs across different cell types#############################################################
library(dplyr)
library(ggplot2)
setwd("D:\\7-Zip\\01_adrenal\\11_DEG\\")
phon <- read.csv('OM_YM_deg_0.35.csv',header=T)
phon$threshold <- as.factor(ifelse(abs(phon$avg_log2FC) >=0,ifelse(phon$avg_log2FC >=0 ,'Up regulated','Down regulated'),'Unchanged')) 
df<-phon%>%
  group_by(celltype,threshold)%>%
  summarise(Counts=n())
celltype_order <- as.data.frame(table(phon$celltype))
colnames(celltype_order) <-c("celltype","freq")
celltype_order <- celltype_order[rev(order(celltype_order$freq)),]
celltype_order <- celltype_order$celltype
df$celltype <- factor(df$celltype,levels=celltype_order)

pdf('adrenal_number_DEG0.35.pdf',width=10,height=6)
ggplot(df,
       aes(x=celltype,y=Counts,fill=threshold))+ 
  geom_bar(stat="identity",alpha = 0.7,width = 0.6)+
  scale_fill_manual(values=c('SteelBlue', 'DarkRed'))+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,color="black",size=15,face="plain"),
        plot.title=element_text(face="bold.italic", size=14, color="black"),
        axis.title=element_text(face="bold.italic",
                                size=14, color="black"),
        axis.text=element_text(face="bold", size=9,
                               color="black"),panel.background=element_blank())+
  labs(x="celltype",
       y="The numeber of up/down-regulated genes")
dev.off()




###Adrenocortical region###
#install.packages("UpSetR")
library(UpSetR)
library(RColorBrewer)
library(ggplot2)
####上调的upset图###
DEG <- read.csv("D:/7-Zip/01_adrenal/11_DEG/OM_YM_deg_0.35.csv")
DEG_tmp <- subset(DEG,avg_log2FC>0)
ZG_tmp <- subset(DEG_tmp,celltype=="ZG")$gene
ZF_tmp <- subset(DEG_tmp,celltype=="ZF")$gene
T_tmp <- subset(DEG_tmp,celltype=="Transitional_cell")$gene 
ZR_tmp <- subset(DEG_tmp,celltype=="ZR")$gene

upset_list <- list(ZG_tmp, T_tmp,ZF_tmp,ZR_tmp)   # 制作Upset图搜所需要的列表文件
names(upset_list) <- c("ZG","Transition cell","ZF", "ZR")
p <- upset(fromList(upset_list),
      nsets = length(upset_list), # 显示数据集的所有数据, nsets = 数值调整可视化数据集数量
      nintersects = 15, # 显示前多少个
      sets = rev(c("ZG","Transition cell","ZF", "ZR")), 
      keep.order = TRUE, # 指定集合或用keep.order = TRUE保持集合按输入的顺序排序
      number.angles = 0, # 交互集合柱状图的柱标倾角
      point.size = 4, # 图中点的大小
      line.size = 1, # 图中连接线粗细
      mainbar.y.label = "Intersection size", # y轴的标签
      main.bar.color = 'black', # y轴柱状图颜色
      matrix.color = "black", # x轴点的颜色
      sets.x.label = "Set size", # x轴的标签
      sets.bar.color=c("ZG"='#F2CD5C',"Transitional cell"="#F8766D","ZF"="#C77CFF","ZR"= "#BA1C96"), # x轴柱状图的颜色; Set1中只有9个颜色，Set3中有12个颜色，Paired中有12个颜色
      mb.ratio = c(0.7, 0.3), # bar plot和matrix plot图形高度的占比
      order.by = "degree", # y轴矩阵排序,如"freq"频率，"degree"程度
      text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1), # 6个参数intersection size title（y标题大小）,intersection size tick labels（y刻度标签大小）, set size title（set标题大小）, set size tick labels（set刻度标签大小）, set names（set 分类标签大小）, numbers above bars（柱数字大小）的设置
      shade.color = "#12507B", # 图中阴影部分的颜色
      #queries=list(list(query = intersects, params = list("votes"), color = "purple", active = T), # 设置自己想要展示的特定组的交集，通过queries参数进行设置，需要展示几个关注组合的颜色，就展示几个
      #             list(query = intersects, params = list("votes","length"), color = "orange", active = T))
)


pdf("D:/7-Zip/01_adrenal/11_DEG/cortex_diff_0.35_up.pdf",width =8, height =5)
print(p)
dev.off()
####下调的upset图###
DEG <- read.csv("D:/7-Zip/01_adrenal/11_DEG/OM_YM_deg_0.35.csv")
DEG_tmp <- subset(DEG,avg_log2FC<0)
ZG_tmp <- subset(DEG_tmp,celltype=="ZG")$gene
ZF_tmp <- subset(DEG_tmp,celltype=="ZF")$gene
T_tmp <- subset(DEG_tmp,celltype=="Transitional_cell")$gene 
ZR_tmp <- subset(DEG_tmp,celltype=="ZR")$gene

upset_list <- list(ZG_tmp, T_tmp,ZF_tmp,ZR_tmp)   # 制作Upset图搜所需要的列表文件
names(upset_list) <- c("ZG","Transition cell","ZF", "ZR")
#names(upset_list) <- colnames(upset_dat[1:9]) 

p <- upset(fromList(upset_list),
      nsets = length(upset_list), # 显示数据集的所有数据, nsets = 数值调整可视化数据集数量
      nintersects = 15, # 显示前多少个
      sets = rev(c("ZG","Transition cell","ZF", "ZR")), 
      keep.order = TRUE, # 指定集合或用keep.order = TRUE保持集合按输入的顺序排序
      number.angles = 0, # 交互集合柱状图的柱标倾角
      point.size = 4, # 图中点的大小
      line.size = 1, # 图中连接线粗细
      mainbar.y.label = "Intersection size", # y轴的标签
      main.bar.color = 'black', # y轴柱状图颜色
      matrix.color = "black", # x轴点的颜色
      sets.x.label = "Set size", # x轴的标签
      sets.bar.color=c("ZG"='#F2CD5C',"Transitional cell"="#F8766D","ZF"="#C77CFF","ZR"= "#BA1C96"), # x轴柱状图的颜色; Set1中只有9个颜色，Set3中有12个颜色，Paired中有12个颜色
      mb.ratio = c(0.7, 0.3), # bar plot和matrix plot图形高度的占比
      order.by = "degree", # y轴矩阵排序,如"freq"频率，"degree"程度
      text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1), # 6个参数intersection size title（y标题大小）,intersection size tick labels（y刻度标签大小）, set size title（set标题大小）, set size tick labels（set刻度标签大小）, set names（set 分类标签大小）, numbers above bars（柱数字大小）的设置
      shade.color = "#12507B", # 图中阴影部分的颜色
      #queries=list(list(query = intersects, params = list("votes"), color = "purple", active = T), # 设置自己想要展示的特定组的交集，通过queries参数进行设置，需要展示几个关注组合的颜色，就展示几个
      #             list(query = intersects, params = list("votes","length"), color = "orange", active = T))
)


pdf("D:/7-Zip/01_adrenal/11_DEG/cortex_diff_0.35_down.pdf",width =8, height =5)
print(p)
dev.off()