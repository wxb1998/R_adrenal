

library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(SeuratData)
library(patchwork)
library(Rserve)
library(reshape2)
library(ggpubr)
library(CellChat)
CT.col =scPalette(20)
set.seed(2)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2)#100000MB~=100G
CT.col = c("ZG"='#F2CD5C',"Transitional cell"="#F8766D","ZF"="#C77CFF","ZR"= "#BA1C96", "Chromaffin cell"="#B05E27", 
"Neuron"="#E3B7A0","FB"="#00BA38" ,"SMC"= "#00887d","EC"= "#c72e29","Adi"="#C1A3A3","TC"= "#016392", "Mac"="#76c0c1")
adrenal <- readRDS(file = "/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/v5/m_adrenal_umap-final.rds")
DefaultAssay(adrenal) <- "RNA"
adrenal <- NormalizeData(object = adrenal, normalization.method = "LogNormalize")
adrenal <-subset(adrenal,seurat_clusters!=35)
###CT.col <- c('#d5695d','#00887d',"#F0E68C","#00A087B2",'#59a77f','#fb832d','#65a479','#adadad',"#E64B35B2","#4DBBD5B2", "#3C5488B2") 
adrenal <- RenameIdents(adrenal,c("0"="ZF", "1"="ZF", "2"="ZG", "3"="ZF","4"= "ZF", "5"="ZF","6"="ZG", "7"="ZG","8"= "ZG","9"= "ZF",
                                  "10"= "ZG", "11"="ZG", "12"="ZG", "13"="Chromaffin cell", "14"="ZF", "15"="ZG","16"= "Chromaffin cell", "17"="ZF", "18"="EC", "19"="FB",
                                  "20"="EC", "21"="Mac", "22"="ZF", "23"="ZR","24"= "Transitional cell", "25"="ZF","26"="ZF", "27"="TC","28"= "FB","29"= "ZG",
                                  "30"="FB", "31"="EC", "32"="Chromaffin cell", "33"="FB","34"= "ZF","36"="FB", "37"="EC","38"= "SMC","39"= "TC",
                                  "40"="ZF", "41"= "ZF" , "42"="Neuron", "43"="Mac", "44"="FB", "45"="Adi", "46"="Chromaffin cell"))
celltypes <- c("ZG","Transitional cell","ZF", "ZR","Chromaffin cell","Neuron", "FB", "SMC", "EC","Adi", "TC", "Mac")
Idents(adrenal) <- factor(Idents(adrenal), levels=celltypes)
UMAP <- DimPlot(adrenal, reduction = "umap",  label = TRUE, pt.size = 0.3,label.size = 4.0,
                cols = CT.col)
ggsave("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/adrenal_celltype_final.pdf", 
       plot = UMAP, width = 14, height = 10)

Idents(adrenal) <- factor(Idents(adrenal), levels=celltypes)
adrenal$celltype <- Idents(adrenal)

markers <- c("LEF1",'VSNL1',"CYP11B2","FDX1","CYP11A1","CYP17A1","CYP51A1","CYB5A","SULT2A1",'DBH', "CHGB", "PMP22", 'SOX10',"COL1A2", 'DCN','MYH11',"ACTC2", "TAGLN", 'FLT1',"PECAM1",'LEPR','ADIPOQ',"TRAF3IP3",'CD247','CD163','MSR1')
DotPlot<-DotPlot(
  adrenal,col.min = 0,
  col.max = 2,
  features= markers, cols=c('grey90', '#C63C3C'))+
  RotatedAxis()
ggsave("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/adrenal_celltype_DotPlot_final.pdf", plot = DotPlot, width = 9, height = 5)
saveRDS(adrenal, file = "/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_adrenal_celltype_final.rds")

##Analysis############################################################################
adrenal@meta.data$celltype<-adrenal@active.ident
write.csv(adrenal@meta.data,"/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_adrenal_Combination_final.csv", row.names =FALSE)




###Read the RDS file
adrenal <- readRDS(file = "/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_adrenal_celltype_final.rds")
CT.col = c("ZG"='#F2CD5C',"Transitional cell"="#F8766D","ZF"="#C77CFF","ZR"= "#BA1C96", "Chromaffin cell"="#B05E27", 
"Neuron"="#E3B7A0","FB"="#00BA38" ,"SMC"= "#00887d","EC"= "#c72e29","Adi"="#C1A3A3","TC"= "#016392", "Mac"="#76c0c1")
  p1 <- VlnPlot(adrenal, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3,cols=CT.col)
pdf(paste0("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/MT_per.pdf"),width = 10, height = 6)
print(p1)
dev.off()

####markers###This step takes a little longer
marker <- FindAllMarkers(adrenal)
# Screen for unannotated genes in cynomolgus monkeys
marker <- subset(marker, p_val_adj<0.05 )
write.csv(marker,"/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/marker_all.csv", row.names =FALSE)
marker <- read.csv("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/marker_all.csv")

top_marker <- subset(marker,p_val_adj<0.05 & avg_log2FC>0.58)
write.csv(top_marker,"/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/top_av0.58_gene.csv", row.names =FALSE)

top50 <- top_marker %>% group_by(cluster) %>% top_n(n=50, wt=avg_log2FC)
write.csv(top50,"/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/top50_gene.csv", row.names =FALSE)

###top_cprtex <- subset(top_1,cluster %in% c("ZG","ZF","ZR"))
top50 <- read.csv("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/top50_gene.csv")
cluster.averages <- AverageExpression(adrenal, return.seurat = TRUE)
cluster.averages
pdf('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/top_gene_heatmap.pdf',width=8,height=7)
DoHeatmap(cluster.averages, features=top50$gene, size = 3, draw.lines = FALSE,group.colors =CT.col)+scale_fill_gradientn(colors = c("#8AB6D6", "white", "#DC143C"))
dev.off()
pdf('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/top_cortex_gene_heatmap_high.pdf',width=10,height=30)
DoHeatmap(cluster.averages, features=top50$gene, size = 3, draw.lines = FALSE,group.colors =CT.col)+scale_fill_gradientn(colors = c("#8AB6D6", "white", "#DC143C"))
dev.off()
# cell proportion #######Plot the proportion of cells for each cell type and each sample
samples <-c('YM1','YM2','YM3','OM1','OM2','OM3','OM4')
adrenal$sample <- factor(adrenal$sample, levels=samples)
cell_number=table(Idents(adrenal),adrenal$sample)
write.csv(cell_number,"/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_adrenal_cell_number_all_final.csv")

pdf('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/adrenal_umap_bygroup_final.pdf',width=10,height=6)
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "group", cols= c('#d3ba68','#014d64'))
DimPlot(adrenal, reduction = "umap", pt.size = 0.2, group.by = "group",split.by = "group",cols=c('#d3ba68','#014d64'),ncol = 2)
dev.off()



######Calculate the cell ratio according to the ratio number
cell.prop = data.frame(table(adrenal$celltype,adrenal$sample))###Cells were counted by cell type and sample
colnames(cell.prop) = c('celltype', 'sample','number')###Name the extracted data frame variable
cell.prop = cell.prop %>% group_by(sample) %>% mutate(sample.sum = sum(number))##Count the number of cells in each summation sample
cell.prop$sample.prop = cell.prop$number / cell.prop$sample.sum * 100###Calculate the proportion of each cell type in each sample
cell.prop$sample <- factor(cell.prop$sample,
                           levels=samples, ordered=TRUE)
cell.prop = cell.prop %>% group_by(celltype) %>% mutate(celltype.sum = sum(sample.prop))##The sum of the proportions for each cell type is counted
cell.prop$prop = cell.prop$sample.prop / cell.prop$celltype.sum * 100

write.csv(cell.prop,"/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_adrenal_cell_prop_all_final.csv")
###Plot the proportion of cells for each cell type and each sample
p=ggplot(cell.prop, aes(celltype, prop, fill=cell.prop$sample)) +
  geom_bar(stat='identity', position = 'dodge', width=0.9) +
  #scale_y_continuous(expand = c(0,0.005)) +
  #scale_fill_manual(values = c( 'dodgerblue1', 'orange2')) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'),
        axis.text = element_text(color='black'),axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_cell_prop_final.pdf", plot = p, width = 8, height = 6)


library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
CT.col = c("ZG"='#F2CD5C',"Transitional cell"="#F8766D","ZF"="#C77CFF","ZR"= "#BA1C96", "Chromaffin cell"="#B05E27", 
"Neuron"="#E3B7A0","FB"="#00BA38" ,"SMC"= "#00887d","EC"= "#c72e29","Adi"="#C1A3A3","TC"= "#016392", "Mac"="#76c0c1")
####Plot the proportion of cells for each cell type and each sample
cell.prop = data.frame(table(adrenal$celltype,adrenal$sample))###Cells were counted by cell type and sample
colnames(cell.prop) = c('celltype', 'sample','number')###Name the extracted data frame variable
cell.prop = cell.prop %>% group_by(sample) %>% mutate(sample.sum = sum(number))##Count the number of cells in each summation sample
cell.prop$sample.prop = cell.prop$number / cell.prop$sample.sum * 100###Calculate the proportion of each cell type in each sample
cell.prop$group <- str_sub(cell.prop$sample,1,2)###Extract the first two letters of the string as the group name
cell.prop$group = factor(cell.prop$group, levels = c('YM', 'OM'))###Order of presentation


cell.prop = cell.prop %>% group_by(group) %>% mutate(group.sum = sum(sample.prop))
cell.prop$prop = cell.prop$sample.prop / cell.prop$group.sum * 100
cell.prop$group = factor(cell.prop$group, levels = c('YM', 'OM'))

pdf('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_cell_prop2.pdf', height=6, width = 4.5)
ggplot(cell.prop, aes(group, prop, fill = celltype)) +
  geom_bar(stat='identity', position = 'fill', width=0.5) +
  scale_fill_manual(values = CT.col) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))

cell.prop1 = cell.prop %>% group_by(celltype,group) %>% mutate(fin_prop = sum(prop))
cell.prop1 <- as.data.frame(unique(cell.prop1[,c("group","celltype","fin_prop")]))
ggplot(cell.prop1, aes(x=group, y=fin_prop, fill = celltype,
                  stratum=celltype, alluvium=celltype)) + 
  geom_col(width = 0.5, color='black')+
  geom_alluvium( width = 0.5,alpha = 0.8,knot.pos = 0)+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  scale_fill_manual(values =CT.col)
dev.off()

pdf('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_cell_prop3.pdf', height=6, width = 8)
ggplot(cell.prop, aes(celltype, prop, fill = group)) +
  geom_bar(stat='identity', position = 'fill', width=0.8) +
  scale_fill_manual(values =CT.col1) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))
dev.off()
####Draw a lollipop chart to show cell proportions
library(ggpubr)
library(ggplot2)
set.seed(1234) 
aa <- as.data.frame(table(adrenal$celltype))

#col =c('#AAAAAA','#8AB6D6',"#d5aca9", "#c6ac8f", "#74d3ae", "#b8dbd9" , "#a6c48a",
      # "#f6e7cb", "#f4e285", "#c6d2ed","#FFB2A7","#98c9a3")
pdf("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/bangbang.pdf",width=8,height=6)
ggdotchart(aa, x="Var1", y="Freq", color = "Var1",
           palette = CT.col,     
           legend = "right",     
           sorting = "descending",   
           add = "segments",   
           rotate = TRUE,       
           dot.size = 10,        
           label = round(aa$Freq),   
           font.label = list(color="black",size=8, vjust=0.5),   
           ggtheme = theme_pubr())+
  theme(axis.text.y = element_text(color="#252422", size=9),
        axis.text.x = element_text(color="#252422", size=9))+    
  #scale_fill_manual(values = col)+ 
  labs(title="The number of different cell types",
       x="celltype", y="counts", color="celltype")
dev.off()

