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
set.seed(2)
adrenal <- readRDS(file = "/data/home/quj_lab/wangxuebao/01_results/01_adrenal/04_int/adrenal_pca.rds")
###
plan("multicore", workers = 40)
options(future.globals.maxSize = 100000 * 1024^2)#100000MB~=100G
DefaultAssay(adrenal) <- "integrated"
adrenal <- RunPCA(adrenal)
adrenal <- RunUMAP(adrenal, reduction = "pca",dims = 1:16)
adrenal <- FindNeighbors(adrenal, reduction = "pca", dims = 1:16)
adrenal <- FindClusters(adrenal, resolution = 2.4)
pdf("/data/home/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/m_adrenal_umap-final.pdf",width=16,height=12)
DimPlot(adrenal, reduction = "umap", label = TRUE, repel = TRUE,raster=FALSE) + NoLegend()
DimPlot(adrenal, reduction = "umap", label = TRUE, repel = TRUE,raster=FALSE)
DimPlot(adrenal, reduction = "umap", split.by = "sample",pt.size = 0.4,label.size = 6,ncol=3,raster=FALSE)
dev.off()
saveRDS(adrenal, file = "/data/home/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/m_adrenal_umap-final.rds")

###
pdf('/data/home/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/adrenal_umap_bysample_final.pdf',width=16,height=12)
#DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "sample",raster=FALSE, cols= c('#d3ba68','#014d64','#d5695d','#01a2d9','#95e1d3','#00887d','#adadad','#00F5FF', '#54FF9F', '#000080','#458B00','#FFD700', '#A0522D', '#9932CC'))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "sample", cols= c('#d3ba68','#014d64','#d5695d','#01a2d9','#95e1d3','#00887d','#adadad'))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "sample", cols= c('#d3ba68',"NA", "NA", "NA", "NA", "NA", "NA"))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "sample", cols= c("NA",'#014d64', "NA", "NA", "NA", "NA", "NA"))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "sample", cols= c("NA","NA", '#d5695d', "NA", "NA", "NA", "NA"))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "sample", cols= c("NA","NA", "NA",'#01a2d9', "NA", "NA", "NA"))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "sample", cols= c("NA","NA", "NA", "NA", '#95e1d3', "NA", "NA"))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "sample", cols= c("NA","NA", "NA", "NA", "NA", '#00887d', "NA"))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "sample", cols= c("NA","NA", "NA", "NA", "NA","NA",'#adadad'))
dev.off()
####
pdf('/data/home/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/adrenal_umap_byage_final.pdf',width=12,height=9)
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "age",raster=FALSE, cols= c('#014d64','#d5695d'))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "age",raster=FALSE, cols= c('#014d64',"NA"))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "age",raster=FALSE, cols= c("NA",'#d5695d'))
dev.off()

####
pdf('/data/home/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/adrenal_umap_bygroup_final.pdf',width=12,height=9)
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "group",raster=FALSE, cols= c('#d3ba68','#014d64','#d5695d','#9932CC'))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "group",raster=FALSE, cols= c("#d3ba68","NA", "NA", "NA"))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "group",raster=FALSE, cols= c("NA","#014d64", "NA", "NA"))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "group",raster=FALSE, cols= c("NA","NA", "#d5695d", "NA"))
DimPlot(adrenal, reduction = "umap", pt.size = 0.3, group.by = "group",raster=FALSE, cols= c("NA","NA", "NA", "#9932CC"))
dev.off()

####
Idents(adrenal) <- "seurat_clusters"
pdf('/data/home/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/violin_plot_final.pdf',width=20,height=10)
VlnPlot(adrenal, features = "nFeature_RNA",pt.size = 0)
VlnPlot(adrenal, features = "percent.MT",pt.size = 0)
dev.off()

###
pdf("/data/home/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/barplot_qc_final.pdf",width=20,height=6)
cell.prop = data.frame(table(adrenal$seurat_clusters,adrenal$sample))
colnames(cell.prop) = c('seurat_clusters', 'sample','number')
cell.prop = cell.prop %>% group_by(sample) %>% mutate(sample.sum = sum(number))
cell.prop$sample.prop = cell.prop$number / cell.prop$sample.sum * 100
#cell.prop$sample <- factor(cell.prop$sample,
                           #levels=samples, ordered=TRUE)
cell.prop = cell.prop %>% group_by(seurat_clusters) %>% mutate(seurat_clusters.sum = sum(sample.prop))
cell.prop$prop = cell.prop$sample.prop / cell.prop$seurat_clusters.sum * 100
ggplot(data =cell.prop, mapping = aes(x = seurat_clusters,y=prop,fill=sample))+
  scale_fill_manual(values = c('#d3ba68','#014d64','#d5695d','#01a2d9','#95e1d3','#00887d','#adadad','#00F5FF', '#54FF9F', '#000080','#458B00','#FFD700', '#A0522D', '#9932CC'))+
  geom_bar(stat = 'identity',position = "dodge",width=0.6) 
dev.off()


pdf("/data/home/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/barplot_qc8_final.pdf",width=20,height=6)
cell.prop = data.frame(table(adrenal$seurat_clusters, adrenal$sample))
colnames(cell.prop) = c('seurat_clusters', 'sample', 'number')
cell.prop$group <- str_sub(cell.prop$sample,1,2)
cell.prop$group = factor(cell.prop$group, levels = c('YM', 'OM'))
cell.prop = cell.prop %>% group_by(group,seurat_clusters) %>% mutate(group.sum = sum(number))
cell.prop$prop = cell.prop$number / cell.prop$group.sum * 100
ggplot(data =cell.prop, mapping = aes(x = seurat_clusters,y=prop,fill=sample))+
  scale_fill_manual(values = c('#d3ba68','#014d64','#d5695d','#01a2d9','#95e1d3','#00887d','#adadad','#00F5FF', '#54FF9F', '#000080','#458B00','#FFD700', '#A0522D', '#9932CC'))+
  geom_bar(stat = 'identity',position = "dodge",width=0.6) 
dev.off()



