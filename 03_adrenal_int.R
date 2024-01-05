library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(patchwork)
library(sctransform)
set.seed(5)
dd <-c('YM1','YM2','YM3','OM1','OM2','OM3','OM4')
for (i in 1:length(dd)){
  tmp <- readRDS(paste0("/data/home/quj_lab/wangxuebao/01_results/01_adrenal/03_qc/",dd[i], "_final_qc.rds"))
  tmp <- SCTransform(tmp, verbose = FALSE)
  tmp$sample <- tmp$orig.ident
  tmp$sex <- as.factor(ifelse(substr(tmp$sample,2,2)=='M','Male','Female'))
  tmp$age <- as.factor(ifelse(substr(tmp$sample,1,1)=='Y','Young','Old'))
  tmp$group <- paste0(tmp$age,'_',tmp$sex)
  assign(dd[i],tmp)
}


int.list <- list(get("YM1"), get("YM2"), get("YM3"), get("OM1"), get("OM2"), get("OM3"), get("OM4"))
int.features <- SelectIntegrationFeatures(object.list = int.list, nfeatures = 3000)
int.list <- PrepSCTIntegration(object.list = int.list, anchor.features = int.features)
int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT", anchor.features = int.features)
adrenal <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT")
DefaultAssay(adrenal) <- "integrated" 
saveRDS(adrenal, file = '/data/home/quj_lab/wangxuebao/01_results/01_adrenal/04_int/adrenal_final_int.rds')


adrenal <- RunPCA(adrenal)
#adrenal <- JackStraw(adrenal, dims = 30)
#adrenal <- ScoreJackStraw(adrenal, dims = 1:30)
adrenal <- RunUMAP(adrenal, reduction = "pca", dims = 1:30)
pdf("/data/home/quj_lab/wangxuebao/01_results/01_adrenal/04_int/adrenal_pca.pdf",width=16,height=12)
#JackStrawPlot(adrenal, dims = 1:30)
ElbowPlot(adrenal, ndims=40)
DimHeatmap(adrenal, dims = 1:30, cells = 500, balanced = TRUE)
DimPlot(adrenal, reduction = "umap", label = TRUE)
DimPlot(adrenal, reduction = "umap", pt.size = 0.4, group.by = "sample")
dev.off()
pdf("/data/home/quj_lab/wangxuebao/01_results/01_adrenal/04_int/adrenal_heatmap.pdf",width=20,height=60)
DimHeatmap(adrenal, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()
saveRDS(adrenal, file = '/data/home/quj_lab/wangxuebao/01_results/01_adrenal/04_int/adrenal_pca.rds')

  
