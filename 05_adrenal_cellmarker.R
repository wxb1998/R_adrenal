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
###??ȡ????
adrenal <- readRDS("/data/home/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/m_adrenal_umap-final.rds")
DefaultAssay(adrenal) <- "RNA"
plan("multicore", workers = 28)
options(future.globals.maxSize = 100000 * 1024^2)#100000MB~=100G
adrenal <- NormalizeData(object =adrenal, normalization.method = "LogNormalize")
pdf("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/adrenal_fenxing_cortex-final.pdf", width=24,height=18)
###casu
FeaturePlot(adrenal,raster=FALSE, features = c("RSPO3","GLI1","WNT4","SHH","SF1","DAX1"))
VlnPlot(adrenal, features = c("RSPO3","GLI1","WNT4","SHH","SF1","DAX1"), pt.size = 0.5)
## ZG
FeaturePlot(adrenal,raster=FALSE, features = c("VSNL1","CYP11B2"))
VlnPlot(adrenal, features = c("VSNL1",  "CYP11B2"), pt.size = 0)
# ZF
FeaturePlot(adrenal, raster=FALSE,features = c("CYP11B1", "CYP17A1",'FDX1','MEG3'))
VlnPlot(adrenal, features = c("CYP11B1", "CYP17A1",'FDX1','MEG3'), pt.size = 0)
# ZR
FeaturePlot(adrenal,raster=FALSE, features = c("SULT2A1", "CYB5A"))
VlnPlot(adrenal, features = c("SULT2A1", "CYB5A"), pt.size = 0)
dev.off()

pdf("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/adrenal_fenxing_immune-final.pdf", width=24,height=18)
# Macrophages
FeaturePlot(adrenal, raster=FALSE,features = c("MSR1", "MRC1", "RGS1", "CD68","CD163", "VSIG4","CD86"))
VlnPlot(adrenal, features = c("MSR1", "MRC1", "RGS1", "CD68","CD163", "VSIG4","CD86"), pt.size = 0)
# T
FeaturePlot(adrenal, raster=FALSE,features = c("CD247", "CD69","TRAF3IP3",'CD3D','CD3E', 'CD8A'))
VlnPlot(adrenal, features = c("CD247", "CD69","TRAF3IP3",'CD3D','CD3E', 'CD8A'), pt.size = 0)
# B
FeaturePlot(adrenal, raster=FALSE,features = c('CD79A', 'CD79B', 'MZB1', 'IGHG2', 'MS4A1','CD27','PTPRC'))
VlnPlot(adrenal, features = c('CD79A', 'CD79B', 'MZB1', 'IGHG2', 'MS4A1','CD27','PTPRC'), pt.size = 0)
# NK&Neutrophil
#FeaturePlot(adrenal, features = c('SI00A8', 'SI00A9','GZMB', 'GNLY', 'KLRF'))
#VlnPlot(adrenal, features = c('SI00A8', 'SI00A9','GZMB', 'GNLY', 'KLRF'), pt.size = 0)
dev.off()

pdf("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/adrenal_fenxing_medulla-final.pdf",width=24,height=18)
# Medulla cell
FeaturePlot(adrenal,raster=FALSE, features = c("CHGB", "CHGA", "DLK1", "PNMT", "DBH", "TH"))
VlnPlot(adrenal, features = c("CHGB", "CHGA", "DLK1", "PNMT", "DBH", "TH"), pt.size = 0)
# Neural 
FeaturePlot(adrenal,raster=FALSE, features = c("PMP22",  "SOX10"))
VlnPlot(adrenal, features = c("PMP22", "SOX10"), pt.size = 0)
# sympathoblast
FeaturePlot(adrenal, raster=FALSE,features = c('FOXD3','SOX10','ELAVL3', 'ELAVL4','ISL1','TH','PNMT','STMN2','STMN4'))
VlnPlot(adrenal, features = c('FOXD3','SOX10','ELAVL3', 'ELAVL4','ISL1','TH','PNMT','STMN2','STMN4'), pt.size = 0)
#Chromaffin cell
FeaturePlot(adrenal,raster=FALSE, features = c('CHGB','CHGA','NPY', 'SCG2','DLK1'))
VlnPlot(adrenal, features = c('CHGB','CHGA','NPY', 'SCG2','DLK1'), pt.size = 0)
# glia cell
FeaturePlot(adrenal, raster=FALSE,features = c('SPP1','PLP1','CRYAB'))
VlnPlot(adrenal, features = c('SPP1','PLP1','CRYAB'), pt.size = 0)
dev.off()
pdf("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/adrenal_fenxing-final.pdf",width=24,height=18)
# VSMC
FeaturePlot(adrenal,raster=FALSE, features = c("MYH11", "ACTA2", "TAGLN"))
VlnPlot(adrenal, features = c("MYH11", "ACTA2", "TAGLN"), pt.size = 0)
# EC
FeaturePlot(adrenal, raster=FALSE,features = c("PLVAP", "EGFL7", "ADGRL4", "FLT1", "CDH5", "TM4SF1", "PECAM1"))
VlnPlot(adrenal, features = c("PLVAP", "EGFL7", "ADGRL4", "FLT1", "CDH5", "TM4SF1", "PECAM1"), pt.size = 0)
# Fibroblast
FeaturePlot(adrenal, raster=FALSE,features = c("COL1A2", "COL3A1", "DCN"))
VlnPlot(adrenal, features = c("COL1A2", "COL3A1", "DCN"), pt.size = 0)
# pericyte
FeaturePlot(adrenal, raster=FALSE,features = c("CSPG4", "MCAM", "TRPC6","PDGFRB", "RGS5"))
VlnPlot(adrenal, features = c("CSPG4", "MCAM", "TRPC6","PDGFRB", "RGS5"), pt.size = 0)
# Adipocyte
FeaturePlot(adrenal, raster=FALSE,features = c("LEPR", "ADIPOQ","KLB", "CAR3"))
VlnPlot(adrenal, features = c("LEPR", "ADIPOQ","KLB", "CAR3"), pt.size = 0)
dev.off()

pdf("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/adrenal_ct_cellmarker.pdf", width=24,height=18)
FeaturePlot(adrenal,raster=FALSE, features = c("LEF1",'VSNL1',"CYP11B2","FDX1"),cols = c("lightgrey", "#DC3535"))
FeaturePlot(adrenal,raster=FALSE, features = c("CYP11A1","CYP17A1","CYP51A1","SULT2A1"),cols = c("lightgrey", "#DC3535"), min.cutoff = 1, max.cutoff = 3)
FeaturePlot(adrenal,raster=FALSE, features = c('DBH', "CHGB", "PMP22", 'SOX10'),cols = c("lightgrey", "#DC3535"))
FeaturePlot(adrenal,raster=FALSE, features = c("COL1A2", 'DCN','MYH11',"ACTC2"),cols = c("lightgrey", "#DC3535"))
FeaturePlot(adrenal,raster=FALSE, features = c("TAGLN", 'FLT1',"PECAM1",'LEPR'),cols = c("lightgrey", "#DC3535"))
FeaturePlot(adrenal,raster=FALSE, features = c('ADIPOQ',"TRAF3IP3",'CD247','CD163'),cols = c("lightgrey", "#DC3535"))
FeaturePlot(adrenal,raster=FALSE, features = c('MSR1'),cols = c("lightgrey", "#DC3535"))
dev.off()

png("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/adrenal_ZR_cellmarker.png", width=12, height=9,units = "in",res=600)
FeaturePlot(adrenal,raster=FALSE, features = c("SULT2A1"),cols = c("lightgrey", "#DC3535"), min.cutoff = 0, max.cutoff = 1)
dev.off()

png("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/adrenal_ZR2_cellmarker.png", width=12, height=10,units = "in",res=600)
FeaturePlot(adrenal,raster=FALSE, features = c("CYP51A1"),cols = c("lightgrey", "#DC3535"), min.cutoff = 1.8, max.cutoff = 2)
dev.off()
###celltype marker genes###
pdf("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/05_umap/celltype_marker.pdf", width=24,height=18)
FeaturePlot(adrenal,raster=FALSE, features = c("LEF1",'VSNL1',"CYP11B2","FDX1","CYP11A1","CYP17A1","CYP51A1","SULT2A1",'DBH', "CHGB","COL1A2", 'DCN','MYH11',"ACTC2", "TAGLN", 'FLT1',"PECAM1","TRAF3IP3",'CD247','CD163','MSR1', "PMP22", 'SOX10','LEPR','ADIPOQ'))
VlnPlot(adrenal, features = c("LEF1",'VSNL1',"CYP11B2","FDX1","CYP11A1","CYP17A1","CYP51A1","SULT2A1",'DBH', "CHGB","COL1A2", 'DCN','MYH11',"ACTC2", "TAGLN", 'FLT1',"PECAM1","TRAF3IP3",'CD247','CD163','MSR1', "PMP22", 'SOX10','LEPR','ADIPOQ'), pt.size = 0)
dev.off()
