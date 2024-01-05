library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(patchwork)
set.seed(5)

sample <-c('OM1','OM2','OM3','OM4','YM1','YM2','YM3')#,,'OF1','OF3','OM1','OM3','OF2','OF4','OM2','OM4','YF1','YF2','YF4','YM1','YM2','YM3'
rate <- c(0.08, 0.08,0.08, 0.069, 0.08,0.08, 0.054)
doublet.prop <- data.frame(matrix(nrow=0, ncol=3))
colnames(doublet.prop) <- c("Sample", "Number","Doublet_prop")
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 100000 * 1024^2)#100000MB~=100G

for (i in 1:length(sample)){
  tmp <- readRDS(paste0("/data/home/quj_lab/wangxuebao/01_results/01_adrenal/03_qc/", sample[i], "_before_qc.rds"))
  
  ## pK Identification (ground-truth) -----------
  sweep.res.list_tmp <- paramSweep_v3(tmp, PCs = 1:20, sct = T)
  sweep.stats_tmp <- summarizeSweep(sweep.res.list_tmp, GT = FALSE)
  bcmvn <- find.pK(sweep.stats_tmp)
  pk_v <- as.numeric(as.character(bcmvn$pK))
  pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]###��ȡ����pkֵ
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- tmp@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(rate[i]*ncol(tmp))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  tmp <- doubletFinder_v3(tmp, PCs = 1:20, pN = 0.25, pK = pk_good, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
  DF.name = colnames(tmp@meta.data)[grepl("DF.classification", colnames(tmp@meta.data))]
  p5.dimplot=cowplot::plot_grid(ncol = 2, DimPlot(tmp, group.by = "orig.ident") + NoAxes(), 
                                DimPlot(tmp, group.by = DF.name) + NoAxes())
  ggsave(filename=paste0("/data/home/quj_lab/wangxuebao/01_results/01_adrenal/03_qc/", sample[i], "_doublet_dimplot.pdf"),plot=p5.dimplot, width = 10, height = 6)
  
  
  # find doublet
  tmp_qc=tmp[, tmp@meta.data[, DF.name] == "Singlet"]
  saveRDS(tmp_qc, file = paste0("/data/home/quj_lab/wangxuebao/01_results/01_adrenal/03_qc/",sample[i], "_final_qc.rds"))
  tmp[["doublet"]] <- tmp[[paste("DF.classifications_0.25", pk_good, nExp_poi.adj, sep="_")]]
  prop <- nExp_poi.adj/length(tmp@meta.data$cluster)
  prop.tmp <- data.frame(Sample=sample[i], Number=nExp_poi.adj, Doublet_prop=prop)
  doublet.prop <- rbind(doublet.prop, prop.tmp)
  pdf(paste0("/data/home/quj_lab/wangxuebao/01_results/01_adrenal/03_qc/", sample[i],"_after_plot.pdf"),width = 10, height = 6)
  p=VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
  print(p)
  dev.off()
}
write.table(doublet.prop,file ="/data/home/quj_lab/wangxuebao/01_results/01_adrenal/03_qc/doublet.prop.txt", sep = "\t" )
