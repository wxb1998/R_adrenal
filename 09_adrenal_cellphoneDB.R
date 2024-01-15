
library(Seurat)

sm_integrated<-readRDS("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_adrenal_celltype_final.rds")
sm_integrated@meta.data$celltype=sm_integrated@active.ident
sm_integrated_o<-subset(sm_integrated,age=="Old")
sm_integrated_y<-subset(sm_integrated,age=="Young")

dir.create('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/11_cellphoneDB/Old')
setwd('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/11_cellphoneDB/Old')
###old###
write.table(as.matrix(sm_integrated_o@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sm_integrated_o@meta.data), sm_integrated_o@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)


###young###
dir.create('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/11_cellphoneDB/Young')
setwd('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/11_cellphoneDB/Young')
write.table(as.matrix(sm_integrated_y@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sm_integrated_y@meta.data), sm_integrated_y@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown"
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)


####查找年轻年老特异的分子对
library(reshape)
library(ggplot2)
library(dplyr)

read_in <- function(group){
  Mean <- read.table(paste0( '/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/11_cellphoneDB/',group,'/out/significant_means.txt'), sep='\t', header = T)
  Mean <- Mean[, -c(3:4,10:12)]
  Mean <- melt(Mean, id=c('id_cp_interaction', 'interacting_pair', 'gene_a', 'gene_b','secreted', 'receptor_a', 'receptor_b'))
  Mean <- plyr::rename(Mean, c(variable='cell_pair', value='mean'))
  Mean <- Mean[!is.na(Mean$mean),]
  
  pvalue <- read.table(paste0('/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/11_cellphoneDB/',group,'/out/pvalues.txt'), sep='\t', header = T)
  pvalue <- pvalue[, -c(3:11)]
  pvalue <- melt(pvalue, id=c('id_cp_interaction', 'interacting_pair'))
  pvalue <- plyr::rename(pvalue, c(variable='cell_pair', value='pvalue'))
  
  df <- left_join(Mean, pvalue[], by=c('id_cp_interaction', 'interacting_pair', 'cell_pair')) %>% mutate(Group=group)
  df[df$pvalue==0,]$pvalue <- 0.001
  
  return(df)
}

get_group_unique <- function(df1, df2){
  df1$ident <- paste0(df1$id_cp_interaction, '_', df1$cell_pair)
  df2$ident <- paste0(df2$id_cp_interaction, '_', df2$cell_pair)
  overlap_id <- intersect(df1$ident, df2$ident)
  shared <- subset(df1, ident %in% overlap_id)
  unique1 <- subset(df1, ident %in%
                      setdiff(as.character(df1$ident), overlap_id))
  unique1$logmean <- log2(unique1$mean)*(-1)
  unique1$logP <- log10(unique1$pvalue)*(-1)
  unique2 <- subset(df2, ident %in%
                      setdiff(as.character(df2$ident), overlap_id))
  unique2$logmean <- log2(unique2$mean)*(-1)
   unique2$logP <- log10(unique2$pvalue)*(-1)
  share <- inner_join(df1, df2, by='ident')
 
  share_1 <- share[c(1:12)]
  share_2 <- share[c(13:23,12)]
  colnames(share_1) <- colnames(unique2)[1:12]
  colnames(share_2) <- colnames(unique2)[1:12]
  share <- rbind(share_1,share_2)
  share$logmean <- log2(share$mean)*(-1)
   share$logP <- log10(share$pvalue)*(-1)
  unique1$DE <- "down"
  unique2$DE <- "up"
  share$DE <- "share"
  return(list(unique1, unique2, share))
}


###开始读取数据###
  Young = read_in('Young')
  Old = read_in('Old')

  O_Y <- get_group_unique(Young, Old)


  YC_unique_1 <- O_Y[[1]]
  OC_unique_1 <- O_Y[[2]]
  Aging_share <- O_Y[[3]]
 
  tmp.last <- rbind(YC_unique_1,OC_unique_1,Aging_share)
  tmp.last$Tissue <- "Adrenal"
write.csv(tmp.last,paste0("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/11_cellphoneDB/adrneal_cellphoneDB_pairs.csv"),quote=F,row.names=F)


