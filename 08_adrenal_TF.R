
library(Seurat)
library(tidyverse)
library(future)
library(patchwork)
library(SCENIC)
library(RcisTarget)
library(AUCell)
library(arrow)

####将全部差异基因来计算转录因子####
m_adrenal <- readRDS(file = "/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_adrenal_celltype_final.rds")
plan("multicore", workers = 20)
options(future.globals.maxSize = 100000 * 1024^2)#100000MB~=100G
###全部基因转录因子预测
deg = read.csv(paste("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/08_DEG/","OM_YM_deg_0.35.csv",sep=""))
##########按照细胞类型跑转录因子
#####"ZG","Transitional","ZF", "ZR","Chromaffin cell", "Fibroblast", "VSMC","T_cell", "EC", "Neuron", "Macrophage","Adipocyte"
dir.create(paste0("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/09_TF/log2fc0.35/"))
for (i in c("ZG","Transitional cell","ZF", "ZR","Chromaffin cell","Neuron", "FB", "SMC", "EC","Adi", "TC", "Mac"))
{
  tmp <- subset(m_adrenal, celltype==i)
  exp.mat <- GetAssayData(tmp, slot='data')
  deg_tmp=subset(deg, celltype==i)####筛选差异基因
  genes <- unique(as.character(deg_tmp$gene))
  exp.mat <- as.matrix(exp.mat)
  exp.mat <- as.matrix(exp.mat[genes,])
  ##########设置工作目录
  dir.create(paste0("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/09_TF/log2fc0.35/all_",i,"/"))
  setwd(paste0("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/09_TF/log2fc0.35/all_",i,"/"))
  #初始化SCENIC环境，设置分析环境
  org <- "hgnc"
  dbDir <- "/data/home/quj_lab/wangxuebao/02_software/SCENIC/hg19_cisTarget_databases" # RcisTarget databases location
  myDatasetTitle <- "all degs on m_adrenal" # choose a name for your analysis
  data(defaultDbNames)
  dbs <- defaultDbNames[[org]]
  scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)  ##nCores为10，开启十个线程
  saveRDS(scenicOptions, file="int/scenicOptions.Rds")
  
  # Gene filter
  genesKept <- geneFiltering(exp.mat, scenicOptions=scenicOptions,
                             minCountsPerGene=3*.01*ncol(exp.mat),
                             minSamples=ncol(exp.mat)*.01)
  exprMat_filtered <- exp.mat[genesKept, ]
  dim(exprMat_filtered)
  
  
  judgement = tryCatch({
    
    # Run Genie3
    runCorrelation(exprMat_filtered, scenicOptions)
    #correct pipeline
    runGenie3(exprMat_filtered, scenicOptions, nParts = 10)
    # Run the remaining
    scenicOptions@settings$verbose <- TRUE
    scenicOptions@settings$nCores <- 1
    scenicOptions@settings$seed <- 123
    runSCENIC_1_coexNetwork2modules(scenicOptions)
    runSCENIC_2_createRegulons(scenicOptions)
    runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
    
    print(paste0('----------------------',i,'_Genie3 finished','----------------------'))
  },  error=function(e) e
  )
  if(inherits(judgement, 'error')) {
    print(paste0(i,'_error'))
    next}
}
