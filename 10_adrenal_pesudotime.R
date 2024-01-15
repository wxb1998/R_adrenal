library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(ggplot2)
library(future)
library(RColorBrewer)
library(tidyverse)
library(htmlwidgets)  
library(ggpubr)       
library(colorRamps)  
library(magrittr)  
library(monocle)
library(pheatmap)
library(viridis)
set.seed(1998)
setwd("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/13_pesudotime/ZG_ZR")
###创建CellDataSet
m_adrenal <- readRDS(file = "/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/07_celltype/m_adrenal_celltype_final.rds")
m_adrenal <- subset(m_adrenal,celltype == 'ZG' | celltype == 'ZF'| celltype == 'Transition cell'| celltype == 'ZR')
m_adrenal <- subset(m_adrenal,downsample=1100)####每个细胞类型取1100cells，节省计算空间，加快速度，同时要保证每个细胞类型都能取到足够的细胞数，
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix <- as(as.matrix(m_adrenal@assays$RNA@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 
p_data <- m_adrenal@meta.data  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
p_data$celltype <- m_adrenal@active.ident
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(m_adrenal),row.names = row.names(m_adrenal))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)
dim(p_data)
dim(f_data)
dim(expr_matrix)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

#构建CDS对象
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
####估计size factor和离散度
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)



###过滤低质量的细胞和基因
cds <- detectGenes(cds, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed
#print(head(fData(cds)))#此时有个基因
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 44)) 
##################选择定义过程的基因，三种都是无监督分析方法，细胞发育轨迹生成完全不受人工干预
#####用seurat筛选出高变基因来做
m_adrenal= FindVariableFeatures(m_adrenal)
express_genes1 = VariableFeatures(m_adrenal)
##各个细胞类型marker基因来做
deg.cluster <- FindAllMarkers(m_adrenal)
express_genes2 <- subset(deg.cluster,p_val_adj<0.05)$gene
##monocle计算的基因
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id

##
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds,max_components = 2,
                       reduction_method = 'DDRTree')
cds <- orderCells(cds)

###输出文件
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
cols <- viridis_pal(option = "D")(4)
pdf("all-1.pdf",width=8,height=6)
p=plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points = FALSE,cell_size = 1.0)+ scale_color_gradientn(colors=viridis_pal(option = "D")(4))
print(p)
dev.off()
pdf("all-2.pdf",width=8,height=6)
p=plot_cell_trajectory(cds, color_by = "celltype",show_branch_points = FALSE,cell_size = 1.0)+scale_color_manual(values= c('#F2CD5C',"#F8766D","#C77CFF","#BA1C96"))
print(p)
dev.off()
pdf("all-3.pdf",width=10,height=6)
p=plot_cell_trajectory(cds, color_by = "celltype",show_branch_points = FALSE,cell_size = 1.0)+scale_color_manual(values= c('#F2CD5C',"#F8766D","#C77CFF","#BA1C96"))+facet_wrap(~celltype, nrow = 2)
print(p)
dev.off()
pdf("all-14.pdf",width=8,height=6)
p=plot_cell_trajectory(cds, color_by = "celltype",show_branch_points = FALSE,cell_size = 0.8)+scale_color_manual(values= c('#F2CD5C',"#F8766D","#C77CFF","#BA1C96"))+facet_wrap(~group, nrow = 1)
print(p)
dev.off()
pdf("all-15.pdf",width=4,height=6)
p=plot_cell_trajectory(cds, color_by = "group",show_branch_points = FALSE,cell_size = 0.8)+scale_color_manual(values= rev(getPalette(2)))+facet_wrap(~group, nrow = 2)
print(p)
dev.off()

cds$group <- factor(cds$group,levels = c("Young_Male","Old_Male"))
pdf("all-ct_group.pdf",width=6,height=6)
p=plot_cell_trajectory(cds, color_by = "celltype",show_branch_points = FALSE,cell_size = 1.0)+scale_color_manual(values= c('#F2CD5C',"#F8766D","#C77CFF","#BA1C96"))+facet_grid(vars(celltype), vars(group))+theme_bw()
print(p)
dev.off()

pdf("all-6.pdf",width=4,height=6)
p=plot_cell_trajectory(cds, color_by = "seurat_clusters",show_branch_points = FALSE,cell_size = 1.0)+scale_color_manual(values= getPalette(25))
print(p)
dev.off()
pdf("all-7.pdf",width=10,height=6)
p=plot_cell_trajectory(cds, color_by = "seurat_clusters",show_branch_points = FALSE,cell_size = 1.0)+scale_color_manual(values= getPalette(25))+facet_wrap(~group, nrow = 2)
print(p)
dev.off()



#保存结果
pdata <- Biobase::pData(cds)

write.csv(pData(cds), "pseudotime.csv")
save(cds, file = "1100cell_cds.rda")




#####
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(ggplot2)
library(future)
library(RColorBrewer)
library(tidyverse)
library(htmlwidgets)  
library(ggpubr)      
library(colorRamps)   
library(magrittr)  
library(monocle)
library(pheatmap)
library(viridis)
set.seed(1998)
setwd("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/13_pesudotime/ZG_ZR")
cds1 <- load("1100cell_cds.rda")
ordergene <- ordering_genes
#指定基因的可视化
##选择前4个top基因并将其对象取出
keygenes <- head(ordergene,4)
cds_subset <- cds[keygenes,]
##可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3
ggsave("Genes_pseudotimeplot.pdf", plot = plotc, width = 16, height = 8)
#指定基因
s.genes <- intersect(c("SHH","GLI1","WNT4", "MC2R", "PKA", "MRAP","SF1","CYP11B2","CYP11B1","CYB5","CC3","DLK1","ERK1","ERK2","CYB5A"),ordergene)
p1 <- plot_genes_jitter(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p2 <- plot_genes_violin(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "celltype")
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot.pdf", plot = plotc, width = 16, height = 8)

###cluster相关基因
s.genes <- intersect(c("CACNA1D","BCL2","KCNK3","NR4A1","PRKD1","ZNRF3"),rownames(cds))
p1 <- plot_genes_jitter(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p2 <- plot_genes_violin(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "celltype")
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot_cluster2.pdf", plot = plotc, width = 16, height = 8)

###cluster相关基因
s.genes <- intersect(c("LDLR","NR1H4","SCARB1","CYP17A1","PBX1","STAR"),rownames(cds))
p1 <- plot_genes_jitter(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p2 <- plot_genes_violin(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "celltype")
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot_cluster3.pdf", plot = plotc, width = 16, height = 8)

###cluster相关基因
s.genes <- intersect(c("MAP3K5","SOD1","SOD2","SIRT1","NCOA7","ARNT"),rownames(cds))
p1 <- plot_genes_jitter(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p2 <- plot_genes_violin(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "celltype")
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot_cluster4.pdf", plot = plotc, width = 16, height = 8)

###cluster相关基因
s.genes <- intersect(c("SESN2","CYP51A1","SULT2A1","CYB5R3","HMGCR","MVD"),rownames(cds))
p1 <- plot_genes_jitter(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p2 <- plot_genes_violin(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "celltype")
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot_cluster5.pdf", plot = plotc, width = 16, height = 8)

####marker genes
p1 <- plot_genes_in_pseudotime(cds[c("LEF1","SHH","WNT4"),], color_by = "celltype",ncol = 3)
p2 <- plot_genes_in_pseudotime(cds[c("CACNA1D","KCNK3","NR4A1"),], color_by = "celltype",ncol = 3)
p3 <- plot_genes_in_pseudotime(cds[c("CYP17A1","APOE","APOC1"),], color_by = "celltype",ncol = 3)
p4 <- plot_genes_in_pseudotime(cds[c("MAP3K5","NCOA7","SIRT1"),], color_by = "celltype",ncol = 3)
p5 <- plot_genes_in_pseudotime(cds[c("HMGCR","CYP51A1","SULT2A1"),], color_by = "celltype",ncol = 3)
plots <- p1/p2/p3/p4/p5
ggsave("Genes_Jitterplot_cluster1-5.pdf", plot = plots, width = 8, height = 12)

###marker_gene####
marker.genes <- c("VSNL1","CYP11A1","CYP11B2","CYP17A1","CYP51A1","SULT2A1")
p3 <- plot_genes_in_pseudotime(cds[marker.genes,], color_by = "celltype",ncol = 3)
ggsave("marker_Genes_Jitterplot.pdf", plot = p3, width = 3, height = 8)



###寻找拟时差异基因（qvalue体现基因与拟时的密切程度
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 44))
expressed_genes <- diff_test_res$gene_short_name
Time_diff <- differentialGeneTest(cds[expressed_genes,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
#Time_diff <- Time_diff[,c(4,3,2,5,1)] #把gene放前面，也可以不改
write.csv(Time_diff, "Time_diff_all.csv", row.names = F)
Time_genes <- row.names(Time_diff)[order(Time_diff$qval)][1:1000]
write.csv(Time_genes, "Time_diff.csv", row.names = F)
####添加列坐标####
Binner <- function(cells_subset){
  df <- data.frame(pData(cds))
  df <- df[,c("Pseudotime", "celltype")]
  df <- df[order(df$Pseudotime, decreasing = F),]
  len <- length(df$Pseudotime)
  bin<-round(len/4400)
 celltype <- c()
  value <- c()
  for(i in 0:4399){
    if(i < 4399){
      start <- 1+(bin*i)
      stop <- bin+(bin*i)
      getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
      value <- getmode(as.vector(df$celltype[c(start:stop)]))
      celltype <- c(celltype, value)
    }
    else{
      celltype <- c(celltype, value)
    }
  }
  return(as.data.frame(celltype))
}
bin <- Binner(cds)

p3 = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=3, show_rownames=F, return_heatmap=T)
p4 = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=F, return_heatmap=T)
p5 = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=5, show_rownames=F, return_heatmap=T)
p6 = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=6, show_rownames=F, return_heatmap=T)
p7 = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=7, show_rownames=F, return_heatmap=T)

ggsave("Time_heatmapAll3.pdf", p3, width = 4, height = 10)
ggsave("Time_heatmapAll4.pdf", p4, width = 4, height = 10)
ggsave("Time_heatmapAll5.pdf", p5, width = 4, height = 10)
ggsave("Time_heatmapAll6.pdf", p6, width = 4, height = 10)
ggsave("Time_heatmapAll7.pdf", p7, width = 4, height = 10)


p5$tree_row # Call: # hclust(d = d, method = method) 
clusters <- cutree(p5$tree_row, k = 5) 
clustering <- data.frame(clusters) 
write.csv(clustering, "Time_clustering_all.csv", row.names = T)



####重新封装函数####
plot_pse <- function(cds_subset, cluster_rows = TRUE, hclust_method = "ward.D2", 
    num_clusters = 5, hmcols = NULL, add_annotation_row = NULL, 
    add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE, 
    norm_method = c("log", "vstExprs"), scale_max = 3, scale_min = -3, 
    trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE, 
    cores = 1) 
{
    num_clusters <- min(num_clusters, nrow(cds_subset))
    pseudocount <- 1
    newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
        max(pData(cds_subset)$Pseudotime), length.out = 4400))
    m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula, 
        relative_expr = T, new_data = newdata)
    m = m[!apply(m, 1, sum) == 0, ]
    norm_method <- match.arg(norm_method)
    if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == 
        FALSE) {
        m = vstExprs(cds_subset, expr_matrix = m)
    }
    else if (norm_method == "log") {
        m = log10(m + pseudocount)
    }
    m = m[!apply(m, 1, sd) == 0, ]
    m = Matrix::t(scale(Matrix::t(m), center = TRUE))
    m = m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] = 0
    m[m > scale_max] = scale_max
    m[m < scale_min] = scale_min
    heatmap_matrix <- m
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    if (is.null(hmcols)) {
        bks <- seq(-3.1, 3.1, by = 0.1)
        hmcols <- blue2green2red(length(bks) - 1)
    }
    else {
        bks <- seq(-3.1, 3.1, length.out = length(hmcols))
    }
    ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
        cluster_rows = cluster_rows, show_rownames = F, show_colnames = F, 
        clustering_distance_rows = row_dist, clustering_method = hclust_method, 
        cutree_rows = num_clusters, silent = TRUE, filename = NA, 
        breaks = bks, border_color = NA, color = hmcols)
    if (cluster_rows) {
        annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
            num_clusters)))
    }
    else {
        annotation_row <- NULL
    }
    if (!is.null(add_annotation_row)) {
        old_colnames_length <- ncol(annotation_row)
        annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
            ])
        colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
    }
    if (!is.null(add_annotation_col)) {
        if (nrow(add_annotation_col) != 4400) {
            stop("add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!")
        }
        annotation_col <- add_annotation_col
    }
    else {
        annotation_col <- NA
    }
    if (use_gene_short_name == TRUE) {
        if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
            feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 
                "gene_short_name"])
            feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
            row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 
                "gene_short_name"])
            row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
        }
        else {
            feature_label <- row.names(heatmap_matrix)
            row_ann_labels <- row.names(annotation_row)
        }
    }
    else {
        feature_label <- row.names(heatmap_matrix)
        if (!is.null(annotation_row)) 
            row_ann_labels <- row.names(annotation_row)
    }
    row.names(heatmap_matrix) <- feature_label
    if (!is.null(annotation_row)) 
        row.names(annotation_row) <- row_ann_labels
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
        cluster_rows = cluster_rows, show_rownames = show_rownames, 
        show_colnames = F, clustering_distance_rows = row_dist, 
        clustering_method = hclust_method, cutree_rows = num_clusters, 
        annotation_row = annotation_row, annotation_col = annotation_col, 
        treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols, 
        border_color = NA, silent = TRUE, filename = NA)
    grid::grid.rect(gp = grid::gpar("fill", col = NA))
    grid::grid.draw(ph_res$gtable)
    if (return_heatmap) {
        return(ph_res)
    }
}

p=plot_pse(cds[Time_genes,], num_clusters=5, show_rownames=F,add_annotation_col = bin,
                                                                     return_heatmap=T) 
ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)
######拟合成曲线####
#####曲线图######
time_select_gene <- as.data.frame(Time_genes)
colnames(time_select_gene) <- "gene"
my_cds <-  cds
length.out <- nrow(pData(my_cds))
times.df <- seq(0.1,100,length.out = length.out)
pData(my_cds)$age_time <- times.df


hclust_method = "ward.D2"
num_clusters = 5 #### Set the number of clusters
scale_max = 3
scale_min = -3
cores = 20

num_clusters <- min(num_clusters, nrow(my_cds))
pseudocount <- 1
newdata <- data.frame(Pseudotime = seq(min(pData(my_cds)$Pseudotime), 
        max(pData(my_cds)$Pseudotime), length.out =length.out))
### Fit smooth spline curves and return the response matrix
smooth_mat <- genSmoothCurves(my_cds, cores = cores, trend_formula = "~sm.ns(Pseudotime,df=3)",
                     relative_expr = T, new_data = newdata)
smooth_mat = smooth_mat[!apply(smooth_mat, 1, sum) == 0, ]
smooth_mat = vstExprs(my_cds, expr_matrix = smooth_mat)
smooth_mat = smooth_mat[!apply(smooth_mat, 1, sd) == 0, ]
smooth_mat = Matrix::t(scale(Matrix::t(smooth_mat), center = TRUE))
smooth_mat = smooth_mat[is.na(row.names(smooth_mat)) == FALSE, ]
smooth_mat[is.nan(smooth_mat)] = 0
smooth_mat[smooth_mat > scale_max] = scale_max
smooth_mat[smooth_mat < scale_min] = scale_min

heatmap_matrix <- smooth_mat

library(reshape2)
for (i in c(1:5)){
 tmp <- subset(clustering,clusters==i)
 gene <- rownames(tmp)
 newdata <- heatmap_matrix[gene,]
 df <- melt(newdata)
 colnames(df) <- c("gene","order",'expression')
 
 p <- ggplot(df,aes(x=order, y=expression))+
  geom_smooth(formula = y~poly(x,15),method = lm,size = 1,alpha = 1,color="#EE3B3B",se=FALSE)+ 
  #stat_smooth(aes(group=gene),formula = y~poly(x,15),method= lm, size = 0.02,linetype=2, alpha = 0.05,color="#6495ED",se=FALSE)+
  #scale_x_continuous(limits=c(1,1500), breaks=seq(250, 1250, 500),label = c("Ctrl","A_TB","F_TB"))+
  #geom_vline(xintercept=c(1,500,1000,1500),lty=4,col="#000000",lwd=0.6)+
  theme(panel.grid=element_blank())+
  theme_bw()
 pdf(paste0("/dellstorage09/quj_lab/wangxuebao/01_results/01_adrenal/13_pesudotime/ZG_ZR/cluster_",i,".pdf"),height = 3,width = 3.5)
 print(p)
 dev.off()
}























#####不同过滤方法####

###可以将四种细胞类型拟合成一条线
cds <- detectGenes(cds, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 44)) #过滤掉在小于44个细胞中表达的基因，还剩个基因。
diff_test_res <- differentialGeneTest(cds[expressed_genes,],
              fullModelFormulaStr = "~celltype")
#diff_test_res <- subset(diff_test_res,str_sub(gene_short_name,1,7) != "ENSMFAG" )
ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:2000]




