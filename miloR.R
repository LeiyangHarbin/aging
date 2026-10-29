#########################################################
library(miloR)
library(scater)
library(patchwork)
library(ggplot2)
library(ggbeeswarm)
library(draw)
library(dplyr)
load("D:/aging/2.T/Tsce.rda")
setwd("D:/aging/2.T")
Idents(Tsce) <- Tsce$AgeGroup
############################################
sce <- subset(Tsce,idents=c("Young","Middle"))
sce <- as.SingleCellExperiment(sce)
sce_milo <- Milo(sce)
reducedDim(sce_milo,"UMAP") <- reducedDim(sce,"UMAP")
sce_milo <- buildGraph(sce_milo, k = 50, d = 50,reduced.dim = "PCA")
sce_milo <- makeNhoods(sce_milo, prop = 0.2, k = 50, d = 50, 
                       refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(sce_milo)
sce_milo <- countCells(sce_milo, 
                       meta.data = as.data.frame(colData(sce_milo)),
                       sample="orig.ident") # sample代表样本
# colData(scRNA) 相当于代表了metadata的信息
# sample代表了样本数
head(nhoodCounts(sce_milo))
# 这个步骤能够给scRNA数据中增加一个n*m的matrix
# n代表不同的neighborhood数量,m代表实验样本
traj_design <- data.frame(colData(sce_milo))[,c("orig.ident", "AgeGroup")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$orig.ident
#重新排序行名以匹配nhoodCounts（milo）的列
traj_design <- traj_design[colnames(nhoodCounts(sce_milo)), , drop=FALSE]
traj_design$AgeGroup <- factor(traj_design$AgeGroup, level = c("Young", "Middle"))
traj_design
# Milo使用cydar引入的Spatial FDR校正的变体，该校正考虑了邻域之间的重叠。具体来说，每个假设检验的P值都通过k个最近邻距离的倒数进行加权。要使用这个统计量，首先需要在Milo对象中存储最近邻之间的距离。
sce_milo <- calcNhoodDistance(sce_milo, d=50, reduced.dim = "PCA")
da_results <- testNhoods(sce_milo,design = ~ AgeGroup, design.df = traj_design)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 
sce_milo <- buildNhoodGraph(sce_milo)
## Plot single-cell UMAP
umap_pl <- plotReducedDim(sce_milo, dimred = "UMAP", 
                          colour_by="AgeGroup", text_by = "celltype", 
                          text_size = 3, point_size=0.5) + guides(fill="none")
umap_pl
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(sce_milo, da_results,
                                layout="UMAP",alpha=0.9)  #alpha默认0.1
umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
drawExport("young-middle-umap.pdf")

da_results <- annotateNhoods(sce_milo, 
                             da_results, 
                             coldata_col = "celltype")
sig_nhoods <- sum(da_results$SpatialFDR < 0.05, na.rm = TRUE)

da_results$celltype <- factor(da_results$celltype,
                                     levels= c("CD4-CCR6", "CD4-CCR7","CD4-FOXP3", "CD4-GZMK", "CD4-TCF7", "CD8-CCR7","CD8-CXCR1","CD8-FGFBP2","CD8-GZMK","CD8-HAVCR2","Exhausted T","Naive_T","NKT-GNLY"))
plotDAbeeswarm(da_results, group.by = "celltype", alpha=1)+geom_boxplot(aes(group = celltype), width = 0.3, outlier.shape = NA, fill = NA)+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  geom_hline(yintercept=0, linetype="dashed")
drawExport("young-middle-DAbeeswarm.pdf")
#################################################################
sce <- subset(Tsce,idents=c("Young","Old"))
sce <- as.SingleCellExperiment(sce)
sce_milo <- Milo(sce)
reducedDim(sce_milo,"UMAP") <- reducedDim(sce,"UMAP")
sce_milo <- buildGraph(sce_milo, k = 50, d = 50,reduced.dim = "PCA")
sce_milo <- makeNhoods(sce_milo, prop = 0.2, k = 50, d = 50, 
                       refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(sce_milo)
sce_milo <- countCells(sce_milo, 
                       meta.data = as.data.frame(colData(sce_milo)),
                       sample="orig.ident") # sample代表样本
# colData(scRNA) 相当于代表了metadata的信息
# sample代表了样本数
head(nhoodCounts(sce_milo))
# 这个步骤能够给scRNA数据中增加一个n*m的matrix
# n代表不同的neighborhood数量,m代表实验样本
traj_design <- data.frame(colData(sce_milo))[,c("orig.ident", "AgeGroup")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$orig.ident
#重新排序行名以匹配nhoodCounts（milo）的列
traj_design <- traj_design[colnames(nhoodCounts(sce_milo)), , drop=FALSE]
traj_design$AgeGroup <- factor(traj_design$AgeGroup, level = c("Young","Old"))
traj_design
# Milo使用cydar引入的Spatial FDR校正的变体，该校正考虑了邻域之间的重叠。具体来说，每个假设检验的P值都通过k个最近邻距离的倒数进行加权。要使用这个统计量，首先需要在Milo对象中存储最近邻之间的距离。
sce_milo <- calcNhoodDistance(sce_milo, d=50, reduced.dim = "PCA")
da_results <- testNhoods(sce_milo,design = ~ AgeGroup, design.df = traj_design)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 
sce_milo <- buildNhoodGraph(sce_milo)
## Plot single-cell UMAP
umap_pl <- plotReducedDim(sce_milo, dimred = "UMAP", 
                          colour_by="AgeGroup", text_by = "celltype", 
                          text_size = 3, point_size=0.5) + guides(fill="none")
umap_pl
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(sce_milo, da_results,
                                layout="UMAP",alpha=0.9)  #alpha默认0.1
umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
drawExport("young-old-umap.pdf")
da_results <- annotateNhoods(sce_milo, 
                             da_results, 
                             coldata_col = "celltype")
da_results$celltype <- factor(da_results$celltype,
                              levels= c("CD4-CCR6", "CD4-CCR7","CD4-FOXP3", "CD4-GZMK", "CD4-TCF7", "CD8-CCR7","CD8-CXCR1","CD8-FGFBP2","CD8-GZMK","CD8-HAVCR2","Exhausted T","Naive_T","NKT-GNLY"))
plotDAbeeswarm(da_results, group.by = "celltype", alpha=1)+geom_boxplot(aes(group = celltype), width = 0.3, outlier.shape = NA, fill = NA)+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  geom_hline(yintercept=0, linetype="dashed")
drawExport("young-old-DAbeeswarm.pdf")
#################################################################
sce <- subset(Tsce,idents=c("Middle","Old"))
sce <- as.SingleCellExperiment(sce)
sce_milo <- Milo(sce)
reducedDim(sce_milo,"UMAP") <- reducedDim(sce,"UMAP")
sce_milo <- buildGraph(sce_milo, k = 50, d = 50,reduced.dim = "PCA")
sce_milo <- makeNhoods(sce_milo, prop = 0.2, k = 50, d = 50, 
                       refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(sce_milo)
sce_milo <- countCells(sce_milo, 
                       meta.data = as.data.frame(colData(sce_milo)),
                       sample="orig.ident") # sample代表样本
# colData(scRNA) 相当于代表了metadata的信息
# sample代表了样本数
head(nhoodCounts(sce_milo))
# 这个步骤能够给scRNA数据中增加一个n*m的matrix
# n代表不同的neighborhood数量,m代表实验样本
traj_design <- data.frame(colData(sce_milo))[,c("orig.ident", "AgeGroup")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$orig.ident
#重新排序行名以匹配nhoodCounts（milo）的列
traj_design <- traj_design[colnames(nhoodCounts(sce_milo)), , drop=FALSE]
traj_design$AgeGroup <- factor(traj_design$AgeGroup, level = c("Middle","Old"))
traj_design
# Milo使用cydar引入的Spatial FDR校正的变体，该校正考虑了邻域之间的重叠。具体来说，每个假设检验的P值都通过k个最近邻距离的倒数进行加权。要使用这个统计量，首先需要在Milo对象中存储最近邻之间的距离。
sce_milo <- calcNhoodDistance(sce_milo, d=50, reduced.dim = "PCA")
da_results <- testNhoods(sce_milo,design = ~ AgeGroup, design.df = traj_design)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 
sce_milo <- buildNhoodGraph(sce_milo)
## Plot single-cell UMAP
umap_pl <- plotReducedDim(sce_milo, dimred = "UMAP", 
                          colour_by="AgeGroup", text_by = "celltype", 
                          text_size = 3, point_size=0.5) + guides(fill="none")
umap_pl
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(sce_milo, da_results,
                                layout="UMAP",alpha=0.9)  #alpha默认0.1
umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
drawExport("Middle-old-umap.pdf")
da_results <- annotateNhoods(sce_milo, 
                             da_results, 
                             coldata_col = "celltype")
da_results$celltype <- factor(da_results$celltype,
                              levels= c("CD4-CCR6", "CD4-CCR7","CD4-FOXP3", "CD4-GZMK", "CD4-TCF7", "CD8-CCR7","CD8-CXCR1","CD8-FGFBP2","CD8-GZMK","CD8-HAVCR2","Exhausted T","Naive_T","NKT-GNLY"))
plotDAbeeswarm(da_results, group.by = "celltype", alpha=1)+geom_boxplot(aes(group = celltype), width = 0.3, outlier.shape = NA, fill = NA)+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  geom_hline(yintercept=0, linetype="dashed")
drawExport("middle-old-DAbeeswarm.pdf")