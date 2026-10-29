library(Seurat)
library(harmony)
library(draw)
library(GeneNMF)
library(scplotter)
library(COSG)
library(tibble)
library(ggplot2)
library(RColorBrewer)
setwd("D:/aging")
dir='data/GSE213516' 
samples=list.files(dir)
samples 
sceList_lap = lapply(samples,function(pro){ 
  #pro=samples[1] 
  print(pro)
  sce=CreateSeuratObject(counts =  Read10X(file.path(dir,pro) ) ,
                         project =  pro ,
                         min.cells = 5,
                         min.features = 300,)
  return(sce)
})
names(sceList_lap) =  samples
sce <- merge(x = sceList_lap[[1]], 
                   y= sceList_lap[ -1 ] ,
                   add.cell.ids =  samples) 
################################################
setwd("D:/aging/1.sc")
sce$mitoPercent <- PercentageFeatureSet(sce, pattern = "^MT-")
pdf("Feature_Count_VlnPlot.pdf",width =23,height = 10)
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA","mitoPercent"),pt.size = 0, ncol = 3)
VlnPlot(sce, features = c("nFeature_RNA"),pt.size = 0, ncol = 1)
VlnPlot(sce, features = c("nCount_RNA"),pt.size = 0, ncol = 1)
VlnPlot(sce, features = c("mitoPercent"),pt.size = 0, ncol = 1)
dev.off()
sce <- subset(sce,subset=nCount_RNA>800 &
                nFeature_RNA>200 &nFeature_RNA<4000&
                mitoPercent<20)
pdf(file = "FeatureScatter.pdf",width =23,height = 10)
VlnPlot(sce, features = c("nFeature_RNA"),pt.size = 0, ncol = 1)
VlnPlot(sce, features = c("nCount_RNA"),pt.size = 0, ncol = 1)
VlnPlot(sce, features = c("mitoPercent"),pt.size = 0, ncol = 1)
dev.off()
#######################################
sce <- JoinLayers(sce)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000)
sce <- ScaleData(object = sce)
sce = RunPCA(sce, npcs = 50, verbose = FALSE)
sce <- RunHarmony(sce,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
sce<-FindNeighbors(sce,reduction="harmony",
                   dims=1:50)
sce<-FindClusters(sce,resolution=seq(from=0.5,to=1.5,by=0.1))#计算多个分辨率library(clustree)clustree(scRNAc)
scplotter::ClustreePlot(sce, prefix = "RNA_snn_res.") 
drawExport("clustertree.pdf",height = 16,width = 16)
sce<-RunUMAP(sce,dims=1:50,
             reduction="harmony")
sce<-RunTSNE(sce,dims=1:50,
             reduction="harmony")
save(sce,file = "sce.rda")
load(file  = "sce.rda")
CellDimPlot(sce,
            group_by = "RNA_snn_res.1", reduction = "umap",
            label = TRUE, label_insitu = TRUE, label_repel = TRUE,
            label_segment_color = "red")
Idents(sce) <- sce$RNA_snn_res.1
Idents(sce) <- factor(Idents(sce), levels =  c( 0:31))
genes_to_check=list(Tcell=c("CD3D"),
                    B=c("CD79A"),
                    NK=c("NKG7"),
                    #Monocyte=c("CD14","FCGR3A","FCN1"),
                    Platelet=c("PPBP","PF4","GP1BA"),
                    Monocytes = c('IFITM3','FCN1'),
                    #DC=c("IL3RA"),
                    #HSPC=c("CD34","CD38","PROM1"),
                    #Myeloid=c("LYZ","CST3","CD14","CD68"),
                    #Neutrophil=c("S100A8","S100A9","CSF3R"),
                    pDCs=c('IL3RA','CLEC4C','LILRB4'),    
                    'cDC1/cDC2' = c('CD1C','CD1E','FCER1A','XCR1','CLEC9A','IDO1'),
                    #Neutrophils=c('CSF3R','S100A8','S100A9'),
                    Erythroblast=c( "HBA1", "HBB"))
DotPlot(sce,features=genes_to_check,cols="RdYlBu")+
  RotatedAxis()
COSG_markers<-cosg(
  sce,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=50)
a <- COSG_markers$names
celltype=data.frame(ClusterID=0:31,
                    celltype= 0:31) 
celltype[celltype$ClusterID %in% c(0,5,18,20,27),2]='NK'
celltype[celltype$ClusterID %in% c(1,7,11,17),2]='Monocyte'   
celltype[celltype$ClusterID %in% c(2,3,4,8,9,12,13,14,19,26,30),2]='T'  
celltype[celltype$ClusterID %in% c(6,10,22,28,29),2]='B'
celltype[celltype$ClusterID %in% c(15,16,23),2]='Platelet'
celltype[celltype$ClusterID %in% c(21),2]='cDC1/2'
celltype[celltype$ClusterID %in% c(24),2]='Erythroblast'
celltype[celltype$ClusterID %in% c(25),2]='pDCs'
celltype[celltype$ClusterID %in% c(31),2]='HSPC'
sce@meta.data$celltype <- sapply(sce@meta.data$RNA_snn_res.1,function(x){celltype[x,2]})
Idents(sce) <- sce$celltype
CellDimPlot(sce, group_by = "celltype", reduction = "umap")
drawExport("anno.pdf",width = 10,height = 5)
save(sce,file = "sce.rda")
####################################################################ROGUE
library(ROGUE)
library(ggplot2)
library(ggridges)
expr <- GetAssayData(JoinLayers(sce))
options(future.globals.maxSize = 100 * 1024^3)
ent.res <- SE_fun(expr)
head(ent.res)
SEplot(ent.res)
ggsave("SEplot.pdf")
rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value
rogue.res <- rogue(expr, labels = sce$celltype, samples = sce$orig.ident, platform = "UMI",
                   span = 0.6)
rogue.res
save(rogue.res,file = "rogue.rda")
library(tidyr)
long_df <- rogue.res %>%
  rownames_to_column(var = "Sample") %>%  
  pivot_longer(
    cols = -Sample,                     
    names_to = "CellType",              
    values_to = "Rogue Score"            
  )
mycol3 <- c(brewer.pal(n=9,name="Set1"),brewer.pal(n=8,name="Set2")
            ,brewer.pal(n=12,name="Set3"))
ggplot(data = long_df, aes(x = `Rogue Score`, y = CellType, fill = CellType)) +
  geom_density_ridges(alpha = 0.8,
                      color= 'white',
                      rel_min_height= 0.01, #尾部修剪，数值越大修剪程度越高
                      scale= 1.8, #山脊重叠程度调整，scale = 1时刚好触及基线，数值越大重叠度越高
                      quantile_lines= TRUE, #显示分位数线
                      quantiles= 2 #仅显示中位数线
  ) +
  scale_fill_manual(values = mycol3) +
  theme_classic()
ggsave("rogue_boxplot.pdf")
###########################################################
sample_table <- as.data.frame(table(sce$celltype,sce$orig.ident))
#定义列名
names(sample_table ) <- c("celltype","Group","CellNumber")
#指定颜色
colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C")
#横轴为样本,纵轴为细胞比例
ggplot(sample_table,aes(x=Group,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values = colors) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)  )+labs(y="Percentage")+RotatedAxis()
save(sce,file = "sce.rda")
drawExport("Percent.pdf")
#####################################################################
features <- c(
  "CCR7","TCF7","LEF1","SELL", #Naive
  "GZMK","IFNG","GZMH","GNLY","PRF1","NKG7","GZMA","GZMB", #Cytotoxicity
  "LAG3","PDCD1","CTLA4","TIGIT","BTLA","HAVCR2",#Inhibitory
  "IL2RA","FOXP3","IKZF2",#Treg
  "MKI67","CDK1","STMN1"#Proliferating    
)
rows_data <- data.frame(
  features = features,
  group = c(
    "Naive", "Naive", "Naive", "Naive", "Cytotoxicity","Cytotoxicity","Cytotoxicity","Cytotoxicity","Cytotoxicity",
    "Cytotoxicity","Cytotoxicity","Cytotoxicity", "Inhibitory", "Inhibitory", "Inhibitory","Inhibitory","Inhibitory","Inhibitory",
    "Treg","Treg","Treg","Proliferating","Proliferating","Proliferating")
  #TF = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE,TRUE, TRUE, TRUE),
  #CSPA = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE,FALSE, FALSE, FALSE)
)

FeatureStatPlot(sce, features = features, ident = "cell_type", cell_type = "dot",
                plot_type = "heatmap", name = "Expression Level", dot_size = function(x) sum(x > 0) / length(x),
                dot_size_name = "Percent Expressed", add_bg = TRUE,
                rows_data = rows_data, 
                show_row_names = TRUE, 
                rows_split_by = "group", 
                cluster_rows = FALSE,
                palette = "YlOrRd",
                #row_annotation_side = "right"
                )
drawExport("FeatureStatPlot.pdf",height = 10,width = 10)
############################################################
meta <- sce@meta.data
meta$age <- substring(meta$orig.ident,2,3)
meta$age <- as.numeric(meta$age)
sce@meta.data <- meta
save(sce,file = "sce.rda")
################################################组织偏好性RO/E
library(Startrac)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
data <- sce@meta.data[,c(1,17,18)]
data$AgeGroup <- 0
data[which(data$age<=45),4] <- "Young"
data[which(data$age>45&data$age<=60),4] <- "Middle"
data[which(data$age<=80&data$age>60),4] <- "Old"
R_oe <- calTissueDist(data,          
                      byPatient = F,   
                      colname.cluster = "celltype",   
                      #colname.patient ="orig.ident",# 不同样本
                      colname.tissue = "AgeGroup",       
                      method = "chisq",             
                      min.rowSum = 0) 
R_oe
Roe <- as.matrix(R_oe)
col_fun <- colorRamp2(c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)),  
                                     c("#91d1c2b2", "white", "#dc0000b2"))
plotdata <- R_oe[,c("Young","Middle","Old")]
Heatmap(as.matrix(plotdata),  
        show_heatmap_legend = TRUE,  
        cluster_rows = F,       
        cluster_columns = F,       
        row_names_side = 'right',     
        show_column_names = TRUE,        
        show_row_names = TRUE,     
        col = col_fun,       
        row_names_gp = gpar(fontsize = 10), 
        column_names_gp = gpar(fontsize = 10),  
        heatmap_legend_param = list(     
          title = "Ro/e Index",  # 自定义图注名称 
          at = seq(0.5, 1.5, by = 0.5), # 例刻度的位置/自己的数据必须修改一下！ 
          labels = seq(0.5, 1.5, by = 0.5) # 每个刻度的标签/自己的数据必须修改一下！ 
          )
          )
drawExport("Roe.pdf",height = 6,width = 8)
 ########################################################
library(miloR)
scRNA_pre<-as.SingleCellExperiment(sce)
#构建miloR对象
scRNA_pre<-Milo(scRNA_pre)
scRNA_pre
#######################################################
library(dplyr)
library(MetBrewer)
library(ggplot2)
library(ggalluvial)
meta <- sce@meta.data %>% 
  dplyr::select(celltype, age) %>% 
  dplyr::count(celltype, age, name = "number")  
meta$group <- 0
meta[which(meta$age<=45),4] <- 0
meta[which(meta$age>45&meta$age<=60),4] <- 1
meta[which(meta$age<=80&meta$age>60),4] <- 2
result <- meta %>%
  group_by(group) %>%  
  mutate(total = sum(number)) %>%  
  mutate(prop = number / total) %>%  
  ungroup() 


colors <- met.brewer("Cross")

ggplot(result, aes(x = group, y=prop, fill = celltype, stratum = celltype, alluvium = celltype))+
  geom_stratum(width = 0.5, color='white')+
  geom_alluvium(alpha = 0.5, width = 0.5, curve_type = "linear") +
  #theme_void(base_size = 25) +
  theme(legend.position = "top") +
  scale_fill_manual(values = colors)+
  labs(
    x = "age",  # 设置 x 轴标签（可根据需求改名，比如“隐藏状态”）
    title = "stratum"     # 添加大标题
  )
drawExport("stratum.pdf",height = 5,width = 5)
###########################################################各细胞类型差异基因
library(dplyr)
library(openxlsx)

Idents(sce) <- sce$celltype

# 每个细胞类型 vs 其余所有细胞
markers <- FindAllMarkers(
  sce,
  assay           = "RNA",
  only.pos        = TRUE,   # 只保留上调基因
  min.pct         = 0.25,   # 至少在 25% 的细胞中表达
  logfc.threshold = 0.25,
  test.use        = "wilcox", # 已安装 presto，速度较快
  verbose         = TRUE
)

# 显著上调基因：校正 p < 0.05 且 avg_log2FC > 0.5
markers_sig <- markers |>
  filter(p_val_adj < 0.05, avg_log2FC > 0.5)

# 每个细胞类型差异最显著（按 log2FC 排序）的前 50 个基因
markers_top50 <- markers_sig |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 50) |>
  ungroup()

save(markers, markers_sig, markers_top50, file = "markers.rda")
write.csv(markers, "markers_all.csv", row.names = FALSE)

# 输出 excel：汇总表 + 每个细胞类型一个 sheet
safe_sheet <- function(x) substr(gsub("[^A-Za-z0-9_.-]", "_", x), 1, 31)
wb <- createWorkbook()

addWorksheet(wb, "All_markers")
writeData(wb, "All_markers", markers)
freezePane(wb, "All_markers", firstRow = TRUE)
setColWidths(wb, "All_markers", cols = 1:ncol(markers), widths = "auto")

addWorksheet(wb, "Top50_per_celltype")
writeData(wb, "Top50_per_celltype", markers_top50)
freezePane(wb, "Top50_per_celltype", firstRow = TRUE)

for (ct in levels(Idents(sce))) {
  sheet <- safe_sheet(paste0("Sig_", ct))
  addWorksheet(wb, sheet)
  writeData(wb, sheet, filter(markers_sig, cluster == ct))
  freezePane(wb, sheet, firstRow = TRUE)
}
saveWorkbook(wb, "markers_by_celltype.xlsx", overwrite = TRUE)