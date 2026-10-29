library(msigdbr)
library(GSVA)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(readr)
setwd("D:/aging/13.bio/1.gsva")
pd <- read_tsv("D:/aging/13.bio/4.pd/TCGA-LAML.clinical.tsv")
load("d:/aging/12.1getx/TCGAexp.rda")
HMM_res <- read.csv("D:/aging/12.1getx/data_with_states.csv")
HMM_states <- HMM_res[,c("X","hidden_state")]

genelist <- readGMT("h.all.v2025.1.Hs.symbols.gmt")
gsvapar <- gsvaParam(as.matrix(exp), genelist, maxDiff=TRUE)

gsva_es <- gsva(gsvapar)

anno <- merge(pd,HMM_states,by.x="sample",by.y="X")
anno2 <- anno[,c("sample","hidden_state","gender.demographic","race.demographic","prior_malignancy.diagnoses","prior_treatment.diagnoses")]
rownames(gsva_es) <- gsub("^HALLMARK_", "", rownames(gsva_es))

anno2 <- column_to_rownames(anno2,"sample")

sort_order <- order(anno2$hidden_state)
anno2_sorted <- anno2[sort_order, ]
gsva_es_sorted <- gsva_es[, rownames(anno2_sorted)]

hm_class <- c(
  ADIPOGENESIS = "Development",
  ALLOGRAFT_REJECTION = "Immune",
  ANDROGEN_RESPONSE = "Signaling",
  ANGIOGENESIS = "Development",
  APICAL_JUNCTION = "Development",
  APICAL_SURFACE = "Development",
  APOPTOSIS = "Stress",
  BILE_ACID_METABOLISM = "Metabolism",
  CHOLESTEROL_HOMEOSTASIS = "Metabolism",
  COAGULATION = "Development",
  COMPLEMENT = "Immune",
  DNA_REPAIR = "DNA_damage",
  E2F_TARGETS = "Proliferation",
  EPITHELIAL_MESENCHYMAL_TRANSITION = "Development",
  ESTROGEN_RESPONSE_EARLY = "Signaling",
  ESTROGEN_RESPONSE_LATE = "Signaling",
  FATTY_ACID_METABOLISM = "Metabolism",
  G2M_CHECKPOINT = "Proliferation",
  GLYCOLYSIS = "Metabolism",
  HEDGEHOG_SIGNALING = "Signaling",
  HEME_METABOLISM = "Metabolism",
  HYPOXIA = "Stress",
  IL2_STAT5_SIGNALING = "Immune",
  IL6_JAK_STAT3_SIGNALING = "Immune",
  INFLAMMATORY_RESPONSE = "Immune",
  INTERFERON_ALPHA_RESPONSE = "Immune",
  INTERFERON_GAMMA_RESPONSE = "Immune",
  KRAS_SIGNALING_DN = "Signaling",
  KRAS_SIGNALING_UP = "Signaling",
  MITOTIC_SPINDLE = "Proliferation",
  MTORC1_SIGNALING = "Signaling",
  MYC_TARGETS_V1 = "Proliferation",
  MYC_TARGETS_V2 = "Proliferation",
  MYOGENESIS = "Development",
  NOTCH_SIGNALING = "Signaling",
  OXIDATIVE_PHOSPHORYLATION = "Metabolism",
  P53_PATHWAY = "DNA_damage",
  PANCREAS_BETA_CELLS = "Development",
  PEROXISOME = "Metabolism",
  PI3K_AKT_MTOR_SIGNALING = "Signaling",
  PROTEIN_SECRETION = "Stress",
  REACTIVE_OXYGEN_SPECIES_PATHWAY = "Stress",
  SPERMATOGENESIS = "Development",
  TGF_BETA_SIGNALING = "Signaling",
  TNFA_SIGNALING_VIA_NFKB = "Immune",
  UNFOLDED_PROTEIN_RESPONSE = "Stress",
  UV_RESPONSE_DN = "Stress",
  UV_RESPONSE_UP = "Stress",
  WNT_BETA_CATENIN_SIGNALING = "Signaling",
  XENOBIOTIC_METABOLISM = "Metabolism"
)
hm_class_sorted <- sort(hm_class)
hm_class_sorted <- as.data.frame(hm_class_sorted)
gsva_es_sorted <- gsva_es_sorted[rownames(hm_class_sorted), ]
colnames(hm_class_sorted) <- c("Biological function")

anno2_sorted$hidden_state <- as.factor(anno2_sorted$hidden_state)
colnames(anno2_sorted) <- c("HS","Gender","Race","Malignancy","Treatment")
annotation_colors <- list(
  HS = c(
    "0" = "#2E9FDF",
    "1" = "#004529",
    "2" = "#7A142C",
    "3"="#6A51A3"
  ),Gender=c("female"="#B9B458","male"="#FDA542"),
  Race=c("asian"="#FE8214","black or african american"="#FB9684","not reported"="black","white"="#EC4D4E"),
  Malignancy=c("no"="#DA4C59","yes"="#C3AAD2"),
  Treatment=c("No"="#DA4C59","Yes"="#C3AAD2"),
  `Biological function` = c(
    Proliferation = "#DA4C59",
    Metabolism = "#74cb71",
    Immune = "#7D54A5",
    Signaling = "#B9A499",
    Stress = "#EAD27A",
    Development = "#FFFF33",
    DNA_damage = "#A65628"
  )
)

pheatmap(gsva_es_sorted,cluster_rows = F,cluster_cols = F, column_split = anno2_sorted$HS,row_split = hm_class_sorted$`Biological function`, 
          color= colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#DA4C59")), 
          show_colnames = F,scale = "row",show_rownames =T,#border_color = "black",
          annotation_col = anno2_sorted,annotation_row = hm_class_sorted,annotation_colors  = annotation_colors)
drawExport("heatmap.pdf",height = 10,width = 10)
###################################

########################LASSO
library(glmnet)
data <- merge(HMM_states,t(gsva_es_sorted),by.x="X",by.y="row.names")
data$hidden_state <- as.factor(data$hidden_state)
x_matrix <- as.matrix(data %>% select(-X, -hidden_state))
y_factor <- data$hidden_state


alphas <- seq(0.1, 1, by = 0.1) 
# alphas <- c(0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1)

cv_errors <- data.frame(alpha = alphas, cv_error = NA)

set.seed(123)
for (i in 1:length(alphas)) {
  # 对每个 alpha 进行交叉验证
  cv_fit <- cv.glmnet(
    x = x_matrix, 
    y = y_factor, 
    family = "multinomial",
    type.multinomial = "grouped",
    alpha = alphas[i],
    nfolds = 10
  )
  cv_errors$cv_error[i] <- min(cv_fit$cvm)
  cat(sprintf("Alpha = %.2f, CV Error = %.4f\n", 
              alphas[i], cv_errors$cv_error[i]))
}

# 选择最优 alpha（误差最小的）
best_alpha <- cv_errors$alpha[which.min(cv_errors$cv_error)]

final_cv_fit <- cv.glmnet(
  x = x_matrix, 
  y = y_factor, 
  family = "multinomial",
  type.multinomial = "grouped",
  alpha = best_alpha,
  nfolds = 10
)

coef_list <- coef(final_cv_fit, s = "lambda.min")
selected_features <- c()
for (i in 1:length(coef_list)) {
  non_zero_feats <- rownames(coef_list[[i]])[coef_list[[i]][, 1] != 0]
  non_zero_feats <- non_zero_feats[non_zero_feats != "(Intercept)"]
  selected_features <- union(selected_features, non_zero_feats)
}
print(selected_features)
##################################
library(tidyverse)
library(ggrepel)
library(FactoMineR)
library(magrittr)
library(factoextra)
library(RColorBrewer)
library(ggExtra)
library(patchwork)

pca <- t(gsva_es_sorted[selected_features,]) %>% 
  prcomp(., scale. = TRUE)  # 进行标准化 PCA 计算
var_explained <- pca$sdev^2 / sum(pca$sdev^2)

fviz_pca_ind(pca, axes = c(1,2),
             geom.ind = "point",
             pointshape = 20,
             pointsize = 4,
             addEllipses = TRUE,
             ellipse.level = 0.95,
             col.ind = anno2_sorted$HS,
             repel = FALSE) +  # 样本点是否需要避免重叠
  scale_color_manual(values = c("#2E9FDF","#004529","#7A142C","#6A51A3")) +
  labs(x = paste0("(PC1: ", round(var_explained[1] * 100, 2), "%)"),
       y = paste0("(PC2: ", round(var_explained[2] * 100, 2), "%)")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = "black", size = 12, margin = margin(t = 12)),
        axis.title.y = element_text(colour = "black", size = 12, margin = margin(r = 12)),
        axis.text = element_text(color = "black"),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 9),
        legend.position = c(1, 0),
        legend.justification = c(1, 0))


library(Rtsne)
iris_unique <- t(gsva_es_sorted[selected_features,])
tsne <- Rtsne(iris_unique, 
              initial_dims = 30,
              dims = 3,
              perplexity = 30)
tsne_data <- data.frame(tsne$Y,anno2_sorted$HS)
colnames(tsne_data)<-c("tSNE_1","tSNE_2",'tSNE_3',"HS")
iris_tsne<-Rtsne(iris_unique,pca=F)# Run TSNE
tsne<-iris_tsne[["Y"]]%>%as.data.frame()
tsne$HS<-anno2_sorted$HS

ggplot(tsne,aes(x=V1,y=V2,color=HS))+
    geom_point()+
    theme_classic()+
    labs(x='tSNE_1',y='tSNE_2')
    #stat_ellipse(aes(color = HS),level = 0.95, linetype = 2, show.legend = F)

library(Rtsne)
library(ggplot2)
library(ggsci)
library(readxl)
library(plotly)

p <- plotly::plot_ly(data  = tsne_data,
                x = ~tSNE_1,
                y = ~tSNE_2, 
                z = ~tSNE_3,
                marker = list(size = 4),
                color = ~HS,
                colors = c("#2E9FDF","#004529","#7A142C","#6A51A3")) 
plotly::orca(p, file = "tsne.pdf", width = 800, height = 600)
save_image(p,file ="tsne.pdf" )


library(umap)
library(ggplot2)
library(factoextra)
library(plotly)
set.seed(123)

# 自定义 UMAP 参数
custom_config <- umap.defaults
custom_config$n_neighbors <- 15      # 邻居数，通常选 5-50 之间
custom_config$min_dist <- 0.1        # 控制簇的紧凑程度，越小簇越集中
custom_config$metric <- "euclidean"  # 距离度量

# 执行降维
umap_result <- umap(iris_unique, config = custom_config)
umap_coords <- as.data.frame(umap_result$layout)

colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
umap_coords$HS <- anno2_sorted$HS
ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = factor(HS))) +
  geom_point(size = 1.5, alpha = 0.8) +
  #stat_ellipse(geom = "polygon",alpha = 0.15,aes(fill = factor(HS)),type = "norm",level = 0.95) +
  scale_color_manual(values = c("#2E9FDF","#004529","#7A142C","#6A51A3")) +
  scale_fill_manual(values = c("#2E9FDF","#004529","#7A142C","#6A51A3")) +
  labs(title = "2D UMAP Projection of Pathway Activities",
       x = "UMAP_1", y = "UMAP_2", color = "HS", fill = "HS") +
  theme_minimal()
drawExport("UMAP.pdf",height = 5,width = 5)

##############################################################
# 将GSVA结果(gsva_es)导出为Excel
library(openxlsx)
gsva_es_out <- data.frame(Pathway = rownames(gsva_es), gsva_es, check.names = FALSE)
write.xlsx(gsva_es_out, file = "gsva_es.xlsx", rowNames = FALSE)
