library(dplyr)
library(factoextra)
library(FNN)
library(ggplot2) 
library(ggsci)
library(glmnet)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
####################
setwd("D:/aging/11.merge")
norm_expr_data<-function(FM,pseudo_expr){
  FM2 <- FM/monocle:::estimateSizeFactorsForMatrix(FM)
  if (is.null(pseudo_expr)) 
    pseudo_expr <- 1
  FM2 <- FM2 + pseudo_expr
  FM2 <- log2(FM2)
  return(FM2)
}
genes<-cds@featureData@data[,c(1,2)]
expr <- Biobase::exprs(cds)
expr2<-data.frame(genes=genes,as.matrix(expr))
pData <- pData(cds)
cell <- rownames(pData)
pData$cell <-cell 
pData2 <- pData[,c("cell","State","Pseudotime")]
norm <- norm_expr_data(expr2[,-c(1:2)],1)
norm2 <- data.frame(expr2[,c(1:2)],norm)
###########################################
exp <- read.csv("D:/aging/data/bulk/Intergrated Multi-omics Sequencing Reveals Metabolic Reprograms in the Progression of ccRCC/bulkRNA_matrix_TPM.csv")
x <- exp$X
Tumor <- exp[,grep("T$", colnames(exp))]
Tumor <- cbind(x,Tumor)
Normal <- exp[,grep("N$", colnames(exp))]
Normal <- cbind(x,Normal)
rownames(Tumor) <- Tumor$x
rownames(Normal) <- Normal$x
####################################单细胞和bulk merge 
exp <- Tumor
input_mock <- merge(exp,norm2[,-2],by.x="x",by.y="genes.gene_short_name")
gene <- input_mock$x
input <- input_mock[,-1]
rownames(input) <- gene
#input <- input_mock[, -which(names(input_mock) == "genes.use_for_ordering")]
batch<-factor(c(rep(1,ncol(exp[,-1])),rep(2,ncol(norm2[,-c(1:2)]))))
####################################combat去批次
combat = ComBat(dat=as.matrix(input), batch=batch)
combat2<-data.frame(genes=gene,combat)
combat2[,1]<-as.character(combat2[,1])
combat2[combat2<0]<-0
vardata<-apply(combat2[,-1],1,var)
combat2$vardata<-vardata
combat2<-combat2[order(combat2$vardata,decreasing=T),]
combat_filter_T<- combat2[1:round((nrow(combat2)*30)/100),-ncol(combat2)]
save(combat_filter_T,file = "combat_filter_T.rda")
#####################################################
knnProject<-function(FM,knn_k = 10, num_dim = 25,...){
  message("Performing dimensionality reduction with PCA")
  irlba_res <- prcomp(t(combat_filter_T[,-1]))
  irlba_pca_res <- irlba_res$x[,c(1:num_dim)]
  rownames(irlba_pca_res)<-colnames(FM[,-1])
  scrna_pca<-irlba_pca_res[grep("^.+_.+-[0-9]+$", rownames(irlba_pca_res)),]
  tcga_pca <- irlba_pca_res[grep("^.+_.+-[0-9]+$", rownames(irlba_pca_res), invert=T), ]
  knn.list <- list()
  for(patient in 1:length(rownames(tcga_pca))){
    print(patient)
    patient2<-rownames(tcga_pca)[patient]
    tcga.mtx<-tcga_pca[row.names(tcga_pca)==patient2,]
    res<-get.knnx(scrna_pca, t(data.frame(tcga.mtx)), k=knn_k, algorithm="kd_tree")
    knn.list[[patient2]]<-res
    row.names(knn.list[[patient2]]$nn.index)<-patient2
    colnames(knn.list[[patient2]]$nn.index)<-rownames(scrna_pca)[res$nn.index]
    row.names(knn.list[[patient2]]$nn.dist) <- patient2
    leninner<-length(grep(rownames(scrna_pca)[res$nn.index],pattern='inner'))
    lenouter<-length(grep(rownames(scrna_pca)[res$nn.index],pattern='outer'))
    
  }
  return(knn.list)
}
combat_knn_T<-knnProject(combat_filter_T,num_dim=25)
save(combat_knn_T,file = "combat_knn_T.rda")

knn_mock_pseudospace <- NULL
ratio_mock_outer_inner<- NULL 
knn_mock_corr<-NULL
for(patient in names(combat_knn)){
  scRNAseq<-gsub(colnames(combat_knn[[patient]]$nn.index),pattern='\\.',replacement='-')
  knn_pseudospace<-mean(pData2[pData2[,1]%in%scRNAseq,'Pseudotime'])
  knn_mock_pseudospace<-c(knn_mock_pseudospace,knn_pseudospace)
  tab<-as.data.frame(table(pData2[pData2[,1]%in%scRNAseq,'State']))
  a <- as.character(tab[which.max(tab$Freq),1])
  ratio_mock_outer_inner <- c(ratio_mock_outer_inner,a)
  idx_tcga<-grep(colnames(combat_filter),pattern=patient)
  exp_profile_tcga<-combat_filter[,idx_tcga]
  scRNAseq<-colnames(combat_knn[[patient]]$nn.index)
  idx_sc<-grep(colnames(combat_filter),pattern=paste(scRNAseq,collapse='|'))
  exp_profile_sc<-combat_filter[,idx_sc]
  knn_mock_corr<-c(knn_mock_corr,mean(cor(exp_profile_tcga,exp_profile_sc,method='spearman')))
}
#tumors<-sapply(strsplit(names(combat_knn),split='\\.'),'[[',1)
knn_df<-data.frame(patients=names(combat_knn),pseudospace=knn_mock_pseudospace,ratio=ratio_mock_outer_inner)
knn_df$patients <- gsub(pattern = "\\.", replacement = "-", x = knn_df$patients)
save(knn_df,file = "knn_df_T.rda")
###########################################################################
pseudospace_input_sort <- knn_df[order(knn_df$pseudospace,decreasing=F),]
pseudospace_input_sort$patients <- gsub(pattern = "\\-", replacement = ".", x = pseudospace_input_sort$patients)
pseudospace_input_sort<-pseudospace_input_sort[which(pseudospace_input_sort$patients%in%colnames(exp)),]
rownames(pseudospace_input_sort) <- pseudospace_input_sort$patients
exp3 <- as.data.frame(t(Tumor[,-1]))
exp3$id <- rownames(exp3)
EXP_ALL_sort <- merge(pseudospace_input_sort,exp3,by.x="patients",by.y="id")
rownames(EXP_ALL_sort) <- EXP_ALL_sort$patients
EXP_ALL_sort <- EXP_ALL_sort[,-1]
EXP_ALL_sort <- EXP_ALL_sort[, colSums(is.na(EXP_ALL_sort)) == 0]
sample <- rownames(EXP_ALL_sort)
EXP_ALL_sort <- apply(EXP_ALL_sort, 2, as.numeric)
rownames(EXP_ALL_sort) <- sample
save(EXP_ALL_sort,file = "EXP_ALL_sort_T.rda")
###################################lasso
library(glmnet)
set.seed(123)
#选择基因表达0的比例40%之下的基因
EXP_ALL_qc <- EXP_ALL_sort[, colMeans(EXP_ALL_sort == 0, na.rm = TRUE) <= 0.4]
cvfit <-cv.glmnet(x =EXP_ALL_qc[,-c(1:2)], y = EXP_ALL_qc[,1])
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
drawExport("lassoT.pdf")
dfcoef<-as.matrix(coef(cvfit,s= "lambda.min"))
genes_to_use<- names(dfcoef[dfcoef[,1]!=0,])[-1]
#genes_to_use2<- names(dfcoef[dfcoef[,1]<=-0.001|dfcoef[,1]>=0.001,])[-1]
#删除基因中带有-和：的基因
filtered_strings <- genes_to_use[!grepl("-", genes_to_use)]
filtered_strings <- filtered_strings[!grepl(":", filtered_strings)]
EXP_ALL_sort2 <- as.data.frame(EXP_ALL_qc)
HMM_data <- EXP_ALL_sort2[order(EXP_ALL_sort2$pseudospace,decreasing=F),]
data <- HMM_data[,filtered_strings]
#data <- data[, colMeans(data == 0, na.rm = TRUE) <= 0.4]
write.csv(data,file = "HMM_data_T.csv")
#########################################转python




