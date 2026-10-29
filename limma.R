library(limma)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggVolcano)
library(ggplot2)
library(ggsankey)
library(dplyr)
library(patchwork)
setwd("D:/aging/13.bio/1.gsva")
HMM_res <- read.csv("D:/aging/12.1getx/data_with_states.csv")
#exp <- read_tsv("D:/aging/12.BLCA/TCGA-BLCA.star_fpkm.tsv")
HMM_states <- HMM_res[,c("X","hidden_state")]

HMM_states[which(HMM_states$hidden_state==0),2] <- "State0"
HMM_states[which(HMM_states$hidden_state==1),2] <- "State1"
HMM_states[which(HMM_states$hidden_state==2),2] <- "State2"
HMM_states[which(HMM_states$hidden_state==3),2] <- "State3"
HMM_states$hidden_state <- factor(HMM_states$hidden_state, levels = c("State0", "State1", "State2", "State3"))

groups <- HMM_states$hidden_state
names(groups) <- HMM_states$X

# 确保GSVA得分矩阵列名与分组匹配
exp <- exp[, HMM_states$X]

design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)
rownames(design) <- HMM_states$X

contrasts <- makeContrasts(
  State0_vs_Others = State0 - (State1 + State2 + State3)/3,
  State1_vs_Others = State1 - (State0 + State2 + State3)/3,
  State2_vs_Others = State2 - (State0 + State1 + State3)/3,
  State3_vs_Others = State3 - (State0 + State1 + State2)/3,
  levels = design
)
fit <- lmFit(exp, design)
fit_contrasts <- contrasts.fit(fit, contrasts)
fit_ebayes <- eBayes(fit_contrasts)

results <- list()
for (i in 1:ncol(contrasts)) {
  contrast_name <- colnames(contrasts)[i]
  top_table <- topTable(fit_ebayes, coef = i, number = Inf, adjust.method = "BH")
  # 添加通路名作为列
  top_table <- top_table %>%
    rownames_to_column("Pathway") %>%
    mutate(Contrast = contrast_name)
  results[[contrast_name]] <- top_table
}

all_results <- bind_rows(results)

sig_results <- all_results %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 0.6)  
table(sig_results$Contrast)
write.csv(sig_results, "significant_genes.csv", row.names = FALSE)

df <- sig_results %>% filter(Contrast == "State0_vs_Others")

gradual_volcano(df,
                x="logFC",
                y="adj.P.Val",
                label="Pathway",
                label_number=20,
                output= FALSE)
drawExport("S0.pdf",height = 10,width = 10)
df <- all_results %>% filter(Contrast == "State1_vs_Others")
gradual_volcano(df,
                x="logFC",
                y="adj.P.Val",
                label="Pathway",
                label_number=10,
                output= FALSE)
drawExport("S1.pdf",height = 10,width = 10)

df <- all_results %>% filter(Contrast == "State2_vs_Others")
gradual_volcano(df,
                x="logFC",
                y="adj.P.Val",
                label="Pathway",
                label_number=10,
                output= FALSE)
drawExport("S2.pdf",height = 10,width = 10)

df <- all_results %>% filter(Contrast == "State3_vs_Others")
gradual_volcano(df,
                x="logFC",
                y="adj.P.Val",
                label="Pathway",
                label_number=10,
                output= FALSE)
drawExport("S3.pdf",height = 10,width = 10)
#############################################
cluster <- readxl::read_xlsx("NIHMS911030-supplement-8.xlsx",sheet = 2)
HMM_states$clean_x <- sub("-.*$", "", HMM_states$X)
HMM_states <- HMM_states %>%
  mutate(clean_x = str_sub(X, end = -5))
plotdat <- merge(HMM_states,cluster[,c("Case ID","mRNA cluster")],by.x="clean_x",by.y = "Case ID")
df <- plotdat %>%
  dplyr::select(hidden_state, `mRNA cluster`) %>%
  make_long(hidden_state, `mRNA cluster`) %>%
  mutate(stage = ifelse(x == "hidden_state", "HS", "mRNA cluster"))
dagg <- df %>%
  group_by(stage, node) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(stage) %>%
  mutate(pct = round(n / sum(n) * 100, 2))

df2 <- df %>%
  left_join(dagg, by = c("stage", "node"))

ggplot(
  df2, 
  aes(
    x = x, 
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = node,
    label = paste0(node, "\nn=", n, " (", pct, "%)")
  )) +
  geom_sankey(
    flow.alpha = 0.6,
    node.color = "black", 
    node.width = 15,
    width = 0.17,
    smooth = 8,
    show.legend = FALSE
  ) +
  geom_sankey_text(
    size = 3,
    color = "black",
    hjust = 0.45,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      '0' = "#FFB6C1",
      '1' = "#AFEEEE",
      '2' = "#FFDAB9",
      '3' = "#D8BFD8",
      "Neuronal"="#EEE0CC",
      "Basal_squamous"="#dd8c1b",
      "Luminal"="#4b036d",
      "Luminal_infiltrated"="#aa98e9",
      "Luminal_papillary"="#94aa96"
    )) +
  scale_x_discrete(
    labels = c("HS", "mRNA cluster"),
    expand = expansion(add = 0.3)
  ) +
  theme_bw() +
  theme(panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold")
  )
drawExport("sangji.pdf",height = 6,width = 12)
