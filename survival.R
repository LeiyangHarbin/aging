library(survival)
library(survminer)
setwd("D:/aging/12.1GETX")
HMM_res <- read.csv("data_with_states.csv")

HMM_states <- HMM_res[,c("X","hidden_state")]
sur <- read_tsv("TCGA-LAML.survival.tsv")
data <- merge(HMM_states,sur,by.x="X",by.y = "sample")
data$OS.time <- data$OS.time/30
surv_obj <- Surv(time=data$OS.time, event=data$OS)
fit <- survfit(surv_obj ~ hidden_state, data=data)
MaleMod <- coxph(surv_obj ~ hidden_state, data = data)
print(MaleMod)
ggsurvplot(fit,
           data = data,
           conf.int= TRUE,pval = TRUE,
           surv.median.line ="hv", # 添加中位生存时间线
           risk.table = TRUE,
           risk.table.col ="strata",# 根据分层更改风险表颜色
           xlab ="Follow up time(month)",# 指定x轴标签
           legend = c(0.8,0.8),# 指定图例位置
           legend.title ="Kaplan-Meier",
           break.x.by=10,
           break.y.by=0.1,
           palette = c("#2E9FDF","#004529","#7A142C","#6A51A3"),
           ggtheme = theme_bw())
drawExport("KM.pdf",height = 10,width = 10)
