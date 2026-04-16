#!/usr/bin/env Rscript
# =============================================================================
# 📚 迁移体相关胃癌生信文献复刻 | Day3 (风险模型构建与验证)
# 基于 LASSO 选出的 8 个预后基因，构建风险评分模型
# 分组方式：最佳截止值 (optimal cut-off by surv_cutpoint)
# 输出:
#   - 风险评分分布图 + 生存状态散点图
#   - KM 生存曲线 (高低风险组)
#   - 时间依赖 ROC 曲线 (1, 3, 5 年) 及 AUC
#   - 风险评分系数表及样本分组表
# 方法: glmnet, survival, timeROC, ggplot2
# 数据: TCGA-STAD 肿瘤样本 (n=375)
# 作者: 郭俊
# 日期: 2026-04-16
# =============================================================================

rm(list = ls())
PROJECT_DIR <- "/data/nas1/guojun_OD/project/05_paper_reproduction/Migrasome_STAD_Day3_RiskModel"
if (!dir.exists(PROJECT_DIR)) dir.create(PROJECT_DIR, recursive = TRUE)
setwd(PROJECT_DIR)
for (d in c("data", "results", "scripts")) {
  if (!dir.exists(d)) dir.create(d)
}

# -------------------- 1. 加载包 --------------------
required_packages <- c(
  "tidyverse", "survival", "survminer", "glmnet", "timeROC", 
  "ggplot2", "patchwork", "gridExtra"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("survminer", "timeROC")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# -------------------- 2. 读取 Day2 输出 --------------------
# 读取预后基因列表
prognostic_genes <- readLines("../Migrasome_STAD_Day2_Prognosis/results/prognostic_genes.txt")
cat("预后基因 (", length(prognostic_genes), "个): ", paste(prognostic_genes, collapse = ", "), "\n")

# 读取 LASSO 系数（文章 lambda 对应的系数）
lasso_coef <- read.csv("../Migrasome_STAD_Day2_Prognosis/results/LASSO_selected_genes.csv", stringsAsFactors = FALSE)
# 确保系数顺序与基因列表一致
lasso_coef <- lasso_coef[match(prognostic_genes, lasso_coef$Gene), ]
cat("LASSO 系数:\n")
print(lasso_coef)

# 读取表达矩阵和生存数据（需要从 Day2 保存的 RDS 加载）
# 如果 Day2 没有保存，可以重新构建，但为了效率，我们加载 Day2 的工作空间
# 这里假设 Day2 保存了 expr_tumor 和 surv_tumor 对象
if (file.exists("../Migrasome_STAD_Day2_Prognosis/results/Day2_workspace.RData")) {
  load("../Migrasome_STAD_Day2_Prognosis/results/Day2_workspace.RData")
  # 确保 expr_tumor 和 surv_tumor 存在
} else {
  # 备用：重新从原始数据构建（简化版，需要确保路径正确）
  cat("未找到 Day2 workspace，尝试重新构建表达矩阵和生存数据...\n")
  # 读取 Day1 的表达矩阵
  old_data <- readRDS("../Migrasome_STAD_Day1_AllSamples/data/STAD_all_samples.rds")
  counts <- old_data$counts
  sample_type <- old_data$sample_type
  # 重新运行 DESeq2 获取 VST
  library(DESeq2)
  colData <- data.frame(condition = factor(sample_type, levels = c("Normal", "Tumor")))
  rownames(colData) <- colnames(counts)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = FALSE)
  expr_vst <- assay(vsd)
  # 提取肿瘤样本
  tumor_samples <- colnames(counts)[sample_type == "Tumor"]
  expr_tumor <- expr_vst[, intersect(tumor_samples, colnames(expr_vst))]
  # 生存数据
  old_raw <- readRDS("../Migrasome_STAD_Day1/data/TCGA_STAD_raw.rds")
  sample_info <- old_raw$sample_info
  surv_df <- data.frame(barcode = colnames(counts), sample_type = sample_type)
  surv_df$OS.time <- ifelse(!is.na(sample_info$days_to_death), 
                            sample_info$days_to_death,
                            sample_info$days_to_last_follow_up)
  surv_df$OS.event <- ifelse(sample_info$vital_status == "Dead", 1, 0)
  surv_tumor <- surv_df[surv_df$sample_type == "Tumor", ]
  rownames(surv_tumor) <- surv_tumor$barcode
  # 清理生存数据
  surv_tumor <- surv_tumor[!is.na(surv_tumor$OS.time) & surv_tumor$OS.time > 0, ]
  surv_tumor <- surv_tumor[!is.na(surv_tumor$OS.event), ]
  surv_tumor$OS.event <- as.numeric(surv_tumor$OS.event)
  # 对齐样本
  common_samples <- intersect(colnames(expr_tumor), rownames(surv_tumor))
  expr_tumor <- expr_tumor[, common_samples]
  surv_tumor <- surv_tumor[common_samples, ]
  cat("样本数用于风险模型:", ncol(expr_tumor), "\n")
}

# 提取 8 个预后基因的表达矩阵（样本 x 基因）
expr_prog <- t(expr_tumor[prognostic_genes, ])  # 样本 x 基因
# 确保表达矩阵与生存数据样本顺序一致
expr_prog <- expr_prog[rownames(surv_tumor), ]

# -------------------- 3. 计算风险评分 --------------------
# 风险评分公式: RiskScore = sum(coef_i * expr_i)
coef_vector <- lasso_coef$Coefficient
names(coef_vector) <- lasso_coef$Gene
risk_score <- as.numeric(expr_prog %*% coef_vector)
surv_tumor$RiskScore <- risk_score

# -------------------- 4. 计算最佳截止值 (optimal cut-off) --------------------
cat("\n🔍 计算最佳截止值 (基于 log-rank 最大化)...\n")
library(survminer)

cut_data <- data.frame(
  risk_score = risk_score,
  time = surv_tumor$OS.time,
  event = surv_tumor$OS.event
)

# 寻找最佳切点（每组至少 30% 样本）
cutpoint_obj <- surv_cutpoint(
  data = cut_data,
  time = "time",
  event = "event",
  variables = "risk_score",
  minprop = 0.3,
  progressbar = TRUE
)

optimal_cut <- cutpoint_obj$cutpoint[1, "cutpoint"]
cat("最佳截止值 (optimal cut-off):", optimal_cut, "\n")

# 根据最佳截止值分组
surv_tumor$RiskGroup <- ifelse(risk_score > optimal_cut, "High", "Low")
surv_tumor$RiskGroup <- factor(surv_tumor$RiskGroup, levels = c("Low", "High"))
cat("高风险组样本数:", sum(surv_tumor$RiskGroup == "High"), "\n")
cat("低风险组样本数:", sum(surv_tumor$RiskGroup == "Low"), "\n")

# 保存最佳截止值和分组信息
write.csv(data.frame(cutoff = optimal_cut), "results/optimal_cutoff_value.csv", row.names = FALSE)
write.csv(surv_tumor[, c("RiskScore", "RiskGroup")], "results/risk_scores_optimal.csv", row.names = TRUE)

# (可选) 原中位数分组代码，供对比，注释掉
# median_risk <- median(risk_score)
# surv_tumor$RiskGroup_median <- ifelse(risk_score > median_risk, "High", "Low")
# surv_tumor$RiskGroup_median <- factor(surv_tumor$RiskGroup_median, levels = c("Low", "High"))
# cat("中位数截止值:", median_risk, "\n")
# cat("中位数分组 - 高风险组:", sum(surv_tumor$RiskGroup_median == "High"), "\n")
# cat("中位数分组 - 低风险组:", sum(surv_tumor$RiskGroup_median == "Low"), "\n")

# -------------------- 5. 风险评分分布图和生存状态散点图 --------------------
# 按风险评分排序
surv_tumor_sorted <- surv_tumor[order(risk_score), ]
surv_tumor_sorted$Index <- 1:nrow(surv_tumor_sorted)

# 风险评分分布图
p1 <- ggplot(surv_tumor_sorted, aes(x = Index, y = RiskScore, fill = RiskGroup)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("Low" = "#1F77B4", "High" = "#D62728")) +
  labs(title = "Risk Score Distribution (optimal cut-off)", x = "Patients (increasing risk score)", y = "Risk Score") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

# 生存状态散点图（绿色=存活，红色=死亡）
p2 <- ggplot(surv_tumor_sorted, aes(x = Index, y = OS.time, color = as.factor(OS.event))) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c("0" = "#2CA02C", "1" = "#D62728"), labels = c("Alive", "Dead")) +
  labs(title = "Survival Status", x = "Patients (increasing risk score)", y = "Overall survival (days)", color = "Status") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

# 拼合两个图
combined_plot <- p1 / p2 + plot_annotation(
  title = "Risk Model Performance in TCGA-STAD (Optimal Cut-off)",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)
ggsave("results/RiskScore_distribution_survival_optimal.pdf", combined_plot, width = 10, height = 8)
ggsave("results/RiskScore_distribution_survival_optimal.png", combined_plot, width = 10, height = 8, dpi = 300)

# -------------------- 6. KM 生存曲线 --------------------
fit_km <- survfit(Surv(OS.time, OS.event) ~ RiskGroup, data = surv_tumor)
# 计算 log-rank p 值
logrank_p <- surv_pvalue(fit_km)$pval
logrank_p_text <- ifelse(logrank_p < 0.001, "p < 0.001", paste0("p = ", round(logrank_p, 4)))

km_plot <- ggsurvplot(
  fit_km, data = surv_tumor, pval = logrank_p_text, conf.int = TRUE,
  risk.table = TRUE, risk.table.col = "strata", palette = c("#1F77B4", "#D62728"),
  xlab = "Time (days)", ylab = "Overall Survival Probability",
  title = paste0("Kaplan-Meier Curves (optimal cut-off = ", round(optimal_cut, 4), ")"),
  legend.title = "Risk Group", legend.labs = c("Low", "High"),
  ggtheme = theme_minimal()
)

# 保存 KM 曲线（PDF + PNG）
pdf("results/KM_curve_optimal.pdf", width = 7, height = 7)
print(km_plot)
dev.off()
png("results/KM_curve_optimal.png", width = 7, height = 7, units = "in", res = 300)
print(km_plot)
dev.off()

# 额外保存风险表格单独文件（可选）
risk_table <- km_plot$table
ggsave("results/KM_risk_table_optimal.pdf", risk_table, width = 7, height = 3)
ggsave("results/KM_risk_table_optimal.png", risk_table, width = 7, height = 3, dpi = 300)

# -------------------- 7. 时间依赖 ROC 曲线 (1, 3, 5 年) --------------------
# 使用 timeROC 包
roc_obj <- timeROC(
  T = surv_tumor$OS.time,
  delta = surv_tumor$OS.event,
  marker = surv_tumor$RiskScore,
  cause = 1,
  times = c(365, 1095, 1825),  # 1, 3, 5 年 (天数)
  iid = TRUE
)

# 提取 AUC 值
auc_1 <- round(roc_obj$AUC[1], 3)
auc_3 <- round(roc_obj$AUC[2], 3)
auc_5 <- round(roc_obj$AUC[3], 3)
cat("AUC at 1 year:", auc_1, "\n")
cat("AUC at 3 years:", auc_3, "\n")
cat("AUC at 5 years:", auc_5, "\n")

# 绘制 ROC 曲线
roc_data <- data.frame(
  tpr_1 = roc_obj$TP[, 1], fpr_1 = roc_obj$FP[, 1],
  tpr_3 = roc_obj$TP[, 2], fpr_3 = roc_obj$FP[, 2],
  tpr_5 = roc_obj$TP[, 3], fpr_5 = roc_obj$FP[, 3]
)

roc_plot <- ggplot() +
  geom_line(data = roc_data, aes(x = fpr_1, y = tpr_1, color = "1-year"), size = 1.2) +
  geom_line(data = roc_data, aes(x = fpr_3, y = tpr_3, color = "3-year"), size = 1.2) +
  geom_line(data = roc_data, aes(x = fpr_5, y = tpr_5, color = "5-year"), size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("1-year" = "#1F77B4", "3-year" = "#2CA02C", "5-year" = "#D62728"),
                     name = "Time") +
  labs(title = "Time-Dependent ROC Curves (Optimal Cut-off)",
       x = "1 - Specificity", y = "Sensitivity",
       subtitle = paste0("AUC: 1-year = ", auc_1, ", 3-year = ", auc_3, ", 5-year = ", auc_5)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave("results/TimeROC_curves_optimal.pdf", roc_plot, width = 7, height = 7)
ggsave("results/TimeROC_curves_optimal.png", roc_plot, width = 7, height = 7, dpi = 300)

# 保存 AUC 值
write.csv(data.frame(Time = c("1-year", "3-year", "5-year"), AUC = c(auc_1, auc_3, auc_5)),
          "results/ROC_AUC_values_optimal.csv", row.names = FALSE)

# -------------------- 8. 保存所有对象和工作空间 --------------------
save.image("results/Day3_workspace_optimal.RData")
cat("\n🎉 Day3 分析完成！所有结果保存在 results/ 目录。\n")
# ============================================================================
# 使用文章系数重新构建风险模型
# ============================================================================

# 文章中的系数（与基因顺序对应）
article_coef <- c(
  BMP1 = 0.089332494,
  CPQ = 0.0003596,
  PDGFD = 0.074722162,
  TSPAN5 = 0.026623357,
  TSPAN7 = 0.047167114,
  TGFB2 = 0.023629132,
  WNT11 = 0.010867593,
  LEFTY1 = -0.03409899
)

prognostic_genes_article <- names(article_coef)

# 提取表达矩阵（确保基因存在）
expr_prog_article <- t(expr_tumor[prognostic_genes_article, ])
expr_prog_article <- expr_prog_article[rownames(surv_tumor), ]  # 对齐样本

# 计算风险评分
risk_score_article <- as.numeric(expr_prog_article %*% article_coef)
surv_tumor$RiskScore_article <- risk_score_article

# 计算最佳截止值（使用 surv_cutpoint）
library(survminer)
cut_data_article <- data.frame(
  risk_score = risk_score_article,
  time = surv_tumor$OS.time,
  event = surv_tumor$OS.event
)
cutpoint_article <- surv_cutpoint(cut_data_article, time = "time", event = "event", 
                                  variables = "risk_score", minprop = 0.3)
optimal_cut_article <- cutpoint_article$cutpoint[1, "cutpoint"]
cat("文章系数模型的最佳截止值:", optimal_cut_article, "\n")

# 分组
surv_tumor$RiskGroup_article <- ifelse(risk_score_article > optimal_cut_article, "High", "Low")
surv_tumor$RiskGroup_article <- factor(surv_tumor$RiskGroup_article, levels = c("Low", "High"))
table(surv_tumor$RiskGroup_article)

# KM 曲线
fit_km_article <- survfit(Surv(OS.time, OS.event) ~ RiskGroup_article, data = surv_tumor)
logrank_p_article <- surv_pvalue(fit_km_article)$pval
cat("Log-rank p value:", logrank_p_article, "\n")

# 时间依赖 ROC (1,3,5年)
library(timeROC)
roc_article <- timeROC(T = surv_tumor$OS.time, delta = surv_tumor$OS.event,
                       marker = risk_score_article, cause = 1,
                       times = c(365, 1095, 1825), iid = TRUE)
auc_1_article <- round(roc_article$AUC[1], 3)
auc_3_article <- round(roc_article$AUC[2], 3)
auc_5_article <- round(roc_article$AUC[3], 3)
cat("AUC (1,3,5年):", auc_1_article, auc_3_article, auc_5_article, "\n")

# 可选：保存结果
write.csv(surv_tumor[, c("RiskScore_article", "RiskGroup_article")], 
          "results/risk_scores_article_coef.csv", row.names = TRUE)
