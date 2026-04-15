#!/usr/bin/env Rscript
# =============================================================================
# 📚 迁移体相关胃癌生信文献复刻 | Day2 (Figure 2A–F)
# 复刻内容:
#   2A: 单变量 Cox 回归森林图 (P < 0.2)
#   2B: PH 假设检验 (P > 0.05)
#   2C/D: LASSO 回归 (10折交叉验证) + 系数路径图
#   2E: 8个预后基因在肿瘤 vs 正常中的箱线图
#   2F: 8个预后基因在肿瘤样本中的相关性热图
# 方法: survival, glmnet, ggplot2, ComplexHeatmap
# 数据: TCGA-STAD (基于 Day1 的表达矩阵和生存数据)
# 作者: 郭俊
# 日期: 2026-04-13
# =============================================================================

rm(list = ls())
PROJECT_DIR <- "/data/nas1/guojun_OD/project/05_paper_reproduction/Migrasome_STAD_Day2_Prognosis"
if (!dir.exists(PROJECT_DIR)) dir.create(PROJECT_DIR, recursive = TRUE)
setwd(PROJECT_DIR)
for (d in c("data", "results", "scripts")) {
  if (!dir.exists(d)) dir.create(d)
}

# -------------------- 1. 加载包 --------------------
required_packages <- c(
  "DESeq2", "tidyverse", "ggplot2", "ggrepel", "survival", "survminer",
  "glmnet", "ComplexHeatmap", "circlize", "ggpubr", "patchwork"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("DESeq2", "ComplexHeatmap", "survminer")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# -------------------- 2. 读取 Day1 结果 --------------------
# 读取候选基因列表
candidate_genes <- readLines("../Migrasome_STAD_Day1_AllSamples/results/candidate_genes.txt")
cat("候选基因数量:", length(candidate_genes), "\n")
cat("候选基因:", paste(candidate_genes, collapse = ", "), "\n")

# 读取 Day1 保存的表达矩阵 (VST 转换后的数据，包含所有样本)
old_data <- readRDS("../Migrasome_STAD_Day1_AllSamples/data/STAD_all_samples.rds")
counts <- old_data$counts
sample_type <- old_data$sample_type

# 重新运行 DESeq2 获取 VST 表达矩阵 (或直接从 Day1 读取 vst 对象)
# 简便方法：重新构建 dds 并 vst
cat("\n🔧 准备表达数据...\n")
colData <- data.frame(condition = factor(sample_type, levels = c("Normal", "Tumor")))
rownames(colData) <- colnames(counts)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
expr_vst <- assay(vsd)  # 基因 x 样本

# 提取肿瘤样本的表达矩阵 (用于后续相关性和 LASSO)
tumor_samples <- colnames(counts)[sample_type == "Tumor"]
expr_tumor <- expr_vst[, intersect(tumor_samples, colnames(expr_vst))]

# -------------------- 3. 获取 TCGA 生存数据 --------------------
# 从 UCSC Xena 或之前保存的临床数据中读取
# 假设我们已有 survival 数据: 从 Day1 的原始数据中提取 sample_info
old_raw <- readRDS("../Migrasome_STAD_Day1/data/TCGA_STAD_raw.rds")
sample_info <- old_raw$sample_info
# 构建生存数据框
surv_df <- data.frame(
  barcode = colnames(counts),
  sample_type = sample_type,
  stringsAsFactors = FALSE
)
# 匹配临床信息 (需要从 sample_info 中提取 days_to_last_follow_up, days_to_death, vital_status)
# 注意: sample_info 是 colData 对象，需要正确提取
if (is(sample_info, "DataFrame")) {
  surv_df$OS.time <- ifelse(!is.na(sample_info$days_to_death), 
                            sample_info$days_to_death,
                            sample_info$days_to_last_follow_up)
  surv_df$OS.event <- ifelse(sample_info$vital_status == "Dead", 1, 0)
} else {
  # 备用: 直接使用之前保存的生存数据文件 (需手动准备)
  # 这里提供一个示例：从 TCGA 下载的临床文件中读取
  cat("⚠️ 未找到生存数据，请确保 sample_info 包含 days_to_death, days_to_last_follow_up, vital_status\n")
  # 临时: 生成模拟生存数据 (仅用于测试，实际使用时需替换)
  # surv_df$OS.time <- runif(nrow(surv_df), 100, 3000)
  # surv_df$OS.event <- rbinom(nrow(surv_df), 1, 0.4)
  stop("需要正确的生存数据。请检查 sample_info 结构。")
}

# 仅保留肿瘤样本的生存数据
surv_tumor <- surv_df[surv_df$sample_type == "Tumor", ]
rownames(surv_tumor) <- surv_tumor$barcode

# 清理生存数据：排除 OS.time 为 NA、0 或负数的样本
surv_tumor <- surv_tumor[!is.na(surv_tumor$OS.time) & surv_tumor$OS.time > 0, ]
# 同时排除 OS.event 为 NA 的样本
surv_tumor <- surv_tumor[!is.na(surv_tumor$OS.event), ]
# 可选：确保 OS.event 为 0/1（如果不是，需要转换）
surv_tumor$OS.event <- as.numeric(surv_tumor$OS.event)

# 确保表达矩阵和生存数据样本一致
common_samples <- intersect(colnames(expr_tumor), surv_tumor$barcode)
expr_tumor <- expr_tumor[, common_samples]
surv_tumor <- surv_tumor[common_samples, ]
cat("肿瘤样本数用于生存分析:", ncol(expr_tumor), "\n")

# -------------------- 4. 单变量 Cox 回归 (Figure 2A) --------------------
cat("\n📊 单变量 Cox 回归分析 (P < 0.2)...\n")
cox_results <- data.frame()
for (gene in candidate_genes) {
  if (!gene %in% rownames(expr_tumor)) next
  expr_gene <- expr_tumor[gene, ]
  cox_data <- data.frame(
    time = surv_tumor$OS.time,
    status = surv_tumor$OS.event,
    expression = as.numeric(expr_gene)
  )
  cox_model <- coxph(Surv(time, status) ~ expression, data = cox_data)
  cox_sum <- summary(cox_model)
  hr <- round(cox_sum$coefficients[1, "exp(coef)"], 3)
  hr_lower <- round(cox_sum$conf.int[1, "lower .95"], 3)
  hr_upper <- round(cox_sum$conf.int[1, "upper .95"], 3)
  pval <- cox_sum$coefficients[1, "Pr(>|z|)"]
  cox_results <- rbind(cox_results, data.frame(
    Gene = gene,
    HR = hr,
    HR_lower = hr_lower,
    HR_upper = hr_upper,
    P_value = pval,
    stringsAsFactors = FALSE
  ))
}
# 筛选 P < 0.2
cox_sig <- cox_results[cox_results$P_value < 0.2, ]
cox_sig <- cox_sig[order(cox_sig$P_value), ]
cat("单变量 Cox 显著基因 (P<0.2):", nrow(cox_sig), "\n")
write.csv(cox_sig, "results/univariate_Cox_P0.2.csv", row.names = FALSE)

# 森林图 (Figure 2A)
forest_data <- cox_sig
forest_data$Gene <- factor(forest_data$Gene, levels = rev(forest_data$Gene))

# 绘图
forest_plot <- ggplot(forest_data, aes(x = HR, y = Gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(xmin = HR_lower, xmax = HR_upper, y = Gene),
                width = 0.2, orientation = "y",   # 注意 orientation = "y"
                color = "#1F77B4", linewidth = 0.8)+# 注意 linewidth
  geom_point(aes(color = ifelse(HR > 1, "Risk", "Protective")), size = 3) +
  scale_color_manual(values = c("Risk" = "#D62728", "Protective" = "#2CA02C")) +
  labs(title = "Univariate Cox Regression (P < 0.2)",
       x = "Hazard Ratio (95% CI)", y = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/Figure2A_univariate_forest.pdf", forest_plot, width = 8, height = 6)
ggsave("results/Figure2A_univariate_forest.png", forest_plot, width = 8, height = 6, dpi = 300)
# -------------------- 5. PH 假设检验 (Figure 2B) --------------------
# 加载必要的包
library(tidyverse)
library(survival)
library(patchwork)

cat("\n🔍 进行 PH 假设检验并准备残差数据...\n")

# 存储结果
ph_results <- data.frame()
ph_plots <- list()

for (gene in cox_sig$Gene) {
  cat("  Processing", gene, "...\n")
  
  # 提取表达量（确保与生存数据样本匹配）
  expr_gene <- expr_tumor[gene, rownames(surv_tumor)]
  
  # 构建 Cox 模型（只拟合一次）
  cox_data <- data.frame(
    time = surv_tumor$OS.time,
    status = surv_tumor$OS.event,
    expression = as.numeric(expr_gene)
  )
  cox_model <- coxph(Surv(time, status) ~ expression, data = cox_data)
  
  # 进行 PH 检验
  ph_test <- cox.zph(cox_model)
  ph_p <- ph_test$table[1, "p"]
  
  # 保存检验结果
  ph_results <- rbind(ph_results, data.frame(Gene = gene, PH_pvalue = ph_p))
  
  # 提取残差数据用于绘图
  resid_vals <- as.numeric(ph_test$y[, 1])
  time_vals <- as.numeric(rownames(ph_test$y))
  residuals_df <- data.frame(time = time_vals, residuals = resid_vals)
  
  # 绘图（使用 geom_smooth 自动拟合 LOESS 曲线和置信区间）
  p <- ggplot(residuals_df, aes(x = time, y = residuals)) +
    geom_point(color = "#D62728", alpha = 0.6, size = 1.5) +
    geom_smooth(method = "loess", se = TRUE, color = "black", 
                fill = "gray70", alpha = 0.4, size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", alpha = 0.7) +
    labs(title = paste0(gene, " (p = ", round(ph_p, 4), ")"),
         x = "Time (days)", y = "Schoenfeld residuals") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          panel.grid.minor = element_blank())
  
  ph_plots[[gene]] <- p
}

# -------------------- 2. 保存 PH 检验结果表格 --------------------
ph_results$PH_pass <- ph_results$PH_pvalue > 0.05
write.csv(ph_results, "results/PH_assumption_test.csv", row.names = FALSE)
cat("PH 假设通过 (P>0.05) 的基因数量:", sum(ph_results$PH_pass), "\n")

# 合并单变量 Cox 结果与 PH 结果，保留 PH 通过的基因
cox_ph <- merge(cox_sig, ph_results, by = "Gene")
cox_ph <- cox_ph[cox_ph$PH_pass == TRUE, ]
cat("通过单变量 Cox 和 PH 假设的基因数量:", nrow(cox_ph), "\n")
write.csv(cox_ph, "results/univariate_Cox_PH_passed.csv", row.names = FALSE)

# -------------------- 3. 拼合图形并保存 --------------------
cat("\n🎨 拼合 Schoenfeld 残差图...\n")
n_plots <- length(ph_plots)
ncol <- 3
nrow <- ceiling(n_plots / ncol)
combined_plot <- wrap_plots(ph_plots, ncol = ncol, nrow = nrow) +
  plot_annotation(
    title = "Proportional Hazards Assumption Test (Schoenfeld Residuals)",
    subtitle = "Each panel shows residuals over time with fitted curve (black) and 95% CI (gray). p-value > 0.05 indicates PH assumption satisfied.",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 8))
  )

ggsave("results/Figure2B_PH_test.pdf", combined_plot, width = 12, height = 3 * nrow, dpi = 300)
ggsave("results/Figure2B_PH_test.png", combined_plot, width = 12, height = 3 * nrow, dpi = 300)

cat("\n✅ Figure 2B 已保存至 results/Figure2B_PH_test.pdf/png\n")

# -------------------- 6. LASSO 回归 (Figure 2C/D) --------------------
cat("\n🔧 LASSO 回归 (10折交叉验证)...\n")
# 准备 LASSO 输入: X = 表达矩阵 (基因 x 样本), y = 生存时间+状态
x_matrix <- t(expr_tumor[cox_ph$Gene, ])  # 样本 x 基因
y_surv <- Surv(surv_tumor$OS.time, surv_tumor$OS.event)

# 10 折交叉验证
set.seed(42)
cv_fit <- cv.glmnet(x_matrix, y_surv, family = "cox", alpha = 1, nfolds = 10)

# 文章指定的 lambda (log(lambda.min) = -4.2257)
lambda_article <- exp(-4.2257)
cat("文章指定的 lambda =", lambda_article, "\n")

# 使用文章 lambda 提取系数
coef_article <- coef(cv_fit, s = lambda_article)
selected_genes <- rownames(coef_article)[which(coef_article != 0)]
cat("LASSO 选中的基因 (文章 lambda):", length(selected_genes), "\n")
cat("选中的基因:", paste(selected_genes, collapse = ", "), "\n")

# 保存 LASSO 结果（包括交叉验证对象和选中的基因系数）
saveRDS(cv_fit, "results/LASSO_cv_fit.rds")
write.csv(data.frame(Gene = selected_genes, Coefficient = as.numeric(coef_article[which(coef_article != 0)])),
          "results/LASSO_selected_genes.csv", row.names = FALSE)

# -------------------- 绘制 LASSO 系数路径图 (Figure 2C) --------------------
# 文章指定的 lambda 值
lambda_article <- 0.0146151  # exp(-4.2257)

pdf("results/Figure2C_LASSO_path.pdf", width = 8, height = 6)
plot(cv_fit$glmnet.fit, xvar = "lambda", label = FALSE, main = "")  # 不显示默认标题
title(main = "LASSO Coefficient Path", line = 2)  # 标题上移2行
abline(v = log(lambda_article), col = "#D62728", lty = 2, lwd = 1.5)
# 添加图例（右上角）
legend("topright", legend = paste0("Article lambda = ", lambda_article), 
       col = "#D62728", lty = 2, bty = "n", cex = 0.8)
dev.off()

png("results/Figure2C_LASSO_path.png", width = 8, height = 6, units = "in", res = 300)
plot(cv_fit$glmnet.fit, xvar = "lambda", label = FALSE, main = "")
title(main = "LASSO Coefficient Path", line = 2)
abline(v = log(lambda_article), col = "#D62728", lty = 2, lwd = 1.5)
legend("topright", legend = paste0("Article lambda = ", lambda_article), 
       col = "#D62728", lty = 2, bty = "n", cex = 0.8)
dev.off()


# 绘制交叉验证曲线 (Figure 2D)，并标注文章 lambda 位置
# 文章指定的 lambda 值
lambda_article <- 0.0146151  # exp(-4.2257)

pdf("results/Figure2D_LASSO_cv.pdf", width = 8, height = 6)
plot(cv_fit, main = "")  # 不显示默认标题
title(main = "10-Fold Cross-Validation for LASSO", line = 2)  # 标题上移2行
abline(v = log(lambda_article), col = "#D62728", lty = 2, lwd = 1.5)
legend("topright", legend = paste0("Article lambda = ", lambda_article), 
       col = "#D62728", lty = 2, bty = "n", cex = 0.8)
dev.off()

png("results/Figure2D_LASSO_cv.png", width = 8, height = 6, units = "in", res = 300)
plot(cv_fit, main = "")
title(main = "10-Fold Cross-Validation for LASSO", line = 2)
abline(v = log(lambda_article), col = "#D62728", lty = 2, lwd = 1.5)
legend("topright", legend = paste0("Article lambda = ", lambda_article), 
       col = "#D62728", lty = 2, bty = "n", cex = 0.8)
dev.off()

# 最终预后基因
prognostic_genes <- selected_genes
cat("最终预后基因:", paste(prognostic_genes, collapse = ", "), "\n")
writeLines(prognostic_genes, "results/prognostic_genes.txt")
# -------------------- 7. 箱线图 (Figure 2E) --------------------
# -------------------- 7. 箱线图 (Figure 2E) --------------------
cat("\n📊 绘制 8 个预后基因的箱线图 (肿瘤 vs 正常)...\n")

# 使用所有样本的表达矩阵（包括正常）
expr_all <- expr_vst[prognostic_genes, , drop = FALSE]

# 从 dds 的 colData 中获取有行名的样本类型
sample_type_named <- colData(dds)$condition
names(sample_type_named) <- rownames(colData(dds))
sample_type_all <- sample_type_named[colnames(expr_all)]

# 构建数据框
df_box <- as.data.frame(t(expr_all))
df_box$Type <- sample_type_all
df_long <- pivot_longer(df_box, cols = -Type, names_to = "Gene", values_to = "Expression")
df_long$Type <- factor(df_long$Type, levels = c("Normal", "Tumor"))

# 计算 Wilcoxon 检验 p 值，并只保留显著的基因
p_values <- data.frame()
for (gene in prognostic_genes) {
  sub <- df_long[df_long$Gene == gene, ]
  test <- wilcox.test(Expression ~ Type, data = sub)
  p_val <- test$p.value
  if (p_val < 0.05) {  # 只标记显著的
    p_sig <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", "*"))
    p_values <- rbind(p_values, data.frame(Gene = gene, sig = p_sig))
  }
}

# 获取每个基因的最大表达量
max_expr <- df_long %>%
  group_by(Gene) %>%
  summarise(max_exp = max(Expression, na.rm = TRUE))

# 合并标签数据：只保留有显著性的基因
label_df <- merge(p_values, max_expr, by = "Gene")
if (nrow(label_df) > 0) {
  label_df$y_pos <- label_df$max_exp + label_df$max_exp * 0.05
  label_df$x_pos <- 1.5  # 两组之间的位置
}

# 箱线图（基础部分不变）
boxplot <- ggplot(df_long, aes(x = Type, y = Expression, fill = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = Type), width = 0.2, size = 0.8, alpha = 0.4) +
  scale_fill_manual(values = c("Normal" = "#1F77B4", "Tumor" = "#D62728")) +
  scale_color_manual(values = c("Normal" = "#1F77B4", "Tumor" = "#D62728")) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 2) +
  labs(title = "Expression of 8 Prognostic Genes in TCGA-STAD",
       x = "", y = "VST normalized expression") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# 如果有显著的基因，添加标记（关键：inherit.aes = FALSE）
if (nrow(label_df) > 0) {
  boxplot <- boxplot + 
    geom_text(data = label_df, aes(x = x_pos, y = y_pos, label = sig),
              inherit.aes = FALSE, size = 5, fontface = "bold")
}

ggsave("results/Figure2E_boxplot.pdf", boxplot, width = 12, height = 8)
ggsave("results/Figure2E_boxplot.png", boxplot, width = 12, height = 8, dpi = 300)
# -------------------- 8. 相关性热图 (Figure 2F) --------------------
cat("\n🔗 计算预后基因在肿瘤样本中的 Spearman 相关性...\n")
expr_prog <- expr_tumor[prognostic_genes, ]
cor_matrix <- cor(t(expr_prog), method = "spearman")
# 计算 p 值矩阵（需要 corrplot 包）
if (!require(corrplot)) install.packages("corrplot")
library(corrplot)
cor_p <- cor.mtest(t(expr_prog), method = "spearman")$p

# 热图颜色
col_fun <- colorRamp2(c(-1, 0, 1), c("#1F77B4", "white", "#D62728"))

# 创建热图对象（使用 column_title 代替 main）
heatmap <- Heatmap(cor_matrix, name = "Spearman ρ", col = col_fun,
                   cluster_rows = TRUE, cluster_columns = TRUE,
                   row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                   column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if (abs(cor_matrix[i, j]) > 0.3 & cor_p[i, j] < 0.05) {
                       grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y, 
                                 gp = gpar(fontsize = 8))
                     }
                   },
                   column_title = "Spearman Correlation among Prognostic Genes",
                   column_title_gp = gpar(fontsize = 12, fontface = "bold"))

# 保存
pdf("results/Figure2F_correlation_heatmap.pdf", width = 8, height = 7)
draw(heatmap)
dev.off()
png("results/Figure2F_correlation_heatmap.png", width = 8, height = 7, units = "in", res = 300)
draw(heatmap)
dev.off()

# -------------------- 9. 保存所有结果 --------------------
write.csv(cor_matrix, "results/prognostic_genes_correlation.csv")
write.csv(as.data.frame(cor_p), "results/prognostic_genes_correlation_pvalues.csv")
save.image("results/Day2_workspace.RData")

cat("\n🎉 Day2 全部完成！结果保存在 results/ 目录。\n")
