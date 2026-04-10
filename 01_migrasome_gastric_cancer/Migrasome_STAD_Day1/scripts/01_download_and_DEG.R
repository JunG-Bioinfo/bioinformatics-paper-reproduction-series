#!/usr/bin/env Rscript
# ============================================================================
# 📚 迁移体相关胃癌生信文献复刻 | Day 1
# 目标：迁移体相关基因在胃癌中的差异表达分析
# 复刻文献：Migrasome-Related Prognostic Genes in Gastric Cancer (Onco Targets Ther, 2025)
# 输出：火山图 + 热图（SCI标准）
# 作者：郭俊
# 日期：2026-04-10
# ============================================================================

# -------------------- 0. 环境配置 --------------------
rm(list = ls())

# 设置工作目录（⚠️ 请修改为你自己的项目路径）
PROJECT_DIR <- "/data/nas1/guojun_OD/project/05_paper_reproduction/Migrasome_STAD_Day1"
if (!dir.exists(PROJECT_DIR)) dir.create(PROJECT_DIR, recursive = TRUE)
setwd(PROJECT_DIR)

# 创建目录
dirs <- c("data", "results", "scripts")
for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# -------------------- 1. 加载R包 --------------------
# 若缺失则自动安装
required_packages <- c(
  "TCGAbiolinks", "tidyverse", "ggplot2", "ggrepel", "pheatmap",
  "SummarizedExperiment", "limma", "edgeR", "RColorBrewer"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("TCGAbiolinks", "SummarizedExperiment")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

cat("✅ 所有R包加载完成\n")

# -------------------- 2. 定义关键变量 --------------------
# 目标基因集：35个迁移体相关基因（MRGs）
# 来源：文献Methods部分 "35 migrasome-related genes (MRGs)"
migrasome_genes <- c(
  "BMP1", "CPQ", "PDGFD", "TSPAN5", "TSPAN7", "TGFB2", "WNT11", "LEFTY1",
  "TSPAN4", "TSPAN9", "TSPAN14", "TSPAN15", "TSPAN18", "TSPAN31",
  "ITGA5", "ITGAV", "ITGB1", "ITGB3", "ITGB5",
  "NDST1", "NDST2", "NDST3", "NDST4",
  "PIGK", "PIGO", "PIGS", "PIGT", "PIGU", "PIGV", "PIGW", "PIGX", "PIGY", "PIGZ",
  "EVA1A", "EVA1B", "EVA1C"
)

# 已验证的预后基因（文献Figure 5/6结果）
validated_genes <- c("BMP1", "CPQ", "PDGFD", "TSPAN5", "TSPAN7", "TGFB2", "WNT11", "LEFTY1")

# 癌症项目
cancer_project <- "TCGA-STAD"

cat("🎯 目标：", length(migrasome_genes), "个迁移体相关基因\n")
cat("🔬 重点关注：", paste(validated_genes, collapse = ", "), "\n")

# -------------------- 3. 下载TCGA-STAD RNA-seq数据 --------------------
cat("\n📥 步骤1：下载TCGA-STAD RNA-seq数据...\n")

# 查询数据：原发性肿瘤 + 实体正常组织
query <- GDCquery(
  project = cancer_project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# 下载数据（使用缓存避免重复下载）
GDCdownload(query, method = "api", files.per.chunk = 10)

# 准备数据
data <- GDCprepare(query, save = FALSE, summarizedExperiment = TRUE)

cat("✅ 数据下载完成！\n")
cat("  总样本数：", ncol(data), "\n")
cat("  总基因数：", nrow(data), "\n")

# 获取表达矩阵（counts）
count_matrix <- assay(data)

# 获取基因名
gene_names <- rowData(data)$gene_name
rownames(count_matrix) <- make.unique(gene_names)

# 获取样本信息
sample_info <- colData(data)
sample_type <- sample_info$shortLetterCode
sample_type <- ifelse(sample_type == "TP", "Tumor", 
                      ifelse(sample_type == "NT", "Normal", NA))

# 保存原始数据
saveRDS(list(counts = count_matrix, sample_info = sample_info, 
             sample_type = sample_type), 
        file = file.path("data", "TCGA_STAD_raw.rds"))

cat("✅ 原始数据已保存至 data/TCGA_STAD_raw.rds\n")

# -------------------- 4. 数据预处理 --------------------
cat("\n🔧 步骤2：数据预处理...\n")

# 筛选有效样本
valid_idx <- !is.na(sample_type)
counts_filtered <- count_matrix[, valid_idx]
sample_type_filtered <- sample_type[valid_idx]

cat("  有效样本：", ncol(counts_filtered), 
    "（肿瘤：", sum(sample_type_filtered == "Tumor"),
    "，正常：", sum(sample_type_filtered == "Normal"), "）\n")

# 创建DGEList对象用于edgeR标准化
group <- factor(sample_type_filtered, levels = c("Normal", "Tumor"))
dge <- DGEList(counts = counts_filtered, group = group)

# 过滤低表达基因（CPM > 1 至少在50%样本中）
keep <- rowSums(cpm(dge) > 1) >= ncol(dge) * 0.5
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat("  过滤后基因数：", nrow(dge), "\n")

# TMM标准化
dge <- calcNormFactors(dge)

# 获取log2(CPM+1)表达矩阵（用于差异分析和可视化）
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# 提取迁移体相关基因的表达
migrasome_logcpm <- logcpm[rownames(logcpm) %in% migrasome_genes, , drop = FALSE]
cat("  迁移体相关基因在表达矩阵中存在：", nrow(migrasome_logcpm), "/", length(migrasome_genes), "\n")

# 保存标准化后的表达矩阵
write.csv(logcpm, file = file.path("data", "TCGA_STAD_logcpm.csv"))
write.csv(migrasome_logcpm, file = file.path("data", "migrasome_genes_logcpm.csv"))

# -------------------- 5. 差异表达分析（limma-voom） --------------------
cat("\n📊 步骤3：差异表达分析（肿瘤 vs 正常）...\n")

# voom转换
v <- voom(dge, plot = TRUE)
# 保存voom图
png(file.path("results", "voom_plot.png"), width = 1200, height = 1000, res = 150)
voom_plot <- plotMD(v, main = "voom: Mean-Difference Plot")
dev.off()

# 构建设计矩阵
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 线性拟合
fit <- lmFit(v, design)

# 对比矩阵：Tumor vs Normal
cont_matrix <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

# 提取差异结果
deg_results <- topTable(fit2, number = Inf, sort.by = "P", adjust.method = "fdr")

# 添加基因名
deg_results$gene <- rownames(deg_results)

# 定义显著差异阈值
deg_results$regulation <- "Not significant"
deg_results$regulation[deg_results$logFC > 1 & deg_results$adj.P.Val < 0.05] <- "Up"
deg_results$regulation[deg_results$logFC < -1 & deg_results$adj.P.Val < 0.05] <- "Down"

# 筛选迁移体相关基因的差异结果
migrasome_deg <- deg_results[deg_results$gene %in% migrasome_genes, ]

# 标记8个预后基因
deg_results$is_validated <- ifelse(deg_results$gene %in% validated_genes, "Validated", "Other")

# 保存差异结果
write.csv(deg_results, file = file.path("results", "STAD_all_genes_DEG.csv"), row.names = FALSE)
write.csv(migrasome_deg, file = file.path("results", "STAD_migrasome_genes_DEG.csv"), row.names = FALSE)

cat("✅ 差异分析完成！\n")
cat("  上调基因数：", sum(deg_results$regulation == "Up"), "\n")
cat("  下调基因数：", sum(deg_results$regulation == "Down"), "\n")
cat("  迁移体相关基因差异：", nrow(migrasome_deg), "\n")

# -------------------- 6. 火山图（SCI标准） --------------------
cat("\n🎨 步骤4：生成火山图...\n")

# 准备标签数据（只标记8个预后基因）
label_data <- deg_results[deg_results$gene %in% validated_genes, ]

# 火山图
volcano <- ggplot(deg_results, aes(x = logFC, y = -log10(P.Value), color = regulation)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(
    values = c("Up" = "#D62728", "Down" = "#1F77B4", "Not significant" = "#D3D3D3"),
    name = "Regulation"
  ) +
  # 标记8个预后基因
  geom_point(data = label_data, aes(x = logFC, y = -log10(P.Value)), 
             color = "#FFD700", size = 4, shape = 18, alpha = 0.9) +
  geom_text_repel(data = label_data, aes(x = logFC, y = -log10(P.Value), label = gene),
                  size = 4, fontface = "bold", color = "black",
                  box.padding = 0.5, point.padding = 0.3, max.overlaps = 15) +
  # 阈值线
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  # 标签
  labs(
    title = "Gastric Cancer (TCGA-STAD): Differential Expression of Migrasome-Related Genes",
    subtitle = paste0("Tumor (n=", sum(sample_type_filtered == "Tumor"), 
                      ") vs Normal (n=", sum(sample_type_filtered == "Normal"), ")"),
    x = expression(log[2] ~ "Fold Change"),
    y = expression(-log[10] ~ "P-value"),
    caption = "Red: Upregulated, Blue: Downregulated | Dashed lines: |logFC| = 1, P = 0.05\nGold diamonds: 8 validated prognostic genes (BMP1, CPQ, PDGFD, TSPAN5, TSPAN7, TGFB2, WNT11, LEFTY1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    plot.caption = element_text(size = 8, color = "gray50", hjust = 0),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.5)
  )

# 保存火山图（PDF + PNG）
ggsave(file.path("results", "STAD_migrasome_volcano.pdf"), 
       plot = volcano, width = 10, height = 8, device = "pdf", dpi = 300)
ggsave(file.path("results", "STAD_migrasome_volcano.png"), 
       plot = volcano, width = 10, height = 8, device = "png", dpi = 300)

cat("✅ 火山图已保存至 results/\n")

# -------------------- 7. 热图（迁移体相关基因） --------------------
cat("\n🎨 步骤5：生成热图...\n")

# 准备热图数据：迁移体相关基因在肿瘤vs正常中的表达
heatmap_data <- migrasome_logcpm

# 按样本类型排序
sample_order <- order(sample_type_filtered)
heatmap_data <- heatmap_data[, sample_order]
annotation_col <- data.frame(SampleType = sample_type_filtered[sample_order])
rownames(annotation_col) <- colnames(heatmap_data)

# 配色方案
ann_colors <- list(SampleType = c(Tumor = "#D62728", Normal = "#1F77B4"))
heatmap_colors <- colorRampPalette(c("#1F77B4", "white", "#D62728"))(100)

# 设置热图文件名
pdf_file <- file.path("results", "STAD_migrasome_genes_heatmap.pdf")
png_file <- file.path("results", "STAD_migrasome_genes_heatmap.png")

# 生成热图（PDF）
pheatmap(
  heatmap_data,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heatmap_colors,
  main = "Migrasome-Related Genes Expression in Gastric Cancer",
  fontsize_row = 8,
  fontsize_col = 6,
  show_colnames = FALSE,
  border_color = NA,
  filename = pdf_file,
  width = 10,
  height = 8
)

# 生成热图（PNG）
pheatmap(
  heatmap_data,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heatmap_colors,
  main = "Migrasome-Related Genes Expression in Gastric Cancer",
  fontsize_row = 8,
  fontsize_col = 6,
  show_colnames = FALSE,
  border_color = NA,
  filename = png_file,
  width = 10,
  height = 8,
  res = 300
)

cat("✅ 热图已保存至 results/\n")

# -------------------- 8. 8个预后基因的箱线图 --------------------
cat("\n🎨 步骤6：生成8个预后基因的箱线图...\n")

# 提取8个预后基因的表达数据
validated_expr <- logcpm[rownames(logcpm) %in% validated_genes, , drop = FALSE]
validated_df <- as.data.frame(t(validated_expr))
validated_df$SampleType <- sample_type_filtered

# 转换为长格式
validated_long <- validated_df %>%
  pivot_longer(cols = -SampleType, names_to = "Gene", values_to = "Expression")

# 箱线图
boxplot_validated <- ggplot(validated_long, aes(x = SampleType, y = Expression, fill = SampleType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = SampleType), width = 0.2, size = 0.8, alpha = 0.4) +
  scale_fill_manual(values = c(Tumor = "#D62728", Normal = "#1F77B4")) +
  scale_color_manual(values = c(Tumor = "#D62728", Normal = "#1F77B4")) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 2) +
  labs(
    title = "Expression of 8 Validated Prognostic Migrasome-Related Genes",
    subtitle = "TCGA-STAD: Tumor vs Normal",
    x = "", y = expression(log[2] ~ CPM ~ "(normalized)"),
    caption = "Wilcoxon test: * p < 0.05, ** p < 0.01, *** p < 0.001"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    plot.caption = element_text(size = 8, color = "gray50", hjust = 0),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95", color = NA),
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# 保存箱线图
ggsave(file.path("results", "STAD_validated_genes_boxplot.pdf"), 
       plot = boxplot_validated, width = 12, height = 8, device = "pdf", dpi = 300)
ggsave(file.path("results", "STAD_validated_genes_boxplot.png"), 
       plot = boxplot_validated, width = 12, height = 8, device = "png", dpi = 300)

cat("✅ 箱线图已保存至 results/\n")

# -------------------- 9. 统计分析：8个预后基因的Wilcoxon检验 --------------------
cat("\n📈 步骤7：统计分析...\n")

wilcox_results <- data.frame(
  Gene = character(),
  log2FC = numeric(),
  P_value = numeric(),
  Regulation = character(),
  stringsAsFactors = FALSE
)

for (gene in validated_genes) {
  if (gene %in% rownames(logcpm)) {
    tumor_expr <- logcpm[gene, sample_type_filtered == "Tumor"]
    normal_expr <- logcpm[gene, sample_type_filtered == "Normal"]
    
    # Wilcoxon秩和检验
    test_result <- wilcox.test(tumor_expr, normal_expr)
    log2fc <- median(tumor_expr) - median(normal_expr)
    
    regulation <- ifelse(log2fc > 0, "Up", ifelse(log2fc < 0, "Down", "No change"))
    
    wilcox_results <- rbind(wilcox_results, data.frame(
      Gene = gene,
      log2FC = round(log2fc, 3),
      P_value = formatC(test_result$p.value, format = "e", digits = 3),
      Regulation = regulation,
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(wilcox_results, file = file.path("results", "validated_genes_wilcox_test.csv"), row.names = FALSE)
cat("✅ 统计结果已保存至 results/validated_genes_wilcox_test.csv\n")

# -------------------- 10. 完成报告 --------------------
cat("\n" + paste(rep("=", 60), collapse = "") + "\n")
cat("🎉 Day1 分析完成！\n")
cat("\n📁 输出文件：\n")
cat("  - results/STAD_all_genes_DEG.csv              (全部基因差异结果)\n")
cat("  - results/STAD_migrasome_genes_DEG.csv        (迁移体相关基因差异结果)\n")
cat("  - results/STAD_migrasome_volcano.pdf/png      (火山图)\n")
cat("  - results/STAD_migrasome_genes_heatmap.pdf/png (热图)\n")
cat("  - results/STAD_validated_genes_boxplot.pdf/png (箱线图)\n")
cat("  - results/validated_genes_wilcox_test.csv     (统计结果)\n")
cat("\n📊 关键结果摘要：\n")
cat("  迁移体相关基因中上调：", sum(migrasome_deg$logFC > 0 & migrasome_deg$adj.P.Val < 0.05), "\n")
cat("  迁移体相关基因中下调：", sum(migrasome_deg$logFC < 0 & migrasome_deg$adj.P.Val < 0.05), "\n")
cat("  8个预后基因在肿瘤中的表达趋势：\n")
for (i in 1:nrow(wilcox_results)) {
  cat("    -", wilcox_results$Gene[i], ":", wilcox_results$Regulation[i], 
      "(log2FC =", wilcox_results$log2FC[i], ", P =", wilcox_results$P_value[i], ")\n")
}
cat(paste(rep("=", 60), collapse = "") + "\n")

# -------------------- 可选：GEO验证数据下载提示 --------------------
cat("\n💡 提示：Day2将整合GEO验证数据集（GSE84437、GSE15459等）\n")
cat("  如需提前准备，可手动下载：\n")
cat("  - GSE84437: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84437\n")
cat("  - GSE15459: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15459\n")