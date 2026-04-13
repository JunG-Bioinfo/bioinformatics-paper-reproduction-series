#!/usr/bin/env Rscript
# =============================================================================
# 📚 迁移体相关胃癌生信文献复刻 | Day1 (全部样本，精确复刻 Methods)
# 复刻 Figure 1:
#   A: 火山图 (标注 top10 up/down)
#   B: 热图 (top DEGs 表达模式)
#   C: 韦恩图 (DEGs ∩ MRGs)
#   D: GO 富集气泡图 (BP/CC/MF 各 top5)
#   E: KEGG 富集气泡图 (top5)
#   F: PPI 网络 (STRING 导出，用于 Cytoscape)
# 方法: DESeq2 (|log2FC|>0.5, P<0.05)
# 数据: TCGA-STAD (实际 375 tumor + 41 normal)
# 作者: 郭俊
# 日期: 2026-04-11
# =============================================================================

rm(list = ls())
PROJECT_DIR <- "/data/nas1/guojun_OD/project/05_paper_reproduction/Migrasome_STAD_Day1_AllSamples"
if (!dir.exists(PROJECT_DIR)) dir.create(PROJECT_DIR, recursive = TRUE)
setwd(PROJECT_DIR)
for (d in c("data", "results", "scripts", "PPI")) {
  if (!dir.exists(d)) dir.create(d)
}

# -------------------- 1. 加载包 --------------------
required_packages <- c(
  "DESeq2", "tidyverse", "ggplot2", "ggrepel", "ComplexHeatmap", 
  "circlize", "ggvenn", "clusterProfiler", "org.Hs.eg.db", 
  "enrichplot", "STRINGdb", "igraph"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("DESeq2", "ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db", "STRINGdb")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# -------------------- 2. 定义 MRGs (35个) --------------------
migrasome_genes <- c(
  "BMP1", "CPQ", "PDGFD", "TSPAN5", "TSPAN7", "TGFB2", "WNT11", "LEFTY1",
  "TSPAN4", "TSPAN9", "TSPAN14", "TSPAN15", "TSPAN18", "TSPAN31",
  "ITGA5", "ITGAV", "ITGB1", "ITGB3", "ITGB5",
  "NDST1", "NDST2", "NDST3", "NDST4",
  "PIGK", "PIGO", "PIGS", "PIGT", "PIGU", "PIGV", "PIGW", "PIGX", "PIGY", "PIGZ",
  "EVA1A", "EVA1B", "EVA1C"
)

# -------------------- 3. 读取已有数据（全部样本）--------------------
cat("\n📂 读取已有 TCGA-STAD 数据（全部样本）...\n")
old_raw <- readRDS("../Migrasome_STAD_Day1/data/TCGA_STAD_raw.rds")
counts_full <- old_raw$counts
sample_type_full <- old_raw$sample_type

# 直接使用所有有效样本（无筛选）
valid <- !is.na(sample_type_full)
counts <- counts_full[, valid]
sample_type <- sample_type_full[valid]
cat("使用样本: 肿瘤", sum(sample_type == "Tumor"), "正常", sum(sample_type == "Normal"), "\n")

# 保存整理后的数据
saveRDS(list(counts = counts, sample_type = sample_type), file = "data/STAD_all_samples.rds")

# -------------------- 4. DESeq2 差异分析 --------------------
cat("\n🔬 运行 DESeq2 (|log2FC|>0.5, P<0.05)...\n")
colData <- data.frame(condition = factor(sample_type, levels = c("Normal", "Tumor")))
rownames(colData) <- colnames(counts)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
# 保持默认过滤（至少3个样本count>10）
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# 标记上下调 (|log2FC|>0.5 且 P<0.05)
res_df$regulation <- "Not significant"
idx_up <- which(res_df$log2FoldChange > 0.5 & res_df$pvalue < 0.05)
idx_down <- which(res_df$log2FoldChange < -0.5 & res_df$pvalue < 0.05)
res_df$regulation[idx_up] <- "Up"
res_df$regulation[idx_down] <- "Down"

deg_up <- res_df$gene[idx_up]
deg_down <- res_df$gene[idx_down]
cat("上调基因:", length(deg_up), "\n下调基因:", length(deg_down), 
    "\n总DEGs:", length(deg_up)+length(deg_down), "\n")
write.csv(res_df, "results/STAD_DESeq2_all_genes.csv", row.names = FALSE)

# -------------------- 5. 火山图 (Figure 1A) 标注 top10 up/down --------------------
top_up <- res_df[idx_up, ] %>% arrange(desc(log2FoldChange)) %>% head(10) %>% pull(gene)
top_down <- res_df[idx_down, ] %>% arrange(log2FoldChange) %>% head(10) %>% pull(gene)
top20_genes <- c(top_up, top_down)
label_df <- res_df[res_df$gene %in% top20_genes, ]

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Up" = "#D62728", "Down" = "#1F77B4", "Not significant" = "#D3D3D3")) +
  geom_point(data = label_df, aes(x = log2FoldChange, y = -log10(pvalue)), 
             color = "#FFD700", size = 2.5, shape = 18) +
  geom_text_repel(data = label_df, aes(label = gene), size = 3.5, fontface = "bold",
                  box.padding = 0.5, max.overlaps = 20) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  labs(title = "TCGA-STAD: Differential Expression (DESeq2)",
       subtitle = paste0("Up: ", length(deg_up), ", Down: ", length(deg_down), 
                         " | |log2FC|>0.5, P<0.05\nGold: top10 up/down"),
       x = expression(log[2] ~ "Fold Change"), y = expression(-log[10] ~ "P-value")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9))
ggsave("results/Figure1A_volcano.pdf", volcano, width = 9, height = 7)
ggsave("results/Figure1A_volcano.png", volcano, width = 9, height = 7, dpi = 300)

# -------------------- 6. 热图 (Figure 1B) 使用 ComplexHeatmap --------------------
# 取 top 20 DEGs (按 pvalue 排序)
top20_deg <- res_df %>% filter(regulation != "Not significant") %>% 
  arrange(pvalue) %>% head(20) %>% pull(gene)
# 获取 VST 表达矩阵
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[rownames(assay(vsd)) %in% top20_deg, ]
# 按样本类型排序
sample_order <- order(sample_type)
mat <- mat[, sample_order]
annotation_col <- data.frame(Type = sample_type[sample_order])
rownames(annotation_col) <- colnames(mat)

# 动态设置颜色断点（使用数据的5%, 50%, 95%分位数）
library(circlize)
col_fun <- colorRamp2(
  breaks = quantile(mat, c(0.05, 0.5, 0.95), na.rm = TRUE),
  colors = c("#1F77B4", "white", "#D62728")
)

# 绘制热图
pdf("results/Figure1B_heatmap.pdf", width = 10, height = 10)
draw(Heatmap(mat, name = "Expression", col = col_fun, 
             cluster_rows = TRUE, cluster_columns = FALSE,
             show_column_names = FALSE, row_names_gp = gpar(fontsize = 8),
             column_title = "Top 20 DEGs in TCGA-STAD",
             top_annotation = HeatmapAnnotation(Type = sample_type[sample_order],
                                                col = list(Type = c("Tumor" = "#D62728", "Normal" = "#1F77B4")),
                                                annotation_name_gp = gpar(fontsize = 10))))
dev.off()

png("results/Figure1B_heatmap.png", width = 10, height = 10, units = "in", res = 300)
draw(Heatmap(mat, name = "Expression", col = col_fun, 
             cluster_rows = TRUE, cluster_columns = FALSE,
             show_column_names = FALSE, row_names_gp = gpar(fontsize = 8),
             column_title = "Top 20 DEGs in TCGA-STAD",
             top_annotation = HeatmapAnnotation(Type = sample_type[sample_order],
                                                col = list(Type = c("Tumor" = "#D62728", "Normal" = "#1F77B4")),
                                                annotation_name_gp = gpar(fontsize = 10))))
dev.off()
# -------------------- 7. 韦恩图 (Figure 1C) --------------------
deg_genes <- res_df$gene[res_df$regulation != "Not significant"]
candidate_genes <- intersect(deg_genes, migrasome_genes)
cat("候选基因 (DEGs ∩ MRGs):", length(candidate_genes), "\n")
writeLines(candidate_genes, "results/candidate_genes.txt")

venn <- ggvenn(
  list(DEGs = deg_genes, MRGs = migrasome_genes),
  fill_color = c("#D62728", "#1F77B4"), stroke_size = 0.5, set_name_size = 4,
  text_size = 3
) + ggtitle("Overlap between DEGs and MRGs") + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/Figure1C_venn.pdf", venn, width = 6, height = 6)
ggsave("results/Figure1C_venn.png", venn, width = 6, height = 6, dpi = 300)

# -------------------- 8. GO 富集分析 (Figure 1D) --------------------
cat("\n🧬 GO 富集分析...\n")
gene_entrez <- bitr(candidate_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- gene_entrez$ENTREZID
go_bp <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_cc <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_mf <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# 取各分类 top5
plot_go <- function(go_obj, title) {
  if (is.null(go_obj) || nrow(go_obj) == 0) return(NULL)
  df <- as.data.frame(go_obj)[1:min(5, nrow(go_obj)), ]
  ggplot(df, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
    geom_bar(stat = "identity") + scale_fill_gradient(low = "#D62728", high = "#1F77B4") +
    labs(x = "Gene Count", y = "", title = title) + theme_minimal()
}
p_bp <- plot_go(go_bp, "GO Biological Process")
p_cc <- plot_go(go_cc, "GO Cellular Component")
p_mf <- plot_go(go_mf, "GO Molecular Function")

if (!is.null(p_bp)) ggsave("results/Figure1D_GO_BP.pdf", p_bp, width = 8, height = 5)
if (!is.null(p_cc)) ggsave("results/Figure1D_GO_CC.pdf", p_cc, width = 8, height = 5)
if (!is.null(p_mf)) ggsave("results/Figure1D_GO_MF.pdf", p_mf, width = 8, height = 5)

# -------------------- 9. KEGG 富集分析 (Figure 1E) --------------------
kegg <- enrichKEGG(gene = entrez_ids, organism = "hsa", pvalueCutoff = 0.05)
if (!is.null(kegg) && nrow(kegg) > 0) {
  kegg_top5 <- as.data.frame(kegg)[1:min(5, nrow(kegg)), ]
  p_kegg <- ggplot(kegg_top5, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
    geom_bar(stat = "identity") + scale_fill_gradient(low = "#D62728", high = "#1F77B4") +
    labs(x = "Gene Count", y = "", title = "KEGG Pathways") + theme_minimal()
  ggsave("results/Figure1E_KEGG.pdf", p_kegg, width = 8, height = 5)
  ggsave("results/Figure1E_KEGG.png", p_kegg, width = 8, height = 5, dpi = 300)
} else {
  cat("⚠️ KEGG 无显著富集结果。\n")
}

# -------------------- 10. PPI 网络 (Figure 1F) 导出 STRING 数据 --------------------
cat("\n🔗 构建 PPI 网络 (STRING)...\n")
# 使用 STRINGdb 包查询
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400, input_directory = "PPI")
# 映射基因
genes_mapped <- string_db$map(data.frame(gene = candidate_genes), "gene", removeUnmappedRows = TRUE)

if (nrow(genes_mapped) > 1) {
  # 获取相互作用
  interactions <- string_db$get_interactions(genes_mapped$STRING_id)
  # 保存网络文件用于 Cytoscape
  write.csv(interactions, "PPI/STRING_interactions.csv", row.names = FALSE)
  writeLines(candidate_genes, "PPI/candidate_genes_for_Cytoscape.txt")
  cat("PPI 网络文件已保存至 PPI/ 目录，可用 Cytoscape 打开。\n")
  # 可选：在 R 中绘制简单网络图
  if (nrow(interactions) > 0) {
    g <- graph_from_data_frame(interactions[, c("from", "to")], directed = FALSE)
    pdf("results/Figure1F_PPI_network.pdf", width = 8, height = 8)
    plot(g, vertex.label = V(g)$name, vertex.size = 10, vertex.label.cex = 0.8, 
         main = "PPI Network of Candidate Genes")
    dev.off()
    png("results/Figure1F_PPI_network.png", width = 8, height = 8, units = "in", res = 300)
    plot(g, vertex.label = V(g)$name, vertex.size = 10, vertex.label.cex = 0.8,
         main = "PPI Network of Candidate Genes")
    dev.off()
  }
} else {
  cat("⚠️ 候选基因在 STRING 中映射不足，无法构建网络。\n")
}

# -------------------- 11. 保存候选基因列表及 session info --------------------
write.csv(data.frame(Gene = candidate_genes), "results/candidate_gene_list.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
cat("\n🎉 Day1 全部完成！结果保存在 results/ 和 PPI/ 目录。\n")
