# 01 - 迁移体相关基因在胃癌中的预后模型复刻

> 原文：Identification of migrasome-related prognostic genes and construction of a prognostic model in gastric cancer. *Onco Targets and Therapy*, 2025.

---

## 📊 复刻进度（21天系列）

| Day | 分析内容 | 脚本 | 输出图 |
|:---:|---------|------|--------|
| 1 | 数据下载 + 差异表达 | `01_download_and_DEG.R` | 火山图、热图、箱线图 |
| 2 | GEO验证集整合 | `02_GEO_validation.R` | 一致性热图、相关性散点图 |
| 3 | LASSO特征筛选 | `03_LASSO.R` | 系数路径图、交叉验证图 |
| 4 | 多变量Cox回归 | `04_multivariate_Cox.R` | 森林图 |
| ... | ... | ... | ... |
| 21 | 完整总结 | `21_summary.Rmd` | 全部图表合集 |

（完整21天列表请见 [系列大纲](./SERIES_OUTLINE.md)）

---

## 🔧 运行环境

- R 4.3+
- 主要R包：TCGAbiolinks, tidyverse, limma, survival, glmnet, pheatmap, ggplot2

安装依赖：
```r
install.packages(c("tidyverse", "glmnet", "survival", "ggplot2", "pheatmap"))
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "limma", "survival"))
