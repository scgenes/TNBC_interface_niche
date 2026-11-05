# =============================================================================
# TNBC interface-niche code repository
# Original source: project_singlecell/TNBC_result/00_data_inventory/shared_style.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

## ============================================================
## shared_style.R — TNBC 肿瘤免疫边界 统一配色/主题
## 所有 figure 脚本 source 本文件，保证跨图故事线一致
## ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

## ---- 大群 9 类配色（Fig1/2/3/6 复用）----
major_colors <- c(
  "Tumor Epithelial cells" = "#C51B7D",
  "T cells"                 = "#4DAF4A",
  "Macrophages"             = "#FF7F00",
  "Fibroblasts"             = "#377EB8",
  "Plasma cells"            = "#984EA3",
  "B cells"                 = "#FFD92F",
  "Pericytes"               = "#A65628",
  "Endothelial cells"       = "#F781BF",
  "Normal Epithelial cells" = "#66C2A5"
)

## ---- 肿瘤-免疫边界三区配色（Fig2-6 复用）----
region_colors <- c(
  "Tumor Core" = "#C51B7D",
  "Interface"  = "#4DAF4A",
  "Stroma"     = "#377EB8"
)

## ---- inferCNV 恶性/非恶性 ----
malignant_colors <- c("Malignant" = "#E41A1C", "Non-malignant" = "#D9D9D9")

## ---- 统一 UMAP 主题（Nature/Cell 级）----
umap_theme <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(2.5, "mm"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0))

## ---- 表达热图梯度 ----
expr_gradient <- scale_color_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027")

## ---- 亚群 35 色（按谱系色系扩展）----
subtype_palette <- function(subs) {
  setNames(scales::hue_pal(l = 55, c = 90)(length(subs)), subs)
}
