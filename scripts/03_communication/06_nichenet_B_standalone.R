# =============================================================================
# TNBC interface-niche code repository
# Original source: FigureAI/Figure5/scripts/04_nichenet_B_standalone.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

#!/usr/bin/env Rscript
# =============================================================================
# 独立图 fig_nichenet_B_activity_with_targets.{pdf,png}
# 在原图（34 配体 × 3 分区 + 25 靶基因）的基础上，删除 9 个纤维/间质结构基因，
# 靶基因侧 25 -> 16 列；其余（配体行集、分区列、配色、字号、图例、画布尺寸）不变。
#
#   删除：COL1A1, COL1A2, COL3A1, ACTA2, LMNA, COL12A1, POSTN, AEBP1, CDH11
#   理由：这些基因与 MMP11 CAF 丰度高度相关，属发送方转录本混入，不是受体细胞响应
#
# 关键约束：相对权重仍按**删列前 25 列**内该配体的最大值＝1 计算，只删列、不二次归一化，
#           因此保留的每一格数值（因而颜色）与删列前完全一致。
#
# 运行：/path/to/envs/cellchat/bin/Rscript 04_nichenet_B_standalone.R
# =============================================================================
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(grid); library(dplyr); library(tidyr)
})

NICHE <- "/path/to/project/project_spatial/FigureAI/results/nichenet_C1Q_MMP11CAF"
DROP <- c("COL1A1","COL1A2","COL3A1","ACTA2","LMNA",   # 胶原/间质结构
           "COL12A1","POSTN","AEBP1","CDH11")           # 其余纤维/间质标志基因

act <- read.csv(file.path(NICHE,"tables/ligand_activity_B_MMP11_to_C1Q.csv"), stringsAsFactors = FALSE)
ligs <- act %>% group_by(zone) %>% arrange(desc(aupr)) %>% slice_head(n = 18) %>%
  ungroup() %>% distinct(test_ligand) %>% pull(test_ligand)
mat_act <- act %>% select(test_ligand, zone, aupr) %>%
  pivot_wider(names_from = zone, values_from = aupr) %>%
  column_to_rownames("test_ligand") %>% as.matrix()
mat_act <- mat_act[ligs, c("Tumor Core","Interface","Stroma"), drop = FALSE]

lt <- read.csv(file.path(NICHE,"tables/ligand_target_receiver_top.csv"), stringsAsFactors = FALSE)
freq <- tapply(lt$target_frequency, lt$target, max)
targets25 <- names(sort(freq, decreasing = TRUE))[1:25]
w <- lt[lt$ligand %in% ligs & lt$target %in% targets25, c("ligand","target","raw_weight")]
mat <- matrix(0, nrow = length(ligs), ncol = length(targets25), dimnames = list(ligs, targets25))
mat[cbind(match(w$ligand, ligs), match(w$target, targets25))] <- w$raw_weight
rmax <- apply(mat, 1, max); rmax[rmax <= 0] <- 1
mat25 <- mat / rmax                                        # 原 25 列基准
stopifnot(all(DROP %in% colnames(mat25)))
matK <- mat25[, setdiff(colnames(mat25), DROP), drop = FALSE]   # 只删列

col_act <- colorRamp2(seq(0.005, max(mat_act, na.rm = TRUE), length.out = 3),
                      c("#F7F7F7", "#E8789A", "#C51B7D"))
col_tgt <- colorRamp2(c(0, 0.05, 0.3, 1), c("#EDEDED", "#FDDBC7", "#E8789A", "#B2182B"))

ht_act <- Heatmap(mat_act, name = "AUPR", col = col_act, na_col = "grey90",
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_order = seq_along(ligs), column_order = seq_len(ncol(mat_act)),
  row_names_side = "left", row_names_gp = gpar(fontsize = 6.5),
  column_names_rot = 30, column_names_gp = gpar(fontsize = 8),
  column_title = "Ligand activity (AUPR)",
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  width = unit(3.2, "cm"), rect_gp = gpar(col = "white", lwd = 0.5),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 7)))

ht_tgt <- Heatmap(matK, name = "Relative target weight", col = col_tgt,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_order = seq_along(ligs), column_order = seq_len(ncol(matK)),
  show_row_names = FALSE, column_names_rot = 90, column_names_gp = gpar(fontsize = 6.5),
  column_title = "Predicted targets (collagen/stromal genes removed)",
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  width = unit(13.5, "cm"), rect_gp = gpar(col = "white", lwd = 0.4),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 7)))

ht_list <- ht_act + ht_tgt
ttl <- "MMP11 CAF -> C1Q TAM: ligand activity and predicted targets"
for (dev in c("png","pdf")) {
  f <- file.path(NICHE, paste0("figures/fig_nichenet_B_activity_with_targets.", dev))
  if (dev == "png") png(f, width = 10.5, height = 9.5, units = "in", res = 600, type = "cairo", bg = "white")
  else cairo_pdf(f, width = 10.5, height = 9.5)
  suppressWarnings(draw(ht_list, ht_gap = unit(4, "mm"), column_title = ttl,
       column_title_gp = gpar(fontsize = 11, fontface = "bold"),
       padding = unit(c(4, 4, 4, 4), "mm")))
  dev.off()
}

write.csv(as.data.frame(as.table(matK)) %>% rename(ligand = Var1, target = Var2, rel_weight = Freq),
          file.path(NICHE,"tables/ligand_target_plotted_after_drop.csv"), row.names = FALSE)
cat("[Done] standalone: ligands", length(ligs), "| targets", ncol(matK),
    "| dropped:", paste(DROP, collapse = ","), "\n")
