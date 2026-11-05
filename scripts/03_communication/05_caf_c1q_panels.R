# =============================================================================
# TNBC interface-niche code repository
# Original source: FigureAI/Figure5/scripts/01_panels_abc.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

#!/usr/bin/env Rscript
# =============================================================================
# Figure 5 / Figure S5 —— 面板 a、b、c 重绘（R）
#
# 前提：只读取既有结果表，不重算任何分析。
#   a = NicheNet 方向 B（MMP11 CAF -> C1Q TAM）配体活性 + 预测靶基因   98 x 101 mm（暖色）
#   b = CellChat 空间分区 20 对 LR 点图                                76 x 101 mm（双方向）
#   c = NicheNet 方向 A（C1Q TAM -> MMP11 CAF）配体活性 + 受体集靶基因 130 x 112 mm（冷色）
#
# 运行环境：/path/to/envs/cellchat/bin/Rscript 01_panels_abc.R
#   R 4.3.2 / ggplot2 3.5.1 / ComplexHeatmap 2.18.0 / circlize 0.4.16 / viridisLite 0.4.2
#
# 设计约束（与最终图一致）
#   * 所有文字 >= 6 pt
#   * 方向配色全图统一：暖 #E64B35 = MMP11 CAF -> C1Q TAM；冷 #3C5488 = C1Q TAM -> MMP11 CAF
#   * 连续色阶一律 viridis；色条只标实测两端（3 位有效数字，纯小数，不用科学计数）
#   * 面板 b 的方向表头竖排（与面板 a 的列标签同规则），点径 4 档（≤25/≤50/≤75/100%）
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(ComplexHeatmap); library(circlize); library(grid); library(viridisLite)
})

ROOT  <- "/path/to/project/project_spatial/FigureAI"
OUT5  <- file.path(ROOT, "Figure5")
NICHE <- file.path(ROOT, "results/nichenet_C1Q_MMP11CAF")
CCR   <- file.path(ROOT, "results/cellchat_by_region")
dir.create(file.path(OUT5, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT5, "tables"), recursive = TRUE, showWarnings = FALSE)

NPG_COLORS <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000")
WARM <- NPG_COLORS[1]; COOL <- NPG_COLORS[4]
FORCED <- c("TGFB1","CSF1","APOE","SPP1")
FONT <- "sans"

act_B  <- read.csv(file.path(NICHE,"tables/ligand_activity_B_MMP11_to_C1Q.csv"), stringsAsFactors = FALSE)
lt_B   <- read.csv(file.path(NICHE,"tables/ligand_target_receiver_top.csv"), stringsAsFactors = FALSE)
act_A  <- read.csv(file.path(NICHE,"tables/ligand_activity_A_C1Q_to_MMP11.csv"), stringsAsFactors = FALSE)
lt_A   <- read.csv(file.path(NICHE,"tables/ligand_target_receiver_A.csv"), stringsAsFactors = FALSE)
sel20  <- read.csv(file.path(CCR,"tables/cellchat_LR_selected_for_plotting.csv"), stringsAsFactors = FALSE)
lr_lvl <- read.csv(file.path(CCR,"tables/cellchat_LR_level_by_region.csv"), stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# 面板 a / c 共用的双热图函数（保证两面板完全同规格）
# ---------------------------------------------------------------------------
build_nichenet_panel <- function(act, lt_tab, zones, forced, n_rows, n_tgt,
                                 dir_title, dir_col, w_act_mm, w_tgt_mm,
                                 out_pdf, out_png, W_mm, H_mm,
                                 drop_tgt = character(0), drop_lig = character(0)) {
  rank <- act[!is.na(act$aupr), ] %>% group_by(test_ligand) %>%
    summarise(max_aupr = max(aupr), .groups = "drop") %>% arrange(desc(max_aupr))
  all_ligs <- rank$test_ligand
  sel <- head(all_ligs, n_rows)
  add <- setdiff(forced, sel); if (length(add) > 0) sel <- c(head(sel, n_rows - length(add)), add)
  sel <- all_ligs[all_ligs %in% sel]
  stopifnot(length(sel) == n_rows, all(forced %in% sel))
  if (length(drop_lig) > 0) {   # 整行删除指定配体（不影响其它行数值）
    stopifnot(all(drop_lig %in% sel))
    sel <- sel[!sel %in% drop_lig]
  }

  mat_act <- act[act$test_ligand %in% sel, c("test_ligand","zone","aupr")] %>%
    pivot_wider(names_from = zone, values_from = aupr) %>%
    column_to_rownames("test_ligand") %>% as.matrix()
  mat_act <- mat_act[sel, zones, drop = FALSE]

  tgt_rank <- lt_tab %>% group_by(target) %>%
    summarise(freq = max(target_frequency), sw = sum(raw_weight), .groups = "drop") %>%
    arrange(desc(freq), desc(sw))
  tgt <- head(tgt_rank$target, n_tgt)
  w <- lt_tab[lt_tab$ligand %in% sel & lt_tab$target %in% tgt, c("ligand","target","raw_weight")]
  mat_tgt <- matrix(0, nrow = length(sel), ncol = length(tgt), dimnames = list(sel, tgt))
  mat_tgt[cbind(match(w$ligand, sel), match(w$target, tgt))] <- w$raw_weight
  rmax <- apply(mat_tgt, 1, max); rmax[rmax <= 0] <- 1
  mat_tgt <- mat_tgt / rmax
  if (length(drop_tgt) > 0) {   # 只删列，不重新归一化（保留格子的数值与颜色不变）
    stopifnot(all(drop_tgt %in% colnames(mat_tgt)))
    mat_tgt <- mat_tgt[, setdiff(colnames(mat_tgt), drop_tgt), drop = FALSE]
  }

  lim <- signif(range(mat_act, na.rm = TRUE), 3)
  lab <- formatC(lim, format = "g", digits = 3, flag = "#")
  col_act <- colorRamp2(seq(lim[1], lim[2], length.out = 64), viridisLite::viridis(64))
  col_tgt <- colorRamp2(c(0, 1), viridisLite::viridis(2))

  ht_act <- Heatmap(mat_act, name = "AUPR", col = col_act, na_col = "grey90",
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_order = seq_along(sel), column_order = seq_along(zones),
    row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily = FONT),
    column_names_rot = 90, column_names_gp = gpar(fontsize = 6, fontfamily = FONT),
    column_title = "Ligand activity (AUPR)",
    column_title_gp = gpar(fontsize = 6.5, fontface = "bold", fontfamily = FONT),
    width = unit(w_act_mm, "mm"), rect_gp = gpar(col = "white", lwd = 0.35),
    heatmap_legend_param = list(title = "AUPR", at = lim, labels = lab,
      title_gp = gpar(fontsize = 6, fontface = "plain", fontfamily = FONT),
      labels_gp = gpar(fontsize = 6, fontfamily = FONT),
      legend_direction = "horizontal", legend_width = unit(28, "mm"),
      title_position = "topcenter", grid_height = unit(2.2, "mm")))

  ht_tgt <- Heatmap(mat_tgt, name = "Relative weight", col = col_tgt,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_order = seq_along(sel), column_order = seq_along(tgt), show_row_names = FALSE,
    column_names_rot = 90, column_names_gp = gpar(fontsize = 6, fontfamily = FONT),
    column_title = "Predicted target genes",
    column_title_gp = gpar(fontsize = 6.5, fontface = "bold", fontfamily = FONT),
    width = unit(w_tgt_mm, "mm"), rect_gp = gpar(col = "white", lwd = 0.35),
    heatmap_legend_param = list(title = "Relative weight", at = c(0, 0.5, 1),
      title_gp = gpar(fontsize = 6, fontface = "plain", fontfamily = FONT),
      labels_gp = gpar(fontsize = 6, fontfamily = FONT),
      legend_direction = "horizontal", legend_width = unit(28, "mm"),
      title_position = "topcenter", grid_height = unit(2.2, "mm")))

  ht <- ht_act + ht_tgt
  for (dev in c("pdf","png")) {
    f <- if (dev == "pdf") out_pdf else out_png
    if (dev == "pdf") cairo_pdf(f, width = W_mm/25.4, height = H_mm/25.4, family = FONT)
    else png(f, width = W_mm/25.4, height = H_mm/25.4, units = "in", res = 600, type = "cairo", bg = "white")
    # ComplexHeatmap 在内部度量时会触发无害的 mbcs 警告；落盘字体为 DejaVu Sans（见 README 校验）
    suppressWarnings(draw(ht, ht_gap = unit(3, "mm"), legend_gap = unit(12, "mm"),
         column_title = dir_title,
         column_title_gp = gpar(fontsize = 7, fontface = "bold", col = dir_col, fontfamily = FONT),
         heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
         padding = unit(c(4, 4, 4, 6), "mm")))
    dev.off()
  }
  list(mat_act = mat_act, mat_tgt = mat_tgt, ligs = sel, targets = tgt)
}

pan_a <- build_nichenet_panel(
  act = act_B, lt_tab = lt_B, zones = c("Tumor Core","Interface","Stroma"),
  forced = FORCED, n_rows = 20, n_tgt = 15,
  dir_title = "MMP11 CAF → C1Q TAM", dir_col = WARM,
  w_act_mm = 18, w_tgt_mm = 58,
  out_pdf = file.path(OUT5,"panels/panel_a.pdf"), out_png = file.path(OUT5,"panels/panel_a.png"),
  W_mm = 98, H_mm = 101,
  drop_tgt = c("COL1A1","COL1A2","COL3A1","ACTA2","LMNA"))

pan_c <- build_nichenet_panel(
  act = act_A, lt_tab = lt_A, zones = c("Interface","Stroma"),
  forced = FORCED, n_rows = 20, n_tgt = 15,
  dir_title = "C1Q TAM → MMP11 CAF", dir_col = COOL,
  w_act_mm = 16, w_tgt_mm = 80,
  out_pdf = file.path(OUT5,"panels/panel_c.pdf"), out_png = file.path(OUT5,"panels/panel_c.png"),
  W_mm = 130, H_mm = 112, drop_lig = "SPP1")

write.csv(as.data.frame(as.table(pan_a$mat_act)) %>% rename(ligand = Var1, zone = Var2, aupr = Freq),
          file.path(OUT5,"tables/panel_a_activity_plotted.csv"), row.names = FALSE)
write.csv(as.data.frame(as.table(pan_a$mat_tgt)) %>% rename(ligand = Var1, target = Var2, rel_weight = Freq),
          file.path(OUT5,"tables/panel_a_targets_plotted.csv"), row.names = FALSE)
write.csv(as.data.frame(as.table(pan_c$mat_act)) %>% rename(ligand = Var1, zone = Var2, aupr = Freq),
          file.path(OUT5,"tables/panel_c_activity_plotted.csv"), row.names = FALSE)
write.csv(as.data.frame(as.table(pan_c$mat_tgt)) %>% rename(ligand = Var1, target = Var2, rel_weight = Freq),
          file.path(OUT5,"tables/panel_c_targets_plotted.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 面板 b：20 对 LR 的分区点图（76 x 101 mm）
# ---------------------------------------------------------------------------
COL_ORDER <- c("Tumor Core | MMP11 CAF -> C1Q TAM", "Tumor Core | C1Q TAM -> MMP11 CAF",
               "Interface | MMP11 CAF -> C1Q TAM",  "Interface | C1Q TAM -> MMP11 CAF")
DIR_COL   <- c(WARM, COOL, WARM, COOL)
row20  <- sel20$LR
lvl20  <- data.frame(LR = row20, pathway = sel20$pathway, stringsAsFactors = FALSE)
lvl20$rank_top <- ave(seq_along(lvl20$LR), lvl20$pathway, FUN = seq_along)
lvl20$k        <- ave(seq_along(lvl20$LR), lvl20$pathway, FUN = length)
lvl20$y        <- length(row20) - seq_along(lvl20$LR) + 1

lr_lvl$key <- paste(lr_lvl$region, lr_lvl$direction, sep = " | ")
g20 <- expand.grid(LR = row20, col = COL_ORDER, stringsAsFactors = FALSE)
g20 <- merge(g20, lr_lvl[, c("LR","key","median_prob","n_sig_units","n_units_tested")],
             by.x = c("LR","col"), by.y = c("LR","key"), all.x = TRUE, sort = FALSE)
g20$tested <- !is.na(g20$n_units_tested)
g20$repro  <- ifelse(g20$tested, g20$n_sig_units / g20$n_units_tested, NA)
stopifnot(nrow(g20) == 80, sum(g20$tested) == 74)

LVL  <- c("≤25%", "≤50%", "≤75%", "100%")
d20 <- merge(g20, lvl20[, c("LR","pathway","y","k","rank_top")], by = "LR")
d20$xi <- match(d20$col, COL_ORDER)
d20$log10p <- ifelse(d20$tested & d20$median_prob > 0, log10(d20$median_prob), NA_real_)
d20$pct_bin <- factor(
  ifelse(d20$repro <= 0.25, "≤25%",
  ifelse(d20$repro <= 0.50, "≤50%",
  ifelse(d20$repro <= 0.75, "≤75%", "100%"))), levels = LVL)
dd_pos  <- d20[d20$tested & d20$median_prob > 0, ]
dd_zero <- d20[d20$tested & d20$median_prob == 0, ]
blk <- do.call(rbind, lapply(split(lvl20, lvl20$pathway), function(s)
  data.frame(pathway = s$pathway[1], y_top = max(s$y), y_bot = min(s$y), k = nrow(s))))
blk <- blk[match(unique(lvl20$pathway), blk$pathway), ]
seps <- blk$y_bot[-nrow(blk)] - 0.5
ytop <- max(lvl20$y); ybot <- min(lvl20$y)
HDR_DIR <- ytop + 2.45; HDR_REG <- ytop + 4.95; YTOP <- ytop + 5.40
LIMS <- range(d20$log10p, na.rm = TRUE)
ARSZ <- setNames(c(0.75, 1.10, 1.45, 1.80), LVL)
FLW <- 0.1175; FSW <- 0.282; TJ <- 2.11      # 0.33 pt 细线 / 0.8 pt 点描边 / 6.0 pt 字号

p_dot <- ggplot() +
  geom_segment(aes(x = 0.45, xend = 0.45, y = ybot - 0.5, yend = ytop + 0.5), colour = "grey25", linewidth = FLW) +
  geom_segment(aes(x = 0.45, xend = 4.55, y = ybot - 0.5, yend = ybot - 0.5), colour = "grey25", linewidth = FLW) +
  geom_segment(data = data.frame(y = seps), aes(x = 0.45, xend = 4.55, y = y, yend = y),
               colour = "grey82", linewidth = FLW) +
  geom_point(data = dd_pos, aes(x = xi, y = y, size = pct_bin, fill = log10p),
             shape = 21, colour = "white", stroke = FSW) +
  geom_point(data = dd_zero, aes(x = xi, y = y, size = pct_bin), shape = 21,
             fill = "grey60", colour = "white", stroke = FSW) +
  annotate("text", x = 0.36, y = blk$y_top, label = blk$pathway, hjust = 1, vjust = 0.5,
           size = TJ, colour = "grey10", family = "sans") +
  annotate("text", x = 1:4, y = HDR_DIR, label = c("CAF→C1Q","C1Q→CAF","CAF→C1Q","C1Q→CAF"),
           size = TJ, colour = DIR_COL, vjust = 0.5, angle = 90, family = "sans") +
  annotate("text", x = c(1.5, 3.5), y = HDR_REG, label = c("Core", "Interface"),
           size = TJ, colour = "black", vjust = 0.5, family = "sans") +
  scale_fill_viridis_c(option = "viridis", name = "log10(median prob)",
                       limits = LIMS, oob = scales::squish, na.value = "grey50",
                       breaks = c(-3, -2), labels = c("1e-3", "1e-2"),
                       guide = guide_colourbar(title.position = "top", barheight = unit(2.6, "mm"),
                                               barwidth = unit(15, "mm"), frame.colour = NA,
                                               ticks.colour = "grey20", ticks.linewidth = 0.1175,
                                               display = "gradient",
                                               title.theme = element_text(size = 6.0, family = "sans"),
                                               label.theme = element_text(size = 6.0, family = "sans"))) +
  scale_size_manual(values = ARSZ, name = "Reproduced / tested units",
                    breaks = LVL, drop = FALSE,
                    guide = guide_legend(order = 1, title.position = "top", nrow = 1,
                                         override.aes = list(fill = "grey35"),
                                         title.theme = element_text(size = 6.0, family = "sans"),
                                         label.theme = element_text(size = 6.0, family = "sans"))) +
  scale_x_continuous(limits = c(0.20, 4.80), expand = c(0, 0)) +
  scale_y_continuous(limits = c(ybot - 0.6, YTOP), expand = c(0, 0),
                     breaks = lvl20$y, labels = lvl20$LR) +
  coord_cartesian(xlim = c(0.45, 4.55), ylim = c(ybot - 0.5, YTOP), clip = "off") +
  labs(title = "Ligand-receptor pairs", x = NULL, y = NULL) +
  theme_void(base_size = 6.0) +
  theme(text = element_text(family = "sans", size = 6.0),
        axis.text.y = element_text(size = 6.0, colour = "black", hjust = 1, margin = margin(r = 38)),
        axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank(), plot.background = element_blank(),
        plot.title = element_text(size = 7.0, face = "bold", hjust = 0.5, family = "sans", margin = margin(b = 3)),
        plot.title.position = "plot",
        legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal",
        legend.justification = "right", legend.box.just = "right",
        legend.background = element_blank(), legend.key = element_blank(),
        legend.key.size = unit(3.2, "mm"), legend.spacing.y = unit(1.0, "mm"),
        legend.spacing.x = unit(1.6, "mm"), legend.margin = margin(1, 1, 0, 1, unit = "mm"),
        plot.margin = margin(5, 2.0, 2, 1.5, unit = "mm"))
W <- 76/25.4; H <- 101/25.4
ggsave(file.path(OUT5,"panels/panel_b.png"), p_dot, width = W, height = H, dpi = 600, bg = "white")
ggsave(file.path(OUT5,"panels/panel_b.pdf"), p_dot, width = W, height = H, bg = "white", device = cairo_pdf)
write.csv(d20[, c("LR","pathway","col","n_units_tested","n_sig_units","repro","median_prob","log10p","pct_bin","y")],
          file.path(OUT5,"tables/panel_b_dotplot_plotted.csv"), row.names = FALSE)

cat("[Done] panel_a (", length(pan_a$ligs), "ligands x", length(pan_a$targets), "targets )",
    "| panel_b (", nrow(d20), "cells |", sum(d20$tested), "tested )",
    "| panel_c (", length(pan_c$ligs), "ligands )\n")
