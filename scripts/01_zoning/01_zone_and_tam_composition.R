# =============================================================================
# TNBC interface-niche code repository
# Original source: project_singlecell/TNBC_result/_archive_legacy/figure2/analysis.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

## ============================================================
## Figure 2 — 空间去卷积 + 肿瘤免疫边界定义（6 panel）
## A 肿瘤比例空间图(12样本) | B 三区空间图(12样本) | C 肿瘤比例直方图+阈值
## D 每样本区域占比 | E 去卷积质量QC | F 距离肿瘤核心的边界距离梯度
## 只读 region_assignment.rds 已有三区，不重跑 RCTD
## ============================================================
source(file.path("scripts", "lib", "shared_style.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

set.seed(42)
OUT <- "/path/to/project/project_singlecell/TNBC_result/figure2"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

d <- readRDS(file.path(OUT, "region_assignment.rds"))
d$region <- factor(d$region, levels = names(region_colors))
cat("[Data] spots:", nrow(d), "; samples:", length(unique(d$sample)), "\n")

non_tumor <- c("T cells","Macrophages","Fibroblasts","Plasma cells","B cells",
               "Pericytes","Endothelial cells","Normal Epithelial cells")
all9 <- c("Tumor Epithelial cells", non_tumor)

## 坐标归一（每样本 x/y 到 0-1，保证 12 样本统一排版）
d <- d %>% group_by(sample) %>%
  mutate(xn = (x - min(x)) / (max(x) - min(x)),
         yn = (y - min(y)) / (max(y) - min(y))) %>% ungroup()

## ---------- A 肿瘤比例空间图（12 样本 grid） ----------
pA <- ggplot(d, aes(x = xn, y = yn, color = tumor_prop)) +
  geom_point(size = 0.8, stroke = 0) +
  facet_wrap(~ sample, ncol = 4) +
  scale_color_gradientn(colors = c("#D9D9D9", "#FDDBC7", "#C51B7D"),
                        name = "Tumor\nproportion", limits = c(0, 1)) +
  coord_fixed() +
  theme_void(base_size = 7) +
  theme(strip.text = element_text(size = 5, face = "bold"),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 5),
        plot.title = element_text(size = 8, face = "bold", hjust = 0)) +
  ggtitle("RCTD tumor proportion across 12 TNBC samples")
ggsave(file.path(OUT, "figure2A.pdf"), pA, width = 7.0, height = 5.4, device = cairo_pdf)
ggsave(file.path(OUT, "figure2A.png"), pA, width = 7.0, height = 5.4, dpi = 300)

## ---------- B 三区空间图（12 样本 grid） ----------
pB <- ggplot(d, aes(x = xn, y = yn, color = region)) +
  geom_point(size = 0.8, stroke = 0) +
  facet_wrap(~ sample, ncol = 4) +
  scale_color_manual(values = region_colors, name = "Region") +
  coord_fixed() +
  theme_void(base_size = 7) +
  theme(strip.text = element_text(size = 5, face = "bold"),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 5),
        plot.title = element_text(size = 8, face = "bold", hjust = 0)) +
  ggtitle("Tumor Core / Interface / Stroma zones")
ggsave(file.path(OUT, "figure2B.pdf"), pB, width = 7.0, height = 5.4, device = cairo_pdf)
ggsave(file.path(OUT, "figure2B.png"), pB, width = 7.0, height = 5.4, dpi = 300)

## ---------- C 肿瘤比例直方图 + 阈值 ----------
pC <- ggplot(d, aes(x = tumor_prop)) +
  geom_histogram(bins = 60, fill = "#C51B7D", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = c(0.2, 0.5), linetype = "dashed", color = "grey30", linewidth = 0.35) +
  annotate("text", x = 0.1, y = Inf, vjust = 1.6, hjust = 0, size = 2.2,
           label = "Stroma\n(tumor < 0.2)", color = "#377EB8", fontface = "bold") +
  annotate("text", x = 0.35, y = Inf, vjust = 1.6, hjust = 0.5, size = 2.2,
           label = "Interface\n(0.2–0.5)", color = "#4DAF4A", fontface = "bold") +
  annotate("text", x = 0.75, y = Inf, vjust = 1.6, hjust = 0.5, size = 2.2,
           label = "Core\n(> 0.5)", color = "#C51B7D", fontface = "bold") +
  labs(x = "Tumor proportion (RCTD)", y = "Spots",
       title = "Boundary definition by tumor-fraction gradient") +
  theme_bw(base_size = 8) + theme(panel.grid = element_blank())
ggsave(file.path(OUT, "figure2C.pdf"), pC, width = 3.4, height = 2.6, device = cairo_pdf)
ggsave(file.path(OUT, "figure2C.png"), pC, width = 3.4, height = 2.6, dpi = 300)

## ---------- D 每样本区域占比 ----------
regfrac <- d %>% count(sample, region) %>%
  group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup()
regfrac$sample <- factor(regfrac$sample,
                         levels = regfrac %>% filter(region == "Tumor Core") %>%
                           arrange(frac) %>% pull(sample))
pD <- ggplot(regfrac, aes(x = sample, y = frac, fill = region)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = region_colors, name = "Region") +
  labs(x = NULL, y = "Fraction of spots", title = "Zone composition per sample") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        panel.grid = element_blank(), legend.key.size = unit(2.5, "mm"))
ggsave(file.path(OUT, "figure2D.pdf"), pD, width = 4.6, height = 2.8, device = cairo_pdf)
ggsave(file.path(OUT, "figure2D.png"), pD, width = 4.6, height = 2.8, dpi = 300)
write.csv(regfrac, file.path(OUT, "region_fraction_by_sample.csv"), row.names = FALSE)

## ---------- E 去卷积质量 QC（9 类比例之和） ----------
d$prop_sum <- rowSums(d[, all9])
pE <- ggplot(d, aes(x = prop_sum)) +
  geom_histogram(bins = 60, fill = "#4DAF4A", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30", linewidth = 0.35) +
  labs(x = "Sum of 9 deconvolved fractions", y = "Spots",
       title = "RCTD deconvolution QC (per-spot sum ≈ 1)") +
  theme_bw(base_size = 8) + theme(panel.grid = element_blank())
ggsave(file.path(OUT, "figure2E.pdf"), pE, width = 3.4, height = 2.6, device = cairo_pdf)
ggsave(file.path(OUT, "figure2E.png"), pE, width = 3.4, height = 2.6, dpi = 300)

## ---------- F 距离肿瘤核心的边界距离梯度（2 代表样本） ----------
## 每 spot 到最近 Core spot 的欧氏距离 → 物理侵袭前沿梯度
dist_to_core <- d %>% group_by(sample) %>% group_modify(~ {
  core <- .x[.x$region == "Tumor Core", c("x","y")]
  if (nrow(core) < 3) { .x$dist_core <- NA_real_; return(.x) }
  dd <- as.matrix(dist(rbind(as.matrix(.x[, c("x","y")]), as.matrix(core))))
  nd <- dd[1:nrow(.x), (nrow(.x)+1):(nrow(.x)+nrow(core)), drop = FALSE]
  .x$dist_core <- apply(nd, 1, min)
  .x
}) %>% ungroup()
d <- merge(d, dist_to_core[, c("barcode","dist_core")], by = "barcode")

rep_samples <- c("GSM6433589_093C", "GSM6433593_094C")
pF <- ggplot(d %>% filter(sample %in% rep_samples), aes(x = xn, y = yn, color = dist_core)) +
  geom_point(size = 1.0, stroke = 0) +
  facet_wrap(~ sample, ncol = 2) +
  scale_color_gradientn(colors = c("#C51B7D", "#FDDBC7", "#377EB8"),
                        name = "Dist.\nto core") +
  coord_fixed() +
  theme_void(base_size = 7) +
  theme(strip.text = element_text(size = 6, face = "bold"),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 5),
        plot.title = element_text(size = 8, face = "bold", hjust = 0)) +
  ggtitle("Distance-to-tumor-core gradient (invasive front)")
ggsave(file.path(OUT, "figure2F.pdf"), pF, width = 5.0, height = 2.7, device = cairo_pdf)
ggsave(file.path(OUT, "figure2F.png"), pF, width = 5.0, height = 2.7, dpi = 300)

## ---------- 保存数据 ----------
write.csv(d %>% select(barcode, sample, region, tumor_prop, prop_sum, dist_core, x, y, all_of(all9)),
          file.path(OUT, "region_assignment_extended.csv"), row.names = FALSE)
saveRDS(d, file.path(OUT, "region_assignment_extended.rds"))
cat("[Done] Figure 2 完成 (6 panel)\n")
