# =============================================================================
# TNBC interface-niche code repository
# Original source: project_singlecell/TNBC_result/figure19/analysis.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

## ============================================================
## Figure 19 — BANKSY 空间域识别 vs RCTD 三区 + 非空间 kmeans 基线
##   核心问题:
##   1) BANKSY 空间域与 RCTD 三区(Core/Interface/Stroma)是否一致 (ARI)
##   2) BANKSY 域能否复现 纤维化区室(MMP11+CAF) 与 免疫生态位(SPP1+/C1Q+TAM) 空间分离
##   3) 空间感知(BANKSY) vs 非空间(kmeans) 谁更清晰区分两区室
##   (预注册: k_geom=18, lambda=0.2, npcs=20, leiden res=0.8; kmeans k=BANKSY 域数)
## ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Banksy)
  library(SingleCellExperiment)
  library(scater)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(mclust)
  library(Matrix)
  library(scales)
})

set.seed(42)
OUT   <- "/path/to/project/project_singlecell/TNBC_result/figure19"
OUT20 <- "/path/to/project/project_singlecell/TNBC_result/figure20"
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(OUT20, showWarnings = FALSE, recursive = TRUE)
DATA <- "/path/to/data/other_ans_project/三阴性乳腺癌分析/data"
FIG2 <- "/path/to/project/project_singlecell/TNBC_result/figure2"

K_GEOM <- 18L
LAMBDA <- 0.2
NPCS   <- 20L
RES    <- 0.8

region_colors <- c("Tumor Core" = "#C51B7D", "Interface" = "#4DAF4A", "Stroma" = "#377EB8")

## ---------- 1. SCT 表达矩阵 (genes x spots) ----------
cat("[1] loading SCT expression ...\n")
obj <- readRDS(file.path(DATA, "integrated_seurat_12.rds"))
X <- GetAssayData(obj, assay = "SCT", slot = "data")
rm(obj); gc()
cells <- colnames(X)
cat(sprintf("[1] SCT: %d genes x %d spots\n", nrow(X), ncol(X)))

## ---------- 2. 三区标签 + 坐标 ----------
cat("[2] loading region assignment ...\n")
ra <- read.csv(file.path(FIG2, "region_assignment.csv"),
               stringsAsFactors = FALSE, check.names = FALSE)
ra <- ra[match(cells, ra$barcode), ]
stopifnot(all(ra$barcode == cells))
sample_id <- ra$sample
region    <- ra$region
xy        <- ra[, c("x", "y")]

## ---------- 3. 35 亚群去卷积比例 ----------
cat("[3] loading 35-subcluster proportions ...\n")
sub <- readRDS(file.path(DATA, "Visium_spatial_deconvolution_亚群.rds"))
submd <- sub@meta.data
rm(sub); gc()
meta_cols <- c("orig.ident","nCount_Spatial","nFeature_Spatial","sample","percent.mt",
               "nCount_SCT","nFeature_SCT","SCT_snn_res.0.5","seurat_clusters")
subcols <- colnames(submd)[!(colnames(submd) %in% meta_cols)]
sub <- as.matrix(submd[, subcols, drop = FALSE])
if (all(cells %in% rownames(sub))) {
  sub <- sub[cells, , drop = FALSE]
} else {
  base <- function(b) sub("_\\d+$", "", b)
  sub <- sub[match(base(cells), base(rownames(sub))), , drop = FALSE]
}
stopifnot(nrow(sub) == length(cells))
rownames(sub) <- cells
cat(sprintf("[3] subcluster proportions: %d spots x %d subclusters\n", nrow(sub), ncol(sub)))

## 关键亚群 (与任务一致)
stopifnot(all(c("MMP11 CAFs","SPP1 TAMs","C1Q TAMs","CCL18 TAMs") %in% subcols))
CD8_subs <- grep("^CD8_", subcols, value = TRUE)
CD4_subs <- grep("^CD4_", subcols, value = TRUE)

## ---------- 保存 mistyR 输入 (base RDS, 可被 Seurat 4.4.0 / MISTy 读取) ----------
saveRDS(as.data.frame(sub), file.path(OUT20, "subproportions.rds"))
saveRDS(data.frame(barcode = cells, sample = sample_id, region = region,
                   x = xy$x, y = xy$y, stringsAsFactors = FALSE),
        file.path(OUT20, "spatial_meta.rds"))
cat("[3] saved mistyR inputs -> figure20/subproportions.rds, spatial_meta.rds\n")

## ---------- 4. 每样本 BANKSY 空间域 + 非空间 kmeans (带缓存) ----------
samples <- unique(sample_id)
dom_all <- rep(NA_character_, length(cells)); names(dom_all) <- cells
km_all  <- rep(NA_character_, length(cells)); names(km_all)  <- cells

cache_file <- file.path(OUT, "banksy_loop_outputs.rds")
if (file.exists(cache_file)) {
  cat("[4] loading cached BANKSY loop outputs ...\n")
  ll <- readRDS(cache_file)
  dom_all <- ll$dom_all; km_all <- ll$km_all
  dom_records <- ll$dom_records; ari_rows <- ll$ari_rows
} else {
  cat("[4] running BANKSY + kmeans per sample ...\n")
  dom_records <- list()
  ari_rows    <- list()
  for (s in samples) {
    idx    <- which(sample_id == s)
    cs     <- cells[idx]
    expr_s <- X[, idx, drop = FALSE]
    cd     <- data.frame(x = xy$x[idx], y = xy$y[idx], row.names = cs)
    sce    <- SingleCellExperiment(assays = list(SCT = expr_s), colData = cd)

    sce <- computeBanksy(sce, assay_name = "SCT", coord_names = c("x", "y"),
                         k_geom = K_GEOM, compute_agf = TRUE)
    sce <- runBanksyPCA(sce, use_agf = TRUE, lambda = LAMBDA, npcs = NPCS)
    sce <- runBanksyPCA(sce, use_agf = FALSE, npcs = NPCS)

    ## 空间域 (leiden)
    cn0 <- colnames(colData(sce))
    sce <- clusterBanksy(sce, use_agf = TRUE, lambda = LAMBDA, resolution = RES,
                         algo = "leiden", npcs = NPCS, k_neighbors = 50L)
    dom_col <- setdiff(colnames(colData(sce)), cn0)
    stopifnot(length(dom_col) == 1)
    dom     <- as.character(colData(sce)[[dom_col]])
    ndom    <- length(unique(dom))

    ## 非空间 kmeans (k = BANKSY 域数)
    cn1 <- colnames(colData(sce))
    sce <- clusterBanksy(sce, use_agf = FALSE, algo = "kmeans",
                         kmeans.centers = ndom, npcs = NPCS)
    km_col <- setdiff(colnames(colData(sce)), cn1)
    stopifnot(length(km_col) == 1)
    km     <- as.character(colData(sce)[[km_col]])

    dom_all[idx] <- paste0(s, "_D", dom)
    km_all[idx]  <- paste0(s, "_K", km)

    ari_d <- mclust::adjustedRandIndex(dom, region[idx])
    ari_k <- mclust::adjustedRandIndex(km,  region[idx])
    ari_rows[[length(ari_rows) + 1]] <- data.frame(
      sample = s, n_spots = length(idx), n_domains = ndom,
      ARI_banksy_vs_region = ari_d, ARI_kmeans_vs_region = ari_k,
      stringsAsFactors = FALSE)

    ## 域组成
    for (d in unique(dom)) {
      didx <- idx[dom == d]
      dom_records[[length(dom_records) + 1]] <- data.frame(
        sample = s, domain = d, n_spots = length(didx),
        MMP11_CAFs = mean(sub[didx, "MMP11 CAFs"]),
        SPP1_TAMs  = mean(sub[didx, "SPP1 TAMs"]),
        C1Q_TAMs   = mean(sub[didx, "C1Q TAMs"]),
        CCL18_TAMs = mean(sub[didx, "CCL18 TAMs"]),
        CD8_T      = mean(rowSums(sub[didx, CD8_subs, drop = FALSE])),
        CD4_T      = mean(rowSums(sub[didx, CD4_subs, drop = FALSE])),
        frac_Core      = mean(region[didx] == "Tumor Core"),
        frac_Interface = mean(region[didx] == "Interface"),
        frac_Stroma    = mean(region[didx] == "Stroma"),
        stringsAsFactors = FALSE)
    }
    cat(sprintf("  [%s] n=%d domains=%d  ARI(BANKSY)=%.3f ARI(kmeans)=%.3f\n",
                s, length(idx), ndom, ari_d, ari_k))
    rm(sce, expr_s); gc()
  }
  saveRDS(list(dom_all = dom_all, km_all = km_all,
               dom_records = dom_records, ari_rows = ari_rows), cache_file)
  cat("[4] done. cached -> banksy_loop_outputs.rds\n")
}
ari_tab <- do.call(rbind, ari_rows)
cat("[4] total domains:", length(unique(dom_all)), "\n")

## ---------- 5. ARI 一致性统计 (每样本 + 合并) ----------
ari_pooled_d <- mclust::adjustedRandIndex(dom_all, region)
ari_pooled_k <- mclust::adjustedRandIndex(km_all,  region)
ari_tab <- rbind(ari_tab, data.frame(
  sample = "POOLED", n_spots = length(cells),
  n_domains = length(unique(dom_all)),
  ARI_banksy_vs_region = ari_pooled_d,
  ARI_kmeans_vs_region = ari_pooled_k, stringsAsFactors = FALSE))
write.csv(ari_tab, file.path(OUT, "ari_consistency.csv"), row.names = FALSE)
cat(sprintf("[5] ARI pooled: BANKSY=%.3f kmeans=%.3f\n", ari_pooled_d, ari_pooled_k))

## ---------- 6. 域组成 + CAF/TAM 区室分离 ----------
dr <- do.call(rbind, dom_records)
## 样本内 z 化 (跨域)
dr <- dr %>% group_by(sample) %>%
  mutate(z_CAF = as.numeric(scale(MMP11_CAFs)),
         z_TAM = as.numeric(scale(SPP1_TAMs + C1Q_TAMs))) %>% ungroup()
dr$CAF_dom   <- dr$z_CAF > 1
dr$TAM_dom   <- dr$z_TAM > 1
dr$mixed     <- dr$CAF_dom & dr$TAM_dom
dr$dominant  <- ifelse(dr$mixed, "mixed",
                ifelse(dr$CAF_dom, "CAF",
                ifelse(dr$TAM_dom, "TAM", "neutral")))
write.csv(dr, file.path(OUT, "domain_enrichment.csv"), row.names = FALSE)

## 区室分离统计 (BANKSY)
cor_per_sample <- sapply(split(dr, dr$sample), function(dd) {
  if (nrow(dd) < 3 || sd(dd$MMP11_CAFs) == 0) return(NA_real_)
  cor(dd$MMP11_CAFs, dd$SPP1_TAMs + dd$C1Q_TAMs)
})
sep_d <- data.frame(
  method = "BANKSY",
  n_domains = nrow(dr),
  n_CAF_dom = sum(dr$CAF_dom),
  n_TAM_dom = sum(dr$TAM_dom),
  n_mixed   = sum(dr$mixed),
  mixed_fraction = mean(dr$mixed),
  CAF_TAM_cor_mean = mean(cor_per_sample, na.rm = TRUE),
  stringsAsFactors = FALSE)

## 空间邻接: CAF-dominant 域 spot 的邻域中是否含 TAM-dominant 域 spot
adj_fracs <- c()
for (s in samples) {
  didx <- which(sample_id == s)
  if (sum(dr$CAF_dom[dr$sample == s]) == 0 || sum(dr$TAM_dom[dr$sample == s]) == 0) next
  sx <- xy$x[didx]; sy <- xy$y[didx]
  dmat <- as.matrix(dist(cbind(sx, sy)))
  pitch <- median(apply(dmat, 1, function(z) min(z[z > 0])))
  thr <- pitch * 1.5
  dd <- dr[dr$sample == s, ]
  ddmap <- setNames(dd$dominant, paste0(s, "_D", dd$domain))
  spot_flags <- ddmap[dom_all[didx]]
  caf_idx <- which(spot_flags == "CAF")
  tam_idx <- which(spot_flags == "TAM")
  if (length(caf_idx) == 0 || length(tam_idx) == 0) next
  adj <- (dmat <= thr)
  has_tam_neigh <- rowSums(adj[caf_idx, tam_idx, drop = FALSE]) > 0
  adj_fracs <- c(adj_fracs, mean(has_tam_neigh))
}
sep_d$CAF_TAM_spatial_adjacency <- if (length(adj_fracs) > 0) mean(adj_fracs) else NA_real_
write.csv(sep_d, file.path(OUT, "compartment_separation.csv"), row.names = FALSE)
cat("[6] BANKSY domains: CAF-dom=", sep_d$n_CAF_dom, " TAM-dom=", sep_d$n_TAM_dom,
    " mixed=", sep_d$n_mixed, " cor=", round(sep_d$CAF_TAM_cor_mean, 3),
    " adj=", round(sep_d$CAF_TAM_spatial_adjacency, 3), "\n", sep = "")

## ---------- 7. 图: 空间域 / kmeans / 三区 空间图 (12 样本) ----------
plot_sp <- function(lab, title, pal = NULL) {
  d <- data.frame(sample = sample_id, x = xy$x, y = xy$y, lab = lab)
  d <- d %>% group_by(sample) %>%
    mutate(xn = (x - min(x)) / (max(x) - min(x)),
           yn = (y - min(y)) / (max(y) - min(y))) %>% ungroup()
  p <- ggplot(d, aes(x = xn, y = yn, color = lab)) +
    geom_point(size = 0.4, stroke = 0) +
    facet_wrap(~ sample, ncol = 4) +
    coord_fixed() + theme_void(base_size = 7) +
    theme(strip.text = element_text(size = 5, face = "bold"),
          legend.key.size = unit(2.5, "mm"),
          legend.text = element_text(size = 5),
          plot.title = element_text(size = 8, face = "bold", hjust = 0)) +
    ggtitle(title)
  if (!is.null(pal)) p <- p + scale_color_manual(values = pal, name = "")
  p
}
p_dom <- plot_sp(dom_all, "BANKSY spatial domains (leiden, k=18, lam=0.2, res=0.8)")
p_km  <- plot_sp(km_all,  "Non-spatial kmeans (k = n BANKSY domains)")
p_reg <- plot_sp(region,  "RCTD three-region assignment") +
  scale_color_manual(values = region_colors, name = "")
ggsave(file.path(OUT, "figure19A_domains.png"), p_dom, width = 7.2, height = 5.6, dpi = 300)
ggsave(file.path(OUT, "figure19B_kmeans.png"),  p_km,  width = 7.2, height = 5.6, dpi = 300)
ggsave(file.path(OUT, "figure19C_regions.png"), p_reg, width = 7.2, height = 5.6, dpi = 300)

## ---------- 8. 域组成热图 ----------
hm_feat <- c("MMP11_CAFs","SPP1_TAMs","C1Q_TAMs","CCL18_TAMs","CD8_T","CD4_T",
             "frac_Core","frac_Interface","frac_Stroma")
hml <- dr %>% select(sample, domain, all_of(hm_feat)) %>%
  group_by(sample) %>%
  mutate(across(all_of(hm_feat), ~ as.numeric(scale(.x)))) %>% ungroup()
hml$row <- paste0(hml$sample, "_D", hml$domain)
hml$dominant <- dr$dominant
hml$row <- factor(hml$row, levels = rev(unique(hml$row)))
hml_long <- pivot_longer(hml, cols = all_of(hm_feat), names_to = "feature", values_to = "z")
hml_long$feature <- factor(hml_long$feature, levels = hm_feat)
p_hm <- ggplot(hml_long, aes(x = feature, y = row, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                       name = "z-score", na.value = "grey90") +
  labs(x = NULL, y = NULL, title = "Domain composition (per-sample z-scored)") +
  theme_bw(base_size = 7) +
  theme(axis.text.y = element_text(size = 3),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        panel.grid = element_blank(),
        plot.title = element_text(size = 8, face = "bold", hjust = 0))
ggsave(file.path(OUT, "figure19D_composition_heatmap.png"), p_hm,
       width = 5.2, height = 9.0, dpi = 300)

## ---------- 9. ARI 对比 (BANKSY vs kmeans) ----------
ari_plot <- ari_tab %>% filter(sample != "POOLED") %>%
  select(sample, ARI_banksy_vs_region, ARI_kmeans_vs_region) %>%
  pivot_longer(-sample, names_to = "method", values_to = "ARI") %>%
  mutate(method = ifelse(grepl("banksy", method), "BANKSY (spatial)", "kmeans (non-spatial)"))
p_ari <- ggplot(ari_plot, aes(x = sample, y = ARI, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = c("BANKSY (spatial)" = "#C51B7D",
                               "kmeans (non-spatial)" = "#999999"), name = NULL) +
  geom_hline(yintercept = ari_pooled_d, linetype = "dashed", color = "#C51B7D", linewidth = 0.4) +
  geom_hline(yintercept = ari_pooled_k, linetype = "dashed", color = "#666666", linewidth = 0.4) +
  labs(x = NULL, y = "Adjusted Rand Index vs RCTD regions",
       title = "Domain-region consistency: spatial vs non-spatial") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        panel.grid = element_blank(), legend.key.size = unit(2.5, "mm"))
ggsave(file.path(OUT, "figure19E_ari_comparison.png"), p_ari, width = 6.0, height = 3.2, dpi = 300)

## ---------- 10. 结论 ----------
cat("\n[结论]\n")
cat(sprintf("  ARI(BANKSY vs region) = %.3f ; ARI(kmeans vs region) = %.3f (pooled)\n",
            ari_pooled_d, ari_pooled_k))
cat(sprintf("  BANKSY 域: %d 个 (CAF-dominant=%d, TAM-dominant=%d, mixed=%d, mixed_frac=%.2f)\n",
            sep_d$n_domains, sep_d$n_CAF_dom, sep_d$n_TAM_dom, sep_d$n_mixed, sep_d$mixed_fraction))
cat(sprintf("  CAF-TAM 域间相关(均值)=%.3f ; CAF/TAM 域空间邻接=%.3f\n",
            sep_d$CAF_TAM_cor_mean, sep_d$CAF_TAM_spatial_adjacency))
cat("[Done] Figure 19 BANKSY 分析完成\n")
