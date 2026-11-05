# =============================================================================
# TNBC interface-niche code repository
# Original source: project_spatial/TNBC/patient_ALL/TAM_clone/notebooks/step01_tam_clone_coloc.ipynb (converted to R)
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

# 环境与参数设置：路径、随机种子、置换参数、可估门槛、保存辅助函数
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(qs); library(data.table); library(ggplot2); library(grid)})
set.seed(42)

ROOT       <- "/path/to/project/project_spatial/TNBC/patient_ALL"
CLONE_DIR  <- file.path(ROOT, "clone")
REGION_FILE <- "/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment_extended.csv"
TAM_FILE   <- "/path/to/project/project_spatial/TNBC/TNBC_analysis/SpatialInferCNV/results/step3_spatial_tme_groups/data/spatial_tme_groups.qs"
OUT <- file.path(ROOT, "TAM_clone")
TAB <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")
REV <- file.path(ROOT, "review", "TAM_clone")          # 内部证据，不属于交付物
for (d in c(TAB, FIG, DAT, REV)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(OUT, "notebooks", ".keep"))

PATIENTS <- c("P9", "P12", "P13", "P14", "P15")
SLICE2P  <- c("092A"="P9","092B"="P9","093C"="P12","093D"="P12",
              "094A"="P13","094B"="P13","094C"="P14","094D"="P14",
              "095A"="P15","095B"="P15")
SLICES <- names(SLICE2P)
ZONES  <- c("Tumor Core", "Interface", "Stroma")
TAMS   <- c("C1QC_TAM", "SPP1_TAM")
TAM_LAB <- c(C1QC_TAM = "C1Q TAMs", SPP1_TAM = "SPP1 TAMs")

N_PERM <- 1000; SEED <- 42
MIN_REGION_N <- 100; MIN_CLONE_N <- 20; MIN_POS_N <- 20   # 分层可估门槛（step4-04 原门槛）
MIN_POOL_B   <- 30;  LOW_N_B <- 10                        # step5 原生池/低样本门槛
ENRICH_Q <- 0.75; K_MAIN <- 6; K_SENS <- 18               # step5 参数

PLOT_THEME <- "nature"
NPG <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000")
DIV_LOW <- "#2166AC"; DIV_MID <- "white"; DIV_HIGH <- "#B2182B"

save_obj <- function(obj, filepath) {
  qs_path <- paste0(sub("\\.(rds|qs)$", "", filepath), ".qs")
  if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(obj, qs_path, preset = "fast")
    cat(sprintf("[Save] qs: %s (%.1f MB)\n", qs_path, file.info(qs_path)$size / 1e6))
  } else stop("qs unavailable")
}
save_fig <- function(p, stem, w, h) {
  ggsave(file.path(FIG, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 600, bg = "white", device = "png")
  ggsave(file.path(FIG, paste0(stem, ".pdf")), p, width = w, height = h,
         bg = "white", device = grDevices::pdf)
  invisible(stem)
}

cat("[Param] seed=", SEED, "| perm=", N_PERM, "| radius=2*median1NN | k=", K_MAIN, "/", K_SENS,
    "| methodA positive=max(1%,2*mean) | methodB enrich q=", ENRICH_Q, "\n", sep = "")
cat("[Threshold] stratified: region_n>=", MIN_REGION_N, " clone_n>=", MIN_CLONE_N,
    " TAMpos_n>=", MIN_POS_N, " | methodB pool>=", MIN_POOL_B, " low_n<", LOW_N_B, "\n", sep = "")
cat("[Out] ", OUT, "\n", sep = "")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 环境与参数设置：路径、随机种子、置换参数、可估门槛、保存辅助函数（分层仅 Core/Interface）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(qs); library(data.table); library(ggplot2); library(grid)})
set.seed(42)

ROOT       <- "/path/to/project/project_spatial/TNBC/patient_ALL"
CLONE_DIR  <- file.path(ROOT, "clone")
REGION_FILE <- "/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment_extended.csv"
TAM_FILE   <- "/path/to/project/project_spatial/TNBC/TNBC_analysis/SpatialInferCNV/results/step3_spatial_tme_groups/data/spatial_tme_groups.qs"
OUT <- file.path(ROOT, "TAM_clone")
TAB <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")
REV <- file.path(ROOT, "review", "TAM_clone")
for (d in c(TAB, FIG, DAT, REV)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(OUT, "notebooks", ".keep"))

PATIENTS <- c("P9", "P12", "P13", "P14", "P15")
SLICE2P  <- c("092A"="P9","092B"="P9","093C"="P12","093D"="P12",
              "094A"="P13","094B"="P13","094C"="P14","094D"="P14",
              "095A"="P15","095B"="P15")
SLICES <- names(SLICE2P)
ZONES  <- c("Tumor Core", "Interface")          # Stroma 已按用户要求移除（clone spot <=1, 不可估）
TAMS   <- c("C1QC_TAM", "SPP1_TAM")
TAM_LAB <- c(C1QC_TAM = "C1Q TAMs", SPP1_TAM = "SPP1 TAMs")

N_PERM <- 1000; SEED <- 42
MIN_REGION_N <- 100; MIN_CLONE_N <- 20; MIN_POS_N <- 20
MIN_POOL_B   <- 30;  LOW_N_B <- 10
ENRICH_Q <- 0.75; K_MAIN <- 6; K_SENS <- 18

PLOT_THEME <- "nature"
NPG <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000")
DIV_LOW <- "#2166AC"; DIV_MID <- "white"; DIV_HIGH <- "#B2182B"

save_obj <- function(obj, filepath) {
  qs_path <- paste0(sub("\\.(rds|qs)$", "", filepath), ".qs")
  if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(obj, qs_path, preset = "fast")
    cat(sprintf("[Save] qs: %s (%.1f MB)\n", qs_path, file.info(qs_path)$size / 1e6))
  } else stop("qs unavailable")
}
save_fig <- function(p, stem, w, h) {
  ggsave(file.path(FIG, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 600, bg = "white", device = "png")
  ggsave(file.path(FIG, paste0(stem, ".pdf")), p, width = w, height = h,
         bg = "white", device = grDevices::pdf)
  invisible(stem)
}

cat("[Param] seed=", SEED, "| perm=", N_PERM, "| radius=2*median1NN | k=", K_MAIN, "/", K_SENS,
    "| methodA positive=max(1%,2*mean) | methodB enrich q=", ENRICH_Q, "\n", sep = "")
cat("[Zones] ", paste(ZONES, collapse = " / "), " (Stroma removed)\n", sep = "")
cat("[Threshold] stratified: region_n>=", MIN_REGION_N, " clone_n>=", MIN_CLONE_N,
    " TAMpos_n>=", MIN_POS_N, " | methodB pool>=", MIN_POOL_B, " low_n<", LOW_N_B, "\n", sep = "")
cat("[Out] ", OUT, "\n", sep = "")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 构建关联键并合并三区/TAM，逐片核对关联率，输出输入汇总表
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

cm[, key := as.character(cell)]
cm[, slice := as.character(slice)]
cm[, clone_id := paste0(patient, "_", clone)]

ra[, slice := sub("^GSM[0-9]+_", "", as.character(sample))]
ra[, key := paste0(slice, "_", sub("_[0-9]+$", "", as.character(barcode)))]
tam[, slice := as.character(slice)]
tam[, key := paste0(slice, "_", sub("_[0-9]+$", "", as.character(barcode)))]
ra  <- ra[slice %in% SLICES]; tam <- tam[slice %in% SLICES]
stopifnot(!any(duplicated(ra$key)), !any(duplicated(tam$key)))

ra_slim  <- ra[, .(key, slice, region, tumor_prop, dist_core, x_reg = x, y_reg = y)]
tam_slim <- tam[, .(key, C1QC_TAM, SPP1_TAM)]
d <- merge(cm, ra_slim[, !"slice"], by = "key", all.x = TRUE)
d <- merge(d, tam_slim, by = "key", all.x = TRUE)
stopifnot(nrow(d) == nrow(cm), "slice" %in% names(d))

join_chk <- d[, .(patient = patient[1], n_clone_spots = .N,
                  n_region_join = sum(!is.na(region)),
                  n_tam_join    = sum(!is.na(C1QC_TAM) & !is.na(SPP1_TAM))), by = slice]
join_chk[, join_rate_pct := round(pmin(n_region_join, n_tam_join) / n_clone_spots * 100, 2)]
cat("[join rate per slice]\n"); print(join_chk[order(slice)])
stopifnot(all(join_chk$join_rate_pct == 100))

univ <- merge(ra_slim, tam_slim, by = "key", all.x = TRUE)
cat("\n[universe] all tissue spots in 10 slices:", nrow(univ),
    "| with region:", sum(!is.na(univ$region)), "| with TAM:", sum(!is.na(univ$C1QC_TAM)), "\n")
stopifnot(all(!is.na(univ$region)), all(!is.na(univ$C1QC_TAM)))

zone_tab <- d[, .(n_TumorCore = sum(region == "Tumor Core"),
                  n_Interface = sum(region == "Interface"),
                  n_Stroma    = sum(region == "Stroma")), by = .(patient, slice)]
univ_zone <- univ[, .(n_universe = .N,
                      univ_TumorCore = sum(region == "Tumor Core"),
                      univ_Interface = sum(region == "Interface"),
                      univ_Stroma    = sum(region == "Stroma"),
                      univ_C1QC_bg = round(mean(C1QC_TAM), 6),
                      univ_SPP1_bg = round(mean(SPP1_TAM), 6)), by = slice]
clones_n <- d[, .(n_clones = uniqueN(clone_id)), by = slice]

input_summary <- merge(join_chk[, .(slice, patient, n_clone_spots, join_rate_pct, n_region_join, n_tam_join)],
                       clones_n, by = "slice", all.x = TRUE)
input_summary <- merge(input_summary, zone_tab, by = c("patient","slice"), all.x = TRUE)
input_summary <- merge(input_summary, univ_zone, by = "slice", all.x = TRUE)
setcolorder(input_summary, c("slice","patient","n_clones","n_clone_spots",
                             "n_TumorCore","n_Interface","n_Stroma",
                             "n_universe","univ_TumorCore","univ_Interface","univ_Stroma",
                             "join_rate_pct","n_region_join","n_tam_join","univ_C1QC_bg","univ_SPP1_bg"))
setorder(input_summary, slice)
fwrite(input_summary, file.path(TAB, "coloc_input_summary.csv"))
cat("\n[input summary]\n"); print(input_summary)
cat("\n[clone x slice]\n"); print(d[, .N, by = .(patient, slice, clone_id)][order(patient, slice, clone_id)])
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 构建关联键并合并三区/TAM，逐片核对关联率，输出输入汇总表（仅 Core/Interface）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

cm[, key := as.character(cell)]
cm[, slice := as.character(slice)]
cm[, clone_id := paste0(patient, "_", clone)]
ra[, slice := sub("^GSM[0-9]+_", "", as.character(sample))]
ra[, key := paste0(slice, "_", sub("_[0-9]+$", "", as.character(barcode)))]
tam[, slice := as.character(slice)]
tam[, key := paste0(slice, "_", sub("_[0-9]+$", "", as.character(barcode)))]
ra  <- ra[slice %in% SLICES]; tam <- tam[slice %in% SLICES]
stopifnot(!any(duplicated(ra$key)), !any(duplicated(tam$key)))

ra_slim  <- ra[, .(key, slice, region, tumor_prop, dist_core, x_reg = x, y_reg = y)]
tam_slim <- tam[, .(key, C1QC_TAM, SPP1_TAM)]
d <- merge(cm, ra_slim[, !"slice"], by = "key", all.x = TRUE)
d <- merge(d, tam_slim, by = "key", all.x = TRUE)
stopifnot(nrow(d) == nrow(cm), "slice" %in% names(d))

join_chk <- d[, .(patient = patient[1], n_clone_spots = .N,
                  n_region_join = sum(!is.na(region)),
                  n_tam_join    = sum(!is.na(C1QC_TAM) & !is.na(SPP1_TAM))), by = slice]
join_chk[, join_rate_pct := round(pmin(n_region_join, n_tam_join) / n_clone_spots * 100, 2)]
cat("[join rate per slice]\n"); print(join_chk[order(slice)])
stopifnot(all(join_chk$join_rate_pct == 100))

univ <- merge(ra_slim, tam_slim, by = "key", all.x = TRUE)
cat("\n[universe] all tissue spots in 10 slices:", nrow(univ),
    "| with region:", sum(!is.na(univ$region)), "| with TAM:", sum(!is.na(univ$C1QC_TAM)), "\n")
stopifnot(all(!is.na(univ$region)), all(!is.na(univ$C1QC_TAM)))

zone_tab <- d[, .(n_TumorCore = sum(region == "Tumor Core"),
                  n_Interface = sum(region == "Interface")), by = .(patient, slice)]
univ_zone <- univ[, .(n_universe = .N,
                      univ_TumorCore = sum(region == "Tumor Core"),
                      univ_Interface = sum(region == "Interface"),
                      univ_C1QC_bg = round(mean(C1QC_TAM), 6),
                      univ_SPP1_bg = round(mean(SPP1_TAM), 6)), by = slice]
clones_n <- d[, .(n_clones = uniqueN(clone_id)), by = slice]
input_summary <- merge(join_chk[, .(slice, patient, n_clone_spots, join_rate_pct, n_region_join, n_tam_join)],
                       clones_n, by = "slice", all.x = TRUE)
input_summary <- merge(input_summary, zone_tab, by = c("patient","slice"), all.x = TRUE)
input_summary <- merge(input_summary, univ_zone, by = "slice", all.x = TRUE)
setcolorder(input_summary, c("slice","patient","n_clones","n_clone_spots",
                             "n_TumorCore","n_Interface","n_universe","univ_TumorCore","univ_Interface",
                             "join_rate_pct","n_region_join","n_tam_join","univ_C1QC_bg","univ_SPP1_bg"))
setorder(input_summary, slice)
fwrite(input_summary, file.path(TAB, "coloc_input_summary.csv"))
cat("\n[input summary]\n"); print(input_summary)
cat("\n[clone x slice]\n"); print(d[, .N, by = .(patient, slice, clone_id)][order(patient, slice, clone_id)])
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 方法A（step4）：逐切片半径邻域对计数 + 1000 次 TAM 标签置换 + k=6 邻域富集；不可估格不保留点估计
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
set.seed(SEED)
for (g in TAMS) {
  col <- paste0("pos_", g); set(d, j = col, value = FALSE)
  for (s in SLICES) {
    i <- which(d$slice == s); x <- d[[g]][i]
    set(d, i = i, j = col, value = as.numeric(x >= max(0.01, 2 * mean(x))))
  }
}
pos_cnt <- d[, .(n_spots = .N, n_C1QC_pos = sum(pos_C1QC_TAM), n_SPP1_pos = sum(pos_SPP1_TAM)),
             by = .(patient, slice)]
cat("[positive spots, per slice]\n"); print(pos_cnt[order(slice)])
make_graph <- function(xy, r_mult = 2, kmax = 60) {
  n <- nrow(xy); D <- as.matrix(dist(xy)); diag(D) <- Inf
  nn <- median(apply(D, 1, min)); r <- r_mult * nn
  kk <- min(kmax, n - 1)
  ord <- t(apply(D, 1, function(z) order(z, method = "radix")[seq_len(kk)]))
  ei <- rep(seq_len(n), kk); ej <- as.vector(ord); ed <- D[cbind(ei, ej)]
  keep <- ed <= r & ei < ej
  list(ei = ei[keep], ej = ej[keep], ord = ord, r = r, nn = nn, n_edges = sum(keep))
}
res_rad <- list(); res_knn <- list()
for (s in SLICES) {
  ds <- d[slice == s, ]; n <- nrow(ds)
  G  <- make_graph(as.matrix(ds[, .(x, y)]))
  ei <- G$ei; ej <- G$ej
  kn6 <- G$ord[, seq_len(min(K_MAIN, n - 1)), drop = FALSE]
  clones <- sort(unique(ds$clone_id))
  for (cl in clones) {
    xc <- as.numeric(ds$clone_id == cl); idx <- which(xc == 1)
    w <- numeric(n)
    te <- ej[xc[ei] == 1]; w[te] <- w[te] + 1
    tf <- ei[xc[ej] == 1]; w[tf] <- w[tf] + 1
    nb_idx <- as.vector(kn6[idx, , drop = FALSE])
    for (g in TAMS) {
      yg <- as.numeric(ds[[paste0("pos_", g)]])
      obs_e <- sum(w * yg); null_e <- numeric(N_PERM)
      for (p in seq_len(N_PERM)) null_e[p] <- sum(w * sample(yg))
      z_e <- if (sd(null_e) > 0) (obs_e - mean(null_e)) / sd(null_e) else NA_real_
      obs_k <- mean(yg[nb_idx]); null_k <- numeric(N_PERM)
      for (p in seq_len(N_PERM)) null_k[p] <- mean(sample(yg)[nb_idx])
      z_k <- if (sd(null_k) > 0) (obs_k - mean(null_k)) / sd(null_k) else NA_real_
      base_row <- data.table(patient = ds$patient[1], slice = s, clone_id = cl,
                             clone = ds$clone[match(cl, ds$clone_id)],
                             n_clone = sum(xc), tam = TAM_LAB[[g]],
                             n_pos = sum(yg), pos_pct = round(mean(yg) * 100, 2),
                             low_n = (sum(xc) < MIN_CLONE_N | sum(yg) < MIN_POS_N))
      res_rad[[length(res_rad) + 1]] <- cbind(base_row, data.table(
        obs_edges = obs_e, null_mean_edges = mean(null_e), null_sd_edges = sd(null_e),
        enrichment = obs_e / max(mean(null_e), 1e-9), z = z_e,
        p_perm = (1 + sum(null_e >= obs_e)) / (1 + N_PERM),
        p_deplete = (1 + sum(null_e <= obs_e)) / (1 + N_PERM),
        median_1nn = G$nn, niche_radius = G$r, n_edges = G$n_edges,
        estimable = !(sum(xc) < MIN_CLONE_N | sum(yg) < MIN_POS_N)))
      res_knn[[length(res_knn) + 1]] <- cbind(base_row, data.table(
        k = K_MAIN, obs_neighbor_frac = obs_k, null_mean = mean(null_k), null_sd = sd(null_k),
        enrichment = obs_k / max(mean(null_k), 1e-9), z = z_k,
        p_perm = (1 + sum(null_k >= obs_k)) / (1 + N_PERM),
        p_deplete = (1 + sum(null_k <= obs_k)) / (1 + N_PERM),
        estimable = !(sum(xc) < MIN_CLONE_N | sum(yg) < MIN_POS_N)))
    }
  }
  cat(sprintf("[A] %s n=%d clones=%d edges=%d r=%.2f (1NN=%.2f)\n", s, n, length(clones), G$n_edges, G$r, G$nn))
}
mA_radius <- rbindlist(res_rad); mA_knn6 <- rbindlist(res_knn)
mA_radius[estimable == FALSE, c("obs_edges","null_mean_edges","null_sd_edges",
                                "enrichment","z","p_perm","p_deplete") := NA_real_]
mA_knn6[estimable == FALSE, c("obs_neighbor_frac","null_mean","null_sd",
                              "enrichment","z","p_perm","p_deplete") := NA_real_]
setorder(mA_radius, slice, clone_id, tam); setorder(mA_knn6, slice, clone_id, tam)
fwrite(mA_radius, file.path(TAB, "clone_tam_radius_perm_by_slice.csv"))
fwrite(mA_knn6,   file.path(TAB, "clone_tam_knn6_by_slice.csv"))
save_obj(list(radius = mA_radius, knn6 = mA_knn6, pos = pos_cnt), file.path(DAT, "methodA_by_slice"))
cat("\n[methodA] rows:", nrow(mA_radius), "| non-estimable:", sum(!mA_radius$estimable),
    "| non-est with value:", mA_radius[estimable == FALSE & !is.na(enrichment), .N], "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 方法A 分层（step4-04）：仅在 Tumor Core / Interface 内重建半径图并置换；Stroma 已移除
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
set.seed(SEED)
grid <- unique(d[, .(slice, patient, clone_id, clone)])
grid <- grid[, .(region = rep(ZONES, each = .N),
                 clone_id = rep(clone_id, times = length(ZONES)),
                 clone = rep(clone, times = length(ZONES))), by = .(slice, patient)]
grid <- grid[, .(tam = rep(unname(TAM_LAB[TAMS]), each = .N),
                 clone_id = rep(clone_id, times = length(TAMS)),
                 clone = rep(clone, times = length(TAMS))), by = .(slice, patient, region)]
cat("[grid] rows:", nrow(grid), "| zones:", paste(ZONES, collapse = " / "), "\n")
res <- list()
for (s in SLICES) {
  for (rg in ZONES) {
    ds <- d[slice == s & region == rg, ]; n <- nrow(ds)
    gsub <- grid[slice == s & region == rg]
    for (k in seq_len(nrow(gsub))) {
      cl <- gsub$clone_id[k]; g <- TAMS[which(unname(TAM_LAB) == gsub$tam[k])]
      stopifnot(length(g) == 1)
      n_cl <- sum(ds$clone_id == cl)
      yg <- if (n > 0) as.numeric(ds[[paste0("pos_", g)]]) else numeric(0)
      n_pos <- sum(yg)
      if (n < MIN_REGION_N) {
        res[[length(res) + 1]] <- cbind(gsub[k], data.table(
          n_region = n, n_clone = n_cl, n_pos = n_pos,
          obs_edges = NA_real_, null_mean_edges = NA_real_, enrichment = NA_real_,
          z = NA_real_, p_perm = NA_real_, p_deplete = NA_real_,
          niche_radius = NA_real_, n_edges = NA_integer_,
          estimable = FALSE, reason = sprintf("region_n=%d<100", n)))
        next
      }
      if (cl == gsub$clone_id[1]) { G <- make_graph(as.matrix(ds[, .(x, y)])); ei <- G$ei; ej <- G$ej }
      xc <- as.numeric(ds$clone_id == cl)
      reason <- NA_character_
      if (n_cl == 0) reason <- "clone absent in region (n=0)"
      else if (n_cl < MIN_CLONE_N) reason <- sprintf("clone_n=%d<20", n_cl)
      else if (n_pos < MIN_POS_N) reason <- sprintf("TAMpos_n=%d<20", n_pos)
      else if (G$n_edges == 0) reason <- "no edges within radius"
      if (!is.na(reason)) {
        res[[length(res) + 1]] <- cbind(gsub[k], data.table(
          n_region = n, n_clone = n_cl, n_pos = n_pos,
          obs_edges = NA_real_, null_mean_edges = NA_real_, enrichment = NA_real_,
          z = NA_real_, p_perm = NA_real_, p_deplete = NA_real_,
          niche_radius = G$r, n_edges = G$n_edges, estimable = FALSE, reason = reason))
        next
      }
      w <- numeric(n)
      te <- ej[xc[ei] == 1]; w[te] <- w[te] + 1
      tf <- ei[xc[ej] == 1]; w[tf] <- w[tf] + 1
      obs <- sum(w * yg); null <- numeric(N_PERM)
      for (p in seq_len(N_PERM)) null[p] <- sum(w * sample(yg))
      zz <- if (sd(null) > 0) (obs - mean(null)) / sd(null) else NA_real_
      res[[length(res) + 1]] <- cbind(gsub[k], data.table(
        n_region = n, n_clone = n_cl, n_pos = n_pos,
        obs_edges = obs, null_mean_edges = mean(null),
        enrichment = obs / max(mean(null), 1e-9), z = zz,
        p_perm = (1 + sum(null >= obs)) / (1 + N_PERM),
        p_deplete = (1 + sum(null <= obs)) / (1 + N_PERM),
        niche_radius = G$r, n_edges = G$n_edges, estimable = TRUE, reason = NA_character_))
    }
    cat(sprintf("[A-strat] %s %-10s region_n=%d edges=%d\n", s, rg, n, if (n >= MIN_REGION_N) G$n_edges else NA))
  }
}
mA_strat <- rbindlist(res, fill = TRUE)
setorder(mA_strat, slice, region, clone_id, tam)
fwrite(mA_strat, file.path(TAB, "clone_tam_radius_perm_region_stratified.csv"))
save_obj(mA_strat, file.path(DAT, "methodA_region_stratified"))
cat("\n[stratified] rows:", nrow(mA_strat), "| estimable:", sum(mA_strat$estimable), "\n")
cat("[estimable by region]\n"); print(mA_strat[, .(rows = .N, estimable = sum(estimable)), by = region])
cat("[not estimable reasons]\n"); print(mA_strat[estimable == FALSE, .N, by = reason][order(-N)])
cat("\n[significant estimable p<0.05]\n")
print(mA_strat[estimable == TRUE & p_perm < 0.05,
               .(slice, region, clone_id, tam, n_clone, n_pos, enrichment, z, p_perm)][order(p_perm)])
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 方法B（step5）：top25% TAM 富集 + 同spot/kNN 指标 + 层内 clone 标签置换（仅 All/Core/Interface）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
set.seed(SEED)
cl_map <- unique(d[, .(key, clone_id, clone)])
res_same <- list(); res_nbr <- list(); blk_meta <- list()
METRICS <- c("samespot_C1QC","samespot_SPP1","neighbor_C1QC_k6","neighbor_SPP1_k6",
             "neighbor_C1QC_k18","neighbor_SPP1_k18")
for (s in SLICES) {
  us <- univ[slice == s]; setorder(us, key); n <- nrow(us)
  D <- as.matrix(dist(as.matrix(us[, .(x_reg, y_reg)]))); diag(D) <- Inf
  ord <- t(apply(D, 1, function(z) order(z, method = "radix")[seq_len(max(K_MAIN, K_SENS))]))
  rm(D); invisible(gc(FALSE))
  fl <- list(); thr <- setNames(numeric(length(TAMS)), TAMS)
  for (tm in TAMS) {
    q <- unname(quantile(us[[tm]], ENRICH_Q, na.rm = TRUE, type = 7))
    thr[tm] <- q; fl[[tm]] <- as.numeric(us[[tm]] >= q)
  }
  V <- cbind(samespot_C1QC    = fl[["C1QC_TAM"]],
             samespot_SPP1    = fl[["SPP1_TAM"]],
             neighbor_C1QC_k6 = rowMeans(matrix(fl[["C1QC_TAM"]][ord[, seq_len(K_MAIN)]], nrow = n)),
             neighbor_SPP1_k6 = rowMeans(matrix(fl[["SPP1_TAM"]][ord[, seq_len(K_MAIN)]], nrow = n)),
             neighbor_C1QC_k18 = rowMeans(matrix(fl[["C1QC_TAM"]][ord[, seq_len(K_SENS)]], nrow = n)),
             neighbor_SPP1_k18 = rowMeans(matrix(fl[["SPP1_TAM"]][ord[, seq_len(K_SENS)]], nrow = n)))
  is_obs <- us$key %in% d$key
  cid <- cl_map$clone_id[match(us$key, cl_map$key)]
  lev_all <- sort(unique(d[slice == s, clone_id]))
  for (zs in c("All", ZONES)) {
    pool <- is_obs & if (zs == "All") rep(TRUE, n) else us$region == zs
    n_pool <- sum(pool)
    levs <- lev_all[lev_all %in% cid[pool]]
    sizes <- as.numeric(table(factor(cid[pool], levels = levs)))
    npool_tam <- sapply(TAMS, function(tm) sum(fl[[tm]][pool]))
    nuniv_tam <- sapply(TAMS, function(tm) if (zs == "All") sum(fl[[tm]]) else sum(fl[[tm]][us$region == zs]))
    blk_meta[[length(blk_meta) + 1]] <- data.table(
      slice = s, zone = zs, n_pool = n_pool, n_clones = length(levs),
      n_C1QC_pos_pool = npool_tam[["C1QC_TAM"]], n_SPP1_pos_pool = npool_tam[["SPP1_TAM"]])
    stat <- NULL
    if (n_pool >= MIN_POOL_B && length(levs) > 0) {
      Vp <- V[pool, , drop = FALSE]; cl_fac <- factor(cid[pool], levels = levs)
      obs <- rowsum(Vp, cl_fac)[levs, , drop = FALSE] / sizes
      bg  <- colMeans(Vp)
      null <- array(NA_real_, dim = c(N_PERM, length(levs), ncol(V)))
      for (p in seq_len(N_PERM)) {
        pl <- sample(cl_fac)
        null[p, , ] <- rowsum(Vp, pl)[levs, , drop = FALSE] / sizes
      }
      stat <- array(NA_real_, dim = c(length(levs), ncol(V), 6),
                    dimnames = list(levs, METRICS, c("obs","bg","roe","z","p_enrich","p_deplete")))
      for (ci in seq_along(levs)) for (j in seq_along(METRICS)) {
        o <- obs[ci, j]; nl <- null[, ci, j]
        stat[ci, j, "obs"] <- o; stat[ci, j, "bg"] <- bg[j]
        stat[ci, j, "roe"] <- o / bg[j]
        stat[ci, j, "z"] <- if (sd(nl) > 0) (o - mean(nl)) / sd(nl) else NA_real_
        stat[ci, j, "p_enrich"]  <- (sum(nl >= o) + 1) / (N_PERM + 1)
        stat[ci, j, "p_deplete"] <- (sum(nl <= o) + 1) / (N_PERM + 1)
      }
    }
    for (cl in lev_all) {
      for (tm in TAMS) {
        n_cl <- sum(cid[pool] == cl, na.rm = TRUE)
        npos_p <- npool_tam[[tm]]; npos_u <- nuniv_tam[[tm]]
        base_reason <- NA_character_
        if (n_pool < MIN_POOL_B) base_reason <- sprintf("stratum_pool=%d<30", n_pool)
        else if (n_pool < MIN_REGION_N) base_reason <- sprintf("stratum_pool=%d<100", n_pool)
        else if (n_cl == 0) base_reason <- "clone absent in stratum (n=0)"
        else if (n_cl < MIN_CLONE_N) base_reason <- sprintf("clone_n=%d<20", n_cl)
        else if (npos_p < MIN_POS_N) base_reason <- sprintf("TAMpos_n=%d<20", npos_p)
        est <- is.na(base_reason)
        j_same <- if (tm == "C1QC_TAM") "samespot_C1QC" else "samespot_SPP1"
        for (metric in c(j_same, sub("samespot", "neighbor", j_same) |> (\(x) c(paste0(x,"_k6"), paste0(x,"_k18")))())) {
          st <- if (!is.null(stat) && cl %in% levs) stat[cl, metric, ] else
                setNames(rep(NA_real_, 6), c("obs","bg","roe","z","p_enrich","p_deplete"))
          row <- data.table(patient = SLICE2P[[s]], slice = s, clone_id = cl,
                            clone = d$clone[match(cl, d$clone_id)], tam = TAM_LAB[[tm]],
                            zone = zs, n_spots_clone = n_cl, n_pool = n_pool,
                            n_pos_pool = npos_p, n_pos_universe = npos_u,
                            roe = st[["roe"]], z = st[["z"]],
                            p_enrich = st[["p_enrich"]], p_deplete = st[["p_deplete"]],
                            low_n = (n_cl < LOW_N_B), estimable = est, reason = base_reason)
          if (grepl("^samespot", metric)) {
            row[, `:=`(obs_frac = st[["obs"]], bg_frac = st[["bg"]])]
            res_same[[length(res_same) + 1]] <- row
          } else {
            row[, `:=`(k = as.integer(sub(".*_k", "", metric)),
                       obs_density = st[["obs"]], bg_density = st[["bg"]])]
            res_nbr[[length(res_nbr) + 1]] <- row
          }
        }
      }
    }
    cat(sprintf("[B] %s %-10s n_pool=%d clones=%d\n", s, zs, n_pool, length(levs)))
  }
}
mB_same <- rbindlist(res_same); mB_nbr <- rbindlist(res_nbr); blk <- rbindlist(blk_meta)
mB_same[estimable == FALSE, c("obs_frac","roe","z","p_enrich","p_deplete") := NA_real_]
mB_nbr[estimable == FALSE,  c("obs_density","roe","z","p_enrich","p_deplete") := NA_real_]
setorder(mB_same, slice, clone_id, tam, zone); setorder(mB_nbr, k, slice, clone_id, tam, zone)
fwrite(mB_same, file.path(TAB, "clone_tam_coloc_samespot.csv"))
fwrite(mB_nbr,  file.path(TAB, "clone_tam_coloc_neighbor.csv"))
save_obj(list(same = mB_same, nbr = mB_nbr, block = blk), file.path(DAT, "methodB_coloc"))
cat("\n[methodB] samespot rows:", nrow(mB_same), "| neighbor rows:", nrow(mB_nbr),
    "| zones:", paste(unique(mB_same$zone), collapse=" / "),
    "| non-est with roe:", mB_same[estimable == FALSE & !is.na(roe), .N], "\n")
cat("[not estimable reasons]\n"); print(mB_same[estimable == FALSE, .N, by = reason][order(-N)])
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 汇总方法B逐片结果、阈值 QC、两方法逐 clone 方向对照（仅 All/Core/Interface）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
tam_key_of <- function(x) ifelse(x == "C1Q TAMs", "C1QC_TAM", "SPP1_TAM")
dir_of <- function(est, roe, p_enrich, p_deplete) {
  out <- rep("not estimable", length(est))
  for (i in which(est)) {
    p_dir <- if (!is.na(roe[i]) && roe[i] > 1) p_enrich[i] else p_deplete[i]
    out[i] <- if (!is.na(p_dir) && p_dir < 0.05) (if (roe[i] > 1) "enriched" else "depleted") else "not significant"
  }
  out
}
enr_qc <- rbindlist(lapply(SLICES, function(s) {
  us <- univ[slice == s]
  rbindlist(lapply(TAMS, function(tm) {
    q <- unname(quantile(us[[tm]], ENRICH_Q, na.rm = TRUE, type = 7))
    data.table(slice = s, patient = SLICE2P[[s]], tam = TAM_LAB[[tm]], tam_key = tm,
               n_universe = nrow(us), threshold_q75 = q,
               n_enriched = sum(us[[tm]] >= q), enriched_pct = round(mean(us[[tm]] >= q) * 100, 2),
               bg_frac_universe = mean(us[[tm]] >= q))
  }))
}))
sm <- mB_same[zone == "All", .(samespot_bg_frac = bg_frac[1], n_est_same = sum(estimable)),
              by = .(slice, tam_key = tam_key_of(tam))]
nb6 <- mB_nbr[zone == "All" & k == K_MAIN, .(
  k6_bg_density = bg_density[1],
  n_clone_enriched_p05 = sum(estimable & p_enrich < 0.05),
  n_clone_depleted_p05 = sum(estimable & p_deplete < 0.05),
  roe_min_k6 = suppressWarnings(min(roe[estimable], na.rm = TRUE)),
  roe_max_k6 = suppressWarnings(max(roe[estimable], na.rm = TRUE))),
  by = .(slice, tam_key = tam_key_of(tam))]
nb6[!is.finite(roe_min_k6), roe_min_k6 := NA_real_]
nb6[!is.finite(roe_max_k6), roe_max_k6 := NA_real_]
cons <- mB_nbr[zone == "All" & k == K_MAIN & estimable == TRUE, .(
  direction_consensus = if (all(roe > 1)) "all_enriched" else if (all(roe < 1)) "all_depleted" else "mixed"),
  by = .(slice, tam_key = tam_key_of(tam))]
zone_est <- mB_nbr[zone != "All", .(n_estimable_zone_cells = sum(estimable),
                                    estimable_zones = paste(sort(unique(zone[estimable])), collapse = ";")),
                   by = .(slice, tam_key = tam_key_of(tam))]
summ <- merge(enr_qc, sm,  by = c("slice","tam_key"), all.x = TRUE)
summ <- merge(summ, nb6, by = c("slice","tam_key"), all.x = TRUE)
summ <- merge(summ, cons, by = c("slice","tam_key"), all.x = TRUE)
summ <- merge(summ, zone_est, by = c("slice","tam_key"), all.x = TRUE)
summ[is.na(direction_consensus), direction_consensus := "not estimable"]
summ[, enriched_pct_flag := enriched_pct > 30]
setcolorder(summ, c("slice","patient","tam","n_universe","threshold_q75","n_enriched","enriched_pct",
                    "enriched_pct_flag","bg_frac_universe","samespot_bg_frac","k6_bg_density",
                    "n_clone_enriched_p05","n_clone_depleted_p05","roe_min_k6","roe_max_k6",
                    "direction_consensus","n_est_same","n_estimable_zone_cells","estimable_zones"))
setorder(summ, slice, tam)
fwrite(summ, file.path(TAB, "coloc_summary_by_slice.csv"))
cat("[coloc_summary_by_slice]\n"); print(summ)

A  <- mA_radius[, .(slice, clone_id, tam, n_clone, n_pos, roe_A = enrichment, z_A = z,
                    pA_perm = p_perm, pA_dep = p_deplete, est_A = estimable)]
B1 <- mB_same[zone == "All", .(slice, clone_id, tam, roe_Bs = roe, z_Bs = z,
                               pBs_enr = p_enrich, pBs_dep = p_deplete, est_Bs = estimable,
                               n_pos_pool, n_pool)]
B2 <- mB_nbr[zone == "All" & k == K_MAIN, .(slice, clone_id, tam, roe_B6 = roe, z_B6 = z,
                                            pB6_enr = p_enrich, pB6_dep = p_deplete, est_B6 = estimable)]
cmp <- Reduce(function(a, b) merge(a, b, by = c("slice","clone_id","tam"), all = TRUE), list(A, B1, B2))
cmp[, dir_A  := dir_of(est_A,  roe_A,  pA_perm,  pA_dep)]
cmp[, dir_Bs := dir_of(est_Bs, roe_Bs, pBs_enr, pBs_dep)]
cmp[, dir_B6 := dir_of(est_B6, roe_B6, pB6_enr, pB6_dep)]
cmp[, sign_agree_A_B6 := ifelse(est_A & est_B6, sign(roe_A) == sign(roe_B6), NA)]
cmp[, sign_agree_A_Bs := ifelse(est_A & est_Bs, sign(roe_A) == sign(roe_Bs), NA)]
cmp[, dir_agree_A_B6  := ifelse(est_A & est_B6, dir_A == dir_B6, NA)]
cmp[, dir_agree_A_Bs  := ifelse(est_A & est_Bs, dir_A == dir_Bs, NA)]
cmp[, sig_agree_A_B6  := ifelse(dir_A %in% c("enriched","depleted") & dir_B6 %in% c("enriched","depleted"),
                                dir_A == dir_B6, NA)]
setorder(cmp, slice, clone_id, tam)
fwrite(cmp, file.path(TAB, "clone_tam_direction_comparison.csv"))
cmp_stats <- list(n_pairs = nrow(cmp),
  comparable_A_B6 = sum(!is.na(cmp$sign_agree_A_B6)), agree_A_B6 = sum(cmp$sign_agree_A_B6, na.rm = TRUE),
  rate_A_B6 = round(mean(cmp$sign_agree_A_B6, na.rm = TRUE) * 100, 1),
  comparable_A_Bs = sum(!is.na(cmp$sign_agree_A_Bs)), agree_A_Bs = sum(cmp$sign_agree_A_Bs, na.rm = TRUE),
  rate_A_Bs = round(mean(cmp$sign_agree_A_Bs, na.rm = TRUE) * 100, 1),
  dir_agree_A_B6 = sum(cmp$dir_agree_A_B6, na.rm = TRUE), dir_agree_A_Bs = sum(cmp$dir_agree_A_Bs, na.rm = TRUE),
  both_sig_A_B6 = sum(!is.na(cmp$sig_agree_A_B6)), same_dir_sig = sum(cmp$sig_agree_A_B6, na.rm = TRUE))
print(cmp_stats)
save_obj(list(summary = summ, comparison = cmp, cmp_stats = cmp_stats), file.path(DAT, "methodB_summary_comparison"))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 生成每个切片 clone x TAM 热图（All + Tumor Core + Interface；Stroma 已移除）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
ZONE_LEVELS <- c("All", "Tumor Core", "Interface")
clone_levels <- function(x) { u <- unique(as.character(x)); u[order(as.numeric(gsub("\\D", "", u)))] }
pdA <- rbind(
  mA_radius[, .(slice, region = "All", clone_id, clone, tam, roe = enrichment,
                estimable, p_enrich = p_perm, p_deplete)],
  mA_strat[, .(slice, region, clone_id, clone, tam, roe = enrichment,
               estimable, p_enrich = p_perm, p_deplete)])
pdB <- mB_nbr[k == K_MAIN, .(slice, region = zone, clone_id, clone, tam, roe,
                             estimable, p_enrich, p_deplete)]
prep_heat <- function(dat) {
  dat <- copy(dat)
  dat[, region := factor(region, levels = ZONE_LEVELS)]
  dat[, `:=`(dir_p = ifelse(!is.na(roe) & roe > 1, p_enrich, p_deplete))]
  dat[, `:=`(sig = estimable & !is.na(dir_p) & dir_p < 0.05)]
  dat[, fill := pmax(pmin(log2(pmax(roe, 1e-3)), 2), -2)]
  dat[is.na(roe) | !estimable, fill := NA_real_]
  dat[, lab := fifelse(is.na(roe) | !estimable, "n/a",
                       sprintf("%.2f%s", roe, fifelse(sig, "*", "")))]
  dat[, clone := factor(clone, levels = clone_levels(clone))]
  dat[]
}
heat_plot <- function(dat, title) {
  ggplot(dat, aes(clone, tam, fill = fill)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = lab), size = 2.1, colour = "grey10") +
    scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH, midpoint = 0,
                         limits = c(-2, 2), oob = scales::squish,
                         na.value = "grey85", name = "log2(Ro/e)") +
    facet_grid(region ~ ., scales = "free_y", space = "free_y") +
    labs(title = title, subtitle = "grey n/a = not estimable; * one-sided permutation p < 0.05; Stroma removed",
         x = "Clone", y = "TAM subtype") +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.text = element_text(size = 7, colour = "black"), axis.ticks = element_blank(),
          strip.text = element_text(size = 7, face = "bold", hjust = 0.5),
          plot.title = element_text(hjust = 0.5, size = 8),
          plot.subtitle = element_text(hjust = 0.5, size = 6.5, colour = "grey35"),
          legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
          legend.key.height = unit(4, "mm"))
}
pdA <- prep_heat(pdA); pdB <- prep_heat(pdB)
figA <- character(0); figB <- character(0)
for (s in SLICES) {
  ncl <- uniqueN(pdA[slice == s, clone_id]); w <- max(4, 0.9 * ncl + 2.4); h <- 3.8
  pA <- heat_plot(pdA[slice == s], sprintf("Clone x TAM co-localization - %s (%s) - Method A radius", s, SLICE2P[[s]]))
  pB <- heat_plot(pdB[slice == s], sprintf("Clone x TAM co-localization - %s (%s) - Method B k=%d", s, SLICE2P[[s]], K_MAIN))
  figA <- c(figA, save_fig(pA, sprintf("heatmap_methodA_%s", s), w, h))
  figB <- c(figB, save_fig(pB, sprintf("heatmap_methodB_%s", s), w, h))
}
save_obj(list(methodA = figA, methodB = figB), file.path(DAT, "figure_stems_per_slice"))
cat("[figures] methodA:", length(figA), "| methodB:", length(figB),
    "| zones:", paste(ZONE_LEVELS, collapse=" / "), "| per-slice size:", w, "x", h, "in\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 094 系列专项面板：仅 Tumor Core / Interface（2 区），4 片 x 2 区 = 8 panel；每片独立 clone 轴
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
S094 <- c("094A","094B","094C","094D")
PANEL_ZONES <- c("Tumor Core", "Interface")
heat_panel <- function(dat, title) {
  dat <- copy(dat)[slice %in% S094 & region %in% PANEL_ZONES]
  dat[, slice := factor(slice, levels = S094)]
  dat[, clone_id := factor(clone_id, levels = sort(unique(clone_id)))]
  dat[, panel_id := factor(paste(as.character(slice), as.character(region), sep = "\n"),
                           levels = do.call(c, lapply(S094, function(s) paste(s, PANEL_ZONES, sep = "\n"))))]
  ggplot() +
    geom_blank(data = dat, aes(clone_id, tam)) +
    geom_tile(data = dat[estimable == TRUE], aes(clone_id, tam, fill = fill),
              color = "white", linewidth = 0.25) +
    geom_tile(data = dat[estimable == FALSE], aes(clone_id, tam),
              fill = "grey85", color = "white", linewidth = 0.25) +
    geom_text(data = dat, aes(clone_id, tam, label = lab), size = 1.9, colour = "grey10") +
    scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH, midpoint = 0,
                         limits = c(-2, 2), oob = scales::squish, name = "log2(Ro/e)") +
    facet_wrap(~ panel_id, ncol = length(PANEL_ZONES), scales = "free_x") +
    labs(title = title,
         subtitle = "Tumor Core / Interface only (Stroma removed); each panel has its own clone axis; grey n/a = not estimable; * p < 0.05",
         x = "Clone (patient-prefixed)", y = "TAM subtype") +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.text = element_text(size = 6.5, colour = "black"), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.text = element_text(size = 7, face = "bold", hjust = 0.5),
          plot.title = element_text(hjust = 0.5, size = 8),
          plot.subtitle = element_text(hjust = 0.5, size = 6.0, colour = "grey35"),
          legend.title = element_text(size = 7), legend.text = element_text(size = 6.5))
}
p94A <- heat_panel(pdA, "094 series - Clone x TAM co-localization - Method A radius")
p94B <- heat_panel(pdB, sprintf("094 series - Clone x TAM co-localization - Method B k=%d", K_MAIN))
fig94 <- c(save_fig(p94A, "panel_094_methodA", 7.2, 7.8),
           save_fig(p94B, "panel_094_methodB", 7.2, 7.8))
save_obj(fig94, file.path(DAT, "figure_stems_094"))
cat("[094 panel] 4 slices x", length(PANEL_ZONES), "zones =", 4*length(PANEL_ZONES), "panels | size 7.2 x 7.8 in\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 图件完整性机械自检（Stroma 移除后）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
pngs <- sort(list.files(FIG, pattern = "\\.png$", full.names = TRUE))
chk <- rbindlist(lapply(pngs, function(f) {
  im <- magick::image_read(f); ii <- magick::image_info(im)
  arr <- as.numeric(magick::image_data(magick::image_resize(im, "800x800"))[1:3, , , drop = FALSE]) / 255
  stem <- sub("\\.png$", "", basename(f))
  data.table(stem = stem, width = ii$width, height = ii$height,
             size_kb = round(file.info(f)$size / 1024, 1),
             sd_rgb = round(sd(arr), 4), nonwhite_frac = round(mean(arr < 0.98), 4),
             pdf_pair = file.exists(file.path(FIG, paste0(stem, ".pdf"))))
}))
chk[, pass := sd_rgb > 0.02 & nonwhite_frac > 0.01 & width >= 1200 & height >= 1200 & pdf_pair]
fwrite(chk, file.path(REV, "figure_pixel_check.csv"))
cat("[Figure Wholeness]", sum(chk$pass), "/", nrow(chk), "figures pixel-verified\n")
print(chk[stem %in% c("heatmap_methodA_092A","heatmap_methodB_094C","panel_094_methodA","panel_094_methodB")])
if (any(!chk$pass)) print(chk[pass == FALSE]) else cat("[PASS] all figures non-blank and paired\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 报告用汇总（Stroma 移除后）：逐片方向、分层可估、一致性
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
dir_tab <- cmp[, .(slice, clone_id, tam, n_clone, roe_A, dir_A, roe_Bs, dir_Bs,
                   roe_B6, dir_B6, sign_agree_A_B6, dir_agree_A_B6)]
setorder(dir_tab, slice, clone_id, tam)
fwrite(dir_tab, file.path(TAB, "clone_tam_direction_by_slice.csv"))
cat("[method A stratified estimable cells]\n")
print(mA_strat[, .(n_rows = .N, n_estimable = sum(estimable)), by = region])
cat("\n[method B zone estimable cells, k=6]\n")
print(mB_nbr[k == K_MAIN & zone != "All", .(n_rows = .N, n_estimable = sum(estimable)), by = zone])
cat("\n[method A whole-slice] enriched:", mA_radius[estimable == TRUE & enrichment > 1 & p_perm < 0.05, .N],
    "| depleted:", mA_radius[estimable == TRUE & enrichment < 1 & p_deplete < 0.05, .N],
    "| not estimable:", mA_radius[estimable == FALSE, .N], "\n")
cat("[method B k6 whole-slice] enriched:", mB_nbr[zone=="All" & k==K_MAIN & estimable==TRUE & roe > 1 & p_enrich < 0.05, .N],
    "| depleted:", mB_nbr[zone=="All" & k==K_MAIN & estimable==TRUE & roe < 1 & p_deplete < 0.05, .N],
    "| not estimable:", mB_nbr[zone=="All" & k==K_MAIN & estimable==FALSE, .N], "\n")
print(cmp_stats)
cat("\n[sign-disagreement A vs B samespot]\n")
print(cmp[!is.na(sign_agree_A_Bs) & sign_agree_A_Bs == FALSE,
          .(slice, clone_id, tam, roe_A, roe_Bs, dir_A, dir_Bs)])
save_obj(list(cmp_stats = cmp_stats, dir_tab = dir_tab), file.path(DAT, "report_numbers"))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 生成两层 HTML 交付文档（结论简报 + 技术说明），图片 base64 内嵌
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

html_table <- function(dt, digits = 3) {
  d <- copy(as.data.frame(dt))
  num <- sapply(d, is.numeric)
  d[num] <- lapply(d[num], function(v) ifelse(is.na(v), "NA", formatC(v, digits = digits, format = "g")))
  hdr <- paste0("<tr>", paste0("<th>", names(d), "</th>", collapse = ""), "</tr>")
  body <- apply(d, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<table>", hdr, paste0(body, collapse = ""), "</table>")
}
img_tag <- function(stem, caption, width = "100%") {
  b64 <- base64enc::base64encode(file.path(FIG, paste0(stem, ".png")))
  sprintf('<figure><img src="data:image/png;base64,%s" style="width:%s"><figcaption>%s</figcaption></figure>',
          b64, width, caption)
}
n_spots <- nrow(d); n_slices <- uniqueN(d$slice); n_clones_u <- uniqueN(d$clone_id)
n_slice_clone <- nrow(unique(d[, .(slice, clone_id)]))
mA_est <- mA_radius[estimable == TRUE]
a_enr <- sum(mA_est$enrichment > 1 & mA_est$p_perm < 0.05)
a_dep <- sum(mA_est$enrichment < 1 & mA_est$p_deplete < 0.05)
a_ns  <- nrow(mA_est) - a_enr - a_dep
bs <- mB_same[zone == "All" & estimable == TRUE]
bs_enr <- sum(bs$roe > 1 & bs$p_enrich < 0.05); bs_dep <- sum(bs$roe < 1 & bs$p_deplete < 0.05)
b6 <- mB_nbr[zone == "All" & k == K_MAIN & estimable == TRUE]
b_enr <- sum(b6$roe > 1 & b6$p_enrich < 0.05); b_dep <- sum(b6$roe < 1 & b6$p_deplete < 0.05)
b6_ns <- nrow(b6) - b_enr - b_dep
strat_A <- mA_strat[, .(estimable = sum(estimable), rows = .N), by = region]
strat_B <- mB_nbr[k == K_MAIN & zone != "All", .(estimable = sum(estimable), rows = .N), by = zone]
inp_show <- input_summary[, .(slice, patient, n_clones, n_clone_spots, n_TumorCore, n_Interface, n_Stroma,
                              n_universe, join_rate_pct)]

lead <- sprintf(paste0(
  "对 10 个切片、5 位病人的 <b>%s 个全量观测 spot</b>（无下采样），用两套独立方法逐片计算 clone × TAM 共定位。",
  "两方法方向高度一致：k=6 邻域 <b>%d/%d（%.1f%%）</b>、同 spot %d/%d（%.1f%%）；27 对双方显著的 clone×TAM 方向 100%% 相同。",
  "clone 层面的 TAM 富集/缺失是<b>切片特异</b>的，未发现跨片一致的单一方向；三区分层中 Stroma 全部不可估（每片 clone spot 0–1 个）。"),
  format(n_spots, big.mark = ","), cmp_stats$agree_A_B6, cmp_stats$comparable_A_B6, cmp_stats$rate_A_B6,
  cmp_stats$agree_A_Bs, cmp_stats$comparable_A_Bs, cmp_stats$rate_A_Bs)

blk_keys <- sprintf(paste0(
  '<ul><li>输入 %s spot / 10 片 / %d 个唯一 clone → <code>coloc_input_summary.csv</code> 与 <code>methods_note.md</code>。</li>',
  '<li>方法A 整片：%d 富集 / %d 缺失 / %d 不显著 / %d 不可估 → <code>clone_tam_radius_perm_by_slice.csv</code>。</li>',
  '<li>方法B k=6 整片：%d 富集 / %d 缺失 / %d 不显著 / %d 不可估 → <code>clone_tam_coloc_neighbor.csv</code>（k=6, zone=All）。</li>',
  '<li>分层可估：方法A Tumor Core 40 / Interface 16 / Stroma 0；方法B k=6 Tumor Core 43 / Interface 24 / Stroma 0。</li>',
  '<li>方向一致率：k=6 39/39（100%%），同 spot 38/39（97.4%%），显著方向 27/27 → <code>clone_tam_direction_comparison.csv</code>。</li></ul>'),
  format(n_spots, big.mark = ","), n_clones_u, a_enr, a_dep, a_ns, sum(!mA_radius$estimable),
  b_enr, b_dep, b6_ns, sum(!mB_nbr[zone == "All" & k == K_MAIN, estimable]))

html <- c(
'<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">',
'<title>clone × TAM 共定位汇总</title>',
'<style>
body{font-family:-apple-system,"Segoe UI","Noto Sans CJK SC",sans-serif;max-width:1180px;margin:32px auto;padding:0 22px;color:#1a1a1a;line-height:1.65;background:#fff}
h1{font-size:24px;border-bottom:2px solid #333;padding-bottom:8px}
h2{font-size:18px;margin-top:30px;border-left:4px solid #3C5488;padding-left:10px}
h3{font-size:15px;margin-top:20px}
table{border-collapse:collapse;margin:12px 0;font-size:12.5px}
th,td{border:1px solid #cfcfcf;padding:4px 7px;text-align:left}
th{background:#f0f2f5;font-weight:600}
.lead{background:#f7f9fc;border-left:4px solid #4DBBD5;padding:10px 14px;font-size:15px}
.note{background:#fff8e6;border-left:4px solid #E6AB02;padding:8px 12px;font-size:13px}
figure{margin:16px 0;text-align:center}
figcaption{font-size:12px;color:#555;margin-top:4px}
hr{margin:36px 0;border:none;border-top:1px dashed #bbb}
code{background:#f3f3f3;padding:1px 4px;border-radius:3px;font-size:12px}
</style></head><body>',
'<h1>clone × TAM 共定位（patient_ALL，全量观测 spot，10 片）</h1>',
sprintf('<p class="lead">%s</p>', lead),
'<h2>① 一句话结论</h2>',
sprintf(paste0('<p>在每个切片内部，clone 与 C1Q/SPP1 TAM 的空间关系主要表现为<b>切片特异的富集或缺失</b>：',
  '方法A 在 42 个可估的整片 clone×TAM 组合中检出 %d 个富集、%d 个缺失；方法B（k=6）在 45 个可估组合中检出 %d 个富集、%d 个缺失。',
  '两方法对方向的判断高度一致，但强度与显著性并不相同，因此分别成表、不合并为单一结论。</p>'),
  a_enr, a_dep, b_enr, b_dep),
'<h2>② 方法要点（方向参照系）</h2>',
'<ul>',
'<li><b>方法A（step4）</b>：clone spot 为分析域；TAM 阳性 = 片内 ≥ max(1%, 2×片内均值)；半径 r = 2×median 1-NN；统计量为半径邻域内 clone–TAM 阳性对的观测边数；置换 <b>TAM 阳性标签</b> 1000 次；另给 k=6 邻域占比。</li>',
'<li><b>方法B（step5）</b>：全部组织 spot 定义 TAM 富集（片内 top 25%）并算 kNN；统计量为 clone 内指标均值 / 层内背景（Ro/e）；置换 <b>clone 标签</b>（保持 clone 大小）1000 次。</li>',
'<li><b>方向参照系</b>：Ro/e &gt; 1 = 该 clone 所在位置的 TAM 比同层背景更富集；&lt; 1 = 更少（缺失）；<code>*</code> = 片内单侧置换 p &lt; 0.05。方法A 的 “富集/缺失” 相对于 clone spot 内部的 TAM 随机分布；方法B 相对于层内 clone 标签随机分布。</li>',
'<li>随机种子 42；所有统计均<b>逐片独立</b>，不合并、不跨片置换；clone 编号带病人前缀，不跨病人比较。</li>',
'</ul>',
'<h2>③ 关键结果</h2>',
'<h3>3.1 输入与关联率（逐片 100%）</h3>',
html_table(inp_show, 2),
sprintf(paste0('<p>clone 文件仅来自 <code>patient_ALL/clone/</code>：P9 2,474、P12 2,164、P13 1,485、P14 1,617、P15 1,938 行，',
  '合计 <b>%s</b> spot；%d 片、%d 个病人、%d 个 “切片×clone” 条目、%d 个唯一 clone。三区来自 project_singlecell V2 底表的 <code>region</code> 列，TAM 比例来自 step3 的 <code>spatial_tme_groups.qs</code>。</p>'),
  format(n_spots, big.mark = ","), n_slices, uniqueN(d$patient), n_slice_clone, n_clones_u),
'<h3>3.2 整片（All）方向计数</h3>',
sprintf('<p><b>方法A 半径邻域</b>：可估 %d / 54，富集 %d、缺失 %d、不显著 %d，不可估 %d（clone 或 TAM 阳性 &lt; 20）。<br><b>方法B 同 spot</b>：可估 %d / 54，富集 %d、缺失 %d、不显著 %d。<br><b>方法B k=6</b>：可估 %d / 54，富集 %d、缺失 %d、不显著 %d。</p>',
  nrow(mA_est), a_enr, a_dep, a_ns, sum(!mA_radius$estimable), nrow(bs), bs_enr, bs_dep, nrow(bs)-bs_enr-bs_dep,
  nrow(b6), b_enr, b_dep, b6_ns),
'<h3>3.3 三区分层可估性</h3>',
html_table(strat_A, 0), '<p style="font-size:12px;color:#555">表：方法A 分层（每区 54 格 = 10 片 × 该片 clone × 2 TAM）。</p>',
html_table(strat_B, 0), '<p style="font-size:12px;color:#555">表：方法B k=6 分层。</p>',
'<div class="note">Stroma 在 10 片中 clone spot 仅 0–1 个（093C 为 1），全部低于 region_n≥100，<b>不可估</b>，表中不填 0。Interface 在 092A/092B/094C/094D 的 clone spot 为 0/4/10/16，亦不可估；其余 6 片可估。</div>',
'<h3>3.4 两方法方向一致性</h3>',
sprintf(paste0('<ul><li>方法A vs 方法B(k=6)：可比 %d 对，方向一致 <b>%d 对（%.1f%%）</b>。</li>',
  '<li>方法A vs 方法B(同 spot)：可比 %d 对，方向一致 %d 对（%.1f%%）。</li>',
  '<li>双方均显著的 %d 对中，方向相同 <b>%d 对</b>。</li></ul>'),
  cmp_stats$comparable_A_B6, cmp_stats$agree_A_B6, cmp_stats$rate_A_B6,
  cmp_stats$comparable_A_Bs, cmp_stats$agree_A_Bs, cmp_stats$rate_A_Bs,
  cmp_stats$both_sig_A_B6, cmp_stats$same_dir_sig),
'<h3>3.5 094 系列专项</h3>',
'<p>094A/094B（P13，2 个 clone）与 094C/094D（P14，5–6 个 clone）显示同一病人内部 clone 间的 TAM 偏好相反：例如 094A_Clone1 富集 C1Q TAM、Clone2 缺失 C1Q TAM；094C_Clone1/Clone3 富集 C1Q TAM，而 Clone2/Clone4 缺失。094D 有 4 个 clone 在该区低于门槛，如实标注不可估。</p>',
img_tag("panel_094_methodA", "094 系列：方法A 半径邻域，log2(Ro/e)；n/a = 不可估（按切片分面，clone 带病人前缀）", "100%"),
img_tag("panel_094_methodB", "094 系列：方法B k=6 邻域，log2(Ro/e)；n/a = 不可估", "100%"),
'<h3>3.6 逐片热图（示例）</h3>',
img_tag("heatmap_methodA_094C", "094C（P14）方法A：All + 三区分面；每格数字为 Ro/e，* 为 p<0.05", "78%"),
img_tag("heatmap_methodB_094C", "094C（P14）方法B k=6：All + 三区分面", "78%"),
'<h2>④ 局限与注意</h2>',
'<ul>',
'<li>所有 p 为<b>片内名义 p</b>，未做跨片 FDR 或 p 值合并；切片/病人之间不做显著性汇总。</li>',
'<li>两方法的空模型与 TAM 阳性定义不同（方法A：≥max(1%,2×均值)；方法B：top 25%），因此强度不可直接互比，只比较方向。</li>',
'<li>方法A 在 TAM 阳性 spot 很少的切片（如 092B C1QC 仅 6 个）会出现极端 Ro/e，已用 n_pos≥20 标记不可估。</li>',
'<li>小 clone（n=3/12/22）在部分层不可估；Stroma 全片不可估；这些格在表中保留并写明原因，不填 0、不画空面板。</li>',
'<li>RCTD 比例全 &gt; 0，故 “富集” 是片内相对分位，不代表某类细胞绝对存在/缺失。</li>',
'<li>clone 由 CNV 定义；clone 与肿瘤上皮比例存在同向混杂（此前分析已报告），本步未重新校正。</li>',
'</ul>',
'<h2>⑤ 图表索引</h2>',
'<ul>',
'<li><code>tables/coloc_input_summary.csv</code> — 切片 × 病人 × clone 数 × spot 数 × 三区 spot 数 × 关联率</li>',
'<li><code>tables/clone_tam_radius_perm_by_slice.csv</code> — 方法A 整片：obs/null/z/p</li>',
'<li><code>tables/clone_tam_knn6_by_slice.csv</code> — 方法A k=6 邻域</li>',
'<li><code>tables/clone_tam_radius_perm_region_stratified.csv</code> — 方法A 三区分层（含 n_spots 与不可估原因）</li>',
'<li><code>tables/clone_tam_coloc_samespot.csv</code> — 方法B 同 spot Ro/e</li>',
'<li><code>tables/clone_tam_coloc_neighbor.csv</code> — 方法B k=6/18 × 切片 × clone × TAM × zone</li>',
'<li><code>tables/coloc_summary_by_slice.csv</code> — 方法B 逐片汇总与阈值 QC</li>',
'<li><code>tables/clone_tam_direction_comparison.csv</code> / <code>clone_tam_direction_by_slice.csv</code> — 两方法逐 clone 方向对照</li>',
'<li><code>tables/methods_note.md</code> — 参数、门槛、两方法差异</li>',
'<li><code>figures/heatmap_methodA_&lt;slice&gt;.{png,pdf}</code>、<code>heatmap_methodB_&lt;slice&gt;.{png,pdf}</code> — 每片两方法各一张</li>',
'<li><code>figures/panel_094_methodA.{png,pdf}</code>、<code>panel_094_methodB.{png,pdf}</code> — 094 系列专项</li>',
'</ul>',
'<hr>',
'<h2>技术说明书</h2>',
'<h3>STEP 01 — clone × TAM 共定位（逐切片、三区分层）</h3>',
'<p>配套 notebook：<code>notebooks/step01_tam_clone_coloc.ipynb</code>（R 4.3.3，singlecell_r）</p>',
'<h3>1. 本步回答的问题</h3>',
'<ol><li>在每个切片内部，各 clone 的 spot 是否比同片/同区背景更靠近 C1Q 或 SPP1 TAM？</li>',
'<li>这种关系在 Tumor Core / Interface / Stroma 内是否依然存在？</li>',
'<li>两套独立方法给出的方向是否一致？</li></ol>',
'<h3>2. 方法流程</h3>',
'<ol>',
'<li>读入 5 个 clone_membership.csv（仅 patient_ALL/clone），逐片与 project_singlecell 三区表、step3 TAM 比例按 key 合并；逐片关联率 100%。</li>',
'<li>方法A：逐片建 clone spot 半径图（r=2×median 1-NN），定义 TAM 阳性，计算观测边数，置换 1000 次；k=6 邻域同样处理；分层版在每区重建图并按区内门槛判定可估性。</li>',
'<li>方法B：逐片在全部组织 spot 上定义 top-25% TAM 富集 spot，算同 spot 与 k=6/18 邻域指标；在 All/三区 pool 内置换 clone 标签 1000 次，计算 Ro/e 与单侧 p。</li>',
'<li>汇总：逐片方向计数、分层可估格数、两方法方向一致性；出每片热图与 094 面板。</li>',
'</ol>',
'<h3>3. 产物逐项说明</h3>',
'<ul>',
'<li><b>clone_tam_radius_perm_by_slice.csv</b>：每行 = 切片 × clone × TAM。obs_edges/null_mean_edges/enrichment/z/p_perm/p_deplete；n_clone、n_pos、niche_radius、n_edges、low_n、estimable。</li>',
'<li><b>clone_tam_knn6_by_slice.csv</b>：同上，但统计量是 clone spot 的 6 近邻中 TAM 阳性占比（obs_neighbor_frac）。</li>',
'<li><b>clone_tam_radius_perm_region_stratified.csv</b>：每行 = 切片 × 区 × clone × TAM，共 162 行；estimable=FALSE 的 reason 写明 region_n/clone_n/TAMpos_n 哪一项不足。</li>',
'<li><b>clone_tam_coloc_samespot.csv</b>：每行 = 切片 × clone × TAM × zone（All/三区）。obs_frac、bg_frac、roe、p_enrich、p_deplete、n_spots_clone、n_pool、n_pos_pool、n_pos_universe、low_n、estimable、reason。</li>',
'<li><b>clone_tam_coloc_neighbor.csv</b>：同上，k=6 与 k=18 各一行；obs_density、bg_density、roe。</li>',
'<li><b>coloc_summary_by_slice.csv</b>：逐片 × TAM 的阈值（q75）、实际富集比例（25.00–25.06%）、背景、显著 clone 数、Ro/e 范围、方向一致性、可估区。</li>',
'<li><b>热图</b>：x=clone，y=TAM，fill=log2(Ro/e) 截断 [-2,2]，蓝-白-红；灰格 “n/a” = 不可估；数字为 Ro/e，* 为 p<0.05。</li>',
'</ul>',
'<h3>4. 关键结论与数字出处</h3>',
blk_keys,
'<h3>5. 使用注意</h3>',
'<ul><li>不要跨切片/跨病人合并 p 值；本步未做 FDR。</li>',
'<li>方法A 与方法B 的阳性定义不同，强度不可直接互比；只比较方向与可估性。</li>',
'<li>标注 “不可估” 的格是样本量/阳性数不足，不是阴性结果。</li>',
'<li>方法A 在极小 clone（n=3）上会出现 Ro/e 极端值，已不可估处理。</li></ul>',
'<h2>交付边界说明</h2>',
'<p>本交付物为 <code>TAM_clone/</code> 目录下的 <code>tables/</code>、<code>figures/</code> 与本 HTML。作图与统计的判据自检记录、像素自检表和内部溯源表保存在交付目录之外的内部目录，不随交付物提供；正文中的结论与数字已在上面完整给出。</p>',
'</body></html>')

out_html <- file.path(OUT, "TAM_clone_汇总.html")
writeLines(html, out_html, useBytes = TRUE)
cat("[Write] HTML:", out_html, "| size:", round(file.info(out_html)$size/1024, 1), "KB\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 终检：措辞修正、doc-disk 对账、复核重绘的 094 面板、写内部溯源表
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

mn_f <- file.path(TAB, "methods_note.md"); html_f <- file.path(OUT, "TAM_clone_汇总.html")
x <- readLines(mn_f, encoding = "UTF-8")
x <- gsub("无任何下采样", "未做任何抽样", x, fixed = TRUE)
writeLines(x, mn_f, useBytes = TRUE)
h <- readLines(html_f, encoding = "UTF-8", warn = FALSE)
h <- gsub("（无下采样）", "（全量，未做任何抽样）", h, fixed = TRUE)
writeLines(h, html_f, useBytes = TRUE)
cat("[wording] 300/slice in methods_note:", any(grepl("300/slice", x, fixed = TRUE)),
    "| in html:", any(grepl("300/slice", h, fixed = TRUE)),
    "| 降采样 in html:", any(grepl("降采样", h, fixed = TRUE)), "\n")

# 复核重绘的两个 094 面板
for (stem in c("panel_094_methodA","panel_094_methodB")) {
  f <- file.path(FIG, paste0(stem, ".png"))
  im <- magick::image_read(f); ii <- magick::image_info(im)
  arr <- as.numeric(magick::image_data(magick::image_resize(im, "800x800"))[1:3,,,drop=FALSE]) / 255
  cat(sprintf("[recheck] %s %dx%d sd=%.4f nonwhite=%.4f pdf=%s\n", stem, ii$width, ii$height,
              sd(arr), mean(arr < 0.98), file.exists(file.path(FIG, paste0(stem, ".pdf")))))
}

# 交付物清单与大小
deliv <- data.table(
  file = c(list.files(TAB, full.names = FALSE), list.files(FIG, pattern = "\\.(png|pdf)$", full.names = FALSE),
           basename(html_f), basename(file.path(OUT, "notebooks/step01_tam_clone_coloc.ipynb"))),
  dir = c(rep("tables", length(list.files(TAB))), rep("figures", length(list.files(FIG, pattern = "\\.(png|pdf)$"))),
          "root", "notebooks"))
deliv[, size_kb := round(file.info(file.path(OUT, dir, file))$size / 1024, 1)]
cat("\n[deliverables] n =", nrow(deliv), "\n"); print(deliv)

# 内部溯源表（不属于交付物）
tr <- c("# TRACEABILITY — TAM_clone 交付数字到内部证据的映射",
"",
"| 交付数字 | 值 | 证据文件 | 列 |",
"|---|---|---|---|",
sprintf("| 总 spot | %s | tables/coloc_input_summary.csv | sum(n_clone_spots) |", format(n_spots, big.mark=",")),
sprintf("| 切片数 | %d | tables/coloc_input_summary.csv | nrow |", n_slices),
sprintf("| 唯一 clone | %d | tables/clone_tam_direction_by_slice.csv | unique(clone_id) |", n_clones_u),
sprintf("| 关联率 | 100%% x 10 | tables/coloc_input_summary.csv | join_rate_pct |"),
sprintf("| 方法A 整片 富集/缺失/不显著/不可估 | %d/%d/%d/%d | tables/clone_tam_radius_perm_by_slice.csv | estimable, enrichment, p_perm, p_deplete |", a_enr, a_dep, a_ns, sum(!mA_radius$estimable)),
sprintf("| 方法B k6 整片 富集/缺失/不显著/不可估 | %d/%d/%d/%d | tables/clone_tam_coloc_neighbor.csv | k==6 & zone==All: estimable, roe, p_enrich, p_deplete |", b_enr, b_dep, b6_ns, sum(!mB_nbr[zone=="All" & k==K_MAIN, estimable])),
sprintf("| 分层可估 方法A Core/Interface/Stroma | 40/16/0 | tables/clone_tam_radius_perm_region_stratified.csv | estimable by region |"),
sprintf("| 分层可估 方法B k6 Core/Interface/Stroma | 43/24/0 | tables/clone_tam_coloc_neighbor.csv | k==6, zone!=All: estimable by zone |"),
sprintf("| 方向一致 A vs B k6 | %d/%d (%.1f%%) | tables/clone_tam_direction_comparison.csv | sign_agree_A_B6 |", cmp_stats$agree_A_B6, cmp_stats$comparable_A_B6, cmp_stats$rate_A_B6),
sprintf("| 方向一致 A vs B 同spot | %d/%d (%.1f%%) | tables/clone_tam_direction_comparison.csv | sign_agree_A_Bs |", cmp_stats$agree_A_Bs, cmp_stats$comparable_A_Bs, cmp_stats$rate_A_Bs),
sprintf("| 显著方向一致 | %d/%d | tables/clone_tam_direction_comparison.csv | sig_agree_A_B6 |", cmp_stats$same_dir_sig, cmp_stats$both_sig_A_B6),
"| 图件像素自检 | 22/22 + 复核 2/2 | review/TAM_clone/figure_pixel_check.csv | pass |",
"",
"## 输入文件行数",
paste0("- ", clone_files, " : ", clone_rows))
writeLines(tr, file.path(REV, "TRACEABILITY.md"), useBytes = TRUE)

# doc <-> disk 机械对账
script <- "/path/to/.claude/scripts/doc_disk_reconcile.py"
if (file.exists(script)) {
  out <- suppressWarnings(system2("python3", c(script, "--all", OUT), stdout = TRUE, stderr = TRUE))
  cat("\n[doc_disk_reconcile]\n"); cat(paste(out, collapse = "\n"), "\n")
  cat("[exit]", attr(out, "status"), "\n")
} else cat("[doc_disk_reconcile] script not found\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 清理汇总表中的 Inf/-Inf（无可估 clone 时）并扫描所有 CSV 是否残留 Inf/NaN 字面量
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

summ[is.infinite(roe_min_k6), roe_min_k6 := NA_real_]
summ[is.infinite(roe_max_k6), roe_max_k6 := NA_real_]
fwrite(summ, file.path(TAB, "coloc_summary_by_slice.csv"))
cat("[fix] Inf -> NA in coloc_summary_by_slice.csv\n")

csvs <- list.files(TAB, pattern = "\\.csv$", full.names = TRUE)
bad <- rbindlist(lapply(csvs, function(f) {
  x <- readLines(f, warn = FALSE)
  data.table(file = basename(f), has_inf = any(grepl("Inf|-Inf", x, fixed = TRUE)),
             has_nan = any(grepl("NaN", x, fixed = TRUE)))
}))
print(bad)
cat("[csv scan] any Inf:", any(bad$has_inf), "| any NaN:", any(bad$has_nan), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 交付文本合规检查：不得出现 review/ 文件名、300/slice、降采样；确认内嵌图数量与关键产物存在
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

h <- paste(readLines(file.path(OUT, "TAM_clone_汇总.html"), warn = FALSE), collapse = "\n")
cat("[html] review/ mentions:", grepl("review/", h, fixed = TRUE),
    "| 300/slice:", grepl("300/slice", h, fixed = TRUE),
    "| 降采样:", grepl("降采样", h, fixed = TRUE),
    "| embedded images:", lengths(regmatches(h, gregexpr("data:image/png;base64,", h, fixed = TRUE))), "\n")
mn <- paste(readLines(file.path(TAB, "methods_note.md")), collapse = "\n")
cat("[methods_note] review/ mentions:", grepl("review/", mn, fixed = TRUE),
    "| 300/slice:", grepl("300/slice", mn, fixed = TRUE),
    "| 降采样:", grepl("降采样", mn, fixed = TRUE), "\n")

need <- c("tables/coloc_input_summary.csv", "tables/clone_tam_radius_perm_by_slice.csv",
          "tables/clone_tam_knn6_by_slice.csv", "tables/clone_tam_radius_perm_region_stratified.csv",
          "tables/clone_tam_coloc_samespot.csv", "tables/clone_tam_coloc_neighbor.csv",
          "tables/coloc_summary_by_slice.csv", "tables/methods_note.md",
          "figures/panel_094_methodA.png", "figures/panel_094_methodB.png",
          "TAM_clone_汇总.html")
cat("[required outputs]\n")
print(data.table(file = need, exists = file.exists(file.path(OUT, need)),
                 size_kb = round(file.info(file.path(OUT, need))$size/1024, 1)))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 内部一致性复核：行数/唯一键/p 值下界/阈值比例/关联率
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

expA <- nrow(unique(d[, .(slice, clone_id, tam = rep(c("C1Q TAMs","SPP1 TAMs"), times = .N))]))
cat("[rows] methodA radius:", nrow(mA_radius), "| expected slice x clone x TAM:", nrow(unique(d[,.(slice,clone_id)])) * 2, "\n")
cat("[dupes] methodA:", sum(duplicated(mA_radius[, .(slice, clone_id, tam)])),
    "| methodA strat:", sum(duplicated(mA_strat[, .(slice, region, clone_id, tam)])),
    "| methodB same:", sum(duplicated(mB_same[, .(slice, clone_id, tam, zone)])),
    "| methodB nbr:", sum(duplicated(mB_nbr[, .(k, slice, clone_id, tam, zone)])), "\n")
cat("[p lower bound] 1/(1+1000) =", 1/1001, "\n")
cat("[p min] A radius:", min(mA_radius$p_perm, na.rm=TRUE), "| B same:", min(mB_same$p_enrich, na.rm=TRUE),
    "| B nbr:", min(mB_nbr$p_enrich, na.rm=TRUE), "\n")
cat("[exclusive p] A both<0.05:", mA_radius[estimable==TRUE & p_perm<0.05 & p_deplete<0.05, .N],
    "| B same both<0.05:", mB_same[estimable==TRUE & p_enrich<0.05 & p_deplete<0.05, .N],
    "| B nbr both<0.05:", mB_nbr[estimable==TRUE & p_enrich<0.05 & p_deplete<0.05, .N], "\n")
cat("[enrich pct] range:", paste(range(summ$enriched_pct), collapse = " .. "), "| any >30:", any(summ$enriched_pct > 30), "\n")
cat("[join] all 100%:", all(input_summary$join_rate_pct == 100), "\n")
cat("[estimable] A rows:", sum(mA_radius$estimable), "/", nrow(mA_radius),
    "| B same:", sum(mB_same$estimable), "/", nrow(mB_same),
    "| B nbr k6 All:", mB_nbr[zone=="All" & k==6, sum(estimable)], "/", mB_nbr[zone=="All" & k==6, .N], "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 补充：方法B 同 spot 的逐片方向计数（报告用）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
cat("[method B same-spot / All]\n")
print(dcast(cmp[, .N, by = .(slice, dir_Bs)], slice ~ dir_Bs, value.var = "N", fill = 0))
cat("\n[method A vs B same-spot sign disagreement rows]\n")
print(cmp[!is.na(sign_agree_A_Bs) & sign_agree_A_Bs == FALSE,
          .(slice, clone_id, tam, roe_A, dir_A, roe_Bs, dir_Bs, roe_B6, dir_B6)])
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 检查 HTML 内嵌图与磁盘图的时序（决定是否需要重新生成 HTML）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
fs <- file.info(file.path(FIG, c("panel_094_methodA.png","panel_094_methodB.png",
                                 "heatmap_methodA_094C.png","heatmap_methodB_094C.png",
                                 "heatmap_methodA_092A.png","heatmap_methodB_092A.png")))
fs$file <- basename(rownames(fs))
html_t <- file.info(file.path(OUT, "TAM_clone_汇总.html"))$mtime
cat("[html mtime]", format(html_t, "%H:%M:%S"), "\n")
print(fs[, c("file","mtime","size")])
cat("[stale embedded images]", any(fs$mtime > html_t), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# HTML 内嵌图刷新（094 面板顺序修复后）：仅替换 4 个 base64 数据 URI
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
html_f <- file.path(OUT, "TAM_clone_汇总.html")
txt <- paste(readLines(html_f, warn = FALSE), collapse = "\n")
stems <- c("panel_094_methodA","panel_094_methodB","heatmap_methodA_094C","heatmap_methodB_094C")
locs <- regmatches(txt, gregexpr("data:image/png;base64,[A-Za-z0-9+/=]+", txt))[[1]]
stopifnot(length(locs) == length(stems))
for (i in seq_along(stems)) {
  new <- paste0("data:image/png;base64,", base64enc::base64encode(file.path(FIG, paste0(stems[i], ".png"))))
  txt <- sub(locs[i], new, txt, fixed = TRUE)
}
writeLines(txt, html_f, useBytes = TRUE)
cat("[refreshed] size:", round(file.info(html_f)$size/1024,1), "KB | mtime:",
    format(file.info(html_f)$mtime, "%H:%M:%S"), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# HTML 刷新后合规复查 + 内嵌图 base64 数量
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
h <- paste(readLines(file.path(OUT, "TAM_clone_汇总.html"), warn = FALSE), collapse = "\n")
cat("[html] review/:", grepl("review/", h, fixed=TRUE),
    "| 300/slice:", grepl("300/slice", h, fixed=TRUE),
    "| 降采样:", grepl("降采样", h, fixed=TRUE),
    "| TODO/placeholder:", grepl("TODO|placeholder", h, ignore.case=TRUE),
    "| img:", lengths(regmatches(h, gregexpr("data:image/png;base64,", h, fixed=TRUE))),
    "| size_kb:", round(nchar(h, type="bytes")/1024,1), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 去除 base64 后再查文本占位符（避免 base64 随机子串误报）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
h <- paste(readLines(file.path(OUT, "TAM_clone_汇总.html"), warn = FALSE), collapse = "\n")
h_txt <- gsub("data:image/png;base64,[A-Za-z0-9+/=]+", "[IMG]", h)
cat("[text-only] contains TODO/placeholder:", grepl("TODO|placeholder", h_txt, ignore.case = TRUE),
    "| contains 'n/a':", grepl("n/a", h_txt, fixed = TRUE), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 诊断 094D 面板为何缺少全 n/a 的 P14_Clone3 列
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
cat("[pdB 094D rows by clone]\n"); print(pdB[slice == "094D", .N, by = clone_id][order(clone_id)])
cat("[pdB 094D clone_id levels]\n"); print(sort(unique(as.character(pdB[slice == "094D"]$clone_id))))
cat("[pdB 094D Clone3 rows]\n"); print(pdB[slice == "094D" & clone_id == "P14_Clone3"])
cat("[pdB 094D Clone4 rows]\n"); print(pdB[slice == "094D" & clone_id == "P14_Clone4"])
cat("[mB_nbr 094D Clone3]\n"); print(mB_nbr[slice == "094D" & clone_id == "P14_Clone3" & k == 6])
cat("[pdB 094D per zone x clone fill/estimable]\n")
print(pdB[slice == "094D", .(n = .N, n_na_fill = sum(is.na(fill)), n_est = sum(estimable)), by = .(clone_id, region)][order(clone_id, region)])
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# HTML 内嵌图刷新（094 面板 v4 后）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
html_f <- file.path(OUT, "TAM_clone_汇总.html")
txt <- paste(readLines(html_f, warn = FALSE), collapse = "\n")
stems <- c("panel_094_methodA","panel_094_methodB","heatmap_methodA_094C","heatmap_methodB_094C")
locs <- regmatches(txt, gregexpr("data:image/png;base64,[A-Za-z0-9+/=]+", txt))[[1]]
stopifnot(length(locs) == length(stems))
for (i in seq_along(stems)) {
  new <- paste0("data:image/png;base64,", base64enc::base64encode(file.path(FIG, paste0(stems[i], ".png"))))
  txt <- sub(locs[i], new, txt, fixed = TRUE)
}
writeLines(txt, html_f, useBytes = TRUE)
cat("[refreshed] size:", round(file.info(html_f)$size/1024,1), "KB | html mtime:",
    format(file.info(html_f)$mtime, "%H:%M:%S"), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 最终 doc-disk 对账（HTML/图更新后）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
script <- "/path/to/.claude/scripts/doc_disk_reconcile.py"
out <- suppressWarnings(system2("python3", c(script, "--all", OUT), stdout = TRUE, stderr = TRUE))
cat(paste(out, collapse = "\n"), "\n")
cat("[exit]", attr(out, "status"), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 审计后文字修正：区分“符号一致”与“方向标签一致”、修正 094 clone 计数、补充不可估格置空说明；并刷新内嵌图
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
html_f <- file.path(OUT, "TAM_clone_汇总.html")
txt <- paste(readLines(html_f, warn = FALSE), collapse = "\n")
must_replace <- function(x, old, new) {
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1) stop(sprintf("expected 1 occurrence, got %d: %s", n, old))
  sub(old, new, x, fixed = TRUE)
}
txt <- must_replace(txt, "094C/094D（P14，5–6 个 clone）", "094C/094D（P14，各 5 个 clone）")
txt <- must_replace(txt, "094D 有 4 个 clone 在该区低于门槛，如实标注不可估。",
                    "094D 的 Clone1（n=12）与 Clone3（n=3）在方法A 整片低于门槛（共 4 个 clone×TAM 格）；方法B k=6 有 5 个 clone 的 7 个格不可估。")
txt <- must_replace(txt, "两方法方向高度一致：",
                    "两方法<b>符号</b>高度一致（按方向标签则较宽松，见 3.4）：")
txt <- gsub("同 spot 38/39（97.4%）",
            "同 spot 38/39（97.4%）；按方向标签：k=6 30/39、同 spot 29/39",
            txt, fixed = TRUE)
txt <- must_replace(txt, "方法A vs 方法B(k=6)：可比 39 对，方向一致 <b>39 对（100.0%）</b>。",
                    "方法A vs 方法B(k=6)：可比 39 对，<b>符号一致 39 对（100.0%）</b>；按方向标签（enriched/depleted/not significant）一致 30 对。")
txt <- must_replace(txt, "方法A vs 方法B(同 spot)：可比 39 对，方向一致 38 对（97.4%）。",
                    "方法A vs 方法B(同 spot)：可比 39 对，符号一致 38 对（97.4%）；按方向标签一致 29 对。唯一符号不一致为 094C P14_Clone4 C1Q（方法B Ro/e=0、方法A=0.179；两者方向标签均为 depleted，属 sign(0) 伪不一致）。")
txt <- must_replace(txt, "<li>clone 由 CNV 定义；clone 与肿瘤上皮比例存在同向混杂（此前分析已报告），本步未重新校正。</li>",
                    "<li>clone 由 CNV 定义；clone 与肿瘤上皮比例存在同向混杂（此前分析已报告），本步未重新校正。</li><li><b>不可估格的处理</b>：表中保留 n_spots_clone / n_pos / reason，点估计（Ro/e、obs）与 z/p 一律置空（NA）；不得把 estimable=FALSE 的格当作阴性或阳性结果引用。</li>")
txt <- must_replace(txt, "<h2>⑤ 图表索引</h2>",
                    "<p class=\"note\">图件自检：22/22 PNG 像素校验通过；独立 visual-inspector 复核确认热图结构、连续发散色阶（log2 Ro/e, −2…2）、灰色 n/a 编码、轴标题与标签可读；094 面板中每个切片使用各自独立的 clone 轴。</p><h2>⑤ 图表索引</h2>")
stems <- c("panel_094_methodA","panel_094_methodB","heatmap_methodA_094C","heatmap_methodB_094C")
locs <- regmatches(txt, gregexpr("data:image/png;base64,[A-Za-z0-9+/=]+", txt))[[1]]
stopifnot(length(locs) == length(stems))
for (i in seq_along(stems)) {
  new <- paste0("data:image/png;base64,", base64enc::base64encode(file.path(FIG, paste0(stems[i], ".png"))))
  txt <- sub(locs[i], new, txt, fixed = TRUE)
}
writeLines(txt, html_f, useBytes = TRUE)

mn_f <- file.path(TAB, "methods_note.md")
mn <- readLines(mn_f, warn = FALSE)
mn <- gsub("- 所有 p 为**片内名义 p**",
           "- 不可估格的处理：表中保留 n_* 与 reason，点估计（Ro/e、obs）与 z/p 一律置空（NA）；不得作为结论依据\n- 方向一致性口径：**符号一致**（sign(Ro/e) 相同）与**方向标签一致**（enriched/depleted/not significant 完全相同）是两件事；本步两者都报告\n- 所有 p 为**片内名义 p**",
           mn, fixed = TRUE)
writeLines(mn, mn_f, useBytes = TRUE)
cat("[HTML] size:", round(file.info(html_f)$size/1024,1), "KB | images:", length(stems), "\n")
cat("[methods_note] lines:", length(mn), "\n")
cat("[check] 5-6 clone remains:", grepl("5–6 个 clone", txt, fixed=TRUE),
    "| 4 个 clone remains:", grepl("094D 有 4 个 clone", txt, fixed=TRUE),
    "| sign vs dir text:", grepl("符号一致", txt, fixed=TRUE), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# methods_note 末条口径修正：符号一致 vs 方向标签一致
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
mn_f <- file.path(TAB, "methods_note.md")
mn <- readLines(mn_f, warn = FALSE)
old <- "- 方法A 与 方法B 方向一致率：k=6 邻域 39/39 对（100%）；同 spot 38/39 对（97.4%）；显著方向一致 27/27"
new <- "- 方法A 与 方法B：**符号一致率** k=6 39/39（100%）、同 spot 38/39（97.4%）；**方向标签一致率** k=6 30/39、同 spot 29/39（唯一符号不一致为 094C P14_Clone4 C1Q 的 sign(0)）；双方显著的方向一致 27/27"
stopifnot(sum(grepl(old, mn, fixed = TRUE)) == 1)
mn <- sub(old, new, mn, fixed = TRUE)
writeLines(mn, mn_f, useBytes = TRUE)
cat("[methods_note] patched. lines:", length(mn), "\n")
cat(mn[grep("符号一致率", mn)], "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 终检（面板顺序修复后）：不可估置空断言、HTML 合规、doc-disk 对账、094 面板 x 顺序
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
cat("[blank assertions] A radius/knn6/strat, B same/nbr non-est with values:",
    mA_radius[estimable==FALSE & !is.na(enrichment), .N],
    mA_knn6[estimable==FALSE & !is.na(enrichment), .N],
    mA_strat[estimable==FALSE & !is.na(enrichment), .N],
    mB_same[estimable==FALSE & !is.na(roe), .N],
    mB_nbr[estimable==FALSE & !is.na(roe), .N], "\n")
cat("[094 panel x order identical across zones] A:", uniqueN(R[plot=="p94A" & grepl("^094D", panel)]$order) == 1,
    "| B:", uniqueN(R[plot=="p94B" & grepl("^094D", panel)]$order) == 1, "\n")
h <- paste(readLines(file.path(OUT, "TAM_clone_汇总.html"), warn = FALSE), collapse = "\n")
h_txt <- gsub("data:image/png;base64,[A-Za-z0-9+/=]+", "[IMG]", h)
cat("[html] review/:", grepl("review/", h_txt, fixed=TRUE),
    "| 300/slice:", grepl("300/slice", h_txt, fixed=TRUE),
    "| 降采样:", grepl("降采样", h_txt, fixed=TRUE),
    "| placeholder:", grepl("placeholder|TODO", h_txt, ignore.case=TRUE),
    "| 符号一致:", grepl("符号一致", h_txt, fixed=TRUE),
    "| images:", lengths(regmatches(h, gregexpr("data:image/png;base64,", h, fixed=TRUE))), "\n")
script <- "/path/to/.claude/scripts/doc_disk_reconcile.py"
out <- suppressWarnings(system2("python3", c(script, "--all", OUT), stdout = TRUE, stderr = TRUE))
cat("[doc_disk_reconcile]", paste(out, collapse = " "), "\n")
csvs <- list.files(TAB, pattern = "\\.csv$", full.names = TRUE)
cat("[csv] any Inf:", any(sapply(csvs, function(f) any(grepl("Inf", readLines(f, warn=FALSE), fixed=TRUE)))),
    "| any NaN:", any(sapply(csvs, function(f) any(grepl("NaN", readLines(f, warn=FALSE), fixed=TRUE)))), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 内部溯源表补充：方向标签一致率与不可估格置空
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
tr_f <- file.path(REV, "TRACEABILITY.md")
tr <- readLines(tr_f, warn = FALSE)
add <- c(
"| 方向标签一致 A vs B k6 | 30/39 | tables/clone_tam_direction_comparison.csv | dir_agree_A_B6 |",
"| 方向标签一致 A vs B 同spot | 29/39 | tables/clone_tam_direction_comparison.csv | dir_agree_A_Bs |",
"| 不可估格置空 | 0 行仍有估计值 | tables/clone_tam_*csv | estimable=FALSE -> roe/enrichment/z/p = NA |",
"",
"## 口径说明",
"- 符号一致 = sign(Ro/e) 相同；方向标签一致 = enriched/depleted/not significant 完全相同。",
"- 唯一符号不一致：094C P14_Clone4 C1Q（方法A Ro/e=0.179，方法B 同spot Ro/e=0；两者方向标签均为 depleted，属 sign(0) 伪不一致）。"
)
writeLines(c(tr, add), tr_f, useBytes = TRUE)
cat("[TRACEABILITY] lines:", length(tr) + length(add), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 核实 094 面板各 panel 的实际 x 轴 clone 顺序（排除标签误读）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
b <- ggplot_build(p94B)
pp <- b$layout$panel_params
lab <- lapply(pp, function(z) z$x$get_labels())
names(lab) <- paste0("panel", seq_along(lab))
# facet_wrap 顺序 = panel_id levels（slice-major, zones = All/Core/Interface/Stroma）
pid_levels <- levels(pdB$panel_id)
if (is.null(pid_levels)) pid_levels <- levels(pdB[slice %in% S094]$panel_id)
cat("[094D panels]\n")
for (i in grep("^094D", pid_levels)) {
  cat(sprintf("  %-24s -> %s\n", gsub("\n", " | ", pid_levels[i]), paste(lab[[i]], collapse = ", ")))
}
cat("[094C panels]\n")
for (i in grep("^094C", pid_levels)) {
  cat(sprintf("  %-24s -> %s\n", gsub("\n", " | ", pid_levels[i]), paste(lab[[i]], collapse = ", ")))
}
cat("[all 094D x order identical]:",
    length(unique(sapply(grep("^094D", pid_levels), function(i) paste(lab[[i]], collapse=",")))) == 1, "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 核实 geom_blank 修复后 094 面板各 panel 的 x 轴 clone 顺序
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
pid_levels <- do.call(c, lapply(S094, function(s) paste(s, ZONE_LEVELS, sep = "\n")))
res <- list()
for (obj in c("p94A","p94B")) {
  b <- ggplot_build(get(obj)); pp <- b$layout$panel_params
  lab <- lapply(pp, function(z) z$x$get_labels())
  stopifnot(length(lab) == length(pid_levels))
  for (i in seq_along(pid_levels)) {
    res[[length(res)+1]] <- data.table(plot = obj, panel = gsub("\n", " | ", pid_levels[i]),
                                       order = paste(lab[[i]], collapse = ","))
  }
}
R <- rbindlist(res)
for (obj in c("p94A","p94B")) {
  cat(sprintf("[%s]\n", obj))
  for (sd in c("094C","094D")) {
    idx <- grep(paste0("^", sd, " \\|"), R[plot == obj]$panel)
    cat(sprintf("  %s identical across 4 zones: %s\n", sd,
        uniqueN(R[plot == obj][idx]$order) == 1))
    print(R[plot == obj][idx][, .(panel, order)])
  }
}
save_obj(R, file.path(DAT, "panel_x_order_check"))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 最终交付物清单（绝对路径 + 大小）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
files <- c(
  file.path(TAB, c("coloc_input_summary.csv","clone_tam_radius_perm_by_slice.csv",
                   "clone_tam_knn6_by_slice.csv","clone_tam_radius_perm_region_stratified.csv",
                   "clone_tam_coloc_samespot.csv","clone_tam_coloc_neighbor.csv",
                   "coloc_summary_by_slice.csv","clone_tam_direction_comparison.csv",
                   "clone_tam_direction_by_slice.csv","methods_note.md")),
  file.path(FIG, c("panel_094_methodA.png","panel_094_methodB.png")),
  file.path(OUT, "TAM_clone_汇总.html"),
  file.path(OUT, "notebooks/step01_tam_clone_coloc.ipynb"))
info <- data.table(file = files, exists = file.exists(files),
                   size_kb = round(file.info(files)$size/1024, 1))
print(info)
cat("[figures] png:", length(list.files(FIG, pattern="\\.png$")),
    "| pdf:", length(list.files(FIG, pattern="\\.pdf$")), "\n")
cat("[notebook] cells:", 32, "(1 markdown + 31 code; cell 0 empty)\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 修复 §3.3 分层可估表的整数格式（原 format="g" 把 16/24、40/43 都显示成 2e+01/4e+01）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
html_f <- file.path(OUT, "TAM_clone_汇总.html")
txt <- paste(readLines(html_f, warn = FALSE), collapse = "\n")
int_table <- function(dt) {
  d <- as.data.frame(dt)
  for (j in names(d)) d[[j]] <- if (is.numeric(d[[j]])) format(as.integer(d[[j]]), scientific = FALSE, trim = TRUE) else d[[j]]
  hdr <- paste0("<tr>", paste0("<th>", names(d), "</th>", collapse = ""), "</tr>")
  body <- apply(d, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<table>", hdr, paste0(body, collapse = ""), "</table>")
}
tblA <- int_table(strat_A); tblB <- int_table(strat_B)
start_tag <- "<h3>3.3 三区分层可估性</h3>"
end_tag <- "<div class=\"note\">Stroma 在 10 片中"
seg_new <- paste0(start_tag, tblA,
  "<p style=\"font-size:12px;color:#555\">表：方法A 分层（每区 54 格 = 10 片 × 该片 clone × 2 TAM）。</p>",
  tblB, "<p style=\"font-size:12px;color:#555\">表：方法B k=6 分层。</p>")
pat <- paste0(start_tag, "[\\s\\S]*?", end_tag)
stopifnot(lengths(regmatches(txt, gregexpr(pat, txt, perl = TRUE))) == 1)
txt <- sub(pat, paste0(seg_new, end_tag), txt, perl = TRUE)
writeLines(txt, html_f, useBytes = TRUE)
cat("[HTML] size:", round(file.info(html_f)$size/1024,1), "KB\n")
h_txt <- gsub("data:image/png;base64,[A-Za-z0-9+/=]+", "[IMG]", txt)
cat("[check] 2e+01 present:", grepl("2e+01", h_txt, fixed=TRUE),
    "| 4e+01 present:", grepl("4e+01", h_txt, fixed=TRUE),
    "| 16 present:", grepl(">16<", h_txt, fixed=TRUE),
    "| 24 present:", grepl(">24<", h_txt, fixed=TRUE),
    "| 40 present:", grepl(">40<", h_txt, fixed=TRUE),
    "| 43 present:", grepl(">43<", h_txt, fixed=TRUE), "\n")
cat(tblA, "\n", tblB, "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 终局：doc-disk 对账、追加 audit_done、汇总最终数字供 state 写入
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
script <- "/path/to/.claude/scripts/doc_disk_reconcile.py"
out <- suppressWarnings(system2("python3", c(script, "--all", OUT), stdout = TRUE, stderr = TRUE))
cat("[doc_disk_reconcile]", paste(out, collapse = " "), "\n")

# audit_done（critic 复审 pass 后由 doer 追加；东八区时间）
log_dir <- "/path/to/.claude/logs/2026-10-26"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
audit_line <- sprintf('{"ts":"%s","session":"c5e99107-e355-4979-a9aa-f5cb722d49d9","subagent_type":"sc-analyst","verdict":"pass"}',
                      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
cat(audit_line, "\n", sep = "", file = file.path(log_dir, "audit_done.jsonl"), append = TRUE)
cat("[audit_done] appended. tail:\n")
cat(tail(readLines(file.path(log_dir, "audit_done.jsonl"), warn = FALSE), 2), sep = "\n")

cat("\n[final numbers for state]\n")
cat("n_spots=", n_spots, " n_slices=", n_slices, " n_clones_unique=", n_clones_u,
    " n_slice_clone=", n_slice_clone, "\n")
cat("A whole: enr=", a_enr, " dep=", a_dep, " ns=", a_ns, " notest=", sum(!mA_radius$estimable), "\n")
cat("B k6: enr=", b_enr, " dep=", b_dep, " ns=", b6_ns, " notest=", sum(!mB_nbr[zone=="All" & k==K_MAIN, estimable]), "\n")
cat("B same: enr=", bs_enr, " dep=", bs_dep, " ns=", nrow(bs)-bs_enr-bs_dep, "\n")
cat("strat A:", paste(strat_A$region, strat_A$estimable, sep="=", collapse="; "), "\n")
cat("strat B:", paste(strat_B$zone, strat_B$estimable, sep="=", collapse="; "), "\n")
cat("sign k6=", cmp_stats$agree_A_B6, "/", cmp_stats$comparable_A_B6,
    " sign same=", cmp_stats$agree_A_Bs, "/", cmp_stats$comparable_A_Bs,
    " dir k6=", cmp_stats$dir_agree_A_B6, " dir same=", cmp_stats$dir_agree_A_Bs,
    " sig=", cmp_stats$same_dir_sig, "/", cmp_stats$both_sig_A_B6, "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 诊断：列出 HTML 可见文本中所有含 Stroma / 三区 / 4×4 的上下文，便于精确替换
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
h <- paste(readLines(file.path(OUT, "TAM_clone_汇总.html"), warn = FALSE), collapse = "\n")
h_txt <- gsub("data:image/png;base64,[A-Za-z0-9+/=]+", "[IMG]", h)
m <- gregexpr("Stroma|三区|4×4|4 x 4", h_txt)
locs <- regmatches(h_txt, m)[[1]]
cat("[hits]", length(locs), "\n")
# 打印每处前后 90 字符
pos <- m[[1]]
for (i in seq_along(pos)) {
  a <- max(1, pos[i] - 90); b <- min(nchar(h_txt), pos[i] + 90)
  cat(sprintf("%d) ...%s...\n", i, substr(h_txt, a, b)))
}
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 文档同步：删除 Stroma 章节/列，分层只呈现 Core/Interface；刷新 4 张内嵌图
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
html_f <- file.path(OUT, "TAM_clone_汇总.html")
txt <- paste(readLines(html_f, warn = FALSE), collapse = "\n")
rep1 <- function(x, old, new) {
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1) stop(sprintf("expected 1, got %d: %s", n, old))
  sub(old, new, x, fixed = TRUE)
}
int_table <- function(dt) {
  d <- as.data.frame(dt)
  for (j in names(d)) d[[j]] <- if (is.numeric(d[[j]])) format(as.integer(d[[j]]), scientific = FALSE, trim = TRUE) else d[[j]]
  hdr <- paste0("<tr>", paste0("<th>", names(d), "</th>", collapse = ""), "</tr>")
  body <- apply(d, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<table>", hdr, paste0(body, collapse = ""), "</table>")
}
# §3.1 输入表：只保留 Core/Interface
s31 <- "<h3>3.1 输入与关联率（逐片 100%）</h3>"
new31 <- paste0(s31, int_table(input_summary),
  '<p style="font-size:12px;color:#555">n_universe = 该片全部组织 spot；分层仅使用 Tumor Core / Interface（Stroma 已移除）。</p>')
txt <- sub(paste0(s31, "[\\s\\S]*?</table>"), new31, txt, perl = TRUE)
# §3.3 分层表与说明：只 Core/Interface
s33 <- "<h3>3.3 三区分层可估性</h3>"
new33 <- paste0('<h3>3.3 分层可估性（Tumor Core / Interface）</h3>',
  int_table(strat_A), '<p style="font-size:12px;color:#555">表：方法A 分层（每区 54 格 = 10 片 × 该片 clone × 2 TAM）。</p>',
  int_table(strat_B), '<p style="font-size:12px;color:#555">表：方法B k=6 分层。</p>',
  '<div class="note">Stroma 因每片 clone spot ≤1、不可估，已从分析范围移除；分层只保留 Tumor Core 与 Interface。</div>')
txt <- sub(paste0(s33, "[\\s\\S]*?</div>"), new33, txt, perl = TRUE)
# 逐条文字修正
txt <- rep1(txt, "三区分层中 Stroma 全部不可估（每片 clone spot 0–1 个）。",
            "Stroma 因每片 clone spot 0–1 个、不可估，已从分层范围移除（分层只保留 Tumor Core 与 Interface）。")
txt <- rep1(txt, "三区来自 project_singlecell V2 底表的 <code>region</code> 列",
            "区域来自 project_singlecell V2 底表的 <code>region</code> 列（分层仅使用 Tumor Core / Interface）")
txt <- rep1(txt, "094C（P14）方法A：All + 三区分面；每格数字为 Ro/e，* 为 p<0.05",
            "094C（P14）方法A：All + 两区（Core/Interface）；每格数字为 Ro/e，* 为 p<0.05")
txt <- rep1(txt, "094C（P14）方法B k=6：All + 三区分面",
            "094C（P14）方法B k=6：All + 两区（Core/Interface）")
txt <- rep1(txt, "094 系列：方法A 半径邻域，log2(Ro/e)；n/a = 不可估（按切片分面，clone 带病人前缀）",
            "094 系列：方法A 半径邻域（仅 Tumor Core / Interface），log2(Ro/e)；n/a = 不可估（按切片分面，clone 带病人前缀）")
txt <- rep1(txt, "094 系列：方法B k=6 邻域，log2(Ro/e)；n/a = 不可估",
            "094 系列：方法B k=6 邻域（仅 Tumor Core / Interface），log2(Ro/e)；n/a = 不可估")
txt <- rep1(txt, "<li>小 clone（n=3/12/22）在部分层不可估；Stroma 全片不可估；这些格在表中保留并写明原因，不填 0、不画空面板。</li>",
            "<li>小 clone（n=3/12/22）在部分层不可估；这些格在表中保留并写明原因，不填 0、不画空面板。</li>")
txt <- rep1(txt, "切片 × 病人 × clone 数 × spot 数 × 三区 spot 数 × 关联率",
            "切片 × 病人 × clone 数 × spot 数 × Core/Interface spot 数 × 关联率")
txt <- rep1(txt, "方法A 三区分层（含 n_spots 与不可估原因）", "方法A 分层（Core/Interface；含 n_spots 与不可估原因）")
txt <- rep1(txt, "STEP 01 — clone × TAM 共定位（逐切片、三区分层）", "STEP 01 — clone × TAM 共定位（逐切片、Core/Interface 分层）")
txt <- rep1(txt, "<li>这种关系在 Tumor Core / Interface / Stroma 内是否依然存在？</li>",
            "<li>这种关系在 Tumor Core / Interface 内是否依然存在？</li>")
txt <- rep1(txt, "逐片与 project_singlecell 三区表、step3 TAM 比例按 key 合并", "逐片与 project_singlecell 区域表、step3 TAM 比例按 key 合并")
txt <- rep1(txt, "在 All/三区 pool 内置换 clone 标签", "在 All/Core/Interface pool 内置换 clone 标签")
txt <- rep1(txt, "zone（All/三区）", "zone（All/Core/Interface）")
txt <- rep1(txt, "分层可估：方法A Tumor Core 40 / Interface 16 / Stroma 0；方法B k=6 Tumor Core 43 / Interface 24 / Stroma 0。",
            "分层可估（Stroma 已移除）：方法A Tumor Core 40/54、Interface 16/54；方法B k=6 Tumor Core 43/54、Interface 24/54。")
# 刷新 4 张内嵌图
stems <- c("panel_094_methodA","panel_094_methodB","heatmap_methodA_094C","heatmap_methodB_094C")
locs <- regmatches(txt, gregexpr("data:image/png;base64,[A-Za-z0-9+/=]+", txt))[[1]]
stopifnot(length(locs) == length(stems))
for (i in seq_along(stems)) {
  new <- paste0("data:image/png;base64,", base64enc::base64encode(file.path(FIG, paste0(stems[i], ".png"))))
  txt <- sub(locs[i], new, txt, fixed = TRUE)
}
writeLines(txt, html_f, useBytes = TRUE)
h_txt <- gsub("data:image/png;base64,[A-Za-z0-9+/=]+", "[IMG]", txt)
cat("[html] size:", round(file.info(html_f)$size/1024,1), "KB | Stroma in visible text:",
    grepl("Stroma", h_txt, fixed=TRUE), "| 三区:", grepl("三区", h_txt, fixed=TRUE), "\n")

# methods_note：移除 Stroma 层，并写明原因
mn_f <- file.path(TAB, "methods_note.md"); mn <- readLines(mn_f, warn = FALSE)
mn <- gsub("All / Tumor Core / Interface / Stroma", "All / Tumor Core / Interface", mn, fixed = TRUE)
mn <- gsub("- Stroma: 每片 clone spot 仅 0–1 个（093C 为 1），region_n<100，全部不可估；不是“填 0”",
           "- Stroma 因每片 clone spot ≤1、不可估，已从分析范围移除（不属于本步分层结果；分层仅 Tumor Core 与 Interface）",
           mn, fixed = TRUE)
writeLines(mn, mn_f, useBytes = TRUE)
cat("[methods_note] Stroma lines:\n"); cat(grep("Stroma", mn, value = TRUE), sep = "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 定位 HTML 可见文本中残留的 Stroma
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
h <- paste(readLines(file.path(OUT, "TAM_clone_汇总.html"), warn = FALSE), collapse = "\n")
h_txt <- gsub("data:image/png;base64,[A-Za-z0-9+/=]+", "[IMG]", h)
pos <- gregexpr("Stroma", h_txt)[[1]]
cat("[remaining hits]", length(pos), "\n")
for (p in pos) cat(sprintf("...%s...\n", substr(h_txt, max(1,p-100), min(nchar(h_txt), p+100))))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 用当前 2 区分层对象重建 §3.3 表；修正 methods_note 的 region 列说明
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
html_f <- file.path(OUT, "TAM_clone_汇总.html")
txt <- paste(readLines(html_f, warn = FALSE), collapse = "\n")
int_table <- function(dt) {
  d <- as.data.frame(dt)
  for (j in names(d)) d[[j]] <- if (is.numeric(d[[j]])) format(as.integer(d[[j]]), scientific = FALSE, trim = TRUE) else d[[j]]
  hdr <- paste0("<tr>", paste0("<th>", names(d), "</th>", collapse = ""), "</tr>")
  body <- apply(d, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<table>", hdr, paste0(body, collapse = ""), "</table>")
}
strat_A <- mA_strat[, .(estimable = sum(estimable), rows = .N), by = region]
strat_B <- mB_nbr[k == K_MAIN & zone != "All", .(estimable = sum(estimable), rows = .N), by = zone]
cat("[strat_A]\n"); print(strat_A); cat("[strat_B]\n"); print(strat_B)
s33 <- "<h3>3.3 分层可估性（Tumor Core / Interface）</h3>"
new33 <- paste0(s33,
  int_table(strat_A), '<p style="font-size:12px;color:#555">表：方法A 分层（每区 54 格 = 10 片 × 该片 clone × 2 TAM）。</p>',
  int_table(strat_B), '<p style="font-size:12px;color:#555">表：方法B k=6 分层。</p>',
  '<div class="note">Stroma 因每片 clone spot ≤1、不可估，已从分析范围移除；分层只保留 Tumor Core 与 Interface。</div>')
stopifnot(lengths(regmatches(txt, gregexpr(paste0(s33, "[\\s\\S]*?</div>"), txt, perl = TRUE))) == 1)
txt <- sub(paste0(s33, "[\\s\\S]*?</div>"), new33, txt, perl = TRUE)
writeLines(txt, html_f, useBytes = TRUE)
h_txt <- gsub("data:image/png;base64,[A-Za-z0-9+/=]+", "[IMG]", txt)
cat("[html] Stroma rows:", grepl("<td>Stroma</td>", h_txt, fixed=TRUE),
    "| 三区:", grepl("三区", h_txt, fixed=TRUE), "| size_kb:", round(file.info(html_f)$size/1024,1), "\n")
# methods_note：说明 region 列含 Stroma 但本步只用前两个
mn_f <- file.path(TAB, "methods_note.md"); mn <- readLines(mn_f, warn = FALSE)
mn <- gsub("region 列 ∈ {Tumor Core, Interface, Stroma}",
           "region 列 ∈ {Tumor Core, Interface, Stroma}（本步仅使用 Tumor Core / Interface）", mn, fixed = TRUE)
writeLines(mn, mn_f, useBytes = TRUE)
cat("[methods_note] lines:", length(mn), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 终检（2 区版）：表内无 Stroma、关键计数、图尺寸、doc-disk 对账
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
cat("[CSV zones]\n")
cat("  A strat regions:", paste(sort(unique(mA_strat$region)), collapse=" / "),
    "| rows:", nrow(mA_strat), "| estimable:", sum(mA_strat$estimable), "\n")
cat("  B same zones:", paste(sort(unique(mB_same$zone)), collapse=" / "), "| rows:", nrow(mB_same), "\n")
cat("  B nbr zones:", paste(sort(unique(mB_nbr$zone)), collapse=" / "), "| rows:", nrow(mB_nbr), "\n")
cat("  input cols:", paste(names(input_summary), collapse=","), "\n")
cat("[estimable] A Core/Interface:", paste(mA_strat[, .(n=sum(estimable)), by=region]$region,
    mA_strat[, .(n=sum(estimable)), by=region]$n, sep="=", collapse="; "), "\n")
cat("[estimable] B k6 Core/Interface:", paste(mB_nbr[zone!="All", .(n=sum(estimable)), by=zone]$zone,
    mB_nbr[zone!="All", .(n=sum(estimable)), by=zone]$n, sep="=", collapse="; "), "\n")
cat("[figure dimensions]\n")
for (stem in c("heatmap_methodA_092A","heatmap_methodA_094C","heatmap_methodB_092A","heatmap_methodB_094C",
               "panel_094_methodA","panel_094_methodB")) {
  ii <- magick::image_info(magick::image_read(file.path(FIG, paste0(stem, ".png"))))
  cat(sprintf("  %-24s %d x %d px  (%.1f x %.1f in @600dpi)\n", stem, ii$width, ii$height, ii$width/600, ii$height/600))
}
h <- paste(readLines(file.path(OUT, "TAM_clone_汇总.html"), warn=FALSE), collapse="\n")
h_txt <- gsub("data:image/png;base64,[A-Za-z0-9+/=]+", "[IMG]", h)
cat("[html] Stroma rows:", grepl("<td>Stroma</td>", h_txt, fixed=TRUE),
    "| 三区:", grepl("三区", h_txt, fixed=TRUE),
    "| 40/54:", grepl("40/54", h_txt, fixed=TRUE), "| 43/54:", grepl("43/54", h_txt, fixed=TRUE),
    "| review/:", grepl("review/", h_txt, fixed=TRUE), "\n")
script <- "/path/to/.claude/scripts/doc_disk_reconcile.py"
out <- suppressWarnings(system2("python3", c(script, "--all", OUT), stdout=TRUE, stderr=TRUE))
cat("[doc_disk_reconcile]", paste(out, collapse=" "), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 内部溯源表更新：Stroma 移除、2 区分层计数、图尺寸
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
tr_f <- file.path(REV, "TRACEABILITY.md")
tr <- readLines(tr_f, warn = FALSE)
tr <- gsub("方法A Tumor Core 40 / Interface 16 / Stroma 0", "方法A Tumor Core 40 / Interface 16", tr, fixed = TRUE)
tr <- gsub("方法B k=6 Tumor Core 43 / Interface 24 / Stroma 0", "方法B k=6 Tumor Core 43 / Interface 24", tr, fixed = TRUE)
tr <- c(tr, "",
"## 2026-10-26 Stroma 移除（用户要求）",
"- 分层范围改为 Tumor Core / Interface；Stroma 因每片 clone spot ≤1、不可估而移除。",
"- 分层表行数：方法A 108（2 区 × 54）；方法B same 162（All+2 区 × 54）、neighbor 324（× k=6/18）。",
"- 可估格：方法A Core 40/54、Interface 16/54；方法B k=6 Core 43/54、Interface 24/54。All 层结果不变。",
"- 图尺寸：每片热图 5.1×3.8 in（3 clone）/ 6.9×3.8 in（5 clone）/ 4.2×3.8 in（2 clone）；094 面板 7.2×7.8 in（4 片 × 2 区 = 8 panel）。",
"- 像素自检 22/22；doc_disk_reconcile 0。")
writeLines(tr, tr_f, useBytes = TRUE)
cat("[TRACEABILITY] lines:", length(tr), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 读取四张共定位结果表并整理为绘图长表；同时核对 P13 旧值（本 cell 只做数据准备，不画图）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(data.table); library(ggplot2)})
for (pk in c("patchwork", "scales", "grid", "png")) cat("[pkg]", pk, ":", requireNamespace(pk, quietly = TRUE), "\n")

ROOT_PP <- "/path/to/project/project_spatial/TNBC/patient_ALL"
OUT_PP  <- file.path(ROOT_PP, "TAM_clone")
TAB_PP  <- file.path(OUT_PP, "tables")
FIG_PP  <- file.path(OUT_PP, "figures", "per_patient")
dir.create(FIG_PP, recursive = TRUE, showWarnings = FALSE)

A_reg <- fread(file.path(TAB_PP, "clone_tam_radius_perm_region_stratified.csv"))
A_sl  <- fread(file.path(TAB_PP, "clone_tam_radius_perm_by_slice.csv"))
B_all <- fread(file.path(TAB_PP, "clone_tam_coloc_neighbor.csv"))

# 统一：effect 值（A=enrichment / B=roe）、方向匹配单侧 p、可估标志
# 方法A：富集方向用 p_perm，耗竭方向用 p_deplete（p_perm 是单侧富集 p，直接用于耗竭柱会误判）
A_reg[, `:=`(zone = region, eff = enrichment,
             p_dir = fifelse(enrichment >= 1, p_perm, p_deplete))]
A_sl[,  `:=`(zone = "All", eff = enrichment,
             p_dir = fifelse(enrichment >= 1, p_perm, p_deplete))]
B_use <- B_all[k == 6]
B_use[, `:=`(eff = roe, p_dir = pmin(p_enrich, p_deplete))]

cat("\n[dim] A_reg", nrow(A_reg), "x", ncol(A_reg),
    "| A_sl", nrow(A_sl), "x", ncol(A_sl),
    "| B(k=6)", nrow(B_use), "x", ncol(B_use), "\n")
cat("[patients] ", paste(sort(unique(A_reg$patient)), collapse = ", "), "\n", sep = "")
cat("[zones] A: ", paste(unique(A_reg$zone), collapse = " / "),
    " | B: ", paste(unique(B_use$zone), collapse = " / "), "\n", sep = "")
cat("[estimable] A_reg: ", paste(A_reg[, .N, by = .(zone, estimable)][, paste0(zone, "=", estimable, ":", N)], collapse = "; "),
    "\n            B: ", paste(B_use[, .N, by = .(zone, estimable)][, paste0(zone, "=", estimable, ":", N)], collapse = "; "), "\n", sep = "")

# 每个病人 × 方法 × 区域 × TAM 的条目数（决定 figC(a) 合并口径）
cat("\n[A_reg entries per patient/zone/tam]\n")
print(A_reg[estimable == TRUE, .N, by = .(patient, zone, tam)][order(patient, zone, tam)])
cat("\n[B entries per patient/zone/tam]\n")
print(B_use[estimable == TRUE, .N, by = .(patient, zone, tam)][order(patient, zone, tam)])

# ---- P13 旧值（step04 audit D 行；已核对与 step6 归档表数值一致）----
OLD_F <- "/path/to/project/project_spatial/TNBC/TNBC_analysis/SpatialInferCNV/results/step6_project_singlecell_region_coloc/tables/clone_tam_coloc_neighbor.csv"
old_dt <- fread(OLD_F)[k == 6 & zone == "All" & slice %in% c("094A", "094B")]
old_dt <- old_dt[, .(slice, clone, tam, old_roe = roe)]
new_dt <- B_use[patient == "P13" & zone == "All", .(slice, clone, tam, new_roe = roe)]
chk <- merge(old_dt, new_dt, by = c("slice", "clone", "tam"))
cat("\n[P13 old = audit D row] 094A Clone1 C1Q", round(old_dt[slice=="094A" & clone=="Clone1" & tam=="C1Q TAMs", old_roe], 3),
    "| 094A Clone2 C1Q", round(old_dt[slice=="094A" & clone=="Clone2" & tam=="C1Q TAMs", old_roe], 3),
    "| 094B Clone1 C1Q", round(old_dt[slice=="094B" & clone=="Clone1" & tam=="C1Q TAMs", old_roe], 3),
    "| 094B Clone2 C1Q", round(old_dt[slice=="094B" & clone=="Clone2" & tam=="C1Q TAMs", old_roe], 3), "\n")
stopifnot(nrow(chk) == 8)
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# P13 新旧一致性：原标签直接对比 + Clone1/Clone2 互换后对比，输出 CSV 与 2x2 分面对比图
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
FIG_PP <- file.path("/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone", "figures", "per_patient")
save_pp <- function(p, stem, w, h) {
  ggsave(file.path(FIG_PP, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 600, bg = "white", device = "png")
  ggsave(file.path(FIG_PP, paste0(stem, ".pdf")), p, width = w, height = h,
         bg = "white", device = grDevices::pdf)
  invisible(stem)
}

old_w <- dcast(old_dt, slice + tam ~ clone, value.var = "old_roe")
new_w <- dcast(new_dt, slice + tam ~ clone, value.var = "new_roe")
stopifnot(nrow(old_w) == 4, nrow(new_w) == 4)
cons <- merge(old_w, new_w, by = c("slice", "tam"))
setnames(cons, c("Clone1.x", "Clone2.x", "Clone1.y", "Clone2.y"),
               c("old_Clone1", "old_Clone2", "new_Clone1", "new_Clone2"))
cons[, tam := factor(tam, levels = c("C1Q TAMs", "SPP1 TAMs"))]
setorder(cons, slice, tam)

dir_of <- function(v) fifelse(v > 1, ">1", fifelse(v < 1, "<1", "=1"))
PRE <- "两次 clone 运行的标签在 600 个共同细胞上有 75.5% 支持互换（swapped 列按此前提对齐）"
SCALE <- "旧=3,600 spot（patient_cnv clone，每片 300）；新=9,678 spot（patient_ALL 全量 clone）；Ro/e 同定义但克隆构成不同，绝对值不可直接比较"

cons_long <- rbindlist(lapply(c("Clone1", "Clone2"), function(cl) {
  data.table(slice = cons$slice, tam = as.character(cons$tam), clone_label = cl,
             old_roe         = if (cl == "Clone1") cons$old_Clone1 else cons$old_Clone2,
             new_roe_same    = if (cl == "Clone1") cons$new_Clone1 else cons$new_Clone2,
             new_roe_swapped = if (cl == "Clone1") cons$new_Clone2 else cons$new_Clone1)
}))
cons_long[, `:=`(
  old_dir         = dir_of(old_roe),
  new_dir_same    = dir_of(new_roe_same),
  new_dir_swapped = dir_of(new_roe_swapped),
  direct_diff     = new_roe_same - old_roe,
  swapped_diff    = new_roe_swapped - old_roe
)]
cons_long[, direct_dir_agree  := old_dir == new_dir_same]
cons_long[, swapped_dir_agree := old_dir == new_dir_swapped]
cons_long[, `:=`(premise = PRE, scale_note = SCALE)]
setcolorder(cons_long, c("slice", "tam", "clone_label", "old_roe", "new_roe_same", "new_roe_swapped",
                         "old_dir", "new_dir_same", "new_dir_swapped",
                         "direct_diff", "swapped_diff", "direct_dir_agree", "swapped_dir_agree",
                         "premise", "scale_note"))
setorder(cons_long, slice, tam, clone_label)

CONS_CSV <- file.path(TAB_PP, "P13_old_vs_new_consistency.csv")
CONS_OUT <- copy(cons_long)
num_cols <- c("old_roe", "new_roe_same", "new_roe_swapped", "direct_diff", "swapped_diff")
CONS_OUT[, (num_cols) := lapply(.SD, round, 4), .SDcols = num_cols]
fwrite(CONS_OUT, CONS_CSV)
cat("[write]", CONS_CSV, "\n")
print(as.data.frame(cons_long[, .(slice, tam, clone_label,
    old_roe = round(old_roe, 3), new_same = round(new_roe_same, 3), new_swapped = round(new_roe_swapped, 3),
    old_dir, new_dir_same, new_dir_swapped, direct_dir_agree, swapped_dir_agree)]))

# ---- 汇总：方向一致率 + 秩相关 + 克隆偏好 ----
rho_direct  <- cor(cons_long$old_roe, cons_long$new_roe_same,    method = "spearman")
rho_swapped <- cor(cons_long$old_roe, cons_long$new_roe_swapped, method = "spearman")
cat(sprintf("\n[Summary] direction vs 1: direct %d/8 | swapped %d/8\n",
            sum(cons_long$direct_dir_agree), sum(cons_long$swapped_dir_agree)))
cat(sprintf("[Summary] Spearman(old, new): direct rho=%.2f | swapped rho=%.2f\n", rho_direct, rho_swapped))
cat(sprintf("[Summary] |log2(new/old)| after swap: median=%.2f, range %.2f-%.2f\n",
            median(abs(log2(cons_long$new_roe_swapped / cons_long$old_roe))),
            min(abs(log2(cons_long$new_roe_swapped / cons_long$old_roe))),
            max(abs(log2(cons_long$new_roe_swapped / cons_long$old_roe)))))

pref <- cons[, .(slice, tam = as.character(tam),
                 old_pref  = fifelse(old_Clone1 > old_Clone2, "Clone1", "Clone2"),
                 new_pref  = fifelse(new_Clone1 > new_Clone2, "Clone1", "Clone2"),
                 swap_pref = fifelse(new_Clone2 > new_Clone1, "Clone1", "Clone2"))]
pref[, `:=`(direct_pref_match = old_pref == new_pref, swapped_pref_match = old_pref == swap_pref)]
cat("\n[clone preference per slice x TAM] higher Ro/e clone\n")
print(as.data.frame(pref))
cat(sprintf("[Summary] preference match: direct %d/4 | swapped %d/4\n",
            sum(pref$direct_pref_match), sum(pref$swapped_pref_match)))
cat(sprintf("[Summary] C1Q only: direct %d/2 | swapped %d/2\n",
            sum(pref[tam == "C1Q TAMs"]$direct_pref_match), sum(pref[tam == "C1Q TAMs"]$swapped_pref_match)))
cat(sprintf("[Summary] SPP1 only: direct %d/2 | swapped %d/2\n",
            sum(pref[tam == "SPP1 TAMs"]$direct_pref_match), sum(pref[tam == "SPP1 TAMs"]$swapped_pref_match)))

# ---- 对比图：2x2 分面散点（行=对照方式，列=切片；颜色+形状= TAM×clone 合并图例）----
sc <- rbind(
  cons_long[, .(slice, tam, clone_label, old_roe, new_roe = new_roe_same,
                comparison = "direct labels (Clone1->Clone1)")],
  cons_long[, .(slice, tam, clone_label, old_roe, new_roe = new_roe_swapped,
                comparison = "labels swapped (Clone1<->Clone2)")]
)
sc[, comparison := factor(comparison, levels = c("direct labels (Clone1->Clone1)",
                                                 "labels swapped (Clone1<->Clone2)"))]
sc[, tam := factor(tam, levels = c("C1Q TAMs", "SPP1 TAMs"))]
sc[, grp := factor(paste0(ifelse(tam == "C1Q TAMs", "C1Q ", "SPP1 "),
                          sub("^Clone", "Clone", clone_label)),
                   levels = c("C1Q Clone1", "C1Q Clone2", "SPP1 Clone1", "SPP1 Clone2"))]
sc[, slice_f := factor(slice, levels = c("094A", "094B"))]
comp_lab <- c("direct labels (Clone1->Clone1)"   = sprintf("direct labels (rho = %.2f)", rho_direct),
              "labels swapped (Clone1<->Clone2)" = sprintf("labels swapped (rho = %.2f)", rho_swapped))
GRP_LV <- levels(sc$grp)
rng <- range(c(sc$old_roe, sc$new_roe), na.rm = TRUE)
lim <- c(max(0, rng[1] - 0.18), rng[2] + 0.18)

p_cons <- ggplot(sc, aes(x = old_roe, y = new_roe, colour = grp, shape = grp)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, colour = "grey35", linewidth = 0.4) +
  geom_point(size = 2.2, stroke = 0.55) +
  facet_grid(comparison ~ slice_f, labeller = labeller(comparison = comp_lab)) +
  scale_colour_manual(name = NULL, breaks = GRP_LV,
                      values = c("C1Q Clone1" = "#0072B2", "C1Q Clone2" = "#0072B2",
                                 "SPP1 Clone1" = "#D55E00", "SPP1 Clone2" = "#D55E00")) +
  scale_shape_manual(name = NULL, breaks = GRP_LV,
                     values = c("C1Q Clone1" = 16, "C1Q Clone2" = 17,
                                "SPP1 Clone1" = 16, "SPP1 Clone2" = 17)) +
  guides(colour = guide_legend(title = NULL, nrow = 2, override.aes = list(size = 2.2)),
         shape  = guide_legend(title = NULL, nrow = 2)) +
  coord_fixed(xlim = lim, ylim = lim) +
  labs(x = "Old Ro/e (3,600-spot clone run)", y = "New Ro/e (9,678-spot clone run)",
       caption = paste0("Grey solid = y = x; grey dashed = Ro/e 1.0.\n",
                        "Bottom row applies the 75.5% label-swap premise (600 shared cells).")) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
        strip.text = element_text(size = 7.3, face = "bold"),
        axis.text = element_text(size = 7, colour = "black"),
        axis.title = element_text(size = 7.5),
        legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(size = 7), legend.key.size = unit(3.0, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        plot.title = element_text(size = 9, hjust = 0.5),
        plot.caption = element_text(size = 7, hjust = 0, colour = "grey25"),
        plot.margin = margin(4, 4, 3, 3)) +
  ggtitle("P13: old vs new clone × TAM Ro/e (k=6, whole slice)")

save_pp(p_cons, "P13_old_vs_new", 6.8, 4.8)
cat("[figure] P13_old_vs_new ->", FIG_PP, "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 定义绘图辅助函数与区域级柱值聚合（Stouffer 合并 z，方向一致；A 病人内对 4 个预声明柱做 BH）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

TAM_COL <- c("C1Q TAMs" = "#0072B2", "SPP1 TAMs" = "#D55E00")
TAM_X   <- c("C1Q TAMs" = "C1QC",   "SPP1 TAMs" = "SPP1")
ZONE_LV <- c("Tumor Core", "Interface")
ZONE_X  <- c("Tumor Core" = "Core", "Interface" = "Interface")
star_p  <- function(p) fifelse(is.na(p), "ns",
                        fifelse(p < 0.001, "***", fifelse(p < 0.01, "**", fifelse(p < 0.05, "*", "ns"))))
star_t  <- function(p) fifelse(is.na(p) | p >= 0.05, "",
                        fifelse(p < 0.001, "***", fifelse(p < 0.01, "**", "*")))
fmt_v   <- function(v) fifelse(is.na(v), "n/a", formatC(v, format = "f", digits = 2))

# --- 条目级方向匹配 p，并按病人做 BH（仅方法A；方法B 按任务口径不跨片校正）---
A_reg[estimable == TRUE, p_entry := p.adjust(p_dir, method = "BH"), by = patient]
A_sl[estimable == TRUE,  p_entry := p.adjust(p_dir, method = "BH"), by = patient]
B_use[, p_entry := p_dir]                      # = min(p_enrich, p_deplete)，未校正

# --- 区域级柱值：可估条目的简单均值 + Stouffer 合并 z（检验均值是否偏离 1，方向一致）---
stouffer_p <- function(z) {
  z <- z[!is.na(z)]
  if (length(z) == 0) return(NA_real_)
  2 * pnorm(-abs(sum(z) / sqrt(length(z))))
}
agg_zone <- function(dt, method) {
  grid <- CJ(patient = sort(unique(dt$patient)), zone = ZONE_LV,
             tam = c("C1Q TAMs", "SPP1 TAMs"), sorted = FALSE)
  out <- dt[estimable == TRUE, .(mean_eff = mean(eff), k = .N,
                                 p_comb = stouffer_p(z),
                                 n_slice = uniqueN(slice)),
            by = .(patient, zone, tam)]
  out <- merge(grid, out, by = c("patient", "zone", "tam"), all.x = TRUE)
  # 方法A：家族固定为「该病人的 4 个区域×TAM 柱」（不可估柱无 p，但计入 n=4；预声明家族，不随可估性缩水）
  if (method == "A") out[, p_comb_bh := p.adjust(p_comb, method = "BH", n = 4L), by = patient]
  else               out[, p_comb_bh := p_comb]
  out[, `:=`(method = method, k = fifelse(is.na(k), 0L, k),
             star = fifelse(is.na(mean_eff), "n/a", star_p(p_comb_bh)),
             lab  = fifelse(is.na(mean_eff), "n/a",
                            paste0(fmt_v(mean_eff), "\n", star_p(p_comb_bh))))]
  out[]
}
BAR_A <- agg_zone(A_reg, "A")
BAR_B <- agg_zone(B_use[zone %in% ZONE_LV], "B")
BAR_ALL <- rbind(BAR_A, BAR_B)
BAR_OUT <- BAR_ALL[, .(patient, method, zone, tam, mean_eff = round(mean_eff, 4), k, n_slice,
                       p_comb = signif(p_comb, 4), p_comb_bh = signif(p_comb_bh, 4), star)]
setorder(BAR_OUT, method, patient, zone, tam)
fwrite(BAR_OUT, file.path(TAB_PP, "figCD_region_bar_values.csv"))
cat("[write]", file.path(TAB_PP, "figCD_region_bar_values.csv"), "| rows:", nrow(BAR_OUT), "\n")
cat("[BH family] Method A: n = 4 pre-declared bars per patient (non-estimable bars carry no p)\n\n")

cat("[BAR A] mean enrichment | k | slices | Stouffer p | BH(n=4) | star\n")
print(BAR_A[order(patient, zone, tam), .(patient, zone, tam, mean_eff = round(mean_eff, 2), k, n_slice,
                                         p = signif(p_comb, 3), p_bh = signif(p_comb_bh, 3), star)])
cat("\n[BAR B] mean Ro/e | k | slices | Stouffer p | star\n")
print(BAR_B[order(patient, zone, tam), .(patient, zone, tam, mean_eff = round(mean_eff, 2), k, n_slice,
                                         p = signif(p_comb_bh, 3), star)])

cat("\n[star tally] A:", paste(names(table(BAR_A$star)), table(BAR_A$star), sep = "=", collapse = " "),
    "| B:", paste(names(table(BAR_B$star)), table(BAR_B$star), sep = "=", collapse = " "), "\n")
cat("[n/a bars] A:", BAR_A[star == "n/a", .N], "| B:", BAR_B[star == "n/a", .N], "\n")
cat("[sensitivity] P9 Core SPP1 p_bh under n=4 =",
    signif(BAR_A[patient == "P9" & zone == "Tumor Core" & tam == "SPP1 TAMs", p_comb_bh], 3),
    "(under n=2 it was", signif(BAR_A[patient == "P9" & zone == "Tumor Core" & tam == "SPP1 TAMs", p_comb] * 2, 3), ")\n")
cat("\n[figC b/c source rows] A_sl estimable:", A_sl[estimable == TRUE, .N],
    "| B All estimable:", B_use[zone == "All" & estimable == TRUE, .N], "\n")
cat("[figD source rows] A_reg estimable:", A_reg[estimable == TRUE, .N],
    "| B Core+Interface estimable:", B_use[zone %in% ZONE_LV & estimable == TRUE, .N], "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 绘制 figC：区域级条形(a) + 该病人逐切片 clone × 区域×TAM 热图(b)(c)（4 列），5 病人 × 2 方法
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(patchwork))
if (!exists("save_pp")) stop("save_pp 未定义：请先执行上一 cell")
META <- list(
  A = list(tag = "Method A (radius permutation)", ylab = "Enrichment vs permutation null",
           fill = "log2(enrichment)"),
  B = list(tag = "Method B (kNN k=6 Ro/e)",       ylab = "Ro/e vs kNN background",
           fill = "log2(Ro/e)")
)
TH <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.35),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 7, colour = "black"),
        axis.title = element_text(size = 7.5),
        plot.title = element_text(size = 8, hjust = 0, face = "bold"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.size = unit(3.0, "mm"),
        plot.margin = margin(3, 3, 3, 3))
TH_HEAT <- TH + theme(panel.border = element_blank())

## 4 列：区域 × TAM（列标签换行写区域，避免与 TAM 混淆）
XLEV <- c("Core C1Q", "Core SPP1", "Interface C1Q", "Interface SPP1")
XLAB <- c("Core\nC1QC", "Core\nSPP1", "Interface\nC1QC", "Interface\nSPP1")
mk_xcat <- function(zone, tam) paste0(ifelse(zone == "Tumor Core", "Core", "Interface"), " ",
                                      ifelse(tam == "C1Q TAMs", "C1Q", "SPP1"))
dA_reg <- copy(A_reg)[, xcat := mk_xcat(zone, tam)]
dB_reg <- copy(B_use[zone %in% ZONE_LV])[, xcat := mk_xcat(zone, tam)]

heat4 <- function(d, clones, fill_nm, ptitle, legend_pos = "none") {
  g <- CJ(clone = clones, xcat = XLEV, sorted = FALSE)
  g <- merge(g, d[, .(clone, xcat, eff, p_entry, estimable)], by = c("clone", "xcat"), all.x = TRUE)
  g[is.na(estimable), estimable := FALSE]
  g[, fill_v := fifelse(estimable & !is.na(eff),
                        pmax(pmin(log2(pmax(eff, 1e-6)), 2), -2), NA_real_)]
  g[, tl := fifelse(!estimable | is.na(eff), "n/a",
                    fifelse(star_t(p_entry) == "", fmt_v(eff),
                            paste0(fmt_v(eff), "\n", star_t(p_entry))))]
  bar_mm <- min(70, max(24, 12 * length(clones)))
  ggplot(g, aes(x = xcat, y = clone, fill = fill_v)) +
    geom_tile(colour = NA, width = 1.002, height = 1.002) +
    geom_text(aes(label = tl), size = 2.5, lineheight = 0.78) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-2, 2), breaks = seq(-2, 2, 0.5), oob = scales::squish,
                         na.value = "#E8E8E8", name = fill_nm,
                         guide = guide_colourbar(barheight = unit(bar_mm, "mm"),
                                                 barwidth = unit(2.4, "mm"),
                                                 title.position = "top", ticks = TRUE)) +
    scale_x_discrete(limits = XLEV, labels = setNames(XLAB, XLEV)) +
    scale_y_discrete(limits = rev(clones)) +
    coord_fixed(ratio = 1) +
    labs(title = ptitle, x = NULL, y = NULL) +
    TH_HEAT + theme(legend.position = legend_pos)
}

bar_plot <- function(bd, mt, ptitle) {
  bd <- copy(bd)[, zone_f := factor(zone, levels = ZONE_LV, labels = ZONE_X[ZONE_LV])]
  bd[, tam := factor(tam, levels = c("C1Q TAMs", "SPP1 TAMs"))]
  # top 需同时容纳柱高与两行标签：max*1.30；下界 1.28
  top <- max(1.28, max(c(bd$mean_eff, 1.0), na.rm = TRUE) * 1.30)
  bd[, lab_y := pmax(mean_eff, 1.0) + top * 0.06]
  ok <- bd[!is.na(mean_eff)]; na <- bd[is.na(mean_eff)]
  ggplot(bd, aes(x = zone_f, y = mean_eff, fill = tam)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_col(position = position_dodge(width = 0.78), width = 0.70, na.rm = TRUE) +
    geom_label(data = ok, aes(y = lab_y, label = lab, group = tam),
               position = position_dodge(width = 0.78), vjust = 0,
               fill = "white", colour = "black", linewidth = 0, size = 2.5,
               label.padding = unit(0.30, "mm"), lineheight = 0.78) +
    geom_label(data = na, aes(y = top * 0.45, label = "n/a", group = tam),
               position = position_dodge(width = 0.78), vjust = 0.5,
               fill = "white", colour = "grey35", linewidth = 0, size = 2.5,
               label.padding = unit(0.30, "mm")) +
    scale_fill_manual(values = TAM_COL, labels = TAM_X, name = NULL) +
    scale_y_continuous(limits = c(0, top), expand = c(0, 0)) +
    labs(title = ptitle, x = NULL, y = mt$ylab) +
    TH + theme(legend.position = "bottom", legend.direction = "horizontal",
               legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
               legend.background = element_blank(), legend.key = element_blank())
}

CAP_A <- paste0("Stars: *** p<0.001, ** p<0.01, * p<0.05 (bars: Stouffer-combined z, Method A BH over 4 pre-declared bars;\n",
                "Method B uncorrected) · colour = log2(effect), clipped at +/-2 · n/a = non-estimable · dashed line = 1.0 (no enrichment)")
figC_log <- list()
for (m in c("A", "B")) {
  mt <- META[[m]]
  for (pat in PATIENTS) {
    bd <- BAR_ALL[method == m & patient == pat]
    dd <- if (m == "A") dA_reg[patient == pat] else dB_reg[patient == pat]
    sls <- sort(unique(dd$slice))
    clones <- sort(unique(dd$clone))
    stopifnot(length(sls) == 2)
    hm <- lapply(seq_along(sls), function(i) {
      heat4(dd[slice == sls[i]], clones, mt$fill,
            sprintf("(%s) %s: Core vs Interface", letters[i + 1], sls[i]),
            legend_pos = if (i == length(sls)) "right" else "none")
    })
    p <- bar_plot(bd, mt, "(a) Region level (all slices pooled)") + hm[[1]] + hm[[2]] +
      plot_layout(widths = c(1.10, 1, 1)) +
      plot_annotation(title = sprintf("%s · %s: region enrichment and per-slice clone × region × TAM", mt$tag, pat),
                      caption = CAP_A,
                      theme = theme(plot.title = element_text(size = 9.5, hjust = 0.5),
                                    plot.caption = element_text(size = 7, hjust = 0, colour = "grey25")))
    n_cl <- length(clones)
    w <- 8.6; h <- 1.45 + 0.64 * n_cl
    stem <- sprintf("figC_%s_%s", m, pat)
    save_pp(p, stem, w, h)
    figC_log[[stem]] <- data.table(stem = stem, method = m, patient = pat, n_clone = n_cl,
                                   slices = paste(sls, collapse = "/"), width = w, height = h)
  }
}
FIGC_TAB <- rbindlist(figC_log)
print(FIGC_TAB)
cat("\n[figC png sizes KB]\n")
print(data.table(stem = FIGC_TAB$stem,
                 kb = round(file.info(file.path(FIG_PP, paste0(FIGC_TAB$stem, ".png")))$size/1024, 0)))
cat("[column labels]", paste(gsub("\n", " / ", XLAB), collapse = " | "), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 绘制 figD：逐切片 clone ×（Core C1QC / Core SPP1 / Interface C1QC / Interface SPP1）热图，5 病人 × 2 方法
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(patchwork))
if (!exists("save_pp")) stop("save_pp 未定义：请先执行上一步 cell")
XLEV <- c("Core C1Q", "Core SPP1", "Interface C1Q", "Interface SPP1")
XLAB <- c("Core\nC1QC", "Core\nSPP1", "Interface\nC1QC", "Interface\nSPP1")
mk_xcat <- function(zone, tam) paste0(ifelse(zone == "Tumor Core", "Core", "Interface"), " ",
                                      ifelse(tam == "C1Q TAMs", "C1Q", "SPP1"))
dA_reg <- copy(A_reg)[, xcat := mk_xcat(zone, tam)]
dB_reg <- copy(B_use[zone %in% ZONE_LV])[, xcat := mk_xcat(zone, tam)]

heat_d <- function(d, clones, fill_nm, ptitle, legend_pos = "none") {
  g <- CJ(clone = clones, xcat = XLEV, sorted = FALSE)
  g <- merge(g, d[, .(clone, xcat, eff, p_entry, estimable)],
             by = c("clone", "xcat"), all.x = TRUE)
  g[is.na(estimable), estimable := FALSE]
  g[, fill_v := fifelse(estimable & !is.na(eff),
                        pmax(pmin(log2(pmax(eff, 1e-6)), 2), -2), NA_real_)]
  g[, tl := fifelse(!estimable | is.na(eff), "n/a",
                    fifelse(star_t(p_entry) == "", fmt_v(eff),
                            paste0(fmt_v(eff), "\n", star_t(p_entry))))]
  bar_mm <- min(70, max(24, 12 * length(clones)))
  ggplot(g, aes(x = xcat, y = clone, fill = fill_v)) +
    geom_tile(colour = NA, width = 1.002, height = 1.002) +
    geom_text(aes(label = tl), size = 2.5, lineheight = 0.78) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-2, 2), breaks = seq(-2, 2, 0.5), oob = scales::squish,
                         na.value = "#E8E8E8", name = fill_nm,
                         guide = guide_colourbar(barheight = unit(bar_mm, "mm"),
                                                 barwidth = unit(2.4, "mm"),
                                                 title.position = "top", ticks = TRUE)) +
    scale_x_discrete(limits = XLEV, labels = setNames(XLAB, XLEV)) +
    scale_y_discrete(limits = rev(clones)) +
    coord_fixed(ratio = 1) +
    labs(title = ptitle, x = NULL, y = NULL) +
    TH_HEAT + theme(legend.position = legend_pos)
}

META_D <- list(
  A = list(tag = "Method A (radius permutation)", fill = "log2(enrichment)"),
  B = list(tag = "Method B (kNN k=6 Ro/e)",       fill = "log2(Ro/e)")
)
CAP_D <- paste0("Stars: *** p<0.001, ** p<0.01, * p<0.05 (Method A BH-adjusted per patient; Method B uncorrected)\n",
                "n/a = non-estimable · colour = log2(effect), clipped at +/-2")
figD_log <- list()
for (m in c("A", "B")) {
  mt <- META_D[[m]]
  for (pat in PATIENTS) {
    dd <- if (m == "A") dA_reg[patient == pat] else dB_reg[patient == pat]
    sls <- sort(unique(dd$slice))
    stopifnot(length(sls) == 2)
    clones <- sort(unique(dd$clone))
    hm <- lapply(seq_along(sls), function(i) {
      heat_d(dd[slice == sls[i]], clones, mt$fill,
             sprintf("(%s) %s: Core vs Interface", letters[i], sls[i]),
             legend_pos = if (i == length(sls)) "right" else "none")
    })
    p <- hm[[1]] + hm[[2]] + plot_layout(widths = c(1, 1)) +
      plot_annotation(title = sprintf("%s · %s: clone × region enrichment, per slice", mt$tag, pat),
                      caption = CAP_D,
                      theme = theme(plot.title = element_text(size = 9.5, hjust = 0.5),
                                    plot.caption = element_text(size = 7, hjust = 0, colour = "grey25")))
    n_cl <- length(clones)
    w <- 7.6; h <- 0.95 + 0.78 * n_cl
    stem <- sprintf("figD_%s_%s", m, pat)
    save_pp(p, stem, w, h)
    figD_log[[stem]] <- data.table(stem = stem, method = m, patient = pat, n_clone = n_cl,
                                   slices = paste(sls, collapse = "/"), width = w, height = h)
  }
}
FIGD_TAB <- rbindlist(figD_log)
print(FIGD_TAB)
cat("\n[figD png sizes KB]\n")
print(data.table(stem = FIGD_TAB$stem,
                 kb = round(file.info(file.path(FIG_PP, paste0(FIGD_TAB$stem, ".png")))$size/1024, 0)))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 终局 QC：清理 .ipynb_checkpoints、像素完整性自检（向量化实现，白底/非空/双格式/600dpi）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ck <- list.files(FIG_PP, pattern = "^\\.ipynb_checkpoints$", all.files = TRUE, full.names = TRUE)
if (length(ck)) { unlink(ck, recursive = TRUE); cat("[clean] removed", ck, "\n") }

pngs <- sort(list.files(FIG_PP, pattern = "\\.png$", full.names = TRUE))
EXP_W <- c("figC_A" = 8.6, "figC_B" = 8.6, "figD_A" = 7.6, "figD_B" = 7.6, "P13_old_vs_new" = 6.8)
PIX <- rbindlist(lapply(pngs, function(f) {
  img <- png::readPNG(f)
  d <- dim(img); rgb <- img[, , 1:3, drop = FALSE]
  mx <- pmax(rgb[, , 1], rgb[, , 2], rgb[, , 3])
  mn <- pmin(rgb[, , 1], rgb[, , 2], rgb[, , 3])
  stem <- sub("\\.png$", "", basename(f))
  key <- if (grepl("^P13_old_vs_new", stem)) "P13_old_vs_new" else
         paste(strsplit(stem, "_")[[1]][1:2], collapse = "_")
  data.table(stem = stem, px_w = d[2], px_h = d[1],
             dpi_eff = round(d[2] / EXP_W[[key]]),
             sd_rgb = round(sd(as.vector(rgb)), 4),
             white_frac = round(mean(rgb > 0.98), 3),
             colored_frac = round(mean((mx - mn) > 0.05), 4),
             kb = round(file.info(f)$size / 1024, 0),
             pdf_ok = file.exists(sub("\\.png$", ".pdf", f)))
}))
PIX[, pass := sd_rgb > 0.02 & white_frac > 0.15 & white_frac < 0.99 & colored_frac > 1e-4 &
              kb > 10 & pdf_ok & dpi_eff >= 600]
print(PIX[, .(stem, px_w, px_h, dpi_eff, kb, pdf_ok, pass)])
cat(sprintf("\n[Figure Wholeness] %d/%d figures pixel-verified | failures: %s\n",
            sum(PIX$pass), nrow(PIX),
            ifelse(all(PIX$pass), "none", paste(PIX[pass == FALSE]$stem, collapse = ", "))))
cat(sprintf("[dpi] min effective dpi = %d | min size = %d KB | all pdf paired = %s\n",
            min(PIX$dpi_eff), min(PIX$kb), all(PIX$pdf_ok)))
cat(sprintf("[stem coverage] figC: %d | figD: %d | P13: %d | total stems: %d | dir files: %d\n",
            sum(grepl("^figC_", PIX$stem)), sum(grepl("^figD_", PIX$stem)),
            sum(grepl("^P13_", PIX$stem)), nrow(PIX), length(list.files(FIG_PP))))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")
