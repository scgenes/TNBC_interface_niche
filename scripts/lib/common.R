# =============================================================================
# TNBC interface-niche code repository
# Original source: project_singlecell/TNBC_result/new_analysis/common.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

## ============================================================
## common.R — shared helpers for TNBC tumor-immune boundary validation
##   (A1 spatial zoning / A3 spatial lag cross-correlation / A4 T-cell proximity)
##   Reuses the spatial-lag framework from _archive_legacy/figure12/analysis.R
##   (kNN K=6, within-sample permutation, Fisher-z cross-sample combine)
## ============================================================
suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

DATA_DIR <- "/path/to/project/project_singlecell/TNBC_result"
OUT_DIR  <- file.path(DATA_DIR, "new_analysis")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

region_colors <- c("Tumor Core" = "#C51B7D", "Interface" = "#4DAF4A", "Stroma" = "#377EB8")

## subtype colors (consistent with lineage palette used across figures)
subtype_colors <- c(
  "SPP1 TAMs"  = "#D73027",
  "C1Q TAMs"   = "#4DAF4A",
  "CCL18 TAMs" = "#66C2A5",
  "MMP11 CAFs" = "#4575B4"
)

## ---------- load spot-level subtype proportions + spatial meta ----------
load_data <- function() {
  prop <- readRDS(file.path(DATA_DIR, "figure20", "subproportions.rds"))
  meta <- readRDS(file.path(DATA_DIR, "figure20", "spatial_meta.rds"))
  stopifnot(all(meta$barcode == rownames(prop)))
  list(prop = prop, meta = meta)
}

## ---------- per-spot distance to nearest Tumor Core spot, calibrated to um ----------
## (figure8 method: median nearest-neighbour pixel pitch -> 100 um)
dist_to_core_um <- function(meta) {
  res <- list()
  for (s in unique(meta$sample)) {
    ds <- meta[meta$sample == s, , drop = FALSE]
    xy <- as.matrix(ds[, c("x", "y")])
    dmat_all <- as.matrix(dist(xy))
    diag(dmat_all) <- Inf
    nn_px <- apply(dmat_all, 1, min)
    med_pitch_px <- median(nn_px)
    core <- which(ds$region == "Tumor Core")
    if (length(core) < 3) {
      ds$dist_core_px <- NA_real_
      ds$dist_core_um <- NA_real_
      ds$med_pitch_px <- med_pitch_px
      res[[s]] <- ds; next
    }
    dmat_core <- dmat_all[, core, drop = FALSE]
    ds$dist_core_px <- apply(dmat_core, 1, min)
    ds$dist_core_um <- ds$dist_core_px / med_pitch_px * 100
    ds$med_pitch_px <- med_pitch_px
    res[[s]] <- ds
  }
  do.call(rbind, res)
}

## ---------- kNN (k=6) row-normalized spatial weight matrix, per sample ----------
## (figure12 build_W: symmetric kNN, row-standardized)
build_W <- function(xy, sample_id, k = 6) {
  N <- nrow(xy)
  uniq <- unique(sample_id)
  I <- integer(0); J <- integer(0); V <- numeric(0)
  for (s in uniq) {
    idx <- which(sample_id == s); n <- length(idx)
    d <- as.matrix(dist(xy[idx, , drop = FALSE])); diag(d) <- Inf
    adj <- matrix(0, n, n)
    for (i in seq_len(n)) adj[i, order(d[i, ])[seq_len(min(k, n - 1))]] <- 1
    adj <- (adj + t(adj) > 0) * 1
    rs <- pmax(rowSums(adj), 1)
    adj <- sweep(adj, 1, rs, "/")
    gi <- which(adj != 0, arr.ind = TRUE)
    I <- c(I, idx[gi[, 1]]); J <- c(J, idx[gi[, 2]]); V <- c(V, adj[gi])
  }
  Matrix::sparseMatrix(i = I, j = J, x = V, dims = c(N, N))
}

fisher_combine <- function(rhos) tanh(mean(atanh(rhos), na.rm = TRUE))

## ---------- p-value formatting (G14D threshold notation) ----------
fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 1e-15) return("<1e-15")
  if (p < 1e-10) return(format(p, digits = 2, scientific = TRUE))
  if (p < 0.001) return(format(p, digits = 3, scientific = TRUE))
  format(round(p, 4), nsmall = 4)
}
