# =============================================================================
# TNBC interface-niche code repository
# Original source: project_singlecell/TNBC_result/final_figures/V3/scripts/figure3/data/spatial_niche.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

## ============================================================
## Figure 3 — the spatial niche of C1Q+ TAMs versus SPP1+ TAMs.
##
## Three questions, one framework (kNN K=6 spatial lag, within-sample
## permutation, Fisher-z cross-sample combine — identical to
## new_analysis/A3_spatial_lag.R so the numbers are comparable):
##
##   (1) CAF-subtype specificity of the C1Q+ TAM niche
##       sender  = each of the six CAF subtypes
##       receiver= C1Q+ TAM  and  SPP1+ TAM
##       6 x 2 = 12 tests.  The earlier run tested MMP11 alone; without the
##       other five CAFs "MMP11 is the partner" is not a claim, only a guess.
##
##   (2) TAM -> T-cell coupling, per TAM substate
##       sender  = each of the seven TAM substates (+ all non-SPP1 TAMs pooled)
##       receiver= total T cells
##       8 tests.  The pooled "other TAM" row is the contrast partner: T cells
##       sit next to TAMs in general, and SPP1 is the state that breaks it.
##
##   (3) spot-level export for the representative-section maps
##       (C1Q TAM, MMP11 CAF, SPP1 TAM, T cells on the same coordinate frame)
##
## Everything is written to panels/figure3/, which is where the Python panels
## read from; nothing here draws a figure.
##
##   Rscript spatial_niche.R
## ============================================================
source(file.path("scripts", "lib", "common.R"))
set.seed(42)

OUT <- "/path/to/project/project_singlecell/TNBC_result/final_figures/V3/panels/figure3"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

d <- load_data()
prop <- d$prop
meta <- d$meta
meta$sname <- sub("^GSM[0-9]+_", "", meta$sample)

sample_id <- meta$sample
xy <- as.matrix(meta[, c("x", "y")])
W <- build_W(xy, sample_id, k = 6)
cat(sprintf("[F3] W: %d x %d | nonzeros: %d\n", nrow(W), ncol(W), length(W@x)))

TAM_COLS <- c("FCN1 TAMs", "ISG15 TAMs", "CCL18 TAMs", "C1Q TAMs",
              "SPP1 TAMs", "MKI67 TAMs", "MMP9 TAMs")
CAF_COLS <- c("APOD CAFs", "CD74 CAFs", "CFD CAFs", "IGF1 CAFs",
              "MMP11 CAFs", "STMN1 CAFs")
TCELL_COLS <- grep("T cells$", colnames(prop), value = TRUE)
stopifnot(all(c(TAM_COLS, CAF_COLS) %in% colnames(prop)))
cat("[F3] T-cell columns (", length(TCELL_COLS), "): ",
    paste(TCELL_COLS, collapse = " | "), "\n", sep = "")

tcell_total <- rowSums(prop[, TCELL_COLS, drop = FALSE])
other_tam <- rowSums(prop[, setdiff(TAM_COLS, "SPP1 TAMs"), drop = FALSE])

## ---------- spatial lag + within-sample permutation (A3-compatible) ----------
run_spatial_lag <- function(S, R, sample_id, W, nperm = 1000) {
  SL <- as.numeric(W %*% S)
  uniq <- unique(sample_id)
  idxs <- lapply(uniq, function(s) which(sample_id == s)); names(idxs) <- uniq
  xs <- lapply(idxs, function(i) { x <- rank(SL[i]); x - mean(x) })
  ys <- lapply(idxs, function(i) { y <- rank(R[i]);  y - mean(y) })
  denom <- sapply(seq_along(uniq),
                  function(j) sqrt(sum(xs[[j]]^2)) * sqrt(sum(ys[[j]]^2)))
  valid <- denom > 1e-12
  rho_comb <- function() {
    out <- sapply(seq_along(uniq),
                  function(j) if (valid[j]) sum(xs[[j]] * ys[[j]]) / denom[j] else NA_real_)
    fisher_combine(out)
  }
  rho_obs <- rho_comb()
  per_sample_rho <- sapply(seq_along(uniq), function(j)
    if (valid[j]) cor(SL[idxs[[j]]], R[idxs[[j]]], method = "spearman") else NA_real_)
  null <- numeric(nperm)
  for (p in seq_len(nperm)) {
    for (j in seq_along(uniq)) if (valid[j]) ys[[j]] <- sample(ys[[j]])
    null[p] <- rho_comb()
  }
  list(rho = rho_obs,
       p_hi = (1 + sum(null >= rho_obs)) / (1 + nperm),
       p_lo = (1 + sum(null <= rho_obs)) / (1 + nperm),
       null_mean = mean(null), null_sd = sd(null),
       n_valid = sum(valid), per_sample_rho = per_sample_rho, samples = uniq)
}

## ---------- (1) CAF subtype -> TAM substate ----------
cat("\n[F3] (1) CAF subtype -> TAM substate spatial lag\n")
caf_rows <- list(); caf_ps <- list()
for (caf in CAF_COLS) {
  S <- as.numeric(prop[[caf]])
  for (tam in c("C1Q TAMs", "SPP1 TAMs")) {
    r <- run_spatial_lag(S, as.numeric(prop[[tam]]), sample_id, W, nperm = 1000)
    cat(sprintf("  %-12s -> %-11s rho=% .4f  p_hi=%s\n",
                caf, tam, r$rho, fmt_p(r$p_hi)))
    caf_rows[[length(caf_rows) + 1]] <- data.frame(
      sender = caf, receiver = tam, rho = r$rho,
      perm_p = r$p_hi, perm_p_neg = r$p_lo,
      null_mean = r$null_mean, null_sd = r$null_sd,
      n_valid = r$n_valid,
      n_positive = sum(r$per_sample_rho > 0, na.rm = TRUE),
      stringsAsFactors = FALSE)
    caf_ps[[paste(caf, tam)]] <- data.frame(
      sample = sub("^GSM[0-9]+_", "", r$samples),
      sender = caf, receiver = tam, rho = r$per_sample_rho)
  }
}
caf_res <- do.call(rbind, caf_rows)
caf_res$padj <- p.adjust(caf_res$perm_p, method = "BH")
caf_res$supported <- caf_res$padj < 0.05 & caf_res$rho > 0
write.csv(caf_res, file.path(OUT, "fig3_caf_tam_lag.csv"), row.names = FALSE)
write.csv(do.call(rbind, caf_ps), file.path(OUT, "fig3_caf_tam_lag_per_sample.csv"),
          row.names = FALSE)
cat("\n[F3] CAF x TAM matrix, BH-corrected:\n")
print(caf_res[, c("sender", "receiver", "rho", "perm_p", "padj",
                  "n_positive", "supported")], row.names = FALSE, digits = 3)

## ---------- (1b) paired contrast, per sample: C1Q vs SPP1 in the MMP11 field ----------
caf_ps <- do.call(rbind, caf_ps)
mm <- reshape(caf_ps[caf_ps$sender == "MMP11 CAFs",
                     c("sample", "receiver", "rho")],
              idvar = "sample", timevar = "receiver", direction = "wide")
names(mm) <- sub("^rho\\.", "", names(mm))
mm$contrast <- mm[["C1Q TAMs"]] - mm[["SPP1 TAMs"]]
write.csv(mm, file.path(OUT, "fig3_mmp11_c1q_spp1_contrast.csv"), row.names = FALSE)
npos <- sum(mm$contrast > 0, na.rm = TRUE)
cat(sprintf("\n[F3] MMP11 field, C1Q minus SPP1 (per sample): mean=% .4f, %d/%d positive, sign p=%s\n",
            mean(mm$contrast), npos, nrow(mm),
            fmt_p(binom.test(npos, nrow(mm), alternative = "greater")$p.value)))
cat(sprintf("[F3]   paired Wilcoxon p=%s\n",
            fmt_p(wilcox.test(mm[["C1Q TAMs"]], mm[["SPP1 TAMs"]], paired = TRUE,
                              exact = FALSE)$p.value)))
write.csv(data.frame(
  claim = "MMP11_field_C1Q_minus_SPP1", mean_contrast = mean(mm$contrast),
  n_positive = npos, n_samples = nrow(mm),
  sign_p = binom.test(npos, nrow(mm), alternative = "greater")$p.value,
  wilcox_p = wilcox.test(mm[["C1Q TAMs"]], mm[["SPP1 TAMs"]], paired = TRUE,
                         exact = FALSE)$p.value),
  file.path(OUT, "fig3_mmp11_c1q_spp1_contrast_summary.csv"), row.names = FALSE)

## ---------- (2) TAM substate -> T cells ----------
cat("\n[F3] (2) TAM substate -> total T cells spatial lag\n")
tam_rows <- list(); tam_ps <- list()
for (tam in c(TAM_COLS, "__other__")) {
  S <- if (tam == "__other__") other_tam else as.numeric(prop[[tam]])
  r <- run_spatial_lag(S, tcell_total, sample_id, W, nperm = 1000)
  lab <- if (tam == "__other__") "Other TAMs (pooled)" else tam
  cat(sprintf("  %-22s -> T cells rho=% .4f  p_hi=%s  p_lo=%s\n",
              lab, r$rho, fmt_p(r$p_hi), fmt_p(r$p_lo)))
  tam_rows[[length(tam_rows) + 1]] <- data.frame(
    sender = lab, receiver = "T cells", rho = r$rho,
    perm_p = r$p_hi, perm_p_neg = r$p_lo,
    null_mean = r$null_mean, null_sd = r$null_sd,
    n_valid = r$n_valid,
    n_positive = sum(r$per_sample_rho > 0, na.rm = TRUE),
    stringsAsFactors = FALSE)
  tam_ps[[lab]] <- data.frame(
    sample = sub("^GSM[0-9]+_", "", r$samples),
    sender = lab, rho = r$per_sample_rho)
}
tam_res <- do.call(rbind, tam_rows)
tam_res$padj <- p.adjust(tam_res$perm_p, method = "BH")
tam_res$supported <- tam_res$padj < 0.05 & tam_res$rho > 0
write.csv(tam_res, file.path(OUT, "fig3_tam_tcell_lag.csv"), row.names = FALSE)
ps <- do.call(rbind, tam_ps)
write.csv(ps, file.path(OUT, "fig3_tam_tcell_lag_per_sample.csv"), row.names = FALSE)
cat("\n[F3] TAM -> T cells, BH-corrected:\n")
print(tam_res[, c("sender", "rho", "perm_p", "padj", "n_positive", "supported")],
      row.names = FALSE, digits = 3)

## ---------- (2b) paired contrast: SPP1 vs pooled other TAMs ----------
wide <- reshape(ps[, c("sample", "sender", "rho")], idvar = "sample",
                timevar = "sender", direction = "wide")
names(wide) <- sub("^rho\\.", "", names(wide))
con <- data.frame(sample = wide$sample,
                  rho_other = wide[["Other TAMs (pooled)"]],
                  rho_spp1  = wide[["SPP1 TAMs"]])
con$contrast <- con$rho_other - con$rho_spp1
write.csv(con, file.path(OUT, "fig3_spp1_tcell_contrast.csv"), row.names = FALSE)
npos <- sum(con$contrast > 0, na.rm = TRUE)
sign_p <- binom.test(npos, nrow(con), p = 0.5, alternative = "greater")$p.value
wil <- wilcox.test(con$rho_other, con$rho_spp1, paired = TRUE, exact = FALSE)
cat(sprintf("\n[F3] Paired contrast other-TAM minus SPP1: mean=% .4f, %d/%d positive, sign p=%s, Wilcoxon p=%s\n",
            mean(con$contrast), npos, nrow(con), fmt_p(sign_p), fmt_p(wil$p.value)))
contrast_summary <- data.frame(
  claim = "other_TAM_minus_SPP1_vs_Tcell",
  mean_contrast = mean(con$contrast), median_contrast = median(con$contrast),
  n_positive = npos, n_samples = nrow(con),
  sign_p = sign_p, wilcox_p = wil$p.value)
write.csv(contrast_summary, file.path(OUT, "fig3_spp1_tcell_contrast_summary.csv"),
          row.names = FALSE)

## ---------- (2c) T-cell neighbourhood enrichment (A4 statistic) ----------
## For every T-cell-enriched spot (global Q75), the fraction of its six
## nearest neighbours that are X-enriched; null = permute the X labels within
## each sample, 1000 draws.  Different from the lag above: threshold-based,
## count-based, and it conditions on T-cell spots rather than on the whole
## section.  The full null is exported so the panel can draw it, not just the z.
cat("\n[F3] (2c) T-cell neighbourhood enrichment\n")
K <- 6
thr_t <- unname(quantile(tcell_total, 0.75))
is_t <- tcell_total > thr_t
samp <- list()
for (s in unique(sample_id)) {
  idx <- which(sample_id == s)
  dd <- as.matrix(dist(xy[idx, , drop = FALSE])); diag(dd) <- Inf
  kk <- min(K, length(idx) - 1)
  samp[[s]] <- list(idx = idx, nbr = t(apply(dd, 1, function(r) order(r)[seq_len(kk)])),
                    t_pos = which(is_t[idx]))
}
cat(sprintf("[F3] T-cell-enriched spots (Q75=%.4f): %d\n", thr_t, sum(is_t)))

compute_E <- function(is_X) {
  vals <- numeric(0)
  for (s in names(samp)) {
    sn <- samp[[s]]
    if (length(sn$t_pos) == 0) next
    Xf <- as.numeric(is_X[sn$idx])
    nb <- sn$nbr[sn$t_pos, , drop = FALSE]
    vals <- c(vals, rowMeans(matrix(Xf[nb], nrow = nrow(nb))))
  }
  mean(vals)
}
permute_E <- function(is_X, nperm = 1000) {
  obs <- compute_E(is_X)
  null <- numeric(nperm)
  for (p in seq_len(nperm)) {
    is_Xp <- is_X
    for (s in names(samp)) {
      idx <- samp[[s]]$idx
      is_Xp[idx] <- sample(is_X[idx])
    }
    null[p] <- compute_E(is_Xp)
  }
  list(obs = obs, null = null,
       p = (1 + sum(null >= obs)) / (1 + nperm),
       p_lo = (1 + sum(null <= obs)) / (1 + nperm),
       z = (obs - mean(null)) / sd(null), null_mean = mean(null),
       null_sd = sd(null))
}
knn_types <- c("SPP1 TAMs", "C1Q TAMs", "MMP11 CAFs")
knn_rows <- list(); null_rows <- list()
for (ty in knn_types) {
  is_X <- as.numeric(prop[[ty]]) > unname(quantile(prop[[ty]], 0.75))
  r <- permute_E(is_X)
  cat(sprintf("  %-11s obs=%.4f  null=%.4f  z=% .2f  p_hi=%s  p_lo=%s\n",
              ty, r$obs, r$null_mean, r$z, fmt_p(r$p), fmt_p(r$p_lo)))
  knn_rows[[ty]] <- data.frame(
    type = ty, obs = r$obs, null_mean = r$null_mean, null_sd = r$null_sd,
    z = r$z, perm_p_enriched = r$p, perm_p_depleted = r$p_lo,
    n_enriched_spots = sum(is_X), stringsAsFactors = FALSE)
  null_rows[[ty]] <- data.frame(type = ty, draw = seq_along(r$null), value = r$null)
}
knn_res <- do.call(rbind, knn_rows)
knn_res$padj <- p.adjust(pmin(knn_res$perm_p_enriched, knn_res$perm_p_depleted),
                         method = "BH")
write.csv(knn_res, file.path(OUT, "fig3_knn_tcell_enrichment.csv"), row.names = FALSE)
write.csv(do.call(rbind, null_rows), file.path(OUT, "fig3_knn_tcell_null.csv"),
          row.names = FALSE)

## ---------- (3) spot-level export for the maps ----------
cat("\n[F3] (3) spot-level export\n")
spot <- data.frame(
  barcode = meta$barcode, sname = meta$sname,
  sample = meta$sample, x = meta$x, y = meta$y, region = meta$region,
  c1q = as.numeric(prop[["C1Q TAMs"]]),
  spp1 = as.numeric(prop[["SPP1 TAMs"]]),
  mmp11 = as.numeric(prop[["MMP11 CAFs"]]),
  tcell = tcell_total,
  caf_sum = rowSums(prop[, CAF_COLS, drop = FALSE]),
  tam_sum = rowSums(prop[, TAM_COLS, drop = FALSE]),
  stringsAsFactors = FALSE)
stopifnot(nrow(spot) == nrow(meta), !anyDuplicated(spot$barcode))
write.csv(spot, file.path(OUT, "fig3_spots.csv"), row.names = FALSE)
cat(sprintf("[F3] spots: %d across %d samples\n",
            nrow(spot), length(unique(spot$sname))))
print(sapply(spot[, c("c1q", "spp1", "mmp11", "tcell", "caf_sum", "tam_sum")],
             quantile, probs = c(0.5, 0.75, 0.9, 1.0), na.rm = TRUE))

cat("\n[F3] Done.\n")
