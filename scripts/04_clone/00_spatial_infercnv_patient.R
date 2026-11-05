# =============================================================================
# TNBC interface-niche code repository
# Original source: project_spatial/TNBC/TNBC_analysis/SpatialInferCNV/results/step1_three_zone/notebooks/step01_three_zone.ipynb (converted to R)
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

# 初始化：加载 infercnv（读取 S4 对象必需）、定义路径与 qs/RDS 保存工具
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
suppressPackageStartupMessages({
  library(infercnv)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})
set.seed(20261020L)

BASE      <- "/path/to/project/project_spatial/TNBC/TNBC_analysis/SpatialInferCNV"
CNV_ROOT  <- "/path/to/data/project_spatial_data/TNBC/SpatialInferCNV/patient_cnv/results"
OUT       <- file.path(BASE, "results/step1_three_zone")
REVIEW    <- file.path(BASE, "review/step1_three_zone")
for (d in c(file.path(OUT, "data"), file.path(OUT, "tables"), file.path(OUT, "figures"), REVIEW)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

PATIENTS <- c("P9", "P11", "P12", "P13", "P14", "P15")
SLICE_EXPECT <- list(
  P9  = c("092A", "092B"),
  P11 = c("093B"),
  P12 = c("093C", "093D"),
  P13 = c("094A", "094B"),
  P14 = c("094C", "094D"),
  P15 = c("095A", "095B")
)
SLICES <- unlist(SLICE_EXPECT, use.names = FALSE)
plot_theme <- "nature"   # 项目主题（sc-state plot_theme = nature）

save_obj <- function(obj, filepath) {
  qs_path <- paste0(sub("\\.(rds|qs)$", "", filepath), ".qs")
  if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(obj, qs_path, preset = "fast")
    cat(sprintf("[Save] qs: %s (%.1f MB)\n", qs_path, file.info(qs_path)$size / 1e6))
  } else {
    rds_path <- paste0(sub("\\.qs$", "", qs_path), ".rds")
    saveRDS(obj, rds_path)
    cat(sprintf("[Save] RDS fallback: %s\n", rds_path))
  }
}
load_obj <- function(filepath) {
  qs_path <- paste0(sub("\\.(rds|qs)$", "", filepath), ".qs")
  rds_path <- paste0(sub("\\.qs$", "", qs_path), ".rds")
  if (file.exists(qs_path) && requireNamespace("qs", quietly = TRUE)) {
    cat(sprintf("[Load] qs: %s\n", qs_path)); qs::qread(qs_path)
  } else if (file.exists(rds_path)) {
    cat(sprintf("[Load] RDS: %s\n", rds_path)); readRDS(rds_path)
  } else stop(sprintf("File not found: %s or %s", qs_path, rds_path))
}

cat("R:", R.version.string, "| infercnv:", as.character(packageVersion("infercnv")),
    "| qs:", as.character(packageVersion("qs")), "\n")
cat("患者数:", length(PATIENTS), "| 切片数:", length(SLICES), "\n")
cat("输出目录:", OUT, "\n内部证据:", REVIEW, "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 读取基因 chr/start 注释，按着丝粒位置把每个基因划入染色体臂（CENTRO 复用参考实现）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
GENE_ORD <- file.path(BASE, "input/gene_order_filtered.tsv")
go <- read.delim(GENE_ORD, header = FALSE,
                 col.names = c("gene", "chr", "start", "end"),
                 stringsAsFactors = FALSE)

# 染色体着丝粒位置 (Mb)，与 run_clone_patient.R:243 完全一致
CENTRO <- c(chr1 = 125, chr2 = 93.3, chr3 = 91, chr4 = 50.4, chr5 = 48.4, chr6 = 60.3,
            chr7 = 60.1, chr8 = 45.6, chr9 = 49, chr10 = 40.2, chr11 = 53.7, chr12 = 35.8,
            chr13 = 17.9, chr14 = 17, chr15 = 19, chr16 = 36.8, chr17 = 25.1, chr18 = 18.5,
            chr19 = 26.1, chr20 = 28.1, chr21 = 12, chr22 = 15, chrX = 61, chrY = 10.4)

go <- go[!duplicated(go$gene), , drop = FALSE]
go$chr   <- as.character(go$chr)
go$start <- as.numeric(go$start)
go$arm   <- ifelse(go$chr %in% names(CENTRO),
                   paste0(go$chr, ifelse(go$start / 1e6 < CENTRO[go$chr], "p", "q")),
                   NA_character_)

ARMS <- unique(go$arm[!is.na(go$arm)])
arm_chr <- sub("[pq]$", "", ARMS)
arm_side <- sub("^.*[0-9XY]([pq])$", "\\1", ARMS)
arm_num <- suppressWarnings(as.numeric(sub("^chr", "", arm_chr)))
arm_num[arm_chr == "chrX"] <- 23
arm_num[arm_chr == "chrY"] <- 24
ARMS <- ARMS[order(arm_num, arm_side)]
arm_n <- table(factor(go$arm, levels = ARMS))

cat("基因注释表:", nrow(go), "个唯一基因\n")
cat("可定义染色体臂的基因:", sum(!is.na(go$arm)), "| 无法定义(非主染色体/缺着丝粒):",
    sum(is.na(go$arm)), "\n")
cat("染色体臂数:", length(ARMS), "| 覆盖染色体:",
    paste(unique(arm_chr), collapse = " "), "\n")
cat("每条臂基因数(前 10):", paste(sprintf("%s=%d", names(arm_n)[1:10], arm_n[1:10]), collapse = ", "), "\n")
cat("无着丝粒定义的 chr 取值:", paste(sort(unique(go$chr[is.na(go$arm)])), collapse = ", "), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 预检：读取一个病人的 inferCNV 终态对象，核对 S4 结构、维度、观测/参考 spot 与命名格式
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
f_probe <- file.path(CNV_ROOT, "P15", "run.final.infercnv_obj")
o <- readRDS(f_probe)
cat("class:", paste(class(o), collapse = "/"), "\n")
cat("slots:", paste(slotNames(o), collapse = ", "), "\n")
ex <- o@expr.data
cat("expr.data:", nrow(ex), "genes x", ncol(ex), "cells | class:", class(ex)[1], "\n")
cat("genes head:", paste(head(rownames(ex), 5), collapse = ", "), "\n")
cat("cells head:", paste(head(colnames(ex), 5), collapse = ", "), "\n")
obs_idx <- unlist(o@observation_grouped_cell_indices, use.names = FALSE)
ref_idx <- unlist(o@reference_grouped_cell_indices, use.names = FALSE)
cat("obs groups:", length(o@observation_grouped_cell_indices),
    "| obs spots:", length(obs_idx),
    "| ref groups:", length(o@reference_grouped_cell_indices),
    "| ref spots:", length(ref_idx), "\n")
obs_cells <- colnames(ex)[obs_idx]
tb_pref <- table(sub("_.*$", "", obs_cells))
cat("obs slice prefixes:", paste(sprintf("%s=%d", names(tb_pref), as.integer(tb_pref)), collapse = " | "), "\n")
cat("expr value range (first 50 obs spots):",
    paste(round(range(ex[, obs_idx[seq_len(min(50, length(obs_idx)))]], na.rm = TRUE), 4), collapse = " .. "), "\n")
cat("genes with chr/start annotation:", sum(rownames(ex) %in% go$gene), "/", nrow(ex), "\n")
rm(o, ex); invisible(gc(verbose = FALSE))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 输入一致性核查：逐个病人比对 inferCNV 对象的观测 spot 与 clone_membership 表（不计算 CNV）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
chk <- list()
for (pid in PATIENTS) {
  o <- readRDS(file.path(CNV_ROOT, pid, "run.final.infercnv_obj"))
  ex <- o@expr.data
  obs_cells <- colnames(ex)[unlist(o@observation_grouped_cell_indices, use.names = FALSE)]
  ref_n <- length(unlist(o@reference_grouped_cell_indices, use.names = FALSE))
  mem <- read.csv(file.path(BASE, "patient_clone", pid, "clone_membership.csv"),
                  stringsAsFactors = FALSE)
  tb_obj <- table(sub("_.*$", "", obs_cells))
  tb_mem <- table(mem$slice)
  sls <- union(names(tb_obj), names(tb_mem))
  chk[[pid]] <- data.frame(
    patient = pid, slice = sls,
    n_obs_obj = as.integer(tb_obj[sls]), n_membership = as.integer(tb_mem[sls]),
    n_ref_obj = ref_n, n_gene_obj = nrow(ex),
    n_cell_overlap = vapply(sls, function(s)
      sum(obs_cells[sub("_.*$", "", obs_cells) == s] %in% mem$cell[mem$slice == s]), 0L),
    stringsAsFactors = FALSE)
  rm(o, ex); invisible(gc(verbose = FALSE))
}
chk <- do.call(rbind, chk); rownames(chk) <- NULL
chk$match <- ifelse(chk$n_obs_obj == chk$n_membership & chk$n_cell_overlap == chk$n_membership,
                    "OK", "MISMATCH")
print(chk, row.names = FALSE)
cat("\n总计: obj 观测 spot =", sum(chk$n_obs_obj),
    "| membership spot =", sum(chk$n_membership),
    "| 重叠 =", sum(chk$n_cell_overlap), "\n")
cat("不一致切片数:", sum(chk$match != "OK"), "/", nrow(chk), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 输入路径核对(只读)：检查 samples/ 下逐切片 inferCNV 对象的规模与 clone_membership 是否一致
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
SIF_ROOT    <- "/path/to/data/project_spatial_data/TNBC/SpatialInferCNV"
SAMPLE_ROOT <- file.path(SIF_ROOT, "samples")
sample_rows <- list()
for (pid in PATIENTS) {
  for (sl in SLICE_EXPECT[[pid]]) {
    f <- file.path(SAMPLE_ROOT, sl, "run.final.infercnv_obj")
    if (!file.exists(f)) {
      sample_rows[[length(sample_rows) + 1]] <- data.frame(
        patient = pid, slice = sl, exists = FALSE, n_obs_obj = NA_integer_,
        n_ref_obj = NA_integer_, n_gene_obj = NA_integer_, n_membership = NA_integer_,
        n_cell_overlap = NA_integer_, stringsAsFactors = FALSE)
      next
    }
    o <- readRDS(f)
    ex <- o@expr.data
    obs_cells <- colnames(ex)[unlist(o@observation_grouped_cell_indices, use.names = FALSE)]
    ref_n <- length(unlist(o@reference_grouped_cell_indices, use.names = FALSE))
    mem <- read.csv(file.path(BASE, "patient_clone", pid, "clone_membership.csv"),
                    stringsAsFactors = FALSE)
    mem_sl <- mem[mem$slice == sl, ]
    sample_rows[[length(sample_rows) + 1]] <- data.frame(
      patient = pid, slice = sl, exists = TRUE,
      n_obs_obj = length(obs_cells), n_ref_obj = ref_n, n_gene_obj = nrow(ex),
      n_membership = nrow(mem_sl),
      n_cell_overlap = sum(obs_cells %in% mem_sl$cell),
      stringsAsFactors = FALSE)
    rm(o, ex); invisible(gc(verbose = FALSE))
  }
}
sample_chk <- do.call(rbind, sample_rows); rownames(sample_chk) <- NULL
sample_chk$match <- with(sample_chk,
  ifelse(exists & n_obs_obj == n_membership & n_cell_overlap == n_membership, "OK", "MISMATCH"))
print(sample_chk, row.names = FALSE)
cat("\nsamples/ 路径:", SAMPLE_ROOT, "\n")
cat("samples/ 观测 spot 总数:", sum(sample_chk$n_obs_obj, na.rm = TRUE),
    "| membership 总数:", sum(sample_chk$n_membership, na.rm = TRUE), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 逐切片读取 samples/ 全量 inferCNV 对象，按染色体臂计算 mean(|CNV-1|) 并跨臂取均值
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
SAMPLE_ROOT <- "/path/to/data/project_spatial_data/TNBC/SpatialInferCNV/samples"
SLICE2PAT  <- setNames(rep(names(SLICE_EXPECT), lengths(SLICE_EXPECT)), unlist(SLICE_EXPECT))
gene_arm_ok <- go$gene[!is.na(go$arm)]

score_rows <- list()
read_log   <- list()
for (sl in SLICES) {
  t0 <- Sys.time()
  pid <- SLICE2PAT[[sl]]
  f <- file.path(SAMPLE_ROOT, sl, "run.final.infercnv_obj")
  if (!file.exists(f)) stop("missing inferCNV object: ", f)
  o  <- readRDS(f)
  ex <- o@expr.data
  obs_cells <- colnames(ex)[unlist(o@observation_grouped_cell_indices, use.names = FALSE)]
  ref_n     <- length(unlist(o@reference_grouped_cell_indices, use.names = FALSE))
  pref <- sub("_.*$", "", obs_cells)
  if (any(pref != sl)) stop(sprintf("%s: cell prefix mismatch (%s)", sl, paste(unique(pref), collapse = ",")))
  if (anyDuplicated(obs_cells)) stop(sprintf("%s: duplicated obs cells", sl))

  genes_use <- rownames(ex)[rownames(ex) %in% gene_arm_ok]
  arm_score <- matrix(NA_real_, nrow = length(ARMS), ncol = length(obs_cells),
                      dimnames = list(ARMS, obs_cells))
  for (a in ARMS) {
    ga <- go$gene[go$arm == a & !is.na(go$arm)]
    ga <- intersect(ga, genes_use)
    if (!length(ga)) next
    arm_score[a, ] <- colMeans(abs(ex[ga, obs_cells, drop = FALSE] - 1))
  }
  n_arm_used <- sum(is.finite(arm_score[, 1]))
  score <- colMeans(arm_score, na.rm = TRUE)
  if (any(!is.finite(score))) stop(sprintf("%s: non-finite CNV score", sl))

  score_rows[[sl]] <- data.frame(patient = pid, slice = sl, cell = obs_cells,
                                 cnv_arm_score = as.numeric(score), stringsAsFactors = FALSE)
  read_log[[sl]] <- data.frame(
    patient = pid, slice = sl, n_obs = length(obs_cells), n_ref = ref_n,
    n_gene_object = nrow(ex), n_gene_used = length(genes_use),
    n_arm_used = n_arm_used, stringsAsFactors = FALSE)
  cat(sprintf("%s (%s): obs=%d ref=%d | genes=%d/%d | arms=%d | score min/med/max=%.4f/%.4f/%.4f | %.1fs\n",
              sl, pid, length(obs_cells), ref_n, length(genes_use), nrow(ex), n_arm_used,
              min(score), median(score), max(score),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  rm(o, ex, arm_score); invisible(gc(verbose = FALSE))
  cat(sprintf("  [memory] remaining objects: %d; %.0f MB in use\n", length(score_rows), sum(gc()[, 2])))
}
spot <- do.call(rbind, score_rows); rownames(spot) <- NULL
read_log <- do.call(rbind, read_log); rownames(read_log) <- NULL
cat("\n合计观测 spot:", nrow(spot), "| 切片数:", length(unique(spot$slice)),
    "| 细胞名重复:", sum(duplicated(spot$cell)), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 切片内 min–max 归一化并按 0.2 / 0.5 赋三区（不跨切片合并、不改阈值）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
spot$cnv_norm <- NA_real_
norm_log <- list()
for (sl in SLICES) {
  i  <- which(spot$slice == sl)
  v  <- spot$cnv_arm_score[i]
  rg <- range(v, na.rm = TRUE)
  if (!all(is.finite(rg))) stop(sprintf("%s: non-finite score range", sl))
  if (rg[2] - rg[1] <= 0) {
    warning(sprintf("%s: degenerate score distribution (range=0); cnv_norm left NA and reported", sl))
    norm_log[[sl]] <- data.frame(slice = sl, n_spot = length(i), min = rg[1], max = rg[2],
                                 range = rg[2] - rg[1], degenerate = TRUE, stringsAsFactors = FALSE)
    next
  }
  spot$cnv_norm[i] <- (v - rg[1]) / (rg[2] - rg[1])
  norm_log[[sl]] <- data.frame(slice = sl, n_spot = length(i), min = rg[1], max = rg[2],
                               range = rg[2] - rg[1], degenerate = FALSE, stringsAsFactors = FALSE)
}
spot$zone <- ifelse(is.na(spot$cnv_norm), NA_character_,
             ifelse(spot$cnv_norm < 0.2, "Stroma",
             ifelse(spot$cnv_norm <= 0.5, "Interface", "Core")))
spot$zone <- factor(spot$zone, levels = c("Stroma", "Interface", "Core"))
norm_log <- do.call(rbind, norm_log); rownames(norm_log) <- NULL

zone_by_slice <- as.data.frame.matrix(table(spot$slice, spot$zone))
zone_by_slice <- cbind(slice = rownames(zone_by_slice), zone_by_slice)
rownames(zone_by_slice) <- NULL
zone_by_slice$n_spot <- rowSums(zone_by_slice[, c("Stroma", "Interface", "Core")])
print(zone_by_slice, row.names = FALSE)
cat("\n归一化后分数范围(应全部为 [0,1]):",
    paste(round(range(spot$cnv_norm, na.rm = TRUE), 6), collapse = " .. "), "\n")
cat("三区 spot 合计:", paste(sprintf("%s=%d", levels(spot$zone), as.integer(table(spot$zone))), collapse = " | "), "\n")
cat("空区切片数:", sum(zone_by_slice$Stroma == 0 | zone_by_slice$Interface == 0 | zone_by_slice$Core == 0), "\n")
cat("臂覆盖: 每切片使用", unique(read_log$n_arm_used), "条臂 / 全局定义", length(ARMS), "条\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 合并 clone/坐标信息并做数据完整性校验（行数、重复 cell、切片一致性、缺失值）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
cm <- do.call(rbind, lapply(PATIENTS, function(p)
  read.csv(file.path(BASE, "patient_clone", p, "clone_membership.csv"), stringsAsFactors = FALSE)))
cm$cell <- as.character(cm$cell)
m <- match(spot$cell, cm$cell)
spot$x     <- cm$x[m]
spot$y     <- cm$y[m]
spot$clone <- cm$clone[m]

slice_cmp <- data.frame(
  slice = SLICES,
  n_spot_table = as.integer(table(factor(spot$slice, levels = SLICES))),
  n_membership = as.integer(table(factor(cm$slice, levels = SLICES))),
  stringsAsFactors = FALSE)
slice_cmp$match <- slice_cmp$n_spot_table == slice_cmp$n_membership

checks <- data.frame(
  check = c("total spot rows", "membership rows", "per-slice row counts all equal",
            "duplicated cell in spot table", "duplicated cell in membership",
            "unmatched cells (no clone/x/y)", "slice label mismatch",
            "NA x or y", "NA clone", "all zones non-empty in every slice"),
  value = c(nrow(spot), nrow(cm), sum(!slice_cmp$match),
            sum(duplicated(spot$cell)), sum(duplicated(cm$cell)),
            sum(is.na(m)), sum(!is.na(m) & cm$slice[m] != spot$slice),
            sum(is.na(spot$x) | is.na(spot$y)), sum(is.na(spot$clone)),
            sum(slice_cmp$n_spot_table == 0 | !complete.cases(zone_by_slice[, c("Stroma","Interface","Core")]) |
                zone_by_slice$Stroma == 0 | zone_by_slice$Interface == 0 | zone_by_slice$Core == 0)),
  expected = c(nrow(cm), nrow(cm), 0, 0, 0, 0, 0, 0, 0, 0),
  stringsAsFactors = FALSE)
checks$ok <- checks$value == checks$expected

print(slice_cmp, row.names = FALSE)
cat("\n")
print(checks, row.names = FALSE)
cat("\n全部完整性检查通过:", all(checks$ok), "\n")
cat("每患者 clone 数:", paste(sprintf("%s=%d", PATIENTS,
    vapply(PATIENTS, function(p) length(unique(cm$clone[cm$slice %in% SLICE_EXPECT[[p]]])), 0L)),
    collapse = " | "), "\n")
write.csv(checks, file.path(REVIEW, "qc_integrity_checks.csv"), row.names = FALSE)
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# QC：三区计数、分数分布、与既有 cnv_burden 的 Spearman、空间邻域一致性、臂覆盖
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
spot$cnv_burden <- cm$cnv_burden[m]

# 1) 分数分布（原始 + 归一化）
dist_qc <- do.call(rbind, lapply(SLICES, function(sl) {
  i <- spot$slice == sl
  v <- spot$cnv_arm_score[i]; nv <- spot$cnv_norm[i]
  q <- quantile(v, c(0, .05, .25, .5, .75, .95, 1), names = FALSE)
  qn <- quantile(nv, c(0, .05, .25, .5, .75, .95, 1), names = FALSE)
  data.frame(slice = sl, n = length(v),
             raw_min = q[1], raw_q05 = q[2], raw_median = q[4], raw_q95 = q[6], raw_max = q[7],
             raw_mean = mean(v), raw_sd = sd(v), n_unique_raw = length(unique(round(v, 8))),
             norm_q05 = qn[2], norm_median = qn[4], norm_q95 = qn[6],
             stringsAsFactors = FALSE)
}))

# 2) 与 clone_membership 的 cnv_burden 的方向一致性（诊断，不作结论）
cor_qc <- do.call(rbind, lapply(SLICES, function(sl) {
  i <- spot$slice == sl
  ct <- suppressWarnings(cor.test(spot$cnv_arm_score[i], spot$cnv_burden[i], method = "spearman"))
  data.frame(slice = sl, n = sum(i), spearman_rho = unname(ct$estimate),
             p_value = ct$p.value, stringsAsFactors = FALSE)
}))

# 3) 空间邻域一致性：每 spot 的 6 个最近邻中同区的比例（诊断是否椒盐状）
sp_qc <- list()
for (sl in SLICES) {
  d <- spot[spot$slice == sl & !is.na(spot$zone), ]
  z <- as.character(d$zone)
  dm <- as.matrix(dist(as.matrix(d[, c("x", "y")]))); diag(dm) <- Inf
  k <- min(6L, nrow(d) - 1L)
  nn <- t(apply(dm, 1, function(v) order(v)[seq_len(k)]))
  same <- mean(z[nn] == z)
  p <- as.numeric(table(z)) / length(z)
  sp_qc[[sl]] <- data.frame(slice = sl, n = nrow(d), k = k,
                            same_zone_frac = same, expected_random = sum(p^2),
                            enrichment = same / sum(p^2), stringsAsFactors = FALSE)
  rm(dm); invisible(gc(verbose = FALSE))
}
sp_qc <- do.call(rbind, sp_qc); rownames(sp_qc) <- NULL

# 4) 臂覆盖：逐切片记录使用到的臂数与缺失臂
arm_qc <- list()
for (sl in SLICES) {
  o <- readRDS(file.path(SAMPLE_ROOT, sl, "run.final.infercnv_obj"))
  g <- rownames(o@expr.data)
  ga <- go$arm[match(g, go$gene)]
  tb <- table(factor(ga[!is.na(ga)], levels = ARMS))
  arm_qc[[sl]] <- data.frame(slice = sl, n_arm_used = sum(tb > 0),
                             missing_arms = paste(names(tb)[tb == 0], collapse = ","),
                             stringsAsFactors = FALSE)
  rm(o); invisible(gc(verbose = FALSE))
}
arm_qc <- do.call(rbind, arm_qc); rownames(arm_qc) <- NULL

cat("--- 分数分布 ---\n"); print(dist_qc, row.names = FALSE, digits = 5)
cat("\n--- 与 cnv_burden 的 Spearman（诊断） ---\n"); print(cor_qc, row.names = FALSE, digits = 4)
cat("\n--- 空间邻域一致性（6-NN 同区比例 vs 随机期望） ---\n"); print(sp_qc, row.names = FALSE, digits = 4)
cat("\n--- 臂覆盖 ---\n"); print(arm_qc, row.names = FALSE)

# 落盘内部证据（不随交付）
write.csv(read_log, file.path(REVIEW, "qc_read_log.csv"), row.names = FALSE)
write.csv(zone_by_slice, file.path(REVIEW, "qc_zone_counts_by_slice.csv"), row.names = FALSE)
write.csv(dist_qc, file.path(REVIEW, "qc_score_distribution_by_slice.csv"), row.names = FALSE)
write.csv(cor_qc, file.path(REVIEW, "qc_consistency_vs_cnv_burden.csv"), row.names = FALSE)
write.csv(sp_qc, file.path(REVIEW, "qc_spatial_concordance_by_slice.csv"), row.names = FALSE)
write.csv(arm_qc, file.path(REVIEW, "qc_arm_coverage_by_slice.csv"), row.names = FALSE)
write.csv(checks, file.path(REVIEW, "qc_integrity_checks.csv"), row.names = FALSE)
write.csv(slice_cmp, file.path(REVIEW, "qc_slice_row_counts.csv"), row.names = FALSE)
cat("\n内部 QC 证据已写入:", REVIEW, "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 落盘逐 spot 赋值表（CSV + qs），并回读校验行数/列名/细胞顺序
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
spot_out <- spot[, c("patient", "slice", "cell", "x", "y", "cnv_arm_score", "cnv_norm", "zone", "clone")]
spot_out <- spot_out[order(match(spot_out$patient, PATIENTS),
                           match(spot_out$slice, SLICES), spot_out$cell), , drop = FALSE]
spot_out$zone <- as.character(spot_out$zone)
spot_out$cnv_arm_score <- round(spot_out$cnv_arm_score, 6)
spot_out$cnv_norm      <- round(spot_out$cnv_norm, 6)
rownames(spot_out) <- NULL

csv_path <- file.path(OUT, "tables/zone_spot_assignment.csv")
qs_path  <- file.path(OUT, "data/zone_spot_assignment.qs")
write.csv(spot_out, csv_path, row.names = FALSE)
save_obj(spot_out, qs_path)

back <- read.csv(csv_path, stringsAsFactors = FALSE)
cat("\n回读校验: 行数", nrow(back), "(期望", nrow(spot_out), ") | 列数", ncol(back),
    "| 列名一致:", identical(colnames(back), colnames(spot_out)),
    "| 细胞顺序一致:", identical(as.character(back$cell), as.character(spot_out$cell)), "\n")
cat("CSV 大小:", round(file.info(csv_path)$size / 1e6, 2), "MB\n")

pct_by_slice <- sweep(as.matrix(zone_by_slice[, c("Stroma", "Interface", "Core")]), 1,
                      zone_by_slice$n_spot, "/") * 100
pct_by_slice <- data.frame(slice = zone_by_slice$slice, round(pct_by_slice, 2), check.names = FALSE)
cat("\n--- 每切片三区 spot 占比 (%) ---\n"); print(pct_by_slice, row.names = FALSE)
cat("\n--- 合计 ---\n")
cat(sprintf("Stroma %d (%.2f%%) | Interface %d (%.2f%%) | Core %d (%.2f%%) | 总 %d\n",
            sum(spot_out$zone == "Stroma"), 100 * mean(spot_out$zone == "Stroma"),
            sum(spot_out$zone == "Interface"), 100 * mean(spot_out$zone == "Interface"),
            sum(spot_out$zone == "Core"), 100 * mean(spot_out$zone == "Core"), nrow(spot_out)))
cat("输出:", csv_path, "\n对象:", qs_path, "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 读取 U1 逐 spot 赋值表并校验完整性；定义三区定性配色（nature 主题）与作图参数
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(png)})
ZONE_LEVELS <- c("Stroma", "Interface", "Core")
# nature 主题内的语义顺序配色：低 CNV=深蓝，中=salmon，高 CNV=红
ZONE_COL <- c(Stroma = "#3C5488", Interface = "#F39B7F", Core = "#DC0000")
FIG_DIR <- file.path(OUT, "figures"); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

spot2 <- read.csv(file.path(OUT, "tables/zone_spot_assignment.csv"), stringsAsFactors = FALSE)
spot2$zone <- factor(spot2$zone, levels = ZONE_LEVELS)

chk2 <- data.frame(
  check = c("rows", "duplicated cell", "NA zone", "NA x/y", "NA clone", "n slices", "n patients"),
  value = c(nrow(spot2), sum(duplicated(spot2$cell)), sum(is.na(spot2$zone)),
            sum(is.na(spot2$x) | is.na(spot2$y)), sum(is.na(spot2$clone)),
            length(unique(spot2$slice)), length(unique(spot2$patient))),
  expected = c(10583, 0, 0, 0, 0, 11, 6), stringsAsFactors = FALSE)
chk2$ok <- chk2$value == chk2$expected
print(chk2, row.names = FALSE)
cat("\n逐切片行数:\n"); print(as.data.frame(table(spot2$slice)), row.names = FALSE)
cat("\n三区合计:", paste(sprintf("%s=%d", ZONE_LEVELS, as.integer(table(spot2$zone))), collapse = " | "), "\n")
cat("x 范围:", paste(round(range(spot2$x), 1), collapse = " .. "),
    "| y 范围:", paste(round(range(spot2$y), 1), collapse = " .. "), "\n")
stopifnot(all(chk2$ok))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 生成三区占比（切片/病人两层）与 clone×三区双向百分比表，并落盘
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
PATIENT_ORDER <- c("P9", "P11", "P12", "P13", "P14", "P15")
SLICE_ORDER   <- c("092A", "092B", "093B", "093C", "093D", "094A", "094B", "094C", "094D", "095A", "095B")

# 1) 切片层三区计数与占比
slice_cnt <- spot2 %>%
  count(patient, slice, zone, .drop = FALSE) %>%
  tidyr::complete(tidyr::nesting(patient, slice), zone = ZONE_LEVELS, fill = list(n = 0))
slice_tab <- slice_cnt %>%
  tidyr::pivot_wider(names_from = zone, values_from = n, values_fill = 0) %>%
  mutate(n_spot = Stroma + Interface + Core,
         Stroma_pct = 100 * Stroma / n_spot,
         Interface_pct = 100 * Interface / n_spot,
         Core_pct = 100 * Core / n_spot) %>%
  rename(Stroma_n = Stroma, Interface_n = Interface, Core_n = Core) %>%
  arrange(match(patient, PATIENT_ORDER), match(slice, SLICE_ORDER))

# 2) 病人层三区计数与占比（n_slices = 该患者实际纳入的切片数，distinct slice）
n_slice_by_patient <- table(unique(spot2[, c("patient", "slice")])$patient)
patient_tab <- spot2 %>%
  count(patient, zone, .drop = FALSE) %>%
  tidyr::complete(patient, zone = ZONE_LEVELS, fill = list(n = 0)) %>%
  tidyr::pivot_wider(names_from = zone, values_from = n, values_fill = 0) %>%
  mutate(n_spot = Stroma + Interface + Core,
         Stroma_pct = 100 * Stroma / n_spot,
         Interface_pct = 100 * Interface / n_spot,
         Core_pct = 100 * Core / n_spot,
         n_slices = as.integer(n_slice_by_patient[patient])) %>%
  rename(Stroma_n = Stroma, Interface_n = Interface, Core_n = Core) %>%
  select(patient, n_slices, n_spot, Stroma_n, Interface_n, Core_n,
         Stroma_pct, Interface_pct, Core_pct) %>%
  arrange(match(patient, PATIENT_ORDER))

# 3) clone × 三区计数（长表；仅纳入该切片实际出现的 clone，三区补零）
clone_keys <- spot2 %>% distinct(patient, slice, clone) %>%
  arrange(match(patient, PATIENT_ORDER), match(slice, SLICE_ORDER), clone)
clone_cnt <- clone_keys %>%
  tidyr::expand(tidyr::nesting(patient, slice, clone), zone = ZONE_LEVELS) %>%
  left_join(spot2 %>% count(patient, slice, clone, zone, name = "n_spot"),
            by = c("patient", "slice", "clone", "zone")) %>%
  mutate(n_spot = ifelse(is.na(n_spot), 0L, as.integer(n_spot))) %>%
  arrange(match(patient, PATIENT_ORDER), match(slice, SLICE_ORDER), clone,
          match(zone, ZONE_LEVELS))

# 4) 双向百分比：clone→区（每 clone 三区行和为 100）与 区→clone（每区各 clone 列和为 100）
pct_c2z <- clone_cnt %>% group_by(patient, slice, clone) %>%
  mutate(pct = 100 * n_spot / sum(n_spot)) %>% ungroup() %>% mutate(pct_type = "clone_to_zone")
pct_z2c <- clone_cnt %>% group_by(patient, slice, zone) %>%
  mutate(pct = 100 * n_spot / sum(n_spot)) %>% ungroup() %>% mutate(pct_type = "zone_to_clone")
clone_pct <- bind_rows(pct_c2z, pct_z2c) %>%
  select(patient, slice, clone, zone, pct_type, pct) %>%
  arrange(match(patient, PATIENT_ORDER), match(slice, SLICE_ORDER), pct_type,
          clone, match(zone, ZONE_LEVELS)) %>%
  mutate(pct = round(pct, 4))

# 5) 落盘
write.csv(slice_tab,   file.path(OUT, "tables/zone_proportion_by_slice.csv"),   row.names = FALSE)
write.csv(patient_tab, file.path(OUT, "tables/zone_proportion_by_patient.csv"), row.names = FALSE)
write.csv(clone_cnt,   file.path(OUT, "tables/clone_zone_counts.csv"),          row.names = FALSE)
write.csv(clone_pct,   file.path(OUT, "tables/clone_zone_pct.csv"),             row.names = FALSE)

cat("表行数: slice =", nrow(slice_tab), "| patient =", nrow(patient_tab),
    "| clone counts =", nrow(clone_cnt), "| clone pct =", nrow(clone_pct),
    "| clone-slice 组合 =", nrow(clone_keys), "\n")
cat("\n--- 切片层三区占比 ---\n")
print(slice_tab %>% mutate(across(ends_with("_pct"), ~round(.x, 2))), row.names = FALSE)
cat("\n--- 病人层三区占比 ---\n")
print(patient_tab %>% mutate(across(ends_with("_pct"), ~round(.x, 2))), row.names = FALSE)
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 从磁盘回读四张统计表，校验行数、百分比闭合、计数与百分比互相可反推、n_slices 正确
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
s_tab  <- read.csv(file.path(OUT, "tables/zone_proportion_by_slice.csv"), stringsAsFactors = FALSE)
p_tab  <- read.csv(file.path(OUT, "tables/zone_proportion_by_patient.csv"), stringsAsFactors = FALSE)
cc_tab <- read.csv(file.path(OUT, "tables/clone_zone_counts.csv"), stringsAsFactors = FALSE)
cp_tab <- read.csv(file.path(OUT, "tables/clone_zone_pct.csv"), stringsAsFactors = FALSE)

c2z <- cp_tab %>% filter(pct_type == "clone_to_zone") %>%
  group_by(patient, slice, clone) %>% summarise(s = sum(pct), .groups = "drop")
z2c <- cp_tab %>% filter(pct_type == "zone_to_clone") %>%
  group_by(patient, slice, zone) %>% summarise(s = sum(pct), .groups = "drop")
cc_slice <- cc_tab %>% group_by(patient, slice) %>% summarise(n = sum(n_spot), .groups = "drop")
cc_pat   <- cc_tab %>% group_by(patient) %>% summarise(n = sum(n_spot), .groups = "drop")

re_c2z <- cc_tab %>% group_by(patient, slice, clone) %>%
  mutate(pct = round(100 * n_spot / sum(n_spot), 4)) %>% ungroup() %>%
  select(patient, slice, clone, zone, pct) %>% mutate(pct_type = "clone_to_zone")
re_z2c <- cc_tab %>% group_by(patient, slice, zone) %>%
  mutate(pct = round(100 * n_spot / sum(n_spot), 4)) %>% ungroup() %>%
  select(patient, slice, clone, zone, pct) %>% mutate(pct_type = "zone_to_clone")
cmp <- bind_rows(re_c2z, re_z2c) %>%
  inner_join(cp_tab, by = c("patient","slice","clone","zone","pct_type"), suffix = c("_re","_tab"))

expected_n_slices <- as.integer(table(s_tab$patient)[p_tab$patient])
n_slices_ok <- all(p_tab$n_slices == expected_n_slices)

checks2 <- data.frame(
  check = c("slice table rows", "patient table rows", "clone counts rows", "clone pct rows",
            "patient n_slices = distinct slice count",
            "slice pct max |rowsum-100|", "patient pct max |rowsum-100|",
            "clone_to_zone max |sum-100|", "zone_to_clone max |sum-100|",
            "counts total = 10583", "counts per slice = slice table n_spot",
            "counts per patient = patient table n_spot", "counts vs pct max abs diff"),
  value = c(nrow(s_tab), nrow(p_tab), nrow(cc_tab), nrow(cp_tab),
            as.integer(n_slices_ok),
            max(abs(rowSums(s_tab[, c("Stroma_pct","Interface_pct","Core_pct")]) - 100)),
            max(abs(rowSums(p_tab[, c("Stroma_pct","Interface_pct","Core_pct")]) - 100)),
            max(abs(c2z$s - 100)), max(abs(z2c$s - 100)),
            sum(cc_tab$n_spot),
            max(abs(cc_slice$n - s_tab$n_spot[match(paste(cc_slice$patient, cc_slice$slice),
                                                    paste(s_tab$patient, s_tab$slice))])),
            max(abs(cc_pat$n - p_tab$n_spot[match(cc_pat$patient, p_tab$patient)])),
            max(abs(cmp$pct_re - cmp$pct_tab))),
  expected = c(11, 6, 78, 156, 1, 0, 0, 0, 0, 10583, 0, 0, 0),
  tolerance = c(0, 0, 0, 0, 0, 1e-2, 1e-2, 1e-2, 1e-2, 0, 0, 0, 1e-4),
  stringsAsFactors = FALSE)
checks2$ok <- abs(checks2$value - checks2$expected) <= checks2$tolerance
print(checks2, row.names = FALSE, digits = 6)
cat("\n表间自洽全部通过:", all(checks2$ok), "\n")
cat("patient n_slices:", paste(sprintf("%s=%d", p_tab$patient, p_tab$n_slices), collapse = " | "), "\n")
write.csv(checks2, file.path(REVIEW, "qc_table_consistency.csv"), row.names = FALSE)
stopifnot(all(checks2$ok))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 统一重绘 11 张空间三区图：双行标题、坐标轴标签、固定画幅、1 mm 比例尺、颜色+形状双线索
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
LIM <- c(58, 527)   # 统一到 Visium 捕获区域（约 6.5 mm 视野），并留出点半径边距
SHAPES <- c(Stroma = 16, Interface = 17, Core = 15)   # 圆 / 三角 / 方：灰度下的第二线索
fig_log <- list()
for (sl in SLICE_ORDER) {
  d   <- spot2[spot2$slice == sl, ]
  pid <- unique(d$patient)
  if (length(pid) != 1L) stop("slice maps to multiple patients: ", sl)
  pct <- 100 * as.numeric(table(d$zone)) / nrow(d)
  dm  <- as.matrix(dist(d[, c("x", "y")])); diag(dm) <- Inf
  nn_med  <- median(apply(dm, 1, min))       # 最近邻中位距离 ~ 100 µm（Visium spot pitch）
  bar_px  <- 10 * nn_med                     # 1 mm 对应的像素长度
  x0 <- LIM[1] + 0.04 * diff(LIM); y0 <- LIM[1] + 0.07 * diff(LIM)

  p <- ggplot(d, aes(x = x, y = y, colour = zone, shape = zone)) +
    geom_point(size = 1.1, stroke = 0) +
    annotate("segment", x = x0, xend = x0 + bar_px, y = y0, yend = y0,
             linewidth = 1.4, colour = "black") +
    annotate("text", x = x0 + bar_px / 2, y = y0 + 6, label = "1 mm",
             size = 2.9, vjust = 0) +
    scale_colour_manual(values = ZONE_COL, breaks = ZONE_LEVELS, drop = FALSE, name = "Zone") +
    scale_shape_manual(values = SHAPES, breaks = ZONE_LEVELS, drop = FALSE, name = "Zone") +
    guides(colour = guide_legend(override.aes = list(size = 3)),
           shape  = guide_legend(override.aes = list(size = 3))) +
    coord_fixed(xlim = LIM, ylim = LIM, expand = FALSE) +
    labs(title = sprintf("%s (%s) - spatial zones", sl, pid),
         subtitle = sprintf("Stroma %.1f%% | Interface %.1f%% | Core %.1f%%", pct[1], pct[2], pct[3]),
         x = "Spatial x (px)", y = "Spatial y (px)") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9.5),
          axis.title = element_text(size = 8.5),
          axis.text = element_text(size = 6.5),
          axis.ticks = element_line(linewidth = 0.25),
          legend.position = "right", legend.key = element_rect(fill = "white"),
          legend.title = element_text(face = "bold", size = 9.5),
          legend.text = element_text(size = 9),
          plot.margin = margin(8, 10, 8, 8))

  stem <- file.path(FIG_DIR, sprintf("zone_spatial_%s_%s", pid, sl))
  ggsave(paste0(stem, ".png"), p, width = 7.0, height = 6.5, dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = 7.0, height = 6.5, bg = "white")
  fig_log[[sl]] <- data.frame(patient = pid, slice = sl, n_spot = nrow(d),
                              Stroma_pct = pct[1], Interface_pct = pct[2], Core_pct = pct[3],
                              nn_med_px = nn_med, scale_bar_px = bar_px,
                              stem = stem, stringsAsFactors = FALSE)
  cat(sprintf("%s: n=%d | %.1f/%.1f/%.1f%% | 1mm=%.1f px | png %.0f KB\n", sl, nrow(d),
              pct[1], pct[2], pct[3], bar_px, file.info(paste0(stem, ".png"))$size / 1024))
}
fig_log <- do.call(rbind, fig_log); rownames(fig_log) <- NULL
write.csv(fig_log, file.path(REVIEW, "qc_figure_manifest.csv"), row.names = FALSE)
cat("\n图数量: png =", length(list.files(FIG_DIR, pattern = "\\.png$")),
    "| pdf =", length(list.files(FIG_DIR, pattern = "\\.pdf$")), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 像素完整性自检（最终版）：11 张 PNG 逐张检查尺寸、标准差、非白像素占比、文件大小与 pdf/png 配对
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
px <- list()
for (sl in SLICE_ORDER) {
  pid <- fig_log$patient[fig_log$slice == sl]
  f   <- file.path(FIG_DIR, sprintf("zone_spatial_%s_%s.png", pid, sl))
  img <- png::readPNG(f)
  rgb <- img[, , 1:3, drop = FALSE]
  sdv <- sd(as.vector(rgb))
  ink <- mean(rgb < 0.95)
  sz  <- file.info(f)$size
  px[[sl]] <- data.frame(
    slice = sl, file = basename(f), width = dim(rgb)[2], height = dim(rgb)[1],
    size_kb = round(sz / 1024, 1), sd = round(sdv, 4), frac_ink = round(ink, 5),
    pdf_exists = file.exists(sub("\\.png$", ".pdf", f)),
    ok = (sdv > 0.02) && (ink > 0.001) && (sz > 10240) && file.exists(sub("\\.png$", ".pdf", f)),
    stringsAsFactors = FALSE)
}
px <- do.call(rbind, px); rownames(px) <- NULL
print(px, row.names = FALSE)
cat("\n像素自检通过:", sum(px$ok), "/", nrow(px), "\n")
cat("检查项: sd>0.02, 非白像素>0.1%, PNG>10KB, PDF 成对\n")
write.csv(px, file.path(REVIEW, "qc_figure_pixel_check.csv"), row.names = FALSE)
if (!all(px$ok)) stop("pixel integrity check failed; redraw required")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 将汇总 Markdown 转为单文件自包含 HTML（图片 base64 内嵌），并校验产物
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
md_path   <- file.path(OUT, "step1_three_zone_汇总.md")
html_path <- file.path(OUT, "step1_three_zone_汇总.html")
script    <- path.expand("~/.claude/scripts/report_md_to_html.py")
res <- system2("python3", c(shQuote(script), shQuote(md_path), shQuote(html_path)),
               stdout = TRUE, stderr = TRUE)
cat(res, sep = "\n")
cat("\nHTML exists:", file.exists(html_path),
    "| size MB:", round(file.info(html_path)$size / 1e6, 2), "\n")
html_txt <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
cat("embedded PNG count:", lengths(regmatches(html_txt, gregexpr("data:image/png;base64", html_txt))), "\n")
cat("contains review token:", grepl("review/", html_txt, fixed = TRUE), "\n")
cat("contains absolute home path:", grepl("/path/to", html_txt, fixed = TRUE), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# Unit 3 最终文字修订后：文档-磁盘一致性对账 + 自包含/禁词复核
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
md_path <- file.path(OUT, "step1_three_zone_汇总.md")
cat("md exists:", file.exists(md_path), "\n")
script2 <- path.expand("~/.claude/scripts/doc_disk_reconcile.py")
res2 <- system2("python3", c(shQuote(script2), "--all", shQuote(OUT)), stdout = TRUE, stderr = TRUE)
cat(res2, sep = "\n")
st <- attr(res2, "status"); if (is.null(st)) st <- 0L
cat("\nreconcile status:", st, "\n")
html_path <- file.path(OUT, "step1_three_zone_汇总.html")
html_txt <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
cat("HTML MB:", round(file.info(html_path)$size / 1e6, 2),
    "| base64 PNG:", lengths(regmatches(html_txt, gregexpr("data:image/png;base64", html_txt))),
    "| review token:", grepl("review/", html_txt, fixed = TRUE),
    "| absolute home:", grepl("/path/to", html_txt, fixed = TRUE), "\n")
cat("old 'stacked bar not generated' sentence:",
    grepl("堆叠柱状图、clone×区热图、边界距离梯度曲线均未生成", html_txt, fixed = TRUE), "\n")
cat("new pointer present:", grepl("见下方「补充图」小节", html_txt, fixed = TRUE), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# Unit 3 初始化：读取 spot 表、定义三区/clone 配色、记录既有交付物基线（md5 + mtime）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(png)})
cat("base64enc available:", requireNamespace("base64enc", quietly = TRUE), "\n")

OUT     <- "/path/to/project/project_spatial/TNBC/TNBC_analysis/SpatialInferCNV/results/step1_three_zone"
REVIEW  <- "/path/to/project/project_spatial/TNBC/TNBC_analysis/SpatialInferCNV/review/step1_three_zone"
FIG_DIR <- file.path(OUT, "figures")
ZONE_LEVELS <- c("Stroma", "Interface", "Core")
ZONE_COL    <- c(Stroma = "#3C5488", Interface = "#F39B7F", Core = "#DC0000")
NATURE8     <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000")

spot3 <- read.csv(file.path(OUT, "tables/zone_spot_assignment.csv"), stringsAsFactors = FALSE)
spot3$zone <- factor(spot3$zone, levels = ZONE_LEVELS)
SLICE_ORDER <- c("092A", "092B", "093B", "093C", "093D", "094A", "094B", "094C", "094D", "095A", "095B")
sp_map <- unique(spot3[, c("slice", "patient")])
SLICE2PAT <- setNames(sp_map$patient, sp_map$slice)
clone_levels_by_patient <- lapply(split(spot3$clone, spot3$patient), function(x) sort(unique(x)))
CLONE_COL <- lapply(clone_levels_by_patient, function(lv) setNames(NATURE8[seq_along(lv)], lv))

# 既有交付物基线：4 张表 + qs + 11 组旧图（pdf/png）
old_tables <- file.path(OUT, c("tables/zone_spot_assignment.csv", "tables/zone_proportion_by_slice.csv",
                               "tables/zone_proportion_by_patient.csv", "tables/clone_zone_counts.csv",
                               "tables/clone_zone_pct.csv", "data/zone_spot_assignment.qs"))
old_figs <- unlist(lapply(SLICE_ORDER, function(sl) {
  pid <- SLICE2PAT[[sl]]
  file.path(FIG_DIR, sprintf("zone_spatial_%s_%s.%s", pid, sl, c("pdf", "png")))
}), use.names = FALSE)
old_files <- c(old_tables, old_figs)
base <- data.frame(file = old_files, size = file.info(old_files)$size,
                   mtime = as.character(file.info(old_files)$mtime),
                   md5 = as.character(tools::md5sum(old_files)), stringsAsFactors = FALSE)
write.csv(base, file.path(REVIEW, "qc_old_artifacts_baseline_u3.csv"), row.names = FALSE)

cat("spot rows:", nrow(spot3), "| duplicated cell:", sum(duplicated(spot3$cell)),
    "| zones:", paste(levels(spot3$zone), collapse = "/"), "\n")
cat("baseline files:", nrow(base), "| all exist:", all(file.exists(old_files)), "\n")
cat("clone levels by patient:\n"); print(clone_levels_by_patient)
cat("clone colors by patient:\n"); print(CLONE_COL)

html_path <- file.path(OUT, "step1_three_zone_汇总.html")
html <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
cat("HTML chars:", nchar(html),
    "| tech-spec h2 pos:", regexpr('<h2 id="[^"]+">技术说明书</h2>', html)[1],
    "| TOC nav pos:", regexpr('<nav class="toc" id="toc">', html)[1], "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 图 1/4：空间散点按 cnv_norm 连续着色（viridis + colorbar，不二值化）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
LIM <- c(58, 527)
FIG1_LOG <- list()
for (sl in SLICE_ORDER) {
  d   <- spot3[spot3$slice == sl, ]
  pid <- unique(d$patient)
  dm  <- as.matrix(dist(d[, c("x", "y")])); diag(dm) <- Inf
  bar_px <- 10 * median(apply(dm, 1, min))
  x0 <- LIM[1] + 0.04 * diff(LIM); y0 <- LIM[1] + 0.07 * diff(LIM)

  p <- ggplot(d, aes(x = x, y = y, colour = cnv_norm)) +
    geom_point(size = 1.1, stroke = 0) +
    annotate("segment", x = x0, xend = x0 + bar_px, y = y0, yend = y0,
             linewidth = 1.4, colour = "black") +
    annotate("text", x = x0 + bar_px / 2, y = y0 + 6, label = "1 mm",
             size = 2.9, vjust = 0) +
    scale_colour_viridis_c(option = "viridis", name = "cnv_norm",
                           limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    coord_fixed(xlim = LIM, ylim = LIM, expand = FALSE) +
    labs(title = sprintf("%s (%s) - cnv_norm (within-slice)", sl, pid),
         subtitle = "cnv_norm = within-slice min-max; not absolute tumor content",
         x = "Spatial x (px)", y = "Spatial y (px)") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8.5),
          axis.title = element_text(size = 8.5), axis.text = element_text(size = 6.5),
          axis.ticks = element_line(linewidth = 0.25),
          legend.position = "right", legend.title = element_text(face = "bold", size = 9),
          legend.text = element_text(size = 8),
          plot.margin = margin(8, 10, 8, 8))

  stem <- file.path(FIG_DIR, sprintf("zone_cnv_spatial_%s_%s", pid, sl))
  ggsave(paste0(stem, ".png"), p, width = 7.0, height = 6.5, dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = 7.0, height = 6.5, bg = "white")
  FIG1_LOG[[sl]] <- data.frame(patient = pid, slice = sl, n_spot = nrow(d), stem = stem,
                               stringsAsFactors = FALSE)
  cat(sprintf("%s: n=%d | png %.0f KB\n", sl, nrow(d), file.info(paste0(stem, ".png"))$size / 1024))
}
FIG1_LOG <- do.call(rbind, FIG1_LOG); rownames(FIG1_LOG) <- NULL
cat("type1 count:", nrow(FIG1_LOG), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 图 2/4：三区 cnv_norm 箱线图 + 单 spot 横向抖动点（低饱和，箱体置顶）+ n/中位数
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
FIG2_LOG <- list()
for (sl in SLICE_ORDER) {
  d   <- spot3[spot3$slice == sl, ]
  pid <- unique(d$patient)
  lab <- d %>% group_by(zone) %>% summarise(n = n(), med = median(cnv_norm), .groups = "drop") %>%
    mutate(label = sprintf("n=%d\nmedian=%.3f", n, med), y = 1.06)

  p <- ggplot(d, aes(x = zone, y = cnv_norm)) +
    geom_jitter(aes(colour = zone), width = 0.16, height = 0, size = 0.25,
                alpha = 0.18, show.legend = FALSE) +
    geom_boxplot(aes(fill = zone), width = 0.6, outlier.shape = NA, alpha = 0.80,
                 colour = "black", linewidth = 0.4) +
    geom_text(data = lab, aes(x = zone, y = y, label = label), inherit.aes = FALSE,
              size = 2.7, vjust = 0, lineheight = 0.9) +
    scale_fill_manual(values = ZONE_COL, breaks = ZONE_LEVELS, drop = FALSE, name = "Zone") +
    scale_colour_manual(values = ZONE_COL, breaks = ZONE_LEVELS, drop = FALSE, guide = "none") +
    scale_y_continuous(limits = c(0, 1.18), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    labs(title = sprintf("%s (%s) - cnv_norm by zone", sl, pid),
         subtitle = "box = IQR drawn on top; whiskers = 1.5x IQR; points = individual spots (horizontal jitter only)",
         x = "Zone", y = "cnv_norm (within-slice)") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title = element_text(size = 9), axis.text = element_text(size = 8),
          axis.ticks = element_line(linewidth = 0.25),
          legend.position = "right", legend.title = element_text(face = "bold", size = 9.5),
          legend.text = element_text(size = 9),
          plot.margin = margin(8, 10, 8, 8))

  stem <- file.path(FIG_DIR, sprintf("zone_cnv_box_%s_%s", pid, sl))
  ggsave(paste0(stem, ".png"), p, width = 6.5, height = 6.0, dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = 6.5, height = 6.0, bg = "white")
  FIG2_LOG[[sl]] <- data.frame(patient = pid, slice = sl, n_spot = nrow(d), stem = stem,
                               stringsAsFactors = FALSE)
  cat(sprintf("%s: n=%d | medians S/I/C = %.3f/%.3f/%.3f | png %.0f KB\n",
              sl, nrow(d), lab$med[1], lab$med[2], lab$med[3],
              file.info(paste0(stem, ".png"))$size / 1024))
}
FIG2_LOG <- do.call(rbind, FIG2_LOG); rownames(FIG2_LOG) <- NULL
cat("type2 count:", nrow(FIG2_LOG), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 图 3/4：空间散点按 clone 定性着色（同一病人内跨切片颜色一致；图例含 spot 数）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
LIM <- c(58, 527)
FIG3_LOG <- list()
for (sl in SLICE_ORDER) {
  d   <- spot3[spot3$slice == sl, ]
  pid <- unique(d$patient)
  cc  <- CLONE_COL[[pid]]
  present <- names(cc)[names(cc) %in% unique(d$clone)]
  cnt <- d %>% count(clone) %>% mutate(label = sprintf("%s (n=%d)", clone, n))
  labels <- setNames(cnt$label, cnt$clone)
  dm  <- as.matrix(dist(d[, c("x", "y")])); diag(dm) <- Inf
  bar_px <- 10 * median(apply(dm, 1, min))
  x0 <- LIM[1] + 0.04 * diff(LIM); y0 <- LIM[1] + 0.07 * diff(LIM)

  p <- ggplot(d, aes(x = x, y = y, colour = clone)) +
    geom_point(size = 1.1, stroke = 0) +
    annotate("segment", x = x0, xend = x0 + bar_px, y = y0, yend = y0,
             linewidth = 1.4, colour = "black") +
    annotate("text", x = x0 + bar_px / 2, y = y0 + 6, label = "1 mm",
             size = 2.9, vjust = 0) +
    scale_colour_manual(values = cc, breaks = present, labels = labels[present],
                        drop = FALSE, name = "Clone") +
    coord_fixed(xlim = LIM, ylim = LIM, expand = FALSE) +
    labs(title = sprintf("%s (%s) - clone identity", sl, pid),
         subtitle = "clone colors are patient-level; legend n = spots in this slice",
         x = "Spatial x (px)", y = "Spatial y (px)") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8.5),
          axis.title = element_text(size = 8.5), axis.text = element_text(size = 6.5),
          axis.ticks = element_line(linewidth = 0.25),
          legend.position = "right", legend.title = element_text(face = "bold", size = 9.5),
          legend.text = element_text(size = 8.5),
          plot.margin = margin(8, 10, 8, 8))

  stem <- file.path(FIG_DIR, sprintf("zone_clone_spatial_%s_%s", pid, sl))
  ggsave(paste0(stem, ".png"), p, width = 7.0, height = 6.5, dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = 7.0, height = 6.5, bg = "white")
  FIG3_LOG[[sl]] <- data.frame(patient = pid, slice = sl, n_spot = nrow(d), stem = stem,
                               stringsAsFactors = FALSE)
  cat(sprintf("%s: n=%d | clones: %s | png %.0f KB\n", sl, nrow(d),
              paste(sprintf("%s=%d", cnt$clone, cnt$n), collapse = ", "),
              file.info(paste0(stem, ".png"))$size / 1024))
}
FIG3_LOG <- do.call(rbind, FIG3_LOG); rownames(FIG3_LOG) <- NULL
cat("type3 count:", nrow(FIG3_LOG), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 图 4/4：三区 clone 组成堆叠柱（x=区、y=100%、fill=clone，柱顶标 spot 数）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
FIG4_LOG <- list()
for (sl in SLICE_ORDER) {
  d   <- spot3[spot3$slice == sl, ]
  pid <- unique(d$patient)
  cc  <- CLONE_COL[[pid]]
  df  <- d %>% count(zone, clone) %>% group_by(zone) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
  present <- names(cc)[names(cc) %in% unique(df$clone)]
  cnt <- d %>% count(clone) %>% mutate(label = sprintf("%s (n=%d)", clone, n))
  labels <- setNames(cnt$label, cnt$clone)
  tot <- d %>% count(zone) %>% mutate(label = sprintf("n=%d", n))

  p <- ggplot(df, aes(x = zone, y = pct, fill = clone)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
    geom_text(data = tot, aes(x = zone, y = 104, label = label), inherit.aes = FALSE,
              size = 3.0, vjust = 0) +
    scale_fill_manual(values = cc, breaks = present, labels = labels[present],
                      drop = FALSE, name = "Clone") +
    scale_y_continuous(limits = c(0, 114), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    labs(title = sprintf("%s (%s) - clone composition by zone", sl, pid),
         subtitle = "bar height = 100%; label above bar = spot count; exact pct in tables/clone_zone_pct.csv",
         x = "Zone", y = "Spot proportion (%)") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 7.8),
          axis.title = element_text(size = 9), axis.text = element_text(size = 8),
          axis.ticks = element_line(linewidth = 0.25),
          legend.position = "right", legend.title = element_text(face = "bold", size = 9.5),
          legend.text = element_text(size = 8.5),
          plot.margin = margin(8, 10, 8, 8))

  stem <- file.path(FIG_DIR, sprintf("zone_clone_bar_%s_%s", pid, sl))
  ggsave(paste0(stem, ".png"), p, width = 6.5, height = 6.0, dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = 6.5, height = 6.0, bg = "white")
  FIG4_LOG[[sl]] <- data.frame(patient = pid, slice = sl, n_spot = nrow(d), stem = stem,
                               stringsAsFactors = FALSE)
  cat(sprintf("%s: n=%d | clones=%d | png %.0f KB\n", sl, nrow(d), length(present),
              file.info(paste0(stem, ".png"))$size / 1024))
}
FIG4_LOG <- do.call(rbind, FIG4_LOG); rownames(FIG4_LOG) <- NULL
cat("type4 count:", nrow(FIG4_LOG), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 像素完整性自检 44/44；并核对既有 5 张表 + qs + 11 组旧图的 md5 与 mtime 均未变
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
TYPES <- c("zone_cnv_spatial", "zone_cnv_box", "zone_clone_spatial", "zone_clone_bar")
px <- list()
for (sl in SLICE_ORDER) {
  pid <- SLICE2PAT[[sl]]
  for (tp in TYPES) {
    stem <- file.path(FIG_DIR, sprintf("%s_%s_%s", tp, pid, sl))
    f    <- paste0(stem, ".png")
    img  <- png::readPNG(f); rgb <- img[, , 1:3, drop = FALSE]
    sdv  <- sd(as.vector(rgb)); ink <- mean(rgb < 0.95); sz <- file.info(f)$size
    px[[length(px) + 1]] <- data.frame(
      type = tp, slice = sl, file = basename(f),
      width = dim(rgb)[2], height = dim(rgb)[1], size_kb = round(sz / 1024, 1),
      sd = round(sdv, 4), frac_ink = round(ink, 5),
      pdf_exists = file.exists(paste0(stem, ".pdf")),
      ok = (sdv > 0.02) && (ink > 0.001) && (sz > 10240) && file.exists(paste0(stem, ".pdf")),
      stringsAsFactors = FALSE)
  }
}
px <- do.call(rbind, px); rownames(px) <- NULL
cat("figure count by type:\n"); print(table(px$type))
cat("pixel ok:", sum(px$ok), "/", nrow(px), "\n")
write.csv(px, file.path(REVIEW, "qc_extra_plots_pixel_check.csv"), row.names = FALSE)

now_info <- file.info(base$file); now_md5 <- as.character(tools::md5sum(base$file))
oldchk <- data.frame(file = base$file, md5_before = base$md5, md5_after = now_md5,
                     mtime_before = base$mtime, mtime_after = as.character(now_info$mtime),
                     unchanged = (base$md5 == now_md5) & (base$mtime == as.character(now_info$mtime)),
                     stringsAsFactors = FALSE)
write.csv(oldchk, file.path(REVIEW, "qc_old_artifacts_unchanged_u3.csv"), row.names = FALSE)
cat("old artifacts unchanged (md5 + mtime):", sum(oldchk$unchanged), "/", nrow(oldchk), "\n")
cat("mtime mismatches:", sum(oldchk$mtime_before != oldchk$mtime_after),
    "| md5 mismatches:", sum(oldchk$md5_before != oldchk$md5_after), "\n")
if (!all(px$ok)) stop("extra-plot pixel check not ok")
if (!all(oldchk$unchanged)) stop("old artifacts changed")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 将 Unit 3 的 4 类新图以 base64 内嵌方式增补进 HTML（保留既有内容），并加入 TOC 与图索引
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
html_path <- file.path(OUT, "step1_three_zone_汇总.html")
html <- paste(readLines(html_path, warn = FALSE), collapse = "\n")

img_tag <- function(stem, width = 880) {
  b64 <- gsub("[\r\n]", "", base64enc::base64encode(paste0(stem, ".png")))
  sprintf('<p style="text-align:center"><img src="data:image/png;base64,%s" width="%d"></p>', b64, width)
}
add_type <- function(sec, title, desc, prefix, width = 880) {
  sec <- c(sec, sprintf('<h3>%s</h3>', title), sprintf('<p>%s</p>', desc),
           sprintf('<p>文件：<code>figures/%s_&lt;P&gt;_&lt;slice&gt;.{pdf,png}</code></p>', prefix))
  for (sl in SLICE_ORDER) {
    pid  <- SLICE2PAT[[sl]]
    stem <- file.path(FIG_DIR, sprintf("%s_%s_%s", prefix, pid, sl))
    sec  <- c(sec, sprintf('<p><strong>%s（%s）</strong> · <code>figures/%s_%s_%s.png</code></p>',
                           sl, pid, prefix, pid, sl), img_tag(stem, width))
  }
  sec
}

sec <- c('<h2 id="sec-extra">补充图：4 类逐切片图（Unit 3）</h2>',
         '<p>以下 44 组图（11 切片 × 4 类）均只读 <code>tables/zone_spot_assignment.csv</code> 绘制，未重算 cnv_arm_score / cnv_norm / zone / clone；每类均提供 PDF 与 PNG。空间图（第 1、3 类）使用与旧图一致的方形坐标框（58–527 px）与 1 mm 比例尺；箱线图与堆叠柱为非空间图，不含比例尺。</p>')
sec <- add_type(sec, "1. CNV 值空间图（viridis 连续色阶）",
  "颜色为 <code>cnv_norm</code>（0–1），越亮表示切片内相对 CNV 越高；色条为连续 0–1，未做二值化。<strong>cnv_norm = 切片内 min–max 归一化，非绝对肿瘤含量</strong>。",
  "zone_cnv_spatial")
sec <- add_type(sec, "2. 三区 CNV 值箱线图",
  "箱体为 IQR，须线为 1.5×IQR，点为单个 spot（仅横向抖动，不改变 y 值）；每区标注 n 与中位数；x 轴顺序 Stroma → Interface → Core，配色与旧空间图一致。",
  "zone_cnv_box", width = 820)
sec <- add_type(sec, "3. clone 空间图",
  "clone 定性着色；同一病人内 clone 颜色跨切片一致；图例括号内为该切片中该 clone 的 spot 数。",
  "zone_clone_spatial")
sec <- add_type(sec, "4. 三区 clone 组成堆叠柱",
  "每柱为一个区，柱高 100%，颜色为 clone 组成；柱顶 n 为该区 spot 数；完整百分比见 <code>tables/clone_zone_pct.csv</code>。",
  "zone_clone_bar", width = 820)

tech_pos <- regexpr('<h2 id="[^"]+">技术说明书</h2>', html)[1]
if (tech_pos < 1) stop("tech-spec heading not found")
html_new <- paste0(substr(html, 1, tech_pos - 1), paste(sec, collapse = "\n"), "\n",
                   substr(html, tech_pos, nchar(html)))

nav_pos <- regexpr('<nav class="toc" id="toc">', html_new)[1]
if (nav_pos < 1) stop("TOC nav not found")
rel_ul <- regexpr('</ul>', substr(html_new, nav_pos, nchar(html_new)))[1]
ul_pos <- nav_pos + rel_ul - 1
html_new <- paste0(substr(html_new, 1, ul_pos - 1), '<li><a href="#sec-extra">补充图</a></li>',
                   "\n", substr(html_new, ul_pos, nchar(html_new)))

writeLines(html_new, html_path, useBytes = TRUE)
cat("HTML size MB:", round(file.info(html_path)$size / 1e6, 2), "\n")
cat("base64 PNG count:", lengths(regmatches(html_new, gregexpr("data:image/png;base64", html_new))), "\n")
cat("section present:", grepl("补充图：4 类逐切片图", html_new, fixed = TRUE),
    "| TOC entry:", grepl('href="#sec-extra"', html_new, fixed = TRUE),
    "| tech-spec kept:", grepl("技术说明书", html_new, fixed = TRUE), "\n")
cat("review token:", grepl("review/", html_new, fixed = TRUE),
    "| absolute home:", grepl("/path/to", html_new, fixed = TRUE), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")

# 仅替换 HTML 中第 2 类（箱线图）小节的 11 张 base64 图，避免重复插入整节
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
html_path <- file.path(OUT, "step1_three_zone_汇总.html")
html <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
img_tag <- function(stem, width = 820) {
  b64 <- gsub("[\r\n]", "", base64enc::base64encode(paste0(stem, ".png")))
  sprintf('<p style="text-align:center"><img src="data:image/png;base64,%s" width="%d"></p>', b64, width)
}
add_block <- function(title, desc, prefix, width = 820) {
  sec <- c(sprintf('<h3>%s</h3>', title), sprintf('<p>%s</p>', desc),
           sprintf('<p>文件：<code>figures/%s_&lt;P&gt;_&lt;slice&gt;.{pdf,png}</code></p>', prefix))
  for (sl in SLICE_ORDER) {
    pid  <- SLICE2PAT[[sl]]
    stem <- file.path(FIG_DIR, sprintf("%s_%s_%s", prefix, pid, sl))
    sec  <- c(sec, sprintf('<p><strong>%s（%s）</strong> · <code>figures/%s_%s_%s.png</code></p>',
                           sl, pid, prefix, pid, sl), img_tag(stem, width))
  }
  sec
}
start_marker <- '<h3>2. 三区 CNV 值箱线图</h3>'
end_marker   <- '<h3>3. clone 空间图</h3>'
p1 <- regexpr(start_marker, html, fixed = TRUE)[1]
p2 <- regexpr(end_marker, html, fixed = TRUE)[1]
if (p1 < 1 || p2 <= p1) stop("type-2 block markers not found")
new_block <- add_block("2. 三区 CNV 值箱线图",
  "箱体为 IQR，须线为 1.5×IQR，箱体绘制在点层之上；点为单个 spot（仅横向抖动，低饱和，不遮挡箱体）；每区标注 n 与中位数；x 轴顺序 Stroma → Interface → Core，配色与旧空间图一致。",
  "zone_cnv_box", width = 820)
html_new <- paste0(substr(html, 1, p1 - 1), paste(new_block, collapse = "\n"), "\n",
                   substr(html, p2, nchar(html)))
writeLines(html_new, html_path, useBytes = TRUE)
cat("HTML MB:", round(file.info(html_path)$size / 1e6, 2),
    "| base64 PNG:", lengths(regmatches(html_new, gregexpr("data:image/png;base64", html_new))),
    "| type2 heading count:", lengths(regmatches(html_new, gregexpr(start_marker, html_new, fixed = TRUE))),
    "| type3 heading count:", lengths(regmatches(html_new, gregexpr(end_marker, html_new, fixed = TRUE))), "\n")
cat("review token:", grepl("review/", html_new, fixed = TRUE),
    "| absolute home:", grepl("/path/to", html_new, fixed = TRUE), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start), "=====\n")
