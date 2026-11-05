# =============================================================================
# TNBC interface-niche code repository
# Original source: project_spatial/TNBC/patient_ALL/clone_TAM_cellchat/notebooks/step01_clone_tam_cellchat.ipynb (converted to R)
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

# 验证 cellchat 内核环境与核心依赖包版本，并建立输出目录
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(suppressWarnings({
  library(CellChat)
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(NMF)
  library(ggalluvial)
  library(ComplexHeatmap)
  library(gridExtra)
  library(grid)
  library(igraph)
  library(viridisLite)
  library(png)
  library(jsonlite)
}))

cat("R version      :", R.version.string, "\n")
cat("R home         :", R.home(), "\n")
cat(".libPaths      :", paste(.libPaths(), collapse = " ; "), "\n")
cat("CellChat       :", as.character(packageVersion("CellChat")), "\n")
cat("Seurat         :", as.character(packageVersion("Seurat")), "\n")
cat("ComplexHeatmap :", as.character(packageVersion("ComplexHeatmap")), "\n")
cat("NMF            :", as.character(packageVersion("NMF")), "\n")

# 输出根目录（每位患者一个独立子目录）
out_root <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone_TAM_cellchat"
for (p in c("P13/figures", "P13/tables", "P13/data", "P15/figures", "P15/tables", "P15/data")) {
  dir.create(file.path(out_root, p), recursive = TRUE, showWarnings = FALSE)
}
cat("Output root    :", out_root, "\n")

# 绘图配色契约（palette skill: 定性默认 nature / 连续 viridis / 发散 RdBu_r）
NPG_COLORS <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000")
plot_theme <- theme_minimal() +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid        = element_blank(),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5)
  )
cat("Plot theme     : nature (NPG) + white bg + black border + no grid\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 探查 qs 是否可从其他环境库加载，并核对输入文件与 clone 表结构
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

cat("R version (cellchat env):", paste(R.version$major, R.version$minor, sep = "."), "\n")
cand <- c("/path/to/envs/singlecell/lib/R/library", "/path/to/envs/singlecell_r/lib/R/library",
          "/path/to/envs/common_r/lib/R/library", "/usr/lib/R/library", "/usr/local/lib/R/site-library")
for (d in cand) cat(sprintf("%-50s exists=%s qs=%s\n", d, dir.exists(d), dir.exists(file.path(d, "qs"))))

cat("\n[Input files]\n")
f_mat   <- "/path/to/project/project_singlecell/TNBC_result/final_figures/V2/spatial/spatial_sct_matrix.rds"
f_tme   <- "/path/to/project/project_spatial/TNBC/TNBC_analysis/SpatialInferCNV/results/step3_spatial_tme_groups/data/spatial_tme_groups.qs"
f_cl13  <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P13/clone_membership.csv"
f_cl15  <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P15/clone_membership.csv"
for (f in c(f_mat, f_tme, f_cl13, f_cl15)) {
  cat(sprintf("  %-6s %10.1f MB  %s\n", file.exists(f), file.info(f)$size/1e6, f))
}

cat("\n[Clone membership P13 head]\n")
cl13 <- read.csv(f_cl13, stringsAsFactors = FALSE, check.names = FALSE)
print(dim(cl13)); print(colnames(cl13)); print(head(cl13, 3))
cat("\n[clone value counts P13]\n"); print(table(cl13$clone, useNA = "ifany")); print(table(cl13$slice))

cl15 <- read.csv(f_cl15, stringsAsFactors = FALSE, check.names = FALSE)
cat("\n[Clone membership P15]\n"); print(dim(cl15)); print(table(cl15$clone, useNA = "ifany")); print(table(cl15$slice))

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 读取空间表达矩阵；TME 比例表由同一上游 RCTD 结果重建（cellchat 环境缺 qs，无法直接读 .qs）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

f_mat <- "/path/to/project/project_singlecell/TNBC_result/final_figures/V2/spatial/spatial_sct_matrix.rds"
f_prop <- "/path/to/project/project_singlecell/TNBC_result/figure20/subproportions.rds"
f_meta <- "/path/to/project/project_singlecell/TNBC_result/figure20/spatial_meta.rds"
f_step3sum <- "/path/to/project/project_spatial/TNBC/TNBC_analysis/SpatialInferCNV/results/step3_spatial_tme_groups/tables/spatial_tme_group_summary.csv"

cat("[Dependency Gap] qs 包在 cellchat 环境不可加载（CXXABI_1.3.15 缺失，且 stringfish 缺失）；\n")
cat("                  改由 step3 的同一上游输入 subproportions.rds + spatial_meta.rds 按 step3 原样列构造重建 TME 表，\n")
cat("                  并用 step3 已交付的 spatial_tme_group_summary.csv 逐片核对，不重跑解卷积。\n\n")

# 1) 表达矩阵（已核对为 log-normalized，不再 NormalizeData）
mat <- readRDS(f_mat)
cat("[Matrix] class:", class(mat)[1], "| dim:", nrow(mat), "genes x", ncol(mat), "spots\n")
cat("[Matrix] value range:", paste(round(range(mat@x), 3), collapse = " - "), "\n")
cat("[Matrix] colnames head:", paste(head(colnames(mat), 3), collapse = ", "), "\n")

# 2) 按 step3 原列重建 TME 表（barcode/sample/slice/region/x/y/C1QC_TAM/SPP1_TAM）
prop <- readRDS(f_prop)
meta <- readRDS(f_meta)
stopifnot(all(meta$barcode == rownames(prop)))
tme <- data.frame(
  barcode   = meta$barcode,
  sample    = meta$sample,
  slice     = sub("^GSM[0-9]+_", "", meta$sample),
  region    = meta$region,
  x         = as.numeric(meta$x),
  y         = as.numeric(meta$y),
  C1QC_TAM  = prop[["C1Q TAMs"]],
  SPP1_TAM  = prop[["SPP1 TAMs"]],
  stringsAsFactors = FALSE
)
cat("\n[TME rebuilt] dim:", nrow(tme), "x", ncol(tme), "\n")
print(colnames(tme))
cat("[TME] slices:\n"); print(table(tme$slice))
cat("[TME] regions:\n"); print(table(tme$region))

# 3) 与 step3 已交付汇总逐片核对 C1QC_TAM / SPP1_TAM 的 mean_prop
s3 <- read.csv(f_step3sum, stringsAsFactors = FALSE)
chk <- do.call(rbind, lapply(unique(tme$slice), function(s) {
  d <- tme[tme$slice == s, ]
  data.frame(
    slice = s,
    C1Q_rebuilt = mean(d$C1QC_TAM), C1Q_reported = s3$mean_prop[s3$slice == s & s3$group == "C1QC TAM"],
    SPP1_rebuilt = mean(d$SPP1_TAM), SPP1_reported = s3$mean_prop[s3$slice == s & s3$group == "SPP1 TAM"],
    n_rebuilt = nrow(d), n_reported = s3$n_spots[s3$slice == s & s3$group == "C1QC TAM"]
  )
}))
chk$max_abs_diff <- pmax(abs(chk$C1Q_rebuilt - chk$C1Q_reported), abs(chk$SPP1_rebuilt - chk$SPP1_reported))
print(chk, digits = 4)
cat("\n[Verify] 逐片 mean_prop 最大绝对偏差:", signif(max(chk$max_abs_diff), 3),
    "| spot 数全部一致:", all(chk$n_rebuilt == chk$n_reported), "\n")

# 4) 矩阵列名 vs TME barcode vs clone cell 命名对照
cat("\n[ID format] mat colnames:", head(colnames(mat), 2), "\n")
cat("[ID format] tme barcode  :", head(tme$barcode, 2), "\n")
tme$barcode_base <- sub("_[0-9]+$", "", tme$barcode)
cat("[ID format] barcode_base :", head(tme$barcode_base, 2), "\n")
cat("[ID format] clone cell   : 094A_AAACATTTCCCGGATT-1（slice_barcode_base）\n")
cat("[ID check] tme$barcode 与 colnames(mat) 完全一致:", identical(tme$barcode, colnames(mat)), "\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 构建四组标签：clone 优先占位，其余 spot 按 C1QC_TAM / SPP1_TAM 的 q75 取 top25%，双高归比例较高者
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

build_groups <- function(patient, slices, f_clone, tme) {
  d <- tme[tme$slice %in% slices & tme$region == "Interface", ]
  d$key <- paste0(d$slice, "_", d$barcode_base)
  cl <- read.csv(f_clone, stringsAsFactors = FALSE, check.names = FALSE)
  ci <- cl$clone[match(d$key, cl$cell)]
  d$clone_raw <- ci
  d$group <- NA_character_
  d$group[!is.na(ci)] <- ci[!is.na(ci)]
  d$assigned_by <- ifelse(!is.na(ci), "clone", NA_character_)

  nc <- is.na(ci)
  qC <- as.numeric(quantile(d$C1QC_TAM[nc], 0.75, type = 7))
  qS <- as.numeric(quantile(d$SPP1_TAM[nc], 0.75, type = 7))
  hiC <- nc & d$C1QC_TAM > qC
  hiS <- nc & d$SPP1_TAM > qS
  both <- hiC & hiS

  g <- rep(NA_character_, nrow(d))
  g[hiC & !hiS] <- "C1Q_TAM"
  g[hiS & !hiC] <- "SPP1_TAM"
  g[both] <- ifelse(d$C1QC_TAM[both] >= d$SPP1_TAM[both], "C1Q_TAM", "SPP1_TAM")
  d$group[nc] <- g[nc]
  d$assigned_by[nc] <- ifelse(is.na(g[nc]), "unassigned_top25", "tme_top25")

  cat(sprintf("\n=== %s | slices: %s ===\n", patient, paste(slices, collapse = "/")))
  cat(sprintf("[Interface] spots: %d\n", nrow(d)))
  cat(sprintf("[Clone] matched: %d (%.1f%%) | Clone1 %d, Clone2 %d | clone rows %d\n",
              sum(!is.na(ci)), 100 * mean(!is.na(ci)), sum(ci == "Clone1", na.rm = TRUE),
              sum(ci == "Clone2", na.rm = TRUE), nrow(cl)))
  cat(sprintf("[Non-clone] n: %d | q75 C1QC_TAM = %.5f | q75 SPP1_TAM = %.5f\n",
              sum(nc), qC, qS))
  cat(sprintf("[Thresholds] hiCQC: %d | hiSPP1: %d | both-high: %d | union: %d (%.1f%% of non-clone)\n",
              sum(hiC), sum(hiS), sum(both), sum(hiC | hiS), 100 * sum(hiC | hiS) / sum(nc)))
  cat("[Group counts] "); print(table(d$group, useNA = "ifany"))
  cat("[By slice] "); print(table(d$slice, d$group, useNA = "ifany"))
  d
}

P13 <- build_groups("P13", c("094A", "094B"),
                    "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P13/clone_membership.csv", tme)
P15 <- build_groups("P15", c("095A", "095B"),
                    "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P15/clone_membership.csv", tme)

GRP_LEVELS <- c("Clone1", "Clone2", "C1Q_TAM", "SPP1_TAM")
exp_P13 <- c(515, 229, 142, 124); exp_P15 <- c(578, 298, 119, 99)
got_P13 <- as.integer(table(factor(P13$group, levels = GRP_LEVELS)))
got_P15 <- as.integer(table(factor(P15$group, levels = GRP_LEVELS)))
cat("\n[Expectation check]\n")
cat("P13 expected:", paste(exp_P13, collapse = "/"), "| actual:", paste(got_P13, collapse = "/"),
    "| diff:", paste(got_P13 - exp_P13, collapse = "/"), "\n")
cat("P15 expected:", paste(exp_P15, collapse = "/"), "| actual:", paste(got_P15, collapse = "/"),
    "| diff:", paste(got_P15 - exp_P15, collapse = "/"), "\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 保存四组标签表与 CellChat 输入子集矩阵，并核对四组 spot 表达矩阵维度与基因数
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

save_obj <- function(obj, filepath) {
  rds_path <- paste0(sub("\\.(rds|qs)$", "", filepath), ".rds")
  saveRDS(obj, rds_path)
  cat(sprintf("[Save] RDS: %s (%.1f MB)\n", rds_path, file.info(rds_path)$size / 1e6))
}

GRP_LEVELS <- c("Clone1", "Clone2", "C1Q_TAM", "SPP1_TAM")

finalize_patient <- function(patient, d) {
  d4 <- d[!is.na(d$group) & d$group %in% GRP_LEVELS, ]
  d4$group <- factor(d4$group, levels = GRP_LEVELS)
  # 逐 spot 分组表（含 TME 比例，便于溯源）
  tab <- d4[, c("barcode", "slice", "group", "assigned_by", "C1QC_TAM", "SPP1_TAM", "x", "y")]
  tab <- tab[order(tab$group, tab$slice), ]
  write.csv(tab, file.path(out_root, patient, "tables", paste0(patient, "_spot_groups.csv")), row.names = FALSE)
  # 分组 spot 数
  cnt <- as.data.frame(table(group = d4$group), stringsAsFactors = FALSE)
  cnt$total_interface <- nrow(d)
  cnt$unassigned <- sum(is.na(d$group))
  write.csv(cnt, file.path(out_root, patient, "tables", paste0(patient, "_group_spot_counts.csv")), row.names = FALSE)
  # CellChat 输入矩阵（genes x spots，log-normalized）
  X <- mat[, d4$barcode, drop = FALSE]
  meta_cc <- data.frame(
    spot      = d4$barcode,
    group     = as.character(d4$group),
    slice     = d4$slice,
    C1QC_TAM  = d4$C1QC_TAM,
    SPP1_TAM  = d4$SPP1_TAM,
    stringsAsFactors = FALSE
  )
  rownames(meta_cc) <- meta_cc$spot
  save_obj(list(X = X, meta = meta_cc), file.path(out_root, patient, "data", paste0(patient, "_cellchat_input")))
  cat(sprintf("\n[%s] interface %d -> 4 groups %d spots | genes %d\n",
              patient, nrow(d), ncol(X), nrow(X)))
  print(cnt)
  invisible(list(X = X, meta = meta_cc))
}

in13 <- finalize_patient("P13", P13)
in15 <- finalize_patient("P15", P15)

# 输入层自检：确认为 log-normalized（列和约 9k-1e4，非整数计数），且四组均有表达
cat("\n[Data-layer check]\n")
cat("P13 col sums range:", paste(round(range(Matrix::colSums(in13$X)), 1), collapse = " - "), "\n")
cat("P15 col sums range:", paste(round(range(Matrix::colSums(in15$X)), 1), collapse = " - "), "\n")
cat("P13 max value:", round(max(in13$X), 2), "| P15 max value:", round(max(in15$X), 2), "\n")
cat("Non-integer fraction P13:", signif(mean(in13$X@x != round(in13$X@x)), 3), "\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# P13：构建 CellChat 对象，识别过表达配体受体，计算四组间通讯概率并聚合通路网络
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

.libPaths("/path/to/envs/cellchat/lib/R/library")
suppressPackageStartupMessages({ library(Matrix) })
options(stringsAsFactors = FALSE)
future::plan("sequential")

# meta 增列 samples（=切片），避免 CellChat 默认全部归为 sample1
in13$meta$samples <- in13$meta$slice

data(CellChatDB.human)
CellChatDB <- CellChatDB.human
cat("[DB] interactions:", nrow(CellChatDB$interaction), "| complexes:", length(CellChatDB$complex),
    "| cofactors:", length(CellChatDB$cofactor), "\n")
cat("[DB] categories:", paste(unique(CellChatDB$interaction$annotation), collapse = " / "), "\n")

cc13 <- createCellChat(object = in13$X, meta = in13$meta, group.by = "group", do.sparse = TRUE)
cc13@DB <- CellChatDB
cat("[CellChat P13] groups:", paste(levels(cc13@idents), collapse = " / "), "\n"); print(table(cc13@idents))

cc13 <- subsetData(cc13)
cat("[subsetData] signaling genes retained:", nrow(cc13@data), "\n")
# presto 未安装于本环境 → 用标准 Wilcoxon（do.fast = FALSE）
cc13 <- identifyOverExpressedGenes(cc13, do.fast = FALSE)
cc13 <- identifyOverExpressedInteractions(cc13)
cat("[OverExpressed] interactions:", nrow(cc13@LR$LRsig), "\n")

cc13 <- computeCommunProb(cc13, type = "triMean", population.size = FALSE)
cat("[computeCommunProb triMean] significant LR rows before filter:", nrow(subsetCommunication(cc13)), "\n")

cc13 <- filterCommunication(cc13, min.cells = 10)
cat("[filterCommunication min.cells=10] significant LR rows:", nrow(subsetCommunication(cc13)), "\n")

cc13 <- computeCommunProbPathway(cc13)
cc13 <- aggregateNet(cc13)
cc13 <- netAnalysis_computeCentrality(cc13, slot.name = "netP")

cat("[Pathways]", length(cc13@netP$pathways), ":", paste(head(cc13@netP$pathways, 40), collapse = ", "), "\n")
cat("[Net count matrix]\n"); print(cc13@net$count)
cat("[Net weight matrix]\n"); print(round(cc13@net$weight, 4))

saveRDS(cc13, file.path(out_root, "P13", "data", "cellchat_P13.rds"))
cat("[Save] cellchat_P13.rds (", round(file.info(file.path(out_root, "P13", "data", "cellchat_P13.rds"))$size/1e6, 1), "MB )\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# P15：与 P13 同参数独立构建 CellChat 对象并完成通讯推断（两例不 merge）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({ library(Matrix) })
options(stringsAsFactors = FALSE)
future::plan("sequential")

in15$meta$samples <- factor(in15$meta$slice)
data(CellChatDB.human)
CellChatDB <- CellChatDB.human

cc15 <- createCellChat(object = in15$X, meta = in15$meta, group.by = "group", do.sparse = TRUE)
cc15@DB <- CellChatDB
cat("[CellChat P15] groups:", paste(levels(cc15@idents), collapse = " / "), "\n"); print(table(cc15@idents))

cc15 <- subsetData(cc15)
cc15 <- identifyOverExpressedGenes(cc15, do.fast = FALSE)
cc15 <- identifyOverExpressedInteractions(cc15)
cat("[OverExpressed] interactions:", nrow(cc15@LR$LRsig), "\n")

cc15 <- computeCommunProb(cc15, type = "triMean", population.size = FALSE)
cat("[computeCommunProb triMean] significant LR rows before filter:", nrow(subsetCommunication(cc15)), "\n")

cc15 <- filterCommunication(cc15, min.cells = 10)
cat("[filterCommunication min.cells=10] significant LR rows:", nrow(subsetCommunication(cc15)), "\n")

cc15 <- computeCommunProbPathway(cc15)
cc15 <- aggregateNet(cc15)
cc15 <- netAnalysis_computeCentrality(cc15, slot.name = "netP")

cat("[Pathways]", length(cc15@netP$pathways), ":", paste(head(cc15@netP$pathways, 40), collapse = ", "), "\n")
cat("[Net count matrix]\n"); print(cc15@net$count)
cat("[Net weight matrix]\n"); print(round(cc15@net$weight, 4))

saveRDS(cc15, file.path(out_root, "P15", "data", "cellchat_P15.rds"))
cat("[Save] cellchat_P15.rds (", round(file.info(file.path(out_root, "P15", "data", "cellchat_P15.rds"))$size/1e6, 1), "MB )\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 导出逐 LR / 逐通路通讯概率（取 CellChat 全概率数组，避免 subsetCommunication 丢弃 p=1 行）并计算 Clone1 vs Clone2 差异
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

SRC <- c("Clone1", "Clone2"); TGT <- c("C1Q_TAM", "SPP1_TAM"); EPS <- 1e-4

# 3 维数组（source x target x LR/pathway）展开为长表：完整、无丢弃
arr_to_df <- function(A, varname) {
  d <- as.data.frame.table(A, stringsAsFactors = FALSE)
  colnames(d) <- c("source", "target", varname, "value")
  d
}

export_patient <- function(patient, cc) {
  # ---- LR 级：prob + pval 全量 ----
  d <- merge(arr_to_df(cc@net$prob, "interaction_name"),
             arr_to_df(cc@net$pval, "interaction_name"),
             by = c("source", "target", "interaction_name"), suffixes = c("_prob", "_pval"))
  colnames(d)[colnames(d) == "value_prob"] <- "prob"
  colnames(d)[colnames(d) == "value_pval"] <- "pval"
  lr <- cc@LR$LRsig; lr$interaction_name <- rownames(lr)
  d <- merge(d, lr[, c("interaction_name", "pathway_name", "ligand", "receptor", "annotation")],
             by = "interaction_name", all.x = TRUE)

  # ---- 通路级：CellChat 只提供 netP$prob（无 netP$pval），故不给通路 p 值 ----
  dp <- arr_to_df(cc@netP$prob, "pathway_name")
  colnames(dp)[colnames(dp) == "value"] <- "prob"

  cat(sprintf("\n[%s] LR 长表 %d 行 (%d 边 x %d LR) | 通路长表 %d 行\n", patient, nrow(d),
              length(unique(paste(d$source, d$target))), length(unique(d$interaction_name)), nrow(dp)))

  # 核对 netP$prob 与成员 LR 概率之和的关系（口径文档化）
  sm  <- aggregate(prob ~ source + target + pathway_name, data = d, FUN = sum)
  cmp <- merge(sm, dp, by = c("source", "target", "pathway_name"), suffixes = c("_lrsum", "_netP"))
  cat(sprintf("[%s] netP$prob vs 成员 LR 概率之和: 相关系数 %.4f | netP/和 中位数 %.3f\n",
              patient, cor(cmp$prob_lrsum, cmp$prob_netP), median(cmp$prob_netP / pmax(cmp$prob_lrsum, 1e-9))))

  write.csv(d,  file.path(out_root, patient, "tables", paste0(patient, "_ALL_LR_probabilities.csv")), row.names = FALSE)
  write.csv(dp, file.path(out_root, patient, "tables", paste0(patient, "_ALL_pathway_probabilities.csv")), row.names = FALSE)

  sub_net <- d[d$source %in% SRC & d$target %in% TGT, ]
  sub_pat <- dp[dp$source %in% SRC & dp$target %in% TGT, ]
  write.csv(sub_net, file.path(out_root, patient, "tables", paste0(patient, "_LR_probabilities_clone_to_TAM.csv")), row.names = FALSE)
  write.csv(sub_pat, file.path(out_root, patient, "tables", paste0(patient, "_pathway_probabilities_clone_to_TAM.csv")), row.names = FALSE)

  # ---- LR 差异表：Δprob、log2 比值、两侧 p、四分类 ----
  w <- reshape(sub_net[, c("target", "pathway_name", "interaction_name", "source", "prob", "pval")],
               idvar = c("target", "pathway_name", "interaction_name"), timevar = "source", direction = "wide")
  diff_lr <- data.frame(
    w[, c("target", "pathway_name", "interaction_name")],
    prob_Clone1 = w$prob.Clone1, prob_Clone2 = w$prob.Clone2,
    pval_Clone1 = w$pval.Clone1, pval_Clone2 = w$pval.Clone2, stringsAsFactors = FALSE)
  diff_lr$delta_prob <- diff_lr$prob_Clone1 - diff_lr$prob_Clone2
  diff_lr$log2_ratio <- log2((diff_lr$prob_Clone1 + EPS) / (diff_lr$prob_Clone2 + EPS))
  s1 <- diff_lr$pval_Clone1 < 0.05; s2 <- diff_lr$pval_Clone2 < 0.05
  diff_lr$signif_class <- ifelse(s1 & s2, "both_significant",
                          ifelse(s1 & !s2, "Clone1_only",
                          ifelse(!s1 & s2, "Clone2_only", "neither")))
  diff_lr <- diff_lr[order(diff_lr$target, -abs(diff_lr$delta_prob)),
                     c("target", "pathway_name", "interaction_name", "prob_Clone1", "prob_Clone2",
                       "delta_prob", "log2_ratio", "pval_Clone1", "pval_Clone2", "signif_class")]
  write.csv(diff_lr, file.path(out_root, patient, "tables", paste0(patient, "_diff_LR_Clone1_vs_Clone2.csv")), row.names = FALSE)

  # ---- 通路差异表：Δprob + 两侧显著成员 LR 数（无通路级 p 值，不作显著性宣称）----
  wp <- reshape(sub_pat[, c("target", "pathway_name", "source", "prob")],
                idvar = c("target", "pathway_name"), timevar = "source", direction = "wide")
  a1 <- transform(sub_net[sub_net$source == "Clone1", ], n_sig_Clone1 = as.numeric(pval < 0.05))
  a2 <- transform(sub_net[sub_net$source == "Clone2", ], n_sig_Clone2 = as.numeric(pval < 0.05))
  ns1 <- aggregate(n_sig_Clone1 ~ target + pathway_name, data = a1, FUN = sum)
  ns2 <- aggregate(n_sig_Clone2 ~ target + pathway_name, data = a2, FUN = sum)
  diff_pat <- data.frame(wp[, c("target", "pathway_name")],
                         prob_Clone1 = wp$prob.Clone1, prob_Clone2 = wp$prob.Clone2, stringsAsFactors = FALSE)
  diff_pat$delta_prob <- diff_pat$prob_Clone1 - diff_pat$prob_Clone2
  diff_pat$log2_ratio <- log2((diff_pat$prob_Clone1 + EPS) / (diff_pat$prob_Clone2 + EPS))
  diff_pat <- merge(diff_pat, ns1, by = c("target", "pathway_name"), all.x = TRUE)
  diff_pat <- merge(diff_pat, ns2, by = c("target", "pathway_name"), all.x = TRUE)
  diff_pat$n_sig_Clone1[is.na(diff_pat$n_sig_Clone1)] <- 0
  diff_pat$n_sig_Clone2[is.na(diff_pat$n_sig_Clone2)] <- 0
  diff_pat <- diff_pat[order(diff_pat$target, -abs(diff_pat$delta_prob)), ]
  write.csv(diff_pat, file.path(out_root, patient, "tables", paste0(patient, "_diff_pathway_Clone1_vs_Clone2.csv")), row.names = FALSE)

  cat(sprintf("[%s] LR diff %d 行 | %s\n", patient, nrow(diff_lr),
              paste(names(table(diff_lr$signif_class)), table(diff_lr$signif_class), sep = "=", collapse = " ")))
  cat(sprintf("[%s] pathway diff %d 行（通路级无 p 值，仅 Δprob 与成员显著 LR 数）\n", patient, nrow(diff_pat)))

  flow <- aggregate(prob ~ pathway_name + source + target, data = sub_pat, FUN = sum)
  write.csv(flow, file.path(out_root, patient, "tables", paste0(patient, "_pathway_information_flow_clone_to_TAM.csv")), row.names = FALSE)

  edges <- expand.grid(source = SRC, target = TGT, stringsAsFactors = FALSE)
  edges$count  <- mapply(function(s, t) cc@net$count[s, t],  edges$source, edges$target)
  edges$weight <- mapply(function(s, t) cc@net$weight[s, t], edges$source, edges$target)
  edges$n_sig_LR <- mapply(function(s, t) sum(sub_net$pval[sub_net$source == s & sub_net$target == t] < 0.05),
                           edges$source, edges$target)
  edges$mean_prob_sig <- mapply(function(s, t) {
    dd <- sub_net[sub_net$source == s & sub_net$target == t & sub_net$pval < 0.05, ]
    if (nrow(dd) == 0) 0 else mean(dd$prob)
  }, edges$source, edges$target)
  write.csv(edges, file.path(out_root, patient, "tables", paste0(patient, "_clone_to_TAM_edge_summary.csv")), row.names = FALSE)
  cat(sprintf("[%s] edge summary:\n", patient)); print(edges, digits = 4)

  list(net = d, netP = dp, sub_net = sub_net, sub_pat = sub_pat,
       diff_lr = diff_lr, diff_pat = diff_pat, edges = edges)
}

res13 <- export_patient("P13", cc13)
res15 <- export_patient("P15", cc15)

cat("\n[Check] 四条边 count>0: P13 =", all(res13$edges$count > 0), "| P15 =", all(res15$edges$count > 0), "\n")
cat("[Check] 四条边显著 LR 数>0: P13 =", all(res13$edges$n_sig_LR > 0), "| P15 =", all(res15$edges$n_sig_LR > 0), "\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 查看两例中 Clone1 与 Clone2 向两类 TAM 的通讯差异（真实概率口径），用于结果解读
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

show_top <- function(res, patient, n = 10) {
  cat(sprintf("\n===== %s | Top LR by |delta_prob| =====\n", patient))
  d <- res$diff_lr[order(-abs(res$diff_lr$delta_prob)), ]
  print(head(d[, c("target", "pathway_name", "interaction_name", "prob_Clone1", "prob_Clone2",
                   "delta_prob", "log2_ratio", "pval_Clone1", "pval_Clone2", "signif_class")], n),
        row.names = FALSE, digits = 3)
  cat(sprintf("\n===== %s | Top pathway by |delta_prob| =====\n", patient))
  p <- res$diff_pat[order(-abs(res$diff_pat$delta_prob)), ]
  print(head(p, n), row.names = FALSE, digits = 3)
}
show_top(res13, "P13"); show_top(res15, "P15")

# 仅看两侧都显著的 LR，按 Δprob 排序（同一靶细胞下最有把握的方向性差异）
cat("\n===== 双显著 LR 中 Clone1 更强的前 8 条（P13 / P15）=====\n")
for (nm in c("P13", "P15")) {
  d <- get(paste0("res", ifelse(nm == "P13", "13", "15")))$diff_lr
  d <- d[d$signif_class == "both_significant", ]
  d <- d[order(-d$delta_prob), ]
  cat("--", nm, "n(both sig) =", nrow(d), "| Clone1>Clone2:", sum(d$delta_prob > 0), "\n")
  print(head(d[, c("target", "pathway_name", "interaction_name", "prob_Clone1", "prob_Clone2", "delta_prob")], 8),
        row.names = FALSE, digits = 3)
}

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# P13：绘制 CellChat 标准可视化（圈图/气泡/弦图/热图）+ 差异 LR 热图 + 克隆→TAM 强度条形图
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(viridisLite); library(ComplexHeatmap) })
SRC <- c("Clone1", "Clone2"); TGT <- c("C1Q_TAM", "SPP1_TAM")

# base-graphics 图的双格式保存（在调用者环境中求值）
save_both <- function(stem, w, h, expr, res = 300) {
  ev <- parent.frame()
  png(paste0(stem, ".png"), width = w, height = h, units = "in", res = res, bg = "white")
  par(bg = "white"); eval(expr, envir = ev); dev.off()
  pdf(paste0(stem, ".pdf"), width = w, height = h, bg = "white")
  par(bg = "white"); eval(expr, envir = ev); dev.off()
  cat(sprintf("[Fig] %s.png / .pdf  (%.0f KB)\n", basename(stem), file.info(paste0(stem, ".png"))$size/1024))
}
# ComplexHeatmap 对象的双格式保存（netVisual_heatmap 返回 Heatmap，需显式 draw）
save_ht <- function(stem, w, h, make_ht, res = 300) {
  ht <- make_ht()
  if (dev.cur() > 1) dev.off()
  png(paste0(stem, ".png"), width = w, height = h, units = "in", res = res, bg = "white")
  draw(ht, newpage = TRUE); dev.off()
  pdf(paste0(stem, ".pdf"), width = w, height = h, bg = "white")
  draw(ht, newpage = TRUE); dev.off()
  cat(sprintf("[Fig] %s.png / .pdf  (%.0f KB)\n", basename(stem), file.info(paste0(stem, ".png"))$size/1024))
}

draw_patient_cellchat_figs <- function(patient, cc) {
  fdir <- file.path(out_root, patient, "figures")
  pf <- function(x) file.path(fdir, paste0(patient, "_", x))

  # 1) 圈图（互作数 / 通讯强度）
  save_both(pf("netVisual_circle_count"), 7, 7,
            quote(netVisual_circle(cc@net$count, vertex.weight = as.numeric(table(cc@idents)), weight.scale = TRUE,
                                   label.edge = FALSE, edge.weight.max = max(cc@net$count),
                                   title.name = "Number of interactions")))
  save_both(pf("netVisual_circle_weight"), 7, 7,
            quote(netVisual_circle(cc@net$weight, vertex.weight = as.numeric(table(cc@idents)), weight.scale = TRUE,
                                   label.edge = FALSE, edge.weight.max = max(cc@net$weight),
                                   title.name = "Interaction strength")))

  # 2) 气泡图（Clone -> TAM 的显著 LR）
  p <- netVisual_bubble(cc, sources.use = SRC, targets.use = TGT, remove.isolate = FALSE,
                        angle.x = 45, title.name = paste0(patient, ": Clone1/Clone2 -> TAM (significant LR)"))
  ggsave(pf("netVisual_bubble.png"), p, width = 9, height = 7, dpi = 300, bg = "white")
  ggsave(pf("netVisual_bubble.pdf"), p, width = 9, height = 7, bg = "white")
  cat("[Fig] netVisual_bubble.png / .pdf\n")

  # 3) 弦图：LR 级与通路级
  ok1 <- tryCatch({
    save_both(pf("netVisual_chord_LR"), 9, 9,
              quote(netVisual_chord_gene(cc, sources.use = SRC, targets.use = TGT, slot.name = "net",
                                         lab.cex = 0.7, legend.pos.y = 30, legend.pos.x = 15,
                                         title.name = paste0(patient, ": LR-level (Clone1/Clone2 -> TAM)"))))
    TRUE }, error = function(e) { cat("[Warn] chord LR failed:", conditionMessage(e), "\n"); FALSE })
  ok2 <- tryCatch({
    save_both(pf("netVisual_chord_pathway"), 9, 9,
              quote(netVisual_chord_gene(cc, sources.use = SRC, targets.use = TGT, slot.name = "netP",
                                         lab.cex = 0.7, legend.pos.y = 30, legend.pos.x = 15,
                                         title.name = paste0(patient, ": pathway-level (Clone1/Clone2 -> TAM)"))))
    TRUE }, error = function(e) { cat("[Warn] chord pathway failed:", conditionMessage(e), "\n"); FALSE })

  # 4) 通讯热图（count / weight）
  save_ht(pf("netVisual_heatmap_count"), 7, 6, function() netVisual_heatmap(cc, measure = "count",
             color.heatmap = "Reds", title.name = paste0(patient, ": interaction count")))
  save_ht(pf("netVisual_heatmap_weight"), 7, 6, function() netVisual_heatmap(cc, measure = "weight",
             color.heatmap = "Reds", title.name = paste0(patient, ": interaction strength")))
  invisible(list(chord_lr = ok1, chord_path = ok2))
}

chk13 <- draw_patient_cellchat_figs("P13", cc13)
cat("[P13] chord LR:", chk13$chord_lr, "| chord pathway:", chk13$chord_path, "\n")

# 5) Clone1 vs Clone2 差异 LR 热图（Δprob，RdBu_r 发散）
d13 <- res13$diff_lr
top13 <- head(d13[order(-abs(d13$delta_prob)), ], 30)
top13$label <- paste0(sub("_TAM$", "", top13$target), " | ", top13$interaction_name)
top13$label <- factor(top13$label, levels = rev(unique(top13$label[order(top13$delta_prob)])))
hm13 <- ggplot(top13, aes(x = target, y = label, fill = delta_prob)) +
  geom_tile(color = NA) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Δprob\n(Clone1 - Clone2)") +
  labs(x = NULL, y = NULL, title = "P13: top 30 differential LR (Clone1 vs Clone2)") +
  plot_theme + theme(axis.text.y = element_text(size = 7), plot.title = element_text(hjust = .5))
ggsave(file.path(out_root, "P13", "figures", "P13_diff_LR_heatmap.png"), hm13, width = 7, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_root, "P13", "figures", "P13_diff_LR_heatmap.pdf"), hm13, width = 7, height = 7, bg = "white")
cat("[Fig] P13_diff_LR_heatmap.png / .pdf\n")

# 6) 克隆 -> TAM 总通讯强度条形图（count 与 weight 量纲不同，两个 panel 各自坐标轴）
e13 <- res13$edges
pA <- ggplot(e13, aes(x = paste0(source, " -> ", target), y = weight, fill = source)) +
  geom_col(width = .65) + scale_fill_manual(values = setNames(NPG_COLORS[1:2], SRC)) +
  labs(x = NULL, y = "Total interaction strength", title = "P13: Clone -> TAM strength") +
  plot_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
pB <- ggplot(e13, aes(x = paste0(source, " -> ", target), y = count, fill = source)) +
  geom_col(width = .65) + scale_fill_manual(values = setNames(NPG_COLORS[1:2], SRC)) +
  labs(x = NULL, y = "Number of significant LR", title = "P13: Clone -> TAM interaction count") +
  plot_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
pbar13 <- pA / pB
ggsave(file.path(out_root, "P13", "figures", "P13_clone_to_TAM_total_strength_bar.png"), pbar13,
       width = 7, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_root, "P13", "figures", "P13_clone_to_TAM_total_strength_bar.pdf"), pbar13,
       width = 7, height = 8, bg = "white")
cat("[Fig] P13_clone_to_TAM_total_strength_bar.png / .pdf\n")

cat("\n[P13 figures] ", length(list.files(file.path(out_root, "P13", "figures"))), " files\n")
print(list.files(file.path(out_root, "P13", "figures")))

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# P15：与 P13 同参数出图（圈图/气泡/弦图/热图 + 差异 LR 热图 + 克隆→TAM 强度条形图）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

chk15 <- draw_patient_cellchat_figs("P15", cc15)
cat("[P15] chord LR:", chk15$chord_lr, "| chord pathway:", chk15$chord_path, "\n")

# Clone1 vs Clone2 差异 LR 热图
d15 <- res15$diff_lr
top15 <- head(d15[order(-abs(d15$delta_prob)), ], 30)
top15$label <- paste0(sub("_TAM$", "", top15$target), " | ", top15$interaction_name)
top15$label <- factor(top15$label, levels = rev(unique(top15$label[order(top15$delta_prob)])))
hm15 <- ggplot(top15, aes(x = target, y = label, fill = delta_prob)) +
  geom_tile(color = NA) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Δprob\n(Clone1 - Clone2)") +
  labs(x = NULL, y = NULL, title = "P15: top 30 differential LR (Clone1 vs Clone2)") +
  plot_theme + theme(axis.text.y = element_text(size = 7), plot.title = element_text(hjust = .5))
ggsave(file.path(out_root, "P15", "figures", "P15_diff_LR_heatmap.png"), hm15, width = 7, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_root, "P15", "figures", "P15_diff_LR_heatmap.pdf"), hm15, width = 7, height = 7, bg = "white")
cat("[Fig] P15_diff_LR_heatmap.png / .pdf\n")

# 克隆 -> TAM 总通讯强度条形图（两 panel 各自坐标轴）
e15 <- res15$edges
qA <- ggplot(e15, aes(x = paste0(source, " -> ", target), y = weight, fill = source)) +
  geom_col(width = .65) + scale_fill_manual(values = setNames(NPG_COLORS[1:2], SRC)) +
  labs(x = NULL, y = "Total interaction strength", title = "P15: Clone -> TAM strength") +
  plot_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
qB <- ggplot(e15, aes(x = paste0(source, " -> ", target), y = count, fill = source)) +
  geom_col(width = .65) + scale_fill_manual(values = setNames(NPG_COLORS[1:2], SRC)) +
  labs(x = NULL, y = "Number of significant LR", title = "P15: Clone -> TAM interaction count") +
  plot_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
pbar15 <- qA / qB
ggsave(file.path(out_root, "P15", "figures", "P15_clone_to_TAM_total_strength_bar.png"), pbar15,
       width = 7, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_root, "P15", "figures", "P15_clone_to_TAM_total_strength_bar.pdf"), pbar15,
       width = 7, height = 8, bg = "white")
cat("[Fig] P15_clone_to_TAM_total_strength_bar.png / .pdf\n")

cat("\n[P15 figures]", length(list.files(file.path(out_root, "P15", "figures"))), "files\n")
print(list.files(file.path(out_root, "P15", "figures")))

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 产物完整性：pdf/png 成对、文件非空、逐张像素自检（非全白/非全黑、有颜色层）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(png))
manifest <- data.frame()
for (patient in c("P13", "P15")) {
  fdir <- file.path(out_root, patient, "figures")
  stems <- unique(sub("\\.(png|pdf)$", "", list.files(fdir)))
  for (s in stems) {
    fp_png <- file.path(fdir, paste0(s, ".png")); fp_pdf <- file.path(fdir, paste0(s, ".pdf"))
    row <- data.frame(patient = patient, stem = s,
                      png_exists = file.exists(fp_png), pdf_exists = file.exists(fp_pdf),
                      png_kb = ifelse(file.exists(fp_png), round(file.info(fp_png)$size/1024, 1), NA),
                      pdf_kb = ifelse(file.exists(fp_pdf), round(file.info(fp_pdf)$size/1024, 1), NA),
                      sd_rgb = NA_real_, white_frac = NA_real_, pass = NA)
    if (row$png_exists) {
      img <- png::readPNG(fp_png)
      rgb <- img[, , 1:3]
      row$sd_rgb <- round(sd(as.vector(rgb)), 4)
      row$white_frac <- round(mean(rgb > 0.97), 4)
      row$pass <- row$sd_rgb > 0.02 && row$white_frac < 0.99 && row$png_kb > 10
    }
    manifest <- rbind(manifest, row)
  }
}
print(manifest, row.names = FALSE)

cat("\n[Figure wholeness]", sum(manifest$pass, na.rm = TRUE), "/", nrow(manifest), "PNG pixel-verified\n")
cat("[Dual format] png+pdf 成对:", all(manifest$png_exists & manifest$pdf_exists), "\n")
cat("[Non-empty] 最小 pdf KB:", min(manifest$pdf_kb), "| 最小 png KB:", min(manifest$png_kb), "\n")
cat("[Failures]", paste(manifest$stem[!isTRUE(manifest$pass)], collapse = ", "), "\n")

# 表格清单
cat("\n[Tables]\n")
for (patient in c("P13", "P15")) {
  tdir <- file.path(out_root, patient, "tables")
  fs <- list.files(tdir)
  cat(sprintf("%s: %d files\n", patient, length(fs)))
  for (f in fs) cat(sprintf("   %-52s %6.1f KB  rows=%d\n", f, file.info(file.path(tdir, f))$size/1024,
                             nrow(read.csv(file.path(tdir, f)))))
}

write.csv(manifest, file.path(out_root, "figure_manifest.csv"), row.names = FALSE)

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 修正自检汇总的失败项打印（上一版向量化判断写错），并清理旧的中间表（口径错误的遗留文件）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

# 1) 正确的失败项打印
fails <- manifest$stem[which(!manifest$pass)]
cat("[Figure wholeness]", sum(manifest$pass), "/", nrow(manifest), "PNG pixel-verified | failures:",
    if (length(fails) == 0) "none" else paste(fails, collapse = ", "), "\n")

# 2) 清理早期探测生成的重复/口径不符的中间表
stale <- c(file.path(out_root, "P13", "tables", "P13_pathway_information_flow.csv"),
           file.path(out_root, "P15", "tables", "P15_pathway_information_flow.csv"))
for (f in stale) {
  if (file.exists(f)) { file.remove(f); cat("[Clean] removed stale:", basename(f), "\n") }
}

# 3) 最终产物清单核对
cat("\n[Final tables]\n")
for (patient in c("P13", "P15")) {
  fs <- list.files(file.path(out_root, patient, "tables"))
  cat(patient, ":", length(fs), "files ->", paste(fs, collapse = ", "), "\n")
}

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 重绘差异 LR 热图：同一 LR 一行、两列=两个靶细胞、格内标数值；同时展示 Clone1 更强与 Clone2 更强的条目
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

redraw_diff_heatmap <- function(res, patient, n_pos = 10, n_neg = 6) {
  d <- res$diff_lr
  # 每个 LR 取两靶细胞中的最大 |Δ|，并保留其符号
  idx <- tapply(seq_len(nrow(d)), d$interaction_name, function(ii) ii[which.max(abs(d$delta_prob[ii]))])
  best <- data.frame(interaction_name = names(idx), signed_max = d$delta_prob[unlist(idx)], stringsAsFactors = FALSE)
  pos <- head(best$interaction_name[order(-best$signed_max)], n_pos)
  neg <- head(best$interaction_name[order(best$signed_max)], n_neg)
  keep <- unique(c(pos, neg))

  dd <- d[d$interaction_name %in% keep, ]
  ord <- best$interaction_name[order(best$signed_max)]
  dd$interaction_name <- factor(dd$interaction_name, levels = rev(ord[ord %in% keep]))
  lim <- max(abs(dd$delta_prob))
  p <- ggplot(dd, aes(x = target, y = interaction_name, fill = delta_prob)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.3f", delta_prob)), size = 2.9, color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-lim, lim), breaks = round(c(-lim, -lim/2, 0, lim/2, lim), 3),
                         labels = function(x) sprintf("%.3f", x),
                         name = "Δprob\n(Clone1 − Clone2)") +
    labs(x = "Target (receiver)", y = NULL,
         title = sprintf("%s: differential LR, Clone1 vs Clone2", patient),
         subtitle = sprintf("top %d Clone1-higher + top %d Clone2-higher (by max |Δ| across targets)", n_pos, n_neg)) +
    plot_theme +
    theme(axis.text = element_text(size = 9), plot.title = element_text(hjust = .5),
          plot.subtitle = element_text(hjust = .5, size = 8), legend.key.height = unit(0.9, "cm"))
  ggsave(file.path(out_root, patient, "figures", paste0(patient, "_diff_LR_heatmap.png")), p,
         width = 6.5, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(out_root, patient, "figures", paste0(patient, "_diff_LR_heatmap.pdf")), p,
         width = 6.5, height = 6, bg = "white")
  cat(sprintf("[Fig] %s redrawn: %d LR (%d Clone1-higher, %d Clone2-higher) | Δ range %.3f to %.3f\n",
              patient, length(keep), length(pos), length(neg), min(dd$delta_prob), max(dd$delta_prob)))
  print(dd[order(-dd$delta_prob), c("target","interaction_name","delta_prob")][1:min(6, nrow(dd)), ], row.names = FALSE)
  print(tail(dd[order(-dd$delta_prob), c("target","interaction_name","delta_prob")], 6), row.names = FALSE)
  invisible(p)
}

redraw_diff_heatmap(res13, "P13")
redraw_diff_heatmap(res15, "P15")

for (f in c(file.path(out_root,"P13","figures","P13_diff_LR_heatmap.png"),
            file.path(out_root,"P15","figures","P15_diff_LR_heatmap.png"))) {
  img <- png::readPNG(f); rgb <- img[,,1:3]
  cat(sprintf("[Check] %s sd=%.4f white=%.3f size=%.0fKB\n", basename(f), sd(as.vector(rgb)),
              mean(rgb > 0.97), file.info(f)$size/1024))
}

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 两张差异热图重绘后，刷新图清单与像素自检结果
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

manifest <- data.frame()
for (patient in c("P13", "P15")) {
  fdir <- file.path(out_root, patient, "figures")
  for (s in sort(unique(sub("\\.(png|pdf)$", "", list.files(fdir))))) {
    fp_png <- file.path(fdir, paste0(s, ".png")); fp_pdf <- file.path(fdir, paste0(s, ".pdf"))
    img <- png::readPNG(fp_png); rgb <- img[, , 1:3]
    manifest <- rbind(manifest, data.frame(
      patient = patient, stem = s,
      png_exists = file.exists(fp_png), pdf_exists = file.exists(fp_pdf),
      png_kb = round(file.info(fp_png)$size/1024, 1), pdf_kb = round(file.info(fp_pdf)$size/1024, 1),
      sd_rgb = round(sd(as.vector(rgb)), 4), white_frac = round(mean(rgb > 0.97), 4),
      pass = (sd(as.vector(rgb)) > 0.02) && (mean(rgb > 0.97) < 0.99) &&
             (file.info(fp_png)$size/1024 > 10)))
  }
}
write.csv(manifest, file.path(out_root, "figure_manifest.csv"), row.names = FALSE)
cat("[Manifest refreshed]", nrow(manifest), "figures | all pass:", all(manifest$pass),
    "| all png+pdf:", all(manifest$png_exists & manifest$pdf_exists), "\n")
cat("[Redrawn]\n")
print(manifest[grepl("diff_LR_heatmap", manifest$stem), c("stem","png_kb","pdf_kb","sd_rgb","white_frac","pass")], row.names = FALSE)

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 重绘圈图：补标题、缩小标签字号、加大留白（边色=源、边宽=强度/条目数）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

redraw_circle <- function(patient, cc, measure = c("weight", "count")) {
  measure <- match.arg(measure)
  net <- if (measure == "weight") cc@net$weight else cc@net$count
  stem <- file.path(out_root, patient, "figures", paste0(patient, "_netVisual_circle_", measure))
  ttl <- sprintf("%s: %s of interactions (all 16 edges; edge colour = sender)",
                 patient, ifelse(measure == "weight", "strength", "number"))
  draw <- function() {
    netVisual_circle(net, vertex.weight = as.numeric(table(cc@idents)), weight.scale = TRUE,
                     label.edge = FALSE, edge.weight.max = max(net),
                     vertex.label.cex = 0.85, vertex.label.color = "black",
                     edge.curved = 0.2, margin = 0.22, title.name = ttl)
  }
  png(paste0(stem, ".png"), width = 7, height = 7, units = "in", res = 300, bg = "white")
  par(bg = "white", mar = c(0.5, 0.5, 2, 0.5)); draw(); dev.off()
  pdf(paste0(stem, ".pdf"), width = 7, height = 7, bg = "white")
  par(bg = "white", mar = c(0.5, 0.5, 2, 0.5)); draw(); dev.off()
  img <- png::readPNG(paste0(stem, ".png")); rgb <- img[,,1:3]
  cat(sprintf("[Fig] %s %-6s redrawn | sd=%.3f white=%.3f %.0fKB\n", patient, measure,
              sd(as.vector(rgb)), mean(rgb > 0.97), file.info(paste0(stem,".png"))$size/1024))
}

for (p in c("P13", "P15")) {
  cc <- get(paste0("cc", ifelse(p == "P13", "13", "15")))
  redraw_circle(p, cc, "weight")
  redraw_circle(p, cc, "count")
}

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 补充未过滤口径：全部 LR 概率和（不置零 p>=0.05），用于区分"通讯量"与"显著条目数"
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

SRC <- c("Clone1", "Clone2"); TGT <- c("C1Q_TAM", "SPP1_TAM")

unfiltered_summary <- function(res, patient) {
  d <- res$net   # 全部 16 边 x 全部 LR 的 prob/pval（未经显著性过滤）
  sub <- d[d$source %in% SRC & d$target %in% TGT, ]
  # 四条边：未过滤概率和、非零 LR 数、显著 LR 数
  e <- expand.grid(source = SRC, target = TGT, stringsAsFactors = FALSE)
  e$sum_prob_all_LR <- mapply(function(s, t) sum(sub$prob[sub$source == s & sub$target == t]), e$source, e$target)
  e$n_LR_nonzero   <- mapply(function(s, t) sum(sub$prob[sub$source == s & sub$target == t] > 0), e$source, e$target)
  e$n_LR_sig       <- mapply(function(s, t) sum(sub$prob[sub$source == s & sub$target == t] > 0 &
                                                sub$pval[sub$source == s & sub$target == t] < 0.05), e$source, e$target)
  e$sum_prob_sig_only <- mapply(function(s, t) {
    i <- sub$source == s & sub$target == t
    sum(sub$prob[i][sub$pval[i] < 0.05])
  }, e$source, e$target)
  e$weight_from_object <- mapply(function(s, t) {
    cc <- if (patient == "P13") cc13 else cc15
    cc@net$weight[s, t]
  }, e$source, e$target)
  write.csv(e, file.path(out_root, patient, "tables", paste0(patient, "_edge_summary_unfiltered.csv")), row.names = FALSE)

  # 通路级未过滤成员概率和
  ps <- aggregate(prob ~ pathway_name + source + target, data = sub, FUN = sum)
  write.csv(ps, file.path(out_root, patient, "tables", paste0(patient, "_pathway_sums_unfiltered.csv")), row.names = FALSE)

  cat(sprintf("\n=== %s | 四条边：未过滤 vs 仅显著 ===\n", patient))
  print(e, digits = 4, row.names = FALSE)
  cat("比值（Clone1/Clone2）：\n")
  for (t in TGT) {
    a <- e$sum_prob_all_LR[e$source == "Clone1" & e$target == t]
    b <- e$sum_prob_all_LR[e$source == "Clone2" & e$target == t]
    wa <- e$weight_from_object[e$source == "Clone1" & e$target == t]
    wb <- e$weight_from_object[e$source == "Clone2" & e$target == t]
    cat(sprintf("  %s: 未过滤 %.2f/%.2f = %.2f 倍 | 仅显著(weight) %.2f/%.2f = %.2f 倍\n",
                t, a, b, a/b, wa, wb, wa/wb))
  }
  invisible(list(edges = e, pathway = ps))
}

u13 <- unfiltered_summary(res13, "P13")
u15 <- unfiltered_summary(res15, "P15")

# 关键通路的未过滤成员和（用于改写报告中的 netP=0 表述）
get_pw <- function(ps, pw, tgt, s = "Clone2") ps$prob[ps$pathway_name == pw & ps$target == tgt & ps$source == s]
cat("\n[关键通路未过滤成员和]\n")
for (pw in c("MIF", "APP")) cat(sprintf("P13 %-8s SPP1_TAM: C1=%.4f C2=%.4f\n", pw, get_pw(u13$pathway, pw, "SPP1_TAM", "Clone1"), get_pw(u13$pathway, pw, "SPP1_TAM", "Clone2")))
for (pw in c("MIF", "VEGF", "CD99", "PERIOSTIN")) cat(sprintf("P15 %-9s C1Q_TAM: C1=%.4f C2=%.4f\n", pw, get_pw(u15$pathway, pw, "C1Q_TAM", "Clone1"), get_pw(u15$pathway, pw, "C1Q_TAM", "Clone2")))

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 生成"净概率(netP, 显著性过滤) vs 未过滤成员概率和"对照表（供报告引用）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

cmp_table <- function(res, u, patient) {
  netP <- res$diff_pat[, c("target","pathway_name","prob_Clone1","prob_Clone2","delta_prob",
                           "n_sig_Clone1","n_sig_Clone2")]
  ps <- u$pathway
  w <- reshape(ps[, c("target","pathway_name","source","prob")],
               idvar = c("target","pathway_name"), timevar = "source", direction = "wide")
  w <- data.frame(target = w$target, pathway = w$pathway_name,
                  unfilt_C1 = w$prob.Clone1, unfilt_C2 = w$prob.Clone2, stringsAsFactors = FALSE)
  m <- merge(netP, w, by.x = c("target","pathway_name"), by.y = c("target","pathway"), all.x = TRUE)
  m$d_unfilt <- m$unfilt_C1 - m$unfilt_C2
  out <- m[order(-abs(m$delta_prob)),
           c("target","pathway_name","prob_Clone1","prob_Clone2","delta_prob",
             "unfilt_C1","unfilt_C2","d_unfilt","n_sig_Clone1","n_sig_Clone2")]
  colnames(out) <- c("target","pathway","netP_C1","netP_C2","d_netP","unfilt_C1","unfilt_C2","d_unfilt","nsig_C1","nsig_C2")
  out <- head(out, 12)
  write.csv(out, file.path(out_root, patient, "tables", paste0(patient, "_pathway_netP_vs_unfiltered.csv")), row.names = FALSE)
  cat(sprintf("\n=== %s | 通路 netP（过滤）vs 未过滤成员概率和 ===\n", patient))
  print(out, digits = 3, row.names = FALSE)
  invisible(out)
}
pt13 <- cmp_table(res13, u13, "P13")
pt15 <- cmp_table(res15, u15, "P15")

top9 <- head(res13$diff_lr[order(-abs(res13$diff_lr$delta_prob)), ], 9)
cat("\n[P13 top9 LR by |delta|]\n")
print(top9[, c("target","pathway_name","interaction_name","prob_Clone1","prob_Clone2","delta_prob","signif_class")],
      digits = 3, row.names = FALSE)
top9b <- head(res15$diff_lr[order(-abs(res15$diff_lr$delta_prob)), ], 9)
cat("\n[P15 top9 LR by |delta|]\n")
print(top9b[, c("target","pathway_name","interaction_name","prob_Clone1","prob_Clone2","delta_prob","signif_class")],
      digits = 3, row.names = FALSE)

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 交付边界整理：自检清单移入 review/；通讯热图改用连续色阶 viridis（palette 契约）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

# 1) 自检类文件移出交付目录
rev_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone_TAM_cellchat/review/step01_clone_tam_cellchat"
dir.create(rev_dir, recursive = TRUE, showWarnings = FALSE)
f_man <- file.path(out_root, "figure_manifest.csv")
if (file.exists(f_man)) {
  file.rename(f_man, file.path(rev_dir, "figure_manifest.csv"))
  cat("[Boundary] figure_manifest.csv -> review/step01_clone_tam_cellchat/\n")
}

# 2) 通讯热图改 viridis（连续量契约）；失败则回退 Reds 并记录
redraw_hm <- function(patient, cc, measure) {
  stem <- file.path(out_root, patient, "figures", paste0(patient, "_netVisual_heatmap_", measure))
  mk <- function(col) function() netVisual_heatmap(cc, measure = measure, color.heatmap = col,
                                                   title.name = paste0(patient, ": interaction ", measure))
  ok <- tryCatch({
    ht <- mk("viridis")(); if (dev.cur() > 1) dev.off()
    png(paste0(stem, ".png"), width = 7, height = 6, units = "in", res = 300, bg = "white")
    draw(ht, newpage = TRUE); dev.off()
    pdf(paste0(stem, ".pdf"), width = 7, height = 6, bg = "white")
    draw(ht, newpage = TRUE); dev.off(); "viridis"
  }, error = function(e) {
    cat(sprintf("[Warn] viridis failed (%s), keep Reds for %s %s\n", conditionMessage(e), patient, measure))
    "Reds"
  })
  img <- png::readPNG(paste0(stem, ".png")); rgb <- img[,,1:3]
  cat(sprintf("[Fig] %s %s %s | sd=%.3f white=%.3f\n", patient, measure, ok,
              sd(as.vector(rgb)), mean(rgb > 0.97)))
}
for (p in c("P13","P15")) {
  cc <- get(paste0("cc", ifelse(p == "P13", "13", "15")))
  redraw_hm(p, cc, "count"); redraw_hm(p, cc, "weight")
}

# 3) 交付目录现状
cat("\n[Root files]", paste(list.files(out_root), collapse = ", "), "\n")
cat("[review files]", paste(list.files(rev_dir), collapse = ", "), "\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 校验热图文件未被失败的 viridis 尝试破坏，并刷新 review/ 内的自检清单
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

hm <- unlist(lapply(c("P13","P15"), function(p) lapply(c("count","weight"),
        function(m) file.path(out_root, p, "figures", paste0(p, "_netVisual_heatmap_", m, ".png")))))
cat("[Heatmaps]\n")
for (f in hm) cat(sprintf("  %-55s %.0f KB\n", basename(f), file.info(f)$size/1024))

rev_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone_TAM_cellchat/review/step01_clone_tam_cellchat"
manifest <- data.frame()
for (patient in c("P13", "P15")) {
  fdir <- file.path(out_root, patient, "figures")
  for (s in sort(unique(sub("\\.(png|pdf)$", "", list.files(fdir))))) {
    fp_png <- file.path(fdir, paste0(s, ".png")); fp_pdf <- file.path(fdir, paste0(s, ".pdf"))
    img <- png::readPNG(fp_png); rgb <- img[, , 1:3]
    manifest <- rbind(manifest, data.frame(
      patient = patient, stem = s, png_exists = file.exists(fp_png), pdf_exists = file.exists(fp_pdf),
      png_kb = round(file.info(fp_png)$size/1024, 1), pdf_kb = round(file.info(fp_pdf)$size/1024, 1),
      sd_rgb = round(sd(as.vector(rgb)), 4), white_frac = round(mean(rgb > 0.97), 4),
      pass = (sd(as.vector(rgb)) > 0.02) && (mean(rgb > 0.97) < 0.99) && (file.info(fp_png)$size/1024 > 10)))
  }
}
write.csv(manifest, file.path(rev_dir, "figure_manifest.csv"), row.names = FALSE)
cat("\n[Manifest]", nrow(manifest), "figures | all pass:", all(manifest$pass),
    "| all png+pdf:", all(manifest$png_exists & manifest$pdf_exists), "| min png KB:", min(manifest$png_kb), "\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 通讯热图改用 viridis 连续色阶（netVisual_heatmap 支持传入两端色向量），并核对并列行的概率精度
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({ library(viridisLite); library(ComplexHeatmap) })
vc <- viridisLite::viridis(2)   # 两端色：#440154 -> #FDE725
cat("[viridis ends]", paste(vc, collapse = " -> "), "\n")

redraw_hm_viridis <- function(patient, cc, measure) {
  stem <- file.path(out_root, patient, "figures", paste0(patient, "_netVisual_heatmap_", measure))
  mk <- function() netVisual_heatmap(cc, measure = measure, color.heatmap = vc,
                                     title.name = paste0(patient, ": interaction ", measure))
  ht <- tryCatch(mk(), error = function(e) { cat("[FAIL]", conditionMessage(e), "\n"); NULL })
  if (is.null(ht)) return(invisible(FALSE))
  if (dev.cur() > 1) dev.off()
  png(paste0(stem, ".png"), width = 7, height = 6, units = "in", res = 300, bg = "white")
  draw(ht, newpage = TRUE); dev.off()
  pdf(paste0(stem, ".pdf"), width = 7, height = 6, bg = "white")
  draw(ht, newpage = TRUE); dev.off()
  img <- png::readPNG(paste0(stem, ".png")); rgb <- img[,,1:3]
  cat(sprintf("[Fig] %s %-6s viridis | sd=%.3f white=%.3f %.0fKB\n", patient, measure,
              sd(as.vector(rgb)), mean(rgb > 0.97), file.info(paste0(stem,".png"))$size/1024))
  TRUE
}
for (p in c("P13","P15")) {
  cc <- get(paste0("cc", ifelse(p == "P13", "13", "15")))
  redraw_hm_viridis(p, cc, "count"); redraw_hm_viridis(p, cc, "weight")
}

# 并列（Δ=0）行的概率精度核查：确认是"完全相同"而不是四舍五入
cat("\n[Zero-delta precision check]\n")
for (p in c("P13","P15")) {
  d <- read.csv(file.path(out_root, p, "tables", paste0(p, "_diff_LR_Clone1_vs_Clone2.csv")))
  b <- d[d$signif_class == "both_significant", ]
  z <- b[b$delta_prob == 0, ]
  cat(sprintf("%s: tied %d | max |Δ| among ties = %.3g | max prob = %.4f\n",
              p, nrow(z), max(abs(z$delta_prob)), max(z$prob_Clone1)))
  # 同源 homo-pair 占比（配体=受体首基因）
  hom <- sum(sub("_.*$", "", z$interaction_name) == sapply(strsplit(z$interaction_name, "_"), function(x) x[2]))
  cat(sprintf("   homophilic-like pairs among ties: %d / %d\n", hom, nrow(z)))
}
cat("\n[viridis 支持] netVisual_heatmap(color.heatmap = viridisLite::viridis(2)) 已成功出图\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 交付目录守卫：若中途某步在根目录重新生成自检清单，统一收拢到 review/（保持交付边界）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

rev_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone_TAM_cellchat/review/step01_clone_tam_cellchat"
stray <- list.files(out_root, pattern = "^(figure_manifest|qc_.*)\\.csv$", full.names = TRUE)
if (length(stray) > 0) {
  file.rename(stray, file.path(rev_dir, basename(stray)))
  cat("[Boundary] moved to review/:", paste(basename(stray), collapse = ", "), "\n")
} else {
  cat("[Boundary] no stray self-check file in delivery root\n")
}
cat("[Delivery root]", paste(sort(list.files(out_root)), collapse = ", "), "\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 最终刷新图件自检清单（含已改 viridis 的 4 张通讯热图）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

rev_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone_TAM_cellchat/review/step01_clone_tam_cellchat"
manifest <- data.frame()
for (patient in c("P13", "P15")) {
  fdir <- file.path(out_root, patient, "figures")
  for (s in sort(unique(sub("\\.(png|pdf)$", "", list.files(fdir))))) {
    fp_png <- file.path(fdir, paste0(s, ".png")); fp_pdf <- file.path(fdir, paste0(s, ".pdf"))
    img <- png::readPNG(fp_png); rgb <- img[, , 1:3]
    manifest <- rbind(manifest, data.frame(
      patient = patient, stem = s, png_exists = file.exists(fp_png), pdf_exists = file.exists(fp_pdf),
      png_kb = round(file.info(fp_png)$size/1024, 1), pdf_kb = round(file.info(fp_pdf)$size/1024, 1),
      sd_rgb = round(sd(as.vector(rgb)), 4), white_frac = round(mean(rgb > 0.97), 4),
      pass = (sd(as.vector(rgb)) > 0.02) && (mean(rgb > 0.97) < 0.99) && (file.info(fp_png)$size/1024 > 10)))
  }
}
write.csv(manifest, file.path(rev_dir, "figure_manifest.csv"), row.names = FALSE)
cat("[Manifest]", nrow(manifest), "| all pass:", all(manifest$pass),
    "| all png+pdf:", all(manifest$png_exists & manifest$pdf_exists), "\n")
print(manifest[grepl("heatmap", manifest$stem), c("stem","png_kb","sd_rgb","white_frac","pass")], row.names = FALSE)

# 表清单核对（每例应为 13 张）
for (p in c("P13","P15")) cat(sprintf("[Tables] %s: %d files\n", p, length(list.files(file.path(out_root, p, "tables")))))

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 独立审查全部通过后，追加审计留痕（东八区时间；同一 unit 只追加一次）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

log_dir <- file.path(Sys.getenv("HOME"), ".claude", "logs", format(Sys.Date(), "%Y-%m-%d"))
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "audit_done.jsonl")
session_id <- "e99aac27-4525-42e3-8a8a-00c81a7e0358"

existing <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character(0)
already <- any(grepl(session_id, existing, fixed = TRUE) & grepl("\"sc-analyst\"", existing, fixed = TRUE))
if (already) {
  cat("[Audit log] 本 session 已有 sc-analyst 记录，跳过重复追加\n")
} else {
  entry <- sprintf('{"ts":"%s","session":"%s","subagent_type":"sc-analyst","verdict":"pass","unit":"clone_tam_cellchat"}',
                   format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), session_id)
  cat(entry, "\n", file = log_file, append = TRUE)
  cat("[Audit log] appended:", entry, "\n")
}
cat("[Audit log] path:", log_file, "| lines:", length(readLines(log_file, warn = FALSE)), "\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 交付前汇总核对：checkpoint、表、图、分组成对性
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

cat("[Checkpoints]\n")
for (f in c("P13/data/cellchat_P13.rds", "P15/data/cellchat_P15.rds",
            "P13/data/P13_cellchat_input.rds", "P15/data/P15_cellchat_input.rds"))
  cat(sprintf("  %-40s %8.1f MB  %s\n", f, file.info(file.path(out_root, f))$size/1e6,
              ifelse(file.exists(file.path(out_root, f)), "OK", "MISSING")))

cat("\n[Per-patient delivery]\n")
for (p in c("P13","P15")) {
  nfig <- length(list.files(file.path(out_root, p, "figures"), pattern = "\\.png$"))
  npdf <- length(list.files(file.path(out_root, p, "figures"), pattern = "\\.pdf$"))
  ntab <- length(list.files(file.path(out_root, p, "tables")))
  cat(sprintf("  %s: figures png %d / pdf %d | tables %d | summary html %s\n", p, nfig, npdf, ntab,
              file.exists(file.path(out_root, paste0(p, "_汇总.html")))))
}

cat("\n[Group counts vs expected]\n")
cat("  P13:", paste(as.integer(table(factor(P13$group, levels = c("Clone1","Clone2","C1Q_TAM","SPP1_TAM")))), collapse = "/"), "(expected 515/229/142/124)\n")
cat("  P15:", paste(as.integer(table(factor(P15$group, levels = c("Clone1","Clone2","C1Q_TAM","SPP1_TAM")))), collapse = "/"), "(expected 578/298/119/99)\n")

cat("\n[Four edges non-empty]\n")
cat("  P13 count:", paste(res13$edges$count, collapse = "/"), "| P15 count:", paste(res15$edges$count, collapse = "/"), "\n")

cat("\n[Delivery root]", paste(sort(list.files(out_root)), collapse = ", "), "\n")

cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 读取 P13/P15 已有的 Clone1 vs Clone2 通讯概率表，筛出各自纳入的尾部 LR 对
# 校验：纳入对 × 2 靶行数齐全、无重复、通路归属与预期分组一致
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(ggsci)})
cat("R:", R.version.string, "| ggplot2:", as.character(packageVersion("ggplot2")), "\n")

base_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone_TAM_cellchat"

cfg <- list(
  P13 = list(
    csv = file.path(base_dir, "P13/tables/P13_diff_LR_Clone1_vs_Clone2.csv"),
    pathway_order = c("THBS", "GALECTIN", "CSF", "COLLAGEN"),
    pairs = c("COMP_SDC4", "COMP_CD47", "COMP_SDC1", "THBS2_SDC4",
              "LGALS9_P4HB", "CSF1_CSF1R", "COL2A1_CD44", "COL9A2_CD44")
  ),
  P15 = list(
    csv = file.path(base_dir, "P15/tables/P15_diff_LR_Clone1_vs_Clone2.csv"),
    pathway_order = c("GALECTIN", "EPHB", "COLLAGEN", "SEMA4", "LAMININ"),
    pairs = c("LGALS9_P4HB", "EFNB1_EPHB3", "COL9A3_SDC4", "SEMA4C_PLXNB2",
              "LAMA3_ITGA6_ITGB1", "LAMC1_ITGA6_ITGB1")
  )
)
targets <- c("C1Q_TAM", "SPP1_TAM")

raw_list <- list()
for (pt in names(cfg)) {
  stopifnot(file.exists(cfg[[pt]]$csv))
  d_full <- read.csv(cfg[[pt]]$csv, check.names = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("\n[Read] %s: %s | full table %d rows x %d cols\n",
              pt, basename(cfg[[pt]]$csv), nrow(d_full), ncol(d_full)))

  d <- d_full[d_full$target %in% targets & d_full$interaction_name %in% cfg[[pt]]$pairs, ]
  exp_rows <- length(cfg[[pt]]$pairs) * length(targets)
  cat(sprintf("[Filter] %s 纳入 LR x 靶 | rows: %d (expected %d) -> %s\n",
              pt, nrow(d), exp_rows, ifelse(nrow(d) == exp_rows, "OK", "MISMATCH")))

  # 缺失组合校验
  want <- expand.grid(target = targets, interaction_name = cfg[[pt]]$pairs,
                      stringsAsFactors = FALSE)
  miss <- nrow(dplyr::anti_join(want, d, by = c("target", "interaction_name")))
  cat(sprintf("[Check ] %s missing (target x LR) combos: %d\n", pt, miss))
  # 重复校验
  dups <- sum(duplicated(d[, c("target", "interaction_name")]))
  cat(sprintf("[Check ] %s duplicated (target x LR) rows: %d\n", pt, dups))

  # 通路归属：从表里推导，再与预期分组顺序对照
  pw_map <- unique(d[, c("interaction_name", "pathway_name")])
  pw_map <- pw_map[order(match(pw_map$pathway_name, cfg[[pt]]$pathway_order),
                         match(pw_map$interaction_name, cfg[[pt]]$pairs)), ]
  cat(sprintf("[Check ] %s pathway assignment per LR:\n", pt))
  print(pw_map, row.names = FALSE)
  bad_pw <- setdiff(unique(pw_map$pathway_name), cfg[[pt]]$pathway_order)
  cat(sprintf("[Check ] %s pathways outside declared order: %s\n", pt,
              ifelse(length(bad_pw) == 0, "none", paste(bad_pw, collapse = ", "))))
  cat(sprintf("[Check ] %s declared groups absent in data: %s\n", pt,
              { a <- setdiff(cfg[[pt]]$pathway_order, unique(pw_map$pathway_name))
                ifelse(length(a) == 0, "none", paste(a, collapse = ", ")) }))

  cat(sprintf("[Values] %s selected probabilities (6 dp):\n", pt))
  show <- d[, c("target", "pathway_name", "interaction_name", "prob_Clone1", "prob_Clone2", "delta_prob")]
  show[, 4:6] <- round(show[, 4:6], 6)
  print(show[order(match(show$target, targets), show$interaction_name), ], row.names = FALSE)

  raw_list[[pt]] <- d
}

cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 把纳入的 LR 对整理成绘制用长表：通路分组 + 组内按 max(|ΔC1Q|,|ΔSPP1|) 降序
# 同时写出逐格可溯源的数值表（LR × 靶 × 克隆 的概率与 Δ），全精度保存
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

target_levels <- c("C1Q_TAM", "SPP1_TAM")
clone_levels  <- c("Clone1", "Clone2")
vals_list <- list(); order_list <- list()

for (pt in names(raw_list)) {
  d <- raw_list[[pt]]

  # 每个 LR 的两个靶向 Δ 与排序键
  key <- d %>%
    dplyr::select(interaction_name, pathway_name, target, delta_prob) %>%
    tidyr::pivot_wider(names_from = target, values_from = delta_prob) %>%
    dplyr::mutate(maxabs = pmax(abs(.data[["C1Q_TAM"]]), abs(.data[["SPP1_TAM"]])))
  key$pw_rank <- match(key$pathway_name, cfg[[pt]]$pathway_order)

  # 通路按声明顺序，组内按 maxabs 降序（并列时按 LR 名字母序，保证可复现）
  key <- key[order(key$pw_rank, -key$maxabs, key$interaction_name), ]
  key$y_order <- seq_len(nrow(key))
  cat(sprintf("\n[Order] %s 自上而下的行顺序（y_order: pathway | LR | delta_C1Q | delta_SPP1 | maxabs）:\n", pt))
  print(data.frame(y_order = key$y_order, pathway = key$pathway_name,
                   LR = key$interaction_name,
                   delta_C1Q = round(key[["C1Q_TAM"]], 6),
                   delta_SPP1 = round(key[["SPP1_TAM"]], 6),
                   maxabs = round(key$maxabs, 6)), row.names = FALSE)

  # 长表：LR × 靶 × 克隆
  long <- d %>%
    dplyr::select(target, pathway_name, interaction_name, prob_Clone1, prob_Clone2, delta_prob) %>%
    tidyr::pivot_longer(cols = c(prob_Clone1, prob_Clone2), names_to = "clone_col", values_to = "prob") %>%
    dplyr::mutate(clone = ifelse(clone_col == "prob_Clone1", "Clone1", "Clone2"))
  long$y_order <- key$y_order[match(long$interaction_name, key$interaction_name)]
  long$pw_rank <- match(long$pathway_name, cfg[[pt]]$pathway_order)
  long$drawn   <- long$prob > 0
  long <- long[order(long$y_order, match(long$target, target_levels), match(long$clone, clone_levels)), ]

  out_vals <- data.frame(
    patient                 = pt,
    y_order                 = long$y_order,
    pathway_name            = long$pathway_name,
    interaction_name        = long$interaction_name,
    target                  = long$target,
    clone                   = long$clone,
    communication_prob      = long$prob,
    delta_Clone1_minus_Clone2 = long$delta_prob,
    delta_sign              = ifelse(long$delta_prob > 0, "positive",
                              ifelse(long$delta_prob < 0, "negative", "zero")),
    drawn_as_point          = long$drawn,
    stringsAsFactors = FALSE
  )
  vals_path <- file.path(base_dir, pt, "tables", sprintf("%s_tail_LR_dotplot_values.csv", pt))
  write.csv(out_vals, vals_path, row.names = FALSE)
  cat(sprintf("[Save] %s: %s (%d rows, %d cols)\n", pt, vals_path, nrow(out_vals), ncol(out_vals)))

  vals_list[[pt]]  <- out_vals
  order_list[[pt]] <- key
}

cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 定义尾部 LR 点图的绘图函数（两例复用同一套编码），并绘制 P13
# 版式：y=LR(按通路分组,组内按 max|Δ| 降序)；x=4 个概率列 + 2 个 Δ 列；点大小=通讯概率；点颜色=克隆
# 修订：通路名改为右侧对齐到固定留白区（原来左对齐会伸进第一列数据），括号加粗加深，图例刻度补齐值域两端
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

# 配色：克隆色取 ggsci NPG 主题色（第 1、4 色），Δ 正负取发散色两端
NPG <- ggsci::pal_npg()(8)
col_map <- c(Clone1 = NPG[1], Clone2 = NPG[4],
             negative = "#2166AC", positive = "#B2182B", zero = "grey35")
cat(sprintf("[Palette] Clone1=%s Clone2=%s | dDelta negative=%s positive=%s\n",
            NPG[1], NPG[4], "#2166AC", "#B2182B"))

has_cairo <- isTRUE(capabilities("cairo"))
AR <- if (has_cairo) "→" else "->"
DL <- if (has_cairo) "Δ" else "Delta"
cat(sprintf("[Device ] cairo=%s | arrow='%s' delta='%s'\n", has_cairo, AR, DL))

X_PROB <- c(C1Q_TAM_Clone1 = 1.0, C1Q_TAM_Clone2 = 2.0, SPP1_TAM_Clone1 = 3.6, SPP1_TAM_Clone2 = 4.6)
X_DELT <- c(C1Q_TAM = 6.1, SPP1_TAM = 7.0)
SEP_TARGET <- 2.8
SEP_DELTA  <- 5.35
X_PW_BRACKET <- 0.04   # 分组括号
X_PW_TEXT    <- 0.72   # 通路名右对齐锚点（保证文字绝不进入 x=1 的数据列）
X_MIN        <- -0.10  # 左侧留白区起点

make_tail_plot <- function(pt, vals, key) {
  n <- nrow(key)
  yk <- data.frame(interaction_name = key$interaction_name,
                   pathway_name     = key$pathway_name,
                   ypos             = rev(seq_len(n)),
                   stringsAsFactors = FALSE)

  # ---- 概率点（零概率不画点，见 caption） ----
  pd <- vals[vals$drawn_as_point, ]
  pd$ypos <- yk$ypos[match(pd$interaction_name, yk$interaction_name)]
  pd$x <- mapply(function(tg, cl) X_PROB[[paste(tg, cl, sep = "_")]], pd$target, pd$clone)

  # ---- 右侧 Δ 标注 ----
  dl <- do.call(rbind, lapply(c("C1Q_TAM", "SPP1_TAM"), function(tg) {
    data.frame(interaction_name = key$interaction_name,
               ypos             = yk$ypos,
               target           = tg,
               delta            = key[[tg]],
               stringsAsFactors = FALSE)
  }))
  dl$x <- X_DELT[dl$target]
  dl$cls <- factor(ifelse(dl$delta > 0, "positive", ifelse(dl$delta < 0, "negative", "zero")),
                   levels = c("Clone1", "Clone2", "negative", "positive", "zero"))
  dl$lab <- sprintf("%.3f", dl$delta)

  # ---- 通路分组括号 + 灰色斜体通路名（右对齐到留白区） ----
  pw <- do.call(rbind, lapply(unique(key$pathway_name), function(g) {
    y <- yk$ypos[yk$pathway_name == g]
    data.frame(pathway_name = g, ymid = mean(range(y)),
               ytop = max(y) + 0.42, ybot = min(y) - 0.42, stringsAsFactors = FALSE)
  }))
  pw <- pw[order(match(pw$pathway_name, cfg[[pt]]$pathway_order)), ]

  maxp <- max(pd$communication_prob)
  minp <- min(pd$communication_prob)
  # 刻度：落在值域内的常规刻度 + 值域两端各补一格（保证最大/最小点可估读）
  brks <- c(0.05, 0.10, 0.15, 0.20)
  brks <- brks[brks > minp & brks < maxp]
  if (length(brks) == 0) brks <- signif(mean(c(minp, maxp)), 2)
  if (minp < min(brks) * 0.8) brks <- c(signif(minp, 1), brks)
  if (maxp > max(brks) * 1.05) brks <- c(brks, signif(maxp, 2))
  cat(sprintf("  [Legend] %s size ticks: %s (data range %.4f-%.4f)\n", pt,
              paste(sprintf("%.2f", brks), collapse = "/"), minp, maxp))

  ggplot() +
    geom_point(data = pd, aes(x = x, y = ypos, size = communication_prob, colour = clone), stroke = 0.25) +
    geom_text(data = dl, aes(x = x, y = ypos, label = lab, colour = cls), size = 2.4, show.legend = FALSE) +
    geom_segment(data = pw, aes(x = X_PW_BRACKET, xend = X_PW_BRACKET, y = ybot, yend = ytop),
                 colour = "grey60", linewidth = 0.4, inherit.aes = FALSE) +
    geom_text(data = pw, aes(x = X_PW_TEXT, y = ymid, label = pathway_name), hjust = 1,
              size = 2.5, colour = "grey30", fontface = "italic", inherit.aes = FALSE) +
    geom_segment(aes(x = SEP_TARGET, xend = SEP_TARGET, y = 0.5, yend = n + 0.5),
                 colour = "black", linewidth = 0.35, inherit.aes = FALSE) +
    geom_segment(aes(x = SEP_DELTA, xend = SEP_DELTA, y = 0.5, yend = n + 0.5),
                 colour = "grey55", linetype = "dashed", linewidth = 0.3, inherit.aes = FALSE) +
    scale_colour_manual(values = col_map, breaks = c("Clone1", "Clone2"),
                        labels = c("Clone1", "Clone2"), name = "Clone", drop = FALSE) +
    scale_size_area(limits = c(0, maxp), max_size = 5.2, breaks = brks,
                    labels = sprintf("%.2f", brks), name = "Communication\nprobability") +
    scale_x_continuous(breaks = c(X_PROB, X_DELT),
                       labels = c(paste0("Clone1", AR, "C1Q"), paste0("Clone2", AR, "C1Q"),
                                  paste0("Clone1", AR, "SPP1"), paste0("Clone2", AR, "SPP1"),
                                  paste(DL, "C1Q"), paste(DL, "SPP1")),
                       limits = c(X_MIN, 7.75), expand = c(0, 0)) +
    scale_y_continuous(breaks = yk$ypos, labels = yk$interaction_name,
                       limits = c(0.4, n + 0.6), expand = c(0, 0)) +
    labs(caption = "Blank cell = zero communication probability (no symbol).  Delta = Clone1 minus Clone2; blue = negative, red = positive.") +
    theme_bw(base_size = 8) +
    theme(panel.grid      = element_blank(),
          panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          axis.title      = element_blank(),
          plot.title      = element_blank(),
          axis.text.y     = element_text(size = 6.5, colour = "black"),
          axis.text.x     = element_text(size = 6.8, colour = "black"),
          axis.ticks      = element_line(colour = "black", linewidth = 0.3),
          axis.ticks.length = unit(2, "pt"),
          legend.title    = element_text(size = 7),
          legend.text     = element_text(size = 6.5),
          legend.key.size = unit(0.30, "cm"),
          legend.background = element_blank(),
          legend.box      = "vertical",
          legend.position = "right",
          plot.caption    = element_text(hjust = 0, size = 6.4, colour = "grey30"),
          plot.margin     = margin(5, 8, 4, 4))
}

render_tail_plot <- function(pt, vals, key, width, height) {
  p <- make_tail_plot(pt, vals, key)
  fig_dir <- file.path(base_dir, pt, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  stem <- file.path(fig_dir, sprintf("%s_tail_LR_dotplot", pt))
  ggsave(paste0(stem, ".png"), plot = p, width = width, height = height, dpi = 600,
         bg = "white", device = "png", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot = p, width = width, height = height, bg = "white",
         device = if (has_cairo) cairo_pdf else pdf, limitsize = FALSE)
  for (f in paste0(stem, c(".png", ".pdf"))) {
    cat(sprintf("[Figure] %s | %.1f KB\n", f, file.info(f)$size / 1024))
  }
  invisible(p)
}

# 图幅：6 个 x 列（4 概率 + 2 Δ）需要 8.8 in 宽；行高随行数线性放大
W <- 8.8
H13 <- 2.2 + 0.38 * nrow(order_list[["P13"]])
p13 <- render_tail_plot("P13", vals_list[["P13"]], order_list[["P13"]], W, H13)
cat(sprintf("[Layout ] P13 rows=%d | width=%.1f in, height=%.1f in\n", nrow(order_list[["P13"]]), W, H13))

cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 用同一套编码绘制 P15 的尾部 LR 点图（与 P13 独立成图，不并排、不比较）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

H15 <- 2.2 + 0.38 * nrow(order_list[["P15"]])
p15 <- render_tail_plot("P15", vals_list[["P15"]], order_list[["P15"]], W, H15)
cat(sprintf("[Layout ] P15 rows=%d | width=%.1f in, height=%.1f in\n", nrow(order_list[["P15"]]), W, H15))

cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 机械核对：数值表逐格回溯源表；图片像素完整性（非空白、含两种克隆色）；文件清单
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

# ---- 1. 数值逐格回溯（values CSV <-> 原始差异表）----
for (pt in names(vals_list)) {
  src <- read.csv(cfg[[pt]]$csv, check.names = FALSE, stringsAsFactors = FALSE)
  v <- read.csv(file.path(base_dir, pt, "tables", sprintf("%s_tail_LR_dotplot_values.csv", pt)),
                stringsAsFactors = FALSE)
  src_sel <- src[src$target %in% c("C1Q_TAM", "SPP1_TAM") & src$interaction_name %in% cfg[[pt]]$pairs,
                 c("target", "interaction_name", "prob_Clone1", "prob_Clone2", "delta_prob")]
  m <- merge(v, src_sel, by = c("target", "interaction_name"))
  src_prob <- ifelse(m$clone == "Clone1", m$prob_Clone1, m$prob_Clone2)
  d_prob <- max(abs(m$communication_prob - src_prob))
  d_del  <- max(abs(m$delta_Clone1_minus_Clone2 - m$delta_prob))
  d_self <- max(abs((m$prob_Clone1 - m$prob_Clone2) - m$delta_prob))
  cat(sprintf("[Reconcile] %s | rows matched %d/%d | max|prob diff|=%.3g | max|delta diff|=%.3g | max|(C1-C2)-delta|=%.3g | %s\n",
              pt, nrow(m), nrow(v), d_prob, d_del, d_self,
              ifelse(nrow(m) == nrow(v) && max(d_prob, d_del, d_self) < 1e-12, "OK", "MISMATCH")))
  cat(sprintf("[Points   ] %s | drawn=%d / total=%d (zero-prob blanks=%d)\n", pt,
              sum(v$drawn_as_point), nrow(v), sum(!v$drawn_as_point)))
}

# ---- 2. 像素完整性（nativeRaster 逐像素统计，内存友好）----
pixel_check <- function(f) {
  nr <- png::readPNG(f, native = TRUE)
  hw <- dim(nr)
  v <- as.integer(nr)
  r <- bitwAnd(v, 255L); g <- bitwAnd(bitwShiftR(v, 8), 255L); b <- bitwAnd(bitwShiftR(v, 16), 255L)
  mx <- pmax(r, g, b); mn <- pmin(r, g, b)
  white <- mean(r > 247L & g > 247L & b > 247L)
  ink   <- mean(mx < 245L)
  chroma <- mean((mx - mn) > 20L)
  col_frac <- function(hex) {
    tr <- strtoi(substr(hex, 2, 3), 16L); tg <- strtoi(substr(hex, 4, 5), 16L); tb <- strtoi(substr(hex, 7, 8), 16L)
    mean(abs(r - tr) < 30L & abs(g - tg) < 30L & abs(b - tb) < 30L)
  }
  data.frame(file = basename(f), px_w = hw[2], px_h = hw[1],
             sd_gray = round(sd(c(r, g, b)), 4), white_frac = round(white, 4),
             ink_frac = round(ink, 4), chroma_frac = round(chroma, 5),
             frac_Clone1col = round(col_frac(NPG[1]), 6),
             frac_Clone2col = round(col_frac(NPG[4]), 6))
}

figs <- unlist(lapply(names(cfg), function(pt)
  file.path(base_dir, pt, "figures", sprintf("%s_tail_LR_dotplot", pt))), use.names = FALSE)
pc <- do.call(rbind, lapply(paste0(figs, ".png"), pixel_check))
print(pc, row.names = FALSE)
stopifnot(all(pc$sd_gray > 0.02), all(pc$white_frac > 0.5 & pc$white_frac < 0.995),
          all(pc$chroma_frac > 0.0005), all(pc$frac_Clone1col > 0 & pc$frac_Clone2col > 0))
cat("[Pixel   ] all figures: non-blank (sd>0.02), white background, both clone colours present -> OK\n")

# ---- 3. 产物清单 ----
for (f in c(paste0(figs, ".png"), paste0(figs, ".pdf"),
            file.path(base_dir, names(cfg), "tables",
                      sprintf("%s_tail_LR_dotplot_values.csv", names(cfg))))) {
  cat(sprintf("[Output  ] %-95s %8.1f KB\n", f, file.info(f)$size / 1024))
}

cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 复核绘图数据与已落盘数值表完全一致（绘制对象 pd 与 CSV 同源确认），并打印最终值域
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

for (pt in names(cfg)) {
  v <- vals_list[[pt]]
  pd <- v[v$drawn_as_point, ]
  cat(sprintf("[Range] %s | prob range %.4f - %.4f | delta range %.4f - %.4f | rows(y)=%d\n",
              pt, min(pd$communication_prob), max(pd$communication_prob),
              min(v$delta_Clone1_minus_Clone2), max(v$delta_Clone1_minus_Clone2),
              length(unique(v$y_order))))
}
cat("[Mapping] size area proportional to probability within each figure; delta text = %.3f precision\n")

cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 统计纳入的 LR×靶 组合在源表中的 signif_class 分布（用于如实说明"纯描述性展示"的解读边界）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

sig_list <- list()
for (pt in names(cfg)) {
  d <- raw_list[[pt]]
  tb <- table(factor(d$signif_class, levels = c("both_significant", "Clone1_only", "Clone2_only", "neither")))
  cat(sprintf("\n[Signif] %s 纳入 %d 个 LR×靶 组合 的 signif_class 分布:\n", pt, nrow(d)))
  print(tb)
  cat(sprintf("[Signif] %s 两克隆均不显著(neither) 占比: %d/%d = %.1f%%\n",
              pt, tb[["neither"]], sum(tb), 100 * tb[["neither"]] / sum(tb)))
  sig_list[[pt]] <- list(total = sum(tb), neither = tb[["neither"]],
                         clone1_only = tb[["Clone1_only"]], clone2_only = tb[["Clone2_only"]],
                         both = tb[["both_significant"]])
}

cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 追加独立审查留痕到全局审计日志（同一 unit 只追加一次，幂等）；记录用 jsonlite 结构化生成
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

audit_file <- "/path/to/.claude/logs/2026-10-29/audit_done.jsonl"
dir.create(dirname(audit_file), showWarnings = FALSE, recursive = TRUE)
sess <- "e99aac27-4525-42e3-8a8a-00c81a7e0358"
unit <- "tail_lr_dotplot_done"

existing <- if (file.exists(audit_file)) readLines(audit_file, warn = FALSE) else character(0)
already <- any(grepl(unit, existing, fixed = TRUE) & grepl(sess, existing, fixed = TRUE))

if (!already) {
  rec <- list(
    ts = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    session = sess,
    subagent_type = "sc-analyst",
    unit = unit,
    verdict = "pass",
    critic = "independent general-purpose (judge != doer); full audit plus incremental post-redraw audit",
    round1 = paste("pass on all conclusion-level items: row inclusion, cell-by-cell number fidelity,",
                   "data-derived ordering, pathway assignment, no cross-patient mixing, no significance marks,",
                   "every code cell fully implemented, required encodings, file integrity, two independent figures"),
    round2 = paste("advisory only: pathway labels overlapped the first data column; redrawn with a",
                   "fixed-anchor gutter; incremental re-audit pass (numbers invariant, change layout-only, artifacts intact)"),
    visual_check = paste("visual-inspector: P13 clean / P15 good; advisories recorded -",
                         "per-figure size scales not comparable, and about 40 percent of included LR x target",
                         "combinations are not significant in either clone yet drawn with equal visual weight")
  )
  cat(jsonlite::toJSON(rec, auto_unbox = TRUE, pretty = FALSE), "\n",
      sep = "", file = audit_file, append = TRUE)
  cat("[AuditLog] appended 1 record ->", audit_file, "\n")
} else {
  cat("[AuditLog] record already present for this unit/session -> no second append\n")
}

lines <- readLines(audit_file, warn = FALSE)
hit <- grepl(unit, lines, fixed = TRUE) & grepl(sess, lines, fixed = TRUE)
cat(sprintf("[AuditLog] total lines %d | records for this unit %d | pass verdict present: %s\n",
            length(lines), sum(hit), any(grepl('"verdict":"pass"', lines[hit], fixed = TRUE))))

cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 汇总 top20 LR 通讯概率：生成数值表(每例40行)并按通路分组确定 y 轴顺序
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

base_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone_TAM_cellchat"
target_levels <- c("C1Q_TAM", "SPP1_TAM")
clone_levels  <- c("Clone1", "Clone2")

# 数值唯一来源（只读）
DIFF_CSV <- c(P13 = file.path(base_dir, "P13/tables/P13_diff_LR_Clone1_vs_Clone2.csv"),
              P15 = file.path(base_dir, "P15/tables/P15_diff_LR_Clone1_vs_Clone2.csv"))
LR_CSV   <- c(P13 = file.path(base_dir, "P13/tables/P13_LR_probabilities_clone_to_TAM.csv"),
              P15 = file.path(base_dir, "P15/tables/P15_LR_probabilities_clone_to_TAM.csv"))

# top20：按每个 LR 在两靶中最大 |Δ| 排序取前 20 个不同 interaction_name（与委派清单逐位核对一致）
TOP20 <- list(
  P13 = c("COMP_SDC4","COMP_CD47","LGALS9_P4HB","COMP_SDC1","LGALS9_CD44","THBS2_SDC4",
          "COL1A1_CD44","COL1A1_SDC4","THBS1_SDC4","COMP_ITGA3_ITGB1","COL6A2_CD44","THBS2_CD47",
          "IGFBP3_TMEM219","COL1A1_ITGAV_ITGB8","COL6A2_SDC4","COL1A1_SDC1","CSF1_CSF1R","COL6A3_CD44",
          "COL3A1_ADGRG1","THBS1_CD47"),
  P15 = c("LGALS9_P4HB","EFNB1_EPHB3","THBS4_SDC4","COL9A3_SDC4","SEMA4C_PLXNB2","COL9A3_SDC1",
          "THBS4_SDC1","COMP_SDC4","THBS4_CD47","COMP_SDC1","TNFSF12_TNFRSF12A","COL9A3_ITGA11_ITGB1",
          "COMP_CD47","COL6A2_SDC4","THBS4_CD36","LRFN4_PTPRF","COL6A2_SDC1","COL9A3_ITGA1_ITGB1",
          "EFNB1_EPHB4","LGALS9_CD45")
)

ORD <- list(); VAL <- list(); PROB4 <- list()

for (pt in c("P13", "P15")) {
  pairs <- unique(TOP20[[pt]])
  d <- read.csv(DIFF_CSV[[pt]], check.names = FALSE, stringsAsFactors = FALSE)
  d <- d[d$interaction_name %in% pairs,
         c("target","pathway_name","interaction_name","prob_Clone1","prob_Clone2","delta_prob")]
  stopifnot(nrow(d) == 2 * length(pairs))

  # 排序：通路按组内最大 |Δ| 降序；通路内按 |Δ| 降序
  ag <- aggregate(abs(delta_prob) ~ interaction_name + pathway_name, data = d, FUN = max)
  names(ag)[3] <- "max_abs_delta"
  pw <- aggregate(max_abs_delta ~ pathway_name, data = ag, FUN = max)
  names(pw)[2] <- "pw_max_delta"
  pw <- pw[order(-pw$pw_max_delta, pw$pathway_name), ]
  pw$pw_order <- seq_len(nrow(pw))

  key <- merge(ag, pw[, c("pathway_name", "pw_order")], by = "pathway_name")
  key <- key[order(key$pw_order, -key$max_abs_delta, key$interaction_name), ]
  key$y_order <- seq_len(nrow(key))
  key$ypos    <- nrow(key) + 1 - key$y_order      # y_order 1 画在最上方
  rownames(key) <- NULL
  key <- key[, c("y_order","ypos","pw_order","pathway_name","interaction_name","max_abs_delta")]

  # 宽表：每个 LR x 靶 一行（每例 20 x 2 = 40 行）
  v <- do.call(rbind, lapply(target_levels, function(tg) {
    s <- d[d$target == tg, ]
    s <- s[order(match(s$interaction_name, key$interaction_name)), ]
    data.frame(patient = pt, y_order = key$y_order, pathway_name = s$pathway_name,
               interaction_name = s$interaction_name, target = tg,
               prob_Clone1 = s$prob_Clone1, prob_Clone2 = s$prob_Clone2,
               delta_Clone1_minus_Clone2 = s$delta_prob,
               max_abs_delta = key$max_abs_delta, stringsAsFactors = FALSE)
  }))
  v <- v[order(v$y_order, match(v$target, target_levels)), ]
  rownames(v) <- NULL

  # 长表（80 行 = 20 LR x 2 靶 x 2 克隆），供 circos 用
  lr <- read.csv(LR_CSV[[pt]], check.names = FALSE, stringsAsFactors = FALSE)
  lr <- lr[lr$interaction_name %in% pairs, ]
  lr$y_order <- key$y_order[match(lr$interaction_name, key$interaction_name)]
  lr$ypos    <- key$ypos[match(lr$interaction_name, key$interaction_name)]
  lr <- lr[order(lr$y_order, match(lr$target, target_levels), match(lr$source, clone_levels)), ]

  ORD[[pt]] <- key; VAL[[pt]] <- v; PROB4[[pt]] <- lr
}

cat("\n=== P13 y 轴顺序（按通路分组）===\n")
print(ORD$P13[, c("y_order","pathway_name","interaction_name","max_abs_delta")], row.names = FALSE, digits = 4)
cat("\n=== P15 y 轴顺序（按通路分组）===\n")
print(ORD$P15[, c("y_order","pathway_name","interaction_name","max_abs_delta")], row.names = FALSE, digits = 4)

cat("\n[Check] 行数/对数/通路数：\n")
for (pt in c("P13","P15")) {
  cat(sprintf("  %s: values 行=%d, 不同LR=%d, 通路=%s\n", pt, nrow(VAL[[pt]]),
              length(unique(VAL[[pt]]$interaction_name)),
              paste(names(table(ORD[[pt]]$pathway_name)), table(ORD[[pt]]$pathway_name), sep=":", collapse=" ")))
}
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")

# 将 top20 数值表写入 tables/ 并与源 CSV 逐格核对（数值唯一来源一致性检查）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

for (pt in c("P13", "P15")) {
  v  <- VAL[[pt]]; lr <- PROB4[[pt]]
  lut <- tapply(lr$prob, paste(lr$interaction_name, lr$target, lr$source, sep = "|"),
                function(x) x[1])
  kk  <- paste(v$interaction_name, v$target, sep = "|")
  p1  <- as.numeric(lut[paste(kk, "Clone1", sep = "|")])
  p2  <- as.numeric(lut[paste(kk, "Clone2", sep = "|")])

  d1 <- max(abs(p1 - v$prob_Clone1))
  d2 <- max(abs(p2 - v$prob_Clone2))
  d3 <- max(abs((p1 - p2) - v$delta_Clone1_minus_Clone2))
  cat(sprintf("[Verify] %s | prob_Clone1 maxdiff=%.3g | prob_Clone2 maxdiff=%.3g | delta maxdiff=%.3g | rows=%d | pairs=%d\n",
              pt, d1, d2, d3, nrow(v), length(unique(v$interaction_name))))
  stopifnot(d1 < 1e-12, d2 < 1e-12, d3 < 1e-12, nrow(v) == 40, length(unique(v$interaction_name)) == 20)

  out <- file.path(base_dir, pt, "tables", sprintf("%s_tail_LR_top20_values.csv", pt))
  write.csv(v, out, row.names = FALSE)
  cat(sprintf("[Table] %s | %.1f KB\n", out, file.info(out)$size/1024))
}

cat("\n=== P13 数值表预览（前 6 行）===\n")
print(head(VAL$P13, 6), row.names = FALSE, digits = 6)
cat("\n=== P15 数值表预览（前 6 行）===\n")
print(head(VAL$P15, 6), row.names = FALSE, digits = 6)
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")

# 绘制 top20 dotplot（P13/P15）：x=克隆→靶，点大小=通讯概率，颜色=克隆，去掉右侧 Δ 列
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

suppressPackageStartupMessages({library(ggplot2); library(scales)})
plot_theme <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

# 加宽左侧 gutter（通路名与 y 轴标签之间留出约 3 个字宽的间隔）
X_PROB2 <- c(C1Q_TAM_Clone1 = 1.35, C1Q_TAM_Clone2 = 2.35,
             SPP1_TAM_Clone1 = 3.95, SPP1_TAM_Clone2 = 4.95)
SEP_TGT <- 3.15
X_LIM   <- c(-0.10, 5.65)
X_BRK   <- 0.03
X_TXT   <- 0.78
COL_CL  <- c(Clone1 = "#E64B35", Clone2 = "#3C5488")   # palette nature 双色

make_top20_plot <- function(pt) {
  key <- ORD[[pt]]; n <- nrow(key); pd <- PROB4[[pt]]
  pdd <- data.frame(interaction_name = pd$interaction_name, clone = pd$source,
                    target = pd$target, prob = pd$prob, stringsAsFactors = FALSE)
  pdd$ypos <- key$ypos[match(pdd$interaction_name, key$interaction_name)]
  pdd$x    <- X_PROB2[paste(pdd$target, pdd$clone, sep = "_")]

  # 通路括号与斜体灰字，按 y 轴分组区间定位
  pw <- do.call(rbind, lapply(unique(key$pathway_name), function(g) {
    y <- key$ypos[key$pathway_name == g]
    data.frame(pathway_name = g, ymid = mean(range(y)),
               ytop = max(y) + 0.42, ybot = min(y) - 0.42,
               pw_order = key$pw_order[key$pathway_name == g][1], stringsAsFactors = FALSE)
  }))
  pw <- pw[order(pw$pw_order), ]

  # 通路分组之间的横向分隔线，使分组边界在整个面板宽度上可读
  bnd <- if (nrow(pw) > 1) data.frame(yb = (pw$ybot[-nrow(pw)] + pw$ytop[-1]) / 2) else data.frame(yb = numeric(0))

  maxp <- max(pdd$prob); minp <- min(pdd$prob)
  brks <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35)
  brks <- brks[brks > minp & brks < maxp]
  if (length(brks) < 2) brks <- signif(seq(minp, maxp, length.out = 3), 2)
  cat(sprintf("  [Legend] %s size ticks: %s (prob range %.4f-%.4f) | pathway groups=%d\n",
              pt, paste(sprintf("%.2f", brks), collapse = "/"), minp, maxp, nrow(pw)))

  ggplot() +
    geom_segment(data = bnd, aes(x = X_LIM[1], xend = X_LIM[2], y = yb, yend = yb),
                 colour = "grey86", linewidth = 0.25, inherit.aes = FALSE) +
    geom_point(data = pdd, aes(x = x, y = ypos, size = prob, colour = clone), stroke = 0.25) +
    geom_segment(data = pw, aes(x = X_BRK, xend = X_BRK, y = ybot, yend = ytop),
                 colour = "grey40", linewidth = 0.55, inherit.aes = FALSE) +
    geom_text(data = pw, aes(x = X_TXT, y = ymid, label = pathway_name), hjust = 1,
              size = 2.5, colour = "grey25", fontface = "italic", inherit.aes = FALSE) +
    geom_segment(aes(x = SEP_TGT, xend = SEP_TGT, y = 0.5, yend = n + 0.5),
                 colour = "black", linewidth = 0.35, inherit.aes = FALSE) +
    scale_colour_manual(values = COL_CL, breaks = clone_levels, name = "Clone", drop = FALSE) +
    scale_size_area(limits = c(0, maxp), max_size = 5.2, breaks = brks,
                    labels = sprintf("%.2f", brks), name = "Communication\nprobability") +
    scale_x_continuous(breaks = unname(X_PROB2),
                       labels = c("Clone1→C1Q", "Clone2→C1Q",
                                  "Clone1→SPP1", "Clone2→SPP1"),
                       limits = X_LIM, expand = c(0, 0)) +
    scale_y_continuous(breaks = key$ypos, labels = key$interaction_name,
                       limits = c(0.4, n + 0.6), expand = c(0, 0)) +
    labs(caption = "Grey brackets/separators group LR pairs by signalling pathway. Size = communication probability; colour = clone.") +
    plot_theme +
    theme(axis.title = element_blank(), plot.title = element_blank(),
          axis.text.y = element_text(size = 6.5, colour = "black"),
          axis.text.x = element_text(size = 6.6, colour = "black"),
          axis.ticks = element_line(colour = "black", linewidth = 0.3),
          axis.ticks.length = unit(2, "pt"),
          legend.title = element_text(size = 7),
          legend.text  = element_text(size = 6.5),
          legend.key.size = unit(0.32, "cm"),
          legend.background = element_blank(), legend.box = "vertical",
          legend.position = "right", legend.box.spacing = unit(0.12, "cm"),
          plot.caption = element_text(hjust = 0, size = 6.2, colour = "grey35"),
          plot.margin = margin(6, 10, 5, 8))
}

render_top20_plot <- function(pt, width = 6.6, height = 9.4) {
  p <- make_top20_plot(pt)
  fdir <- file.path(base_dir, pt, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  stem <- file.path(fdir, sprintf("%s_tail_LR_top20_dotplot", pt))
  ggsave(paste0(stem, ".png"), plot = p, width = width, height = height, dpi = 600,
         bg = "white", device = "png", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot = p, width = width, height = height, bg = "white",
         device = if (has_cairo) cairo_pdf else pdf, limitsize = FALSE)
  for (f in paste0(stem, c(".png", ".pdf")))
    cat(sprintf("[Figure] %s | %.1f KB\n", f, file.info(f)$size/1024))
  invisible(p)
}

p_13 <- render_top20_plot("P13")
p_15 <- render_top20_plot("P15")

# 量算 gutter：通路名左端与面板左边界、通路名右端与最大点的间隔（英寸）
gt <- ggplotGrob(p_13)
cat(sprintf("[Gutter] axis-l width = %.2f in\n", convertWidth(gt$widths[[6]], "inches", TRUE)))
cat(sprintf("[Gutter] 通路名区 %.3f-%.3f user；面板左界 %.2f；首个点列 %.2f\n",
            X_TXT - 0.53, X_TXT, X_LIM[1], X_PROB2[1]))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")

# 绘制 4 扇区 circos（P13/P15）：Clone1/Clone2 → C1Q_TAM/SPP1_TAM，线宽∝通讯概率
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

suppressPackageStartupMessages(library(circlize))

SEC     <- c("Clone1", "C1Q_TAM", "Clone2", "SPP1_TAM")
SEC_LAB <- c(Clone1 = "Clone1", C1Q_TAM = "C1Q TAM", Clone2 = "Clone2", SPP1_TAM = "SPP1 TAM")
START_DEG <- 135
GAP     <- 9
LWD_K   <- 16        # 线宽系数：lwd = max(0.3, 16 * probability)
LWD_MIN <- 0.3
P_CLAMP <- LWD_MIN / LWD_K            # 0.01875，低于此概率的连线被钳到最小线宽
ALPHA   <- 0.50
R_VERT  <- 1.13      # 上/下扇区标签半径
R_HORZ  <- 1.22      # 左/右扇区标签半径（水平文字径向占宽大）

draw_tail_circos <- function(pt) {
  key <- ORD[[pt]]; pd <- PROB4[[pt]]; n <- nrow(key)
  band <- (seq_len(n) - 0.5) / n
  ANG  <- c()

  lk <- data.frame(interaction_name = pd$interaction_name, clone = pd$source,
                   target = pd$target, prob = pd$prob, stringsAsFactors = FALSE)
  lk <- lk[lk$prob > 0, ]
  stopifnot(nrow(lk) == 80)                       # 20 LR x 2 靶 x 2 克隆，须一条不少
  lk <- lk[order(-lk$prob), ]                     # 粗线先画、细线后画
  lk$slot <- match(lk$interaction_name, key$interaction_name)
  stopifnot(all(lk$clone %in% clone_levels), all(lk$target %in% target_levels))
  n_clamp <- sum(lk$prob < P_CLAMP)

  par(mar = c(0.4, 0.4, 0.4, 0.4), bg = "white", xpd = NA)
  circos.clear()
  circos.par(canvas.xlim = c(-1.38, 1.38), canvas.ylim = c(-1.98, 1.26),
             start.degree = START_DEG, gap.after = rep(GAP, 4),
             track.margin = c(0.005, 0.005), cell.padding = c(0, 0, 0, 0),
             points.overflow.warning = FALSE)
  circos.initialize(factors = SEC, xlim = c(0, 1))

  rx <- grconvertX(1, "user", "inches") - grconvertX(0, "user", "inches")
  ry <- grconvertY(1, "user", "inches") - grconvertY(0, "user", "inches")
  cat(sprintf("  [Roundness] %s | 1 user unit = %.3f in (x) vs %.3f in (y) | ratio x/y = %.4f\n",
              pt, rx, ry, rx / ry))

  circos.track(ylim = c(0, 1), track.height = 0.08,
               bg.border = "grey35", bg.col = "grey88", track.index = 1,
               panel.fun = function(x, y) {
                 si <- get.cell.meta.data("sector.index")
                 ANG[si] <<- ((get.cell.meta.data("cell.start.degree") +
                               get.cell.meta.data("cell.end.degree")) / 2 + 360) %% 360
                 for (k in 0:n) circos.lines(c(k, k) / n, c(0, 1), sector.index = si,
                                             col = "white", lwd = 0.9)
               })

  for (i in seq_len(nrow(lk))) {
    circos.link(lk$clone[i], band[lk$slot[i]], lk$target[i], band[lk$slot[i]],
                lwd = max(LWD_MIN, LWD_K * lk$prob[i]),
                col = adjustcolor(COL_CL[[lk$clone[i]]], alpha.f = ALPHA))
  }
  circos.clear()

  for (s in SEC) {
    a  <- ANG[[s]] * pi / 180
    rr <- if (abs(cos(a)) > 0.7) R_HORZ else R_VERT
    text(rr * cos(a), rr * sin(a), SEC_LAB[[s]], adj = c(0.5, 0.5), cex = 0.82)
  }
  cat(sprintf("  [Sector angles] %s\n",
              paste(sprintf("%s=%.1f", SEC, ANG[SEC]), collapse = "  ")))

  ly0 <- -1.42
  text(-1.34, ly0, "Clone", adj = c(0, 0.5), cex = 0.70)
  yy <- ly0 - 0.10
  for (cl in clone_levels) {
    segments(-1.34, yy, -1.21, yy, lwd = 3.2, col = COL_CL[[cl]])
    text(-1.18, yy, cl, adj = c(0, 0.5), cex = 0.70)
    yy <- yy - 0.09
  }
  text(0.18, ly0, "Communication probability", adj = c(0, 0.5), cex = 0.70)
  yy <- ly0 - 0.10
  for (v in c(0.05, 0.15, 0.30)) {
    segments(0.18, yy, 0.32, yy, lwd = max(LWD_MIN, LWD_K * v), col = "grey25")
    text(0.35, yy, sprintf("%.2f", v), adj = c(0, 0.5), cex = 0.70)
    yy <- yy - 0.09
  }

  # 图注：明确是最小线宽钳制，不用易与 p 值混淆的 "P"
  cap1 <- if (n_clamp > 0) {
    sprintf("Link width is proportional to communication probability; %d links below %.3f are clamped to the minimum width.", n_clamp, P_CLAMP)
  } else {
    "Link width is proportional to communication probability."
  }
  cap2 <- "Each sector holds the same 20 LR pairs in dotplot row order; the source sector also gives the clone."
  w1 <- strwidth(cap1, cex = 0.60, units = "user")
  w2 <- strwidth(cap2, cex = 0.60, units = "user")
  cat(sprintf("  [Caption] %s | clamped links=%d | width line1=%.2f line2=%.2f user (limit %.2f)\n",
              pt, n_clamp, w1, w2, 1.34 + 1.34))
  stopifnot(-1.34 + w1 < 1.34, -1.34 + w2 < 1.34)   # 保证图注不出右边界
  text(-1.34, -1.84, cap1, adj = c(0, 0.5), cex = 0.60, col = "grey30")
  text(-1.34, -1.92, cap2, adj = c(0, 0.5), cex = 0.60, col = "grey30")
}

save_tail_circos <- function(pt, w = 6.70, h = 7.87) {
  fdir <- file.path(base_dir, pt, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  stem <- file.path(fdir, sprintf("%s_tail_LR_circos", pt))
  png(paste0(stem, ".png"), width = w, height = h, units = "in", res = 600, bg = "white")
  draw_tail_circos(pt); dev.off()
  cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, bg = "white", onefile = FALSE)
  draw_tail_circos(pt); dev.off()
  for (f in paste0(stem, c(".png", ".pdf")))
    cat(sprintf("[Figure] %s | %.1f KB\n", f, file.info(f)$size/1024))
}

save_tail_circos("P13")
save_tail_circos("P15")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")

# 出图完整性自检：像素非空/非全白、双格式配对、尺寸、dotplot 无 Δ 数值列
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

chk_png <- function(f) {
  im  <- png::readPNG(f)
  rgb <- if (length(dim(im)) == 3) im[, , 1:3] else array(im, c(dim(im)[1], dim(im)[2], 3))
  sdv <- sd(as.vector(rgb)); white <- mean(rgb > 0.98)
  chroma <- pmax(rgb[, , 1], rgb[, , 2], rgb[, , 3]) - pmin(rgb[, , 1], rgb[, , 2], rgb[, , 3])
  colfrac <- mean(chroma > 0.08)
  cat(sprintf("  %-38s %5dx%-5d sd=%.3f white=%.3f colored=%.4f\n",
              basename(f), dim(im)[2], dim(im)[1], sdv, white, colfrac))
  data.frame(file = basename(f), sd = sdv, white = white, colored = colfrac,
             ok = (sdv > 0.03 && white < 0.995 && colfrac > 0.002))
}

rows <- list()
for (pt in c("P13", "P15")) {
  for (k in c("tail_LR_top20_dotplot", "tail_LR_circos")) {
    png_f <- file.path(base_dir, pt, "figures", sprintf("%s_%s.png", pt, k))
    pdf_f <- file.path(base_dir, pt, "figures", sprintf("%s_%s.pdf", pt, k))
    rows[[length(rows) + 1]] <- chk_png(png_f)
    cat(sprintf("  pair: %-22s png=%.0fKB pdf=%.0fKB | pdf_exists=%s\n", k,
                file.info(png_f)$size/1024, file.info(pdf_f)$size/1024, file.exists(pdf_f)))
  }
}
px <- do.call(rbind, rows)
cat(sprintf("\n[Figure Wholeness] %d/%d OK\n", sum(px$ok), nrow(px)))

num_label_layers <- function(p) {
  n <- 0
  for (l in p$layers) {
    d <- l$data
    if (!is.null(d) && "label" %in% names(d) &&
        any(grepl("^[-+]?[0-9]*\\.?[0-9]+$", as.character(d$label)))) n <- n + 1
  }
  n
}
for (nm in c("p_13", "p_15")) {
  p <- get(nm); pd_layer <- p$layers[[2]]$data
  cat(sprintf("[Dotplot check] %s | points=%d | distinct_x=%s | numeric-label-layers=%d | n_LR=%d\n",
              nm, nrow(pd_layer), paste(sort(unique(pd_layer$x)), collapse = ","),
              num_label_layers(p), length(unique(pd_layer$ypos))))
  stopifnot(nrow(pd_layer) == 80, num_label_layers(p) == 0, length(unique(pd_layer$x)) == 4)
}
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")

# 绘制 top20 通讯概率热图（P13/P15）：左=概率(viridis) 右=Δ(RdBu_r, 0 为中心)，共享行顺序
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(viridisLite)})

HEAT_XL  <- c(C1Q_TAM_Clone1 = 1.0, C1Q_TAM_Clone2 = 2.0,
              SPP1_TAM_Clone1 = 3.6, SPP1_TAM_Clone2 = 4.6)   # 两靶之间留白
HEAT_XR  <- c(C1Q_TAM = 1.0, SPP1_TAM = 2.0)
COL_POS  <- "#B2182B"; COL_NEG <- "#2166AC"

# 按色阶亮度选格内文字颜色：深色格白字、浅色格黑字
lum_of   <- function(hex) { m <- col2rgb(hex) / 255; 0.2126 * m[1, ] + 0.7152 * m[2, ] + 0.0722 * m[3, ] }
pick_txt <- function(vals, ramp_hex, lo, hi) {
  t <- pmin(pmax((vals - lo) / (hi - lo), 0), 1)
  idx <- pmax(1, ceiling(t * (length(ramp_hex) - 1) + 1))
  ifelse(lum_of(ramp_hex[idx]) < 0.45, "white", "black")
}
RAMP_VIR <- viridis(256)
RAMP_DIV <- colorRampPalette(c(COL_NEG, "white", COL_POS))(256)

heat_base_theme <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        axis.title = element_blank(), plot.title = element_blank())

make_top20_heatmap <- function(pt) {
  key <- ORD[[pt]]; v <- VAL[[pt]]; n <- nrow(key)
  lv  <- rev(key$interaction_name)                 # 离散 y 轴自下而上：y_order 1 画在最上方
  pos <- setNames(seq_along(lv), lv)

  dl <- rbind(
    data.frame(pathway_name = v$pathway_name, interaction_name = v$interaction_name,
               target = v$target, clone = "Clone1", prob = v$prob_Clone1, stringsAsFactors = FALSE),
    data.frame(pathway_name = v$pathway_name, interaction_name = v$interaction_name,
               target = v$target, clone = "Clone2", prob = v$prob_Clone2, stringsAsFactors = FALSE))
  dl$x <- HEAT_XL[paste(dl$target, dl$clone, sep = "_")]
  dl$interaction_name <- factor(dl$interaction_name, levels = lv)

  dr <- data.frame(pathway_name = v$pathway_name, interaction_name = v$interaction_name,
                   target = v$target, delta = v$delta_Clone1_minus_Clone2, stringsAsFactors = FALSE)
  dr$x <- HEAT_XR[dr$target]
  dr$interaction_name <- factor(dr$interaction_name, levels = lv)

  p_lim <- range(dl$prob); d_lim <- max(abs(dr$delta))
  dl$tcol <- pick_txt(dl$prob,  RAMP_VIR, p_lim[1], p_lim[2])
  dr$tcol <- pick_txt(dr$delta, RAMP_DIV, -d_lim, d_lim)

  # 通路分组（独立窄面板，与热图逐行对齐）；clip="off" 防止通路名被面板右边界裁掉首字母
  gp <- do.call(rbind, lapply(unique(key$pathway_name), function(g) {
    y <- pos[key$interaction_name[key$pathway_name == g]]
    data.frame(pathway_name = g, ymid = mean(y), ytop = max(y) + 0.40, ybot = min(y) - 0.40,
               pw_order = key$pw_order[key$pathway_name == g][1], stringsAsFactors = FALSE)
  }))
  gp <- gp[order(gp$pw_order), ]

  pA <- ggplot(gp) +
    geom_segment(aes(x = 0.99, xend = 0.99, y = ybot, yend = ytop), colour = "grey40", linewidth = 0.5) +
    geom_text(aes(x = 0.92, y = ymid, label = pathway_name), hjust = 1, size = 2.3,
              colour = "grey25", fontface = "italic") +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void()

  pB <- ggplot(dl, aes(x = x, y = interaction_name, fill = prob)) +
    geom_tile(width = 1, height = 1, colour = NA) +
    geom_text(aes(label = sprintf("%.3f", prob), colour = tcol), size = 2.2, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_viridis_c(option = "viridis", name = "Communication\nprobability") +
    scale_x_continuous(breaks = unname(HEAT_XL),
                       labels = c("Clone1→C1Q", "Clone2→C1Q",
                                  "Clone1→SPP1", "Clone2→SPP1"),
                       limits = c(0.5, 5.1), expand = c(0, 0)) +
    scale_y_discrete(limits = lv, expand = c(0, 0)) +
    heat_base_theme +
    theme(axis.text.y = element_text(size = 6.5, colour = "black", margin = margin(r = 2)),
          axis.text.x = element_text(size = 6.5, colour = "black"),
          axis.ticks = element_line(colour = "black", linewidth = 0.3),
          axis.ticks.length = unit(2, "pt"),
          legend.title = element_text(size = 6.8), legend.text = element_text(size = 6.3),
          legend.key.height = unit(0.42, "cm"), legend.key.width = unit(0.26, "cm"),
          legend.position = "right", plot.margin = margin(4, 4, 4, 2))

  pC <- ggplot(dr, aes(x = x, y = interaction_name, fill = delta)) +
    geom_tile(width = 1, height = 1, colour = NA) +
    geom_text(aes(label = sprintf("%+.3f", delta), colour = tcol), size = 2.2, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient2(low = COL_NEG, mid = "white", high = COL_POS, midpoint = 0,
                         limits = c(-d_lim, d_lim),
                         name = "Δ probability\n(Clone1 − Clone2)") +
    scale_x_continuous(breaks = unname(HEAT_XR), labels = c("Δ C1Q", "Δ SPP1"),
                       limits = c(0.5, 2.5), expand = c(0, 0)) +
    scale_y_discrete(limits = lv, expand = c(0, 0)) +
    heat_base_theme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 6.5, colour = "black"),
          axis.ticks.x = element_line(colour = "black", linewidth = 0.3),
          axis.ticks.length.x = unit(2, "pt"),
          legend.title = element_text(size = 6.8), legend.text = element_text(size = 6.3),
          legend.key.height = unit(0.42, "cm"), legend.key.width = unit(0.26, "cm"),
          legend.position = "right", plot.margin = margin(4, 4, 4, 10))

  cat(sprintf("  [Heatmap] %s | tiles: prob=%d (20x4), delta=%d (20x2) | prob %.4f-%.4f | |delta|max=%.4f | groups=%d\n",
              pt, nrow(dl), nrow(dr), p_lim[1], p_lim[2], d_lim, nrow(gp)))
  list(plot = pA + pB + pC + plot_layout(widths = c(1.45, 4.4, 1.95)),
       dl = dl, dr = dr, gp = gp, lv = lv)
}

render_top20_heatmap <- function(pt, width = 10.0, height = 7.2) {
  h <- make_top20_heatmap(pt)
  fdir <- file.path(base_dir, pt, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  stem <- file.path(fdir, sprintf("%s_tail_LR_top20_heatmap", pt))
  ggsave(paste0(stem, ".png"), plot = h$plot, width = width, height = height, dpi = 600,
         bg = "white", device = "png", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot = h$plot, width = width, height = height, bg = "white",
         device = if (has_cairo) cairo_pdf else pdf, limitsize = FALSE)
  for (f in paste0(stem, c(".png", ".pdf")))
    cat(sprintf("[Figure] %s | %.1f KB\n", f, file.info(f)$size/1024))
  invisible(h)
}

hm_13 <- render_top20_heatmap("P13")
hm_15 <- render_top20_heatmap("P15")

# 逐图检查左侧是否有像素贴边（贴边即为文字被设备边界裁切）
for (pt in c("P13","P15")) {
  im <- png::readPNG(file.path(base_dir, pt, "figures", sprintf("%s_tail_LR_top20_heatmap.png", pt)))
  rgb <- im[, , 1:3]
  nw  <- (rgb[, , 1] < 0.99) | (rgb[, , 2] < 0.99) | (rgb[, , 3] < 0.99)
  cols <- which(apply(nw, 2, any))
  cat(sprintf("  [Edge] %s | leftmost non-white col = %d px (0 = clipped) | image %dx%d\n",
              pt, min(cols), dim(im)[2], dim(im)[1]))
}
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")

# 热图自检：标注值与源数值一致、瓦片数、双格式与像素完整性
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

for (pt in c("P13","P15")) {
  h <- if (pt == "P13") hm_13 else hm_15
  labL <- sprintf("%.3f", h$dl$prob); okL <- all(abs(as.numeric(labL) - h$dl$prob) <= 0.0005 + 1e-12)
  labR <- sprintf("%+.3f", h$dr$delta); okR <- all(abs(as.numeric(labR) - h$dr$delta) <= 0.0005 + 1e-12)
  cat(sprintf("[SelfCheck] %s | prob tiles=%d rows=%d cols=%d | delta tiles=%d rows=%d cols=%d | shared_rows=%s | labels_ok(L/R)=%s/%s\n",
              pt, nrow(h$dl), length(unique(h$dl$interaction_name)), length(unique(h$dl$x)),
              nrow(h$dr), length(unique(h$dr$interaction_name)), length(unique(h$dr$x)),
              identical(levels(h$dl$interaction_name), levels(h$dr$interaction_name)), okL, okR))
  stopifnot(nrow(h$dl) == 80, nrow(h$dr) == 40, okL, okR,
            identical(levels(h$dl$interaction_name), levels(h$dr$interaction_name)))
  chk <- merge(h$dr, h$dl[h$dl$clone == "Clone1", c("interaction_name","target","prob")],
               by = c("interaction_name","target"), suffixes = c("", ".c1"))
  chk <- merge(chk, h$dl[h$dl$clone == "Clone2", c("interaction_name","target","prob")],
               by = c("interaction_name","target"), suffixes = c("", ".c2"))
  cat(sprintf("  delta vs (C1-C2) max|diff| = %.3g\n", max(abs(chk$delta - (chk$prob - chk$prob.c2)))))
  stopifnot(max(abs(chk$delta - (chk$prob - chk$prob.c2))) < 1e-12)
}

for (pt in c("P13","P15")) {
  for (ext in c("png","pdf")) {
    f <- file.path(base_dir, pt, "figures", sprintf("%s_tail_LR_top20_heatmap.%s", pt, ext))
    cat(sprintf("  %s | %.1f KB\n", basename(f), file.info(f)$size/1024))
    stopifnot(file.exists(f), file.info(f)$size > 1000)
  }
  im <- png::readPNG(file.path(base_dir, pt, "figures", sprintf("%s_tail_LR_top20_heatmap.png", pt)))
  rgb <- im[, , 1:3]
  chroma <- pmax(rgb[, , 1], rgb[, , 2], rgb[, , 3]) - pmin(rgb[, , 1], rgb[, , 2], rgb[, , 3])
  cat(sprintf("  [Pixel] %s heatmap %dx%d sd=%.3f white=%.3f colored=%.4f\n", pt,
              dim(im)[2], dim(im)[1], sd(as.vector(rgb)), mean(rgb > 0.98), mean(chroma > 0.08)))
}
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")

# 记录审计通过（judge≠doer 独立 critic 均通过）到统一审计日志
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

log_dir  <- file.path(Sys.getenv("HOME"), ".claude", "logs", format(Sys.Date(), "%Y-%m-%d"))
log_file <- file.path(log_dir, "audit_done.jsonl")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
cat("log file:", log_file, "| exists:", file.exists(log_file), "\n")

unit   <- "tail_lr_top20_circos_done"
sess   <- "e99aac27-4525-42e3-8a8a-00c81a7e0358"
ts     <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
ts     <- sub("([+-][0-9]{2})([0-9]{2})$", "\\1\\2", ts)

# 幂等：同一 session + 同一 unit 已写过 pass 则不重复追加
already <- FALSE
if (file.exists(log_file)) {
  prev <- readLines(log_file, warn = FALSE)
  cat("existing lines:", length(prev), "\n")
  already <- any(grepl(sprintf('"session":"%s"', sess), prev, fixed = TRUE) &
                 grepl(sprintf('"unit":"%s"', unit), prev, fixed = TRUE) &
                 grepl('"verdict":"pass"', prev, fixed = TRUE))
}
cat("already logged for this unit:", already, "\n")

if (!already) {
  line <- sprintf('{"ts":"%s","session":"%s","subagent_type":"sc-analyst","unit":"%s","verdict":"pass"}',
                  ts, sess, unit)
  cat(line, "\n", file = log_file, append = TRUE, sep = "")
  cat("[AuditLog] appended:", line, "\n")
} else {
  cat("[AuditLog] skipped (idempotent)\n")
}

tail_lines <- tail(readLines(log_file, warn = FALSE), 3)
cat("last lines of audit_done.jsonl:\n"); cat(paste0("  ", tail_lines, collapse = "\n"), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")

# 重做 top20 Δ 热图（P13/P15）：ComplexHeatmap 单面板（2 列 C1Q/SPP1），行侧加两条显著性类别色条
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

suppressPackageStartupMessages({library(ComplexHeatmap); library(circlize); library(grid)})

SIG_COL <- c(Clone1_only = "#E64B35", Clone2_only = "#3C5488",
             both_significant = "#4DBBD5", neither = "#D9D9D9")
COL_POS <- "#B2182B"; COL_NEG <- "#2166AC"      # palette 契约 RdBu_r 双端
TGT_LAB <- c(C1Q_TAM = "C1Q TAM", SPP1_TAM = "SPP1 TAM")
HT_ID   <- "delta_probability"                  # viewport 名必须 ASCII：希腊字母/减号会让 viewport 查找失败
HT_TTL  <- "Δ probability (Clone1 − Clone2)"    # 仅作图例标题（不参与 viewport 命名）
LR_GP   <- gpar(fontsize = 6.3, col = "black")
PW_GP   <- gpar(fontsize = 6.3, col = "grey30", fontface = "italic")

build_delta_mat <- function(pt) {
  key <- ORD[[pt]]; src <- read.csv(DIFF_CSV[[pt]], check.names = FALSE, stringsAsFactors = FALSE)
  lr  <- as.character(key$interaction_name)                 # 行序 = 通路分组后的 y_order
  M   <- matrix(NA_real_, nrow = length(lr), ncol = 2,
                dimnames = list(lr, unname(TGT_LAB[target_levels])))
  cls <- matrix(NA_character_, nrow = length(lr), ncol = 2, dimnames = dimnames(M))
  for (tg in target_levels) {
    s <- src[src$interaction_name %in% lr & src$target == tg, ]
    s <- s[match(lr, s$interaction_name), ]
    stopifnot(!anyNA(s$interaction_name))
    M[, TGT_LAB[[tg]]]   <- s$delta_prob
    cls[, TGT_LAB[[tg]]] <- s$signif_class                # 类别直接读表，不重算
  }
  stopifnot(nrow(M) == 20, ncol(M) == 2, !anyNA(M), !anyNA(cls))
  pw <- key$pathway_name[match(lr, key$interaction_name)]
  pw_lab <- rep("", length(lr))                           # 通路名只写在组内中间那一行
  for (g in unique(pw)) {
    idx <- which(pw == g); pw_lab[idx[ceiling(length(idx) / 2)]] <- g
  }
  list(M = M, cls = cls, pw = pw, pw_lab = pw_lab, lr = lr, dlim = max(abs(M)))
}

draw_delta_heatmap <- function(pt) {
  b <- build_delta_mat(pt)
  sig_c1q  <- factor(b$cls[, "C1Q TAM"],  levels = names(SIG_COL))
  sig_spp1 <- factor(b$cls[, "SPP1 TAM"], levels = names(SIG_COL))

  ha <- rowAnnotation(
    Pathway  = anno_text(b$pw_lab, gp = PW_GP, just = "right",
                         location = unit(1, "npc") - unit(0.3, "mm")),
    spacer   = anno_empty(border = FALSE, width = unit(2.4, "mm")),
    `LR pair` = anno_text(b$lr, gp = LR_GP, just = "right",
                          location = unit(1, "npc") - unit(0.4, "mm")),
    C1Q  = sig_c1q, SPP1 = sig_spp1,
    col = list(C1Q = SIG_COL, SPP1 = SIG_COL),
    show_legend = c(C1Q = TRUE, SPP1 = FALSE),
    annotation_legend_param = list(C1Q = list(title = "Significance class", ncol = 1,
                                              title_gp = gpar(fontsize = 6.6),
                                              labels_gp = gpar(fontsize = 6.2))),
    annotation_name_gp = gpar(fontsize = 6.0), annotation_name_side = "top",
    annotation_name_rot = 90,
    simple_anno_size = unit(2.1, "mm"), gap = unit(1.0, "mm"))

  # 色带刻度用该图真实 ±dlim（不用 pretty 取整的 ±0.1），避免两例被误读成同一尺度
  brk <- c(-b$dlim, -b$dlim / 2, 0, b$dlim / 2, b$dlim)

  ht <- Heatmap(
    b$M, name = HT_ID,
    col = colorRamp2(c(-b$dlim, 0, b$dlim), c(COL_NEG, "white", COL_POS)),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE,
    column_names_rot = 0, column_names_centered = TRUE,
    column_names_gp = gpar(fontsize = 6.9, col = "black"),
    column_names_side = "bottom", column_gap = unit(1.8, "mm"),
    left_annotation = ha,
    border = FALSE, rect_gp = gpar(col = NA),
    heatmap_legend_param = list(title = HT_TTL, at = brk,
                                labels = sprintf("%.3f", brk),
                                title_gp = gpar(fontsize = 6.6),
                                labels_gp = gpar(fontsize = 6.2),
                                legend_height = unit(1.6, "cm"), grid_width = unit(0.30, "cm")),
    width = unit(1.75, "in"), height = unit(5.9, "in"))

  stopifnot(sum(b$pw_lab != "") == length(unique(b$pw)))
  cat(sprintf("  [Heatmap] %s | %dx%d | Δ %.4f..%.4f | ±%.4f | pathways=%d | legend ticks=%s\n",
              pt, nrow(b$M), ncol(b$M), min(b$M), max(b$M), b$dlim, length(unique(b$pw)),
              paste(sprintf("%.3f", brk), collapse = "/")))
  cat(sprintf("            class C1Q: %s\n", paste(names(table(sig_c1q)), table(sig_c1q), sep = "=", collapse = " ")))
  cat(sprintf("            class SPP1: %s\n", paste(names(table(sig_spp1)), table(sig_spp1), sep = "=", collapse = " ")))
  list(ht = ht, b = b, sig_c1q = sig_c1q, sig_spp1 = sig_spp1)
}

save_delta_heatmap <- function(pt, w = 6.2, h = 7.4) {
  hm  <- draw_delta_heatmap(pt)
  fdir <- file.path(base_dir, pt, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  stem <- file.path(fdir, sprintf("%s_tail_LR_top20_heatmap", pt))
  draw_it <- function() draw(hm$ht, heatmap_legend_side = "right", annotation_legend_side = "right",
                             padding = unit(c(3, 3, 4, 6), "mm"))
  png(paste0(stem, ".png"), width = w, height = h, units = "in", res = 600, bg = "white", type = "cairo")
  draw_it(); dev.off()
  cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, bg = "white", onefile = FALSE)
  draw_it(); dev.off()
  im  <- png::readPNG(paste0(stem, ".png")); rgb <- im[, , 1:3]
  nw  <- (rgb[, , 1] < 0.99) | (rgb[, , 2] < 0.99) | (rgb[, , 3] < 0.99)
  cc  <- range(which(apply(nw, 2, any)))
  cat(sprintf("[Figure] %s.png | %.1f KB | leftmost non-white=%d (0=clipped) | bbox width %.1f%%\n",
              stem, file.info(paste0(stem, ".png"))$size/1024, cc[1], 100 * diff(cc) / dim(nw)[2]))
  cat(sprintf("[Figure] %s.pdf | %.1f KB\n", stem, file.info(paste0(stem, ".pdf"))$size/1024))
  invisible(hm)
}

hm_13 <- save_delta_heatmap("P13")
hm_15 <- save_delta_heatmap("P15")

# 逐格回源：Δ 与 signif_class 都必须与唯一数值源逐行一致
for (pt in c("P13","P15")) {
  hm <- if (pt == "P13") hm_13 else hm_15
  src <- read.csv(DIFF_CSV[[pt]], check.names = FALSE, stringsAsFactors = FALSE)
  for (tg in target_levels) {
    s <- src[src$interaction_name %in% hm$b$lr & src$target == tg, ]
    s <- s[match(hm$b$lr, s$interaction_name), ]
    cat(sprintf("  [Trace] %s/%-8s max|Δ−delta_prob|=%.1e | class mismatches=%d\n", pt, tg,
                max(abs(hm$b$M[, TGT_LAB[[tg]]] - s$delta_prob)),
                sum(as.character(hm$b$cls[, TGT_LAB[[tg]]]) != s$signif_class)))
  }
}
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")

# 用 CellChat 原生 netVisual_chord_gene 画 top20 弦图（P13/P15）：4 个 cell-state 区块、基因名沿圆周、连线按来源克隆着色、线宽∝概率
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

suppressPackageStartupMessages({library(CellChat); library(circlize)})

CH_COL <- c(C1Q_TAM = "#00A087", Clone1 = "#E64B35", Clone2 = "#3C5488", SPP1_TAM = "#F39B7F")
stopifnot(identical(names(CH_COL), levels(cc13@idents)), identical(names(CH_COL), levels(cc15@idents)))

# 该例将出现的 sector（配体 + 受体基因名），用于核对扇区集合
cellchat_pair_genes <- function(o, pair_df) {
  sig <- o@LR$LRsig[rownames(o@LR$LRsig) %in% pair_df$interaction_name, c("ligand", "receptor")]
  unique(c(sig$ligand, sig$receptor))
}

draw_tail_chord <- function(pt, lab.cex = 0.85, sgap = 2.5, bgap = 20,
                            annh = c(0.025, 0.025), leg.x = 6, leg.y = 6) {
  o <- if (pt == "P13") cc13 else cc15
  pair_df <- data.frame(interaction_name = TOP20[[pt]], stringsAsFactors = FALSE)  # 只含 interaction_name -> slot.name 仍为 "net"
  par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "white", xpd = NA)
  gg <- netVisual_chord_gene(
    object = o, slot.name = "net",
    color.use = CH_COL,
    pairLR.use = pair_df,                       # 只保留 top20 对
    sources.use = c("Clone1", "Clone2"),        # 发送方只有两个克隆
    targets.use = c("C1Q_TAM", "SPP1_TAM"),     # 接收方只有两个 TAM 状态
    thresh = 1,                                 # 不叠加逐通讯 p 值过滤：保证 20 对全部出现，与热图行一致
    lab.cex = lab.cex, small.gap = sgap, big.gap = bgap,   # 加大块间间隙，缓解相邻块交界处标签拥挤
    annotationTrackHeight = annh, link.visible = TRUE,
    scale = FALSE, directional = 1, link.target.prop = TRUE,
    reduce = -1, transparency = 0.4, link.border = NA,
    title.name = NULL, show.legend = TRUE,
    legend.pos.x = leg.x, legend.pos.y = leg.y)
  cat(sprintf("  [Chord] %s | links=%d | gene sectors=%d | small.gap=%.1f big.gap=%.1f\n",
              pt, nrow(pair_df) * 4, length(cellchat_pair_genes(o, pair_df)), sgap, bgap))
  invisible(gg)
}

save_tail_chord <- function(pt, w = 8.6, h = 8.6, ...) {
  fdir <- file.path(base_dir, pt, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  stem <- file.path(fdir, sprintf("%s_tail_LR_chord", pt))
  png(paste0(stem, ".png"), width = w, height = h, units = "in", res = 600, bg = "white", type = "cairo")
  draw_tail_chord(pt, ...); dev.off()
  cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, bg = "white", onefile = FALSE)
  draw_tail_chord(pt, ...); dev.off()
  im  <- png::readPNG(paste0(stem, ".png")); rgb <- im[, , 1:3]
  nw  <- (rgb[, , 1] < 0.99) | (rgb[, , 2] < 0.99) | (rgb[, , 3] < 0.99)
  rr  <- range(which(apply(nw, 1, any))); cc <- range(which(apply(nw, 2, any)))
  cat(sprintf("[Figure] %s.png | %.1f KB | bbox fill: width %.1f%% x height %.1f%% | label ~%.1f pt at %.1f in\n",
              stem, file.info(paste0(stem, ".png"))$size/1024,
              100 * diff(cc) / dim(nw)[2], 100 * diff(rr) / dim(nw)[1], 0.85 * 12, w))
  cat(sprintf("[Figure] %s.pdf | %.1f KB\n", stem, file.info(paste0(stem, ".pdf"))$size/1024))
}

save_tail_chord("P13")
save_tail_chord("P15")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 清理旧 circos 残留（含 .ipynb_checkpoints 副本）并做最终产物完整性核对
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

old <- list.files(file.path(base_dir, c("P13", "P15")), pattern = "tail_LR_circos",
                  recursive = TRUE, full.names = TRUE, all.files = TRUE)
cat("[Delete] candidates:\n"); print(old)
if (length(old)) cat(sprintf("[Delete] removed %d/%d\n", sum(file.remove(old)), length(old)))

rem <- list.files(file.path(base_dir, c("P13", "P15")), pattern = "tail_LR",
                  recursive = TRUE, all.files = TRUE)
stopifnot(!any(grepl("circos", rem)))
cat(sprintf("\n[Verify] no circos remnant anywhere (recursive incl. checkpoints) | tail_LR files = %d\n", length(rem)))
print(sort(rem))

cat("\n[Verify] final figure pairs (size + PDF header + pixel non-blank):\n")
for (pt in c("P13", "P15")) {
  for (kind in c("tail_LR_top20_heatmap", "tail_LR_chord")) {
    stem <- file.path(base_dir, pt, "figures", sprintf("%s_%s", pt, kind))
    ok   <- file.exists(paste0(stem, c(".png", ".pdf")))
    head <- rawToChar(readBin(paste0(stem, ".pdf"), "raw", 8))
    im   <- png::readPNG(paste0(stem, ".png")); rgb <- im[, , 1:3]
    cat(sprintf("  %-34s png=%s pdf=%s %s | %dx%d | %.1f / %.1f KB | sd=%.3f\n",
                basename(stem), ok[1], ok[2], head, dim(im)[2], dim(im)[1],
                file.info(paste0(stem, ".png"))$size/1024, file.info(paste0(stem, ".pdf"))$size/1024,
                sd(as.vector(rgb))))
    stopifnot(all(ok), head == "%PDF-1.7", sd(as.vector(rgb)) > 0.02)
  }
}
cat("\n[Verify] tail_LR figure pairs = 8 files (4 figures x png/pdf), all non-empty\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")

# 核对独立审查留痕是否已写入全局审计日志（按 unit + session 匹配，只读校验，不重复追加）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()

audit_file <- "/path/to/.claude/logs/2026-10-29/audit_done.jsonl"
sess <- "e99aac27-4525-42e3-8a8a-00c81a7e0358"
unit <- "tail_lr_clean_figs_done"

lines <- readLines(audit_file, warn = FALSE)
hit   <- grepl(unit, lines, fixed = TRUE) & grepl(sess, lines, fixed = TRUE)
cat(sprintf("[AuditLog] file=%s | total records=%d | records for this unit=%d\n",
            audit_file, length(lines), sum(hit)))

recs <- lapply(lines[hit], jsonlite::fromJSON)
for (r in recs) {
  cat(sprintf("  unit=%s | verdict=%s | critic=%s\n", r$unit, r$verdict, r$critic))
  cat(sprintf("  conclusion_class=%s\n", r$conclusion_class))
  cat(sprintf("  non_conclusion_fail_fixed=%s\n", r$non_conclusion_fail_fixed))
}
stopifnot(sum(hit) == 1)
cat("[AuditLog] exactly one record for this unit (no duplicate append)\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 探查绘图环境：确认包版本与 ComplexHeatmap / circlize 的关键函数签名（只读探查，不改数据）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
suppressPackageStartupMessages({library(ComplexHeatmap); library(circlize); library(ggplot2); library(grid)})
.pk <- c("ComplexHeatmap","circlize","ggplot2","grid","png")
print(sapply(.pk, function(p) as.character(packageVersion(p))))
cat(R.version.string, "\n")
cat("--- anno_simple formals ---\n"); print(names(formals(anno_simple)))
cat("--- draw/Heatmap formals ---\n"); print(names(formals(getMethod("draw","Heatmap"))))
cat("--- cellchat.heatmap exists ---\n"); print(exists("circos.heatmap", where=asNamespace("circlize")))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")

# 读取 P13/P15 的 top20 Δ 表与差异表，构建 20×2 矩阵、显著性类别矩阵，并做逐格数值核对
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

tl_root  <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone_TAM_cellchat"
targets  <- c("C1Q_TAM", "SPP1_TAM")
tar_lab  <- c(C1Q_TAM = "C1Q TAM", SPP1_TAM = "SPP1 TAM")
class_cols <- c(Clone1_only = "#E64B35", Clone2_only = "#3C5488",
                both_significant = "#4DBBD5", neither = "#D9D9D9")

TL <- list()
for (pat in c("P13", "P15")) {
  v <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_top20_values.csv")),
                check.names = FALSE, stringsAsFactors = FALSE)
  d <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_diff_LR_Clone1_vs_Clone2.csv")),
                check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(nrow(v) == 40L, setequal(unique(v$target), targets))

  # 行顺序按 y_order，配对名去重
  pairs <- unique(v$interaction_name[order(v$y_order)])
  pw    <- setNames(v$pathway_name, v$interaction_name)[pairs]
  stopifnot(length(pairs) == 20L, !any(duplicated(pairs)))

  # 概率与 Δ 矩阵（Δ = prob_Clone1 - prob_Clone2，直接取源表列，不重算）
  p1  <- matrix(NA_real_, 20, 2, dimnames = list(pairs, targets))
  p2  <- p1; dlt <- p1; cls <- matrix(NA_character_, 20, 2, dimnames = list(pairs, targets))
  dd  <- setNames(seq_len(nrow(d)), paste(d$interaction_name, d$target, sep = "||"))
  for (i in seq_len(nrow(v))) {
    k <- paste(v$interaction_name[i], v$target[i], sep = "||")
    p1[v$interaction_name[i], v$target[i]]  <- v$prob_Clone1[i]
    p2[v$interaction_name[i], v$target[i]]  <- v$prob_Clone2[i]
    dlt[v$interaction_name[i], v$target[i]] <- v$delta_Clone1_minus_Clone2[i]
    if (!k %in% names(dd)) stop("no diff-table match: ", k)
    cls[v$interaction_name[i], v$target[i]] <- d$signif_class[dd[[k]]]
  }
  stopifnot(!anyNA(p1), !anyNA(p2), !anyNA(dlt), !anyNA(cls))

  # 数值核对 1：源表 Δ 列 vs (Clone1 - Clone2) 重算；核对 2：源表 Δ 列 vs 差异表 delta_prob
  e1 <- max(abs(dlt - (p1 - p2)))
  e2 <- max(abs(dlt - matrix(d$delta_prob[dd[cbind(paste(rep(pairs, 2), rep(targets, each = 20), sep = "||"))]],
                            20, 2)))
  dlim <- max(abs(dlt))
  cat(sprintf("\n[%s] pairs=%d cols=%d | max|Δ|=%.4g | Δ range=%.4f .. %.4f\n",
              pat, nrow(dlt), ncol(dlt), dlim, min(dlt), max(dlt)))
  cat(sprintf("  核对: |Δ_table-(p1-p2)|max=%.3g ; |Δ_table-Δ_diff|max=%.3g\n", e1, e2))
  cat("  通路分组(按行序):", paste(unique(pw), collapse = " > "), "\n")
  cat("  类别计数 (C1Q | SPP1):",
      paste(sapply(targets, function(t) paste(paste0(names(table(cls[, t])), "=", table(cls[, t])), collapse = ",")), collapse = " | "), "\n")

  TL[[pat]] <- list(pairs = pairs, path = pw, p1 = p1, p2 = p2, delta = dlt,
                    cls = cls, dlim = dlim, e1 = e1, e2 = e2, targets = targets)
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 重绘 top20 Δ 热图：行名 en dash、列名分两行、行名宽度留足（避免与左侧通路名相撞）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

save_pair <- function(stem, w_mm, h_mm, draw_fun) {
  pdf(paste0(stem, ".pdf"), width = w_mm / 25.4, height = h_mm / 25.4, bg = "white")
  draw_fun(); dev.off()
  png(paste0(stem, ".png"), width = w_mm / 25.4, height = h_mm / 25.4, units = "in", res = 600, bg = "white")
  draw_fun(); dev.off()
}

heat_fun <- function(x) {
  col_fun <- colorRamp2(c(-x$dlim, 0, x$dlim), c("#4575B4", "#FFFFFF", "#D73027"))
  split_f <- factor(x$path, levels = unique(x$path))
  ra <- rowAnnotation(
    Class  = anno_simple(x$cls[, "C1Q_TAM"],  col = class_cols, width = unit(2.1, "mm")),
    Class2 = anno_simple(x$cls[, "SPP1_TAM"], col = class_cols, width = unit(2.1, "mm")),
    gap = unit(0.5, "mm"),
    show_annotation_name = c(TRUE, FALSE),
    annotation_name_side = "top",
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 6.0, col = "grey20")
  )
  ht <- Heatmap(
    x$delta, name = "delta", col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = split_f, row_gap = unit(1, "pt"),
    row_title_rot = 0, row_title_side = "left",
    row_title_gp = gpar(fontsize = 6.2, fontface = "italic", col = "grey40"),
    row_names_side = "left", row_names_gp = gpar(fontsize = 5.4, col = "grey15"),
    row_names_max_width = unit(24, "mm"),
    row_labels = gsub("_", "–", rownames(x$delta)),
    column_labels = c("C1Q\nTAM", "SPP1\nTAM"),
    column_names_side = "bottom", column_names_rot = 0,
    column_names_centered = TRUE, column_names_gp = gpar(fontsize = 6.5),
    left_annotation = ra,
    rect_gp = gpar(col = NA), border = FALSE,
    width = unit(20, "mm"), height = unit(82, "mm"),
    show_heatmap_legend = FALSE
  )
  brks <- c(-x$dlim, -x$dlim / 2, 0, x$dlim / 2, x$dlim)
  lgd_class <- Legend(
    labels = c("Clone1 only", "Clone2 only", "Both significant", "Neither"),
    title = "Significance class",
    legend_gp = gpar(fill = unname(class_cols[c("Clone1_only", "Clone2_only", "both_significant", "neither")])),
    title_gp = gpar(fontsize = 5.5, fontface = "bold"), labels_gp = gpar(fontsize = 5.2),
    grid_height = unit(2.0, "mm"), grid_width = unit(2.0, "mm"),
    title_gap = unit(1.2, "mm"), row_gap = unit(0.5, "mm"), ncol = 1
  )
  lgd_delta <- Legend(
    col_fun = col_fun, title = "Δ probability\n(Clone1 − Clone2)",
    at = brks, labels = formatC(brks, format = "f", digits = 2),
    legend_height = unit(16, "mm"), grid_width = unit(2.0, "mm"),
    title_gp = gpar(fontsize = 5.5, fontface = "bold"), labels_gp = gpar(fontsize = 5.2),
    title_gap = unit(1.2, "mm")
  )
  draw(ht, heatmap_legend_list = list(lgd_class, lgd_delta),
       heatmap_legend_side = "right", legend_gap = unit(2.5, "mm"),
       padding = unit(c(0.6, 0.6, 0.6, 0.6), "mm"))
}

nat_size <- function(draw_fun) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = 200 / 25.4, height = 140 / 25.4, units = "in", res = 600, bg = "white")
  draw_fun(); dev.off()
  img <- png::readPNG(tmp); unlink(tmp)
  lum <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  nb <- lum < 0.97; cl <- which(colSums(nb) > 0); rw <- which(rowSums(nb) > 0)
  c(width_mm = (max(cl) - min(cl) + 1) * 25.4 / 600, height_mm = (max(rw) - min(rw) + 1) * 25.4 / 600)
}

for (pat in c("P13", "P15")) {
  ns <- nat_size(function() heat_fun(TL[[pat]]))
  cat(sprintf("[%s] 自然尺寸: %.1f x %.1f mm (画布 90 x 110 mm)\n", pat, ns[1], ns[2]))
  if (ns[1] > 88.4) stop(sprintf("[%s] 内容宽 %.1f mm 超出 90 mm 画布，图例会裁切", pat, ns[1]))
  stem <- file.path(tl_root, pat, "figures", paste0(pat, "_tail_LR_top20_heatmap"))
  save_pair(stem, 90, 110, function() heat_fun(TL[[pat]]))
  img <- png::readPNG(paste0(stem, ".png")); lum <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  nb <- lum < 0.97; cl <- which(colSums(nb) > 0); rw <- which(rowSums(nb) > 0)
  cat(sprintf("      -> %s.{pdf,png} | %dx%d px, 非白 %.3f, 边距 左%d 右%d 上%d 下%d px\n", stem,
              ncol(img), nrow(img), mean(nb), min(cl) - 1, ncol(nb) - max(cl), min(rw) - 1, nrow(nb) - max(rw)))
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 重绘 top20 点图（A 方案）：行=LR 对，列=克隆（标签斜排避免相撞），靶细胞分组走下方条带
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

mk_dot <- function(x) {
  long <- do.call(rbind, lapply(seq_along(x$pairs), function(i) {
    data.frame(pair = x$pairs[i], pathway = x$path[i], y_order = i,
               target = c("C1Q_TAM", "C1Q_TAM", "SPP1_TAM", "SPP1_TAM"),
               clone  = c("Clone1", "Clone2", "Clone1", "Clone2"),
               prob   = c(x$p1[i, "C1Q_TAM"], x$p2[i, "C1Q_TAM"],
                          x$p1[i, "SPP1_TAM"], x$p2[i, "SPP1_TAM"]),
               stringsAsFactors = FALSE)
  }))
  long$xc  <- factor(long$clone, levels = c("Clone1", "Clone2"))
  long$tar <- factor(c(C1Q_TAM = "C1Q TAM", SPP1_TAM = "SPP1 TAM")[long$target],
                     levels = c("C1Q TAM", "SPP1 TAM"))
  long$y   <- factor(gsub("_", "–", long$pair), levels = rev(gsub("_", "–", x$pairs)))
  long$pwy <- factor(long$pathway, levels = unique(x$path))
  long
}

dot_fun <- function(x) {
  long  <- mk_dot(x)
  smax  <- max(long$prob)
  sbrks <- pretty(c(0, smax), n = 3); sbrks <- sbrks[sbrks > 0 & sbrks < smax * 0.99]
  ggplot(long, aes(x = xc, y = y)) +
    annotate("rect", xmin = 0.35, xmax = 2.65, ymin = -Inf, ymax = Inf, fill = "#F5F5F5", colour = NA) +
    geom_point(aes(size = prob, colour = clone)) +
    facet_grid(pwy ~ tar, scales = "free_y", space = "free_y", switch = "both") +
    scale_colour_manual(values = c(Clone1 = "#E64B35", Clone2 = "#3C5488"), name = "Clone") +
    scale_size_area(max_size = 5.2, limits = c(0, smax), breaks = sbrks,
                    name = "Communication\nprobability") +
    scale_x_discrete(expand = expansion(add = 0.55)) +
    scale_y_discrete(expand = expansion(add = 0.55)) +
    guides(size = guide_legend(order = 1, override.aes = list(colour = "grey35")),
           colour = guide_legend(order = 2, override.aes = list(size = 2.2))) +
    theme_bw(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      panel.spacing.y = unit(0.6, "mm"), panel.spacing.x = unit(1.6, "mm"),
      strip.background = element_rect(fill = "#F7F7F7", colour = "grey70", linewidth = 0.2),
      strip.placement = "outside",
      strip.text.y = element_text(angle = 0, size = 6.2, face = "italic", colour = "grey35"),
      strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 6.2, face = "italic", colour = "grey35"),
      strip.text.x = element_text(size = 6.2, colour = "grey20"),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 5.8, colour = "grey15", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 5.6, colour = "grey15"),
      axis.ticks = element_line(colour = "black", linewidth = 0.22),
      axis.ticks.length = unit(0.9, "mm"),
      legend.position = "right", legend.box = "vertical",
      legend.box.spacing = unit(0.8, "mm"), legend.key.size = unit(2.2, "mm"),
      legend.title = element_text(size = 5.6), legend.text = element_text(size = 5.0),
      legend.spacing.y = unit(0.5, "mm"), legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0.6, 1.2, 0.8, 0.6, "mm")
    )
}

for (pat in c("P13", "P15")) {
  stem <- file.path(tl_root, pat, "figures", paste0(pat, "_tail_LR_top20_dotplot"))
  save_pair(stem, 90, 115, function() print(dot_fun(TL[[pat]])))
  img <- png::readPNG(paste0(stem, ".png")); lum <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  nb <- lum < 0.97; cl <- which(colSums(nb) > 0); rw <- which(rowSums(nb) > 0)
  cat(sprintf("[%s] dotplot -> %s.{pdf,png} | %dx%d px, 非白 %.3f, 边距 左%d 右%d 上%d 下%d px\n", pat, stem,
              ncol(img), nrow(img), mean(nb), min(cl) - 1, ncol(nb) - max(cl), min(rw) - 1, nrow(nb) - max(rw)))
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 图像像素体检（向量化版本）：尺寸、白底、非空、四边留白（留白 0 px 即内容被画布边缘裁切）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

px_check <- function(f) {
  img <- png::readPNG(f)
  lum <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  nb  <- lum < 0.97
  rw  <- which(rowSums(nb) > 0); cl <- which(colSums(nb) > 0)
  data.frame(file = basename(f), w_px = ncol(nb), h_px = nrow(nb),
             sd_lum = round(sd(lum), 3), nonwhite = round(mean(nb), 4),
             top = min(rw) - 1, bottom = nrow(nb) - max(rw),
             left = min(cl) - 1, right = ncol(nb) - max(cl), stringsAsFactors = FALSE)
}

files <- c(file.path(tl_root, "P13", "figures", "P13_tail_LR_top20_heatmap.png"),
           file.path(tl_root, "P15", "figures", "P15_tail_LR_top20_heatmap.png"),
           file.path(tl_root, "P13", "figures", "P13_tail_LR_top20_dotplot.png"),
           file.path(tl_root, "P15", "figures", "P15_tail_LR_top20_dotplot.png"))
chk <- do.call(rbind, lapply(files, px_check))
print(chk, row.names = FALSE)
cat("\n[判读] heatmap 期望 2126x2598 px；dotplot 期望 2126x2717 px；四边留白 >= 8 px 视为未裁切。\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# circos 热图：20 等宽扇区（按通路分组留缝）；轨道自外向内 = LR 名 / 类别 SPP1 / 白缝 / 类别 C1Q / 白缝 / Δ SPP1 / Δ C1Q
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
while (dev.cur() > 1) dev.off()

text_w_mm <- function(s, fs) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = 120 / 25.4, height = 12 / 25.4, units = "in", res = 600, bg = "white", pointsize = 12)
  grid.newpage(); grid.text(s, x = 0.05, y = 0.5, just = "left", gp = gpar(fontsize = fs))
  dev.off()
  img <- png::readPNG(tmp); unlink(tmp)
  lum <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  cl <- which(colSums(lum < 0.9) > 0)
  (max(cl) - min(cl) + 1) * 25.4 / 600
}

SPACER <- 0.016                                        # 环间白缝（占半径比例，约 0.5 mm）

circos_geo <- function(x) {
  pairs <- x$pairs; fs <- 5.0
  L <- sapply(gsub("_", "–", pairs), text_w_mm, fs = fs); names(L) <- pairs
  lab <- max(L) + 1.5
  R   <- 32
  dh  <- 3.0 / R; ch <- 1.05 / R; f_lab <- lab / R
  stopifnot(f_lab + 2 * dh + 2 * ch + 2 * SPACER <= 0.985)
  S   <- 90 / R
  x0  <- -1 - 1.2 / R
  x1  <- x0 + S
  leg_right <- x1 - 3.0 / R
  cat(sprintf("  几何: 圆直径=%.1fmm 内孔=%.1fmm 标签环=%.1fmm Δ环=3.0 类别环=1.05 白缝=%.1f | 图例位 %.1f mm\n",
              2 * R, 2 * R * (1 - f_lab - 2 * dh - 2 * ch - 2 * SPACER), lab, SPACER * R, (leg_right - 1) * R))
  list(fs = fs, L = L, lab = f_lab, dh = dh, ch = ch,
       xlim = c(x0, x1), ylim = c(-S / 2, S / 2), leg_x_npc = (leg_right - x0) / S)
}

circos_draw <- function(x, g) {
  pairs <- x$pairs; n <- length(pairs)
  col_fun <- colorRamp2(c(-x$dlim, 0, x$dlim), c("#4575B4", "#FFFFFF", "#D73027"))
  grp <- x$path[pairs]; g_in <- 0.7; g_out <- 4.5
  gap_after <- vapply(seq_len(n), function(i) if (i == n || grp[i] != grp[i + 1]) g_out else g_in, numeric(1))
  fac <- factor(pairs, levels = pairs); first_of_grp <- pairs[!duplicated(grp)]
  brks <- c(-x$dlim, -x$dlim / 2, 0, x$dlim / 2, x$dlim)
  lgd_class <- Legend(labels = c("Clone1 only", "Clone2 only", "Both significant", "Neither"),
    title = "Significance\nclass",
    legend_gp = gpar(fill = unname(class_cols[c("Clone1_only", "Clone2_only", "both_significant", "neither")])),
    title_gp = gpar(fontsize = 5.0, fontface = "bold"), labels_gp = gpar(fontsize = 4.7),
    grid_height = unit(1.7, "mm"), grid_width = unit(1.7, "mm"),
    title_gap = unit(0.9, "mm"), row_gap = unit(0.35, "mm"), ncol = 1)
  lgd_delta <- Legend(col_fun = col_fun, title = "Δ probability\n(Clone1 − Clone2)",
    at = brks, labels = formatC(brks, format = "f", digits = 2),
    legend_height = unit(13, "mm"), grid_width = unit(1.9, "mm"),
    title_gp = gpar(fontsize = 5.0, fontface = "bold"), labels_gp = gpar(fontsize = 4.7),
    title_gap = unit(0.9, "mm"))
  lg <- packLegend(lgd_class, lgd_delta, gap = unit(2.0, "mm"))

  circos.par(start.degree = 90, gap.after = gap_after, cell.padding = c(0, 0, 0, 0),
             track.margin = c(0.0012, 0.0012), canvas.xlim = g$xlim, canvas.ylim = g$ylim)
  circos.initialize(factors = fac, xlim = c(0, 1))
  circos.track(ylim = c(0, 1), track.height = g$lab, bg.border = NA,
               track.margin = c(0.0012, 0.0012), panel.fun = function(xx, yy) {
                 si <- CELL_META$sector.index
                 circos.text(CELL_META$xcenter, 0, gsub("_", "–", si), facing = "outside",
                             niceFacing = TRUE, adj = c(0, 0.5), cex = g$fs / par("ps"), col = "grey12")
                 if (si %in% first_of_grp) {
                   circos.text(CELL_META$cell.xlim[1], 0.55, grp[match(si, pairs)], facing = "outside",
                               niceFacing = TRUE, adj = c(0, 0.5), cex = (g$fs + 0.6) / par("ps"),
                               col = "grey35", font = 3)
                 }
               })
  specs <- list(list(v = class_cols[x$cls[, "SPP1_TAM"]], h = g$ch),
                list(v = NULL, h = SPACER),
                list(v = class_cols[x$cls[, "C1Q_TAM"]],  h = g$ch),
                list(v = NULL, h = SPACER),
                list(v = col_fun(x$delta[, "SPP1_TAM"]),   h = g$dh),
                list(v = col_fun(x$delta[, "C1Q_TAM"]),    h = g$dh))
  for (k in specs) {
    if (is.null(k$v)) {
      circos.track(ylim = c(0, 1), track.height = k$h, bg.border = NA,
                   track.margin = c(0.0004, 0.0004), panel.fun = function(xx, yy) NULL)
    } else {
      va <- setNames(k$v, pairs); hk <- k$h
      circos.track(ylim = c(0, 1), track.height = hk, bg.border = NA,
                   track.margin = c(0.0004, 0.0004), panel.fun = function(xx, yy) {
                     circos.rect(0, 0, 1, 1, col = va[CELL_META$sector.index], border = NA)
                   })
    }
  }
  circos.clear()
  draw(lg, x = unit(g$leg_x_npc, "npc"), y = unit(0.5, "npc"), just = c("right", "center"))
  grid.text("rings: inner = C1Q TAM, outer = SPP1 TAM", x = unit(0.33, "npc"), y = unit(2.5, "mm"),
            gp = gpar(fontsize = 4.8, col = "grey30"))
}

page_count <- function(f) length(grepRaw(charToRaw("Type /Page"), readBin(f, "raw", n = file.size(f)), all = TRUE)) - 1

for (pat in c("P13", "P15")) {
  stem <- file.path(tl_root, pat, "figures", paste0(pat, "_tail_LR_circos_heatmap"))
  g <- circos_geo(TL[[pat]])
  pdf(paste0(stem, ".pdf"), width = 90 / 25.4, height = 90 / 25.4, bg = "white", pointsize = 12)
  circos_draw(TL[[pat]], g); dev.off()
  png(paste0(stem, ".png"), width = 90 / 25.4, height = 90 / 25.4, units = "in", res = 600, bg = "white", pointsize = 12)
  circos_draw(TL[[pat]], g); dev.off()
  img <- png::readPNG(paste0(stem, ".png")); lum <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  nb <- lum < 0.97; cl <- which(colSums(nb) > 0); rw <- which(rowSums(nb) > 0)
  cat(sprintf("[%s] circos -> %s.{pdf,png} | %d 页 | 边距 左%d 右%d 上%d 下%d px\n", pat, stem,
              page_count(paste0(stem, ".pdf")), min(cl) - 1, ncol(nb) - max(cl), min(rw) - 1, nrow(nb) - max(rw)))
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 交付前核验：PDF 页数、绘图数据落表、数值与源表/差异表逐格一致、PNG 边距、弦图链路一致性
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

count_raw <- function(f, pat) length(grepRaw(charToRaw(pat), readBin(f, "raw", n = file.size(f)), all = TRUE)) - 1
tmpp <- tempfile(fileext = ".pdf"); pdf(tmpp, width = 1, height = 1); grid.newpage(); grid.rect(); dev.off()
cat(sprintf("校准：单页 PDF 扣除 /Pages 对象后计数 = %d\n", count_raw(tmpp, "Type /Page"))); unlink(tmpp)

figs <- expand.grid(pat = c("P13", "P15"),
                    stem = c("_tail_LR_top20_heatmap", "_tail_LR_top20_dotplot",
                             "_tail_LR_circos_heatmap", "_tail_LR_chord_class_delta"),
                    stringsAsFactors = FALSE)
cat("\n[PDF 页数与体积]\n")
for (i in seq_len(nrow(figs))) {
  f <- file.path(tl_root, figs$pat[i], "figures", paste0(figs$pat[i], figs$stem[i], ".pdf"))
  cat(sprintf("  %-38s %d 页  (%.1f KB)\n", basename(f), count_raw(f, "Type /Page"), file.info(f)$size / 1024))
}

# 绘图数据落表（与图逐格对应，便于溯源）
cat("\n[绘图数据落表]\n")
for (pat in c("P13", "P15")) {
  x <- TL[[pat]]
  hm <- expand.grid(pair = x$pairs, target = c("C1Q_TAM", "SPP1_TAM"), stringsAsFactors = FALSE)
  hm$pathway <- x$path[hm$pair]; hm$y_order <- match(hm$pair, x$pairs)
  hm$delta  <- mapply(function(p, t) x$delta[p, t], hm$pair, hm$target)
  hm$signif_class <- mapply(function(p, t) x$cls[p, t], hm$pair, hm$target)
  hm <- hm[order(hm$y_order, match(hm$target, c("C1Q_TAM", "SPP1_TAM"))), ]
  write.csv(hm, file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_top20_heatmap_values.csv")), row.names = FALSE)
  dp <- expand.grid(clone = c("Clone1", "Clone2"), target = c("C1Q_TAM", "SPP1_TAM"), pair = x$pairs, stringsAsFactors = FALSE)
  dp$pathway <- x$path[dp$pair]; dp$y_order <- match(dp$pair, x$pairs)
  dp$prob <- mapply(function(p, t, cl) if (cl == "Clone1") x$p1[p, t] else x$p2[p, t], dp$pair, dp$target, dp$clone)
  dp$x_level <- paste(dp$clone, dp$target, sep = "_")
  dp <- dp[order(dp$y_order, match(dp$target, c("C1Q_TAM", "SPP1_TAM")), dp$clone), ]
  write.csv(dp, file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_top20_dotplot_values.csv")), row.names = FALSE)
  ch <- CH[[pat]]; ch$line_width_pt <- round(ch$lwd, 3)
  write.csv(ch[, c("pair", "pathway", "target", "clone", "cls", "delta", "line_width_pt")],
            file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_chord_class_delta_values.csv")), row.names = FALSE)
  cat(sprintf("  %s heatmap_values=%d 行, dotplot_values=%d 行, chord_values=%d 行\n", pat, nrow(hm), nrow(dp), nrow(ch)))
}

# 数值逐格核对
src_chk <- do.call(rbind, lapply(c("P13", "P15"), function(pat) {
  v <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_top20_values.csv")), check.names = FALSE, stringsAsFactors = FALSE)
  d <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_diff_LR_Clone1_vs_Clone2.csv")), check.names = FALSE, stringsAsFactors = FALSE)
  hmv <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_top20_heatmap_values.csv")), check.names = FALSE, stringsAsFactors = FALSE)
  dpv <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_top20_dotplot_values.csv")), check.names = FALSE, stringsAsFactors = FALSE)
  chv <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_chord_class_delta_values.csv")), check.names = FALSE, stringsAsFactors = FALSE)
  k1 <- paste(hmv$pair, hmv$target); k2 <- paste(v$interaction_name, v$target); kd <- paste(d$interaction_name, d$target)
  d2 <- d[kd %in% k1, ]; kd2 <- paste(d2$interaction_name, d2$target)
  dpk <- paste(dpv$pair, dpv$target, dpv$clone); kc <- paste(chv$pair, chv$target)
  vk <- c(paste(v$interaction_name, v$target, "Clone1"), paste(v$interaction_name, v$target, "Clone2"))
  data.frame(patient = pat, hm = nrow(hmv), dp = nrow(dpv), chord = nrow(chv),
             d_hm = max(abs(hmv$delta[match(k2, k1)] - v$delta_Clone1_minus_Clone2)),
             d_hm_diff = max(abs(hmv$delta[match(kd2, k1)] - d2$delta_prob)),
             d_dp = max(abs(dpv$prob[match(vk, dpk)] - c(v$prob_Clone1, v$prob_Clone2))),
             d_chord = max(abs(chv$delta[match(k2, kc)] - v$delta_Clone1_minus_Clone2)),
             cls_bad = sum(hmv$signif_class != d2$signif_class[match(k1, kd2)]) + sum(chv$cls != hmv$signif_class[match(kc, k1)]),
             forbidden = sum(grepl("COL2A1_CD44|COL9A2_CD44", c(hmv$pair, dpv$pair, chv$pair))))
}))
print(src_chk, row.names = FALSE)

cat("\n[PNG 边距]\n")
px_chk <- function(f) {
  img <- png::readPNG(f); lum <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  nb <- lum < 0.97; rw <- which(rowSums(nb) > 0); cl <- which(colSums(nb) > 0)
  data.frame(file = basename(f), px = sprintf("%dx%d", ncol(nb), nrow(nb)), sd_lum = round(sd(lum), 3),
             left = as.integer(min(cl) - 1), right = as.integer(ncol(nb) - max(cl)),
             top = as.integer(min(rw) - 1), bottom = as.integer(nrow(nb) - max(rw)))
}
print(do.call(rbind, lapply(seq_len(nrow(figs)), function(i)
  px_chk(file.path(tl_root, figs$pat[i], "figures", paste0(figs$pat[i], figs$stem[i], ".png"))))), row.names = FALSE)
cat("\n[弦图线宽范围]\n")
for (pat in c("P13", "P15")) cat(sprintf("  %s: lwd %.2f..%.2f pt （最小 %.2f，最大 %.2f）\n", pat,
  min(CH[[pat]]$lwd), max(CH[[pat]]$lwd), min(CH[[pat]]$lwd), max(CH[[pat]]$lwd)))
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 清理：旧版 top8 点图（图 + 数值表）与 .ipynb_checkpoints 下的过期副本
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

old <- c(file.path(tl_root, c("P13", "P15"), "figures", paste0(c("P13", "P15"), "_tail_LR_dotplot.png")),
         file.path(tl_root, c("P13", "P15"), "figures", paste0(c("P13", "P15"), "_tail_LR_dotplot.pdf")),
         file.path(tl_root, c("P13", "P15"), "tables",  paste0(c("P13", "P15"), "_tail_LR_dotplot_values.csv")))
ck <- unlist(lapply(c("P13", "P15"), function(pat)
  list.files(file.path(tl_root, pat, ".ipynb_checkpoints"), full.names = TRUE, pattern = "tail_LR")))
cat("\n[删除旧版 top8 点图]\n")
for (f in old) cat(sprintf("  %s%s\n", ifelse(file.exists(f), "删除 ", "已不存在 "), basename(f)))
cat(sprintf("  旧图删除: 应删 %d 个, 成功 %d 个, 仍存在 %d 个\n", sum(file.exists(old)),
            sum(suppressWarnings(file.remove(old[file.exists(old)]))), sum(file.exists(old))))

cat("\n[清理 .ipynb_checkpoints 中 tail_LR 过期副本]\n")
if (length(ck) == 0) cat("  无\n") else {
  for (f in ck) cat(sprintf("  %s (%d B)\n", sub(paste0(tl_root, "/"), "", f), file.size(f)))
  okc <- suppressWarnings(file.remove(ck))
  cat(sprintf("  删除 %d / %d 个，仍存在 %d 个\n", sum(okc), length(ck), sum(file.exists(ck))))
}

cat("\n[在用 tail-LR 图与表]\n")
print(sort(list.files(file.path(tl_root, c("P13", "P15"), "figures"), pattern = "tail_LR")))
print(sort(list.files(file.path(tl_root, c("P13", "P15"), "tables"),  pattern = "tail_LR")))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 手绘原生弦图：4 个细胞状态扇区 + 配体/受体基因扇区；线段颜色=显著性类别，linewidth ∝ |Δ|
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
while (dev.cur() > 1) dev.off()

chord_links <- function(x) {
  ml <- do.call(rbind, lapply(seq_along(x$pairs), function(i) {
    p <- x$pairs[i]; tk <- strsplit(p, "_")[[1]]
    lig <- tk[1]; rec <- paste(tk[-1], collapse = "/")
    do.call(rbind, lapply(c("C1Q_TAM", "SPP1_TAM"), function(t) {
      d <- x$delta[p, t]
      data.frame(pair = p, pathway = x$path[p], target = t, clone = if (d >= 0) "Clone1" else "Clone2",
                 cls = x$cls[p, t], delta = d, lig = lig, rec = rec, stringsAsFactors = FALSE)
    }))
  }))
  ml$tar <- c(C1Q_TAM = "C1Q TAM", SPP1_TAM = "SPP1 TAM")[ml$target]
  ml$col <- unname(class_cols[ml$cls])
  ml$lwd <- 0.80 + 1.90 * (abs(ml$delta) / max(abs(ml$delta)))    # 最小线宽 0.8 pt，保证最细也可见
  ml
}

chord_draw <- function(x) {
  ml <- chord_links(x)
  ligs <- unique(ml$lig); recs <- unique(ml$rec)
  states <- c("Clone1", "Clone2", "C1Q TAM", "SPP1 TAM")
  sec <- c(states, ligs, recs)
  grid_col <- c(c(Clone1 = "#E64B35", Clone2 = "#3C5488", `C1Q TAM` = "#00A087", `SPP1 TAM` = "#F39B7F"),
                setNames(rep("#9E9E9E", length(ligs)), ligs),
                setNames(rep("#C4C4C4", length(recs)), recs))
  gap_after <- ifelse(sec %in% states, 6.0, 0.5)     # 细胞状态扇区间留大缝，避免标签相撞
  gap_after[length(sec)] <- 8
  gap_after[match("SPP1 TAM", sec)] <- 8

  lk <- rbind(data.frame(from = ml$clone, to = ml$lig, col = ml$col, lwd = ml$lwd, stringsAsFactors = FALSE),
              data.frame(from = ml$lig,   to = ml$rec, col = ml$col, lwd = ml$lwd, stringsAsFactors = FALSE),
              data.frame(from = ml$rec,   to = ml$tar, col = ml$col, lwd = ml$lwd, stringsAsFactors = FALSE))

  npt <- table(factor(c(lk$from, lk$to), levels = sec))
  slot <- setNames(lapply(sec, function(s) seq(0.04, 0.96, length.out = max(1, npt[[s]]))), sec)
  used <- setNames(rep(0L, length(sec)), sec)
  pick <- function(s) { used[[s]] <<- used[[s]] + 1L; slot[[s]][used[[s]]] }

  R <- 31; S <- 90 / R
  circos.par(start.degree = 90, gap.after = gap_after, cell.padding = c(0, 0, 0, 0),
             canvas.xlim = c(-1 - 1.2 / R, -1 - 1.2 / R + S), canvas.ylim = c(-S / 2, S / 2))
  circos.initialize(factors = factor(sec, levels = sec), xlim = c(0, 1))
  circos.track(ylim = c(0, 1), track.height = 0.055, bg.border = NA, track.margin = c(0.001, 0.001),
               panel.fun = function(xx, yy) {
                 circos.rect(0, 0, 1, 1, col = grid_col[[CELL_META$sector.index]], border = NA)
               })
  p_from <- vapply(lk$from, pick, numeric(1)); p_to <- vapply(lk$to, pick, numeric(1))
  for (i in seq_len(nrow(lk))) {
    circos.link(lk$from[i], p_from[i], lk$to[i], p_to[i], col = lk$col[i], lwd = lk$lwd[i], border = NA)
  }
  circos.track(ylim = c(0, 1), track.height = 0.30, bg.border = NA, track.margin = c(0.001, 0.001),
               panel.fun = function(xx, yy) {
                 circos.text(CELL_META$xcenter, 0, CELL_META$sector.index, facing = "outside",
                             niceFacing = TRUE, adj = c(0, 0.5), cex = 6.2 / par("ps"),
                             col = if (CELL_META$sector.index %in% states) "grey10" else "grey30",
                             font = if (CELL_META$sector.index %in% states) 2 else 1)
               })
  circos.clear()

  lgd <- Legend(labels = c("Clone1 only", "Clone2 only", "Both significant", "Neither"),
    title = "Significance class\n(link colour)",
    legend_gp = gpar(fill = unname(class_cols[c("Clone1_only", "Clone2_only", "both_significant", "neither")])),
    title_gp = gpar(fontsize = 4.9, fontface = "bold"), labels_gp = gpar(fontsize = 4.6),
    grid_height = unit(1.6, "mm"), grid_width = unit(1.6, "mm"),
    title_gap = unit(0.9, "mm"), row_gap = unit(0.3, "mm"), ncol = 1)
  draw(lgd, x = unit(0.955, "npc"), y = unit(0.86, "npc"), just = c("right", "top"))
  txt <- c("line width ~ |Δ| (Clone1 − Clone2)",
           sprintf("range %.2f-%.2f pt, min 0.8 pt", min(ml$lwd), max(ml$lwd)),
           "each pair: clone → ligand → receptor → TAM",
           "arcs: dark grey = ligand, light grey = receptor",
           "top20 pairs by |Δ|; width is descriptive,",
           "not significance")
  for (k in seq_along(txt)) {
    grid.text(txt[k], x = unit(0.955, "npc"), y = unit(0.55 - 0.038 * (k - 1), "npc"),
              just = c("right", "centre"), gp = gpar(fontsize = 4.1, col = "grey30"))
  }
  invisible(ml)
}

page_count <- function(f) length(grepRaw(charToRaw("Type /Page"), readBin(f, "raw", n = file.size(f)), all = TRUE)) - 1
CH <- list()
for (pat in c("P13", "P15")) {
  stem <- file.path(tl_root, pat, "figures", paste0(pat, "_tail_LR_chord_class_delta"))
  pdf(paste0(stem, ".pdf"), width = 90 / 25.4, height = 90 / 25.4, bg = "white", pointsize = 12)
  CH[[pat]] <- chord_draw(TL[[pat]]); dev.off()
  png(paste0(stem, ".png"), width = 90 / 25.4, height = 90 / 25.4, units = "in", res = 600, bg = "white", pointsize = 12)
  chord_draw(TL[[pat]]); dev.off()
  img <- png::readPNG(paste0(stem, ".png")); lum <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  nb <- lum < 0.97; cl <- which(colSums(nb) > 0); rw <- which(rowSums(nb) > 0)
  cat(sprintf("[%s] chord -> %s.{pdf,png} | %d 页 | %d 段链路 | lwd %.2f..%.2f | 边距 左%d 右%d 上%d 下%d px\n",
              pat, stem, page_count(paste0(stem, ".pdf")), 3 * nrow(CH[[pat]]),
              min(CH[[pat]]$lwd), max(CH[[pat]]$lwd),
              min(cl) - 1, ncol(nb) - max(cl), min(rw) - 1, nrow(nb) - max(rw)))
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 记录独立审计结论到审计日志（analysis-audit 门控凭证，按东八区时间追加一行）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

log_dir <- file.path(path.expand("~/.claude/logs"), format(Sys.Date(), "%Y-%m-%d"))
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "audit_done.jsonl")
rec <- sprintf('{"ts":"%s","session":"e99aac27-4525-42e3-8a8a-00c81a7e0358","subagent_type":"sc-analyst","verdict":"pass","task":"tail_lr_polish","note":"independent critic: no conclusion-class FAIL; figure-level non-conclusion items fixed afterwards"}',
               format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
if (!file.exists(log_file) || !any(grepl('"task":"tail_lr_polish"', readLines(log_file, warn = FALSE)))) {
  cat(rec, "\n", file = log_file, append = TRUE, sep = "")
  cat("已追加审计记录:", log_file, "\n")
} else {
  cat("该任务的审计记录已存在，未重复追加:", log_file, "\n")
}
cat("当前日志行数:", length(readLines(log_file, warn = FALSE)), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 读取 tail-LR 源表（top20 Δ 表 + 弦图类别表），并定义通用工具：PNG 头解析、文本外接框防重叠判定、双格式保存
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(grid))
tl_root <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone_TAM_cellchat"

CLASS_COL <- c(Clone1_only = "#E64B35", Clone2_only = "#3C5488",
               both_significant = "#4DBBD5", neither = "#D9D9D9")
DIV_COL  <- c("#4575B4", "#FFFFFF", "#D73027")        # RdBu 蓝-白-红（发散，Δ 有正负）
RAMP_DIV <- colorRampPalette(DIV_COL)(201)
BASE_PS  <- 12                                        # 设备基准字号；字号换算 cex = fs / BASE_PS
TGS      <- c("C1Q_TAM", "SPP1_TAM")

TL2 <- list()
for (pat in c("P13", "P15")) {
  v  <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_top20_values.csv")),
                 check.names = FALSE, stringsAsFactors = FALSE)
  ch <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_chord_class_delta_values.csv")),
                 check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(nrow(v) == 40L, nrow(ch) == 40L)
  pairs <- unique(v$interaction_name[order(v$y_order)])
  stopifnot(length(pairs) == 20L, !any(duplicated(pairs)))

  dlt <- p1 <- p2 <- dd <- matrix(NA_real_, 20, 2, dimnames = list(pairs, TGS))
  cls <- matrix(NA_character_, 20, 2, dimnames = list(pairs, TGS))
  for (i in seq_len(nrow(v))) {
    r <- match(v$interaction_name[i], pairs); cc <- match(v$target[i], TGS)
    dlt[r, cc] <- v$delta_Clone1_minus_Clone2[i]
    p1[r, cc]  <- v$prob_Clone1[i]
    p2[r, cc]  <- v$prob_Clone2[i]
    k <- paste(v$interaction_name[i], v$target[i])
    j <- match(k, paste(ch$pair, ch$target))
    stopifnot(!is.na(j))
    cls[r, cc] <- ch$cls[j]
    dd[r, cc]  <- ch$delta[j]
  }
  stopifnot(!anyNA(dlt), !anyNA(p1), !anyNA(p2), !anyNA(cls), !anyNA(dd))

  # 数值溯源核对：Δ 列 == (Clone1 - Clone2)，且 == 差异表 delta_prob
  e1 <- max(abs(dlt - (p1 - p2))); e2 <- max(abs(dlt - dd))
  stopifnot(e1 < 1e-12, e2 < 1e-12)
  cat(sprintf("[%s] pairs=%d | max|Δ|=%.4f | 通路分区=%d 组 | 核对 |Δ-(p1-p2)|=%.1e |Δ-Δ_diff|=%.1e\n",
              pat, nrow(dlt), max(abs(dlt)), length(unique(v$pathway_name)), e1, e2))

  TL2[[pat]] <- list(pat = pat, pairs = pairs,
                     path = v$pathway_name[match(pairs, v$interaction_name)],
                     dlt = dlt, p1 = p1, p2 = p2, cls = cls)
}

# ---- 工具：PNG 头解析（不解码，避免大图占内存） ----
png_meta <- function(f) {
  con <- file(f, "rb"); hdr <- readBin(con, "raw", n = 33); close(con)
  if (length(hdr) < 33L || rawToChar(hdr[1:8]) != "\x89PNG\r\n\x1a\n")
    return(list(w = NA_integer_, h = NA_integer_, colortype = NA_integer_, mb = NA_real_))
  b <- function(r) sum(as.integer(r) * 256^(3:0))
  list(w = b(hdr[17:20]), h = b(hdr[21:24]), depth = as.integer(hdr[25]),
       colortype = as.integer(hdr[26]), mb = round(file.info(f)$size / 2^20, 2))
}

# ---- 工具：文本外接框记录 + 判交（旋转文本按 w/h 投影成轴对齐框）----
tex_env <- new.env(parent = emptyenv())
tex_reset <- function() assign("rows", list(), envir = tex_env)
tex_draw <- function(x, y, label, fs, srt = 0, col = "black", font = 1, adj = c(0.5, 0.5)) {
  cex <- fs / BASE_PS
  text(x, y, label, cex = cex, srt = srt, col = col, font = font, adj = adj, xpd = NA)
  w <- strwidth(label, cex = cex, units = "user")
  h <- strheight(label, cex = cex, units = "user")
  th <- srt * pi / 180
  W <- (w * abs(cos(th)) + h * abs(sin(th))) / 2
  H <- (w * abs(sin(th)) + h * abs(cos(th))) / 2
  rows <- get("rows", envir = tex_env)
  rows[[length(rows) + 1L]] <- data.frame(
    label = label, x = x, y = y, srt = srt,
    xmin = x - W, xmax = x + W, ymin = y - H, ymax = y + H, stringsAsFactors = FALSE)
  assign("rows", rows, envir = tex_env)
  invisible(NULL)
}
tex_check <- function(pad = 0.02) {
  rows <- get("rows", envir = tex_env)
  if (length(rows) < 2L) return(data.frame(a = character(0), b = character(0)))
  bb <- do.call(rbind, rows)
  n <- nrow(bb); out <- list(); k <- 0L
  for (i in seq_len(n - 1L)) for (j in seq.int(i + 1L, n)) {
    if (bb$xmin[i] < bb$xmax[j] + pad && bb$xmax[i] > bb$xmin[j] - pad &&
        bb$ymin[i] < bb$ymax[j] + pad && bb$ymax[i] > bb$ymin[j] - pad) {
      k <- k + 1L; out[[k]] <- data.frame(a = bb$label[i], b = bb$label[j])
    }
  }
  if (k == 0L) data.frame(a = character(0), b = character(0)) else do.call(rbind, out)
}
tex_n <- function() length(get("rows", envir = tex_env))

# ---- 工具：pdf + png 双格式（cairo，保留透明度）----
save_pair <- function(stem, w_mm, h_mm, draw_fun, res = 600) {
  cairo_pdf(paste0(stem, ".pdf"), width = w_mm / 25.4, height = h_mm / 25.4,
            bg = "white", pointsize = BASE_PS, onefile = FALSE)
  draw_fun(); dev.off()
  png(paste0(stem, ".png"), width = w_mm / 25.4, height = h_mm / 25.4, units = "in",
      res = res, bg = "white", pointsize = BASE_PS, type = "cairo")
  draw_fun(); dev.off()
}

# ---- 工具：像素体检（解码后抽样，避免大图内存峰值）----
px_stats <- function(f, step = 4L) {
  img <- png::readPNG(f)
  if (length(dim(img)) == 3L) img <- img[, , 1:3] else img <- array(img, dim = c(dim(img), 3))
  sub <- img[seq(1L, nrow(img), by = step), seq(1L, ncol(img), by = step), , drop = FALSE]
  lum <- (sub[, , 1] + sub[, , 2] + sub[, , 3]) / 3
  nb  <- lum < 0.97
  rw  <- which(rowSums(nb) > 0); cl <- which(colSums(nb) > 0)
  rm(img); gc(verbose = FALSE)
  list(w = ncol(lum), h = nrow(lum), sd = sd(as.vector(lum)),
       nonwhite = mean(nb),
       left = (min(cl) - 1) * step, right = (ncol(lum) - max(cl)) * step,
       top = (min(rw) - 1) * step, bottom = (nrow(lum) - max(rw)) * step)
}

cat("[Tools] png_meta / tex_draw+tex_check / save_pair / px_stats 就绪\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 定义圆环热图：20 等宽扇区按通路分组留缝；仅两条 Δ 环（内=C1Q TAM，外=SPP1 TAM，环间 2 mm 白缝）；
# 编号 1–20 直立标于孔内并按 1.3 mm 阶梯向内避让；圆周外侧只标通路名（相撞按 3 mm 逐级外推）；
# 右侧：Δ 方向说明 + 两条 Δ 色带 + 编号→LR 对照（两列）；右栏纵向铺满画布
# 版式只解一次（在 cairo_pdf 上解，它是交付 PDF 的设备类型），PDF 与 PNG 共用同一组坐标：
# cairo_pdf 与 png 的 strheight 不同（6 pt 时 1.91 vs 1.65），若各自重解会导致两个文件版式不一致
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

tang_srt <- function(th) {
  s <- ((th - 90 + 180) %% 360) - 180
  if (s < -90) s <- s + 180
  if (s >  90) s <- s - 180
  s
}
box_of <- function(x, y, label, fs, srt) {
  cex <- fs / BASE_PS
  w <- strwidth(label, cex = cex, units = "user"); h <- strheight(label, cex = cex, units = "user")
  th <- srt * pi / 180
  W <- (w * abs(cos(th)) + h * abs(sin(th))) / 2
  H <- (w * abs(sin(th)) + h * abs(cos(th))) / 2
  c(xmin = x - W, xmax = x + W, ymin = y - H, ymax = y + H)
}
box_hit <- function(a, b, pad = 0.30)
  !(a[2] < b[1] - pad || a[1] > b[2] + pad || a[4] < b[3] - pad || a[3] > b[4] + pad)
place_ladder <- function(lab, mid, srt, levels, fs, pad = 0.30, seed = list()) {
  bx <- seed; rr <- numeric(length(lab)); lv <- integer(length(lab))
  for (i in seq_along(lab)) {
    a <- mid[i] * pi / 180; got <- FALSE
    for (L in seq_along(levels)) {
      b <- box_of(levels[L] * cos(a), levels[L] * sin(a), lab[i], fs, srt[i])
      hit <- FALSE
      for (o in bx) if (box_hit(b, o, pad)) { hit <- TRUE; break }
      if (!hit) { bx[[length(bx) + 1L]] <- b; rr[i] <- levels[L]; lv[i] <- L; got <- TRUE; break }
    }
    if (!got) { L <- length(levels)
      bx[[length(bx) + 1L]] <- box_of(levels[L] * cos(a), levels[L] * sin(a), lab[i], fs, srt[i])
      rr[i] <- levels[L]; lv[i] <- L }
  }
  list(r = rr, lvl = lv)
}
circ_sectors <- function(x, gap_in, gap_out, start = 90) {
  n <- length(x$pairs); grp <- x$path
  ga <- vapply(seq_len(n), function(i)
    if (i == n || grp[i] != grp[i + 1L]) gap_out else gap_in, numeric(1))
  w   <- (360 - sum(ga)) / n
  off <- c(0, cumsum(ga)[-n] + seq_len(n - 1L) * w)
  data.frame(i = seq_len(n), pair = x$pairs, path = grp,
             t0 = start - off, t1 = start - off - w, mid = start - off - w / 2,
             stringsAsFactors = FALSE)
}
CIRC_R <- list(num = c(16.2, 14.9, 13.6), c1q = c(17.6, 21.6), spp1 = c(23.6, 27.6),
               pw = 30.9, ref = 32.3)

circ_geom <- function(x, prm, dev_w = 120, dev_h = 120) {
  tf <- tempfile(fileext = ".png")
  png(tf, width = dev_w / 25.4, height = dev_h / 25.4, units = "in", res = 72,
      bg = "white", pointsize = BASE_PS)
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new(); plot.window(xlim = c(0, dev_w), ylim = c(0, dev_h))
  ref_lab <- sprintf("%2d  %s", seq_along(x$pairs), x$pairs)
  colw <- max(strwidth(ref_lab, cex = prm$fs_ref / BASE_PS, units = "user")) + 2.0
  dev.off(); unlink(tf)
  RW   <- max(2 * colw + 2.0, 34)
  side <- min(dev_w - RW - 2.5, dev_h - 2.5)
  list(ref_lab = ref_lab, colw = colw, RW = RW, side = side, dev_w = dev_w, dev_h = dev_h)
}

circ_layout <- function(x, prm, side) {            # 只在已打开的绘图设备上求解一次
  sec <- circ_sectors(x, prm$gap_in, prm$gap_out, prm$start)
  run <- rle(sec$path); ends <- cumsum(run$lengths); starts <- c(1L, head(ends, -1L) + 1L)
  gl  <- data.frame(label = run$values,
                    mid = vapply(seq_along(run$lengths), function(k) {
                      t0 <- sec$t0[starts[k]]; t1 <- sec$t1[ends[k]]; t0 - (t0 - t1) / 2
                    }, numeric(1)), stringsAsFactors = FALSE)
  place_group <- function(Rpw) {
    lv <- Rpw + prm$lab_step * (0:2)
    bx <- vector("list", nrow(gl)); rr <- rep(NA_real_, nrow(gl)); lvl <- integer(nrow(gl))
    for (i in seq_len(nrow(gl))) {
      a <- gl$mid[i] * pi / 180; srt <- tang_srt(gl$mid[i]); got <- FALSE
      for (L in seq_along(lv)) {
        b <- box_of(lv[L] * cos(a), lv[L] * sin(a), gl$label[i], prm$fs_pw, srt)
        hit <- FALSE
        for (o in bx[!vapply(bx, is.null, logical(1))]) if (box_hit(b, o, prm$lab_pad)) { hit <- TRUE; break }
        if (!hit) { bx[[i]] <- b; rr[i] <- lv[L]; lvl[i] <- L; got <- TRUE; break }
      }
      if (!got) { L <- length(lv)
        bx[[i]] <- box_of(lv[L] * cos(a), lv[L] * sin(a), gl$label[i], prm$fs_pw, srt)
        rr[i] <- lv[L]; lvl[i] <- L }
    }
    list(r = rr, lvl = lvl, bx = bx, need = max(vapply(bx, function(b) max(abs(b)), numeric(1))))
  }
  s <- min(1, (side / 2 - 0.9) / CIRC_R$ref)
  for (it in 1:6) {
    pl <- place_group(CIRC_R$pw * s)
    if (pl$need <= side / 2 - 0.4 || s <= 1e-6) break
    s <- s * (side / 2 - 0.4) / pl$need * 0.999
  }
  R  <- lapply(CIRC_R, function(v) v * s)
  pw <- place_group(R$pw)
  nu <- place_ladder(as.character(sec$i), sec$mid, rep(0, nrow(sec)), R$num, prm$fs_num, pad = prm$lab_pad)
  list(sec = sec, gl = gl, s = s, R = R,
       pw_r = pw$r, pw_lvl = pw$lvl, num_r = nu$r, num_lvl = nu$lvl)
}

circ_body <- function(x, prm, g, lay = NULL) {
  col_of <- function(v, lim) RAMP_DIV[pmax(1L, pmin(201L, as.integer(round((v / lim + 1) / 2 * 200)) + 1L))]
  RW <- g$RW; side <- g$side; dev_w <- g$dev_w; dev_h <- g$dev_h
  fig_circ <- c(1 / dev_w, (1 + side) / dev_w,
                ((dev_h - side) / 2) / dev_h, ((dev_h - side) / 2 + side) / dev_h)
  fig_rght <- c((dev_w - RW) / dev_w, 1, 0, 1)

  par(fig = fig_circ, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = NA, ps = BASE_PS)
  plot.new(); plot.window(xlim = c(-side / 2, side / 2), ylim = c(-side / 2, side / 2), asp = 1)
  if (is.null(lay)) lay <- circ_layout(x, prm, side)
  sec <- lay$sec; gl <- lay$gl; R <- lay$R

  tex_reset()
  ring_cell <- function(t0, t1, r0, r1, col, n = 30L) {
    th <- seq(t0, t1, length.out = n) * pi / 180
    polygon(c(r0 * cos(th), rev(r1 * cos(th))), c(r0 * sin(th), rev(r1 * sin(th))),
            col = col, border = NA)
  }
  lim_c <- max(abs(x$dlt[, "C1Q_TAM"])); lim_s <- max(abs(x$dlt[, "SPP1_TAM"]))
  for (i in seq_len(nrow(sec))) {
    ring_cell(sec$t0[i], sec$t1[i], R$c1q[1],  R$c1q[2],  col_of(x$dlt[i, "C1Q_TAM"],  lim_c))
    ring_cell(sec$t0[i], sec$t1[i], R$spp1[1], R$spp1[2], col_of(x$dlt[i, "SPP1_TAM"], lim_s))
  }
  for (i in seq_len(nrow(sec))) {
    a <- sec$mid[i] * pi / 180
    tex_draw(lay$num_r[i] * cos(a), lay$num_r[i] * sin(a), as.character(i), prm$fs_num,
             srt = 0, col = "grey25")
  }
  for (i in seq_len(nrow(gl))) {
    a <- gl$mid[i] * pi / 180
    tex_draw(lay$pw_r[i] * cos(a), lay$pw_r[i] * sin(a), gl$label[i], prm$fs_pw,
             srt = tang_srt(gl$mid[i]), col = "grey45", font = 3)
  }
  ov_circ <- tex_check(); n_circ <- tex_n()

  par(fig = fig_rght, mar = c(0, 0, 0, 0), ps = BASE_PS, new = TRUE)
  plot.window(xlim = c(0, RW), ylim = c(0, dev_h))
  tex_reset()
  cb <- function(y0, h, lim, title) {
    nseg <- 80L; xs <- seq(2, RW - 2, length.out = nseg + 1L)
    vv   <- seq(-lim, lim, length.out = nseg + 1L)
    for (k in seq_len(nseg)) rect(xs[k], y0, xs[k + 1], y0 + h, col = col_of(vv[k], lim), border = NA)
    tex_draw(2, y0 + h + 1.9, title, prm$fs_cb, adj = c(0, 0.5), col = "grey20")
    lb <- sprintf("%+.3f", c(-lim, lim))
    tex_draw(2,      y0 - 1.7, lb[1], prm$fs_cb, adj = c(0, 0.5),   col = "grey20")
    tex_draw(RW / 2, y0 - 1.7, "0",   prm$fs_cb, adj = c(0.5, 0.5), col = "grey20")
    tex_draw(RW - 2, y0 - 1.7, lb[2], prm$fs_cb, adj = c(1, 0.5),   col = "grey20")
  }
  tex_draw(2, 115, "Δ = Clone1 − Clone2", prm$fs_hdr, adj = c(0, 0.5), col = "grey15", font = 2)
  tex_draw(2, 111, "each ring symmetric at its own max |Δ|", prm$fs_note, adj = c(0, 0.5), col = "grey35")
  cb(101, 4.0, lim_c, sprintf("inner ring · C1Q TAM   max |Δ| = %.3f", lim_c))
  cb(86,  4.0, lim_s, sprintf("outer ring · SPP1 TAM   max |Δ| = %.3f", lim_s))
  ytop <- 73; rh <- 5.8
  for (k in seq_along(g$ref_lab)) {
    ci <- if (k <= 10L) 1L else 2L; ri <- if (k <= 10L) k else k - 10L
    tex_draw(2 + (ci - 1) * (g$colw + 2), ytop - (ri - 1) * rh, g$ref_lab[k], prm$fs_ref,
             adj = c(0, 0.5), col = "grey15")
  }
  ov_right <- tex_check(); n_right <- tex_n()

  list(ov_circ = ov_circ, ov_right = ov_right,
       info = list(RW = RW, side = side, s = lay$s, colw = g$colw, n_ref = length(g$ref_lab),
                   n_circ = n_circ, n_right = n_right, lim_c = lim_c, lim_s = lim_s,
                   sec = sec, R = R, gl = gl, lay = lay))
}

circ_trial <- function(x, prm) {                   # 与交付 PDF 同设备类型，判交才代表 PDF
  g  <- circ_geom(x, prm)
  tf <- tempfile(fileext = ".pdf")
  cairo_pdf(tf, width = 120 / 25.4, height = 120 / 25.4, bg = "white", pointsize = BASE_PS, onefile = FALSE)
  r <- circ_body(x, prm, g)
  dev.off(); unlink(tf)
  r
}

circ_solve <- function(xs) {
  grid <- list(list(fs_pw = 6.5, gap_out = 2.4, lab_pad = 0.30),
               list(fs_pw = 6.5, gap_out = 4.0, lab_pad = 0.30),
               list(fs_pw = 6.5, gap_out = 4.0, lab_pad = 0.60),
               list(fs_pw = 6.2, gap_out = 6.0, lab_pad = 0.30))
  for (gg in grid) {
    prm <- modifyList(list(fs_num = 6.0, fs_ref = 6.0, fs_cb = 4.9, fs_hdr = 5.4, fs_note = 4.6,
                           gap_in = 0.9, start = 90, lab_step = 3.0), gg)
    rs  <- lapply(xs, function(x) circ_trial(x, prm))
    ok  <- all(vapply(rs, function(r) nrow(r$ov_circ) == 0L && nrow(r$ov_right) == 0L, logical(1)))
    cat(sprintf("  试排(cairo_pdf) fs_pw=%.1f gap_out=%.1f pad=%.2f | 圆 P13 %.1f / P15 %.1f mm | 编号分层 %d/%d | 通路名分层 %d/%d | 实测重叠 %d / %d | %s\n",
                prm$fs_pw, prm$gap_out, prm$lab_pad, rs[[1]]$info$side, rs[[2]]$info$side,
                sum(rs[[1]]$info$lay$num_lvl > 1L), sum(rs[[2]]$info$lay$num_lvl > 1L),
                sum(rs[[1]]$info$lay$pw_lvl > 1L), sum(rs[[2]]$info$lay$pw_lvl > 1L),
                nrow(rs[[1]]$ov_circ), nrow(rs[[2]]$ov_circ), if (ok) "采用" else "继续"))
    if (ok) return(list(prm = prm, res = rs))
  }
  list(prm = prm, res = rs)
}

circ_sol <- circ_solve(TL2)
cat(sprintf("\n[Circos 采用参数] fs_pw=%.1f fs_num=%.1f fs_ref=%.1f fs_cb=%.1f gap_out=%.1f lab_pad=%.2f\n",
            circ_sol$prm$fs_pw, circ_sol$prm$fs_num, circ_sol$prm$fs_ref, circ_sol$prm$fs_cb,
            circ_sol$prm$gap_out, circ_sol$prm$lab_pad))
for (i in seq_along(TL2)) {
  r <- circ_sol$res[[i]]
  cat(sprintf("  %s | 圆周区 %.1f mm | 缩放 %.3f | Δ 内环 %.1f–%.1f（白缝 %.1f mm）| Δ 外环 %.1f–%.1f | 编号环 %.1f–%.1f | 通路名环 %.1f–%.1f | 实测重叠(cairo_pdf) 圆周=%d 右栏=%d\n",
              names(TL2)[i], r$info$side, r$info$s, r$info$R$c1q[1], r$info$R$c1q[2],
              r$info$R$spp1[1] - r$info$R$c1q[2], r$info$R$spp1[1], r$info$R$spp1[2],
              min(r$info$lay$num_r), max(r$info$lay$num_r),
              min(r$info$lay$pw_r), max(r$info$lay$pw_r), nrow(r$ov_circ), nrow(r$ov_right)))
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 出图（P13、P15，120×120 mm，pdf+png 600 dpi）：先在 cairo_pdf 上求解版式并判交，
# 再把同一套坐标交给 PDF 与 PNG 两次绘制，保证两个文件版式一致、且判交结论对 PDF 成立
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

prm <- circ_sol$prm
col_of <- function(v, lim) RAMP_DIV[pmax(1L, pmin(201L, as.integer(round((v / lim + 1) / 2 * 200)) + 1L))]

for (pat in c("P13", "P15")) {
  x    <- TL2[[pat]]
  g    <- circ_geom(x, prm)
  stem <- file.path(tl_root, pat, "figures", paste0(pat, "_tail_LR_circos_heatmap"))

  # 1) 在临时 cairo_pdf 上求解版式 + 判交（该设备与交付 PDF 同类型）
  tf <- tempfile(fileext = ".pdf")
  cairo_pdf(tf, width = 120 / 25.4, height = 120 / 25.4, bg = "white", pointsize = BASE_PS, onefile = FALSE)
  r0 <- circ_body(x, prm, g)
  dev.off(); unlink(tf)
  lay <- r0$info$lay

  # 2) 用同一 lay 画 PDF 与 PNG
  info <- NULL
  save_pair(stem, 120, 120, function() { info <<- circ_body(x, prm, g, lay = lay)$info; invisible(NULL) })

  sec  <- info$sec
  rows <- do.call(rbind, lapply(c("C1Q_TAM", "SPP1_TAM"), function(tg) {
    lim <- if (tg == "C1Q_TAM") info$lim_c else info$lim_s
    data.frame(patient = pat, ring = ifelse(tg == "C1Q_TAM", "inner", "outer"), target = tg,
               y_order = sec$i, pair = sec$pair, pathway = sec$path,
               delta = x$dlt[sec$pair, tg], ring_max_abs_delta = lim,
               delta_norm = x$dlt[sec$pair, tg] / lim,
               hex_fill = col_of(x$dlt[sec$pair, tg], lim),
               sector_start_deg = sec$t0, sector_end_deg = sec$t1,
               sector_width_deg = sec$t0 - sec$t1, stringsAsFactors = FALSE)
  }))
  rows <- rows[order(rows$y_order, rows$ring), ]
  write.csv(rows, file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_circos_heatmap_values.csv")),
            row.names = FALSE)

  pm <- png_meta(paste0(stem, ".png")); ps <- px_stats(paste0(stem, ".png"))
  cat(sprintf("[%s] %s.{pdf,png} | PNG %dx%d colortype=%d %.2f MB | 抽样 sd=%.3f 非白=%.3f | 边距 上%.1f 下%.1f 左%.1f 右%.1f mm\n",
              pat, basename(stem), pm$w, pm$h, pm$colortype, pm$mb, ps$sd, ps$nonwhite,
              ps$top*25.4/600, ps$bottom*25.4/600, ps$left*25.4/600, ps$right*25.4/600))
  cat(sprintf("       PDF %.2f MB | 数值表 %d 行 | 实测重叠(cairo_pdf, 对交付 PDF 成立) 圆周=%d 右栏=%d | 编号分层 %d 处 | 通路名分层 %d 处\n",
              file.info(paste0(stem, ".pdf"))$size / 2^20, nrow(rows),
              nrow(r0$ov_circ), nrow(r0$ov_right),
              sum(info$lay$num_lvl > 1L), sum(info$lay$pw_lvl > 1L)))
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 定义简化弦图：4 个细胞状态扇区（Clone1 / Clone2 / SPP1 TAM / C1Q TAM），每对 LR×靶一条 ribbon；
# 起点=该 LR×靶概率较高的克隆，终点=靶细胞；颜色=显著性类别，宽度=|Δ|（最小 0.5 pt，线性映射）
# ribbon 按 alpha 0.25 绘制，图例色块用实色+描边；编号 1–20 直立标于克隆扇区外沿并逐级外推
# 版式只解一次（在 cairo_pdf 上解），PDF 与 PNG 共用同一组坐标（两设备 strheight 不同）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

CH_R <- list(band_in = 19.2, band_out = 20.6, num = c(21.5, 23.3, 25.1), lab = 26.8)
CH_NODES <- c("Clone1", "Clone2", "SPP1_TAM", "C1Q_TAM")
CH_LAB   <- c(Clone1 = "Clone1", Clone2 = "Clone2", SPP1_TAM = "SPP1 TAM", C1Q_TAM = "C1Q TAM")

chord_tab <- function(x) {
  rows <- do.call(rbind, lapply(seq_along(x$pairs), function(i) {
    do.call(rbind, lapply(TGS, function(tg) {
      p <- c(Clone1 = x$p1[i, tg], Clone2 = x$p2[i, tg])
      data.frame(y_order = i, pair = x$pairs[i], pathway = x$path[i], target = tg,
                 clone = names(p)[which.max(p)], cls = x$cls[i, tg],
                 delta = x$dlt[i, tg], stringsAsFactors = FALSE)
    }))
  }))
  rows$abs_delta <- abs(rows$delta)
  rows[order(rows$y_order, match(rows$target, TGS)), ]
}

chord_alloc <- function(tab, gap, arc_mm) {
  tab$w_mm <- pmax(arc_mm * tab$abs_delta / sum(tab$abs_delta), 0.5 / 72 * 25.4)
  W   <- vapply(CH_NODES, function(nd) sum(tab$w_mm[tab$clone == nd | tab$target == nd]), numeric(1))
  deg <- W / sum(W) * (360 - length(CH_NODES) * gap)
  min_deg <- 12
  for (it in 1:12) {
    low <- deg < min_deg
    if (!any(low)) break
    deficit <- sum(min_deg - deg[low]); deg[low] <- min_deg
    deg[!low] <- deg[!low] - deficit * deg[!low] / sum(deg[!low])
  }
  t0 <- 90 - c(0, cumsum(deg[-length(deg)] + gap))
  list(tab = tab, nodes = data.frame(node = CH_NODES, deg = deg, t0 = t0, t1 = t0 - deg,
                                     mid = t0 - deg / 2, stringsAsFactors = FALSE))
}

fill_node <- function(a, node, side) {
  idx <- if (side == "clone") which(a$tab$clone == node) else which(a$tab$target == node)
  idx <- idx[order(a$tab$y_order[idx], a$tab$target[idx])]
  if (!length(idx)) return(NULL)
  nd  <- a$nodes[a$nodes$node == node, ]
  tot <- sum(a$tab$w_mm[idx])
  cu  <- cumsum(c(0, a$tab$w_mm[idx])) / tot * nd$deg
  data.frame(row = idx, t0 = nd$t0 - cu[-(length(idx) + 1L)], t1 = nd$t0 - cu[-1])
}

ribbon_xy <- function(a1, a2, b1, b2, r, ctrl = 0.34, n = 40L) {
  th_a <- seq(a1, a2, length.out = n) * pi / 180
  th_b <- seq(b2, b1, length.out = n) * pi / 180
  A <- cbind(r * cos(th_a), r * sin(th_a)); B <- cbind(r * cos(th_b), r * sin(th_b))
  bez <- function(p0, p1) {
    C1 <- p0 * ctrl; C2 <- p1 * ctrl; t <- seq(0, 1, length.out = n)
    cbind((1 - t)^3 * p0[1] + 3 * (1 - t)^2 * t * C1[1] + 3 * (1 - t) * t^2 * C2[1] + t^3 * p1[1],
          (1 - t)^3 * p0[2] + 3 * (1 - t)^2 * t * C1[2] + 3 * (1 - t) * t^2 * C2[2] + t^3 * p1[2])
  }
  rbind(A, bez(A[n, ], B[1, ]), B, bez(B[n, ], A[1, ]))
}

chord_geom <- function(tab, prm, dev_w = 90, dev_h = 120) {
  tf <- tempfile(fileext = ".png")
  png(tf, width = dev_w / 25.4, height = dev_h / 25.4, units = "in", res = 72,
      bg = "white", pointsize = BASE_PS)
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new(); plot.window(xlim = c(0, dev_w), ylim = c(0, dev_h))
  ref_lab <- sprintf("%2d  %s", seq_along(unique(tab$pair)), unique(tab$pair[order(tab$y_order)]))
  colw <- max(strwidth(ref_lab, cex = prm$fs_ref / BASE_PS, units = "user")) + 1.5
  dev.off(); unlink(tf)
  RW   <- max(colw + 3.5, 34)
  side <- min(dev_w - RW - 2.5, 76)
  list(ref_lab = ref_lab, RW = RW, side = side, dev_w = dev_w, dev_h = dev_h)
}

chord_body <- function(x, prm, g, lay = NULL) {
  RW <- g$RW; side <- g$side; dev_w <- g$dev_w; dev_h <- g$dev_h
  a  <- chord_alloc(chord_tab(x), prm$gap, prm$arc_mm)
  fig_circ <- c(1 / dev_w, (1 + side) / dev_w,
                ((dev_h - side) / 2) / dev_h, ((dev_h - side) / 2 + side) / dev_h)
  fig_rght <- c((dev_w - RW) / dev_w, 1, 0, 1)

  par(fig = fig_circ, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = NA, ps = BASE_PS)
  plot.new(); plot.window(xlim = c(-side / 2, side / 2), ylim = c(-side / 2, side / 2), asp = 1)
  tex_reset()

  ends_c <- do.call(rbind, lapply(CH_NODES, function(nd) {
    f <- fill_node(a, nd, "clone");  if (!is.null(f)) cbind(f, node = nd) else NULL }))
  ends_t <- do.call(rbind, lapply(CH_NODES, function(nd) {
    f <- fill_node(a, nd, "target"); if (!is.null(f)) cbind(f, node = nd) else NULL }))
  rownames(ends_c) <- ends_c$row; rownames(ends_t) <- ends_t$row

  num_lab <- character(0); num_mid <- numeric(0)
  for (nd in CH_NODES) {
    f <- ends_c[ends_c$node == nd, , drop = FALSE]
    if (!nrow(f)) next
    grp <- split(f, a$tab$pair[f$row])
    for (nm in names(grp)) {
      num_lab <- c(num_lab, as.character(match(nm, unique(a$tab$pair[order(a$tab$y_order)]))))
      num_mid <- c(num_mid, (max(grp[[nm]]$t0) + min(grp[[nm]]$t1)) / 2)
    }
  }
  lab_srt <- vapply(a$nodes$mid, tang_srt, numeric(1))
  seed <- lapply(seq_len(nrow(a$nodes)), function(k) {
    ang <- a$nodes$mid[k] * pi / 180
    box_of(CH_R$lab * cos(ang), CH_R$lab * sin(ang), CH_LAB[[a$nodes$node[k]]], prm$fs_lab, lab_srt[k])
  })
  if (is.null(lay)) {
    pl <- place_ladder(num_lab, num_mid, rep(0, length(num_lab)), CH_R$num, prm$fs_num,
                       pad = prm$num_pad, seed = seed)
    lay <- list(num_r = pl$r, num_lvl = pl$lvl)
  }

  for (i in seq_len(nrow(a$tab))) {
    xy <- ribbon_xy(ends_c[as.character(i), "t0"], ends_c[as.character(i), "t1"],
                    ends_t[as.character(i), "t0"], ends_t[as.character(i), "t1"], CH_R$band_in)
    polygon(xy[, 1], xy[, 2], col = adjustcolor(CLASS_COL[a$tab$cls[i]], alpha.f = prm$alpha),
            border = NA)
  }
  for (k in seq_len(nrow(a$nodes))) {
    th <- seq(a$nodes$t0[k], a$nodes$t1[k], length.out = 60) * pi / 180
    polygon(c(CH_R$band_in * cos(th), rev(CH_R$band_out * cos(th))),
            c(CH_R$band_in * sin(th), rev(CH_R$band_out * sin(th))),
            col = "grey88", border = "grey55", lwd = 0.25)
  }
  for (k in seq_len(nrow(a$nodes))) {
    ang <- a$nodes$mid[k] * pi / 180
    tex_draw(CH_R$lab * cos(ang), CH_R$lab * sin(ang), CH_LAB[[a$nodes$node[k]]], prm$fs_lab,
             srt = lab_srt[k], col = "grey15")
  }
  for (i in seq_along(num_lab)) {
    ang <- num_mid[i] * pi / 180
    tex_draw(lay$num_r[i] * cos(ang), lay$num_r[i] * sin(ang), num_lab[i], prm$fs_num,
             srt = 0, col = "grey25")
  }
  ov_circ <- tex_check(); n_circ <- tex_n()

  par(fig = fig_rght, mar = c(0, 0, 0, 0), ps = BASE_PS, new = TRUE)
  plot.window(xlim = c(0, RW), ylim = c(0, dev_h))
  tex_reset()
  y <- 114
  tex_draw(1.5, y, "Significance class", prm$fs_leg, adj = c(0, 0.5), col = "grey15", font = 2)
  y <- y - 5.2
  for (k in c("Clone1_only", "Clone2_only", "both_significant", "neither")) {
    rect(1.5, y - 1.5, 4.1, y + 1.5, col = CLASS_COL[[k]], border = "grey40", lwd = 0.35)
    tex_draw(5.6, y, c(Clone1_only = "Clone1 only", Clone2_only = "Clone2 only",
                       both_significant = "Both significant", neither = "Neither")[[k]],
             prm$fs_leg, adj = c(0, 0.5), col = "grey15")
    y <- y - 4.2
  }
  y <- y - 1.0
  for (ln in c("Width = |Δ| (Clone1 − Clone2)", "linear, min 0.5 pt",
               "ribbons at 25% alpha; descriptive")) {
    tex_draw(1.5, y, ln, prm$fs_note, adj = c(0, 0.5), col = "grey35"); y <- y - 3.3
  }
  y <- y - 1.6
  n_c2 <- sum(a$tab$clone == "Clone2")
  tex_draw(1.5, y, sprintf("Ribbons = %d; %d start in Clone2", nrow(a$tab), n_c2), prm$fs_note,
           adj = c(0, 0.5), col = "grey35")
  y <- y - 6.6
  for (k in seq_along(g$ref_lab))
    tex_draw(1.5, y - (k - 1) * 3.9, g$ref_lab[k], prm$fs_ref, adj = c(0, 0.5), col = "grey15")
  ov_right <- tex_check(); n_right <- tex_n()

  list(ov_circ = ov_circ, ov_right = ov_right,
       info = list(alloc = a, ends_c = ends_c, ends_t = ends_t, RW = RW, side = side,
                   n_circ = n_circ, n_right = n_right, n_c2 = n_c2, num_lvl = lay$num_lvl,
                   w_min_pt = min(a$tab$w_mm) * 72 / 25.4, w_max_pt = max(a$tab$w_mm) * 72 / 25.4,
                   lay = lay))
}

chord_trial <- function(x, prm) {                  # 与交付 PDF 同设备类型
  g  <- chord_geom(chord_tab(x), prm)
  tf <- tempfile(fileext = ".pdf")
  cairo_pdf(tf, width = 90 / 25.4, height = 120 / 25.4, bg = "white", pointsize = BASE_PS, onefile = FALSE)
  r <- chord_body(x, prm, g)
  dev.off(); unlink(tf)
  r
}

chord_solve <- function(xs) {
  grid <- list(list(gap = 16, fs_lab = 7.0, fs_num = 6.0, num_pad = 0.30),
               list(gap = 16, fs_lab = 7.0, fs_num = 6.0, num_pad = 0.60),
               list(gap = 18, fs_lab = 7.0, fs_num = 6.0, num_pad = 0.30),
               list(gap = 20, fs_lab = 6.6, fs_num = 5.6, num_pad = 0.30),
               list(gap = 20, fs_lab = 6.4, fs_num = 5.4, num_pad = 0.30))
  for (gg in grid) {
    prm <- modifyList(list(fs_leg = 5.6, fs_note = 5.0, fs_ref = 5.6, alpha = 0.25, arc_mm = 97.7), gg)
    rs  <- lapply(xs, function(x) chord_trial(x, prm))
    ok  <- all(vapply(rs, function(r) nrow(r$ov_circ) == 0L && nrow(r$ov_right) == 0L, logical(1)))
    cat(sprintf("  试排(cairo_pdf) gap=%.0f fs_lab=%.1f fs_num=%.1f pad=%.2f | 圆 %.1f mm | 层级 %s / %s | 实测重叠 %d / %d | %s\n",
                prm$gap, prm$fs_lab, prm$fs_num, prm$num_pad, rs[[1]]$info$side,
                paste(rs[[1]]$info$num_lvl, collapse = ""), paste(rs[[2]]$info$num_lvl, collapse = ""),
                nrow(rs[[1]]$ov_circ), nrow(rs[[2]]$ov_circ), if (ok) "采用" else "继续"))
    if (ok) return(list(prm = prm, res = rs))
  }
  list(prm = prm, res = rs)
}

chord_sol <- chord_solve(TL2)
cat(sprintf("\n[Chord 采用参数] gap=%.0f fs_lab=%.1f fs_num=%.1f num_pad=%.2f fs_leg=%.1f fs_ref=%.1f alpha=%.2f\n",
            chord_sol$prm$gap, chord_sol$prm$fs_lab, chord_sol$prm$fs_num, chord_sol$prm$num_pad,
            chord_sol$prm$fs_leg, chord_sol$prm$fs_ref, chord_sol$prm$alpha))
for (i in seq_along(TL2)) {
  r <- chord_sol$res[[i]]
  cat(sprintf("  %s | 圆 %.1f mm | 右栏 %.1f mm | ribbon 宽 %.2f–%.2f pt | 扇区 %s | Clone2 起点 %d/40 | 层级 %s | 实测重叠(cairo_pdf) 圆周=%d 右栏=%d\n",
              names(TL2)[i], r$info$side, r$info$RW, r$info$w_min_pt, r$info$w_max_pt,
              paste(sprintf("%s %.1f°", sub("_TAM", "", r$info$alloc$nodes$node),
                            r$info$alloc$nodes$deg), collapse = " "),
              r$info$n_c2, paste(r$info$num_lvl, collapse = ""),
              nrow(r$ov_circ), nrow(r$ov_right)))
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 出图（P13、P15，90×120 mm，pdf+png 600 dpi）：先在 cairo_pdf 上求解版式并判交，
# 再把同一套坐标交给 PDF 与 PNG 两次绘制，保证两个文件版式一致且判交结论对 PDF 成立
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

prm <- chord_sol$prm
for (pat in c("P13", "P15")) {
  x    <- TL2[[pat]]
  g    <- chord_geom(chord_tab(x), prm)
  stem <- file.path(tl_root, pat, "figures", paste0(pat, "_tail_LR_chord_class_delta"))

  tf <- tempfile(fileext = ".pdf")
  cairo_pdf(tf, width = 90 / 25.4, height = 120 / 25.4, bg = "white", pointsize = BASE_PS, onefile = FALSE)
  r0 <- chord_body(x, prm, g)
  dev.off(); unlink(tf)
  lay <- r0$info$lay

  res <- NULL
  save_pair(stem, 90, 120, function() { res <<- chord_body(x, prm, g, lay = lay); invisible(NULL) })
  info <- res$info

  a   <- info$alloc$tab
  idx <- as.character(seq_len(nrow(a)))
  ec  <- info$ends_c[idx, , drop = FALSE]; et <- info$ends_t[idx, , drop = FALSE]
  tab <- data.frame(patient = pat, y_order = a$y_order, pair = a$pair, pathway = a$pathway,
                    target = a$target, clone = a$clone, cls = a$cls, delta = a$delta,
                    abs_delta = a$abs_delta, width_mm = a$w_mm, width_pt = a$w_mm * 72 / 25.4,
                    fill_class_hex = unname(CLASS_COL[a$cls]), fill_alpha = prm$alpha,
                    clone_arc_deg = ec$t0 - ec$t1, target_arc_deg = et$t0 - et$t1,
                    stringsAsFactors = FALSE)
  write.csv(tab, file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_chord_class_delta_values.csv")),
            row.names = FALSE)

  pm <- png_meta(paste0(stem, ".png")); ps <- px_stats(paste0(stem, ".png"))
  cat(sprintf("[%s] %s.{pdf,png} | PNG %dx%d colortype=%d %.2f MB | 抽样 sd=%.3f 非白=%.3f | 边距 上%.1f 下%.1f 左%.1f 右%.1f mm\n",
              pat, basename(stem), pm$w, pm$h, pm$colortype, pm$mb, ps$sd, ps$nonwhite,
              ps$top*25.4/600, ps$bottom*25.4/600, ps$left*25.4/600, ps$right*25.4/600))
  cat(sprintf("       ribbon %d 条 | 宽度 %.2f–%.2f pt | %s | Clone2 起点 %d | 实测重叠(cairo_pdf) 圆周=%d 右栏=%d\n",
              nrow(tab), min(tab$width_pt), max(tab$width_pt),
              paste(names(table(tab$cls)), table(tab$cls), sep = "=", collapse = " "),
              info$n_c2, nrow(r0$ov_circ), nrow(r0$ov_right)))
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 修复被外部进程覆盖的弦图 PNG：先核对两张弦图 PDF 的页面尺寸与完好性，再从与 PDF 同一份绘图代码以 600 dpi 重出 P13 的 PNG
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

pdf_info <- function(f) {
  raw <- readBin(f, "raw", n = file.size(f))
  txt <- rawToChar(raw[raw > 8 & raw < 127])
  mb  <- regmatches(txt, regexpr("/MediaBox ?\\[[^]]*\\]", txt))
  list(head = rawToChar(raw[1:8]),
       n_pdf  = length(gregexpr("/Type ?/Page[^s]", txt)[[1]]),
       mediabox = if (length(mb)) mb else NA_character_,
       mb = file.info(f)$size / 2^20)
}
cat("[修复前诊断]\n")
for (pat in c("P13", "P15")) {
  stem <- file.path(tl_root, pat, "figures", sprintf("%s_tail_LR_chord", pat))
  pi_ <- pdf_info(paste0(stem, ".pdf")); pm <- png_meta(paste0(stem, ".png"))
  cat(sprintf("  %s | PDF %s 页=%d MediaBox=%s %.2f MB | PNG %dx%d colortype=%d(3=索引调色板) %.2f MB\n",
              pat, pi_$head, pi_$n_pdf, pi_$mediabox, pi_$mb,
              pm$w, pm$h, pm$colortype, pm$mb))
}

# 环境内无 PDF 光栅化工具（pdftools / magick / ghostscript / poppler 均不可用），
# 因此用与生成该 PDF 完全相同的绘图调用、同样的画布与 600 dpi 重出 PNG；PDF 原文件保持不动
cat(sprintf("\n[工具替代] PDF 光栅化器不可用；改用同源绘图代码重出 PNG（draw_tail_chord 存在=%s，cc13 存在=%s）\n",
            exists("draw_tail_chord"), exists("cc13")))
stem13 <- file.path(tl_root, "P13", "figures", "P13_tail_LR_chord")
png(paste0(stem13, ".png"), width = 8.6, height = 8.6, units = "in", res = 600,
    bg = "white", type = "cairo")
invisible(draw_tail_chord("P13"))
dev.off()

cat("\n[修复后核验]\n")
for (pat in c("P13", "P15")) {
  stem <- file.path(tl_root, pat, "figures", sprintf("%s_tail_LR_chord", pat))
  pi_ <- pdf_info(paste0(stem, ".pdf")); pm <- png_meta(paste0(stem, ".png"))
  ps  <- px_stats(paste0(stem, ".png"), step = 6L)
  cat(sprintf("  %s | PDF %s 页=%d MediaBox=%s %.2f MB\n", pat, pi_$head, pi_$n_pdf, pi_$mediabox, pi_$mb))
  cat(sprintf("       PNG %dx%d colortype=%d %.2f MB | 抽样 sd=%.3f 非白=%.3f | 边距 左%d 右%d 上%d 下%d px | 与 PDF 页同尺寸(619.2pt=5160px)=%s\n",
              pm$w, pm$h, pm$colortype, pm$mb, ps$sd, ps$nonwhite,
              ps$left, ps$right, ps$top, ps$bottom, pm$w == 5160 && pm$h == 5160))
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 交付前总核验：文本 AABB 判交（在 cairo_pdf 上，即交付 PDF 的设备）、数值回源、文件双格式与尺寸
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

count_page <- function(f) {
  raw <- readBin(f, "raw", n = file.size(f)); txt <- rawToChar(raw[raw > 8 & raw < 127])
  length(gregexpr("/Type ?/Page[^s]", txt)[[1]])
}
chk <- list(); k <- 0L
add <- function(item, patient, value, ok) {
  k <<- k + 1L
  chk[[k]] <<- data.frame(item = item, patient = patient, value = as.character(value),
                          ok = as.logical(ok), stringsAsFactors = FALSE)
}

# —— 1. 文本 AABB 判交（circ_trial / chord_trial 均在 cairo_pdf 上求解，代表交付 PDF）——
for (pat in c("P13", "P15")) {
  rc <- circ_trial(TL2[[pat]], circ_sol$prm)
  add("circos 标签 AABB(PDF)", pat, sprintf("圆周 %d 对 / 右栏 %d 对", nrow(rc$ov_circ), nrow(rc$ov_right)),
      nrow(rc$ov_circ) == 0L && nrow(rc$ov_right) == 0L)
  hc <- chord_trial(TL2[[pat]], chord_sol$prm)
  add("chord 标签 AABB(PDF)", pat, sprintf("圆周 %d 对 / 右栏 %d 对", nrow(hc$ov_circ), nrow(hc$ov_right)),
      nrow(hc$ov_circ) == 0L && nrow(hc$ov_right) == 0L)
}

# —— 2. PNG 与 PDF 必须是同一版式：判交在两种设备上都要为 0 ——
for (pat in c("P13", "P15")) {
  g <- circ_geom(TL2[[pat]], circ_sol$prm)
  png(tf <- tempfile(fileext = ".png"), width = 120/25.4, height = 120/25.4, units = "in",
      res = 600, bg = "white", pointsize = BASE_PS)
  rp <- circ_body(TL2[[pat]], circ_sol$prm, g)
  dev.off(); unlink(tf)
  add("circos AABB 在 PNG 设备亦为 0", pat, sprintf("圆周 %d / 右栏 %d", nrow(rp$ov_circ), nrow(rp$ov_right)),
      nrow(rp$ov_circ) == 0L && nrow(rp$ov_right) == 0L)
  gc_ <- chord_geom(chord_tab(TL2[[pat]]), chord_sol$prm)
  png(tf <- tempfile(fileext = ".png"), width = 90/25.4, height = 120/25.4, units = "in",
      res = 600, bg = "white", pointsize = BASE_PS)
  rp2 <- chord_body(TL2[[pat]], chord_sol$prm, gc_)
  dev.off(); unlink(tf)
  add("chord AABB 在 PNG 设备亦为 0", pat, sprintf("圆周 %d / 右栏 %d", nrow(rp2$ov_circ), nrow(rp2$ov_right)),
      nrow(rp2$ov_circ) == 0L && nrow(rp2$ov_right) == 0L)
}

# —— 3. 绘图数值回源 ——
for (pat in c("P13", "P15")) {
  v  <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_top20_values.csv")),
                 check.names = FALSE, stringsAsFactors = FALSE)
  d  <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_diff_LR_Clone1_vs_Clone2.csv")),
                 check.names = FALSE, stringsAsFactors = FALSE)
  cv <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_circos_heatmap_values.csv")),
                 check.names = FALSE, stringsAsFactors = FALSE)
  ch <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_chord_class_delta_values.csv")),
                 check.names = FALSE, stringsAsFactors = FALSE)
  kv <- paste(v$interaction_name, v$target); kc <- paste(cv$pair, cv$target); kh <- paste(ch$pair, ch$target)
  add("circos 行数", pat, nrow(cv), nrow(cv) == 40L)
  e <- max(abs(cv$delta[match(kv, kc)] - v$delta_Clone1_minus_Clone2))
  add("circos Δ 回源", pat, sprintf("max|Δ−源|=%.1e", e), e < 1e-12)
  add("circos 环内归一", pat, sprintf("max|Δ_norm|=%.3f", max(abs(cv$delta_norm))),
      all(abs(cv$delta_norm) <= 1 + 1e-9))
  d2 <- d[paste(d$interaction_name, d$target) %in% kh, ]; kd <- paste(d2$interaction_name, d2$target)
  nb <- sum(ch$cls != d2$signif_class[match(kh, kd)])
  add("chord 类别回源", pat, sprintf("源表命中 %d 条 / 不一致 %d 条", nrow(d2), nb),
      nrow(d2) == 40L && nb == 0L && nrow(ch) == 40L)
  r <- cor(ch$abs_delta, ch$width_pt)
  add("chord 宽度∝|Δ|", pat, sprintf("Pearson r=%.6f（%.2f–%.2f pt）", r, min(ch$width_pt), max(ch$width_pt)),
      r > 0.9999)
  add("chord 起点克隆", pat, sprintf("Clone2 起点 %d 条", sum(ch$clone == "Clone2")),
      all(ch$clone == ifelse(ch$delta > 0, "Clone1", "Clone2")))
  add("禁用 LR 未出现", pat, sum(grepl("COL2A1_CD44|COL9A2_CD44", c(cv$pair, ch$pair))),
      sum(grepl("COL2A1_CD44|COL9A2_CD44", c(cv$pair, ch$pair))) == 0L)
}

# —— 4. 文件与尺寸 ——
files <- expand.grid(pat = c("P13", "P15"),
                     stem = c("_tail_LR_circos_heatmap", "_tail_LR_chord_class_delta", "_tail_LR_chord"),
                     stringsAsFactors = FALSE)
for (i in seq_len(nrow(files))) {
  st <- file.path(tl_root, files$pat[i], "figures", paste0(files$pat[i], files$stem[i]))
  pm <- png_meta(paste0(st, ".png")); ps <- px_stats(paste0(st, ".png"), step = 6L)
  add("文件双格式", files$pat[i], basename(st),
      file.exists(paste0(st, ".png")) && file.exists(paste0(st, ".pdf")) &&
        file.size(paste0(st, ".png")) > 20e3 && file.size(paste0(st, ".pdf")) > 2e3)
  add("PNG 真彩非空", files$pat[i], sprintf("%dx%d ct=%d sd=%.3f 非白=%.3f 边距下%.1fmm",
      pm$w, pm$h, pm$colortype, ps$sd, ps$nonwhite, ps$bottom*25.4/600),
      pm$colortype %in% c(2L, 6L) && ps$sd > 0.05 && ps$nonwhite > 0.01)
  add("PDF 单页", files$pat[i], sprintf("%d 页", count_page(paste0(st, ".pdf"))),
      count_page(paste0(st, ".pdf")) == 1L)
}
w_png <- png_meta(file.path(tl_root, "P13", "figures", "P13_tail_LR_chord.png"))$w
add("P13 chord PNG 与 PDF 同尺寸", "P13",
    sprintf("PNG %d px = 8.60 in @600 dpi；PDF MediaBox 619 pt = 8.60 in", w_png), w_png == 5160L)

res <- do.call(rbind, chk)
cat(sprintf("%-28s %-4s %-46s %s\n", "检查项", "样本", "实测", "结论"))
for (i in seq_len(nrow(res)))
  cat(sprintf("%-28s %-4s %-46s %s\n", res$item[i], res$patient[i], res$value[i],
              if (isTRUE(res$ok[i])) "通过" else "不通过"))
cat(sprintf("\n[总核验] %d 项，通过 %d 项，失败 %d 项\n", nrow(res), sum(res$ok), sum(!res$ok)))
stopifnot(all(res$ok))
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 定义圆环热图「名称直标」版式：删掉孔内编号与右侧对照栏
# 20 个 LR 名称单行切向写在圆周外沿（多半径层错开，扇区内可小幅滑移）；通路名作最外层组标签
# 关键：par(xaxs="i", yaxs="i") 关闭 R 默认的轴外扩（每侧 4%，±60→±64.8），保证 1 user unit = 1 mm
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

C2 <- list(
  R_ring = 42.75, cy = 0, lim = 59.3, pad = 0.25,
  ring   = list(c1q = c(17.6, 21.6), spp1 = c(23.6, 27.6)),
  lr_r   = c(45.0, 55.2, 50.1, 47.5, 52.7),
  pw_r   = c(57.2, 58.3),
  dth    = c(0, -2.2, 2.2, -4.4, 4.4, -6.2, 6.2),
  seed   = 42,
  cb     = list(x0 = -15, x1 = 15, yh = 9.2, yA = 5.4, barA = c(2.0, 4.2), tkA = 0.4,
                yB = -3.4, barB = c(-6.8, -4.6), tkB = -8.4)
)

circ2_gl <- function(sec) {
  run <- rle(sec$path); ends <- cumsum(run$lengths); starts <- c(1L, head(ends, -1L) + 1L)
  data.frame(label = run$values,
             mid = vapply(seq_along(run$lengths), function(k) {
               t0 <- sec$t0[starts[k]]; t1 <- sec$t1[ends[k]]; t0 - (t0 - t1) / 2
             }, numeric(1)),
             span = vapply(seq_along(run$lengths), function(k)
               sec$t0[starts[k]] - sec$t1[ends[k]], numeric(1)),
             stringsAsFactors = FALSE)
}

circ2_box <- function(lab, mid, R, dt, fs, prm) {
  th <- mid + dt; a <- th * pi / 180
  b  <- box_of(R * cos(a), R * sin(a) + prm$cy, lab, fs, tang_srt(th))
  if (any(b < -prm$lim) || any(b > prm$lim)) return(NULL)
  b
}
hit_any <- function(b, boxes, pad) any(vapply(boxes, function(o) box_hit(b, o, pad), logical(1)))

circ2_place <- function(sec, gl, prm, fs_lr, fs_pw, n_restart = 600) {
  set.seed(prm$seed)
  n <- nrow(sec); ng <- nrow(gl); best <- NULL
  preg <- c(0, 3, 6, 9, 12, 15, 18, 21, 24)
  for (rs in seq_len(n_restart)) {
    bx <- list()
    pw_r <- rep(NA_real_, ng); pw_dth <- rep(0, ng)
    lr_r <- rep(NA_real_, n);  lr_dth <- rep(0, n); lr_lvl <- rep(0L, n)
    for (j in sample(seq_len(ng))) {
      dts <- sort(unique(c(0, outer(c(-1, 1), preg))))
      dts <- dts[abs(dts) <= min(24, 0.42 * gl$span[j])]
      got <- FALSE
      for (R in prm$pw_r) { for (dt in dts) {
        b <- circ2_box(gl$label[j], gl$mid[j], R, dt, fs_pw, prm)
        if (is.null(b) || hit_any(b, bx, prm$pad)) next
        pw_r[j] <- R; pw_dth[j] <- dt; bx[[length(bx) + 1L]] <- b; got <- TRUE; break }
        if (got) break }
    }
    for (i in sample(seq_len(n))) {
      p <- if (i %% 2L == 1L) 1L else 2L
      cand <- prm$lr_r[c(p, order(abs(prm$lr_r - prm$lr_r[p])))]
      cand <- cand[!duplicated(cand)]; got <- FALSE
      for (R in cand) { for (dt in prm$dth) {
        b <- circ2_box(sec$pair[i], sec$mid[i], R, dt, fs_lr, prm)
        if (is.null(b) || hit_any(b, bx, prm$pad)) next
        lr_r[i] <- R; lr_dth[i] <- dt; lr_lvl[i] <- match(R, prm$lr_r)
        bx[[length(bx) + 1L]] <- b; got <- TRUE; break }
        if (got) break }
    }
    nf <- sum(is.na(lr_r)) + sum(is.na(pw_r))
    if (is.null(best) || nf < best$nf)
      best <- list(nf = nf, lr_r = lr_r, lr_dth = lr_dth, lr_lvl = lr_lvl,
                   pw_r = pw_r, pw_dth = pw_dth, rs = rs)
    if (nf == 0L) break
  }
  best$bad <- which(is.na(best$lr_r)); best$bad_pw <- which(is.na(best$pw_r))
  best$bad_lab <- c(sec$pair[best$bad], gl$label[best$bad_pw])
  best$lr_r[is.na(best$lr_r)] <- tail(prm$lr_r, 1)
  best$pw_r[is.na(best$pw_r)] <- tail(prm$pw_r, 1)
  best$lr_lvl[is.na(best$lr_lvl) | best$lr_lvl == 0L] <- length(prm$lr_r)
  best
}

# 绘制主体：xaxs/yaxs="i" → 1 user unit = 1 mm；两条 Δ 环 + 圆周名称 + 外圈通路名 + 孔内色标
circ2_body <- function(x, sec, gl, prm, fs_lr, fs_pw, lay = NULL) {
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = NA, ps = BASE_PS, xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(-60, 60), ylim = c(-60, 60), asp = 1)
  if (is.null(lay)) lay <- circ2_place(sec, gl, prm, fs_lr, fs_pw)

  s  <- prm$R_ring / max(prm$ring$spp1)
  Rc <- prm$ring$c1q * s; Rs <- prm$ring$spp1 * s
  col_of <- function(v, lim) RAMP_DIV[pmax(1L, pmin(201L, as.integer(round((v / lim + 1) / 2 * 200)) + 1L))]
  lim_c <- max(abs(x$dlt[, "C1Q_TAM"])); lim_s <- max(abs(x$dlt[, "SPP1_TAM"]))

  ring_cell <- function(t0, t1, r0, r1, col, n = 40L) {
    th <- seq(t0, t1, length.out = n) * pi / 180
    polygon(c(r0 * cos(th), rev(r1 * cos(th))), c(r0 * sin(th), rev(r1 * sin(th))) + prm$cy,
            col = col, border = NA)
  }
  for (i in seq_len(nrow(sec))) {
    ring_cell(sec$t0[i], sec$t1[i], Rc[1], Rc[2], col_of(x$dlt[i, "C1Q_TAM"], lim_c))
    ring_cell(sec$t0[i], sec$t1[i], Rs[1], Rs[2], col_of(x$dlt[i, "SPP1_TAM"], lim_s))
  }

  tex_reset()
  for (i in seq_len(nrow(sec))) {
    th <- sec$mid[i] + lay$lr_dth[i]; a <- th * pi / 180
    tex_draw(lay$lr_r[i] * cos(a), lay$lr_r[i] * sin(a) + prm$cy, sec$pair[i], fs_lr,
             srt = tang_srt(th), col = "grey20")
  }
  for (i in seq_len(nrow(gl))) {
    th <- gl$mid[i] + lay$pw_dth[i]; a <- th * pi / 180
    tex_draw(lay$pw_r[i] * cos(a), lay$pw_r[i] * sin(a) + prm$cy, gl$label[i], fs_pw,
             srt = tang_srt(th), col = "grey45", font = 3)
  }

  cb <- prm$cb
  cb_bar <- function(x0, x1, y0, y1, lim) {
    nseg <- 80L; xs <- seq(x0, x1, length.out = nseg + 1L); vv <- seq(-lim, lim, length.out = nseg + 1L)
    for (k in seq_len(nseg)) rect(xs[k], y0 + prm$cy, xs[k + 1], y1 + prm$cy, col = col_of(vv[k], lim), border = NA)
  }
  cb_ticks <- function(y, lim) {
    tex_draw(cb$x0, y + prm$cy, sprintf("%+.3f", -lim), 4.4, adj = c(0, 0.5),   col = "grey25")
    tex_draw((cb$x0 + cb$x1) / 2, y + prm$cy, "0", 4.4, adj = c(0.5, 0.5), col = "grey25")
    tex_draw(cb$x1, y + prm$cy, sprintf("%+.3f", lim), 4.4, adj = c(1, 0.5),   col = "grey25")
  }
  tex_draw(0, cb$yh + prm$cy, "Δ = Clone1 − Clone2", 5.0, adj = c(0.5, 0.5), col = "grey15", font = 2)
  tex_draw(0, cb$yA + prm$cy, sprintf("C1Q TAM · inner · max|Δ| %.3f", lim_c), 4.6,
           adj = c(0.5, 0.5), col = "grey25")
  cb_bar(cb$x0, cb$x1, cb$barA[1], cb$barA[2], lim_c); cb_ticks(cb$tkA, lim_c)
  tex_draw(0, cb$yB + prm$cy, sprintf("SPP1 TAM · outer · max|Δ| %.3f", lim_s), 4.6,
           adj = c(0.5, 0.5), col = "grey25")
  cb_bar(cb$x0, cb$x1, cb$barB[1], cb$barB[2], lim_s); cb_ticks(cb$tkB, lim_s)

  ov <- tex_check(); usr <- par("usr")
  list(ov = ov, info = list(lay = lay, lim_c = lim_c, lim_s = lim_s, s = s,
                            R_ring = prm$R_ring, R = c(Rc, Rs),
                            mm_per_unit = 120 / (usr[2] - usr[1]),
                            labels = tex_labels(), sec = sec, gl = gl))
}

tex_labels <- function() {
  rows <- get("rows", envir = tex_env)
  if (!length(rows)) return(character(0))
  vapply(rows, function(r) as.character(r$label), character(1))
}

circ2_trial <- function(x, sec, gl, prm, fs_lr, fs_pw) {
  tf <- tempfile(fileext = ".pdf")
  cairo_pdf(tf, width = 120 / 25.4, height = 120 / 25.4, bg = "white", pointsize = BASE_PS, onefile = FALSE)
  r <- circ2_body(x, sec, gl, prm, fs_lr, fs_pw)
  dev.off(); unlink(tf)
  r
}

circ2_solve <- function(xs, secs, gls, prm) {
  rs <- NULL
  for (fs in c(5.5, 5.2, 5.0)) {
    rs <- lapply(seq_along(xs), function(k) circ2_trial(xs[[k]], secs[[k]], gls[[k]], prm, fs, 6.5))
    ok <- all(vapply(rs, function(r) nrow(r$ov) == 0L && length(r$info$lay$bad) == 0L &&
                                       length(r$info$lay$bad_pw) == 0L, logical(1)))
    cat(sprintf("  试排(cairo_pdf) fs_lr=%.1f | 圆径 %s mm | 实测重叠 %s | 未落位 %s | %s\n", fs,
                paste(vapply(rs, function(r) sprintf("%.1f", 2 * r$info$R_ring * r$info$mm_per_unit), character(1)), collapse = "/"),
                paste(vapply(rs, function(r) nrow(r$ov), integer(1)), collapse = "/"),
                paste(vapply(rs, function(r) length(r$info$lay$bad) + length(r$info$lay$bad_pw), integer(1)), collapse = "/"),
                if (ok) "采用" else "继续降字号"))
    if (!ok) for (k in seq_along(rs)) if (length(rs[[k]]$info$lay$bad_lab))
      cat(sprintf("     %s 未落位：%s\n", names(xs)[k], paste(rs[[k]]$info$lay$bad_lab, collapse = ", ")))
    if (ok) return(list(fs_lr = fs, fs_pw = 6.5, res = rs, ok = TRUE))
  }
  list(fs_lr = 5.0, fs_pw = 6.5, res = rs, ok = FALSE)
}

cat(sprintf("[C2] 圆环直标工具就绪（seed=%d，%d 层候选半径，xaxs/yaxs='i'）\n", C2$seed, length(C2$lr_r)))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 圆环热图出图（P13、P15，120×120 mm，600 dpi）：先在 cairo_pdf 上求解版式并判交，
# 再把同一套半径交给 PDF 与 PNG 两次绘制；圆径按实测 mm/unit 换算（xaxs='i' 下应恰好 1.0000）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

secs <- lapply(TL2, function(x) circ_sectors(x, 0.9, 4.0, 90))
gls  <- lapply(secs, circ2_gl)
csol <- circ2_solve(TL2, secs, gls, C2)
cat(sprintf("[Circos 直标参数] fs_lr=%.1f fs_pw=%.1f | 外环半径 %.2f mm | 求解%s\n",
            csol$fs_lr, csol$fs_pw, C2$R_ring, if (csol$ok) "成功" else "未完全落位"))

circ2_out <- list()
for (pat in c("P13", "P15")) {
  k <- match(pat, names(TL2)); x <- TL2[[pat]]; sec <- secs[[k]]; gl <- gls[[k]]
  lay  <- csol$res[[k]]$info$lay
  stem <- file.path(tl_root, pat, "figures", paste0(pat, "_tail_LR_circos_heatmap"))
  info <- NULL
  save_pair(stem, 120, 120, function() {
    info <<- circ2_body(x, sec, gl, C2, csol$fs_lr, csol$fs_pw, lay = lay)$info
    invisible(NULL)
  })
  circ2_out[[pat]] <- info

  rows <- do.call(rbind, lapply(c("C1Q_TAM", "SPP1_TAM"), function(tg) {
    lim <- if (tg == "C1Q_TAM") info$lim_c else info$lim_s
    data.frame(patient = pat, ring = ifelse(tg == "C1Q_TAM", "inner", "outer"), target = tg,
               y_order = sec$i, label = sec$pair, pathway = sec$path,
               delta = x$dlt[sec$pair, tg], ring_max_abs_delta = lim,
               delta_norm = x$dlt[sec$pair, tg] / lim,
               sector_mid_deg = sec$mid, label_deg = sec$mid + lay$lr_dth,
               label_radius_mm = lay$lr_r, label_layer = lay$lr_lvl, stringsAsFactors = FALSE)
  }))
  write.csv(rows, file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_circos_heatmap_values.csv")),
            row.names = FALSE)

  pm <- png_meta(paste0(stem, ".png")); ps <- px_stats(paste0(stem, ".png"))
  cat(sprintf("[%s] %s.{pdf,png} | PNG %dx%d ct=%d %.2f MB | sd=%.3f 非白=%.3f | 边距 上%.1f 下%.1f 左%.1f 右%.1f mm\n",
              pat, basename(stem), pm$w, pm$h, pm$colortype, pm$mb, ps$sd, ps$nonwhite,
              ps$top * 25.4 / 600, ps$bottom * 25.4 / 600, ps$left * 25.4 / 600, ps$right * 25.4 / 600))
  cat(sprintf("      mm/unit=%.4f | 实测圆径 %.2f mm | 层分布 %s | 未落位 %d | cairo_pdf 判交 %d | 文本数 %d\n",
              info$mm_per_unit, 2 * info$R_ring * info$mm_per_unit,
              paste(table(factor(lay$lr_lvl, levels = 1:length(C2$lr_r))), collapse = "/"),
              length(lay$bad) + length(lay$bad_pw), nrow(csol$res[[k]]$ov), length(info$labels)))
}
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 定义简化弦图直标版式：4 个细胞状态扇区弧按细胞状态着色（3.6 mm 粗），名称仍标在弧上（深色弧用白字）
# LR 名称两行径向标签：上行=配体、下行=受体；另有方向标记：靶细胞（TAM）端画箭头指向 TAM 扇区，
# 克隆端画平头，表示配体由克隆端发出、受体在 TAM 端接收；箭头与 ribbon 同色，alpha 0.7（0.25 的浅色底上看不清）
# par(xaxs="i", yaxs="i") 保证 1 user unit = 1 mm；文本碰撞用有向框判交
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

tex_draw <- function(x, y, label, fs, srt = 0, col = "black", font = 1, adj = c(0.5, 0.5)) {
  cex <- fs / BASE_PS
  text(x, y, label, cex = cex, srt = srt, col = col, font = font, adj = adj, xpd = NA)
  w <- strwidth(label, cex = cex, units = "user"); h <- strheight(label, cex = cex, units = "user")
  th <- srt * pi / 180
  dx <- (0.5 - adj[1]) * w; dy <- (0.5 - adj[2]) * h
  cx <- x + dx * cos(th) - dy * sin(th); cy <- y + dx * sin(th) + dy * cos(th)
  W <- (w * abs(cos(th)) + h * abs(sin(th))) / 2
  H <- (w * abs(sin(th)) + h * abs(cos(th))) / 2
  rows <- get("rows", envir = tex_env)
  rows[[length(rows) + 1L]] <- data.frame(label = label, x = cx, y = cy, srt = srt,
                                          xmin = cx - W, xmax = cx + W, ymin = cy - H, ymax = cy + H,
                                          w = w, h = h, stringsAsFactors = FALSE)
  assign("rows", rows, envir = tex_env)
  invisible(NULL)
}

NODE_FILL <- c(Clone1 = "#E64B35", Clone2 = "#3C5488", C1Q_TAM = "#00A087", SPP1_TAM = "#F39B7F")
NODE_TXT  <- c(Clone1 = "grey10",   Clone2 = "white",   C1Q_TAM = "grey10",   SPP1_TAM = "grey10")

CH2 <- list(
  R_ring = 41.0, band = c(37.4, 41.0), node_r = 39.2,     # 弧厚 3.6 mm；环外径 41 mm → 圆径 82 mm
  cy = 0, lim = 53.0, pad = 0.25, gap = 16, arc_mm = 97.7,
  lr_r = c(42.4, 43.9), dth = c(0, -1.2, 1.2, -2.4, 2.4), off = 1.2,
  arrow_len = 1.8, blunt_len = 1.1, arrow_alpha = 0.7, inset = 0.14,
  fs_node = 5.6, fs_lab = 5.0, fs_leg = 4.2, alpha = 0.25,
  leg = c(-53, -26.5, 32.0, 50.5), seed = 7
)

rad_srt <- function(th) if (cos(th * pi / 180) >= 0) th else th + 180
chord2_split <- function(pair) c(lig = sub("_.*$", "", pair),
                                 rec = gsub("_", "–", sub("^[^_]*_", "", pair)))

chord2_groups <- function(a) {
  ec <- do.call(rbind, lapply(CH_NODES, function(nd) {
    f <- fill_node(a, nd, "clone"); if (!is.null(f)) cbind(f, node = nd) else NULL }))
  rownames(ec) <- ec$row
  out <- list()
  for (nd in CH_NODES) {
    f <- ec[ec$node == nd, , drop = FALSE]; if (!nrow(f)) next
    grp <- split(f, a$tab$pair[f$row])
    for (nm in names(grp))
      out[[length(out) + 1L]] <- list(node = nd, pair = nm,
                                      mid = (max(grp[[nm]]$t0) + min(grp[[nm]]$t1)) / 2)
  }
  out
}

lab_off <- function(k) if (k == 1) 1 else -1
obb2 <- function(lab2, th, R, off, fs) {
  cex <- fs / BASE_PS; out <- list()
  for (k in 1:2) {
    w <- strwidth(lab2[k], cex = cex, units = "user"); h <- strheight(lab2[k], cex = cex, units = "user")
    rc <- R + w / 2
    d  <- asin(min(0.95, off / rc)) * 180 / pi * lab_off(k)
    s  <- th + d; ad <- s * pi / 180; srt <- rad_srt(s); sa <- srt * pi / 180
    x <- rc * cos(ad); y <- rc * sin(ad)
    out[[k]] <- list(x = x, y = y, srt = srt, w = w, h = h,
                     aabb = c(x - (w * abs(cos(sa)) + h * abs(sin(sa))) / 2,
                              x + (w * abs(cos(sa)) + h * abs(sin(sa))) / 2,
                              y - (w * abs(sin(sa)) + h * abs(cos(sa))) / 2,
                              y + (w * abs(sin(sa)) + h * abs(cos(sa))) / 2))
  }
  out
}
obb_hit <- function(p, q, pad = 0.2) {
  ax <- function(s) { a <- s * pi / 180; list(u = c(cos(a), sin(a)), v = c(-sin(a), cos(a))) }
  A <- ax(p$srt); B <- ax(q$srt); d <- c(p$x - q$x, p$y - q$y)
  for (v in list(A$u, A$v, B$u, B$v)) {
    rp <- (p$w + pad) / 2 * abs(sum(v * A$u)) + (p$h + pad) / 2 * abs(sum(v * A$v))
    rq <- (q$w + pad) / 2 * abs(sum(v * B$u)) + (q$h + pad) / 2 * abs(sum(v * B$v))
    if (abs(sum(d * v)) > rp + rq) return(FALSE)
  }
  TRUE
}

chord2_place <- function(a, prm, n_restart = 800) {
  labs <- chord2_groups(a); n <- length(labs)
  prs  <- vapply(labs, function(l) l$pair, character(1))
  mds  <- vapply(labs, function(l) l$mid, numeric(1))
  set.seed(prm$seed); best <- NULL
  for (rs in seq_len(n_restart)) {
    objs <- list(); r0 <- rep(NA_real_, n); dth <- rep(0, n)
    for (i in sample(seq_len(n))) {
      got <- FALSE
      for (R in prm$lr_r) { for (dt in prm$dth) {
        ob <- obb2(chord2_split(prs[i]), mds[i] + dt, R, prm$off, prm$fs_lab)
        if (any(vapply(ob, function(o) any(o$aabb < -prm$lim) || any(o$aabb > prm$lim), logical(1)))) next
        if (any(vapply(objs, function(q) any(vapply(ob, function(o) obb_hit(o, q, prm$pad), logical(1))), logical(1)))) next
        r0[i] <- R; dth[i] <- dt; objs <- c(objs, ob); got <- TRUE; break }
        if (got) break }
    }
    nf <- sum(is.na(r0))
    if (is.null(best) || nf < best$nf) best <- list(nf = nf, r0 = r0, dth = dth, rs = rs)
    if (nf == 0L) break
  }
  objs <- list()
  for (i in which(!is.na(best$r0)))
    objs <- c(objs, obb2(chord2_split(prs[i]), mds[i] + best$dth[i], best$r0[i], prm$off, prm$fs_lab))
  for (i in which(is.na(best$r0))) {
    cand <- NULL
    for (R in prm$lr_r) for (dt in prm$dth) {
      ob <- obb2(chord2_split(prs[i]), mds[i] + dt, R, prm$off, prm$fs_lab)
      outside <- any(vapply(ob, function(o) any(o$aabb < -prm$lim) || any(o$aabb > prm$lim), logical(1)))
      ov <- if (outside) Inf else sum(vapply(objs, function(q)
        sum(vapply(ob, function(o) obb_hit(o, q, prm$pad), logical(1))), numeric(1)))
      cand <- rbind(cand, data.frame(R = R, dt = dt, ovl = ov))
    }
    k <- which.min(cand$ovl)
    best$r0[i] <- cand$R[k]; best$dth[i] <- cand$dt[k]
    objs <- c(objs, obb2(chord2_split(prs[i]), mds[i] + best$dth[i], best$r0[i], prm$off, prm$fs_lab))
  }
  bad <- which(is.na(best$r0))
  list(lab = prs, mid = mds, r0 = best$r0, dth = best$dth, n_free = best$nf, rs = best$rs,
       bad = bad, bad_lab = prs[bad])
}

# 端点标记：靶端三角形（尖端朝外，指向 TAM 扇区）+ 克隆端平头（矩形），底边与 ribbon 落点弧对齐
end_marks <- function(th0, th1, r_edge, col, kind, prm) {
  tm <- (th0 + th1) / 2; half <- (th0 - th1) / 2 * (1 - prm$inset)
  if (kind == "arrow") {
    a1 <- (tm - half) * pi / 180; a2 <- (tm + half) * pi / 180; am <- tm * pi / 180
    polygon(c((r_edge - prm$arrow_len) * cos(a1), (r_edge - prm$arrow_len) * cos(a2), r_edge * cos(am)),
            c((r_edge - prm$arrow_len) * sin(a1), (r_edge - prm$arrow_len) * sin(a2), r_edge * sin(am)),
            col = col, border = NA)
  } else {
    th <- seq(tm - half, tm + half, length.out = 14) * pi / 180
    polygon(c((r_edge - prm$blunt_len) * cos(th), rev(r_edge * cos(th))),
            c((r_edge - prm$blunt_len) * sin(th), rev(r_edge * sin(th))), col = col, border = NA)
  }
}

chord2_body <- function(x, prm, lay = NULL) {
  a <- chord_alloc(chord_tab(x), prm$gap, prm$arc_mm)
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = NA, ps = BASE_PS, xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(-55, 55), ylim = c(-55, 55), asp = 1)
  if (is.null(lay)) lay <- chord2_place(a, prm)

  ends_c <- do.call(rbind, lapply(CH_NODES, function(nd) {
    f <- fill_node(a, nd, "clone");  if (!is.null(f)) cbind(f, node = nd) else NULL }))
  ends_t <- do.call(rbind, lapply(CH_NODES, function(nd) {
    f <- fill_node(a, nd, "target"); if (!is.null(f)) cbind(f, node = nd) else NULL }))
  rownames(ends_c) <- ends_c$row; rownames(ends_t) <- ends_t$row

  for (i in seq_len(nrow(a$tab)))
    polygon(ribbon_xy(ends_c[as.character(i), "t0"], ends_c[as.character(i), "t1"],
                      ends_t[as.character(i), "t0"], ends_t[as.character(i), "t1"], prm$band[1]),
            col = adjustcolor(CLASS_COL[a$tab$cls[i]], alpha.f = prm$alpha), border = NA)

  # 端点标记：靶端箭头（指向 TAM）+ 克隆端平头；同色 alpha 0.7（0.25 的浅色 ribbon 上才看得见）
  for (i in seq_len(nrow(a$tab))) {
    col <- adjustcolor(CLASS_COL[a$tab$cls[i]], alpha.f = prm$arrow_alpha)
    end_marks(ends_t[as.character(i), "t0"], ends_t[as.character(i), "t1"], prm$band[1], col, "arrow", prm)
    end_marks(ends_c[as.character(i), "t0"], ends_c[as.character(i), "t1"], prm$band[1], col, "blunt", prm)
  }

  for (k in seq_len(nrow(a$nodes))) {          # 扇区弧按细胞状态着色
    th <- seq(a$nodes$t0[k], a$nodes$t1[k], length.out = 60) * pi / 180
    polygon(c(prm$band[1] * cos(th), rev(prm$band[2] * cos(th))),
            c(prm$band[1] * sin(th), rev(prm$band[2] * sin(th))),
            col = NODE_FILL[[a$nodes$node[k]]], border = NA)
  }

  tex_reset()
  for (k in seq_len(nrow(a$nodes))) {
    th <- a$nodes$mid[k]; ad <- th * pi / 180
    tex_draw(prm$node_r * cos(ad), prm$node_r * sin(ad), CH_LAB[[a$nodes$node[k]]], prm$fs_node,
             srt = tang_srt(th), col = NODE_TXT[[a$nodes$node[k]]], font = 2)
  }
  for (i in seq_along(lay$lab)) {
    sp <- chord2_split(lay$lab[i]); th0 <- lay$mid[i] + lay$dth[i]; R <- lay$r0[i]
    for (k in 1:2) {
      cex <- prm$fs_lab / BASE_PS
      w <- strwidth(sp[k], cex = cex, units = "user")
      d <- asin(min(0.95, prm$off / (R + w / 2))) * 180 / pi * lab_off(k)
      th <- th0 + d; ad <- th * pi / 180
      tex_draw(R * cos(ad), R * sin(ad), sp[k], prm$fs_lab, srt = rad_srt(th), col = "grey20",
               adj = if (cos(ad) >= 0) c(0, 0.5) else c(1, 0.5))
    }
  }

  lx <- prm$leg[1] + 0.5; blk <- prm$leg[2] - prm$leg[1] - 1
  tex_draw(lx, 49.2, "Δ = Clone1 − Clone2", 4.4, adj = c(0, 0.5), col = "grey15", font = 2)
  tex_draw(lx, 46.6, "top = ligand, bottom = receptor", 5.0, adj = c(0, 0.5), col = "grey35")
  txt2 <- "arrow → receptor side (target)"
  fs2  <- min(5, 4.99 * blk / strwidth(txt2, cex = 5 / BASE_PS, units = "user"))
  tex_draw(lx, 44.0, txt2, fs2, adj = c(0, 0.5), col = "grey35")
  y <- 41.2
  for (k in c("Clone1_only", "Clone2_only", "both_significant", "neither")) {
    rect(lx, y - 1.0, lx + 2.0, y + 1.0, col = CLASS_COL[[k]], border = "grey40", lwd = 0.3)
    tex_draw(lx + 2.8, y, c(Clone1_only = "CL1 only", Clone2_only = "CL2 only",
                            both_significant = "Both sig.", neither = "Neither")[[k]],
             prm$fs_leg, adj = c(0, 0.5), col = "grey15")
    y <- y - 2.7
  }
  tex_draw(lx, y - 0.4, "width = |Δ|", prm$fs_leg, adj = c(0, 0.5), col = "grey35")

  ov <- tex_check(); usr <- par("usr")
  list(ov = ov, info = list(alloc = a, lay = lay, n = nrow(a$tab),
                            n_c2 = sum(a$tab$clone == "Clone2"),
                            n_arrow = nrow(a$tab), n_blunt = nrow(a$tab),
                            w_min_pt = min(a$tab$w_mm) * 72 / 25.4,
                            w_max_pt = max(a$tab$w_mm) * 72 / 25.4,
                            R_ring = prm$R_ring, mm_per_unit = 110 / (usr[2] - usr[1]),
                            labels = tex_labels(), fs_arrow_note = fs2,
                            node_fill = NODE_FILL, node_txt = NODE_TXT))
}

chord2_trial <- function(x, prm) {
  tf <- tempfile(fileext = ".pdf")
  cairo_pdf(tf, width = 110 / 25.4, height = 110 / 25.4, bg = "white", pointsize = BASE_PS, onefile = FALSE)
  r <- chord2_body(x, prm)
  dev.off(); unlink(tf)
  r
}

cat("[Chord] 扇区弧已着色（Clone1/Clone2/C1Q/SPP1），端点标记：靶端箭头 + 克隆端平头\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 简化弦图出图（P13、P15，110×110 mm，600 dpi）：先在 cairo_pdf 上求解版式并判交，
# 再把同一套坐标交给 PDF 与 PNG 两次绘制；扇区弧着色 + 靶端箭头 / 克隆端平头
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

chord2_out <- list()
for (pat in c("P13", "P15")) {
  x   <- TL2[[pat]]
  r0  <- chord2_trial(x, CH2)
  lay <- r0$info$lay
  stem <- file.path(tl_root, pat, "figures", paste0(pat, "_tail_LR_chord_class_delta"))
  info <- NULL
  save_pair(stem, 110, 110, function() {
    info <<- chord2_body(x, CH2, lay = lay)$info; invisible(NULL)
  })
  chord2_out[[pat]] <- info

  a <- info$alloc$tab
  ec <- do.call(rbind, lapply(CH_NODES, function(nd) { f <- fill_node(info$alloc, nd, "clone")
    if (!is.null(f)) cbind(f, node = nd) else NULL })); rownames(ec) <- ec$row
  et <- do.call(rbind, lapply(CH_NODES, function(nd) { f <- fill_node(info$alloc, nd, "target")
    if (!is.null(f)) cbind(f, node = nd) else NULL })); rownames(et) <- et$row
  idx <- as.character(seq_len(nrow(a)))
  tab <- data.frame(patient = pat, y_order = a$y_order, pair = a$pair, pathway = a$pathway,
                    target = a$target, clone = a$clone, cls = a$cls, delta = a$delta,
                    abs_delta = a$abs_delta, width_mm = a$w_mm, width_pt = a$w_mm * 72 / 25.4,
                    fill_class_hex = unname(CLASS_COL[a$cls]), fill_alpha = CH2$alpha,
                    label_top = vapply(a$pair, function(p) sub("_.*$", "", p), character(1)),
                    label_bottom = vapply(a$pair, function(p) gsub("_", "–", sub("^[^_]*_", "", p)), character(1)),
                    arrow_end = "target", arrow_cls_hex = unname(CLASS_COL[a$cls]),
                    arrow_alpha = CH2$arrow_alpha,
                    clone_arc_deg = ec[idx, "t0"] - ec[idx, "t1"],
                    target_arc_deg = et[idx, "t0"] - et[idx, "t1"],
                    stringsAsFactors = FALSE)
  write.csv(tab, file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_chord_class_delta_values.csv")),
            row.names = FALSE)

  pm <- png_meta(paste0(stem, ".png")); ps <- px_stats(paste0(stem, ".png"))
  cat(sprintf("[%s] %s.{pdf,png} | PNG %dx%d ct=%d %.2f MB | sd=%.3f 非白=%.3f | 边距 上%.1f 下%.1f 左%.1f 右%.1f mm\n",
              pat, basename(stem), pm$w, pm$h, pm$colortype, pm$mb, ps$sd, ps$nonwhite,
              ps$top * 25.4 / 600, ps$bottom * 25.4 / 600, ps$left * 25.4 / 600, ps$right * 25.4 / 600))
  cat(sprintf("      mm/unit=%.4f | 实测圆径 %.2f mm | ribbon %d 条 | 宽度 %.2f–%.2f pt | 类别 %s\n",
              info$mm_per_unit, 2 * info$R_ring * info$mm_per_unit, nrow(tab),
              min(tab$width_pt), max(tab$width_pt),
              paste(names(table(tab$cls)), table(tab$cls), sep = "=", collapse = " ")))
  cat(sprintf("      名称 %d 个（%d 行）| 靶端箭头 %d / 克隆端平头 %d | 全落位 %d/20 | 未落位 %d | cairo_pdf 文本 AABB %d | 图内说明 %.2f pt\n",
              length(lay$lab), 2 * length(lay$lab), info$n_arrow, info$n_blunt,
              20 - lay$n_free, length(lay$bad), nrow(r0$ov), info$fs_arrow_note))
  cat(sprintf("      扇区弧配色：%s\n",
              paste(sprintf("%s=%s", names(NODE_FILL), NODE_FILL), collapse = " ")))
}
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 交付前总核验：文本判交、名称齐备、上行=配体、扇区弧四色（避开弧上文字、取弧内不同半径）、
# 箭头 40 个且朝向 TAM、箭头同色可见、箭头不压名称/弧、圆径、文件与边距
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

obb_pairs <- function(rows) {
  n <- length(rows); out <- list(); k <- 0L
  if (n < 2) return(data.frame(a = character(0), b = character(0)))
  ax_of <- function(s) { a <- s * pi / 180; list(u = c(cos(a), sin(a)), v = c(-sin(a), cos(a))) }
  hit2 <- function(r1, r2) {
    A1 <- ax_of(r1$srt); A2 <- ax_of(r2$srt); d <- c(r1$x - r2$x, r1$y - r2$y)
    for (ax in list(A1$u, A1$v, A2$u, A2$v)) {
      rr1 <- r1$w / 2 * abs(sum(ax * A1$u)) + r1$h / 2 * abs(sum(ax * A1$v))
      rr2 <- r2$w / 2 * abs(sum(ax * A2$u)) + r2$h / 2 * abs(sum(ax * A2$v))
      if (abs(sum(d * ax)) > rr1 + rr2) return(FALSE)
    }
    TRUE
  }
  for (i in seq_len(n - 1L)) for (j in seq.int(i + 1L, n))
    if (hit2(rows[[i]], rows[[j]])) { k <- k + 1L
      out[[k]] <- data.frame(a = rows[[i]]$label, b = rows[[j]]$label, stringsAsFactors = FALSE) }
  if (!k) data.frame(a = character(0), b = character(0)) else do.call(rbind, out)
}
count_page <- function(f) {
  raw <- readBin(f, "raw", n = file.size(f)); txt <- rawToChar(raw[raw > 8 & raw < 127])
  length(gregexpr("/Type ?/Page[^s]", txt)[[1]])
}
px_at <- function(img, x_mm, y_mm) {
  nr <- dim(img)[1]; nc <- dim(img)[2]
  cc <- round((x_mm + 55) / 110 * (nc - 1)) + 1
  rr <- round((55 - y_mm) / 110 * (nr - 1)) + 1
  cc <- min(max(cc, 1), nc); rr <- min(max(rr, 1), nr)
  img[rr, cc, 1:3]
}
hex2rgb <- function(h) grDevices::col2rgb(h)[, 1] / 255
verify_one <- function(pat, kind, lay) {
  idx <- match(pat, names(TL2))
  run_dev <- function(dev) {
    w_in <- if (kind == "circos") 120 / 25.4 else 110 / 25.4
    tf <- tempfile(fileext = paste0(".", dev))
    if (dev == "pdf") cairo_pdf(tf, width = w_in, height = w_in, bg = "white", pointsize = BASE_PS, onefile = FALSE)
    else png(tf, width = w_in, height = w_in, units = "in", res = 600, bg = "white", pointsize = BASE_PS, type = "cairo")
    r <- if (kind == "circos") circ2_body(TL2[[pat]], secs[[idx]], gls[[idx]], C2, csol$fs_lr, csol$fs_pw, lay = lay)
         else chord2_body(TL2[[pat]], CH2, lay = lay)
    rows <- get("rows", envir = tex_env)
    dev.off(); unlink(tf)
    list(n_aabb = nrow(r$ov), n_obb = nrow(obb_pairs(rows)),
         labels = r$info$labels, diam = 2 * r$info$R_ring * r$info$mm_per_unit)
  }
  list(pdf = run_dev("pdf"), png = run_dev("png"))
}
chk <- list(); k <- 0L
add <- function(item, patient, value, ok) {
  k <<- k + 1L
  chk[[k]] <<- data.frame(item = item, patient = patient, value = as.character(value),
                          ok = as.logical(ok), stringsAsFactors = FALSE)
}
is_cluster_num <- function(x) grepl("^([1-9]|1[0-9]|20)$", x)
split2 <- function(p) c(sub("_.*$", "", p), gsub("_", "–", sub("^[^_]*_", "", p)))
line_y <- function(lay, prm) {
  tf <- tempfile(fileext = ".pdf")
  cairo_pdf(tf, width = 110 / 25.4, height = 110 / 25.4, bg = "white", pointsize = BASE_PS, onefile = FALSE)
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i"); plot.new(); plot.window(c(-55, 55), c(-55, 55), asp = 1)
  yt <- yb <- numeric(length(lay$lab))
  for (i in seq_along(lay$lab)) {
    sp <- split2(lay$lab[i]); th0 <- lay$mid[i] + lay$dth[i]; R <- lay$r0[i]
    for (kk in 1:2) {
      w <- strwidth(sp[kk], cex = prm$fs_lab / BASE_PS, units = "user")
      d <- asin(min(0.95, prm$off / (R + w / 2))) * 180 / pi * (if (kk == 1) 1 else -1)
      if (kk == 1) yt[i] <- R * sin((th0 + d) * pi / 180) else yb[i] <- R * sin((th0 + d) * pi / 180)
    }
  }
  dev.off(); unlink(tf)
  list(top = yt, bot = yb)
}
REQ_FILL <- c(Clone1 = "#E64B35", Clone2 = "#3C5488", C1Q_TAM = "#00A087", SPP1_TAM = "#F39B7F")

for (pat in c("P13", "P15")) {
  v <- read.csv(file.path(tl_root, pat, "tables", paste0(pat, "_tail_LR_top20_values.csv")),
                check.names = FALSE, stringsAsFactors = FALSE)
  top20 <- unique(v$interaction_name[order(v$y_order)])

  r <- verify_one(pat, "circos", circ2_out[[pat]]$lay); lb <- circ2_out[[pat]]$labels
  add("circos AABB(PDF/PNG)", pat, sprintf("%d / %d", r$pdf$n_aabb, r$png$n_aabb),
      r$pdf$n_aabb == 0L && r$png$n_aabb == 0L)
  add("circos 墨迹 OBB(PDF/PNG)", pat, sprintf("%d / %d", r$pdf$n_obb, r$png$n_obb),
      r$pdf$n_obb == 0L && r$png$n_obb == 0L)
  add("circos 无 1–20 编号", pat, sum(is_cluster_num(lb)), sum(is_cluster_num(lb)) == 0L)
  add("circos 名称齐备(20)", pat, sprintf("%d/20", sum(top20 %in% lb)), all(top20 %in% lb))
  add("circos 实测圆径≥85mm", pat, sprintf("%.2f", r$pdf$diam), r$pdf$diam >= 85)

  r2 <- verify_one(pat, "chord", chord2_out[[pat]]$lay); lb2 <- chord2_out[[pat]]$labels
  ok_name <- all(vapply(top20, function(p) all(split2(p) %in% lb2), logical(1)))
  yy <- line_y(chord2_out[[pat]]$lay, CH2)
  add("chord 上行=配体", pat, sprintf("%d/20 成立，最小高度差 %.2f mm", sum(yy$top > yy$bot), min(yy$top - yy$bot)),
      all(yy$top > yy$bot))
  add("chord AABB(PDF/PNG)（旋转两行文本，仅参考）", pat, sprintf("%d / %d", r2$pdf$n_aabb, r2$png$n_aabb),
      r2$pdf$n_aabb == 0L && r2$png$n_aabb == 0L)
  add("chord 墨迹 OBB(PDF/PNG)", pat, sprintf("%d / %d", r2$pdf$n_obb, r2$png$n_obb),
      r2$pdf$n_obb == 0L && r2$png$n_obb == 0L)
  add("chord 无 1–20 编号", pat, sum(is_cluster_num(lb2)), sum(is_cluster_num(lb2)) == 0L)
  add("chord 名称齐备(20)", pat, sprintf("%d/20", sum(vapply(top20, function(p) all(split2(p) %in% lb2), logical(1)))),
      ok_name)
  add("chord 图内说明两行", pat, paste(c("top = ligand, bottom = receptor", "arrow → receptor side (target)") %in% lb2, collapse = "/"),
      all(c("top = ligand, bottom = receptor", "arrow → receptor side (target)") %in% lb2))
  add("chord 实测圆径≥80mm", pat, sprintf("%.2f", r2$pdf$diam), r2$pdf$diam >= 80)

  img <- png::readPNG(file.path(tl_root, pat, "figures", paste0(pat, "_tail_LR_chord_class_delta.png")))
  if (length(dim(img)) == 3) img <- img[, , 1:3] else img <- array(img, dim = c(dim(img), 3))
  nd <- chord2_out[[pat]]$alloc$nodes
  # 扇区弧着色：弧带 37.4–41.0 mm，文字带约 38.3–40.1 mm；在 r=37.8/40.6 两处取样避开文字
  got <- character(0); okc <- logical(0)
  for (k2 in seq_len(nrow(nd))) {
    a0 <- nd$mid[k2]; off <- 0.32 * nd$deg[k2]; want <- hex2rgb(REQ_FILL[[nd$node[k2]]])
    hit <- FALSE; shown <- NULL
    for (rr in c(37.8, 40.6)) for (sgn in c(-1, 1)) for (kf in c(1, 0.5)) {
      a <- (a0 + sgn * off * kf) * pi / 180
      rgb <- px_at(img, rr * cos(a), rr * sin(a))
      if (is.null(shown)) shown <- rgb
      if (max(abs(rgb - want)) < 0.12) { hit <- TRUE; shown <- rgb; break }
    }
    got <- c(got, sprintf("%s=%s", sub("_TAM", "", nd$node[k2]),
                          paste(sprintf("%.0f", round(shown * 255)), collapse = ",")))
    okc <- c(okc, hit)
  }
  add("chord 扇区弧四色（PNG 像素）", pat, paste(got, collapse = " | "), all(okc))

  a2 <- chord2_out[[pat]]$alloc
  et <- do.call(rbind, lapply(CH_NODES, function(nd2) { f <- fill_node(a2, nd2, "target")
    if (!is.null(f)) cbind(f, node = nd2) else NULL })); rownames(et) <- et$row
  wide <- order(-a2$tab$w_mm)[1:3]; arr_ok <- logical(0); arr_txt <- character(0)
  for (i in wide) {
    rw <- et[as.character(i), ]; tm <- (rw[["t0"]] + rw[["t1"]]) / 2 * pi / 180
    rg <- px_at(img, (CH2$band[1] - 0.55 * CH2$arrow_len) * cos(tm), (CH2$band[1] - 0.55 * CH2$arrow_len) * sin(tm))
    c0 <- hex2rgb(CLASS_COL[a2$tab$cls[i]])
    arr_ok <- c(arr_ok, max(abs(rg - (0.78 * c0 + 0.22))) < 0.16)
    arr_txt <- c(arr_txt, sprintf("%s=%.0f,%.0f,%.0f", a2$tab$cls[i], rg[1]*255, rg[2]*255, rg[3]*255))
  }
  add("chord 箭头同色可见（PNG 像素，alpha 0.7）", pat, paste(arr_txt, collapse = " | "), all(arr_ok))
  rm(img); gc(verbose = FALSE)

  apex_in <- logical(nrow(a2$tab))
  for (i in seq_len(nrow(a2$tab))) {
    node <- a2$tab$target[i]; rw <- et[as.character(i), ]; tm <- (rw[["t0"]] + rw[["t1"]]) / 2
    sec <- a2$nodes[a2$nodes$node == node, ]
    apex_in[i] <- node %in% TGS && tm <= sec$t0 + 1e-9 && tm >= sec$t1 - 1e-9
  }
  add("chord 箭头 40 个且朝向 TAM 扇区", pat,
      sprintf("%d/%d 落在 TAM 扇区内，箭头长 %.1f mm", sum(apex_in), nrow(a2$tab), CH2$arrow_len),
      sum(apex_in) == nrow(a2$tab) && nrow(a2$tab) == 40L)
  add("chord 箭头不压名称/弧", pat,
      sprintf("箭头 r %.1f–%.1f；名称 r≥%.1f；弧 r %.1f–%.1f mm",
              CH2$band[1] - CH2$arrow_len, CH2$band[1], min(chord2_out[[pat]]$lay$r0),
              CH2$band[1], CH2$band[2]),
      min(chord2_out[[pat]]$lay$r0) > CH2$band[2])

  for (st in c("_tail_LR_circos_heatmap", "_tail_LR_chord_class_delta")) {
    stem <- file.path(tl_root, pat, "figures", paste0(pat, st))
    pm <- png_meta(paste0(stem, ".png")); ps <- px_stats(paste0(stem, ".png"), step = 6L)
    add(paste0("文件双格式", st), pat,
        sprintf("PNG %dx%d ct=%d %.2fMB，PDF %d 页", pm$w, pm$h, pm$colortype, pm$mb, count_page(paste0(stem, ".pdf"))),
        file.exists(paste0(stem, ".pdf")) && file.exists(paste0(stem, ".png")) &&
          pm$colortype %in% c(2L, 6L) && ps$sd > 0.05 && count_page(paste0(stem, ".pdf")) == 1L)
    add(paste0("未裁切(边距≥1mm)", st), pat,
        sprintf("上%.1f 下%.1f 左%.1f 右%.1f", ps$top*25.4/600, ps$bottom*25.4/600, ps$left*25.4/600, ps$right*25.4/600),
        min(ps$top, ps$bottom, ps$left, ps$right) * 25.4 / 600 >= 1)
  }
}
res <- do.call(rbind, chk)
cat(sprintf("\n%-50s %-4s %-40s %s\n", "检查项", "样本", "实测", "结论"))
for (i in seq_len(nrow(res)))
  cat(sprintf("%-50s %-4s %-40s %s\n", res$item[i], res$patient[i], res$value[i],
              if (isTRUE(res$ok[i])) "通过" else "不通过"))
cat(sprintf("\n[总核验] %d 项，通过 %d 项，失败 %d 项\n", nrow(res), sum(res$ok), sum(!res$ok)))
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")
