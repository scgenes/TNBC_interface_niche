# =============================================================================
# TNBC interface-niche code repository
# Original source: project_spatial/TNBC/TNBC_analysis/SpatialInferCNV/119B_clone_analysis/run_clone_analysis.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

#!/usr/bin/env Rscript
# ============================================================================
# 119B 单样品克隆分析
#   问题1: 这个样品能划分几个 clone?
#   问题2: 不同 clone 之间的 CNV 变化是什么?
#
# ---------------------------------------------------------------------------
# 数据与口径
# ---------------------------------------------------------------------------
# 输入(全部只读, 不修改任何已有文件):
#   samples/119B/run.final.infercnv_obj                      CNV 矩阵 (denoise 后终态)
#   samples/119B/HMM_CNV_predictions.HMMi6...pred_cnv_genes.dat  HMM 基因级状态
#   input/gene_order_filtered.tsv                            gene -> chr/pos
#   RCID/spatial_noHE/spatial_data_cache.rds                 空间坐标(仅为画克隆空间图)
#
# 关键: expr.data 的基线(二倍体) = 参考组中位数 = 1.0036, 不是 1.0。
#       本样品参考组(TME 148 细胞)在**臂级**上完全中性, 但基因级 HMM 同样给出
#       大量 state 2/4 —— 说明 HMM 的 2/4 被高频焦点噪声主导。
#       因此: 基因级 state 只用来做"事件描述", 克隆划分与差异检验用**臂级平均**。
#
# HMM state -> CNV 水平(本样品实测标定, 见 hmm_state_calibration.csv):
#   state 2 = loss   (中位 0.868)
#   state 3 = 中性    (中位 1.0036)
#   state 4 = gain   (中位 1.136)
#   state 5 = gain+  (中位 1.303)
#   state 6 = gain++ (中位 1.579)
#
# ---------------------------------------------------------------------------
# 方法
# ---------------------------------------------------------------------------
# 1) 参考组体检: 确认基线、确认臂级中性
# 2) cell x 42 臂 的 CNV 矩阵 (expr.data 基线校正后按臂平均)
# 3) 克隆数 k 的选择(不硬切):
#      - silhouette(k=2..10) 在两种分辨率下做: 臂级 / 基因级 PCA
#      - 100 次 80% 重抽样的簇稳定性 (Jaccard, 按理论上限归一化)
#      - 与 infercnv 自带 hclust (tumor_subclusters$hc$all_observations) 的一致性
#    k* = silhouette 最大者; 并检查 k*+1 是否带来新的显著 CNV 差异
# 4) 克隆间 CNV 差异: 逐臂 Wilcoxon + BH-FDR; 多克隆用 Kruskal-Wallis
# 5) 空间验证: 克隆是否成片(同克隆近邻比例 vs 打乱), 这是"真亚克隆"的关键证据
# 6) 出图 + 导出 Newick(供 phylogram / ggtree 用)
#
# 输出: 本目录 (119B_clone_analysis/), 共 10 个表 + 5 张图 + 2 个 Newick + 1 个日志
# 回滚: 直接删除本目录即可; 输入文件全程只读
# 用法: Rscript run_clone_analysis.R
# ============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(ggplot2); library(FNN); library(RColorBrewer)
})

SAMPLE   <- Sys.getenv("SAMPLE", "119B")
BASE_DIR <- Sys.getenv("SICNV_BASE", "/path/to/data/project_spatial_data/TNBC/SpatialInferCNV")
SAMP_DIR <- file.path(BASE_DIR, "samples", SAMPLE)
OUT_DIR  <- Sys.getenv("CLONE_OUT", file.path(BASE_DIR, paste0(SAMPLE, "_clone_analysis")))
XY_CACHE <- Sys.getenv("XY_CACHE",
                       "/path/to/data/project_spatial_data/TNBC/RCID/spatial_noHE/spatial_data_cache.rds")
STEPS    <- 2:10
N_BOOT   <- as.integer(Sys.getenv("N_BOOT", "100"))
SEED     <- 20261010L

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT_DIR, "run.log")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG)   # 覆盖上次日志

# 克隆配色: Okabe-Ito (CVD 安全); 顺序固定, 不循环
CLONE_COL <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
               "#56B4E9", "#F0E442", "#999999", "#000000")

say("================ %s 克隆分析 ================", SAMPLE)
say("输出目录: %s", OUT_DIR)

# ==================================================================
# 1. 读对象(只读) + 参考组体检
# ==================================================================
obj_f <- file.path(SAMP_DIR, "run.final.infercnv_obj")
if (!file.exists(obj_f)) stop("找不到 ", obj_f)
obj <- readRDS(obj_f)
ex  <- obj@expr.data
ref <- colnames(ex)[unlist(obj@reference_grouped_cell_indices)]
obs <- colnames(ex)[unlist(obj@observation_grouped_cell_indices)]
BASE_VAL <- median(ex[, ref])                 # 二倍体基线 = 参考组中位数
say("\n[1] 对象: %s", basename(obj_f))
say("    expr.data %d 基因 x %d 细胞 (观测 %d / 参考 %d)",
    nrow(ex), ncol(ex), length(obs), length(ref))
say("    二倍体基线 (参考组中位数) = %.4f", BASE_VAL)

# ==================================================================
# 2. HMM state 标定(用数据本身, 不假设 state 编号含义)
# ==================================================================
pg_fs <- list.files(SAMP_DIR, pattern = "pred_cnv_genes\\.dat$", full.names = TRUE)
# 优先 Pnorm_0.5 版(= BayesNet 精修后的终态调用, 与 run.final 对象同源);
# 17_ 版是 HMM 早阶段产物, state 3 未被发射, 直接用它会把 state 3 误报成不存在。
pg_f <- pg_fs[grepl("Pnorm_0\\.5|HMM_CNV_predictions", pg_fs)][1]
if (is.na(pg_f)) pg_f <- pg_fs[1]
cal <- NULL
if (!is.na(pg_f)) {
  say("\n[2] HMM state 标定来源: %s (mtime %s)", basename(pg_f),
      format(file.info(pg_f)$mtime))
  pg <- read.delim(pg_f, stringsAsFactors = FALSE)
  i <- cbind(match(pg$gene, rownames(ex)), match(pg$cell_group_name, colnames(ex)))
  k <- !is.na(i[, 1]) & !is.na(i[, 2])
  v <- ex[i[k, , drop = FALSE]]; st <- pg$state[k]
  cal <- do.call(rbind, lapply(sort(unique(st)), function(s) {
    x <- v[st == s]
    data.frame(state = s, n_gene_calls = length(x), median = round(median(x), 4),
               q05 = round(quantile(x, .05), 4), q95 = round(quantile(x, .95), 4))
  }))
  write.csv(cal, file.path(OUT_DIR, "hmm_state_calibration.csv"), row.names = FALSE)
  say("    n=%.0f%% 的 gene-call 落在 state 2/4, 故 state 仅用于事件描述, 不用于划分克隆:",
      100 * sum(cal$n_gene_calls[cal$state %in% c(2, 4)]) / sum(cal$n_gene_calls))
  for (j in seq_len(nrow(cal)))
    say("      state %d: n=%9d  中位 %.4f", cal$state[j], cal$n_gene_calls[j], cal$median[j])
}

# ==================================================================
# 3. 基因 -> 染色体臂, 建 cell x arm 矩阵
# ==================================================================
CENTRO <- c(chr1=125, chr2=93.3, chr3=91, chr4=50.4, chr5=48.4, chr6=60.3, chr7=60.1,
            chr8=45.6, chr9=49, chr10=40.2, chr11=53.7, chr12=35.8, chr13=17.9, chr14=17,
            chr15=19, chr16=36.8, chr17=25.1, chr18=18.5, chr19=26.1, chr20=28.1,
            chr21=12, chr22=15, chrX=61, chrY=10.4)
go <- read.delim(file.path(BASE_DIR, "input/gene_order_filtered.tsv"), header = FALSE,
                 col.names = c("gene", "chr", "start", "end"))
go <- go[!duplicated(go$gene), ]
go$arm <- ifelse(is.na(CENTRO[go$chr]), NA_character_,
                 paste0(go$chr, ifelse(go$start / 1e6 < CENTRO[go$chr], "p", "q")))
m  <- match(rownames(ex), go$gene)
ok <- !is.na(m) & !is.na(go$arm[m])
gx <- ex[ok, , drop = FALSE]                  # 基因 x cell (只留有臂归属的基因)
arm_of <- go$arm[m[ok]]
ARMS <- sort(unique(arm_of))
# 按染色体号 + p/q 排 (必须同时去掉 "chr" 前缀, 否则 as.numeric 全 NA -> 静默退化为字典序)
# chrX/chrY 无法数值化, suppressWarnings 后由 NA 自然排到末尾(惯例位置)
ARMS <- ARMS[suppressWarnings(order(as.numeric(sub("^chr", "", sub("[pq]$", "", ARMS))),
                                    sub("^chr[0-9XY]+", "", ARMS)))]
say("\n[3] 基因->臂: %d/%d 基因有臂归属, 覆盖 %d 条臂", sum(ok), nrow(ex), length(ARMS))

arm_mean <- function(cells) {
  M <- vapply(ARMS, function(a) colMeans(gx[arm_of == a, cells, drop = FALSE]),
              numeric(length(cells)))
  M <- matrix(M, nrow = length(cells), dimnames = list(cells, ARMS)); M
}
M    <- arm_mean(obs)                          # 观测细胞的臂级矩阵
Mref <- arm_mean(ref)                          # 参考组同口径(体检)
say("    参考组臂级: |dev|>0.02 的臂 %d 条 (=0 说明参考组臂级中性, 基线可信)",
    sum(abs(colMeans(Mref) - BASE_VAL) > 0.02))

# ==================================================================
# 4. 克隆划分: k 的选择
# ==================================================================
D_arm <- dist(M)                               # 臂级欧氏距离
hc_a  <- hclust(D_arm, method = "ward.D2")

Xg    <- t(gx[, obs, drop = FALSE])            # 409 x n_gene  (gx 已按 ok 过滤)
pca   <- prcomp(Xg, center = TRUE, scale. = FALSE, rank. = 30)
cumv  <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
NPC   <- max(3, min(which(cumv >= 0.80)[1], 15))
D_gn  <- dist(pca$x[, 1:NPC, drop = FALSE])
hc_g  <- hclust(D_gn, method = "ward.D2")
say("\n[4] 基因级 PCA: 保留 %d 个 PC (累计方差 %.1f%%); PC1=%.1f%%",
    NPC, 100 * cumv[NPC], 100 * cumv[1])

sil <- function(cl, d) mean(cluster::silhouette(cl, d)[, 3])
ks  <- data.frame(k = STEPS,
                  sil_arm  = NA_real_, sil_gene = NA_real_,
                  n_min = NA_integer_, n_max = NA_integer_)
for (i in seq_along(STEPS)) {
  ka <- cutree(hc_a, k = STEPS[i]); kg <- cutree(hc_g, k = STEPS[i])
  ks$sil_arm[i]  <- sil(ka, D_arm)
  ks$sil_gene[i] <- sil(kg, D_gn)
  ks$n_min[i] <- min(table(kg)); ks$n_max[i] <- max(table(kg))
}
set.seed(SEED)
# 稳定性: 80% 重抽样 100 次; Jaccard 用"理论上限"归一化
# (子样本只含 80% 细胞, 完美重现也只能拿到 <=0.8 的 Jaccard, 不归一化会系统性低估)
Zg <- pca$x[, 1:NPC, drop = FALSE]
boot <- function(Z, k) {
  n <- nrow(Z)
  ref_cl <- cutree(hclust(dist(Z), method = "ward.D2"), k = k)
  out <- numeric(0); ceil <- numeric(0)
  for (b in seq_len(N_BOOT)) {
    ii <- sample(n, floor(0.8 * n))
    c2 <- cutree(hclust(dist(Z[ii, , drop = FALSE]), method = "ward.D2"), k = k)
    for (g in unique(ref_cl)) {
      A <- which(ref_cl == g)                 # 全集中的索引
      inb <- A[A %in% ii]                     # A 在子样本中的部分
      B <- which(c2 == g)                     # 子样本内索引(相对 ii)
      if (length(B) < 5 || length(inb) < 5) next
      jac <- sapply(unique(c2), function(gg) {
        Bg <- ii[which(c2 == gg)]
        length(intersect(inb, Bg)) / length(union(inb, Bg))
      })
      out  <- c(out, max(jac))
      ceil <- c(ceil, length(inb) / length(A))
    }
  }
  c(median_norm = median(out / ceil, na.rm = TRUE),
    q25_norm = quantile(out / ceil, .25, na.rm = TRUE, names = FALSE),
    n_eval = length(out))
}
ks$stability_gene <- NA_real_
for (k in 2:4) {
  b <- boot(Zg, k)
  ks$stability_gene[ks$k == k] <- round(b[["median_norm"]], 3)
  say("    [boot] k=%d 稳定性(上限归一) 中位 %.3f, q25 %.3f (n=%d)", k,
      b[["median_norm"]], b[["q25_norm"]], b[["n_eval"]])
}
write.csv(data.frame(k = 2:4, median_norm = ks$stability_gene[ks$k %in% 2:4]),
          file.path(OUT_DIR, "boot_stability.csv"), row.names = FALSE)
write.csv(ks, file.path(OUT_DIR, "clone_k_selection.csv"), row.names = FALSE)

say("\n    k  silhouette(臂级)  silhouette(基因级)  稳定性(上限归一)  簇大小")
for (i in seq_len(nrow(ks)))
  say("    %2d   %.3f              %.3f               %s          %d..%d",
      ks$k[i], ks$sil_arm[i], ks$sil_gene[i],
      ifelse(is.na(ks$stability_gene[i]), "  -  ", sprintf("%.3f", ks$stability_gene[i])),
      ks$n_min[i], ks$n_max[i])

K <- ks$k[which.max(ks$sil_gene)]              # 主判据: 基因级 silhouette 最大
say("\n    >>> 选定 k* = %d  (基因级 silhouette %.3f, 臂级 %.3f)",
    K, ks$sil_gene[ks$k == K], ks$sil_arm[ks$k == K])

cl_g <- cutree(hc_g, k = K); cl_a <- cutree(hc_a, k = K)
agree <- mean(cl_g == cl_a)
say("    两种分辨率划分一致率: %.1f%%", 100 * agree)
write.csv(data.frame(k = K, agreement_arm_vs_gene = round(agree, 4)),
          file.path(OUT_DIR, "resolution_agreement.csv"), row.names = FALSE)

# 与 infercnv 自带的 hclust 对照
ts <- tryCatch(obj@tumor_subclusters, error = function(e) NULL)
if (!is.null(ts) && !is.null(ts$hc[["all_observations"]])) {
  hc_i <- ts$hc[["all_observations"]]
  nl <- length(unique(ts$subclusters[["all_observations"]]))
  say("    infercnv 自带: 显著性检验给出 %d 个 subcluster; 其 hclust 树顶切 2 份大小 = %s",
      nl, paste(sort(table(cutree(hc_i, 2)), decreasing = TRUE), collapse = " + "))
  # 顺序对齐后算一致率(取两种标签对应关系中最优)
  ci <- cutree(hc_i, k = K)
  if (length(ci) == length(cl_g)) {
    tab <- table(cl_g, ci)
    ag  <- sum(apply(tab, 1, max)) / sum(tab)   # 标签最优对应下的上限一致率
    say("    infercnv 树切 k=%d 与本次划分的最优标签一致率: %.1f%%", K, 100 * ag)
    write.csv(as.data.frame.matrix(tab), file.path(OUT_DIR, "infercnv_hc_crosstab.csv"))
  }
}

# ==================================================================
# 5. 克隆成员表 / 汇总 / 差异臂
# ==================================================================
CL <- factor(cl_g, levels = 1:K)
# 按大小重命名: C1 = 最大簇
rel <- order(-as.integer(table(CL)))
newlab <- integer(K); newlab[rel] <- seq_len(K)
CL <- factor(newlab[as.integer(CL)], levels = seq_len(K))
names(CL) <- obs
say("\n[5] 克隆大小: %s", paste(sprintf("C%d n=%d (%.1f%%)", seq_len(K),
    as.integer(table(CL)), 100 * as.integer(table(CL)) / length(obs)), collapse = " | "))

burden <- colMeans(abs(gx[, obs, drop = FALSE] - BASE_VAL) > 0.10)
say("    CNV burden (|expr-基线|>0.10 的基因占比) 中位: %s",
    paste(sprintf("C%d %.3f", seq_len(K), tapply(burden, CL, median)), collapse = " | "))

# 克隆 x 臂 平均谱
CM <- do.call(rbind, lapply(seq_len(K), function(g) colMeans(M[CL == g, , drop = FALSE])))
rownames(CM) <- paste0("C", seq_len(K))
write.csv(round(CM, 4), file.path(OUT_DIR, "clone_arm_mean.csv"))
write.csv(round(CM - BASE_VAL, 4), file.path(OUT_DIR, "clone_arm_deviation.csv"))

# ---- 克隆 QC: 排除"低测序深度/低质量 spot"冒充克隆 ----
# 这是本分析最关键的对照: 若某簇 UMI 远低于其它簇, 它的谱会被噪声抹平,
# 看起来像"独立的平坦克隆", 实为技术伪影。
lib   <- colSums(obj@count.data)[obs]
ngene <- colSums(obj@count.data[, obs, drop = FALSE] > 0)
ref_prof <- colMeans(Mref)
qc <- data.frame(clone = paste0("C", seq_len(K)),
                 n_cell = as.integer(table(CL)),
                 umi_median = round(tapply(lib, CL, median)),
                 gene_median = round(tapply(ngene, CL, median)),
                 burden_median = round(tapply(burden, CL, median), 3),
                 dist_to_reference_profile = round(sapply(seq_len(K), function(g)
                   sqrt(sum((CM[g, ] - ref_prof)^2))), 4))
write.csv(qc, file.path(OUT_DIR, "clone_qc.csv"), row.names = FALSE)
say("\n    -- 克隆 QC (artifact 对照) --")
for (g in seq_len(K))
  say("    C%d: UMI 中位 %6.0f | 检出基因 %5.0f | burden %.3f | 谱到参考组距离 %.4f",
      g, qc$umi_median[g], qc$gene_median[g], qc$burden_median[g],
      qc$dist_to_reference_profile[g])
say("    全样品 UMI 中位 %.0f (q25 %.0f) | burden vs log10(UMI) r = %.3f (接近 0 = 非深度驱动)",
    median(lib), quantile(lib, .25), cor(burden, log10(lib)))
if (any(qc$umi_median < 0.25 * median(lib)))
  say("    [警告] 有克隆 UMI 中位 < 全样品 25%%, 该簇可能是低深度伪影, 不可当克隆")

# ---- k+1 诊断: 多切出来的那 1 簇是真克隆还是伪影 ----
if (K + 1 <= max(STEPS)) {
  K1 <- K + 1
  cl1 <- cutree(hc_g, k = K1)
  rl1 <- order(-as.integer(table(cl1))); nl1 <- integer(K1); nl1[rl1] <- seq_len(K1)
  cl1 <- factor(nl1[cl1], levels = seq_len(K1))
  cm1 <- do.call(rbind, lapply(seq_len(K1), function(g) colMeans(M[cl1 == g, , drop = FALSE])))
  say("\n    -- k=%d 诊断 (为什么不多切) -- silhouette = %.3f", K1, sil(as.integer(cl1), D_gn))
  for (g in seq_len(K1))
    say("    C%d: n=%3d | UMI 中位 %6.0f | burden %.3f | 谱到参考组距离 %.4f",
        g, sum(cl1 == g), median(lib[cl1 == g]), median(burden[cl1 == g]),
        sqrt(sum((cm1[g, ] - ref_prof)^2)))
  if (min(tapply(lib, cl1, median)) < 0.25 * median(lib))
    say("    -> 多出的簇 UMI 中位仅 %.0f (全样品 25%% 分位 %.0f), 判定为低深度伪影, 不采纳",
        min(tapply(lib, cl1, median)), quantile(lib, .25))
  write.csv(data.frame(clone = paste0("C", seq_len(K1)), n_cell = as.integer(table(cl1)),
                       umi_median = round(tapply(lib, cl1, median)),
                       burden_median = round(tapply(burden, cl1, median), 3)),
            file.path(OUT_DIR, "clone_kplus1_diagnostic.csv"), row.names = FALSE)
}

# 逐臂检验: K>=3 用 Kruskal-Wallis, 否则 Wilcoxon; 并与最大簇(C1)两两比较
pw <- function(x, y) tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
diff <- do.call(rbind, lapply(ARMS, function(a) {
  v <- M[, a]
  kw <- if (K >= 3) tryCatch(kruskal.test(v, CL)$p.value, error = function(e) NA_real_) else
    pw(v[CL == 1], v[CL == 2])
  row <- data.frame(arm = a, test = if (K >= 3) "kruskal" else "wilcox",
                    p = kw, mean_all = round(mean(v), 4))
  for (g in 2:K) {
    row[[paste0("delta_C", g, "_vs_C1")]] <- round(median(v[CL == g]) - median(v[CL == 1]), 4)
    row[[paste0("p_C", g, "_vs_C1")]]     <- pw(v[CL == g], v[CL == 1])
  }
  row
}))
diff$fdr <- p.adjust(diff$p, "BH")
diff <- diff[order(diff$fdr), ]
write.csv(diff, file.path(OUT_DIR, "clone_diff_arms.csv"), row.names = FALSE)
say("\n    逐臂差异检验(%s): FDR<0.05 的臂 %d / %d", diff$test[1],
    sum(diff$fdr < 0.05, na.rm = TRUE), nrow(diff))

# 每个克隆的 top gain / loss
summ <- do.call(rbind, lapply(seq_len(K), function(g) {
  d <- CM[g, ] - BASE_VAL
  up <- sort(d, decreasing = TRUE)[1:5]; dn <- sort(d)[1:5]
  data.frame(clone = paste0("C", g), n_cell = sum(CL == g),
             pct = round(100 * sum(CL == g) / length(obs), 1),
             burden_median = round(median(burden[CL == g]), 3),
             top_gains = paste(sprintf("%s%+.3f", names(up), up), collapse = " "),
             top_losses = paste(sprintf("%s%+.3f", names(dn), dn), collapse = " "))
}))
write.csv(summ, file.path(OUT_DIR, "clone_summary.csv"), row.names = FALSE)
say("\n    -- 各克隆 CNV 特征 (臂级偏离基线) --")
for (g in seq_len(K)) {
  d <- sort(CM[g, ] - BASE_VAL, decreasing = TRUE)
  say("    C%d (n=%d, %.1f%%): 高 %s", g, sum(CL == g), 100 * sum(CL == g) / length(obs),
      paste(sprintf("%s%+.3f", names(d)[1:5], d[1:5]), collapse = " "))
  say("                    低 %s", paste(sprintf("%s%+.3f", names(d)[(length(d)-4):length(d)],
                                                 d[(length(d)-4):length(d)]), collapse = " "))
}

# 克隆两两相似度
if (K >= 2) {
  cc <- cor(t(CM)); dd <- as.matrix(dist(CM))
  pw_df <- do.call(rbind, lapply(combn(K, 2, simplify = FALSE), function(p) {
    i <- p[1]; j <- p[2]
    # 发生变化的臂(|delta|>0.15)
    dl <- CM[j, ] - CM[i, ]
    oo <- order(abs(dl), decreasing = TRUE)[1:5]   # 注意: 取 order(), 不是 sort() 的数值
    data.frame(clone_a = paste0("C", i), clone_b = paste0("C", j),
               arm_profile_r = round(cc[i, j], 4),
               arm_profile_dist = round(dd[i, j], 4),
               n_arms_delta_gt_0.05 = sum(abs(dl) > 0.05),
               n_arms_changed = sum(abs(dl) > 0.15),
               biggest_changes = paste(sprintf("%s%+.3f", names(dl)[oo], dl[oo]), collapse = " "))
  }))
  write.csv(pw_df, file.path(OUT_DIR, "clone_pairwise.csv"), row.names = FALSE)
  say("\n    克隆两两: %s", paste(sprintf("%s vs %s r=%.3f, |delta|>0.15 的臂 %d 条",
      pw_df$clone_a, pw_df$clone_b, pw_df$arm_profile_r, pw_df$n_arms_changed), collapse = " | "))
}

# ==================================================================
# 6. 空间验证 + 坐标
# ==================================================================
SP <- NULL
if (file.exists(XY_CACHE)) {
  xy <- readRDS(XY_CACHE)
  key <- if ("barcode" %in% colnames(xy)) xy$barcode else xy[[1]]
  b <- sub(paste0("^", SAMPLE, "_"), "", obs)
  mm <- match(b, sub(paste0("^", SAMPLE, "_"), "", as.character(key)))
  if (sum(!is.na(mm)) == length(obs) &&
      all(c("x", "y") %in% colnames(xy)) && "sample" %in% colnames(xy)) {
    SP <- data.frame(cell = obs, x = xy$x[mm], y = xy$y[mm],
                     clone = CL, burden = burden, stringsAsFactors = FALSE)
    P <- as.matrix(SP[, c("x", "y")]); cc <- as.integer(SP$clone)
    kn <- get.knn(P, k = 6)$nn.index
    same <- function(v) mean(sapply(seq_len(nrow(P)), function(j) mean(v[kn[j, ]] == v[j])))
    set.seed(SEED); perm <- replicate(500, same(sample(cc)))
    s0 <- same(cc)
    pv <- (1 + sum(perm >= s0)) / (1 + length(perm))
    say("\n[6] 空间验证: 同克隆近邻比例 %.3f | 打乱后 %.3f | 富集 %.2fx | p=%.4f",
        s0, mean(perm), s0 / mean(perm), pv)
    write.csv(data.frame(metric = c("same_clone_neighbor_frac", "shuffled_mean", "enrichment", "p_perm"),
                         value = round(c(s0, mean(perm), s0 / mean(perm), pv), 4)),
              file.path(OUT_DIR, "spatial_coherence.csv"), row.names = FALSE)
  } else say("\n[6] 坐标未能匹配, 跳过空间图")
}

mem <- data.frame(cell = obs, clone = paste0("C", as.integer(CL)),
                  cnv_burden = round(burden, 4))
if (!is.null(SP)) mem <- cbind(mem, x = SP$x, y = SP$y)
write.csv(mem, file.path(OUT_DIR, "clone_membership.csv"), row.names = FALSE)

# ==================================================================
# 7. 树导出 (Newick, 供 phylogram / ggtree)
# ==================================================================
phy  <- ape::as.phylo(hc_g)
ape::write.tree(phy, file = file.path(OUT_DIR, "tree_allcells.nwk"))
cm_hc <- hclust(dist(CM), method = "ward.D2")   # 克隆层面树
ape::write.tree(ape::as.phylo(cm_hc), file = file.path(OUT_DIR, "tree_clones.nwk"))
say("\n[7] Newick: tree_allcells.nwk (%d tips), tree_clones.nwk (%d tips)",
    length(phy$tip.label), K)

# ==================================================================
# 8. 图 (全部 ASCII 标签 —— R 设备无中文字体)
# ==================================================================
pal_div <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(101)   # 发散: loss=蓝, gain=红

# --- Fig1: k 选择 ---
kdf <- rbind(data.frame(k = ks$k, silhouette = ks$sil_arm,  res = "arm (42 features)"),
             data.frame(k = ks$k, silhouette = ks$sil_gene,
                        res = paste0("gene PCA (", NPC, " PCs)")))
RES_G <- paste0("gene PCA (", NPC, " PCs)")
p1 <- ggplot(kdf, aes(k, silhouette, colour = res)) +
  geom_line(linewidth = .8) + geom_point(size = 2.4) +
  geom_vline(xintercept = K, linetype = 2, colour = "grey35") +
  annotate("text", x = K, y = max(kdf$silhouette), label = paste0("k* = ", K),
           vjust = 1.4, hjust = -.1, size = 3.2) +
  scale_colour_manual(values = setNames(c("#0072B2", "#D55E00"),
                                        c("arm (42 features)", RES_G)), name = NULL) +
  scale_x_continuous(breaks = STEPS) +
  labs(title = sprintf("%s: choosing the number of clones", SAMPLE),
       subtitle = "silhouette width vs k (higher = better separated)",
       x = "k (number of clones)", y = "mean silhouette width") +
  theme_classic(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(OUT_DIR, "Fig1_k_selection.png"), p1, width = 6.2, height = 4, dpi = 300, bg = "white")

# --- Fig2: 细胞树 + 克隆色带 ---
png(file.path(OUT_DIR, "Fig2_clone_tree.png"), width = 2600, height = 1500, res = 260)
par(mar = c(3.2, 2.5, 3, 1))
plot(hc_g, labels = FALSE, hang = -1, main = sprintf(
  "%s: %d epithelial cells, ward.D2 on CNV (%d PCs) - clone colour strip",
  SAMPLE, length(obs), NPC), xlab = "", sub = "", yaxs = "i", cex.main = 1.0)
usr <- par("usr"); yr <- usr[4] - usr[3]
ord <- hc_g$order
cols <- CLONE_COL[as.integer(CL)][ord]
rect(seq_along(ord) - .5, usr[3] - .055 * yr, seq_along(ord) + .5, usr[3] - .015 * yr,
     col = cols, border = NA, xpd = NA)
legend("topright", legend = paste0("C", seq_len(K), " (n=", as.integer(table(CL)), ")"),
       fill = CLONE_COL[seq_len(K)], border = NA, bty = "n", cex = .85)
invisible(dev.off())

# --- Fig3: 克隆 x 臂 热图 ---
hm <- do.call(rbind, lapply(seq_len(K), function(g) data.frame(
  arm = factor(ARMS, levels = ARMS), clone = paste0("C", g), dev = CM[g, ] - BASE_VAL)))
hm$sig <- ""
for (a in diff$arm[diff$fdr < 0.05]) {
  i <- which(diff$arm == a)
  for (g in 2:K) {
    dv <- diff[[paste0("delta_C", g, "_vs_C1")]][i]
    if (!is.na(dv) && abs(dv) > 0.02)
      hm$sig[hm$arm == a & hm$clone == paste0("C", g)] <- "*"
  }
}
lim <- max(abs(hm$dev))
p3 <- ggplot(hm, aes(arm, clone, fill = dev)) +
  geom_tile(colour = "white", linewidth = .35) +
  geom_text(aes(label = sig), size = 4, vjust = .78, colour = "grey15") +
  scale_fill_gradientn(colours = pal_div, limits = c(-lim, lim),
                       name = "CNV level\n(deviation from\ndiploid 1.0036)") +
  labs(title = sprintf("%s: mean CNV by clone and chromosome arm", SAMPLE),
       subtitle = "* = FDR<0.05 vs C1 (Wilcoxon/Kruskal-Wallis, BH)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 8),
        panel.grid = element_blank())
ggsave(file.path(OUT_DIR, "Fig3_clone_arm_heatmap.png"), p3,
       width = 9.5, height = 1.9 + .6 * K, dpi = 300, bg = "white", limitsize = FALSE)

# --- Fig4: 空间图 ---
if (!is.null(SP)) {
  sp <- SP; sp$clone <- factor(paste0("C", sp$clone), levels = paste0("C", seq_len(K)))
  p4a <- ggplot(sp, aes(x, y, colour = clone)) + geom_point(size = .55) +
    scale_colour_manual(values = CLONE_COL[seq_len(K)], name = NULL) +
    coord_fixed() + labs(title = "Clone assignment") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = .5, size = 11, face = "bold"),
          legend.position = "bottom")
  p4b <- ggplot(sp, aes(x, y, colour = burden)) + geom_point(size = .55) +
    scale_colour_viridis_c(option = "inferno", name = "CNV burden") +
    coord_fixed() + labs(title = "CNV burden") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = .5, size = 11, face = "bold"),
          legend.position = "bottom", legend.key.width = unit(1.2, "cm"),
          legend.key.height = unit(.28, "cm"))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p4 <- patchwork::wrap_plots(p4a, p4b, ncol = 2) +
      patchwork::plot_annotation(title = sprintf("%s: clones on tissue", SAMPLE),
                                 theme = theme(plot.title = element_text(hjust = .5, face = "bold")))
    ggsave(file.path(OUT_DIR, "Fig4_clone_spatial.png"), p4, width = 11, height = 5.6,
           dpi = 300, bg = "white")
  } else {
    ggsave(file.path(OUT_DIR, "Fig4_clone_spatial.png"), p4a, width = 6, height = 5.6,
           dpi = 300, bg = "white")
  }
}

# --- Fig5: burden 分布 ---
b5 <- data.frame(clone = factor(paste0("C", as.integer(CL)), levels = paste0("C", seq_len(K))),
                 burden = burden)
p5 <- ggplot(b5, aes(clone, burden, fill = clone)) +
  geom_boxplot(outlier.size = .4, width = .55, colour = "grey25") +
  scale_fill_manual(values = CLONE_COL[seq_len(K)], guide = "none") +
  labs(title = sprintf("%s: CNV burden per cell by clone", SAMPLE),
       subtitle = "fraction of genes with |CNV - diploid| > 0.10",
       x = NULL, y = "CNV burden") +
  theme_classic(base_size = 11)
ggsave(file.path(OUT_DIR, "Fig5_clone_burden.png"), p5, width = 5.2, height = 4, dpi = 300, bg = "white")

say("\n[8] 图: Fig1_k_selection / Fig2_clone_tree / Fig3_clone_arm_heatmap / Fig4_clone_spatial / Fig5_clone_burden")
say("\n表: clone_k_selection / hmm_state_calibration / clone_membership / clone_summary")
say("    clone_arm_mean / clone_arm_deviation / clone_diff_arms / clone_pairwise")
say("树: tree_allcells.nwk / tree_clones.nwk")
say("\n输出目录: %s", OUT_DIR)
say("========== Done ==========")
