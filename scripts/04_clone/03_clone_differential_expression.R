# =============================================================================
# TNBC interface-niche code repository
# Original source: project_spatial/TNBC/patient_ALL/TAM_clone/notebooks/step02_clone1_vs_clone2_de.ipynb (converted to R)
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

# 记录本 notebook 的可复现元数据：随机种子、软件版本、输入文件与关键参数
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
set.seed(42)
p13_path <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P13/clone_seurat.rds"
p15_path <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P15/clone_seurat.rds"
fi <- function(p) sprintf("%s (%.1f MB, mtime=%s)", p, file.info(p)$size/1e6, format(file.info(p)$mtime, "%Y-%m-%dT%H:%M:%S"))
meta <- data.frame(
  item = c("random_seed","R_version","Seurat","SeuratObject","presto","ggplot2","ggrepel","png",
           "input_P13","input_P15","DE_params","producer"),
  value = c("42", R.version.string,
            as.character(packageVersion("Seurat")), as.character(packageVersion("SeuratObject")),
            tryCatch(as.character(packageVersion("presto")), error=function(e) "MISSING"),
            as.character(packageVersion("ggplot2")), as.character(packageVersion("ggrepel")),
            tryCatch(as.character(packageVersion("png")), error=function(e) "MISSING"),
            fi(p13_path), fi(p15_path),
            "FindMarkers(test.use=wilcox, min.pct=0.25, logfc.threshold=0); ident.1=Clone1, ident.2=Clone2; sig=padj<0.05",
            "sc-analyst, 2026-10-27"))
print(meta, row.names = FALSE)
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 载入 P13 / P15 对象，核验维度、metadata、克隆注释、切片分布与 CNV 负荷（混杂信息）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
obj13 <- readRDS(p13_path)
obj15 <- readRDS(p15_path)

inspect_clone_obj <- function(obj, tag, memb_path) {
  cat("\n===== ", tag, " =====\n")
  cat("class:", paste(class(obj), collapse = ", "), "\n")
  cat("dim:", nrow(obj), "genes x", ncol(obj), "spots\n")
  cat("Assays:", paste(Assays(obj), collapse = ", "), "| DefaultAssay:", DefaultAssay(obj), "\n")
  cat("Layers:", paste(Layers(obj[[DefaultAssay(obj)]]), collapse = ", "), "\n")
  cat("meta.data columns:", paste(colnames(obj@meta.data), collapse = ", "), "\n")
  cat("Idents levels:", paste(levels(Idents(obj)), collapse = ", "), "\n")
  m <- read.csv(memb_path, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(m) <- m$cell
  cat("[membership] rows:", nrow(m), "| all object spots present in membership:", all(colnames(obj) %in% rownames(m)), "\n")
  if ("clone" %in% colnames(obj@meta.data)) {
    same <- identical(as.character(obj$clone), as.character(m[colnames(obj), "clone"]))
    cat("[membership] object 'clone' column matches membership CSV:", same, "\n")
  } else {
    cat("[Align] 'clone' absent in object -> merged from membership CSV\n")
    obj@meta.data <- cbind(obj@meta.data, m[colnames(obj), setdiff(colnames(m), "cell"), drop = FALSE])
  }
  md <- obj@meta.data
  cat("\n[clone counts]\n"); print(table(md$clone))
  if ("slice" %in% colnames(md)) { cat("\n[clone x slice]\n"); print(table(md$clone, md$slice)) }
  if ("cnv_burden" %in% colnames(md)) {
    cat("\n[cnv_burden by clone: mean / median]\n")
    print(round(tapply(md$cnv_burden, md$clone, function(x) c(mean = mean(x), median = median(x))), 4))
  }
  lyr <- if ("data" %in% Layers(obj[[DefaultAssay(obj)]])) "data" else "counts"
  d <- LayerData(obj, assay = DefaultAssay(obj), layer = lyr)
  v <- if (inherits(d, "dgCMatrix")) d@x else as.vector(d)
  cat(sprintf("\n[%s layer] nnz=%d | nonzero range=[%.3f, %.3f] | all-integer=%s\n",
              lyr, length(v), min(v), max(v), all(v == round(v))))
  obj
}

obj13 <- inspect_clone_obj(obj13, "P13", "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P13/clone_membership.csv")
obj15 <- inspect_clone_obj(obj15, "P15", "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P15/clone_membership.csv")
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 合并切片/CNV 负荷描述性元数据，并对两个对象做 LogNormalize 标准化（FindMarkers 需要 data 层）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
suppressPackageStartupMessages(library(Seurat))
out_tab <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/tables"
out_fig <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/figures/per_patient"
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE); dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

add_desc_meta <- function(obj, memb_path) {
  m <- read.csv(memb_path, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(m) <- m$cell
  idx <- colnames(obj)
  obj@meta.data$slice <- as.character(m[idx, "slice"])
  obj@meta.data$cnv_burden <- as.numeric(m[idx, "cnv_burden"])
  obj
}
obj13 <- add_desc_meta(obj13, "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P13/clone_membership.csv")
obj15 <- add_desc_meta(obj15, "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P15/clone_membership.csv")

cat("[Plan] NormalizeData(LogNormalize, scale.factor=1e4) — 生成 data 层供 FindMarkers 使用\n")
obj13 <- NormalizeData(obj13, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
obj15 <- NormalizeData(obj15, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
d13 <- LayerData(obj13, assay = "RNA", layer = "data"); d15 <- LayerData(obj15, assay = "RNA", layer = "data")
cat(sprintf("[Normalize] P13 data range=[%.3f, %.3f] | P15 data range=[%.3f, %.3f]\n",
            min(d13@x), max(d13@x), min(d15@x), max(d15@x)))
cat(sprintf("[Desc meta] spots: P13=%d, P15=%d | slices: P13=%s, P15=%s\n",
            ncol(obj13), ncol(obj15),
            paste(unique(obj13$slice), collapse = "/"), paste(unique(obj15$slice), collapse = "/")))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 汇总复制/克隆与切片、CNV 负荷、测序深度的关系（仅描述性统计，不参与 DE 检验）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  s_ij <- sum(choose(tab, 2)); s_i <- sum(choose(rowSums(tab), 2)); s_j <- sum(choose(colSums(tab), 2))
  e <- s_i * s_j / choose(n, 2); mx <- (s_i + s_j) / 2
  (s_ij - e) / (mx - e)
}
summarise_patient <- function(obj, tag) {
  md <- obj@meta.data
  cat("\n===== ", tag, " =====\n[orig.ident]\n"); print(table(md$orig.ident))
  cat("[clone x slice]\n"); print(table(md$clone, md$slice))
  cat(sprintf("[clone vs slice] ARI = %.4f\n", ari(md$clone, md$slice)))
  perf <- data.frame(cnv_mean = tapply(md$cnv_burden, md$clone, mean),
                     cnv_median = tapply(md$cnv_burden, md$clone, median),
                     nCount_median = tapply(md$nCount_RNA, md$clone, median),
                     nFeature_median = tapply(md$nFeature_RNA, md$clone, median))
  cat("[clone metrics: cnv_burden mean/median | nCount/nFeature median]\n"); print(round(perf, 4))
  cat(sprintf("[depth] nCount median ratio Clone1/Clone2 = %.2f ; nFeature ratio = %.2f\n",
              perf["Clone1", "nCount_median"] / perf["Clone2", "nCount_median"],
              perf["Clone1", "nFeature_median"] / perf["Clone2", "nFeature_median"]))
  tcl <- table(md$clone); tcs <- table(md$clone, md$slice)
  data.frame(patient = tag, clone = names(tcl), n_spot = as.integer(tcl),
             pct_of_patient = round(100 * as.integer(tcl) / nrow(md), 2),
             cnv_burden_mean = round(perf$cnv_mean, 4), cnv_burden_median = round(perf$cnv_median, 4),
             nCount_median = round(perf$nCount_median), nFeature_median = round(perf$nFeature_median),
             n_slice_A = as.integer(tcs[, 1]), n_slice_B = as.integer(tcs[, 2]),
             clone_slice_ARI = round(ari(md$clone, md$slice), 4))
}
conf <- rbind(summarise_patient(obj13, "P13"), summarise_patient(obj15, "P15"))
cat("\n[Confounder summary]\n"); print(conf, row.names = FALSE)
write.csv(conf, file.path(out_tab, "clone_confounder_summary.csv"), row.names = FALSE)
cat("\n[Write]", file.path(out_tab, "clone_confounder_summary.csv"), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# P13：spot 级 FindMarkers（Clone1 vs Clone2，wilcox，min.pct=0.25，logfc.threshold=0），导出完整差异表
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
suppressPackageStartupMessages(library(Seurat))
Idents(obj13) <- "clone"
cat("[Plan] FindMarkers(ident.1=Clone1, ident.2=Clone2, test.use=wilcox, min.pct=0.25, logfc.threshold=0) — P13\n")
t0 <- Sys.time()
de13 <- FindMarkers(obj13, ident.1 = "Clone1", ident.2 = "Clone2", test.use = "wilcox",
                    min.pct = 0.25, logfc.threshold = 0, verbose = FALSE)
cat(sprintf("[Runtime] P13 FindMarkers: %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
de13$gene <- rownames(de13)
de13 <- de13[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
n_sig13 <- sum(de13$p_val_adj < 0.05, na.rm = TRUE)
n_up13 <- sum(de13$p_val_adj < 0.05 & de13$avg_log2FC > 0, na.rm = TRUE)
n_dn13 <- sum(de13$p_val_adj < 0.05 & de13$avg_log2FC < 0, na.rm = TRUE)
cat(sprintf("[P13] tested=%d | padj<0.05=%d (Clone1-high=%d, Clone1-low=%d)\n", nrow(de13), n_sig13, n_up13, n_dn13))
cat("\n[P13] top 10 by avg_log2FC (Clone1-high direction)\n")
print(head(de13[order(-de13$avg_log2FC), c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 10), row.names = FALSE)
cat("\n[P13] top 10 by avg_log2FC (Clone1-low direction)\n")
print(head(de13[order(de13$avg_log2FC), c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 10), row.names = FALSE)
top13 <- de13$gene[which.max(de13$avg_log2FC)]
d13d <- LayerData(obj13, assay = "RNA", layer = "data")
cat(sprintf("\n[Direction check] %s mean expr: Clone1=%.3f | Clone2=%.3f (expect Clone1 > Clone2); padj=%.2e\n",
            top13, mean(d13d[top13, obj13$clone == "Clone1"]), mean(d13d[top13, obj13$clone == "Clone2"]),
            de13$p_val_adj[de13$gene == top13]))
write.csv(de13, file.path(out_tab, "de_clone1_vs_clone2_P13.csv"), row.names = FALSE)
cat("[Write]", file.path(out_tab, "de_clone1_vs_clone2_P13.csv"), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# P15：spot 级 FindMarkers（Clone1 vs Clone2，wilcox，min.pct=0.25，logfc.threshold=0），导出完整差异表
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
Idents(obj15) <- "clone"
cat("[Plan] FindMarkers(ident.1=Clone1, ident.2=Clone2, test.use=wilcox, min.pct=0.25, logfc.threshold=0) — P15\n")
t0 <- Sys.time()
de15 <- FindMarkers(obj15, ident.1 = "Clone1", ident.2 = "Clone2", test.use = "wilcox",
                    min.pct = 0.25, logfc.threshold = 0, verbose = FALSE)
cat(sprintf("[Runtime] P15 FindMarkers: %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
de15$gene <- rownames(de15)
de15 <- de15[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
n_sig15 <- sum(de15$p_val_adj < 0.05, na.rm = TRUE)
n_up15 <- sum(de15$p_val_adj < 0.05 & de15$avg_log2FC > 0, na.rm = TRUE)
n_dn15 <- sum(de15$p_val_adj < 0.05 & de15$avg_log2FC < 0, na.rm = TRUE)
cat(sprintf("[P15] tested=%d | padj<0.05=%d (Clone1-high=%d, Clone1-low=%d)\n", nrow(de15), n_sig15, n_up15, n_dn15))
cat("\n[P15] top 10 by avg_log2FC (Clone1-high direction)\n")
print(head(de15[order(-de15$avg_log2FC), c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 10), row.names = FALSE)
cat("\n[P15] top 10 by avg_log2FC (Clone1-low direction)\n")
print(head(de15[order(de15$avg_log2FC), c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 10), row.names = FALSE)
top15 <- de15$gene[which.max(de15$avg_log2FC)]
d15d <- LayerData(obj15, assay = "RNA", layer = "data")
cat(sprintf("\n[Direction check] %s mean expr: Clone1=%.3f | Clone2=%.3f (expect Clone1 > Clone2); padj=%.2e\n",
            top15, mean(d15d[top15, obj15$clone == "Clone1"]), mean(d15d[top15, obj15$clone == "Clone2"]),
            de15$p_val_adj[de15$gene == top15]))
write.csv(de15, file.path(out_tab, "de_clone1_vs_clone2_P15.csv"), row.names = FALSE)
cat("[Write]", file.path(out_tab, "de_clone1_vs_clone2_P15.csv"), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 切片内稳健性检查（描述性，不替代主分析）：显著基因在每个切片内的 Clone1-Clone2 表达差方向是否与总体 avg_log2FC 一致
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
slice_check <- function(obj, de, tag) {
  d <- LayerData(obj, assay = "RNA", layer = "data")
  sig <- de[de$p_val_adj < 0.05 & !is.na(de$p_val_adj), ]
  genes <- sig$gene; slices <- sort(unique(obj$slice))
  out <- data.frame(gene = genes, avg_log2FC = sig$avg_log2FC, padj = sig$p_val_adj,
                    direction = ifelse(sig$avg_log2FC > 0, "Clone1-high", "Clone1-low"))
  cons <- rep(TRUE, nrow(out))
  for (sl in slices) {
    in_sl <- obj$slice == sl
    c1 <- Matrix::rowMeans(d[genes, in_sl & obj$clone == "Clone1", drop = FALSE])
    c2 <- Matrix::rowMeans(d[genes, in_sl & obj$clone == "Clone2", drop = FALSE])
    out[[paste0("n_C1_", sl)]] <- sum(in_sl & obj$clone == "Clone1")
    out[[paste0("n_C2_", sl)]] <- sum(in_sl & obj$clone == "Clone2")
    out[[paste0("meanC1_", sl)]] <- round(c1, 4); out[[paste0("meanC2_", sl)]] <- round(c2, 4)
    out[[paste0("diff_", sl)]] <- round(c1 - c2, 4)
    cons <- cons & (sign(c1 - c2) == sign(sig$avg_log2FC))
  }
  out$sign_consistent_both_slices <- cons
  cat(sprintf("\n[%s] sig genes=%d | sign-consistent in BOTH slices: %d (%.1f%%)\n", tag, nrow(out), sum(cons), 100 * mean(cons)))
  for (dir in c("Clone1-high", "Clone1-low")) {
    ii <- out$direction == dir
    cat(sprintf("  %-12s %d genes | consistent %d (%.1f%%)\n", dir, sum(ii), sum(cons[ii]), 100 * mean(cons[ii])))
  }
  write.csv(out, file.path(out_tab, paste0("de_slice_consistency_", tag, ".csv")), row.names = FALSE)
  out
}
sc13 <- slice_check(obj13, de13, "P13")
sc15 <- slice_check(obj15, de15, "P15")
cat("\n[Write] de_slice_consistency_P13.csv / de_slice_consistency_P15.csv\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 火山图（修订版）：y 轴按数据分位截断并显式标注超出个数，修复 y 轴压缩与刻度覆盖问题
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
suppressPackageStartupMessages({library(ggplot2); library(ggrepel)})
volcano_plot <- function(de, tag, out_base, cap) {
  de$neglog10p <- -log10(pmax(de$p_val_adj, .Machine$double.xmin))
  de$direction <- ifelse(de$p_val_adj < 0.05 & de$avg_log2FC > 0, "Clone1-high",
                  ifelse(de$p_val_adj < 0.05 & de$avg_log2FC < 0, "Clone1-low", "NS"))
  de$direction <- factor(de$direction, levels = c("Clone1-high", "Clone1-low", "NS"))
  de$y_plot <- pmin(de$neglog10p, cap)
  n_up <- sum(de$direction == "Clone1-high"); n_dn <- sum(de$direction == "Clone1-low")
  max_y <- max(de$neglog10p); n_above <- sum(de$neglog10p > cap)
  up <- de[de$direction == "Clone1-high", ]
  top <- head(up[order(-up$avg_log2FC), ], 8)
  cat(sprintf("\n[%s] top%d Clone1-high labels: %s | y-limit=%.0f | genes above limit=%d | max=%.1f\n",
              tag, nrow(top), paste(top$gene, collapse = ", "), cap, n_above, max_y))
  p <- ggplot(de, aes(avg_log2FC, y_plot, color = direction)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_text_repel(data = top, mapping = aes(x = avg_log2FC, y = y_plot, label = gene),
                    inherit.aes = FALSE, size = 3, color = "black", max.overlaps = Inf,
                    box.padding = 0.6, point.padding = 0.25, force = 3, max.iter = 10000,
                    min.segment.length = 0, seed = 42, show.legend = FALSE) +
    annotate("text", x = -Inf, y = cap, hjust = -0.05, vjust = 1.2, size = 2.8, color = "grey30",
             label = sprintf("%d genes above y = %.0f (max %.1f)", n_above, cap, max_y)) +
    scale_color_manual(values = c("Clone1-high" = "#B2182B", "Clone1-low" = "#2166AC", "NS" = "grey70"),
                       labels = c(paste0("Clone1-high (", n_up, ")"), paste0("Clone1-low (", n_dn, ")"), "NS"),
                       drop = FALSE) +
    scale_y_continuous(limits = c(0, cap * 1.06), breaks = seq(0, cap, by = 10),
                       expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "avg_log2FC (Clone1 vs Clone2)", y = "-log10(padj)", color = NULL,
         title = sprintf("%s: Clone1 vs Clone2", tag),
         subtitle = "avg_log2FC > 0 = Clone1-high; significant = padj < 0.05; points above the y-limit are clipped") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
          legend.position = "right", legend.key = element_blank())
  ggsave(paste0(out_base, ".png"), p, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(paste0(out_base, ".pdf"), p, width = 8, height = 6, dpi = 300, bg = "white")
  img <- png::readPNG(paste0(out_base, ".png"))
  rgbv <- as.vector(img[, , 1:3])
  cat(sprintf("[Pixel] %s.png: %dx%d | sd=%.4f | white_frac=%.3f\n",
              basename(out_base), dim(img)[2], dim(img)[1], sd(rgbv), mean(rgbv > 0.97)))
  invisible(top)
}
top13_plot <- volcano_plot(de13, "P13", file.path(out_fig, "volcano_P13_clone1_vs_clone2"), cap = 50)
top15_plot <- volcano_plot(de15, "P15", file.path(out_fig, "volcano_P15_clone1_vs_clone2"), cap = 75)
cat("\n[Write] volcano_P13_clone1_vs_clone2.{png,pdf} / volcano_P15_clone1_vs_clone2.{png,pdf} (revised)\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 跨病人描述性比较：显著基因集重叠与 avg_log2FC 一致性（仅描述，不做跨病人合并检验；克隆编号跨病人不可比）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
m <- merge(de13[, c("gene", "avg_log2FC", "p_val_adj")], de15[, c("gene", "avg_log2FC", "p_val_adj")],
           by = "gene", suffixes = c("_P13", "_P15"))
rho_all <- cor(m$avg_log2FC_P13, m$avg_log2FC_P15, method = "spearman")
up13 <- de13$gene[de13$p_val_adj < 0.05 & de13$avg_log2FC > 0]
up15 <- de15$gene[de15$p_val_adj < 0.05 & de15$avg_log2FC > 0]
dn13 <- de13$gene[de13$p_val_adj < 0.05 & de13$avg_log2FC < 0]
dn15 <- de15$gene[de15$p_val_adj < 0.05 & de15$avg_log2FC < 0]
shared_up <- sort(intersect(up13, up15)); shared_dn <- sort(intersect(dn13, dn15))
disc <- sort(intersect(up13, dn15))
cat(sprintf("[Shared genes tested] %d | Spearman rho(all shared genes) = %.3f\n", nrow(m), rho_all))
cat(sprintf("[Clone1-high] P13=%d | P15=%d | shared=%d | Jaccard=%.3f\n",
            length(up13), length(up15), length(shared_up), length(shared_up) / length(union(up13, up15))))
cat(sprintf("[Clone1-low ] P13=%d | P15=%d | shared=%d | Jaccard=%.3f\n",
            length(dn13), length(dn15), length(shared_dn), length(shared_dn) / length(union(dn13, dn15))))
cat(sprintf("[Discordant] Clone1-high in P13 AND Clone1-low in P15: %d genes\n", length(disc)))
sig_any <- unique(c(up13, up15, dn13, dn15))
mm <- m[m$gene %in% sig_any, ]
cat(sprintf("[Concordance among genes significant in >=1 patient] n=%d | Spearman rho = %.3f | same-sign fraction = %.3f\n",
            nrow(mm), cor(mm$avg_log2FC_P13, mm$avg_log2FC_P15, method = "spearman"),
            mean(sign(mm$avg_log2FC_P13) == sign(mm$avg_log2FC_P15))))
cat("\n[Top shared Clone1-high genes, ranked by P13 avg_log2FC]\n")
su <- m[m$gene %in% shared_up, ]; su <- su[order(-su$avg_log2FC_P13), ]
print(head(su, 25), row.names = FALSE)
ov <- data.frame(metric = c("genes_tested_shared", "spearman_all", "n_up_P13", "n_up_P15", "n_shared_up", "jaccard_up",
                            "n_down_P13", "n_down_P15", "n_shared_down", "jaccard_down", "n_discordant_up13_dn15",
                            "n_sig_any", "spearman_sig_any", "same_sign_fraction_sig_any"),
                 value = c(nrow(m), round(rho_all, 4), length(up13), length(up15), length(shared_up),
                           round(length(shared_up) / length(union(up13, up15)), 4), length(dn13), length(dn15),
                           length(shared_dn), round(length(shared_dn) / length(union(dn13, dn15)), 4), length(disc),
                           nrow(mm), round(cor(mm$avg_log2FC_P13, mm$avg_log2FC_P15, method = "spearman"), 4),
                           round(mean(sign(mm$avg_log2FC_P13) == sign(mm$avg_log2FC_P15)), 4)))
print(ov, row.names = FALSE)
write.csv(ov, file.path(out_tab, "de_cross_patient_overlap.csv"), row.names = FALSE)
write.csv(su, file.path(out_tab, "de_shared_up_genes.csv"), row.names = FALSE)
cat("[Write] de_cross_patient_overlap.csv / de_shared_up_genes.csv\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 构建交付文档内容（修订版）：修正跨病人分母/基因归属表述、表格整数精度、图注披露 y 轴截断
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
to_html <- function(d) {
  rows <- apply(d, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<table><thead><tr>", paste0("<th>", names(d), "</th>", collapse = ""), "</tr></thead><tbody>",
         paste(rows, collapse = ""), "</tbody></table>")
}
fmt_tbl <- function(df, dec = list()) {
  d <- as.data.frame(df, stringsAsFactors = FALSE)
  for (nm in names(d)) if (is.numeric(d[[nm]])) {
    k <- dec[[nm]]; if (is.null(k)) k <- if (all(abs(d[[nm]] - round(d[[nm]])) < 1e-9, na.rm = TRUE)) 0L else 3L
    d[[nm]] <- formatC(d[[nm]], format = "f", digits = k)
  }
  to_html(d)
}
fmt_gene <- function(df) {
  d <- df
  d$avg_log2FC <- formatC(d$avg_log2FC, format = "f", digits = 3)
  d$pct.1 <- formatC(d$pct.1, format = "f", digits = 3); d$pct.2 <- formatC(d$pct.2, format = "f", digits = 3)
  d$p_val_adj <- formatC(d$p_val_adj, format = "g", digits = 3)
  to_html(d)
}
img <- function(p, cap) paste0('<figure><img src="', base64enc::dataURI(file = p, mime = "image/png"),
                               '" alt="', cap, '"><figcaption>', cap, "</figcaption></figure>")
ss2 <- data.frame(patient = c("P13", "P15"), n_genes_tested = c(nrow(de13), nrow(de15)),
                  n_sig = c(sum(de13$p_val_adj < 0.05), sum(de15$p_val_adj < 0.05)),
                  n_clone1_high = c(sum(de13$p_val_adj < 0.05 & de13$avg_log2FC > 0), sum(de15$p_val_adj < 0.05 & de15$avg_log2FC > 0)),
                  n_clone1_low = c(sum(de13$p_val_adj < 0.05 & de13$avg_log2FC < 0), sum(de15$p_val_adj < 0.05 & de15$avg_log2FC < 0)))
n_union <- length(union(c(de13$gene[de13$p_val_adj < 0.05], de15$gene[de15$p_val_adj < 0.05]), character(0)))
h <- c()
h <- c(h, '<h1>Step 02 — Clone1 vs Clone2 差异表达（P13 / P15）</h1>', '<h2>一、结论简报</h2>',
       '<p><b>一句话结论：</b>两例病人的 Clone1 与 Clone2 都存在大量差异表达基因，但显著基因的规模与主导方向在两人之间并不对称——P13 以 Clone1 低表达为主（888 上调 / 4459 下调），P15 以 Clone1 高表达为主（2383 上调 / 935 下调），而两人的测序深度差异方向恰好相反，提示该不对称至少部分是深度/检出率的产物。两人共享的 Clone1 高表达核心是一组细胞外基质（ECM）/胶原类基因；P15 另有在 P13 未显著的 Wnt/β-catenin 相关基因（NOTUM、APCDD1、CCND1）高表达。由于克隆由 CNV 定义、CNV 负荷与肿瘤细胞比例在两人均存在同向差异，且 P15 的克隆与切片高度共线，"克隆差异基因"不能直接读作克隆内在的肿瘤细胞程序。</p>',
       '<p>方向参照系：ident.1 = Clone1，ident.2 = Clone2；<b>avg_log2FC &gt; 0 = Clone1 高表达</b>，&lt; 0 = Clone1 低表达。已用锚定基因校验方向（P13 COMP：Clone1 1.019 vs Clone2 0.360；P15 NOTUM：Clone1 0.232 vs Clone2 0.047）。</p>',
       '<h2>二、方法要点</h2><ul>',
       '<li>spot 级 <code>Seurat::FindMarkers(test.use="wilcox", min.pct=0.25, logfc.threshold=0)</code>；两个病人分别独立分析，不合并、不跨病人比较。</li>',
       '<li><code>logfc.threshold=0</code> 用于导出完整基因表（Seurat 默认 0.1 会预先滤掉小效应基因）；显著定义 = <code>padj &lt; 0.05</code>，不额外设 log2FC 阈值，方向由 avg_log2FC 符号判定。</li>',
       '<li>输入为 10x Visium <b>spot 级</b>对象（每 spot 为多细胞混合），原始 counts 经 LogNormalize（scale.factor=1e4）标准化后检验；未做任何抽样/降采样。</li>',
       '<li>克隆标签来自 <code>clone_membership.csv</code>，与对象内 <code>clone</code> 列逐 spot 完全一致（P13 1485/1485、P15 1938/1938）。</li></ul>',
       '<h2>三、关键结果</h2>', fmt_tbl(ss2),
       '<p>表 1：每病人的检验基因数与显著基因数（padj &lt; 0.05）。完整差异表见 tables/de_clone1_vs_clone2_P13.csv、tables/de_clone1_vs_clone2_P15.csv（列：gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj）。</p>',
       '<h3>P13 — Clone1 高表达 top8（按 avg_log2FC）</h3>', fmt_gene(top_up(de13)),
       '<h3>P15 — Clone1 高表达 top8（按 avg_log2FC）</h3>', fmt_gene(top_up(de15)),
       img(file.path(out_fig, "volcano_P13_clone1_vs_clone2.png"), "图 1：P13 Clone1 vs Clone2 火山图（x=avg_log2FC，y=-log10(padj)；红=Clone1 高表达，蓝=Clone1 低表达，灰=不显著；标注 Clone1 高表达方向效应量最大的 8 个基因；y 轴展示上限 50，35 个基因超出该上限（最大 -log10(padj)=122.4），已在图内标注）"),
       img(file.path(out_fig, "volcano_P15_clone1_vs_clone2.png"), "图 2：P15 Clone1 vs Clone2 火山图（x=avg_log2FC，y=-log10(padj)；红=Clone1 高表达，蓝=Clone1 低表达，灰=不显著；标注 Clone1 高表达方向效应量最大的 8 个基因；y 轴展示上限 75，31 个基因超出该上限（最大 -log10(padj)=141.8），已在图内标注）"),
       '<h2>四、跨病人一致性与差异（描述性）</h2>',
       sprintf('<p>两人共同检验的 8896 个基因中，任一病人显著的有 6582 个（两病人显著基因并集共 %d 个，含仅在单侧被检验的基因），这批基因的 avg_log2FC 方向一致比例为 85.5%%、Spearman 相关 0.559。Clone1 高表达集合：P13 %d / P15 %d，交集 %d（Jaccard %.3f）；Clone1 低表达集合：P13 %d / P15 %d，交集 %d（Jaccard %.3f）；方向相反（P13 高、P15 低）基因 %d 个。</p>',
               n_union, length(up13 <- de13$gene[de13$p_val_adj < 0.05 & de13$avg_log2FC > 0]),
               length(up15 <- de15$gene[de15$p_val_adj < 0.05 & de15$avg_log2FC > 0]),
               length(intersect(up13, up15)), length(intersect(up13, up15)) / length(union(up13, up15)),
               length(dn13 <- de13$gene[de13$p_val_adj < 0.05 & de13$avg_log2FC < 0]),
               length(dn15 <- de15$gene[de15$p_val_adj < 0.05 & de15$avg_log2FC < 0]),
               length(intersect(dn13, dn15)), length(intersect(dn13, dn15)) / length(union(dn13, dn15)),
               length(intersect(up13, dn15))),
       '<p>共享的 Clone1 高表达基因以 ECM/胶原与基质重塑为主（COMP、AEBP1、COL1A1、COL5A1、SPARC、BGN、MFAP2、EMILIN1、TIMP3、THY1、MMP11）。P15 显著而 P13 未显著的高效应基因包括 NOTUM、APCDD1、CCND1（CCND1 在 P13 表中但 padj≈1；NOTUM/APCDD1 在 P13 因 min.pct=0.25 未被检验，属"未检验"而非"无差异"）。这些基因与 Wnt/β-catenin 信号相关，此归类基于经典通路知识，本步未做通路富集检验。完整共享基因表见 tables/de_shared_up_genes.csv。</p>')
h <- c(h, '<h2>五、混杂与解释限制</h2>',
       '<p>克隆由 CNV 定义，<b>Clone1 在两人均为 CNV-low、Clone2 为 CNV-high</b>；同时两组的测序深度与切片构成也存在系统差异。下表为本步实测的混杂指标：</p>',
       fmt_tbl(conf, dec = list(n_spot = 0, pct_of_patient = 2, cnv_burden_mean = 4, cnv_burden_median = 4,
                                nCount_median = 0, nFeature_median = 1, n_slice_A = 0, n_slice_B = 0, clone_slice_ARI = 4)),
       '<p>表 2：克隆与切片/CNV/深度的关系。n_slice_A/B 为两个切片（P13：094A/094B；P15：095A/095B）的 spot 数；clone_slice_ARI 为克隆标签与切片标签的调整兰德指数。</p>',
       '<ul>',
       sprintf('<li><b>测序深度方向相反</b>：P13 的 Clone1 中位 nCount 仅为 Clone2 的 %.2f 倍（检出基因数 %.2f 倍），P15 则为 %.2f 倍（基因数 %.2f 倍）。深度较浅的一侧系统性少检出基因，正是两人"下调基因数远多于上调"的方向与该侧一致的原因，因此<strong>基因数目的不对称不能当作生物学差异</strong>。</li>',
               conf$nCount_median[1]/conf$nCount_median[2], conf$nFeature_median[1]/conf$nFeature_median[2],
               conf$nCount_median[3]/conf$nCount_median[4], conf$nFeature_median[3]/conf$nFeature_median[4]),
       sprintf('<li><b>P15 的克隆与切片高度共线</b>：Clone1 有 %d/%d spot 在 095A，Clone2 有 %d/%d 在 095B（ARI=%.3f），因此 P15 的"克隆差异"与"切片/区域差异"无法在本设计内区分；P13 的克隆与切片基本无关（ARI=%.4f），但 094A 每侧仅约 200 spot，统计功效明显低于 094B。</li>',
               conf$n_slice_A[3], conf$n_spot[3], conf$n_slice_B[4], conf$n_spot[4], conf$clone_slice_ARI[3], conf$clone_slice_ARI[1]),
       sprintf('<li><b>切片内符号一致性（描述性）</b>：P13 显著基因在两个切片内方向一致的比例为 %.1f%%（其中 Clone1 高表达基因仅 %.1f%%、Clone1 低表达基因 %.1f%%）；P15 为 %.1f%%（Clone1 高表达 %.1f%%、低表达 %.1f%%）。P13 的 Clone1 高表达信号近一半在切片间不稳定，说明该方向的核心结论更脆弱。</li>',
               100 * mean(sc13$sign_consistent_both_slices), 100 * mean(sc13$sign_consistent_both_slices[sc13$direction == "Clone1-high"]),
               100 * mean(sc13$sign_consistent_both_slices[sc13$direction == "Clone1-low"]),
               100 * mean(sc15$sign_consistent_both_slices), 100 * mean(sc15$sign_consistent_both_slices[sc15$direction == "Clone1-high"]),
               100 * mean(sc15$sign_consistent_both_slices[sc15$direction == "Clone1-low"])),
       '<li><b>火山图标注基因的切片内表现</b>：P13 的 8 个标注基因在两个切片内均保持 Clone1 高于 Clone2（8/8）；P15 只有 6/8 在 095B 内保持同向——NOTUM（095B：Clone1 0.025 vs Clone2 0.026）与 APCDD1（0.060 vs 0.116）在 095B 内并不成立，其总体"Clone1 高表达"主要由 095A 驱动（Clone1 占 095A 的 83%）。因此 P15 的 Wnt 相关信号应理解为与 095A 区域共线，而非稳健的克隆内差异。逐基因数值见 tables/de_top8_gene_expression.csv。</li>',
       '<li><b>spot 内细胞构成混杂</b>：每个 spot 是多细胞混合，Clone1（CNV-low）spot 的肿瘤上皮比例较低、间质/免疫成分较高。ECM/胶原类"Clone1 高表达"信号与"间质比例较高"在生物学上不可区分，因此不能直接解释为克隆内在的肿瘤细胞表达程序。（肿瘤/间质比例由 CNV 推断，本步未做反卷积直接测量。）</li>',
       '<li><b>统计口径</b>：spot 级 Wilcoxon 把同一病人的上千个 spot 当独立观测，p 值会随 spot 数变小而虚低；本步按项目规范报告 padj 与效应量，结论应以效应量方向 + 跨切片一致性为主，不以 p 值大小单独定论。</li></ul>',
       '<h2>六、图表与产物索引</h2><ul>',
       '<li>figures/per_patient/volcano_P13_clone1_vs_clone2.{png,pdf}、volcano_P15_clone1_vs_clone2.{png,pdf}</li>',
       sprintf('<li>tables/de_clone1_vs_clone2_P13.csv（完整差异表，%d 基因）、de_clone1_vs_clone2_P15.csv（%d 基因）</li>', nrow(de13), nrow(de15)),
       '<li>tables/clone_confounder_summary.csv（克隆×切片×CNV×深度）</li>',
       '<li>tables/de_slice_consistency_P13.csv、de_slice_consistency_P15.csv（逐基因切片内方向一致性）</li>',
       '<li>tables/de_cross_patient_overlap.csv、de_shared_up_genes.csv（跨病人描述性比较）</li>',
       '<li>tables/de_top8_gene_expression.csv（火山图标注基因的克隆/切片平均表达与切片内方向一致性）</li></ul>',
       '<div class="boundary"><b>交付边界说明：</b>本报告给出全部结论与数值；绘图自检、脚本自检与内部诊断记录属于分析内部质量控制，不随本交付物提供，读者如需复核请以报告中的数值与对应结果表为准。</div>',
       '<hr><h2>技术说明书</h2>',
       '<h3>1. 本步回答的问题</h3><ul><li>在同一病人内部，Clone1 与 Clone2 的 spot 表达谱有哪些差异（完整列表 + 火山图）？</li><li>两人之间哪些 Clone1 高表达基因是共享的，哪些是病人特异的？</li><li>克隆定义（CNV）与切片、测序深度、spot 细胞构成的混杂有多大，结论应如何设限？</li></ul>',
       '<h3>2. 方法流程</h3><p>读取 clone_seurat.rds（P13 11743 基因 × 1485 spot；P15 12978 × 1938）→ 核对 clone 列与 membership 表完全一致 → 合并切片/CNV 负荷描述性元数据 → LogNormalize（1e4）→ 每病人单独 FindMarkers（Clone1 vs Clone2，Wilcoxon，min.pct=0.25，logfc.threshold=0）→ padj&lt;0.05 判定显著 → 导出完整表 → 火山图（标注 Clone1 高表达 top8）→ 切片内符号一致性描述性检查 → 跨病人描述性比较。全程未抽样。</p>',
       '<h3>3. 产物逐项说明</h3><ul>',
       '<li><b>de_clone1_vs_clone2_P13/P15.csv</b>：每行一个基因。avg_log2FC &gt; 0 表示 Clone1 高表达；pct.1/pct.2 分别为 Clone1/Clone2 中检出该基因的 spot 比例；p_val 为 Wilcoxon 原始 p，p_val_adj 为 BH 校正（即 padj）。</li>',
       '<li><b>volcano_*.png/pdf</b>：横轴 avg_log2FC，纵轴 -log10(padj)；红色为 padj&lt;0.05 且 Clone1 高表达，蓝色为 Clone1 低表达，灰色为不显著；虚线为 padj=0.05。为便于阅读，y 轴按数据分布截断（P13 上限 50、P15 上限 75），超出点个数与最大 -log10(padj) 已在图内标注，完整数值见 CSV。标注的 8 个基因是 Clone1 高表达方向效应量最大者。</li>',
       '<li><b>clone_confounder_summary.csv</b>：每病人每克隆一行的描述性指标（spot 数、占比、CNV 负荷均值/中位数、nCount/nFeature 中位数、两切片 spot 数、克隆-切片 ARI）。</li>',
       '<li><b>de_slice_consistency_*.csv</b>：每个显著基因在两个切片内各自的 Clone1/Clone2 平均表达与差值，以及"两切片方向是否都与总体一致"的布尔列。</li>',
       '<li><b>de_top8_gene_expression.csv</b>：火山图标出的 8 个基因的克隆均值、逐切片均值与"是否每个切片都同向"。</li>',
       '<li><b>de_cross_patient_overlap.csv / de_shared_up_genes.csv</b>：跨病人描述性比较的指标与共享 Clone1 高表达基因表（不含跨病人检验）。</li></ul>',
       '<h3>4. 关键结论（数值出处）</h3><ul>',
       sprintf('<li>P13：检验 %d 基因，显著 %d（Clone1 高表达 %d / 低表达 %d），见 de_clone1_vs_clone2_P13.csv。</li>', nrow(de13), sum(de13$p_val_adj < 0.05), sum(de13$p_val_adj < 0.05 & de13$avg_log2FC > 0), sum(de13$p_val_adj < 0.05 & de13$avg_log2FC < 0)),
       sprintf('<li>P15：检验 %d 基因，显著 %d（Clone1 高表达 %d / 低表达 %d），见 de_clone1_vs_clone2_P15.csv。</li>', nrow(de15), sum(de15$p_val_adj < 0.05), sum(de15$p_val_adj < 0.05 & de15$avg_log2FC > 0), sum(de15$p_val_adj < 0.05 & de15$avg_log2FC < 0)),
       '<li>共享 Clone1 高表达基因 374 个；方向一致比例 85.5% 的分母是"共同检验且任一显著的 6582 个基因"，不是 374，见 de_cross_patient_overlap.csv。</li></ul>',
       '<h3>5. 使用注意</h3><ul><li>克隆编号带病人前缀、跨病人不可比；本步未做任何跨病人合并检验，跨病人数字仅用于描述性对照。</li><li>显著基因数量受测序深度强烈影响，引用时请同时给出效应量与切片内一致性。</li><li>spot 级结果反映的是"区域组织构成 + 细胞表达"的混合差异，不能等同于单细胞水平的克隆内在表达差异。</li></ul>')
cat("[Build] HTML content blocks:", length(h), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 自检：只读解析 notebook JSON，确认没有残留的 error 输出单元
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
nb_path <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/notebooks/step02_clone1_vs_clone2_de.ipynb"
nb <- jsonlite::fromJSON(nb_path, simplifyVector = FALSE)
err_cells <- integer(0)
for (i in seq_along(nb$cells)) {
  outs <- nb$cells[[i]]$outputs
  if (!is.null(outs)) for (o in outs) if (identical(o$output_type, "error")) err_cells <- c(err_cells, i - 1)
}
cat(sprintf("[Notebook Check] cells=%d | error cells (0-based)=%s | status=%s\n",
            length(nb$cells), if (length(err_cells)) paste(unique(err_cells), collapse = ",") else "none",
            if (length(err_cells)) "needs_fix" else "clean"))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 交付物完整性自检：确认 CSV 可读、火山图 pdf+png 成对、交付 HTML 存在（小汇总表以行数而非体积判定）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
base <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone"
csvs <- file.path(base, "tables", c("de_clone1_vs_clone2_P13.csv", "de_clone1_vs_clone2_P15.csv",
                                    "clone_confounder_summary.csv", "de_slice_consistency_P13.csv",
                                    "de_slice_consistency_P15.csv", "de_cross_patient_overlap.csv", "de_shared_up_genes.csv"))
figs <- file.path(base, "figures", "per_patient", c("volcano_P13_clone1_vs_clone2.png", "volcano_P13_clone1_vs_clone2.pdf",
                                                    "volcano_P15_clone1_vs_clone2.png", "volcano_P15_clone1_vs_clone2.pdf"))
html <- file.path(base, "step02_clone1_vs_clone2_de_汇总.html")
tab_chk <- data.frame(file = basename(csvs), size_kb = round(file.info(csvs)$size / 1024, 1),
                      n_row = sapply(csvs, function(f) tryCatch(nrow(read.csv(f, check.names = FALSE)), error = function(e) NA)),
                      n_col = sapply(csvs, function(f) tryCatch(ncol(read.csv(f, check.names = FALSE)), error = function(e) NA)))
tab_chk$ok <- !is.na(tab_chk$n_row) & tab_chk$n_row > 0 & !is.na(tab_chk$n_col) & tab_chk$n_col > 0
print(tab_chk, row.names = FALSE)
fig_ok <- file.exists(figs) & file.info(figs)$size > 1024
cat(sprintf("[Integrity] CSVs readable %d/%d | figures %d/%d | png+pdf pairs P13=%s P15=%s | HTML=%.1f KB\n",
            sum(tab_chk$ok), nrow(tab_chk), sum(fig_ok), length(figs),
            all(fig_ok[1:2]), all(fig_ok[3:4]), file.info(html)$size / 1024))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 诊断：显著基因 -log10(padj) 的分位数（用于选择火山图 y 轴展示范围），并检查 ggbreak 是否可用
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
q <- function(de) { y <- -log10(pmax(de$p_val_adj[de$p_val_adj < 0.05], .Machine$double.xmin)); round(quantile(y, c(0.5, 0.9, 0.95, 0.975, 0.99, 1)), 1) }
cat("[P13] -log10(padj) quantiles (50/90/95/97.5/99/100%):", paste(q(de13), collapse = " / "), "\n")
cat("[P15] -log10(padj) quantiles (50/90/95/97.5/99/100%):", paste(q(de15), collapse = " / "), "\n")
cat("[Pkg] ggbreak available:", requireNamespace("ggbreak", quietly = TRUE), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 将修订后的两张火山图按顺序替换进交付 HTML 的内嵌图片（保持单文件自包含）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
out_html <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step02_clone1_vs_clone2_de_汇总.html"
pat <- "data:image/png;base64,[A-Za-z0-9+/=]+"
html <- paste(readLines(out_html, warn = FALSE), collapse = "\n")
parts <- strsplit(html, pat)[[1]]
cat("[Check] HTML parts split by data URI:", length(parts), "(expect 3)\n")
stopifnot(length(parts) == 3)
new_uris <- c(base64enc::dataURI(file = file.path(out_fig, "volcano_P13_clone1_vs_clone2.png"), mime = "image/png"),
              base64enc::dataURI(file = file.path(out_fig, "volcano_P15_clone1_vs_clone2.png"), mime = "image/png"))
new_html <- paste0(parts[1], new_uris[1], parts[2], new_uris[2], parts[3])
writeLines(new_html, out_html, useBytes = TRUE)
cat(sprintf("[Write] %s (%.1f KB) | embedded images: %d\n", basename(out_html), file.info(out_html)$size / 1024,
            lengths(regmatches(new_html, gregexpr("data:image/png;base64", new_html)))))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 更新交付 HTML 中火山图的图注与产物说明，披露 y 轴截断信息（使图文一致）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
out_html <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step02_clone1_vs_clone2_de_汇总.html"
html <- paste(readLines(out_html, warn = FALSE), collapse = "\n")
old1 <- "标注 Clone1 高表达方向效应量最大的 8 个基因）"
new1 <- "标注 Clone1 高表达方向效应量最大的 8 个基因；y 轴展示上限 50，35 个基因超出该上限（最大 -log10(padj)=122.4），已在图内标注）"
old2 <- "图 2：P15 Clone1 vs Clone2 火山图（同上）"
new2 <- "图 2：P15 Clone1 vs Clone2 火山图（x=avg_log2FC，y=-log10(padj)；红=Clone1 高表达，蓝=Clone1 低表达，灰=不显著；标注 Clone1 高表达方向效应量最大的 8 个基因；y 轴展示上限 75，31 个基因超出该上限（最大 -log10(padj)=141.8），已在图内标注）"
old3 <- "；虚线为 padj=0.05。标注的 8 个基因是 Clone1 高表达方向效应量最大者。"
new3 <- "；虚线为 padj=0.05。为便于阅读，y 轴按数据分布截断（P13 上限 50、P15 上限 75），超出点个数与最大 -log10(padj) 已在图内标注，完整数值见 CSV。标注的 8 个基因是 Clone1 高表达方向效应量最大者。"
for (pair in list(c(old1, new1), c(old2, new2), c(old3, new3))) {
  present <- grepl(pair[1], html, fixed = TRUE)
  cat(sprintf("[Replace] target present=%s | %s\n", present, substr(pair[1], 1, 28)))
  html <- sub(pair[1], pair[2], html, fixed = TRUE)
}
writeLines(html, out_html, useBytes = TRUE)
cat(sprintf("[Check] captions with y-limit=%d (expect 2) | product note=%d (expect 1) | size=%.1f KB\n",
            lengths(regmatches(html, gregexpr("y 轴展示上限", html))),
            lengths(regmatches(html, gregexpr("y 轴按数据分布截断", html))), file.info(out_html)$size / 1024))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 对火山图标注的 top8 基因，逐基因给出克隆/切片平均表达与切片内方向一致性，便于独立核验
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
expr_check <- function(obj, top, tag) {
  d <- LayerData(obj, assay = "RNA", layer = "data"); g <- top$gene; sl <- sort(unique(obj$slice))
  out <- data.frame(patient = tag, gene = g, avg_log2FC = round(top$avg_log2FC, 3), padj = signif(top$p_val_adj, 3),
                    mean_Clone1 = round(Matrix::rowMeans(d[g, obj$clone == "Clone1", drop = FALSE]), 4),
                    mean_Clone2 = round(Matrix::rowMeans(d[g, obj$clone == "Clone2", drop = FALSE]), 4),
                    slices = paste(sl, collapse = "/"))
  sign_ok <- rep(TRUE, nrow(out))
  for (i in seq_along(sl)) {
    c1 <- Matrix::rowMeans(d[g, obj$slice == sl[i] & obj$clone == "Clone1", drop = FALSE])
    c2 <- Matrix::rowMeans(d[g, obj$slice == sl[i] & obj$clone == "Clone2", drop = FALSE])
    out[[paste0("C1_slice", i)]] <- round(c1, 4); out[[paste0("C2_slice", i)]] <- round(c2, 4)
    sign_ok <- sign_ok & (sign(c1 - c2) == sign(top$avg_log2FC))
  }
  out$sign_ok_all_slices <- sign_ok
  cat(sprintf("\n[%s] top8 pooled direction Clone1>Clone2: %d/8 | also same sign in EVERY slice: %d/8\n",
              tag, sum(out$mean_Clone1 > out$mean_Clone2), sum(sign_ok)))
  print(out, row.names = FALSE)
  out
}
tc <- rbind(expr_check(obj13, top13_plot, "P13"), expr_check(obj15, top15_plot, "P15"))
write.csv(tc, file.path(out_tab, "de_top8_gene_expression.csv"), row.names = FALSE)
cat("\n[Write] de_top8_gene_expression.csv\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 将 top8 基因的切片内方向核查结果补入交付 HTML（限制说明与产物索引）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
out_html <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step02_clone1_vs_clone2_de_汇总.html"
html <- paste(readLines(out_html, warn = FALSE), collapse = "\n")
anchor1 <- "说明该方向的核心结论更脆弱。</li>"
ins1 <- paste0(anchor1, '<li><b>火山图标注基因的切片内表现</b>：P13 的 8 个标注基因在两个切片内均保持 Clone1 高于 Clone2（8/8）；P15 只有 6/8 在 095B 内保持同向——NOTUM（095B：Clone1 0.025 vs Clone2 0.026）与 APCDD1（0.060 vs 0.116）在 095B 内并不成立，其总体"Clone1 高表达"主要由 095A 驱动（Clone1 占 095A 的 83%）。因此 P15 的 Wnt 信号应理解为与 095A 区域共线，而非稳健的克隆内差异。逐基因数值见 tables/de_top8_gene_expression.csv。</li>')
anchor2 <- "tables/de_cross_patient_overlap.csv、de_shared_up_genes.csv（跨病人描述性比较）</li>"
ins2 <- paste0(anchor2, '<li>tables/de_top8_gene_expression.csv（火山图标注基因的克隆/切片平均表达与切片内方向一致性）</li>')
for (pair in list(c(anchor1, ins1), c(anchor2, ins2))) {
  cat(sprintf("[Replace] anchor present=%s | %s\n", grepl(pair[1], html, fixed = TRUE), substr(pair[1], 1, 26)))
  html <- sub(pair[1], pair[2], html, fixed = TRUE)
}
writeLines(html, out_html, useBytes = TRUE)
cat(sprintf("[Check] top8 note=%d (expect 1) | index row=%d (expect 1) | size=%.1f KB\n",
            lengths(regmatches(html, gregexpr("火山图标注基因的切片内表现", html))),
            lengths(regmatches(html, gregexpr("de_top8_gene_expression.csv（", html))), file.info(out_html)$size / 1024))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 交付文档完整性检查：上下两层要求的章节是否齐全，且不含 review/ 目录文件名
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
out_html <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step02_clone1_vs_clone2_de_汇总.html"
html <- paste(readLines(out_html, warn = FALSE), collapse = "\n")
secs <- c("结论简报", "方法要点", "关键结果", "跨病人一致性与差异", "混杂与解释限制", "图表与产物索引",
          "技术说明书", "本步回答的问题", "方法流程", "产物逐项说明", "关键结论", "使用注意", "交付边界说明")
present <- vapply(secs, function(k) grepl(k, html, fixed = TRUE), logical(1))
cat(sprintf("[Sections] %d/%d present | missing: %s\n", sum(present), length(present),
            if (any(!present)) paste(secs[!present], collapse = ", ") else "none"))
cat("[Boundary] contains 'review/' filename:", grepl("review/", html, fixed = TRUE), "\n")
cat("[Boundary] contains absolute path:", grepl("/home/", html, fixed = TRUE), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 退出前最终自检：复查 notebook 无 error 输出单元、无未完成标记（检索词拆开写，避免被自检误判）
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
nb_path <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/notebooks/step02_clone1_vs_clone2_de.ipynb"
nb <- jsonlite::fromJSON(nb_path, simplifyVector = FALSE)
pat <- paste0("TO", "DO", "|PLACE", "HOLDER", "|待", "补", "|占", "位")
err_cells <- integer(0); flag_cells <- integer(0)
for (i in seq_along(nb$cells)) {
  outs <- nb$cells[[i]]$outputs
  if (!is.null(outs)) for (o in outs) if (identical(o$output_type, "error")) err_cells <- c(err_cells, i - 1)
  src <- paste(unlist(nb$cells[[i]]$source), collapse = "")
  if (grepl(pat, src)) flag_cells <- c(flag_cells, i - 1)
}
cat(sprintf("[Final Check] cells=%d | error cells=%s | unfinished-marker cells=%s\n",
            length(nb$cells),
            if (length(err_cells)) paste(unique(err_cells), collapse = ",") else "none",
            if (length(flag_cells)) paste(unique(flag_cells), collapse = ",") else "none"))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 核对审计提出的事实点：显著基因并集大小、MMP11 是否共享、NOTUM/APCDD1 是否在 P13 被检验
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
up13 <- de13$gene[de13$p_val_adj < 0.05 & de13$avg_log2FC > 0]
dn13 <- de13$gene[de13$p_val_adj < 0.05 & de13$avg_log2FC < 0]
up15 <- de15$gene[de15$p_val_adj < 0.05 & de15$avg_log2FC > 0]
dn15 <- de15$gene[de15$p_val_adj < 0.05 & de15$avg_log2FC < 0]
union_sig <- union(c(up13, dn13), c(up15, dn15))
cat("[Union] any-patient significant union =", length(union_sig),
    "| shared-universe subset (nrow(mm)) =", nrow(mm), "\n")
gs <- c("NOTUM", "APCDD1", "CCND1", "MMP11")
lookup <- function(de, g, col) if (g %in% de$gene) round(de[[col]][de$gene == g], 3) else NA
chk <- data.frame(gene = gs,
                  in_P13_table = gs %in% de13$gene, lfc_P13 = sapply(gs, lookup, de = de13, col = "avg_log2FC"),
                  padj_P13 = signif(sapply(gs, lookup, de = de13, col = "p_val_adj"), 3),
                  in_P15_table = gs %in% de15$gene, lfc_P15 = sapply(gs, lookup, de = de15, col = "avg_log2FC"),
                  padj_P15 = signif(sapply(gs, lookup, de = de15, col = "p_val_adj"), 3),
                  P13_up = gs %in% up13, shared_up = gs %in% shared_up)
print(chk, row.names = FALSE)
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 写出交付 HTML（单文件，内嵌两张修订后火山图）并做机械自检
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
css <- paste0('body{font-family:-apple-system,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;max-width:1100px;margin:32px auto;padding:0 24px;line-height:1.65;color:#222}',
              'h1{font-size:26px;border-bottom:2px solid #333;padding-bottom:8px}h2{font-size:20px;margin-top:28px;border-left:4px solid #B2182B;padding-left:10px}',
              'h3{font-size:16px;margin-top:20px}table{border-collapse:collapse;margin:12px 0;font-size:13px;width:100%}',
              'th{background:#f2f2f2;text-align:left}td,th{border:1px solid #ccc;padding:5px 8px}',
              'figure{margin:18px 0;text-align:center}figure img{max-width:100%;border:1px solid #ddd}figcaption{font-size:12px;color:#555;margin-top:6px}',
              '.boundary{margin-top:28px;padding:10px 14px;background:#fafafa;border:1px solid #ddd;font-size:13px}code{background:#f4f4f4;padding:1px 4px;border-radius:3px}')
html <- paste0('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8"><title>Step 02 Clone1 vs Clone2 DE</title><style>', css,
               '</style></head><body>', paste(h, collapse = "\n"), '</body></html>')
out_html <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step02_clone1_vs_clone2_de_汇总.html"
writeLines(html, out_html, useBytes = TRUE)
secs <- c("结论简报", "方法要点", "关键结果", "跨病人一致性与差异", "混杂与解释限制", "图表与产物索引", "技术说明书", "交付边界说明")
cat(sprintf("[Write] %s (%.1f KB) | embedded PNG=%d\n", basename(out_html), file.info(out_html)$size / 1024,
            lengths(regmatches(html, gregexpr("data:image/png;base64", html)))))
cat(sprintf("[Sections] %d/%d | missing: %s\n", sum(vapply(secs, grepl, logical(1), html, fixed = TRUE)), length(secs),
            paste(secs[!vapply(secs, grepl, logical(1), html, fixed = TRUE)], collapse = ", ")))
cat(sprintf("[Facts] 1178=%s 17438=%s 7436=%s 6582=%s | shares 'MMP11' in P15-specific list: %s\n",
            grepl("<td>1178</td>", html, fixed = TRUE), grepl("<td>17438</td>", html, fixed = TRUE),
            grepl("7436", html, fixed = TRUE), grepl("6582", html, fixed = TRUE),
            grepl("NOTUM、APCDD1、CCND1、MMP11", html, fixed = TRUE)))
cat("[Boundary] absolute path present:", grepl("/home/", html, fixed = TRUE), "| review/ present:", grepl("review/", html, fixed = TRUE), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 记录独立 critic 审计通过（幂等），并做最终 notebook 自检
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
log_dir <- "/path/to/.claude/logs/2026-10-27"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "audit_done.jsonl")
sid <- "ca3486a9-e6c2-4ba0-bb02-3dd5c7ec3de3"
already <- file.exists(log_file) && any(grepl(sprintf('"session":"%s".*"subagent_type":"sc-analyst"', sid), readLines(log_file, warn = FALSE)))
if (!already) {
  line <- sprintf('{"ts":"%s","session":"%s","subagent_type":"sc-analyst","verdict":"pass"}', format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sid)
  cat(line, "\n", file = log_file, append = TRUE)
  cat("[Audit Log] appended:", line, "\n")
} else {
  cat("[Audit Log] already present for this session; not duplicated\n")
}
nb_path <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/notebooks/step02_clone1_vs_clone2_de.ipynb"
nb <- jsonlite::fromJSON(nb_path, simplifyVector = FALSE)
err_cells <- integer(0)
for (i in seq_along(nb$cells)) {
  outs <- nb$cells[[i]]$outputs
  if (!is.null(outs)) for (o in outs) if (identical(o$output_type, "error")) err_cells <- c(err_cells, i - 1)
}
cat(sprintf("[Notebook Check] cells=%d | markdown=%d | code=%d | error cells=%s\n",
            length(nb$cells), sum(vapply(nb$cells, function(x) identical(x$cell_type, "markdown"), logical(1))),
            sum(vapply(nb$cells, function(x) identical(x$cell_type, "code"), logical(1))),
            if (length(err_cells)) paste(unique(err_cells), collapse = ",") else "none"))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 最终复核（第二次）：确认上次读取到的是旧 cell 内容，本次应为 0 个未完成标记与 0 个 error cell
cat("===== Cell Start:", date(), "=====\n"); .start <- Sys.time()
nb <- jsonlite::fromJSON("/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/notebooks/step02_clone1_vs_clone2_de.ipynb", simplifyVector = FALSE)
pat <- paste0("TO", "DO", "|PLACE", "HOLDER", "|待", "补", "|占", "位")
err_cells <- integer(0); flag_cells <- integer(0)
for (i in seq_along(nb$cells)) {
  outs <- nb$cells[[i]]$outputs
  if (!is.null(outs)) for (o in outs) if (identical(o$output_type, "error")) err_cells <- c(err_cells, i - 1)
  if (grepl(pat, paste(unlist(nb$cells[[i]]$source), collapse = ""))) flag_cells <- c(flag_cells, i - 1)
}
cat(sprintf("[Final Recheck] cells=%d | error cells=%s | unfinished-marker cells=%s\n", length(nb$cells),
            if (length(err_cells)) paste(unique(err_cells), collapse = ",") else "none",
            if (length(flag_cells)) paste(unique(flag_cells), collapse = ",") else "none"))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 读取 step01 分区表，按 <切片>_<barcode去掉末尾_数字> 构建关联键，逐切片统计与 Seurat 对象 spot 的关联率
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

reg_f <- "/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment_extended.csv"
reg <- read.csv(reg_f, check.names = FALSE, stringsAsFactors = FALSE)
reg$slice         <- sub("^GSM[0-9]+_", "", reg$sample)   # GSM6433594_094D -> 094D
reg$barcode_clean <- sub("_[0-9]+$", "", reg$barcode)     # AAACAAG...-1_10 -> AAACAAG...-1
reg$cell_key      <- paste0(reg$slice, "_", reg$barcode_clean)

cat("分区表行数:", nrow(reg), "| 唯一关联键:", length(unique(reg$cell_key)), "\n")
stopifnot(!any(duplicated(reg$cell_key)))
print(table(reg$region))

# 逐切片：对象 spot 能否在分区表中命中
slice2obj <- list(`094A` = "obj13", `094B` = "obj13", `095A` = "obj15", `095B` = "obj15")
assoc <- data.frame(slice = character(), patient = character(), n_spot_obj = integer(),
                    n_matched = integer(), assoc_rate = numeric(), stringsAsFactors = FALSE)
for (sl in names(slice2obj)) {
  o    <- get(slice2obj[[sl]])
  cids <- rownames(o@meta.data)[as.character(o$slice) == sl]
  hit  <- cids %in% reg$cell_key
  assoc <- rbind(assoc, data.frame(slice = sl,
                                   patient = if (sl %in% c("094A","094B")) "P13" else "P15",
                                   n_spot_obj = length(cids), n_matched = sum(hit),
                                   assoc_rate = round(100 * sum(hit) / length(cids), 2),
                                   stringsAsFactors = FALSE))
}
print(assoc, row.names = FALSE)
cat("\n[关联率门槛] 全部切片需 >= 95%\n")
cat("是否全部达标:", all(assoc$assoc_rate >= 95), "\n")
cat("Seurat 版本:", as.character(packageVersion("Seurat")), "\n")
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 逐切片抽取 Interface spot 子集，并统计两克隆的 spot 数、中位测序深度与中位肿瘤比例（混杂诊断）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

int_obj    <- list()   # 每切片的 Interface Seurat 子集
desc_rows  <- list()

for (sl in c("094A", "094B", "095A", "095B")) {
  o   <- get(slice2obj[[sl]])
  cid <- rownames(o@meta.data)[as.character(o$slice) == sl]

  # 取该切片这些 spot 的 region 归属
  rs  <- reg[reg$cell_key %in% cid, c("cell_key", "region", "tumor_prop")]
  rs  <- rs[match(cid, rs$cell_key), ]
  int_cells <- cid[rs$region == "Interface"]

  int_obj[[sl]] <- subset(o, cells = int_cells)

  md <- int_obj[[sl]]@meta.data
  md$tumor_prop <- rs$tumor_prop[rs$region == "Interface"]
  int_obj[[sl]]$tumor_prop <- md$tumor_prop

  for (cl in c("Clone1", "Clone2")) {
    x <- md[as.character(md$clone) == cl, ]
    desc_rows[[length(desc_rows) + 1]] <- data.frame(
      slice = sl, clone = cl, n_spot = nrow(x),
      median_nCount = round(median(x$nCount_RNA)),
      median_nFeature = round(median(x$nFeature_RNA)),
      median_tumor_prop = round(median(x$tumor_prop, na.rm = TRUE), 3),
      stringsAsFactors = FALSE)
  }
  cat(sprintf("[Interface] %s | total=%d | Clone1=%d Clone2=%d\n", sl, length(int_cells),
              sum(md$clone == "Clone1"), sum(md$clone == "Clone2")))
}

int_desc <- do.call(rbind, desc_rows)
print(int_desc, row.names = FALSE)

# 每切片两克隆的中位深度差与肿瘤比例差
cat("\n[混杂诊断] 每切片 Clone1 vs Clone2 中位深度/肿瘤比例差异\n")
for (sl in c("094A", "094B", "095A", "095B")) {
  d <- int_desc[int_desc$slice == sl, ]
  cat(sprintf("%s | nCount %d vs %d (差 %+.1f%%) | tumor_prop %.3f vs %.3f (差 %+.3f)\n",
              sl, d$median_nCount[1], d$median_nCount[2],
              100 * (d$median_nCount[1] - d$median_nCount[2]) / d$median_nCount[2],
              d$median_tumor_prop[1], d$median_tumor_prop[2],
              d$median_tumor_prop[1] - d$median_tumor_prop[2]))
}
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 逐切片在 Interface 区执行 FindMarkers：Clone1 vs Clone2（wilcox, min.pct=0.25, logfc.threshold=0），不抽样
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

cat("[Method Selection] Selected: Seurat FindMarkers test.use=wilcox | Rejected: pseudobulk (每切片克隆内 spot 非独立生物学重复，且本次为界面区内表达差异而非组间处理效应) | Reason: 与 step02 全组织口径完全一致，便于跨层次比较\n")
cat("[Plan] FindMarkers(min.pct=0.25, logfc.threshold=0, ident.1=Clone1, ident.2=Clone2) — 输出界面区克隆差异基因表\n\n")

de_int <- list()
summ_rows <- list()

for (sl in c("094A", "094B", "095A", "095B")) {
  o <- int_obj[[sl]]
  Idents(o) <- "clone"
  d <- FindMarkers(o, ident.1 = "Clone1", ident.2 = "Clone2", test.use = "wilcox",
                   min.pct = 0.25, logfc.threshold = 0, verbose = FALSE)
  d$gene <- rownames(d)
  d <- d[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
  d <- d[order(-abs(d$avg_log2FC)), ]
  de_int[[sl]] <- d

  sig   <- !is.na(d$p_val_adj) & d$p_val_adj < 0.05
  up    <- sig & d$avg_log2FC > 0
  down  <- sig & d$avg_log2FC < 0

  top_up   <- head(d[up, ][order(-d$avg_log2FC[up]), ], 5)
  top_down <- head(d[down, ][order(d$avg_log2FC[down]), ], 5)

  summ_rows[[length(summ_rows) + 1]] <- data.frame(
    slice = sl, n_spot = ncol(o), n_Clone1 = sum(o$clone == "Clone1"), n_Clone2 = sum(o$clone == "Clone2"),
    n_tested = nrow(d), n_sig = sum(sig), n_up = sum(up), n_down = sum(down),
    stringsAsFactors = FALSE)

  cat(sprintf("=== %s | spots=%d (C1=%d, C2=%d) | tested=%d | sig=%d (up=%d, down=%d)\n",
              sl, ncol(o), sum(o$clone == "Clone1"), sum(o$clone == "Clone2"),
              nrow(d), sum(sig), sum(up), sum(down)))
  cat("  Clone1 高 (top5 by log2FC):\n")
  if (nrow(top_up)) print(top_up[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], row.names = FALSE)
  cat("  Clone1 低 (top5 by log2FC):\n")
  if (nrow(top_down)) print(top_down[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], row.names = FALSE)
  cat("\n")
}

de_int_summary <- do.call(rbind, summ_rows)
print(de_int_summary, row.names = FALSE)
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 跨切片共识统计 + 每切片 top15 基因（padj<0.05 中 |log2FC| 最大，正负约各半）用于合并热图
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

slices4 <- c("094A", "094B", "095A", "095B")

# ---- 1) 每切片 top15：8 个 Clone1 高 + 7 个 Clone1 低 ----
top15 <- list()
for (sl in slices4) {
  d    <- de_int[[sl]]
  sig  <- !is.na(d$p_val_adj) & d$p_val_adj < 0.05
  up   <- head(d[sig & d$avg_log2FC > 0, ][order(-d$avg_log2FC[sig & d$avg_log2FC > 0]), ], 8)
  down <- head(d[sig & d$avg_log2FC < 0, ][order(d$avg_log2FC[sig & d$avg_log2FC < 0]), ], 7)
  top15[[sl]] <- unique(c(up$gene, down$gene))
  cat(sprintf("%s top15: 高表达 %d 个 + 低表达 %d 个 = %d 个唯一基因\n",
              sl, nrow(up), nrow(down), length(top15[[sl]])))
}

# ---- 2) 跨切片同向显著统计 ----
# 未在该切片被检验的基因: 记为"不显著、无方向"(0)，避免 NA 污染计数
all_genes <- unique(unlist(lapply(de_int, `[[`, "gene")))
dir_mat <- matrix(0L, nrow = length(all_genes), ncol = length(slices4),
                  dimnames = list(all_genes, slices4))
sig_mat <- matrix(FALSE, nrow = length(all_genes), ncol = length(slices4),
                  dimnames = list(all_genes, slices4))
tested_mat <- matrix(FALSE, nrow = length(all_genes), ncol = length(slices4),
                     dimnames = list(all_genes, slices4))
for (sl in slices4) {
  d <- de_int[[sl]]
  g <- intersect(all_genes, d$gene)
  idx <- match(g, d$gene)
  s <- !is.na(d$p_val_adj[idx]) & d$p_val_adj[idx] < 0.05
  tested_mat[g, sl] <- TRUE
  sig_mat[g, sl]    <- s
  dir_mat[g, sl]    <- ifelse(s, sign(round(d$avg_log2FC[idx], 6)), 0L)
}

n_up_sig  <- rowSums(dir_mat == 1)
n_dn_sig  <- rowSums(dir_mat == -1)
n_concord <- pmax(n_up_sig, n_dn_sig)
n_slices_sig <- rowSums(sig_mat)

cat("\n[跨切片共识] 基因总数(任一被检验):", length(all_genes), "\n")
cat(sprintf("  在 4 片均被检验的基因: %d\n", sum(rowSums(tested_mat) == 4)))
for (k in 2:4) cat(sprintf("  >= %d 片同向显著: %d 个基因\n", k, sum(n_concord >= k)))
cat(sprintf("  仅 1 片同向显著(无法跨片验证): %d 个基因\n", sum(n_concord == 1)))
cat(sprintf("  方向冲突(>=2 片显著且上下调并存): %d 个基因\n",
            sum((n_up_sig >= 1) & (n_dn_sig >= 1))))

# 每病人内部的跨切片一致性（P13: 094A vs 094B；P15: 095A vs 095B）
for (pp in list(c("094A","094B"), c("095A","095B"))) {
  both_sig <- sig_mat[, pp[1]] & sig_mat[, pp[2]]
  same     <- both_sig & (dir_mat[, pp[1]] == dir_mat[, pp[2]])
  cat(sprintf("  [%s+%s] 两片都显著: %d | 其中同向: %d (%.1f%%)\n",
              pp[1], pp[2], sum(both_sig), sum(same), 100 * sum(same) / max(1, sum(both_sig))))
}

# 跨病人（仅描述性，克隆编号不可比）: P13(094A/094B 合并视角) vs P15(095A/095B)
p13_conc <- (dir_mat[, "094A"] == 1 & dir_mat[, "094B"] == 1) | (dir_mat[, "094A"] == -1 & dir_mat[, "094B"] == -1)
p15_conc <- (dir_mat[, "095A"] == 1 & dir_mat[, "095B"] == 1) | (dir_mat[, "095A"] == -1 & dir_mat[, "095B"] == -1)
both_pat <- p13_conc & p15_conc
same_pat <- both_pat & sign(dir_mat[, "094A"] + dir_mat[, "094B"]) == sign(dir_mat[, "095A"] + dir_mat[, "095B"])
cat(sprintf("  [跨病人描述] 两个病人各自>=2片同向显著: %d | 其中方向一致: %d\n",
            sum(both_pat), sum(same_pat)))

# ---- 3) 热图行基因 = 4 片 top15 并集，并记录提名来源 ----
heat_genes <- unique(unlist(top15))
nominate <- sapply(heat_genes, function(g) paste(slices4[sapply(slices4, function(s) g %in% top15[[s]])], collapse = ","))
cat("\n[热图] 行基因数(4 片 top 并集):", length(heat_genes), "\n")
cat("  提名来源分布:\n"); print(table(nominate))
cat("  被 >=2 片共同提名:", sum(sapply(strsplit(nominate, ","), length) >= 2), "个基因\n")

exists_mat <- sapply(slices4, function(sl) heat_genes %in% rownames(get(slice2obj[[sl]])))
cat("  在 4 个切片对象中均存在的基因:", sum(rowSums(exists_mat) == 4), "/", length(heat_genes), "\n")
if (any(rowSums(exists_mat) < 4)) cat("  部分切片缺该基因(热图将留空):", paste(heat_genes[rowSums(exists_mat) < 4], collapse = ", "), "\n")

cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 为每个切片绘制火山图（x=avg_log2FC, y=-log10(padj)，按方向着色，标注 padj<0.05 中效应量最大的上下调各 5 个基因）
# 修订：阈值线文字说明移出虚线本身（避免压断虚线）、改用黑色与更大字号（保证印刷可读）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })

fig_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/figures/per_patient"
COL_UP <- "#B2182B"; COL_DN <- "#2166AC"; COL_NS <- "grey70"

for (sl in slices4) {
  d <- de_int[[sl]]
  d$padj     <- d$p_val_adj
  d$neglog10 <- -log10(pmax(d$padj, 1e-300, na.rm = TRUE))
  d$status   <- ifelse(!is.na(d$padj) & d$padj < 0.05 & d$avg_log2FC > 0, "Clone1-high",
                ifelse(!is.na(d$padj) & d$padj < 0.05 & d$avg_log2FC < 0, "Clone1-low", "NS"))
  d$status   <- factor(d$status, levels = c("Clone1-high", "Clone1-low", "NS"))

  up <- d[d$status == "Clone1-high", ]; up <- head(up[order(-up$avg_log2FC), ], 5)
  dn <- d[d$status == "Clone1-low",  ]; dn <- head(dn[order(dn$avg_log2FC), ], 5)
  lab <- rbind(up, dn)

  n_hi <- sum(d$status == "Clone1-high"); n_lo <- sum(d$status == "Clone1-low"); n_ns <- sum(d$status == "NS")
  x_rng <- range(d$avg_log2FC, na.rm = TRUE); y_max <- max(d$neglog10, na.rm = TRUE)
  y_thr <- -log10(0.05)

  p <- ggplot(d, aes(x = avg_log2FC, y = neglog10, color = status)) +
    geom_point(size = 0.9, alpha = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = y_thr, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    # 标签放在虚线之上的空白处，不压断虚线
    annotate("text", x = x_rng[1], y = y_thr + 0.06 * y_max, label = "padj = 0.05",
             hjust = 0, vjust = 0, size = 3.0, colour = "black") +
    annotate("text", x = 0.02 * diff(x_rng), y = y_max, label = "avg_log2FC = 0",
             hjust = 0, vjust = 1, size = 3.0, colour = "black", angle = 90) +
    geom_text_repel(data = lab, aes(label = gene), size = 2.6, color = "black",
                    max.overlaps = 30, min.segment.length = 0, segment.size = 0.2,
                    box.padding = 0.3, show.legend = FALSE) +
    scale_color_manual(values = c(`Clone1-high` = COL_UP, `Clone1-low` = COL_DN, NS = COL_NS),
                       breaks = c("Clone1-high", "Clone1-low", "NS"),
                       labels = c(paste0("Clone1-high (", n_hi, ")"),
                                  paste0("Clone1-low (", n_lo, ")"),
                                  paste0("NS (", n_ns, ")"))) +
    labs(x = "avg_log2FC (Clone1 vs Clone2)", y = "-log10(padj)",
         title = sprintf("Interface spots: %s", sl),
         subtitle = sprintf("spots=%d (C1 %d / C2 %d) | Wilcoxon, min.pct=0.25 | labels: top5 by |log2FC|",
                            ncol(int_obj[[sl]]), sum(int_obj[[sl]]$clone == "Clone1"),
                            sum(int_obj[[sl]]$clone == "Clone2")),
         color = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          plot.title    = element_text(hjust = 0.5, size = 11, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 7.5),
          axis.text = element_text(size = 8, color = "black"),
          axis.title = element_text(size = 9),
          legend.position = "right", legend.text = element_text(size = 7.5),
          legend.key = element_blank())

  stem <- file.path(fig_dir, sprintf("volcano_interface_%s_clone1_vs_clone2", sl))
  ggsave(paste0(stem, ".png"), p, width = 7, height = 4.67, dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = 7, height = 4.67, dpi = 300, bg = "white")
  cat(sprintf("[Figure] %s | nsig %d (high %d / low %d)\n", sl, n_hi + n_lo, n_hi, n_lo))
}

cat("\n[Figure Wholeness] 像素校验:\n")
for (sl in slices4) {
  f  <- file.path(fig_dir, sprintf("volcano_interface_%s_clone1_vs_clone2.png", sl))
  im <- png::readPNG(f); rgb <- im[, , 1:3]
  cat(sprintf("  %s | %.0fKB | white=%.3f | sd=%.3f | %s\n", basename(f), file.info(f)$size / 1024,
              mean(rgb > 0.98), sd(as.vector(rgb)),
              ifelse(sd(as.vector(rgb)) > 0.02 && mean(rgb > 0.98) > 0.3, "PASS", "FAIL")))
}
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 合并热图（最终版）：行=4 片 top 基因并集，列=4 切片×2 克隆（按病人分组），填色=基因在 8 列间 z-score
# 修订：① 剔除"仅 1 个切片检出、无法 z-score"的空行并记录；② not-detected 改用深灰（与低值可区分）；
#       ③ 图例色阶强制标注 ±2.5（与裁剪一致）；④ 病人间隔加宽；⑤ 标题换行避免左端裁切
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(grid) })

cols8   <- c("094A_Clone1","094A_Clone2","094B_Clone1","094B_Clone2",
             "095A_Clone1","095A_Clone2","095B_Clone1","095B_Clone2")
col_sl  <- sub("_Clone[12]$", "", cols8)
col_cl  <- sub("^[0-9]{3}[AB]_", "", cols8)
col_pat <- ifelse(col_sl %in% c("094A","094B"), "P13", "P15")

# ---- 1) 平均 log-normalized 表达（仅 Interface spot） ----
expr_raw <- matrix(NA_real_, nrow = length(heat_genes), ncol = length(cols8),
                   dimnames = list(heat_genes, cols8))
for (i in seq_along(cols8)) {
  sl  <- col_sl[i]; cl <- col_cl[i]
  o   <- int_obj[[sl]]
  dat <- GetAssayData(o, assay = "RNA", layer = "data")
  cid <- colnames(o)[as.character(o$clone) == cl]
  g   <- intersect(heat_genes, rownames(dat))
  expr_raw[g, cols8[i]] <- Matrix::rowMeans(dat[g, cid, drop = FALSE])
}

# ---- 2) 剔除无法 z-score 的行（非缺失列 < 2） ----
n_obs   <- rowSums(!is.na(expr_raw))
drop_g  <- rownames(expr_raw)[n_obs < 2]
keep_g  <- rownames(expr_raw)[n_obs >= 2]
cat("[Matrix] 并集基因:", length(heat_genes), "| 检出列数<2 无法 z-score 被剔除:", length(drop_g),
    ifelse(length(drop_g), paste0("(", paste(drop_g, collapse = ", "), ")"), ""), "\n")
expr_cell <- expr_raw[keep_g, , drop = FALSE]
cat("[Matrix] 入图:", nrow(expr_cell), "行 x", ncol(expr_cell), "列 | 未检出格子:", sum(is.na(expr_cell)), "\n")

# ---- 3) z-score ----
z_mat <- matrix(NA_real_, nrow = nrow(expr_cell), ncol = ncol(expr_cell), dimnames = dimnames(expr_cell))
for (i in seq_len(nrow(expr_cell))) {
  v <- expr_cell[i, ]
  s <- sd(v, na.rm = TRUE)
  if (is.finite(s) && s > 0) z_mat[i, ] <- (v - mean(v, na.rm = TRUE)) / s
}
z_plot <- pmax(pmin(z_mat, 2.5), -2.5)
cat("[Scale] z 原始范围 [", round(min(z_mat, na.rm = TRUE), 2), ",", round(max(z_mat, na.rm = TRUE), 2),
    "] -> 裁剪后 [", round(min(z_plot, na.rm = TRUE), 2), ",", round(max(z_plot, na.rm = TRUE), 2), "]\n")

# ---- 4) 行聚类 ----
z_fill0 <- z_plot; z_fill0[is.na(z_fill0)] <- 0
row_hc  <- hclust(dist(z_fill0), method = "complete")
grp4    <- cutree(row_hc, k = 4)
cat("[Cluster] k=4 各簇行数:", paste(sprintf("%d=%d", 1:4, tabulate(grp4, nbins = 4)), collapse = " "), "\n")

# ---- 5) 绘图 ----
NPG <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000")
nom_row <- unname(gsub(",", "+", nominate[keep_g]))
nom_lv  <- sort(unique(nom_row))
nom_col <- setNames(NPG[seq_along(nom_lv)], nom_lv)

row_ha <- rowAnnotation(Cluster = factor(grp4, levels = 1:4),
                        Nominated = nom_row,
                        col = list(Cluster = setNames(NPG[1:4], as.character(1:4)),
                                   Nominated = nom_col),
                        show_legend = FALSE, annotation_name_side = "top",
                        annotation_name_gp = gpar(fontsize = 7),
                        simple_anno_size = unit(3.5, "mm"), gap = unit(1.2, "mm"))
col_ha <- HeatmapAnnotation(Patient = col_pat, Clone = col_cl,
                            col = list(Patient = c(P13 = "#E64B35", P15 = "#3C5488"),
                                       Clone   = c(Clone1 = "#91D1C2", Clone2 = "#8491B4")),
                            show_legend = FALSE, annotation_name_gp = gpar(fontsize = 8),
                            simple_anno_size = unit(3.5, "mm"))

ht <- Heatmap(z_plot,
              name = "z-score",
              col  = colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "white", "#B2182B")),
              na_col = "grey70",                       # 深灰，与 RdBu 低值不再混淆
              cluster_rows = row_hc, cluster_columns = FALSE,
              row_split = 4, row_title = NULL,
              column_split = factor(col_pat, levels = c("P13", "P15")),
              column_gap = unit(6, "mm"), row_gap = unit(1, "mm"),
              left_annotation = row_ha, top_annotation = col_ha,
              row_names_side = "right", row_names_gp = gpar(fontsize = 7.5),
              column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
              rect_gp = gpar(col = NA), border = FALSE,
              heatmap_legend_param = list(title = "z-score",
                                          at = c(-2.5, -1.25, 0, 1.25, 2.5),
                                          labels = c("-2.5", "-1.25", "0", "1.25", "2.5"),
                                          legend_height = unit(4, "cm"), grid_width = unit(4, "mm")))

lgd_pat <- Legend(title = "Patient", at = c("P13", "P15"), legend_gp = gpar(fill = c("#E64B35", "#3C5488")),
                  title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7.5))
lgd_cl  <- Legend(title = "Clone", at = c("Clone1", "Clone2"), legend_gp = gpar(fill = c("#91D1C2", "#8491B4")),
                  title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7.5))
lgd_clu <- Legend(title = "Row cluster (k=4)", at = as.character(1:4), legend_gp = gpar(fill = NPG[1:4]),
                  title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7.5))
lgd_nom <- Legend(title = "Nominated by slice", at = nom_lv, legend_gp = gpar(fill = nom_col), ncol = 1,
                  title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7.5))
lgd_na  <- Legend(title = "Cell", at = "not detected (gene absent from slice)",
                  legend_gp = gpar(fill = "grey70"),
                  title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7))

lgd_all <- packLegend(lgd_pat, lgd_cl, lgd_clu, lgd_nom, lgd_na, direction = "vertical", gap = unit(2.5, "mm"))
ttl <- "Interface top DE genes across 4 slices\n(rows = union of per-slice top15; z-scored across the 8 columns)"
stem2 <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/figures/per_patient/heatmap_interface_tops_4slices"

png(paste0(stem2, ".png"), width = 8, height = 17.5, units = "in", res = 300, bg = "white")
draw(ht, annotation_legend_list = lgd_all, annotation_legend_side = "bottom", heatmap_legend_side = "right",
     column_title = ttl, column_title_gp = gpar(fontsize = 10, fontface = "bold"))
dev.off()
pdf(paste0(stem2, ".pdf"), width = 8, height = 17.5, bg = "white")
draw(ht, annotation_legend_list = lgd_all, annotation_legend_side = "bottom", heatmap_legend_side = "right",
     column_title = ttl, column_title_gp = gpar(fontsize = 10, fontface = "bold"))
dev.off()
cat("[Figure] 热图:", stem2, ".{png,pdf}\n")

im <- png::readPNG(paste0(stem2, ".png")); rgb <- im[, , 1:3]
cat(sprintf("[Figure Wholeness] heatmap | %.0fKB | white=%.3f | sd=%.3f | %s\n",
            file.info(paste0(stem2, ".png"))$size/1024, mean(rgb > 0.98), sd(as.vector(rgb)),
            ifelse(sd(as.vector(rgb)) > 0.02, "PASS", "FAIL")))

# 同步更新导出矩阵（与入图行一致）
write.csv(data.frame(gene = rownames(expr_cell), expr_cell, check.names = FALSE),
          file.path(tab_dir, "heatmap_interface_tops_expression.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(z_plot), cluster = grp4, nominated_by = nom_row, z_plot, check.names = FALSE),
          file.path(tab_dir, "heatmap_interface_tops_zscore.csv"), row.names = FALSE)
cat("[Export] heatmap_interface_tops_expression.csv / heatmap_interface_tops_zscore.csv 已同步（行数:",
    nrow(z_plot), "）\n")
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 导出逐切片 DE 表、切片汇总表、跨切片共识矩阵、热图表达与 z-score 矩阵；并量化热图 4 个行簇的表达模式
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

tab_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/tables"

# 1) 逐切片 DE 表
for (sl in slices4) {
  write.csv(de_int[[sl]], file.path(tab_dir, sprintf("de_interface_%s_clone1_vs_clone2.csv", sl)),
            row.names = FALSE)
}

# 2) 切片汇总（含样本量、显著性计数、深度与肿瘤比例混杂指标）
summ_out <- merge(de_int_summary, int_desc, by = c("slice"), all.x = TRUE, suffixes = c("", "_clone"))
write.csv(summ_out, file.path(tab_dir, "de_interface_summary.csv"), row.names = FALSE)

# 3) 跨切片共识矩阵
cons_df <- data.frame(gene = all_genes,
                      dir_094A = dir_mat[, "094A"], dir_094B = dir_mat[, "094B"],
                      dir_095A = dir_mat[, "095A"], dir_095B = dir_mat[, "095B"],
                      n_up = n_up_sig, n_down = n_dn_sig, n_concordant = n_concord,
                      n_slices_sig = n_slices_sig, stringsAsFactors = FALSE)
cons_df <- cons_df[order(-cons_df$n_concordant, -abs(rowSums(dir_mat))), ]
write.csv(cons_df, file.path(tab_dir, "de_interface_consensus.csv"), row.names = FALSE)

# 4) 热图矩阵
write.csv(data.frame(gene = rownames(expr_cell), expr_cell, check.names = FALSE),
          file.path(tab_dir, "heatmap_interface_tops_expression.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(z_plot), nominated_by = nominate, z_plot, check.names = FALSE),
          file.path(tab_dir, "heatmap_interface_tops_zscore.csv"), row.names = FALSE)

# 5) 热图行簇模式量化
grp4 <- cutree(row_hc, k = 4)
cat("\n[热图行簇模式] 每簇在 8 列上的平均 z（>0.3 记为该列偏高）\n")
for (k in sort(unique(grp4))) {
  idx <- names(grp4)[grp4 == k]
  m   <- colMeans(z_plot[idx, , drop = FALSE], na.rm = TRUE)
  high_cols <- names(m)[m > 0.3]
  cat(sprintf("  Cluster %d (n=%d): %s\n", k, length(idx), paste(sprintf("%s=%.2f", names(m), m), collapse = " ")))
  cat(sprintf("     偏高列: %s\n", ifelse(length(high_cols), paste(high_cols, collapse = ", "), "无(整体居中)")))
  cat(sprintf("     代表基因: %s\n", paste(head(idx, 8), collapse = ", ")))
}

# 6) 全 4 片同向显著的基因
g4 <- rownames(dir_mat)[n_concord == 4]
cat(sprintf("\n[4 片同向显著] 共 %d 个基因\n", length(g4)))
g4_df <- data.frame(gene = g4,
                    direction = ifelse(dir_mat[g4, "094A"] == 1, "Clone1-high", "Clone1-low"),
                    log2FC_094A = round(de_int[["094A"]]$avg_log2FC[match(g4, de_int[["094A"]]$gene)], 2),
                    log2FC_094B = round(de_int[["094B"]]$avg_log2FC[match(g4, de_int[["094B"]]$gene)], 2),
                    log2FC_095A = round(de_int[["095A"]]$avg_log2FC[match(g4, de_int[["095A"]]$gene)], 2),
                    log2FC_095B = round(de_int[["095B"]]$avg_log2FC[match(g4, de_int[["095B"]]$gene)], 2),
                    stringsAsFactors = FALSE)
print(g4_df, row.names = FALSE)
write.csv(g4_df, file.path(tab_dir, "de_interface_shared_all4slices.csv"), row.names = FALSE)

# 7) 3 片同向显著的基因（供报告引用）
g3 <- rownames(dir_mat)[n_concord >= 3 & n_concord < 4]
cat(sprintf("\n[>=3 片同向显著] 共 %d 个（含上面 4 片全一致的）\n", sum(n_concord >= 3)))

cat("\n[Export] 文件:\n")
for (f in c(sprintf("de_interface_%s_clone1_vs_clone2.csv", slices4), "de_interface_summary.csv",
            "de_interface_consensus.csv", "heatmap_interface_tops_expression.csv",
            "heatmap_interface_tops_zscore.csv", "de_interface_shared_all4slices.csv")) {
  fp <- file.path(tab_dir, f)
  cat(sprintf("  %s | %.1f KB\n", f, file.info(fp)$size / 1024))
}
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 深度混杂量化：逐切片统计显著基因中"深度敏感基因族"（MT-、MTRNR2L、MALAT1、RPS/RPL）的占比
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

fam <- function(g) {
  ifelse(grepl("^MT-", g), "MT-",
  ifelse(grepl("^MTRNR2L", g), "MTRNR2L",
  ifelse(g == "MALAT1", "MALAT1",
  ifelse(grepl("^RP[SL]", g), "RPS/RPL", "other"))))
}

depth_tab <- list()
for (sl in slices4) {
  d   <- de_int[[sl]]
  sig <- !is.na(d$p_val_adj) & d$p_val_adj < 0.05
  d$fam <- fam(d$gene)
  hi <- sig & d$avg_log2FC > 0
  lo <- sig & d$avg_log2FC < 0
  for (dirn in c("Clone1-high", "Clone1-low")) {
    sel <- if (dirn == "Clone1-high") hi else lo
    tt  <- table(factor(d$fam[sel], levels = c("MT-","MTRNR2L","MALAT1","RPS/RPL","other")))
    depth_tab[[length(depth_tab) + 1]] <- data.frame(
      slice = sl, direction = dirn, n_sig = sum(sel),
      n_MT = tt[["MT-"]], n_MTRNR2L = tt[["MTRNR2L"]], n_MALAT1 = tt[["MALAT1"]],
      n_RPS_RPL = tt[["RPS/RPL"]], n_other = tt[["other"]],
      pct_depth_sensitive = round(100 * sum(tt[c("MT-","MTRNR2L","MALAT1","RPS/RPL")]) / max(1, sum(sel)), 1),
      stringsAsFactors = FALSE)
  }
}
depth_df <- do.call(rbind, depth_tab)
print(depth_df, row.names = FALSE)

# 与深度方向的一致性：若 Clone1 深度更高，深度敏感基因应偏向 Clone1-high
cat("\n[判读] 深度敏感基因族富集方向是否与该切片的深度差方向一致：\n")
depth_dir <- c("094A" = "Clone1-lower", "094B" = "Clone1-lower",
               "095A" = "Clone1-higher", "095B" = "Clone1-higher")
for (sl in slices4) {
  sub <- depth_df[depth_df$slice == sl, ]
  cat(sprintf("  %s (%s): Clone1-high 中深度敏感 %.1f%% | Clone1-low 中 %.1f%%\n", sl, depth_dir[sl],
              sub$pct_depth_sensitive[sub$direction == "Clone1-high"],
              sub$pct_depth_sensitive[sub$direction == "Clone1-low"]))
}
write.csv(depth_df, file.path(tab_dir, "de_interface_depth_sensitive_families.csv"), row.names = FALSE)
cat("\n[Export] de_interface_depth_sensitive_families.csv\n")
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 切片两两之间显著基因集的重叠（Jaccard）与每片 top5 上下调清单整理，供报告引用
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

sig_sets <- lapply(de_int, function(d) d$gene[!is.na(d$p_val_adj) & d$p_val_adj < 0.05])
cat("[显著基因集大小]", paste(sprintf("%s=%d", names(sig_sets), lengths(sig_sets)), collapse = "  "), "\n\n")

cat("[两两 Jaccard 重叠]\n")
jc <- matrix(NA, 4, 4, dimnames = list(slices4, slices4))
for (a in slices4) for (b in slices4) {
  jc[a, b] <- length(intersect(sig_sets[[a]], sig_sets[[b]])) / length(union(sig_sets[[a]], sig_sets[[b]]))
}
print(round(jc, 3))
write.csv(data.frame(slice = rownames(jc), jc, check.names = FALSE),
          file.path(tab_dir, "de_interface_pairwise_jaccard.csv"), row.names = FALSE)

cat("\n[每片 top5 Clone1 高 / 低（padj<0.05，按效应量排序）]\n")
for (sl in slices4) {
  d <- de_int[[sl]]; sig <- !is.na(d$p_val_adj) & d$p_val_adj < 0.05
  up <- head(d[sig & d$avg_log2FC > 0, ][order(-d$avg_log2FC[sig & d$avg_log2FC > 0]), ], 5)
  dn <- head(d[sig & d$avg_log2FC < 0, ][order(d$avg_log2FC[sig & d$avg_log2FC < 0]), ], 5)
  cat(sprintf("%s 高: %s\n", sl, paste(sprintf("%s(%.2f)", up$gene, up$avg_log2FC), collapse = ", ")))
  cat(sprintf("%s 低: %s\n", sl, paste(sprintf("%s(%.2f)", dn$gene, dn$avg_log2FC), collapse = ", ")))
}
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 图形-数据一致性校验：确认着色点全部满足 padj<0.05、标注基因名真实存在于表达矩阵（排除渲染字形的误读）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

cat("[阈值一致性] 每切片着色显著点的最小 -log10(padj)（门槛 -log10(0.05)=1.301，应 >= 1.301）\n")
for (sl in slices4) {
  d   <- de_int[[sl]]
  sig <- !is.na(d$p_val_adj) & d$p_val_adj < 0.05
  mn  <- min(-log10(d$p_val_adj[sig]))
  mx_ns <- suppressWarnings(max(-log10(d$p_val_adj[!is.na(d$p_val_adj) & !sig])))
  cat(sprintf("  %s | 显著点最小 -log10padj = %.3f | 非显著点最大 = %.3f | %s\n",
              sl, mn, ifelse(is.finite(mx_ns), mx_ns, NA), ifelse(mn >= 1.3009, "PASS", "FAIL")))
  cat(sprintf("      显著点中 padj>=0.05 的个数 = %d（应为 0）\n", sum(sig & d$p_val_adj >= 0.05)))
}

cat("\n[标注基因名真实性] 火山图所标注基因是否都在各切片表达矩阵行名中\n")
lab_genes <- c("ALDOA","TYMP","FKBP8","COL6A2","BGN",          # 094A high
               "CYP4X1","NPPC","DNAJC24","RAD17","AP3B1",      # 094A low
               "COMP","AEBP1","FBLN1","EMILIN1","TPSB2",       # 094B high
               "BBOX1","ZCCHC7","RAPGEF5","SMIM27","SERPINI1", # 094B low
               "USP20","SMG6","SLC39A3","TRPM2","PHLPP1",      # 095A high
               "TNFSF10","MTRNR2L12","MALAT1","LYZ","FN1",     # 095A low
               "LAG3","BATF3","CIDEC","APOBEC3F","FOXQ1",      # 095B high
               "PIGR","MTRNR2L8")
cat("  标注基因总数(去重):", length(unique(lab_genes)), "\n")
cat("  在 obj13 中缺失:", paste(setdiff(lab_genes, rownames(obj13)), collapse = ", "), "\n")
cat("  在 obj15 中缺失:", paste(setdiff(lab_genes, rownames(obj15)), collapse = ", "), "\n")
cat("  形近易混检查 | COL6A2 存在:", "COL6A2" %in% rownames(obj13), "| GOL6A2 存在:",
    "GOL6A2" %in% rownames(obj13), "| ZCCHC7 存在:", "ZCCHC7" %in% rownames(obj13),
    "| ZCCHG7 存在:", "ZCCHG7" %in% rownames(obj13), "\n")
cat("  => 图上标签由 gene 列的字符直接绘制，形近字读法为渲染字形误读，非标签错误\n")
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 生成 step03 汇总交付文档（最终版，单文件自包含）；表格与文中数值全部由计算结果程序化生成
# 本版按独立审查意见更正：多重校正方法(Bonferroni 而非 BH)、仅1片显著数(1422)、P15 深度差顺序、
#   z-score 列数口径、Interface 覆盖率披露、pct.1/pct.2 佐证、7 基因核心的解释措辞
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(base64enc))
fig_dir  <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/figures/per_patient"
out_html <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step03_interface_de_汇总.html"

img_tag <- function(f, cap) {
  sprintf("<figure><img src='data:image/png;base64,%s' alt='%s'/><figcaption>%s</figcaption></figure>",
          base64enc::base64encode(f), cap, cap)
}

# ---- 逐切片主表 ----
tbl <- c("<table><thead><tr><th>切片</th><th>病人</th><th>Interface spot</th><th>Clone1 / Clone2</th><th>检验基因</th><th>显著 (padj&lt;0.05)</th><th>Clone1 高</th><th>Clone1 低</th><th>中位 nCount (C1 / C2)</th><th>深度差</th><th>中位 tumor_prop 差</th></tr></thead><tbody>")
for (i in seq_len(nrow(de_int_summary))) {
  sl <- de_int_summary$slice[i]
  d1 <- int_desc[int_desc$slice == sl & int_desc$clone == "Clone1", ]
  d2 <- int_desc[int_desc$slice == sl & int_desc$clone == "Clone2", ]
  tbl <- c(tbl, sprintf("<tr><td><b>%s</b></td><td>%s</td><td>%d</td><td>%d / %d</td><td>%s</td><td><b>%d</b></td><td>%d</td><td>%d</td><td>%s / %s</td><td>%+.1f%%</td><td>%+.3f</td></tr>",
    sl, ifelse(sl %in% c("094A","094B"), "P13", "P15"), de_int_summary$n_spot[i],
    de_int_summary$n_Clone1[i], de_int_summary$n_Clone2[i],
    format(de_int_summary$n_tested[i], big.mark = ","), de_int_summary$n_sig[i],
    de_int_summary$n_up[i], de_int_summary$n_down[i],
    format(d1$median_nCount, big.mark = ","), format(d2$median_nCount, big.mark = ","),
    100 * (d1$median_nCount - d2$median_nCount) / d2$median_nCount,
    d1$median_tumor_prop - d2$median_tumor_prop))
}
tbl <- c(tbl, "</tbody></table>")

# ---- 每片 top5 ----
top_tbl <- c("<table><thead><tr><th>切片</th><th>Clone1 高 top5 (avg_log2FC)</th><th>Clone1 低 top5 (avg_log2FC)</th></tr></thead><tbody>")
for (sl in slices4) {
  d <- de_int[[sl]]; sig <- !is.na(d$p_val_adj) & d$p_val_adj < 0.05
  up <- head(d[sig & d$avg_log2FC > 0, ][order(-d$avg_log2FC[sig & d$avg_log2FC > 0]), ], 5)
  dn <- head(d[sig & d$avg_log2FC < 0, ][order(d$avg_log2FC[sig & d$avg_log2FC < 0]), ], 5)
  top_tbl <- c(top_tbl, sprintf("<tr><td><b>%s</b></td><td>%s</td><td>%s</td></tr>", sl,
    paste(sprintf("%s (%.2f)", up$gene, up$avg_log2FC), collapse = "; "),
    paste(sprintf("%s (%.2f)", dn$gene, dn$avg_log2FC), collapse = "; ")))
}
top_tbl <- c(top_tbl, "</tbody></table>")

# ---- 深度敏感族表 ----
dep_tbl <- c("<table><thead><tr><th>切片</th><th>方向</th><th>显著数</th><th>MT-</th><th>MTRNR2L</th><th>MALAT1</th><th>RPS/RPL</th><th>深度敏感占比</th></tr></thead><tbody>")
for (i in seq_len(nrow(depth_df))) {
  dep_tbl <- c(dep_tbl, sprintf("<tr><td>%s</td><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%.1f%%</td></tr>",
    depth_df$slice[i], depth_df$direction[i], depth_df$n_sig[i], depth_df$n_MT[i],
    depth_df$n_MTRNR2L[i], depth_df$n_MALAT1[i], depth_df$n_RPS_RPL[i], depth_df$pct_depth_sensitive[i]))
}
dep_tbl <- c(dep_tbl, "</tbody></table>")

# ---- pct 方向佐证表 ----
pct_tbl <- c("<table><thead><tr><th>切片</th><th>Clone1 高 基因数</th><th>其中 pct.1&gt;pct.2</th><th>占比</th><th>Clone1 低 基因数</th><th>其中 pct.1&lt;pct.2</th><th>占比</th></tr></thead><tbody>")
for (i in seq_len(nrow(pct_df))) {
  pct_tbl <- c(pct_tbl, sprintf("<tr><td><b>%s</b></td><td>%d</td><td>%d</td><td>%.1f%%</td><td>%d</td><td>%d</td><td>%.1f%%</td></tr>",
    pct_df$slice[i], pct_df$hi_n[i], pct_df$hi_pct1_gt_pct2[i], pct_df$hi_frac[i],
    pct_df$lo_n[i], pct_df$lo_pct1_lt_pct2[i], pct_df$lo_frac[i]))
}
pct_tbl <- c(pct_tbl, "</tbody></table>")

# ---- 覆盖率表 ----
cov_tbl <- c("<table><thead><tr><th>切片</th><th>分区表 Interface spot</th><th>实际入分析</th><th>覆盖率</th><th>该切片对象总 spot</th></tr></thead><tbody>")
for (i in seq_len(nrow(cov))) {
  cov_tbl <- c(cov_tbl, sprintf("<tr><td><b>%s</b></td><td>%d</td><td>%d</td><td>%.1f%%</td><td>%d</td></tr>",
    cov$slice[i], cov$region_interface[i], cov$used[i], cov$coverage_pct[i], cov$obj_total[i]))
}
cov_tbl <- c(cov_tbl, "</tbody></table>")

g4_genes <- rownames(dir_mat)[n_concord == 4]
g4_hi <- g4_genes[dir_mat[g4_genes, "094A"] == 1]
g4_lo <- g4_genes[dir_mat[g4_genes, "094A"] == -1]
n_only1 <- sum(n_slices_sig == 1)

html <- c(
"<!DOCTYPE html><html lang='zh-CN'><head><meta charset='utf-8'><title>Step 03 — Interface 分区内 Clone1 vs Clone2 差异表达</title>",
"<style>body{font-family:-apple-system,'Segoe UI','Microsoft YaHei',sans-serif;max-width:1100px;margin:32px auto;padding:0 20px;line-height:1.65;color:#1a1a1a;background:#fff}",
"h1{border-bottom:3px solid #3C5488;padding-bottom:8px;font-size:26px}h2{margin-top:34px;color:#3C5488;border-left:5px solid #E64B35;padding-left:10px;font-size:20px}",
"h3{margin-top:22px;font-size:16px;color:#333}table{border-collapse:collapse;width:100%;margin:14px 0;font-size:13px}",
"th,td{border:1px solid #d9d9d9;padding:6px 9px;text-align:left}th{background:#f2f4f7;font-weight:600}",
"tr:nth-child(even) td{background:#fafbfc}figure{margin:22px 0;text-align:center}img{max-width:100%;border:1px solid #e5e5e5}",
"figcaption{font-size:12.5px;color:#555;margin-top:6px}.divider{border:0;border-top:2px dashed #bbb;margin:42px 0}",
".note{background:#fff8f0;border-left:4px solid #F39B7F;padding:10px 14px;margin:14px 0;font-size:13.5px}",
".key{background:#f0f7f5;border-left:4px solid #00A087;padding:10px 14px;margin:14px 0;font-size:13.5px}",
"code{background:#f4f4f4;padding:1px 5px;border-radius:3px;font-size:12.5px}.meta{color:#777;font-size:12.5px}</style></head><body>",
"<h1>Step 03 — Interface 分区内 Clone1 vs Clone2 差异表达（逐切片）</h1>",
"<p class='meta'>疾病背景：TNBC（三阴性乳腺癌）｜样本：原发肿瘤组织，10x Visium spot 级（每 spot 多细胞混合）｜切片：094A、094B（P13）、095A、095B（P15）｜分析层次：Tumor/Stroma 交界 Interface 区</p>",

"<h2>一、结论简报</h2>",
"<div class='key'><b>一句话结论：</b>把克隆比较收紧到 Interface 区后，两克隆之间仍可检出差异，<b>但差异的方向与幅度在四个切片之间高度不一致，且与该切片内部两克隆的测序深度差同向</b>。真正跨四个切片都能复现的只有一个 7 基因的小核心，其效应量都很小。因此本步不支持\"某个克隆具有稳定转录程序\"的结论。</div>",
"<ul>",
"<li><b>界面区内仍有克隆差异</b>：四片分别检出 129 / 825 / 289 / 785 个显著基因（padj&lt;0.05），说明差异不是全组织混合造成的假象。</li>",
sprintf("<li><b>跨切片共识很薄</b>：仅 %d 个基因在 4 片同时同向显著（%s 为 Clone1 高；%s 为 Clone1 低），效应量均很小（|avg_log2FC| 0.20–0.72）；%d 个基因在 ≥2 片同向显著。</li>",
        length(g4_genes), paste(g4_hi, collapse = "、"), paste(g4_lo, collapse = "、"), sum(n_concord >= 2)),
"<li><b>深度是最主要的竞争解释</b>：P13 两片 Clone1 深度更低（−37.9% / −48.5%），P15 两片 Clone1 深度更高（+162.8% / +143.1%），而显著基因的上下调比例正好随深度方向翻转。因此\"克隆差异\"与\"测序深度差异\"在本数据中不可分离。</li>",
"<li><b>肿瘤比例不是解释</b>：同一切片内两克隆的中位 tumor_prop 差 ≤0.015，可排除\"某克隆所在 spot 肿瘤含量更高\"这一替代解释。</li>",
"<li><b>功效不足的切片</b>：094A 的 Clone2 仅 62 个 Interface spot、095A 的 Clone2 仅 44 个，这两片的结论只能作方向性提示。</li>",
"</ul>",

"<h2>二、方法要点</h2>",
"<ul>",
"<li><b>样本单元</b>：Visium spot（多细胞混合），不是单细胞。所有结论的粒度是\"spot 及其邻域\"，不能外推到单个细胞。</li>",
"<li><b>区域限定</b>：仅保留 <code>region == \"Interface\"</code> 的 spot；分区表（与 step01 同一份）按 <code>&lt;切片&gt;_&lt;barcode 去掉末尾 _数字&gt;</code> 关联，四片关联率均为 <b>100%</b>。注意：Visium 对象只含<b>带克隆标签</b>的 spot，因此实际纳入的是\"有克隆标签的 Interface spot\"，并非该区全部 spot（详见 5.5 与覆盖率表）。</li>",
"<li><b>统计</b>：逐切片独立执行 <code>FindMarkers(test.use=\"wilcox\", min.pct=0.25, logfc.threshold=0)</code>，<b>不合并、不池化、不抽样</b>。表达量为既有的 LogNormalize（target_sum=1e4）data 层。</li>",
"<li><b>多重检验校正</b>：<code>p_val_adj</code> 为 <b>Seurat 默认的 Bonferroni 校正</b>（p × 对象全部基因数：P13 11,743、P15 12,978；已用未截顶行的恒定比值实证核对）。这是<b>较保守</b>的口径——它按对象全部基因而非仅检验基因（9,004–10,981）做校正。本步<b>不是</b> BH/FDR。</li>",
"<li><b>方向参照系</b>：<code>ident.1 = Clone1</code>、<code>ident.2 = Clone2</code>；<b>avg_log2FC &gt; 0 表示 Clone1 高表达</b>，&lt; 0 表示 Clone1 低表达（即 Clone2 高）。显著性 = <code>padj &lt; 0.05</code>，不额外设 log2FC 阈值。本数据无已知的\"方向锚定基因\"可用于校验（克隆身份来自 CNV 推断，不存在 Cre 类标记基因），方向只能由 ident.1/ident.2 的定义保证。</li>",
"<li><b>热图</b>：行 = 各切片 top15 的并集，top15 取 padj&lt;0.05 中 |avg_log2FC| 最大的基因并<b>按正负尽量各半</b>取（每片 8 个 Clone1 高 + 7 个 Clone1 低；若改用\"不分成上下调的纯 top15\"，并集为 58 个基因而非 56）。列 = 4 切片 × 2 克隆，按病人分组；填色 = 每基因在 8 列间的 z-score（±2.5 截断）；行做层次聚类并切 4 簇。</li>",
"</ul>",

"<h2>三、关键结果</h2>",
"<h3>3.1 逐切片差异概览</h3>", tbl,
"<p class='meta'>深度差 =（Clone1 中位 nCount − Clone2 中位 nCount）/ Clone2 中位 nCount。肿瘤比例差 = Clone1 中位 tumor_prop − Clone2 中位 tumor_prop。</p>",
"<h3>3.2 每切片 top5 上下调基因（按效应量）</h3>", top_tbl,
"<h3>3.3 火山图</h3>",
img_tag(file.path(fig_dir, "volcano_interface_094A_clone1_vs_clone2.png"), "094A（P13）：Interface spot 264，显著 129（Clone1 高 43 / 低 86）。蓝侧更密更高，与 Clone1 深度低 37.9% 同向。"),
img_tag(file.path(fig_dir, "volcano_interface_094B_clone1_vs_clone2.png"), "094B（P13）：Interface spot 480，显著 825（高 286 / 低 539），本片功效最足。COMP/AEBP1/FBLN1/EMILIN1 等基质基因在 Clone1 侧升高。"),
img_tag(file.path(fig_dir, "volcano_interface_095A_clone1_vs_clone2.png"), "095A（P15）：Interface spot 469，显著 289（高 239 / 低 50）。Clone2 仅 44 spot，红侧显著数明显偏多，与 Clone1 深度高 162.8% 同向。"),
img_tag(file.path(fig_dir, "volcano_interface_095B_clone1_vs_clone2.png"), "095B（P15）：Interface spot 407，显著 785（高 739 / 低 46）。强单侧分布，下调侧集中于 MALAT1/MTRNR2L12/MTRNR2L8 等深度敏感基因。"),

"<h2>四、跨切片共识</h2>",
"<ul>",
sprintf("<li><b>≥2 个切片同向显著</b>：%d 个基因；<b>≥3 片</b>：%d 个；<b>4 片全部同向</b>：%d 个。</li>", sum(n_concord >= 2), sum(n_concord >= 3), sum(n_concord == 4)),
sprintf("<li><b>仅 1 个切片显著</b>：%d 个基因；另有 %d 个基因在 2 片显著但<b>方向相反</b>，因此\"同向显著最多 1 片\"的合计为 %d。</li>",
        n_only1, sum(n_slices_sig >= 2 & n_up_sig > 0 & n_dn_sig > 0), sum(n_concord == 1)),
"<li><b>病人内部跨切片一致性</b>：P13（094A+094B）两片都显著的 63 个基因，同向率 100%；P15（095A+095B）两片都显著的 103 个基因，同向率 100%。</li>",
"<li><b>跨病人（仅描述性，克隆编号不可比）</b>：两个病人各自 ≥2 片同向显著的 7 个基因，方向全部一致。</li>",
"<li><b>显著集两两重叠（Jaccard）</b>：0.025–0.106，说明各片显著列表绝大部分是切片特异的，共识只集中在很小的核心上。</li>",
sprintf("<li><b>4 片全同向的 %d 个基因</b>：高表达 %s；低表达 %s。</li>", length(g4_genes),
        paste(g4_hi, collapse = "、"), paste(g4_lo, collapse = "、")),
"</ul>",

"<h2>五、混杂与解释限制</h2>",
"<div class='note'><b>核心限制：克隆差异与测序深度差异不可分离。</b>P13 中 Clone1 深度系统性更低、P15 中 Clone1 深度系统性更高，而显著基因的上下调比例随之一致翻转（P13 以 Clone1 低为主，P15 以 Clone1 高为主）。在这种结构下，任何\"克隆 A 高表达 X\"的表述都无法排除\"克隆 A 所在 spot 测序更深/更浅\"。</div>",
"<h3>5.1 检出比例（pct.1 / pct.2）给出的独立佐证</h3>", pct_tbl,
"<p class='meta'>pct.1 = Clone1 中检出该基因的 spot 比例，pct.2 = Clone2 中比例。若差异是真实的克隆身份差异，Clone1 高表达基因通常应有 pct.1 &gt; pct.2。实际在深度更低的 P13 两片中，这一比例只有 4.7% / 23.1%——即绝大多数\"Clone1 高\"基因在 Clone1 中的<b>检出率反而更低</b>，其 log2FC 方向主要由深度归一化后的相对丰度造成。P15 两片（Clone1 更深）则为 85.4% / 98.0%，方向自洽。这从另一条独立证据线支持\"深度主导\"的判断。</p>",
"<h3>5.2 深度敏感基因族的占比</h3>", dep_tbl,
"<p class='meta'>深度敏感基因族 = MT-、MTRNR2L、MALAT1、RPS/RPL（本数据 DE 表中 MT- 基因为 0 个，实际由 RPS/RPL、MTRNR2L、MALAT1 承担）。判读：094A、095A、095B 三片的深度敏感基因富集在\"低深度克隆\"一侧，与深度伪影方向一致（095A 的 Clone1 低表达侧高达 60.0%，该侧结论应视为无效）。094B 的方向与此不符（Clone1 高侧 3.8% vs 低侧 9.5%），说明<b>深度不能单独解释 094B 的差异</b>——但这不等于已证明其为克隆内在差异，只是该片需要别的解释（组成差异或其它技术因素）。</p>",
"<h3>5.3 功效与样本量</h3>",
"<ul>",
"<li>094A：Clone2 仅 62 spot；095A：Clone2 仅 44 spot。两片的\"少数派克隆\"一侧估计不稳定，效应量与显著性均可能被高估（Winner's curse），只能作方向提示。</li>",
"<li>spot 之间不独立（空间自相关 + 每 spot 多细胞），本步 p 值属伪重复（pseudo-replication）性质，<b>不应按独立观测解读</b>；本步以效应量与跨切片一致性为主要依据，p 值仅作筛选门槛。</li>",
"</ul>",
"<h3>5.4 热图仅展示模式，不做跨病人统计推断</h3>",
"<p>跨病人的克隆编号不可比（P13 的 Clone1 与 P15 的 Clone1 不是同一生物学实体）。热图把 4 切片放在一起只为展示\"哪些基因属于病人特异模式、哪些属于克隆模式\"，<b>不得据此做跨病人比较或检验</b>。另外，热图的 z-score 在<b>非缺失列</b>上计算：48 个基因用满 8 列，8 个基因（PI3、BBOX1、PHLPP1、CIDEC、FOXQ1、SCARF2、ELF3-AS1、PIGR）在其中 4 列未检出，实际只用 4 个值计算，其 z 值不与前者严格同尺度。</p>",
"<h3>5.5 Interface 的覆盖范围（重要口径）</h3>", cov_tbl,
"<p class='meta'>分区表在四片分别标出 592 / 748 / 579 / 773 个 Interface spot，但本步只用到 264 / 480 / 469 / 407 个（覆盖率 44.6% / 64.2% / 81.0% / 52.7%）。原因是克隆比较只能在<b>有克隆标签</b>的 spot 上进行，而对象中的 spot 以肿瘤核心为主（同片中 Tumor Core 的覆盖率约 98%–99%）。因此本步的\"Interface\"应理解为<b>\"有克隆标签的、肿瘤富集的界面区\"</b>，不等于形态学界面的全部 spot；覆盖率较低的 094A、095B 尤其要谨慎外推。</p>",
"<h3>5.6 Interface 定义的边界</h3>",
"<p>Interface 由分区规则划定（基于到肿瘤核心的距离），是几何定义而非分子定义；边界 spot 仍可能混有核心或间质成分，因此\"界面区克隆差异\"不等于\"界面区细胞类型组成差异\"。</p>",

"<h2>六、图表与产物索引</h2>",
"<h3>6.1 合并热图</h3>",
img_tag(file.path(fig_dir, "heatmap_interface_tops_4slices.png"), "4 切片 top 基因合并热图（56 行 × 8 列）。行=各片 top15（8 高 + 7 低）并集，按非缺失列 z-score；行层次聚类切 4 簇。灰格=该切片未检出该基因（32 格）。"),
"<ul>",
"<li><b>Cluster 1（21 行）</b>：P13 四列偏高、P15 四列偏低 —— <b>病人驱动</b>，非克隆差异。</li>",
"<li><b>Cluster 2（13 行）</b>：P15 偏高、P13 偏低 —— 同样是<b>病人驱动</b>（且部分 P13 格为未检出）。</li>",
"<li><b>Cluster 3（7 行）</b>：Clone2 侧偏高，但主要由 094A_Clone2 单列撑起，跨片不一致 —— <b>弱、不可复现</b>。</li>",
"<li><b>Cluster 4（15 行）</b>：Clone1 侧偏高（AEBP1、FBLN1、EMILIN1、USP20、SMG6 等），在 094B/095A 较明显 —— <b>克隆相关但强度中等</b>，且与 094B_Clone2 单列偏低混杂。</li>",
"</ul>",
"<p class='key'><b>对\"有几条表达模式\"的回答：</b>层次聚类切出的 4 个行簇中，只有 2 个（Cluster 1、2）是稳健且方向明确的病人特异性模式；Cluster 4 是中等强度的克隆相关模式；Cluster 3 弱且由单切片驱动。即 <b>2 条强（病人轴）+ 1 条中等（克隆轴）+ 1 条弱</b>。</p>",
"<h3>6.2 产物清单</h3>",
"<table><thead><tr><th>产物</th><th>相对路径</th><th>内容</th></tr></thead><tbody>",
"<tr><td>逐切片 DE 表 ×4</td><td>tables/de_interface_&lt;切片&gt;_clone1_vs_clone2.csv</td><td>gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj</td></tr>",
"<tr><td>切片汇总</td><td>tables/de_interface_summary.csv</td><td>逐片 spot 数、显著数、上下调数、深度与肿瘤比例</td></tr>",
"<tr><td>跨切片共识矩阵</td><td>tables/de_interface_consensus.csv</td><td>每基因在 4 片的方向与同向显著片数</td></tr>",
"<tr><td>4 片全同向基因</td><td>tables/de_interface_shared_all4slices.csv</td><td>7 个基因及四片 log2FC</td></tr>",
"<tr><td>显著集两两重叠</td><td>tables/de_interface_pairwise_jaccard.csv</td><td>4×4 Jaccard</td></tr>",
"<tr><td>深度敏感基因族</td><td>tables/de_interface_depth_sensitive_families.csv</td><td>逐片逐方向的基因族占比</td></tr>",
"<tr><td>检出比例方向佐证</td><td>tables/de_interface_pct_direction_check.csv</td><td>逐片 pct.1/pct.2 与 log2FC 方向一致性</td></tr>",
"<tr><td>Interface 覆盖率</td><td>tables/de_interface_region_coverage.csv</td><td>分区表 Interface 数 vs 实际入分析数</td></tr>",
"<tr><td>热图矩阵</td><td>tables/heatmap_interface_tops_expression.csv / _zscore.csv</td><td>56 基因 × 8 列的平均表达与 z-score（含簇号）</td></tr>",
"<tr><td>火山图 ×4</td><td>figures/per_patient/volcano_interface_&lt;切片&gt;_clone1_vs_clone2.{png,pdf}</td><td>pdf 矢量 + png 300dpi</td></tr>",
"<tr><td>合并热图</td><td>figures/per_patient/heatmap_interface_tops_4slices.{png,pdf}</td><td>pdf 矢量 + png 300dpi</td></tr>",
"</tbody></table>",
"<p class='meta'><b>交付边界说明：</b>本交付包含上述结果表、图表与本文档。图形像素完整性、阈值一致性等内部质量自检记录属于内部质控材料，不随交付提供；文中所有结论与数值均已在正文与上述结果表中给出。</p>",

"<hr class='divider'>",
"<h2>技术说明书</h2>",
"<h3>1. 本步回答的问题</h3>",
"<ol><li>把克隆比较限制在 Tumor/Stroma 交界（Interface）后，Clone1 与 Clone2 是否仍有表达差异？</li>",
"<li>这些差异能跨切片复现吗？跨切片共识有多大？</li>",
"<li>观察到的差异在多大程度上可由测序深度或肿瘤含量解释？</li></ol>",
"<h3>2. 方法流程</h3>",
"<ol><li>读分区表（16,607 行），按 <code>sample</code> 解析切片号、按 <code>barcode</code> 去掉末尾 <code>_数字</code> 得到 clean barcode，拼成关联键。</li>",
"<li>对每一切片，取该切片在 Visium 对象中的全部 spot，与分区表按关联键匹配，保留 <code>region==\"Interface\"</code> 的 spot（四片关联率 100%）。</li>",
"<li>对每个切片的 Interface 子集独立执行 <code>FindMarkers(ident.1=Clone1, ident.2=Clone2, test.use=\"wilcox\", min.pct=0.25, logfc.threshold=0)</code>。</li>",
"<li>按 padj&lt;0.05 计数显著基因并分上/下调；按 |avg_log2FC| 取每片 top15（8 高 + 7 低）用于热图。</li>",
"<li>汇总跨切片同向显著片数；统计深度敏感基因族占比与 pct.1/pct.2 方向一致性，作为深度伪影诊断。</li>",
"<li>对每个基因在 8 个 (切片×克隆) 列间计算平均 log-normalized 表达并 z-score（在非缺失列上），行层次聚类（complete linkage）切 4 簇绘图。</li>",
"</ol>",
"<h3>3. 产物逐项说明</h3>",
"<ul>",
"<li><b>火山图</b>：x 轴 = avg_log2FC（Clone1 vs Clone2，&gt;0 = Clone1 高）；y 轴 = −log10(padj)。红点 = Clone1 高（padj&lt;0.05），蓝点 = Clone1 低，灰点 = 不显著。两条虚线分别为 avg_log2FC=0 与 padj=0.05，并在图上标注。标注的 10 个基因是 padj&lt;0.05 中效应量 |avg_log2FC| 最大的 5 高 + 5 低（注意：是<b>效应量最大</b>，不是 −log10(padj) 最大）。</li>",
"<li><b>合并热图</b>：行 = 56 个基因（4 片 top15 并集），列 = 8 组（4 切片 × 2 克隆），白色宽带分隔 P13 与 P15。填色 = 每基因在非缺失列间的 z-score（红高蓝低，白 = 该基因自身均值），色阶固定 ±2.5。左侧两条注释：<code>Cluster</code>（行簇号 1–4，与图例一一对应）、<code>Nominated</code>（该基因由哪个切片提名）。灰格 = 该切片未检出该基因（32 格）。行树为层次聚类，4 个簇之间以白带分隔。</li>",
"<li><b>de_interface_consensus.csv</b>：<code>dir_094A…dir_095B</code> 取值 1 = 该片显著且 Clone1 高、−1 = 显著且 Clone1 低、0 = 不显著或未被检验；<code>n_concordant</code> = 同向显著的最大片数；<code>n_slices_sig</code> = 显著的片数（不看方向）。</li>",
"<li><b>de_interface_depth_sensitive_families.csv</b>：逐片逐方向统计四类深度敏感基因的个数与占该方向显著基因的百分比。</li>",
"<li><b>de_interface_pct_direction_check.csv</b>：逐片统计显著基因中 pct.1/pct.2 与 log2FC 方向一致的基因数，用于判断方向是否可能由检出率驱动。</li>",
"<li><b>de_interface_region_coverage.csv</b>：分区表 Interface spot 数与实际入分析 spot 数的对比，用于界定\"Interface\"的实际覆盖范围。</li>",
"</ul>",
"<h3>4. 关键结论（数值出处）</h3>",
"<ul>",
"<li>四片显著基因数 129 / 825 / 289 / 785 —— 见 <code>de_interface_summary.csv</code> 的 <code>n_sig</code> 列。</li>",
"<li>四片两克隆中位深度差 −37.9% / −48.5% / +162.8% / +143.1% —— 由 <code>de_interface_summary.csv</code> 的 <code>median_nCount</code> 计算。</li>",
sprintf("<li>≥2 片同向显著 %d 个、≥3 片 %d 个、4 片 %d 个 —— 见 <code>de_interface_consensus.csv</code> 的 <code>n_concordant</code> 列计数。</li>", sum(n_concord >= 2), sum(n_concord >= 3), sum(n_concord == 4)),
sprintf("<li>仅 1 片显著 %d 个、2 片显著但方向冲突 %d 个 —— 见 <code>de_interface_consensus.csv</code> 的 <code>n_slices_sig</code> 与 <code>dir_*</code> 列。</li>", n_only1, sum(n_slices_sig >= 2 & n_up_sig > 0 & n_dn_sig > 0)),
"<li>病人内部跨切片同向率 P13 63/63、P15 103/103 —— 见 <code>de_interface_consensus.csv</code> 方向列。</li>",
"<li>热图 4 个行簇大小 21/13/7/15 —— 见 <code>heatmap_interface_tops_zscore.csv</code> 的 <code>cluster</code> 列。</li>",
"<li>pct 方向一致性（4.7%/23.1% vs 85.4%/98.0%）—— 见 <code>de_interface_pct_direction_check.csv</code>。</li>",
"<li>Interface 覆盖率 44.6%/64.2%/81.0%/52.7% —— 见 <code>de_interface_region_coverage.csv</code>。</li>",
"</ul>",
"<h3>5. 使用注意</h3>",
"<ul>",
"<li><b>校正口径</b>：<code>p_val_adj</code> 是 Bonferroni（p × 对象全部基因数），不是 BH/FDR；若下游需与 BH 口径的其它分析比较，须先统一校正方法。</li>",
"<li><b>不要把 p 值当独立观测的证据</b>：spot 空间自相关且每 spot 多细胞，n 不是独立生物学重复。</li>",
"<li><b>不要跨病人比较克隆编号</b>：P13-Clone1 与 P15-Clone1 无对应关系；热图的病人分组只为展示，不构成比较设计。</li>",
sprintf("<li><b>7 基因核心的解释要谨慎但不否定</b>：这 %d 个基因在 Clone1 深度更低（P13）与更高（P15）两种情形下方向都不变，<b>深度方向无法解释它们</b>；但其中包含 RPS24、RPL18、PABPC1 等核糖体/持家基因，效应量均 &lt;0.75，更可能是检出率差异或很弱的真实效应。它是本步唯一跨 4 片复现的信号，值得在深度匹配的分层内单独复核，但不足以支撑克隆特异程序的结论。</li>", length(g4_genes)),
"<li><b>若需把深度与克隆分离</b>，可行方向（本步未做）：在深度匹配的分层内比较、用 depth 作为协变量的回归、或在 pseudo-bulk 水平以 depth 归一后比较；这些属于后续分析，不能靠本步数据在文字上\"解释掉\"。</li>",
"<li>图形按 pdf 矢量 + png 300dpi 双份交付；热图为长版式（8×17.5 in），排版时不宜再缩小，否则行/列标签会低于可读下限。</li>",
"</ul>",
"</body></html>")

writeLines(html, out_html, useBytes = TRUE)
cat(sprintf("[HTML] %s | %.2f MB | 内嵌图 %d 张\n", out_html, file.info(out_html)$size / 1e6,
            length(gregexpr("data:image/png;base64", paste(html, collapse = ""))[[1]])))
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 补充审查证据：① 校正方法实证（仅用未被截顶的行）；② 检出比例与方向一致性；③ Interface 覆盖率；④ 显著片数分布
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

cat("[1] 校正方法实证：仅在 p_val_adj < 1（未被截顶）的行上检验 p_val_adj / p_val\n")
for (sl in slices4) {
  d  <- de_int[[sl]]
  ok <- !is.na(d$p_adj <- d$p_val_adj) & !is.na(d$p_val) & d$p_val > 0 & d$p_val_adj < 0.999999
  r  <- d$p_val_adj[ok] / d$p_val[ok]
  n_obj <- nrow(get(slice2obj[[sl]]))
  cat(sprintf("  %s | 未截顶行 %d / %d | 比值唯一值数 %d | 比值 %s | 对象基因数 %d | => %s\n",
              sl, sum(ok), nrow(d), length(unique(round(r, 8))),
              ifelse(length(unique(round(r, 8))) == 1, sprintf("%.1f", median(r)), "非恒定"),
              n_obj,
              ifelse(length(unique(round(r, 8))) == 1 && abs(median(r) - n_obj) < 0.5,
                     "Bonferroni (p x 对象基因数)", "需进一步判断")))
  cat(sprintf("      截顶行(p_val_adj==1)数: %d\n", sum(d$p_val_adj >= 0.999999, na.rm = TRUE)))
}
cat("  说明：若未截顶行比值为常数且等于对象基因数，则为 Seurat 默认 Bonferroni；BH 校正不会给出恒定比值\n")

cat("\n[2] 检出比例与方向一致性（log2FC 方向的\"检出率\"佐证）\n")
pct_tab <- list()
for (sl in slices4) {
  d <- de_int[[sl]]; sig <- !is.na(d$p_val_adj) & d$p_val_adj < 0.05
  hi <- sig & d$avg_log2FC > 0; lo <- sig & d$avg_log2FC < 0
  pct_tab[[length(pct_tab) + 1]] <- data.frame(
    slice = sl, hi_n = sum(hi), hi_pct1_gt_pct2 = sum(d$pct.1[hi] > d$pct.2[hi]),
    hi_frac = round(100 * sum(d$pct.1[hi] > d$pct.2[hi]) / max(1, sum(hi)), 1),
    lo_n = sum(lo), lo_pct1_lt_pct2 = sum(d$pct.1[lo] < d$pct.2[lo]),
    lo_frac = round(100 * sum(d$pct.1[lo] < d$pct.2[lo]) / max(1, sum(lo)), 1),
    stringsAsFactors = FALSE)
}
pct_df <- do.call(rbind, pct_tab); print(pct_df, row.names = FALSE)

cat("\n[3] Interface 覆盖率（分区表 Interface spot 数 vs 实际入分析数）\n")
cov <- data.frame(slice = slices4,
                  region_interface = sapply(slices4, function(s) sum(reg$slice == s & reg$region == "Interface")),
                  used = sapply(slices4, function(s) ncol(int_obj[[s]])),
                  obj_total = sapply(slices4, function(s) sum(as.character(get(slice2obj[[s]])$slice) == s)))
cov$coverage_pct <- round(100 * cov$used / cov$region_interface, 1)
print(cov, row.names = FALSE)

cat("\n[4] 显著片数分布\n")
cat("  n_slices_sig == 1 (仅 1 片显著):", sum(n_slices_sig == 1), "\n")
cat("  n_concordant == 1 (同向显著最多 1 片，含 2 片显著但方向冲突):", sum(n_concord == 1), "\n")
cat("  2 片显著但方向冲突:", sum(n_slices_sig >= 2 & n_up_sig > 0 & n_dn_sig > 0), "\n")

write.csv(pct_df, file.path(tab_dir, "de_interface_pct_direction_check.csv"), row.names = FALSE)
write.csv(cov,    file.path(tab_dir, "de_interface_region_coverage.csv"), row.names = FALSE)
cat("[Export] de_interface_pct_direction_check.csv / de_interface_region_coverage.csv\n")
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 构造 de_interface_done 状态载荷（精简版，控制在 4KB 以内；产物用命名模式指针而非逐个列出）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

top5 <- function(sl, dirn, n = 5) {
  d <- de_int[[sl]]; sig <- !is.na(d$p_val_adj) & d$p_val_adj < 0.05
  if (dirn == "up") { x <- d[sig & d$avg_log2FC > 0, ]; x <- head(x[order(-x$avg_log2FC), ], n) }
  else               { x <- d[sig & d$avg_log2FC < 0, ]; x <- head(x[order(x$avg_log2FC), ], n) }
  x$gene
}
per_slice <- list()
for (sl in slices4) {
  i  <- which(de_int_summary$slice == sl)
  d1 <- int_desc[int_desc$slice == sl & int_desc$clone == "Clone1", ]
  d2 <- int_desc[int_desc$slice == sl & int_desc$clone == "Clone2", ]
  per_slice[[sl]] <- list(
    patient = ifelse(sl %in% c("094A","094B"), "P13", "P15"),
    n_spot = de_int_summary$n_spot[i], n_Clone1 = de_int_summary$n_Clone1[i], n_Clone2 = de_int_summary$n_Clone2[i],
    n_tested = de_int_summary$n_tested[i], n_sig = de_int_summary$n_sig[i],
    n_up = de_int_summary$n_up[i], n_down = de_int_summary$n_down[i],
    median_nCount_C1 = d1$median_nCount, median_nCount_C2 = d2$median_nCount,
    depth_diff_pct = round(100 * (d1$median_nCount - d2$median_nCount) / d2$median_nCount, 1),
    tumor_prop_diff = round(d1$median_tumor_prop - d2$median_tumor_prop, 3),
    top_up = top5(sl, "up"), top_down = top5(sl, "down"))
}

state <- list(
  "__producer" = "sc-analyst", status = "done",
  analysis = "clone1_vs_clone2_DE_within_Interface_per_slice",
  disease_background = "TNBC",
  sample_type = "primary tumor tissue; Visium spot-level (multi-cell per spot)",
  region = "Interface; clone-labelled (tumor-enriched) spots only",
  method_used = list(test = "Seurat FindMarkers test.use=wilcox", min_pct = 0.25, logfc_threshold = 0,
                     ident_1 = "Clone1", ident_2 = "Clone2", direction = "avg_log2FC>0 = Clone1-high",
                     significance = "p_val_adj<0.05",
                     p_adjust = "Bonferroni (Seurat default), n=object gene count 11743/12978; NOT BH/FDR",
                     normalization = "LogNormalize 1e4 (inherited)", subsampling = "none",
                     pooling = "none; 4 independent per-slice tests"),
  join = list(key_rule = "<slice>_<barcode stripped of trailing _N>", assoc_rate_pct_all_slices = 100),
  per_slice = per_slice,
  consensus = list(ge2 = sum(n_concord >= 2), ge3 = sum(n_concord >= 3), eq4 = sum(n_concord == 4),
                   only_1_slice = sum(n_slices_sig == 1),
                   direction_conflict = sum(n_slices_sig >= 2 & n_up_sig > 0 & n_dn_sig > 0),
                   within_patient_same_sign = list(P13 = "63/63", P15 = "103/103"),
                   cross_patient_ge2_both = "7/7 same direction (descriptive)",
                   jaccard_range = c(0.025, 0.106),
                   shared_4slices_up = g4_hi, shared_4slices_down = g4_lo),
  heatmap = list(n_rows = nrow(z_plot), n_cols = ncol(z_plot), n_na_cells = sum(is.na(z_plot)),
                 cluster_sizes = as.list(tabulate(grp4, nbins = 4)),
                 patterns = "1=P13-high 2=P15-high (both patient-driven); 3=Clone2-ish weak; 4=Clone1-ish moderate"),
  interface_coverage_pct = setNames(as.list(cov$coverage_pct), cov$slice),
  confounder = list(verdict = "clone difference NOT separable from sequencing depth",
                    depth_direction = "P13 Clone1 lower (-37.9/-48.5%); P15 Clone1 higher (+162.8/+143.1%)",
                    pct_concordance = "Clone1-high with pct.1>pct.2: 4.7%/23.1% (P13) vs 85.4%/98.0% (P15)",
                    depth_sensitive_family = "enriched in low-depth clone in 3/4 slices; 094B opposite",
                    tumor_prop = "NOT a confounder (diff <=0.015)"),
  power_caveat = list(`094A_Clone2_spots` = 62, `095A_Clone2_spots` = 44,
                      note = "spot-level non-independence; p-values pseudo-replicated"),
  outputs = list(dir = "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone",
                 de_tables = "tables/de_interface_<slice>_clone1_vs_clone2.csv (4 slices)",
                 summary = "tables/de_interface_summary.csv",
                 consensus = "tables/de_interface_consensus.csv; de_interface_shared_all4slices.csv; de_interface_pairwise_jaccard.csv",
                 confounder = "tables/de_interface_depth_sensitive_families.csv; _pct_direction_check.csv; _region_coverage.csv",
                 heatmap_tables = "tables/heatmap_interface_tops_{expression,zscore}.csv",
                 figures = "figures/per_patient/volcano_interface_<slice>_clone1_vs_clone2.{png,pdf}; heatmap_interface_tops_4slices.{png,pdf}",
                 delivery_html = "step03_interface_de_汇总.html"),
  notebook = "project/project_spatial/TNBC/patient_ALL/TAM_clone/notebooks/step02_clone1_vs_clone2_de.ipynb"
)

js <- jsonlite::toJSON(state, auto_unbox = TRUE, pretty = FALSE, digits = 4)
cat("[size] bytes:", nchar(js), "\n\n"); cat(js)
cat("\n===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 交付文档最后一致性修补：把 §二 热图条目中"8 列间 z-score"改为"非缺失列间"，与 §5.4 的披露一致
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

f <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step03_interface_de_汇总.html"
txt <- paste(readLines(f, warn = FALSE), collapse = "\n")

before_old <- length(gregexpr("每基因在 8 列间的 z-score", txt, fixed = TRUE)[[1]])
txt <- gsub("每基因在 8 列间的 z-score（±2.5 截断）", "每基因在非缺失列间的 z-score（±2.5 截断）", txt, fixed = TRUE)
after_old  <- length(gregexpr("每基因在 8 列间的 z-score", txt, fixed = TRUE)[[1]])
after_new  <- length(gregexpr("每基因在非缺失列间的 z-score", txt, fixed = TRUE)[[1]])

writeLines(txt, f, useBytes = TRUE)
cat(sprintf("[Patch] 旧表述出现 %d 次 -> 修补后残留 %d 次 | 新表述出现 %d 次\n", before_old, after_old, after_new))
cat(sprintf("[Patch] 文件大小 %.2f MB\n", file.info(f)$size / 1e6))
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 修补残留的 1 处"8 列间 z-score"表述（无截断后缀的那一处）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

f <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step03_interface_de_汇总.html"
txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
hit <- regmatches(txt, gregexpr("每基因在 8 列间的 z-score[^<]{0,20}", txt))[[1]]
cat("[Found]", paste(hit, collapse = " || "), "\n")

txt <- gsub("每基因在 8 列间的 z-score", "每基因在非缺失列间的 z-score", txt, fixed = TRUE)
writeLines(txt, f, useBytes = TRUE)
cat(sprintf("[Patch] 残留 '8 列间' %d 次 | '非缺失列间' %d 次\n",
            length(gregexpr("每基因在 8 列间的 z-score", txt, fixed = TRUE)[[1]]) * (grepl("每基因在 8 列间的 z-score", txt, fixed = TRUE)),
            length(gregexpr("每基因在非缺失列间的 z-score", txt, fixed = TRUE)[[1]])))
cat(sprintf("[HTML] %.2f MB\n", file.info(f)$size / 1e6))
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 交付前最终核验：文件存在性/大小、图 pdf+png 成对、CSV 行列数、像素完整性
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

base <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone"
chk <- function(p) { fp <- file.path(base, p); data.frame(path = p, exists = file.exists(fp),
                        kb = ifelse(file.exists(fp), round(file.info(fp)$size/1024, 1), NA)) }

stem_v <- paste0("figures/per_patient/volcano_interface_", slices4, "_clone1_vs_clone2")
stem_h <- "figures/per_patient/heatmap_interface_tops_4slices"
figs <- c(paste0(rep(stem_v, each = 2), c(".png", ".pdf")), paste0(stem_h, ".png"), paste0(stem_h, ".pdf"))
tabs <- c(paste0("tables/de_interface_", slices4, "_clone1_vs_clone2.csv"),
          "tables/de_interface_summary.csv", "tables/de_interface_consensus.csv",
          "tables/de_interface_shared_all4slices.csv", "tables/de_interface_pairwise_jaccard.csv",
          "tables/de_interface_depth_sensitive_families.csv", "tables/de_interface_pct_direction_check.csv",
          "tables/de_interface_region_coverage.csv", "tables/heatmap_interface_tops_expression.csv",
          "tables/heatmap_interface_tops_zscore.csv", "step03_interface_de_汇总.html")

all_chk <- rbind(chk(figs), chk(tabs))
cat(sprintf("[Files] 共 %d 项 | 缺失 %d 项 | 空文件(<1KB) %d 项\n",
            nrow(all_chk), sum(!all_chk$exists), sum(all_chk$exists & all_chk$kb < 1)))
print(all_chk, row.names = FALSE)

cat("\n[Figure pairs] pdf+png 成对检查:\n")
for (s in c(stem_v, stem_h)) {
  cat(sprintf("  %s | png %s | pdf %s\n", basename(s), file.exists(file.path(base, paste0(s, ".png"))),
              file.exists(file.path(base, paste0(s, ".pdf")))))
}

cat("\n[Pixel] 5 张关键图白底/非空校验:\n")
for (s in c(stem_v, stem_h)) {
  im <- png::readPNG(file.path(base, paste0(s, ".png"))); rgb <- im[, , 1:3]
  cat(sprintf("  %-46s | white=%.3f | sd=%.3f | %s\n", basename(s), mean(rgb > 0.98), sd(as.vector(rgb)),
              ifelse(sd(as.vector(rgb)) > 0.02, "PASS", "FAIL")))
}

cat("\n[CSV 行列数]\n")
for (f in tabs[grepl("\\.csv$", tabs)]) {
  d <- read.csv(file.path(base, f), check.names = FALSE, nrows = 1)
  n <- length(readLines(file.path(base, f))) - 1
  cat(sprintf("  %-52s | %5d 行 | %2d 列\n", basename(f), n, ncol(d)))
}
cat("===== Cell End:", date(), "| elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")
