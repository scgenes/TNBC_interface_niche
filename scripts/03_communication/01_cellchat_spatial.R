# =============================================================================
# TNBC interface-niche code repository
# Original source: project_singlecell/TNBC_result/final_figures/V3/rerun_results/fig4_barrier/cellchat/01_cellchat_compute.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

## ============================================================
## TNBC 图4 细胞通讯 — 步骤1：CellChat 计算（全局配体-受体扫描）
## 层级A：Celltype_V2（大群，~9类）→ 全局通讯网络
## 层级B：Celltype_sub（亚群，~35类，含 7 TAM + 6 CAF）→ 解析两条故事轴的主导亚群
## 完整可复现代码（createCellChat → computeCommunProb → computeCommunProbPathway → aggregateNet）
## ============================================================
suppressPackageStartupMessages({library(CellChat); library(Seurat); library(dplyr)})

OUT <- "/path/to/project/project_singlecell/TNBC_result/results/fig4_cellchat"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "01_cellchat_compute.log"), split = TRUE)

run_cellchat <- function(obj, group.by, tag, population.size = FALSE, nboot = 100) {
  cat("\n############################################################\n")
  cat(sprintf("## CellChat: group.by=%s (tag=%s)\n", group.by, tag))
  cat("############################################################\n")
  cat("[Groups]:\n"); print(sort(table(obj@meta.data[[group.by]]), decreasing = TRUE))

  t0 <- Sys.time()
  cc <- createCellChat(object = obj, group.by = group.by, assay = "RNA")
  cc@DB <- CellChatDB.human
  cat(sprintf("[DB] %d interactions (full CellChatDB.human: Secreted + Cell-Cell Contact + ECM-Receptor)\n",
              nrow(cc@DB$interaction)))

  cc <- subsetData(cc)                        # keep signaling genes present in DB
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)

  t1 <- Sys.time()
  cc <- computeCommunProb(cc, type = "triMean", population.size = population.size,
                          nboot = nboot, seed.use = 1)
  cat(sprintf("[computeCommunProb] %.1f min (nboot=%d, population.size=%s)\n",
              as.numeric(difftime(Sys.time(), t1, units = "mins")), nboot, population.size))

  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)

  saveRDS(cc, file.path(OUT, sprintf("cellchat_%s.rds", tag)))
  cat(sprintf("[Save] cellchat_%s.rds (total %.1f min)\n", tag,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  cc
}

write_tables <- function(cc, tag) {
  ## 全量 LR 表
  comm <- subsetCommunication(cc)
  comm$pval <- as.numeric(comm$pval)
  write.csv(comm, file.path(OUT, sprintf("%s_all_interactions.csv", tag)), row.names = FALSE)
  ## 显著 (permutation p < 0.05)
  sig <- comm[!is.na(comm$pval) & comm$pval < 0.05, ]
  write.csv(sig, file.path(OUT, sprintf("%s_significant_interactions.csv", tag)), row.names = FALSE)
  cat(sprintf("[%s] all interactions: %d | significant (p<0.05): %d\n", tag, nrow(comm), nrow(sig)))

  ## 网络 count / weight（细胞类型对级别）
  write.csv(as.data.frame(cc@net$count),  file.path(OUT, sprintf("%s_net_count.csv", tag)))
  write.csv(as.data.frame(cc@net$weight), file.path(OUT, sprintf("%s_net_weight.csv", tag)))

  ## 通路级 netP 汇总
  netP <- cc@netP$prob
  pw <- dimnames(netP)[[3]]
  path_df <- data.frame(pathway = pw, total_prob = as.numeric(apply(netP, 3, sum)))
  path_df <- path_df[order(-path_df$total_prob), ]
  write.csv(path_df, file.path(OUT, sprintf("%s_pathway_summary.csv", tag)), row.names = FALSE)
  cat(sprintf("[%s] pathways: %d | top 15:\n", tag, length(pw)))
  print(head(path_df, 15) %>% select(pathway, total_prob))
  invisible(path_df)
}

## ---- 载入 lean 对象 ----
t0 <- Sys.time()
obj <- readRDS(file.path(OUT, "sc_lean.rds"))
cat(sprintf("[Load sc_lean.rds] %.1f min | %d cells x %d genes\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), ncol(obj), nrow(obj)))

## ============================================================
## 层级A：大群 Celltype_V2
## ============================================================
cc_major <- tryCatch(run_cellchat(obj, "Celltype_V2", "major"),
                     error = function(e) { cat("[ERROR major]", conditionMessage(e), "\n"); NULL })
if (!is.null(cc_major)) {
  write_tables(cc_major, "major")
  ## 全局网络图：count / weight 热图
  png(file.path(OUT, "fig4_global_net_count.png"), width = 2000, height = 1700, res = 150)
  tryCatch(netVisual_heatmap(cc_major, measure = "count", color.heatmap = "Reds",
                             title.name = "Number of interactions (Celltype_V2)"),
           error = function(e) cat("[plot err count]", conditionMessage(e), "\n"))
  dev.off()
  png(file.path(OUT, "fig4_global_net_weight.png"), width = 2000, height = 1700, res = 150)
  tryCatch(netVisual_heatmap(cc_major, measure = "weight", color.heatmap = "Blues",
                             title.name = "Interaction strength (Celltype_V2)"),
           error = function(e) cat("[plot err weight]", conditionMessage(e), "\n"))
  dev.off()
  rm(cc_major); gc()
}

## ============================================================
## 层级B：亚群 Celltype_sub
## ============================================================
cc_sub <- tryCatch(run_cellchat(obj, "Celltype_sub", "subtype"),
                   error = function(e) { cat("[ERROR subtype]", conditionMessage(e), "\n"); NULL })
if (!is.null(cc_sub)) {
  write_tables(cc_sub, "subtype")
  rm(cc_sub); gc()
}

sink()
cat("[Done] 01_cellchat_compute\n")
