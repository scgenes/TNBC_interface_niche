# =============================================================================
# TNBC interface-niche code repository
# Original source: FigureAI/results/cellchat_by_region/scripts/01_cellchat_by_region.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

#!/usr/bin/env Rscript
# =============================================================================
# CellChat 按区推断 C1Q TAM <-> MMP11 CAF 通讯（Visium 空间数据，逐样本 x 分区）
#
# 运行（cellchat 环境，与 notebook 完全同参）：
#   /path/to/envs/cellchat/bin/Rscript scripts/03_communication/04_cellchat_CAF_C1Q_by_region.R
#
# 设计要点
#   * 分析单元：12 样本 x 3 分区，逐单元单独建 CellChat 网络，再按区汇总（与 COMMOT 一致）
#   * 空间模式：distance.use = TRUE，interaction.range = 500（与 COMMOT dis_thr=500 同一坐标系）
#   * 参数口径：CellChatDB.human；computeCommunProb(type="triMean", nboot=100, seed.use=1)；
#               filterCommunication(min.cells=10)
#   * spot 归属（主口径，用户指定）：C1Q TAM >= 2%、MMP11 CAF >= 5%；
#       两组阈值同时满足 -> 归给比例较高的一组；都不满足 -> "Other"（保留，不剔除）
#   * 敏感性口径：两组阈值同时满足 -> 归给更稀有的 C1Q TAM 组（与主口径相反），仅用于检验
#       归属规则对结论的影响，不参与主结论
#   * 缓存：data/cellchat_by_region_cache.rds 存在时只重画图，不重算网络
# =============================================================================

suppressPackageStartupMessages({
  library(CellChat); library(Seurat); library(ggplot2)
  library(dplyr); library(tidyr); library(patchwork)
})

set.seed(1)

# ---- 路径与参数 ----
OUT     <- "/path/to/project/project_spatial/FigureAI/results/cellchat_by_region"
DIR_TAB <- file.path(OUT, "tables"); DIR_FIG <- file.path(OUT, "figures")
DIR_DAT <- file.path(OUT, "data")
for (d in c(OUT, DIR_TAB, DIR_FIG, DIR_DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

IN_OBJ <- "/path/to/data/other_ans_project/三阴性乳腺癌分析/data/integrated_seurat_12.rds"
IN_SUB <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/subcluster_prop.csv.gz"
IN_RA  <- "/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment.csv"

C1Q_COL <- "C1Q TAMs"; MMP_COL <- "MMP11 CAFs"
THR_C1Q <- 0.02; THR_MMP <- 0.05
MIN_GROUP <- 10          # 与 filterCommunication(min.cells) 对齐
INTERACTION_RANGE <- 500; SCALE_DISTANCE <- 0.01
NBoot <- 100; SEED_USE <- 1L; MIN_CELLS <- 10
GP_MAIN <- c("C1Q TAM", "MMP11 CAF", "Other")
REGIONS <- c("Tumor Core", "Interface", "Stroma")
FOCUS <- c("C1Q TAM -> MMP11 CAF", "MMP11 CAF -> C1Q TAM")

# ---- 1. 输入与分组 ----
obj <- readRDS(IN_OBJ)
sub <- read.csv(gzfile(IN_SUB), check.names = FALSE)
ra  <- read.csv(IN_RA, check.names = FALSE)

stopifnot(identical(colnames(obj), sub$barcode))
sp <- data.frame(barcode = colnames(obj), sample = obj$orig.ident,
                 c1q = sub[[C1Q_COL]], mmp11 = sub[[MMP_COL]], stringsAsFactors = FALSE)
sp <- merge(sp, ra[, c("barcode", "region", "x", "y")], by = "barcode", all.x = TRUE, sort = FALSE)
sp <- sp[match(colnames(obj), sp$barcode), ]
stopifnot(!any(is.na(sp$region)))

q_c1q <- sp$c1q >= THR_C1Q; q_mmp <- sp$mmp11 >= THR_MMP
sp$group <- "Other"
sp$group[q_c1q & !q_mmp] <- "C1Q TAM"
sp$group[q_mmp & !q_c1q] <- "MMP11 CAF"
dual <- q_c1q & q_mmp
sp$group[dual] <- ifelse(sp$c1q[dual] >= sp$mmp11[dual], "C1Q TAM", "MMP11 CAF")  # 主口径
sp$group_alt <- "Other"
sp$group_alt[q_c1q] <- "C1Q TAM"
sp$group_alt[q_mmp & !q_c1q] <- "MMP11 CAF"                                       # 敏感性口径

# ---- 2. 可估性 ----
units <- expand.grid(sample = sort(unique(sp$sample)), region = REGIONS, stringsAsFactors = FALSE)
est <- do.call(rbind, lapply(seq_len(nrow(units)), function(i) {
  d <- sp[sp$sample == units$sample[i] & sp$region == units$region[i], ]
  tm <- table(factor(d$group, levels = GP_MAIN)); ta <- table(factor(d$group_alt, levels = GP_MAIN))
  data.frame(sample = units$sample[i], region = units$region[i], n_spots = nrow(d),
             main_C1Q = as.integer(tm["C1Q TAM"]), main_MMP11 = as.integer(tm["MMP11 CAF"]),
             main_Other = as.integer(tm["Other"]),
             alt_C1Q = as.integer(ta["C1Q TAM"]), alt_MMP11 = as.integer(ta["MMP11 CAF"]),
             alt_Other = as.integer(ta["Other"]), stringsAsFactors = FALSE)
}))
est <- est[est$n_spots > 0, ]
est$main_estimable <- est$main_C1Q >= MIN_GROUP & est$main_MMP11 >= MIN_GROUP
est$alt_estimable  <- est$alt_C1Q  >= MIN_GROUP & est$alt_MMP11  >= MIN_GROUP
write.csv(est, file.path(DIR_TAB, "cellchat_region_estimability.csv"), row.names = FALSE)

# ---- 3. 单单元 CellChat ----
DAT_SCT <- GetAssayData(obj, assay = "SCT", layer = "data")
data("CellChatDB.human", package = "CellChat"); CCDB <- CellChatDB.human

build_cc_unit <- function(d, group_levels, spatial = TRUE) {
  bc <- as.character(d$barcode)
  lab <- droplevels(factor(d$group, levels = group_levels))
  meta <- data.frame(labels = lab, samples = factor(rep("unit", length(bc))), row.names = bc)
  coords <- cbind(x = as.numeric(d$x), y = as.numeric(d$y)); rownames(coords) <- bc
  cc <- createCellChat(object = DAT_SCT[, bc, drop = FALSE], meta = meta, group.by = "labels",
                       datatype = if (spatial) "spatial" else "RNA", coordinates = coords,
                       spatial.factors = data.frame(ratio = 1, tol = 0))
  cc@DB <- CCDB
  cc <- setIdent(cc, ident.use = "labels")
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)   # cellchat 环境无 presto
  cc <- identifyOverExpressedInteractions(cc)
  cc <- if (spatial) {
    computeCommunProb(cc, type = "triMean", trim = 0.1, population.size = FALSE, raw.use = TRUE,
                      distance.use = TRUE, interaction.range = INTERACTION_RANGE,
                      scale.distance = SCALE_DISTANCE, k.min = 10, contact.dependent = TRUE,
                      contact.range = INTERACTION_RANGE, nboot = NBoot, seed.use = SEED_USE)
  } else {
    computeCommunProb(cc, type = "triMean", trim = 0.1, population.size = FALSE,
                      raw.use = TRUE, nboot = NBoot, seed.use = SEED_USE)
  }
  cc <- filterCommunication(cc, min.cells = MIN_CELLS)
  cc <- computeCommunProbPathway(cc)
  aggregateNet(cc)
}

run_batch <- function(units_df, rule) {
  gcol <- if (rule == "main") "group" else "group_alt"
  E <- list(); L <- list(); P <- list(); OBJ <- list()
  for (i in seq_len(nrow(units_df))) {
    s <- units_df$sample[i]; rg <- units_df$region[i]; uid <- paste(s, rg, sep = "|")
    d <- sp[sp$sample == s & sp$region == rg, ]
    d2 <- data.frame(barcode = d$barcode, group = d[[gcol]], x = d$x, y = d$y, stringsAsFactors = FALSE)
    .logf <- tempfile(); .con <- file(.logf, "wt")
    sink(.con, type = "output"); sink(.con, type = "message")
    t0 <- Sys.time(); ok <- TRUE; err <- NA_character_
    cc <- tryCatch(build_cc_unit(d2, GP_MAIN), error = function(e) { ok <<- FALSE; err <<- conditionMessage(e); NULL })
    sink(type = "message"); sink(type = "output"); close(.con)
    rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!ok) { cat(sprintf("[FAIL] %s | %s\n", uid, err)); next }
    w <- if (!is.null(cc@net$weight)) cc@net$weight else apply(cc@net$prob, c(1, 2), sum)
    for (src in rownames(w)) for (tgt in colnames(w)) {
      E[[length(E) + 1]] <- data.frame(rule = rule, unit = uid, sample = s, region = rg,
        source = src, target = tgt, weight = w[src, tgt], n_LR_sig = cc@net$count[src, tgt],
        n_cells_source = sum(cc@idents == src), n_cells_target = sum(cc@idents == tgt),
        n_spots = nrow(d2), spatial_mode = TRUE, runtime_s = round(rt, 1), stringsAsFactors = FALSE)
    }
    lr <- tryCatch(subsetCommunication(cc, slot.name = "net"), error = function(e) NULL)
    if (!is.null(lr)) {
      lr <- lr[lr$pval < 0.05 & ((lr$source == "C1Q TAM" & lr$target == "MMP11 CAF") |
                                 (lr$source == "MMP11 CAF" & lr$target == "C1Q TAM")), ]
      if (nrow(lr) > 0) { lr$rule <- rule; lr$unit <- uid; lr$sample <- s; lr$region <- rg
                          L[[length(L) + 1]] <- lr }
    }
    np <- cc@netP$prob
    for (src in dimnames(np)[[1]]) for (tgt in dimnames(np)[[2]]) {
      if (!((src == "C1Q TAM" & tgt == "MMP11 CAF") | (src == "MMP11 CAF" & tgt == "C1Q TAM"))) next
      pv <- np[src, tgt, ]; pv <- pv[pv > 0]
      if (length(pv)) P[[length(P) + 1]] <- data.frame(rule = rule, unit = uid, sample = s, region = rg,
        source = src, target = tgt, pathway = names(pv), prob = as.numeric(pv), stringsAsFactors = FALSE)
    }
    OBJ[[uid]] <- list(unit = uid, sample = s, region = rg, rule = rule, weight = w,
                       count = cc@net$count, groups = levels(cc@idents),
                       n_cells = as.integer(table(cc@idents)), n_spots = nrow(d2),
                       prob = cc@net$prob, pval = cc@net$pval, pathways = cc@netP$pathways,
                       LR_sig = if (is.null(lr)) NULL else lr, parameter = cc@options$parameter)
    cat(sprintf("[OK] %s | %s | %.1fs\n", uid, rule, rt))
  }
  list(edges = if (length(E)) do.call(rbind, E), lr = if (length(L)) do.call(rbind, L),
       path = if (length(P)) do.call(rbind, P), objs = OBJ)
}

f_cache <- file.path(DIR_DAT, "cellchat_by_region_cache.rds")
if (file.exists(f_cache)) {
  cat("[Cache] 命中，跳过网络计算：", f_cache, "\n")
  cache <- readRDS(f_cache)
} else {
  res_main <- run_batch(est[est$main_estimable, c("sample", "region")], "main")
  res_alt  <- run_batch(est[est$alt_estimable,  c("sample", "region")], "alt")
  # 单元键用 "rule|unit"：主口径与敏感性口径的单元 id 相同，若直接按 unit 命名会互相覆盖
  u_all <- c(res_main$objs, res_alt$objs)
  names(u_all) <- vapply(u_all, function(z) paste(z$rule, z$unit, sep = "|"), character(1))
  cache <- list(units = u_all,
                edges = rbind(res_main$edges, res_alt$edges),
                lr = rbind(res_main$lr, res_alt$lr), path = rbind(res_main$path, res_alt$path),
                estimability = est)
  saveRDS(cache, f_cache, compress = "xz")
}

# ---- 4. 表 ----
dir_lab <- function(s, t) paste0(s, " -> ", t)
E_all <- cache$edges;  E_all$direction <- dir_lab(E_all$source, E_all$target)
L_all <- cache$lr;     L_all$direction <- dir_lab(L_all$source, L_all$target)
P_all <- cache$path;   P_all$direction <- dir_lab(P_all$source, P_all$target)
E_f <- E_all[E_all$direction %in% FOCUS, ]; L_f <- L_all[L_all$direction %in% FOCUS, ]
P_f <- P_all[P_all$direction %in% FOCUS, ]

desc <- E_f %>% group_by(rule, region, direction) %>%
  summarise(n_units = n(), weight_mean = mean(weight), weight_median = median(weight),
            weight_sd = sd(weight), weight_min = min(weight), weight_max = max(weight),
            nLR_mean = mean(n_LR_sig), nLR_median = median(n_LR_sig), .groups = "drop")
write.csv(as.data.frame(desc), file.path(DIR_TAB, "cellchat_focus_edge_by_region.csv"), row.names = FALSE)

paired_test <- function(df, d, rule) {
  a <- df[df$rule == rule & df$direction == d & df$region == "Tumor Core", ]
  b <- df[df$rule == rule & df$direction == d & df$region == "Interface", ]
  s <- intersect(a$sample, b$sample)
  if (length(s) < 3) return(data.frame(rule = rule, direction = d, n_pairs = length(s), p = NA))
  x <- a$weight[match(s, a$sample)]; y <- b$weight[match(s, b$sample)]
  data.frame(rule = rule, direction = d, n_pairs = length(s), median_core = median(x),
             median_interface = median(y), median_diff = median(x - y),
             V = as.numeric(wilcox.test(x, y, paired = TRUE, exact = TRUE)$statistic),
             p = wilcox.test(x, y, paired = TRUE, exact = TRUE)$p.value)
}
pair_tab <- bind_rows(lapply(c("main", "alt"), function(r)
  bind_rows(lapply(FOCUS, function(d) paired_test(E_f, d, r)))))
write.csv(pair_tab, file.path(DIR_TAB, "cellchat_core_vs_interface_paired_test.csv"), row.names = FALSE)

asym <- E_f %>% select(rule, region, sample, direction, weight) %>%
  pivot_wider(names_from = direction, values_from = weight) %>% group_by(rule, region) %>%
  summarise(n = sum(!is.na(.data[[FOCUS[1]]]) & !is.na(.data[[FOCUS[2]]])),
            median_A_to_B = median(.data[[FOCUS[1]]], na.rm = TRUE),
            median_B_to_A = median(.data[[FOCUS[2]]], na.rm = TRUE),
            ratio = median_B_to_A / median_A_to_B,
            p_signedrank = tryCatch(wilcox.test(.data[[FOCUS[2]]], .data[[FOCUS[1]]],
                                                paired = TRUE)$p.value, error = function(e) NA),
            .groups = "drop")
write.csv(as.data.frame(asym), file.path(DIR_TAB, "cellchat_direction_asymmetry.csv"), row.names = FALSE)

n_units_tab <- E_f %>% distinct(rule, region, unit) %>% count(rule, region, name = "n_units")
lr_rec <- L_f %>% group_by(rule, region, direction, ligand, receptor, pathway_name) %>%
  summarise(n_sig_units = n_distinct(unit), mean_prob = mean(prob), median_prob = median(prob),
            .groups = "drop") %>%
  left_join(n_units_tab, by = c("rule", "region")) %>% mutate(frac_units = n_sig_units / n_units)
pw_rec <- P_f %>% group_by(rule, region, direction, pathway) %>%
  summarise(n_units = n_distinct(unit), mean_prob = mean(prob), median_prob = median(prob), .groups = "drop")

write.csv(E_all, file.path(DIR_TAB, "cellchat_edges_by_unit.csv"), row.names = FALSE)
write.csv(L_all, file.path(DIR_TAB, "cellchat_significant_LR_by_unit.csv"), row.names = FALSE)
write.csv(P_all, file.path(DIR_TAB, "cellchat_pathway_by_unit.csv"), row.names = FALSE)
write.csv(as.data.frame(lr_rec), file.path(DIR_TAB, "cellchat_LR_recurrence_by_region.csv"), row.names = FALSE)
write.csv(as.data.frame(pw_rec), file.path(DIR_TAB, "cellchat_pathway_recurrence_by_region.csv"), row.names = FALSE)

cat("[Done] tables written to", DIR_TAB, "\n")
cat("[Note] 图与时序统计见 notebook：notebooks/cellchat_by_region.ipynb（本脚本复现网络计算与核心表）\n")
