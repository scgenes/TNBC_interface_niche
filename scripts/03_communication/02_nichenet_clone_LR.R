# =============================================================================
# TNBC interface-niche code repository
# Original source: project_spatial/TNBC/patient_ALL/TAM_clone/notebooks/step03_clone_LR.ipynb (converted to R)
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

cat("===== Cell Start:", date(), "=====\n")
# 环境探针：确认 R 版本、库路径、qs/data.table 可用性
cat("R:", R.version.string, "\n")
cat("libPaths:\n"); print(.libPaths())
cat("qs:", requireNamespace("qs", quietly = TRUE),
    "| data.table:", requireNamespace("data.table", quietly = TRUE),
    "| jsonlite:", requireNamespace("jsonlite", quietly = TRUE), "\n")
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 只读探查：region 文本表与两张 gz 表的表头/规模（gz 统一经 zcat 命令通道读取，避开 R.utils 依赖）
p_reg  <- "/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment_extended.csv"
p_sub  <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/subcluster_prop.csv.gz"
p_expr <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/lr_expr.csv.gz"
for (p in c(p_reg, p_sub, p_expr)) cat(sprintf("%s\n  exists=%s | size=%.1f MB\n", p, file.exists(p), file.info(p)$size/1e6))

r <- data.table::fread(p_reg, nrows = 3)
cat("\n[region] ncol =", ncol(r), "\n cols =", paste(names(r), collapse = " | "), "\n")
print(r[, 1:8, with = FALSE])

s <- data.table::fread(cmd = paste("zcat", shQuote(p_sub)), nrows = 3)
cat("\n[subcluster_prop] ncol =", ncol(s), "\n cols =", paste(names(s), collapse = " | "), "\n")
print(s[, 1:4, with = FALSE])

e <- data.table::fread(cmd = paste("zcat", shQuote(p_expr)), nrows = 3)
cat("\n[lr_expr] ncol =", ncol(e), "\n cols =", paste(names(e), collapse = " | "), "\n")
print(e[, 1:6, with = FALSE])
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 探查两张 gz 表结构（用 zcat 命令通道绕开 R.utils 依赖），并确认 region 表的 sample→slice 与 barcode 去后缀规则
p_sub  <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/subcluster_prop.csv.gz"
p_expr <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/lr_expr.csv.gz"
p_reg  <- "/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment_extended.csv"

cat("\n--- subcluster_prop.csv.gz ---\n")
s <- data.table::fread(cmd = paste("zcat", shQuote(p_sub)), nrows = 3)
cat("ncol =", ncol(s), "\n")
cat("cols (all) =", paste(names(s), collapse = " | "), "\n")
print(s[, 1:min(6, ncol(s)), with = FALSE])

cat("\n--- lr_expr.csv.gz ---\n")
e <- data.table::fread(cmd = paste("zcat", shQuote(p_expr)), nrows = 3)
cat("ncol =", ncol(e), "\n")
cat("cols =", paste(names(e), collapse = " | "), "\n")
print(e[, 1:min(8, ncol(e)), with = FALSE])

cat("\n--- region: sample -> slice, barcode suffix check ---\n")
r <- data.table::fread(p_reg, select = c("barcode", "sample", "region"))
r[, slice := sub("^GSM[0-9]+_", "", sample)]
r[, bc_nosuf := sub("_[0-9]+$", "", barcode)]
cat("unique slices:", paste(sort(unique(r$slice)), collapse = ", "), "\n")
cat("suffix example per slice:\n")
print(unique(r[, .(slice, suffix = sub("^.*_", "", barcode))])[order(slice, suffix)][1:20])
cat("dup key (slice+bc_nosuf):", sum(duplicated(r[, paste0(slice, "_", bc_nosuf)])), "\n")
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 核对：region 表 sample→slice 是否唯一映射、suffix 是否 1:1、三张 RCTD 表的 barcode suffix 分布与总行数
p_reg  <- "/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment_extended.csv"
p_sub  <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/subcluster_prop.csv.gz"
p_expr <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/lr_expr.csv.gz"

r <- data.table::fread(p_reg, select = c("barcode", "sample", "region"))
cat("region nrow =", nrow(r), "\n")
cat("region: NA sample =", sum(is.na(r$sample) | r$sample == ""), "| NA barcode =", sum(is.na(r$barcode) | r$barcode == ""), "\n")
r[, slice := sub("^GSM[0-9]+_", "", sample)]
cat("\n--- sample -> slice mapping (unique pairs) ---\n")
print(unique(r[, .(sample, slice)])[order(slice)])
cat("每个 slice 对应的 suffix 数:\n")
r[, suffix := sub("^.*_", "", barcode)]
print(r[, .(n = .N, suffixes = paste(sort(unique(suffix)), collapse = ",")), by = slice][order(slice)])

cat("\n--- subcluster_prop: 行数 + suffix 分布 ---\n")
s_key <- data.table::fread(cmd = paste("zcat", shQuote(p_sub), "| cut -d, -f1"), header = TRUE)
cat("subcluster_prop nrow =", nrow(s_key), "| 首行值 =\n"); print(head(s_key, 3))
s_key[, suffix := sub("^.*_", "", barcode)]
print(s_key[, .N, by = suffix][order(as.integer(suffix))])

cat("\n--- lr_expr: 行数 + suffix 分布 ---\n")
e_key <- data.table::fread(cmd = paste("zcat", shQuote(p_expr), "| cut -d, -f1"), header = TRUE)
cat("lr_expr nrow =", nrow(e_key), "| 首行值 =\n"); print(head(e_key, 3))
e_key[, suffix := sub("^.*_", "", barcode)]
print(e_key[, .N, by = suffix][order(as.integer(suffix))])
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 核对：克隆 spot 是否全部落入 RCTD 三表的 spot 全集（逐切片），并统计 RCTD 各切片 spot 数
p_reg  <- "/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment_extended.csv"
p_sub  <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/subcluster_prop.csv.gz"
p_expr <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/lr_expr.csv.gz"
c13 <- data.table::fread("/path/to/project/project_spatial/TNBC/patient_ALL/clone/P13/clone_membership.csv")
c15 <- data.table::fread("/path/to/project/project_spatial/TNBC/patient_ALL/clone/P15/clone_membership.csv")
cat("P13 rows =", nrow(c13), "| slices:", paste(names(table(c13$slice)), table(c13$slice), collapse = " / "), "\n")
cat("P15 rows =", nrow(c15), "| slices:", paste(names(table(c15$slice)), table(c15$slice), collapse = " / "), "\n")
cat("P13 dup cell key:", sum(duplicated(c13$cell)), "| P15 dup cell key:", sum(duplicated(c15$cell)), "\n")
cat("P13 clones:", paste(sort(unique(c13$clone)), collapse = ","),
    "| P15 clones:", paste(sort(unique(c15$clone)), collapse = ","), "\n")

suf2slice <- c("1"="092A","2"="092B","3"="093A","4"="093B","5"="093C","6"="093D",
               "7"="094A","8"="094B","9"="094C","10"="094D","11"="095A","12"="095B")
mk_key <- function(bc) paste0(suf2slice[sub("^.*_", "", bc)], "_", sub("_[0-9]+$", "", bc))
r_key <- mk_key(data.table::fread(p_reg, select = "barcode")$barcode)
s_key <- mk_key(data.table::fread(cmd = paste("zcat", shQuote(p_sub), "| cut -d, -f1"), header = TRUE)$barcode)
e_key <- mk_key(data.table::fread(cmd = paste("zcat", shQuote(p_expr), "| cut -d, -f1"), header = TRUE)$barcode)
cat("\nRCTD key dup:", sum(duplicated(r_key)), "| sub dup:", sum(duplicated(s_key)),
    "| expr dup:", sum(duplicated(e_key)), "\n")
cat("sub vs expr key 集合一致:", setequal(s_key, e_key), "| region vs sub 集合一致:", setequal(r_key, s_key), "\n")

cat("\n--- 逐切片重叠（clone -> region / sub / expr）---\n")
for (nm in c("P13", "P15")) {
  cc <- if (nm == "P13") c13 else c15
  for (sl in sort(unique(cc$slice))) {
    k <- cc$cell[cc$slice == sl]
    cat(sprintf("%s %s | clone=%d | in_region=%d (%.2f%%) | in_sub=%d (%.2f%%) | in_expr=%d (%.2f%%)\n",
                nm, sl, length(k), sum(k %in% r_key), 100*mean(k %in% r_key),
                sum(k %in% s_key), 100*mean(k %in% s_key), sum(k %in% e_key), 100*mean(k %in% e_key)))
  }
}
cat("\n--- RCTD 全集逐切片 spot 数 ---\n")
print(data.table::data.table(slice = sub("_.*$", "", r_key))[, .N, by = slice][order(slice)])
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 诊断：非克隆 spot（RCTD 有、克隆表无）落在哪些区域？决定底表 universe 是否必须保留它们
p_reg  <- "/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment_extended.csv"
p_sub  <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/subcluster_prop.csv.gz"
c13 <- data.table::fread("/path/to/project/project_spatial/TNBC/patient_ALL/clone/P13/clone_membership.csv")
c15 <- data.table::fread("/path/to/project/project_spatial/TNBC/patient_ALL/clone/P15/clone_membership.csv")
suf2slice <- c("1"="092A","2"="092B","3"="093A","4"="093B","5"="093C","6"="093D",
               "7"="094A","8"="094B","9"="094C","10"="094D","11"="095A","12"="095B")

reg2 <- data.table::fread(p_reg)
reg2[, slice := sub("^GSM[0-9]+_", "", sample)]
reg2[, cell  := paste0(slice, "_", sub("_[0-9]+$", "", barcode))]
cat("region cell dup:", sum(duplicated(reg2$cell)), "\n")

c1q <- data.table::fread(cmd = paste("zcat", shQuote(p_sub), "| cut -d, -f1,4"), header = TRUE)
data.table::setnames(c1q, c("barcode", "C1Q_TAM_prop"))
c1q[, cell := paste0(suf2slice[sub("^.*_", "", barcode)], "_", sub("_[0-9]+$", "", barcode))]

clone_use <- rbind(c13[, .(cell, clone)], c15[, .(cell, clone)])
sl_use <- c("094A","094B","095A","095B")
d <- merge(reg2[slice %in% sl_use, .(cell, slice, region)], c1q[, .(cell, C1Q_TAM_prop)], by = "cell", all.x = TRUE)
d <- merge(d, clone_use, by = "cell", all.x = TRUE)
d[, is_clone := !is.na(clone)]
cat("合并后行数:", nrow(d), "| C1Q_TAM_prop 缺失:", sum(is.na(d$C1Q_TAM_prop)), "\n")
cat("\n--- 病人两切片内：克隆 / 非克隆 spot 的区域构成与 C1Q TAM 阳性数 ---\n")
print(d[, .(total = .N, clone = sum(is_clone), nonclone = sum(!is_clone),
            c1q_pos = sum(C1Q_TAM_prop >= 0.02),
            c1q_pos_nonclone = sum(C1Q_TAM_prop >= 0.02 & !is_clone)),
        by = .(slice, region)][order(slice, region)])
cat("\n非克隆 spot 区域合计：\n")
print(d[is_clone == FALSE, .N, by = region][order(-N)])
cat("非克隆 spot 中的 C1Q TAM 阳性合计 =", d[is_clone == FALSE & C1Q_TAM_prop >= 0.02, .N], "\n")
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 正式构建：四输入按 <slice>_<barcode> 对齐，产出 P13/P15 病人内两切片合并底表 + 可估性格子表
suppressPackageStartupMessages(library(data.table))

## 0) 路径与输出目录
p13_path <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P13/clone_membership.csv"
p15_path <- "/path/to/project/project_spatial/TNBC/patient_ALL/clone/P15/clone_membership.csv"
p_reg    <- "/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment_extended.csv"
p_sub    <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/subcluster_prop.csv.gz"
p_expr   <- "/path/to/project/project_singlecell/TNBC_result/figure18/intermediates/lr_expr.csv.gz"
out_dir  <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step03_clone_LR"
dir.create(file.path(out_dir, "data", "checkpoint"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

## 1) 读入四个输入
c13 <- fread(p13_path); c15 <- fread(p15_path)
reg <- fread(p_reg); sub <- fread(cmd = paste("zcat", shQuote(p_sub)))
expr <- fread(cmd = paste("zcat", shQuote(p_expr)))
cat(sprintf("[Load] P13 clone=%d, P15 clone=%d, region=%d, subcluster=%d, expr=%d\n",
            nrow(c13), nrow(c15), nrow(reg), nrow(sub), nrow(expr)))

## 2) 统一 key：<slice>_<barcode>（region 用 sample 派生 slice；其余用 section 后缀映射）
suf2slice <- c("1"="092A","2"="092B","3"="093A","4"="093B","5"="093C","6"="093D",
               "7"="094A","8"="094B","9"="094C","10"="094D","11"="095A","12"="095B")
reg[, slice := sub("^GSM[0-9]+_", "", sample)]
reg[, cell  := paste0(slice, "_", sub("_[0-9]+$", "", barcode))]
sub[, slice := suf2slice[sub("^.*_", "", barcode)]][, cell := paste0(slice, "_", sub("_[0-9]+$", "", barcode))]
expr[, slice := suf2slice[sub("^.*_", "", barcode)]][, cell := paste0(slice, "_", sub("_[0-9]+$", "", barcode))]

## 3) 入口断言：key 唯一、三张 RCTD 表 spot 集合完全一致
stopifnot(!anyDuplicated(reg$cell), !anyDuplicated(sub$cell), !anyDuplicated(expr$cell),
          !anyDuplicated(c13$cell), !anyDuplicated(c15$cell))
stopifnot(setequal(reg$cell, sub$cell), setequal(reg$cell, expr$cell),
          nrow(reg) == nrow(sub), nrow(sub) == nrow(expr))
cat("[Check] key 唯一 + 三张 RCTD 表 spot 集合完全一致 (n =", nrow(reg), ")\n")

## 4) 逐切片关联率（四输入两两）
clone_all <- rbind(c13[, .(cell, clone, clone_int, tree_node, cnv_burden, slice, clone_x = x, clone_y = y)],
                   c15[, .(cell, clone, clone_int, tree_node, cnv_burden, slice, clone_x = x, clone_y = y)])
sl_use <- c("094A","094B","095A","095B")
assoc <- rbindlist(lapply(sort(unique(clone_all$slice)), function(sl) {
  k <- clone_all$cell[clone_all$slice == sl]
  rk <- reg$cell[reg$slice == sl]; sk <- sub$cell[sub$slice == sl]; ek <- expr$cell[expr$slice == sl]
  data.table(slice = sl, n_clone = length(k), n_rctd = length(rk),
             clone_in_region = mean(k %in% rk), clone_in_sub = mean(k %in% sk), clone_in_expr = mean(k %in% ek),
             rctd_in_clone = mean(rk %in% k),
             region_in_sub = mean(rk %in% sk), region_in_expr = mean(rk %in% ek), sub_in_expr = mean(sk %in% ek))
}))
cat("\n--- 逐切片四输入关联率 ---\n"); print(assoc)
stopifnot(all(assoc$clone_in_region == 1), all(assoc$clone_in_sub == 1), all(assoc$clone_in_expr == 1),
          all(assoc$region_in_sub == 1), all(assoc$region_in_expr == 1), all(assoc$sub_in_expr == 1))

## 5) 底表 universe：病人两切片内的全部 RCTD spot（克隆 spot 是其中子集，非克隆 spot 保留 clone=NA）
universe <- reg[slice %in% sl_use, .(cell, slice, barcode, rctd_region = region,
                                     tumor_prop, prop_sum, dist_core, x, y)]
universe[, patient := fifelse(slice %in% c("094A","094B"), "P13", "P15")]
universe <- merge(universe, clone_all[, .(cell, clone, clone_int, tree_node, cnv_burden, clone_x, clone_y)],
                  by = "cell", all.x = TRUE, sort = FALSE)
universe <- merge(universe, sub[, .(cell, c1q_tam_prop = `C1Q TAMs`)], by = "cell", all.x = TRUE, sort = FALSE)
universe <- merge(universe, expr[, c("cell", "COL1A1","FN1","COL3A1","HLA-A","HLA-B","HLA-C",
                                     "HLA-DRA","HLA-DRB1","SPP1","CD44","CD8A","CD4","ITGAV","ITGB1"), with = FALSE],
                  by = "cell", all.x = TRUE, sort = FALSE)
universe[, region := fifelse(rctd_region == "Stroma", NA_character_, rctd_region)]   # 分区只用 Tumor Core / Interface
universe[, c1q_tam_pos := !is.na(c1q_tam_prop) & c1q_tam_prop >= 0.02]

cat("\n--- 底表完整性 ---\n")
cat("universe rows:", nrow(universe), "| dup cell:", anyDuplicated(universe$cell),
    "| clone 缺失(非克隆 spot):", sum(is.na(universe$clone)),
    "| c1q 缺失:", sum(is.na(universe$c1q_tam_prop)), "\n")
print(universe[, .(n = .N, clone_spots = sum(!is.na(clone)), c1q_pos = sum(c1q_tam_pos)),
               by = .(patient, slice, rctd_region)][order(patient, slice, rctd_region)])
cat("克隆 spot 全部保留: P13 =", sum(!is.na(universe$clone) & universe$patient == "P13"),
    "/", nrow(c13), "| P15 =", sum(!is.na(universe$clone) & universe$patient == "P15"), "/", nrow(c15), "\n")

cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 生成可估性格子表（病人 × 区域 × 克隆），并落盘 qs 底表与 estimability.csv
suppressPackageStartupMessages(library(data.table))
out_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step03_clone_LR"

## 1) 空间坐标一致性（决定是否保留克隆表自带 x/y）
xydiff <- universe[!is.na(clone), max(abs(x - clone_x), na.rm = TRUE)]
cat("[Coord] max|x - clone_x| =", xydiff, "| max|y - clone_y| =",
    universe[!is.na(clone), max(abs(y - clone_y), na.rm = TRUE)], "\n")
if (xydiff < 1e-6) universe[, c("clone_x", "clone_y") := NULL]

## 2) 可估性格子：病人 x 区域(Tumor Core/Interface/Stroma) x 克隆
regions_all <- c("Tumor Core", "Interface", "Stroma")
clones_by_patient <- universe[!is.na(clone), .(clones = list(sort(unique(clone)))), by = patient]
est <- rbindlist(lapply(seq_len(nrow(clones_by_patient)), function(i) {
  pt <- clones_by_patient$patient[i]
  rbindlist(lapply(regions_all, function(rg) {
    d <- universe[patient == pt & rctd_region == rg]
    n_rg <- nrow(d); n_recv <- sum(d$c1q_tam_pos)
    rbindlist(lapply(clones_by_patient$clones[[i]], function(cl) {
      n_send <- sum(!is.na(d$clone) & d$clone == cl)
      region_ok <- rg %in% c("Tumor Core", "Interface")
      fail <- c(if (!region_ok) "region_N/A_stroma_excluded",
                if (n_rg < 100) "n_region_spots<100",
                if (n_send < 20) "n_sender_pos<20",
                if (n_recv < 20) "n_receiver_pos<20")
      data.table(patient = pt, region = rg, clone = cl,
                 n_region_spots = n_rg, n_sender_pos = n_send, n_receiver_pos = n_recv,
                 estimable = length(fail) == 0,
                 reason = if (length(fail) == 0) "ok" else paste(fail, collapse = ";"))
    }))
  }))
}))
setcolorder(est, c("patient","region","clone","n_region_spots","n_sender_pos","n_receiver_pos","estimable","reason"))
setorder(est, patient, region, clone)

cat("\n===== estimability.csv =====")
print(est)
cat("\n格子数 =", nrow(est), "| 可估 =", sum(est$estimable), "| 不可估 =", sum(!est$estimable), "\n")
cat("每个病人覆盖的区域:", paste(regions_all, collapse = "/"), "| 每个病人克隆数:",
    paste(clones_by_patient$patient, sapply(clones_by_patient$clones, length), sep = "=", collapse = ", "), "\n")

## 3) 落盘
base_p13 <- universe[patient == "P13"]; base_p15 <- universe[patient == "P15"]
meta <- list(
  built_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
  inputs = list(P13_clone = p13_path, P15_clone = p15_path, rctd_region = p_reg,
                rctd_subcluster = p_sub, lr_expr = p_expr),
  input_sizes_bytes = file.info(c(p13_path, p15_path, p_reg, p_sub, p_expr))$size,
  params = list(c1q_tam_positive_threshold = 0.02, regions_used = c("Tumor Core","Interface"),
                stroma_handling = "N/A (excluded from region-based analysis, spots retained with region=NA)",
                min_region_spots = 100, min_sender_pos = 20, min_receiver_pos = 20,
                subsampling = FALSE, unit = "patient x region x clone",
                patient_slice_map = list(P13 = c("094A","094B"), P15 = c("095A","095B"))),
  key_rule = "<slice>_<barcode>; RCTD tables keyed by section suffix 1..12 -> 092A..095B",
  universe = "all RCTD spots of the patient's two slices; clone=NA for non-clone spots",
  n_spots = list(P13 = nrow(base_p13), P15 = nrow(base_p15)),
  n_clone_spots = list(P13 = sum(!is.na(base_p13$clone)), P15 = sum(!is.na(base_p15$clone))),
  association_rates_by_slice = assoc,
  reason_legend = c(ok = "all criteria met",
                    "region_N/A_stroma_excluded" = "Stroma not used per analysis spec",
                    "n_region_spots<100" = "fewer than 100 spots in region",
                    "n_sender_pos<20" = "fewer than 20 sender-positive spots for this clone",
                    "n_receiver_pos<20" = "fewer than 20 C1Q-TAM-positive spots in region")
)
save_obj <- function(obj, filepath) {
  if (requireNamespace("qs", quietly = TRUE)) { qs::qsave(obj, filepath, preset = "fast")
    cat(sprintf("[Save] qs: %s (%.1f MB)\n", filepath, file.info(filepath)$size/1e6))
  } else { saveRDS(obj, sub("\\.qs$", ".rds", filepath)) }
}
qs_path <- file.path(out_dir, "data", "checkpoint", "clone_lr_input.qs")
save_obj(list(P13 = base_p13, P15 = base_p15, meta = meta), qs_path)
est_path <- file.path(out_dir, "tables", "estimability.csv")
fwrite(est, est_path)
cat("[Save] csv:", est_path, "|", file.info(est_path)$size, "bytes | rows =", nrow(est), "\n")

## 4) 落盘后回读校验
back <- qs::qread(qs_path)
cat("[Verify] qs 回读: P13 =", nrow(back$P13), "| P15 =", nrow(back$P15),
    "| 列数 =", ncol(back$P13), "/", ncol(back$P15),
    "| P13 切片 =", paste(sort(unique(back$P13$slice)), collapse = ","),
    "| P15 切片 =", paste(sort(unique(back$P15$slice)), collapse = ","), "\n")
cat("[Verify] estimability 回读行数 =", nrow(fread(est_path)), "\n")
cat("底表列名:\n"); print(names(back$P13))
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 独立复核：换一条代码路径重算关键格子；检查表达量尺度/缺失；补 receiver 口径说明后重存
suppressPackageStartupMessages(library(data.table))
out_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step03_clone_LR"
qs_path <- file.path(out_dir, "data", "checkpoint", "clone_lr_input.qs")
back <- qs::qread(qs_path)

## 1) 表达量尺度与缺失（G03/G07 数据尺度审计）
genes <- c("COL1A1","FN1","COL3A1","HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","SPP1","CD44","CD8A","CD4","ITGAV","ITGB1")
gstats <- rbindlist(lapply(c("P13","P15"), function(pt) {
  d <- back[[pt]]
  rbindlist(lapply(genes, function(g) data.table(patient = pt, gene = g, min = min(d[[g]]), max = max(d[[g]]),
                                                  mean = mean(d[[g]]), n_na = sum(is.na(d[[g]])))))
}))
cat("\n--- 表达量尺度审计（14 基因 x 2 病人）---\n")
print(rbind(gstats[, .(min = min(min), max = max(max), mean = mean(mean), total_na = sum(n_na)), by = patient]))
cat("负值数:", sum(gstats$min < 0), "| log1p 量级(<20) 的基因数:", sum(gstats$max < 20), "/", nrow(gstats), "\n")

## 2) 克隆字段完整性
cat("\n--- 克隆字段完整性 ---\n")
print(rbind(back$P13[!is.na(clone), .(patient="P13", n=.N, na_clone_int=sum(is.na(clone_int)),
          na_tree=sum(is.na(tree_node)), na_burden=sum(is.na(cnv_burden)))],
            back$P15[!is.na(clone), .(patient="P15", n=.N, na_clone_int=sum(is.na(clone_int)),
          na_tree=sum(is.na(tree_node)), na_burden=sum(is.na(cnv_burden)))]))

## 3) 独立复算（直接用原始三表，不走 universe 构建路径）
c13 <- fread("/path/to/project/project_spatial/TNBC/patient_ALL/clone/P13/clone_membership.csv")
reg <- fread("/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment_extended.csv")
sub <- fread(cmd = "zcat /path/to/project/project_singlecell/TNBC_result/figure18/intermediates/subcluster_prop.csv.gz")
reg[, slice := sub("^GSM[0-9]+_", "", sample)][, cell := paste0(slice, "_", sub("_[0-9]+$", "", barcode))]
sub[, cell := paste0(c("1"="092A","2"="092B","3"="093A","4"="093B","5"="093C","6"="093D","7"="094A",
                       "8"="094B","9"="094C","10"="094D","11"="095A","12"="095B")[sub("^.*_", "", barcode)],
                     "_", sub("_[0-9]+$", "", barcode))]
chk <- merge(reg[slice %in% c("094A","094B"), .(cell, region)], sub[, .(cell, c1q = `C1Q TAMs`)], by = "cell")
chk <- merge(chk, c13[, .(cell, clone)], by = "cell", all.x = TRUE)
cat("\n--- 独立复算 P13 Interface（不走底表构建路径）---\n")
cat("n_region_spots =", chk[region=="Interface", .N],
    "| Clone1 sender =", chk[region=="Interface" & !is.na(clone) & clone=="Clone1", .N],
    "| Clone2 sender =", chk[region=="Interface" & !is.na(clone) & clone=="Clone2", .N],
    "| receiver =", chk[region=="Interface" & c1q >= 0.02, .N], "\n")
cat("--- 独立复算 P15 Tumor Core ---\n")
reg15 <- reg[slice %in% c("095A","095B")]
chk15 <- merge(reg15[, .(cell, region)], sub[, .(cell, c1q = `C1Q TAMs`)], by = "cell")
c15 <- fread("/path/to/project/project_spatial/TNBC/patient_ALL/clone/P15/clone_membership.csv")
chk15 <- merge(chk15, c15[, .(cell, clone)], by = "cell", all.x = TRUE)
cat("n_region_spots =", chk15[region=="Tumor Core", .N],
    "| Clone1 sender =", chk15[region=="Tumor Core" & !is.na(clone) & clone=="Clone1", .N],
    "| Clone2 sender =", chk15[region=="Tumor Core" & !is.na(clone) & clone=="Clone2", .N],
    "| receiver =", chk15[region=="Tumor Core" & c1q >= 0.02, .N], "\n")

## 4) 敏感性：若改为「仅在克隆 spot 内」计算的对照（不改变主交付口径）
sens <- rbindlist(lapply(c("P13","P15"), function(pt) {
  d <- back[[pt]][!is.na(clone)]
  rbindlist(lapply(c("Tumor Core","Interface"), function(rg) {
    rbindlist(lapply(sort(unique(d$clone)), function(cl) {
      dd <- d[rctd_region == rg]
      data.table(patient = pt, region = rg, clone = cl,
                 n_clone_region_spots = nrow(dd), n_sender = sum(dd$clone == cl),
                 n_receiver_within_clone = sum(dd$clone == cl & dd$c1q_tam_pos))
    }))
  }))
}))
cat("\n--- 敏感性对照：克隆 spot 内口径下的发送/受端计数（仅供裁决，不改变主交付）---\n")
print(sens)

## 5) 补记 receiver 口径并重存 qs
back$meta$params$receiver_scope <- "region-level: C1Q TAM positive spots in the region (independent of tumour clone)"
back$meta$params$sender_scope <- "clone-specific: spots in the region assigned to that clone"
back$meta$sensitivity_clone_anchored <- sens
back$meta$association_note <- paste("clone->{region,subcluster,expr} = 100% in every slice;",
  "region/subcluster/expr are the same 16607-spot set;",
  "rctd_in_clone < 100% reflects that the clone table lists only CNV-assigned tumour-clone spots,",
  "so non-clone spots are retained with clone=NA rather than dropped (1072 of 1105 non-clone spots are Interface)")
qs::qsave(back, qs_path, preset = "fast")
cat("\n[Save] qs 更新:", qs_path, "|", file.info(qs_path)$size, "bytes\n")
cat("[Verify] meta.params$receiver_scope =", qs::qread(qs_path)$meta$params$receiver_scope, "\n")
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 收尾自检：notebook 错误 cell 扫描、交付文件存在性、会话标识（供审计日志）
suppressPackageStartupMessages(library(data.table))
nb <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/notebooks/step03_clone_LR.ipynb"
out_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step03_clone_LR"
cat("[Notebook] exists =", file.exists(nb), "| size =", file.info(nb)$size, "bytes\n")
j <- jsonlite::fromJSON(nb, simplifyVector = FALSE)
n_code <- 0; err_cells <- integer(0)
for (i in seq_along(j$cells)) {
  cl <- j$cells[[i]]
  if (identical(cl$cell_type, "code")) {
    n_code <- n_code + 1
    if (any(vapply(cl$outputs, function(o) identical(o$output_type, "error"), logical(1)))) err_cells <- c(err_cells, i)
  }
}
cat("[Notebook] 总 cell =", length(j$cells), "| code cell =", n_code,
    "| 含 error 的 cell 索引 =", if (length(err_cells)) paste(err_cells, collapse = ",") else "none",
    "| status =", if (length(err_cells)) "needs_fix" else "clean", "\n")

cat("\n--- 交付文件 ---\n")
for (f in c(file.path(out_dir, "data", "checkpoint", "clone_lr_input.qs"),
            file.path(out_dir, "tables", "estimability.csv"))) {
  cat(sprintf("%s | exists=%s | size=%s bytes\n", f, file.exists(f), file.info(f)$size))
}
cat("\n--- estimability.csv 原样回读 ---\n")
print(fread(file.path(out_dir, "tables", "estimability.csv")))

cat("\n--- 会话/环境标识 ---\n")
cat("R:", R.version.string, "| qs:", as.character(packageVersion("qs")),
    "| data.table:", as.character(packageVersion("data.table")), "| jsonlite:", as.character(packageVersion("jsonlite")), "\n")
envs <- Sys.getenv()
idx <- grep("SESSION|CLAUDE|JUPYTER", names(envs), value = TRUE)
print(envs[idx[!idx %in% c("JUPYTER_TOKEN")]])
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 将 P13 / P15 两张病人底表从 checkpoint 导出为 csv（显式 na="NA"，保证下游读回为真正的缺失值）
suppressPackageStartupMessages(library(data.table))
out_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step03_clone_LR"
qs_path <- file.path(out_dir, "data", "checkpoint", "clone_lr_input.qs")
x <- qs::qread(qs_path)
p13_csv <- file.path(out_dir, "tables", "clone_lr_input_P13.csv")
p15_csv <- file.path(out_dir, "tables", "clone_lr_input_P15.csv")
fwrite(x$P13, p13_csv, na = "NA"); fwrite(x$P15, p15_csv, na = "NA")
cat(sprintf("[Save] %s | %d rows x %d cols | %.1f MB\n", p13_csv, nrow(x$P13), ncol(x$P13), file.info(p13_csv)$size/1e6))
cat(sprintf("[Save] %s | %d rows x %d cols | %.1f MB\n", p15_csv, nrow(x$P15), ncol(x$P15), file.info(p15_csv)$size/1e6))

# 回读校验：行数/列数/切片/克隆与受端计数与 qs 完全一致
b13 <- fread(p13_csv); b15 <- fread(p15_csv)
stopifnot(nrow(b13) == nrow(x$P13), nrow(b15) == nrow(x$P15),
          setequal(unique(b13$slice), c("094A","094B")), setequal(unique(b15$slice), c("095A","095B")),
          sum(!is.na(b13$clone)) == sum(!is.na(x$P13$clone)),
          sum(!is.na(b15$clone)) == sum(!is.na(x$P15$clone)),
          sum(b13$c1q_tam_pos) == sum(x$P13$c1q_tam_pos), sum(b15$c1q_tam_pos) == sum(x$P15$c1q_tam_pos),
          identical(names(b13), names(x$P13)), identical(names(b15), names(x$P15)))
cat("[Verify] csv 与 qs 一致 | P13: rows =", nrow(b13), "clone spots =", sum(!is.na(b13$clone)),
    "c1q_pos =", sum(b13$c1q_tam_pos), "\n")
cat("[Verify] P15: rows =", nrow(b15), "clone spots =", sum(!is.na(b15$clone)),
    "c1q_pos =", sum(b15$c1q_tam_pos), "\n")
cat("[Verify] 缺失值读回: P13 clone NA =", sum(is.na(b13$clone)), "| 空串 =", sum(b13$clone == "", na.rm = TRUE), "\n")
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 重新扫描 notebook 错误 cell（cell 11 已就地修复）并列出全部交付文件
nb <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/notebooks/step03_clone_LR.ipynb"
j <- jsonlite::fromJSON(nb, simplifyVector = FALSE)
n_code <- 0; err_cells <- integer(0)
for (i in seq_along(j$cells)) {
  cl <- j$cells[[i]]
  if (identical(cl$cell_type, "code")) {
    n_code <- n_code + 1
    if (any(vapply(cl$outputs, function(o) identical(o$output_type, "error"), logical(1)))) err_cells <- c(err_cells, i)
  }
}
cat("[Notebook] 总 cell =", length(j$cells), "| code cell =", n_code,
    "| 含 error 的 cell 索引 =", if (length(err_cells)) paste(err_cells, collapse = ",") else "none",
    "| status =", if (length(err_cells)) "needs_fix" else "clean", "\n")
cat("\n--- 输出目录文件清单 ---\n")
out_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step03_clone_LR"
fs <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
for (f in fs) cat(sprintf("%-70s %10.1f KB\n", sub(paste0(out_dir, "/"), "", f), file.info(f)$size/1024))
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 收尾对账（只读）：格子表与底表的分母/分子是否自洽
suppressPackageStartupMessages(library(data.table))
out_dir <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step03_clone_LR"
x <- qs::qread(file.path(out_dir, "data", "checkpoint", "clone_lr_input.qs"))
est <- fread(file.path(out_dir, "tables", "estimability.csv"))

for (pt in c("P13","P15")) {
  d <- x[[pt]]
  cat("\n[", pt, "] 底表 row =", nrow(d), "| clone spots =", sum(!is.na(d$clone)),
      "| c1q_pos =", sum(d$c1q_tam_pos), "\n")
  cat("  区域点数（rctd_region 原始口径）:\n")
  print(d[, .(spots = .N, clone_spots = sum(!is.na(clone)), c1q_pos = sum(c1q_tam_pos)),
          by = rctd_region][order(rctd_region)])
  e <- est[patient == pt]
  cat("  格子表按区域求和 vs 底表:\n")
  print(e[, .(sum_n_region_spots = sum(n_region_spots), sum_n_sender = sum(n_sender_pos),
              n_receiver_pos = unique(n_receiver_pos)), by = region][order(region)])
}
# 断言：每区域 n_region_spots 在克隆间重复，求和会重复计数，故用 unique 校验
for (pt in c("P13","P15")) {
  d <- x[[pt]]; e <- est[patient == pt]
  for (rg in unique(d$rctd_region)) {
    stopifnot(unique(e$n_region_spots[e$region == rg]) == nrow(d[rctd_region == rg]))
    stopifnot(unique(e$n_receiver_pos[e$region == rg]) == sum(d[rctd_region == rg]$c1q_tam_pos))
    stopifnot(sum(e$n_sender_pos[e$region == rg]) == sum(!is.na(d[rctd_region == rg]$clone)))
  }
}
cat("\n[Verify] 全部区域的分母(点数)、受端计数、发送端合计与底表逐一吻合\n")
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 删除空 cell 后重新扫描：错误 cell / 空 cell 计数
nb <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/notebooks/step03_clone_LR.ipynb"
j <- jsonlite::fromJSON(nb, simplifyVector = FALSE)
n_code <- 0; err <- integer(0); empt <- integer(0)
for (i in seq_along(j$cells)) {
  cl <- j$cells[[i]]
  src <- if (is.list(cl$source)) paste(unlist(cl$source), collapse = "") else cl$source
  if (identical(cl$cell_type, "code")) {
    n_code <- n_code + 1
    if (trimws(src) == "") empt <- c(empt, i)
    if (any(vapply(cl$outputs, function(o) identical(o$output_type, "error"), logical(1)))) err <- c(err, i)
  }
}
cat("[Notebook] total =", length(j$cells), "| code =", n_code,
    "| error cells =", if (length(err)) paste(err, collapse = ",") else "none",
    "| empty code cells =", if (length(empt)) paste(empt, collapse = ",") else "none", "\n")
cat("===== Cell End:", date(), "=====\n")

cat("===== Cell Start:", date(), "=====\n")
# 独立 critic 全部 PASS 后：追加审计凭证到统一日志（东八区时间），并打印 kernel 变量清单
tz <- format(Sys.time(), "%z")
logdir <- file.path("/path/to/.claude/logs", format(Sys.Date(), "%Y-%m-%d"))
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
line <- sprintf('{"ts":"%s","session":"sc-analyst:step03_clone_LR","subagent_type":"sc-analyst","verdict":"pass"}',
                format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
cat("[Audit Log] tz =", tz, "| file =", file.path(logdir, "audit_done.jsonl"), "\n")
cat("[Audit Log] line =", line, "\n")
cat(line, "\n", file = file.path(logdir, "audit_done.jsonl"), append = TRUE, sep = "")
cat("[Audit Log] 写入后最后一行 =", tail(readLines(file.path(logdir, "audit_done.jsonl")), 1), "\n")

# kernel_variables 用：登记内存中可复用的数据对象
objs <- list()
for (v in c("universe", "est", "x", "assoc")) {
  if (exists(v)) {
    o <- get(v)
    objs[[v]] <- list(class = paste(class(o), collapse = "/"),
                      dim = if (!is.null(dim(o))) paste(dim(o), collapse = "x") else length(o))
  }
}
cat("\n[Kernel Variables]\n"); print(objs)
cat("===== Cell End:", date(), "=====\n")

# 读取空间表达矩阵（只读参考），按 <切片>_<barcode> 与底表对齐，构建 NicheNet 输入快照
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

nn_root    <- "/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone/step03_clone_LR"
nn_fig_dir <- file.path(nn_root, "figures")
nn_tab_dir <- file.path(nn_root, "tables")
nn_dat_dir <- file.path(nn_root, "data", "checkpoint")
for (dd in c(nn_fig_dir, nn_tab_dir, nn_dat_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

nn_spatial_path <- "/path/to/project/project_singlecell/TNBC_result/final_figures/V2/spatial/spatial_sct_matrix.rds"
nn_lr_path      <- "/path/to/data/other_ans_project/SeekSpace_code/lr_network.rds"
nn_ltm_path     <- "/path/to/data/other_ans_project/SeekSpace_code/ligand_target_matrix.rds"
nn_smap <- c("092A","092B","093A","093B","093C","093D","094A","094B","094C","094D","095A","095B")

# 表达矩阵列名形如 <barcode>_<section 1..12>，按 section 映射回切片名
nn_expr <- readRDS(nn_spatial_path)
nn_sect <- as.integer(sub("^.*_", "", colnames(nn_expr)))
stopifnot(all(nn_sect %in% 1:12))
colnames(nn_expr) <- paste(nn_smap[nn_sect], sub("_[0-9]+$", "", colnames(nn_expr)), sep = "_")
stopifnot(!anyDuplicated(colnames(nn_expr)))
nn_missing <- setdiff(universe$cell, colnames(nn_expr))
stopifnot(length(nn_missing) == 0)
nn_expr <- nn_expr[, universe$cell]
stopifnot(identical(colnames(nn_expr), universe$cell))

nn_lr  <- unique(as.data.frame(readRDS(nn_lr_path))[, c("from", "to")])
nn_ltm <- readRDS(nn_ltm_path)

cat("[Input] 表达矩阵:", nrow(nn_expr), "基因 x", ncol(nn_expr), "spot\n")
cat("[Input] LR 网络配对:", nrow(nn_lr), "| 配体数:", length(unique(nn_lr$from)), "\n")
cat("[Input] ligand-target 矩阵:", nrow(nn_ltm), "基因 x", ncol(nn_ltm), "配体\n")
cat("[Input] 与底表 spot 完全对齐:", identical(colnames(nn_expr), universe$cell), "| 缺失:", length(nn_missing), "\n")
cat("[Input] 表达值范围:", paste(round(range(nn_expr@x), 3), collapse = " ~ "), "(log 归一化)\n")

qs::qsave(list(expr = nn_expr, cells = universe$cell, smap = nn_smap, source = nn_spatial_path),
          file.path(nn_dat_dir, "nichenet_input.qs"), preset = "fast")
cat("[Save]", file.path(nn_dat_dir, "nichenet_input.qs"),
    round(file.info(file.path(nn_dat_dir, "nichenet_input.qs"))$size / 1e6, 1), "MB\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 定义三个基础函数：受端基因集（Spearman）、候选配体筛选、NicheNet 配体活性打分
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

nn_uni <- as.data.frame(universe)
nn_expr_frac <- 0.10   # 表达率阈值：背景基因 / 配体 / 受体统一 10%
nn_n_geneset  <- 200   # 受端基因集固定 200 个，两克隆必须同规模才可比

# 1) 受端程序：受端 spot 内，基因表达 vs 该克隆指示变量的 Spearman 相关，取 ρ 前 200 正相关基因
nn_receiver_program <- function(pat, rg, cl) {
  d  <- nn_uni[!is.na(nn_uni$region) & nn_uni$region == rg & nn_uni$patient == pat, , drop = FALSE]
  rec <- d[which(d$c1q_tam_pos), , drop = FALSE]
  ind <- as.numeric(!is.na(rec$clone) & rec$clone == cl)

  Xr <- nn_expr[, rec$cell, drop = FALSE]
  rate_r <- Matrix::rowMeans(Xr > 0)

  # Spearman = 秩上的 Pearson；matrixStats::rowRanks 逐基因取平均秩
  rk  <- matrixStats::rowRanks(as.matrix(Xr), ties.method = "average")
  ri  <- rank(ind, ties.method = "average")
  rkc <- rk - rowMeans(rk)
  ric <- ri - mean(ri)
  den <- sqrt(rowSums(rkc^2)) * sqrt(sum(ric^2))
  rho <- as.vector(rkc %*% ric) / den
  rho[!is.finite(rho)] <- NA
  names(rho) <- rownames(Xr)

  # 候选基因限定在 ligand-target 先验模型覆盖的基因内（模型外基因无法参与活性预测，不可避免）
  cand <- intersect(names(rho), rownames(nn_ltm))
  cand <- cand[is.finite(rho[cand])]
  geneset <- utils::head(cand[order(-rho[cand])], nn_n_geneset)
  geneset <- geneset[rho[geneset] > 0]

  background <- names(rate_r)[rate_r >= nn_expr_frac]
  list(geneset = geneset, rho = rho, background = background,
       n_sender = sum(!is.na(d$clone) & d$clone == cl),
       n_receiver = nrow(rec), n_receiver_pos = sum(ind),
       n_candidate_genes = length(cand), candidate_min_rho = min(rho[cand]))
}

# 2) 候选配体：发送端表达率 ≥10%，且受端存在表达率 ≥10% 的受体
nn_potential_ligands <- function(sender_cells, receiver_cells) {
  sender_rate <- Matrix::rowMeans(nn_expr[, sender_cells, drop = FALSE] > 0)
  rec_rate    <- Matrix::rowMeans(nn_expr[, receiver_cells, drop = FALSE] > 0)
  expressed_ligands   <- unique(nn_lr$from[nn_lr$from %in% names(sender_rate)[sender_rate >= nn_expr_frac]])
  expressed_receptors <- names(rec_rate)[rec_rate >= nn_expr_frac]
  lr_f <- nn_lr[nn_lr$from %in% expressed_ligands & nn_lr$to %in% expressed_receptors, , drop = FALSE]
  list(ligands = unique(lr_f$from), sender_rate = sender_rate, n_lr_pairs = nrow(lr_f))
}

# 3) 配体活性：pearson 活性值 + aupr/auroc，rank 按 pearson 降序
nn_activity <- function(geneset, background, potential_ligands) {
  act <- as.data.frame(nichenetr::predict_ligand_activities(
    geneset = geneset, background_expressed_genes = background,
    ligand_target_matrix = nn_ltm, potential_ligands = potential_ligands))
  act <- act[order(-act$pearson), , drop = FALSE]
  act$rank <- seq_len(nrow(act))
  rownames(act) <- NULL
  act
}

cat("[Function] nn_receiver_program / nn_potential_ligands / nn_activity 已定义\n")
cat("[Params] 表达率阈值 =", nn_expr_frac, "| 基因集规模 =", nn_n_geneset, "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 逐格运行 8 个可估格（病人 x 区域 x 克隆）：计算受端基因集、候选配体与配体活性
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

nn_cell_grid <- as.data.frame(est)
nn_cell_grid <- nn_cell_grid[nn_cell_grid$estimable == TRUE, c("patient", "region", "clone")]
stopifnot(nrow(nn_cell_grid) == 8, !any(nn_cell_grid$region == "Stroma"))

nn_prog <- list(); nn_lig <- list(); nn_act <- list()
nn_cell_stats <- list()

for (i in seq_len(nrow(nn_cell_grid))) {
  pat <- nn_cell_grid$patient[i]; rg <- nn_cell_grid$region[i]; cl <- nn_cell_grid$clone[i]
  key <- paste(pat, rg, cl, sep = "|")
  t1 <- Sys.time()

  d  <- nn_uni[!is.na(nn_uni$region) & nn_uni$region == rg & nn_uni$patient == pat, , drop = FALSE]
  sender_cells <- d$cell[!is.na(d$clone) & d$clone == cl]
  rec_cells    <- d$cell[which(d$c1q_tam_pos)]

  prog <- nn_receiver_program(pat, rg, cl)
  lig  <- nn_potential_ligands(sender_cells, rec_cells)
  act  <- nn_activity(prog$geneset, prog$background, lig$ligands)
  act$sender_expr_rate <- as.numeric(lig$sender_rate[act$test_ligand])

  nn_prog[[key]] <- prog
  nn_lig[[key]]  <- lig
  nn_act[[key]]  <- act
  nn_cell_stats[[key]] <- data.frame(
    patient = pat, region = rg, clone = cl,
    n_region_spots = nrow(d), n_sender = length(sender_cells),
    n_receiver = length(rec_cells), n_receiver_pos = prog$n_receiver_pos,
    n_candidate_genes = prog$n_candidate_genes,
    n_geneset = length(prog$geneset), geneset_min_rho = min(prog$rho[prog$geneset]),
    n_background = length(prog$background),
    n_geneset_in_background = sum(prog$geneset %in% prog$background),
    n_potential_ligands = length(lig$ligands),
    top_ligand = act$test_ligand[1], top_ligand_pearson = round(act$pearson[1], 4),
    top_ligand_sender_rate = round(act$sender_expr_rate[1], 3),
    runtime_s = round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1))
  cat(sprintf("[Cell] %s %s %s | sender=%d recv=%d recv_pos=%d geneset=%d bg=%d ligands=%d | top: %s (r=%.3f)\n",
              pat, rg, cl, length(sender_cells), length(rec_cells), prog$n_receiver_pos,
              length(prog$geneset), length(prog$background), length(lig$ligands),
              act$test_ligand[1], act$pearson[1]))
}

nn_cell_stats <- do.call(rbind, nn_cell_stats)
rownames(nn_cell_stats) <- NULL
cat("\n[Summary] 可估格数:", nrow(nn_cell_stats), "| 全部基因集规模:", paste(unique(nn_cell_stats$n_geneset), collapse = ","), "\n")
print(nn_cell_stats[, c("patient","region","clone","n_sender","n_receiver","n_receiver_pos",
                        "n_geneset","n_background","n_potential_ligands","top_ligand")])
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 敏感性分析：同一病人同一区域内，用同一受端基因集与同一配体集合比较两克隆（排除基因集差异带来的不可比）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

nn_regions <- unique(nn_cell_grid[, c("patient", "region")])
rownames(nn_regions) <- NULL
nn_shared <- list(); nn_geneset_overlap <- list()

for (i in seq_len(nrow(nn_regions))) {
  pat <- nn_regions$patient[i]; rg <- nn_regions$region[i]
  k1 <- paste(pat, rg, "Clone1", sep = "|"); k2 <- paste(pat, rg, "Clone2", sep = "|")

  union_lig <- sort(union(nn_lig[[k1]]$ligands, nn_lig[[k2]]$ligands))
  gs1 <- nn_prog[[k1]]$geneset; gs2 <- nn_prog[[k2]]$geneset
  bg  <- nn_prog[[k1]]$background   # 同区域两克隆受端集合相同，背景一致

  nn_geneset_overlap[[paste(pat, rg, sep = "|")]] <- data.frame(
    patient = pat, region = rg,
    n_geneset_clone1 = length(gs1), n_geneset_clone2 = length(gs2),
    n_shared = length(intersect(gs1, gs2)),
    jaccard = round(length(intersect(gs1, gs2)) / length(union(gs1, gs2)), 4))

  for (src in c("Clone1", "Clone2")) {
    gs <- if (src == "Clone1") gs1 else gs2
    act <- nn_activity(gs, bg, union_lig)
    act$geneset_source <- src
    nn_shared[[paste(pat, rg, src, sep = "|")]] <- act
    cat(sprintf("[Shared] %s %s | gene set=%s | ligands=%d | top: %s (r=%.3f)\n",
                pat, rg, src, nrow(act), act$test_ligand[1], act$pearson[1]))
  }
}

nn_geneset_overlap <- do.call(rbind, nn_geneset_overlap)
rownames(nn_geneset_overlap) <- NULL
cat("\n[Gene set overlap] 同区域内两克隆 200 基因集的重叠:\n")
print(nn_geneset_overlap)
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 汇总对比表：主表（各克隆自身受端基因集下的活性差与排序变化）+ 共享受端口径的敏感性表
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

# ---- 1. 长表：8 个可估格的配体活性 ----
nn_act_long <- do.call(rbind, lapply(names(nn_act), function(k) {
  a <- nn_act[[k]]; p <- strsplit(k, "\\|")[[1]]
  data.frame(patient = p[1], region = p[2], clone = p[3], ligand = a$test_ligand,
             pearson = round(a$pearson, 5), aupr_corrected = round(a$aupr_corrected, 5),
             auroc = round(a$auroc, 5), rank = a$rank,
             sender_expr_rate = round(a$sender_expr_rate, 4), stringsAsFactors = FALSE)
}))
rownames(nn_act_long) <- NULL

# ---- 2. 受端基因集明细 ----
nn_gs_tab <- do.call(rbind, lapply(names(nn_prog), function(k) {
  p <- strsplit(k, "\\|")[[1]]; gs <- nn_prog[[k]]$geneset
  data.frame(patient = p[1], region = p[2], clone = p[3], gene = gs,
             rho = round(as.numeric(nn_prog[[k]]$rho[gs]), 4),
             rank_in_geneset = seq_along(gs), stringsAsFactors = FALSE)
}))
rownames(nn_gs_tab) <- NULL

# ---- 3. 主表：同病人同区域内 Clone1 vs Clone2（各自自身基因集/配体集合）----
nn_main <- list()
for (i in seq_len(nrow(nn_regions))) {
  pat <- nn_regions$patient[i]; rg <- nn_regions$region[i]
  k1 <- paste(pat, rg, "Clone1", sep = "|"); k2 <- paste(pat, rg, "Clone2", sep = "|")
  a1 <- nn_act[[k1]]; a2 <- nn_act[[k2]]
  ligs <- sort(union(a1$test_ligand, a2$test_ligand))
  out <- data.frame(patient = pat, region = rg, ligand = ligs, stringsAsFactors = FALSE)
  out$in_clone1_sender  <- ligs %in% a1$test_ligand
  out$in_clone2_sender  <- ligs %in% a2$test_ligand
  out$clone1_pearson <- a1$pearson[match(ligs, a1$test_ligand)]
  out$clone2_pearson <- a2$pearson[match(ligs, a2$test_ligand)]
  out$clone1_rank    <- a1$rank[match(ligs, a1$test_ligand)]
  out$clone2_rank    <- a2$rank[match(ligs, a2$test_ligand)]
  out$delta_pearson  <- round(out$clone1_pearson - out$clone2_pearson, 5)
  out$rank_shift_C1_minus_C2 <- out$clone2_rank - out$clone1_rank
  out$n_potential_clone1 <- nrow(a1); out$n_potential_clone2 <- nrow(a2)
  nn_main[[paste(pat, rg, sep = "|")]] <- out
}
nn_main <- do.call(rbind, nn_main); rownames(nn_main) <- NULL

# ---- 4. 共享受端口径（敏感性）：同一基因集下两克隆配体成员资格的对比 ----
nn_shared_long <- do.call(rbind, lapply(names(nn_shared), function(k) {
  p <- strsplit(k, "\\|")[[1]]; a <- nn_shared[[k]]
  k1 <- paste(p[1], p[2], "Clone1", sep = "|"); k2 <- paste(p[1], p[2], "Clone2", sep = "|")
  data.frame(patient = p[1], region = p[2], geneset_source = p[3], ligand = a$test_ligand,
             pearson = round(a$pearson, 5), rank = a$rank,
             in_clone1_sender = a$test_ligand %in% nn_act[[k1]]$test_ligand,
             in_clone2_sender = a$test_ligand %in% nn_act[[k2]]$test_ligand,
             stringsAsFactors = FALSE)
}))
rownames(nn_shared_long) <- NULL

nn_shared_summary <- list()
for (pk in unique(nn_shared_long$patient)) for (rg in unique(nn_shared_long$region)) {
  for (src in c("Clone1", "Clone2")) {
    s <- nn_shared_long[nn_shared_long$patient == pk & nn_shared_long$region == rg & nn_shared_long$geneset_source == src, ]
    if (nrow(s) == 0) next
    c1only <- s[s$in_clone1_sender & !s$in_clone2_sender, ]
    c2only <- s[s$in_clone2_sender & !s$in_clone1_sender, ]
    both   <- s[s$in_clone1_sender & s$in_clone2_sender, ]
    top <- function(x) if (nrow(x) == 0) c(NA_character_, NA_real_) else c(x$ligand[which.max(x$pearson)], max(x$pearson))
    t1 <- top(c1only); t2 <- top(c2only); tb <- top(both)
    nn_shared_summary[[paste(pk, rg, src, sep = "|")]] <- data.frame(
      patient = pk, region = rg, geneset_source = src,
      n_clone1_only = nrow(c1only), n_clone2_only = nrow(c2only), n_shared = nrow(both),
      top_clone1_only_ligand = t1[1], top_clone1_only_pearson = as.numeric(t1[2]),
      top_clone2_only_ligand = t2[1], top_clone2_only_pearson = as.numeric(t2[2]),
      top_shared_ligand = tb[1], top_shared_pearson = as.numeric(tb[2]), stringsAsFactors = FALSE)
  }
}
nn_shared_summary <- do.call(rbind, nn_shared_summary); rownames(nn_shared_summary) <- NULL

# ---- 5. 落盘 ----
fwrite_safe <- function(x, f) data.table::fwrite(x, file.path(nn_tab_dir, f), na = "NA")
fwrite_safe(nn_cell_stats,      "nichenet_cell_summary.csv")
fwrite_safe(nn_act_long,        "nichenet_ligand_activity_long.csv")
fwrite_safe(nn_gs_tab,          "nichenet_receiver_genesets.csv")
fwrite_safe(nn_main,            "nichenet_clone1_vs_clone2.csv")
fwrite_safe(nn_shared_long,     "nichenet_shared_readout.csv")
fwrite_safe(nn_shared_summary,  "nichenet_shared_readout_summary.csv")
fwrite_safe(nn_geneset_overlap, "nichenet_geneset_overlap.csv")

cat("[Save] tables:\n")
for (f in c("nichenet_cell_summary.csv","nichenet_ligand_activity_long.csv","nichenet_receiver_genesets.csv",
            "nichenet_clone1_vs_clone2.csv","nichenet_shared_readout.csv","nichenet_shared_readout_summary.csv",
            "nichenet_geneset_overlap.csv")) {
  cat(sprintf("  %-40s %8.1f KB  rows=%d\n", f,
              file.info(file.path(nn_tab_dir, f))$size/1e3,
              nrow(data.table::fread(file.path(nn_tab_dir, f)))))
}

cat("\n[Main] 每区域内活性最高的 3 个配体（各自基因集）:\n")
for (pk in unique(nn_main$patient)) for (rg in unique(nn_main$region)) {
  s <- nn_main[nn_main$patient == pk & nn_main$region == rg, ]
  a1 <- s[order(-s$clone1_pearson), ][1:3, ]; a2 <- s[order(-s$clone2_pearson), ][1:3, ]
  cat(sprintf("  %s %s | Clone1: %s | Clone2: %s\n", pk, rg,
              paste0(a1$ligand, "(", round(a1$clone1_pearson,3), ")", collapse=", "),
              paste0(a2$ligand, "(", round(a2$clone2_pearson,3), ")", collapse=", ")))
}
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 图1：配体活性热图（病人 x 区域 x 克隆），行为各格 top 配体的并集，颜色为 NicheNet pearson 活性
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
suppressPackageStartupMessages({ library(ggplot2); library(grid) })

nn_theme <- "nature"
nn_pal_qual <- if (requireNamespace("ggsci", quietly = TRUE)) ggsci::pal_npg()(8) else
  c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000")

# 每格取 top 5 配体的并集（按最高活性排序，最多 30 行）
nn_top5 <- nn_act_long[nn_act_long$rank <= 5, ]
nn_top_ligs <- names(sort(tapply(nn_top5$pearson, nn_top5$ligand, max), decreasing = TRUE))[1:30]

nn_hm <- nn_act_long[nn_act_long$ligand %in% nn_top_ligs, ]
nn_hm$cell_label <- factor(paste0(nn_hm$patient, " ", ifelse(nn_hm$region == "Tumor Core", "TC", "IF"),
                                  "\n", nn_hm$clone),
                           levels = c("P13 TC\nClone1","P13 TC\nClone2","P13 IF\nClone1","P13 IF\nClone2",
                                      "P15 TC\nClone1","P15 TC\nClone2","P15 IF\nClone1","P15 IF\nClone2"))
nn_hm$ligand <- factor(nn_hm$ligand, levels = rev(nn_top_ligs))

nn_p1 <- ggplot(nn_hm, aes(x = cell_label, y = ligand, fill = pearson)) +
  geom_tile(color = NA) +
  scale_fill_viridis_c(option = "viridis", na.value = "grey92", name = "配体活性\n(Pearson r)") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right", legend.key.height = unit(1.1, "cm")) +
  labs(title = "NicheNet 配体活性：克隆特异性发送端 → 区域级 C1QC TAM 受端",
       subtitle = "每格使用该克隆自身的受端基因集；灰色 = 该配体不在该格的候选配体集合内",
       x = NULL, y = NULL)

nn_w <- max(8, 0.55 * 8 + 3); nn_h <- max(6, 0.28 * length(nn_top_ligs) + 3)
ggsave(file.path(nn_fig_dir, "nichenet_ligand_activity_heatmap.pdf"), nn_p1, width = nn_w, height = nn_h, bg = "white")
ggsave(file.path(nn_fig_dir, "nichenet_ligand_activity_heatmap.png"), nn_p1, width = nn_w, height = nn_h, dpi = 300, bg = "white")

cat("[Figure] nichenet_ligand_activity_heatmap | rows(ligands)=", length(nn_top_ligs),
    "| cols=", length(unique(nn_hm$cell_label)), "| size=", nn_w, "x", nn_h, "\n")
cat("[Check] PNG size:", round(file.info(file.path(nn_fig_dir, "nichenet_ligand_activity_heatmap.png"))$size/1e3, 1), "KB\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 图2：Clone1 vs Clone2 配体活性散点（同一病人同一区域，两克隆各自基因集；对角线为活性相等）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

nn_sc <- nn_main[nn_main$in_clone1_sender & nn_main$in_clone2_sender, ]
nn_sc$facet <- factor(paste0(nn_sc$patient, " ", ifelse(nn_sc$region == "Tumor Core", "Tumor Core", "Interface")),
                      levels = c("P13 Tumor Core","P13 Interface","P15 Tumor Core","P15 Interface"))
nn_sc$abs_delta <- abs(nn_sc$delta_pearson)
nn_lab <- do.call(rbind, lapply(split(nn_sc, nn_sc$facet), function(x) x[order(-x$abs_delta), ][1:3, ]))

nn_p2 <- ggplot(nn_sc, aes(x = clone2_pearson, y = clone1_pearson)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = "dotted") +
  geom_point(size = 1.1, alpha = 0.55, color = nn_pal_qual[1]) +
  geom_point(data = nn_lab, size = 2.2, color = nn_pal_qual[4]) +
  geom_text(data = nn_lab, aes(label = ligand), size = 2.6, vjust = -0.6, check_overlap = TRUE) +
  facet_wrap(~ facet, nrow = 2, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(hjust = 0.5, color = "black"),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  labs(title = "同一病人同一区域内 Clone1 与 Clone2 的配体活性对比",
       subtitle = "点 = 两克隆均为候选配体的配体；虚线上方 = Clone1 发送端活性更高；标注 = 每panel |Δ| 最大的 3 个配体",
       x = "Clone2 发送端活性 (Pearson r)", y = "Clone1 发送端活性 (Pearson r)")

ggsave(file.path(nn_fig_dir, "nichenet_clone1_vs_clone2_scatter.pdf"), nn_p2, width = 8, height = 7, bg = "white")
ggsave(file.path(nn_fig_dir, "nichenet_clone1_vs_clone2_scatter.png"), nn_p2, width = 8, height = 7, dpi = 300, bg = "white")
while (dev.cur() > 1) dev.off()

cat("[Figure] nichenet_clone1_vs_clone2_scatter | points=", nrow(nn_sc), "| size=8x7\n")
cat("[Check] PNG size:", round(file.info(file.path(nn_fig_dir, "nichenet_clone1_vs_clone2_scatter.png"))$size/1e3, 1), "KB\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 图3：Clone1 - Clone2 配体活性差（各自基因集下的主对比），每panel 取 |Δ| 最大的 5 个正向与 5 个负向配体
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

nn_d <- nn_main[nn_main$in_clone1_sender & nn_main$in_clone2_sender, ]
nn_d$facet <- factor(paste0(nn_d$patient, " ", ifelse(nn_d$region == "Tumor Core", "Tumor Core", "Interface")),
                     levels = c("P13 Tumor Core","P13 Interface","P15 Tumor Core","P15 Interface"))
nn_d_bar <- do.call(rbind, lapply(split(nn_d, nn_d$facet), function(x) {
  pos <- x[order(-x$delta_pearson), ][1:5, ]
  neg <- x[order(x$delta_pearson), ][1:5, ]
  rbind(pos, neg)
}))
nn_d_bar$ligand <- factor(nn_d_bar$ligand, levels = unique(nn_d_bar$ligand[order(nn_d_bar$facet, nn_d_bar$delta_pearson)]))

nn_p3 <- ggplot(nn_d_bar, aes(x = delta_pearson, y = ligand, fill = delta_pearson > 0)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Clone1 更高", `FALSE` = "Clone2 更高"), name = NULL) +
  facet_wrap(~ facet, nrow = 2, scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(hjust = 0.5, color = "black"),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  labs(title = "配体活性差：Clone1 发送端 − Clone2 发送端",
       subtitle = "每格使用该克隆自身受端基因集；仅统计两克隆均为候选配体的配体；正值(红)=Clone1 更强",
       x = "Δ 配体活性 (Pearson r)", y = NULL)

ggsave(file.path(nn_fig_dir, "nichenet_clone_delta_barplot.pdf"), nn_p3, width = 12, height = 8, bg = "white")
ggsave(file.path(nn_fig_dir, "nichenet_clone_delta_barplot.png"), nn_p3, width = 12, height = 8, dpi = 300, bg = "white")
while (dev.cur() > 1) dev.off()

cat("[Figure] nichenet_clone_delta_barplot | bars=", nrow(nn_d_bar), "| size=12x8\n")
cat("[Check] PNG size:", round(file.info(file.path(nn_fig_dir, "nichenet_clone_delta_barplot.png"))$size/1e3, 1), "KB\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")

# 图4：共享受端口径敏感性——同一基因集下，各配体的活性与其发送端归属（Clone1 独有 / Clone2 独有 / 共有）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

nn_sr <- nn_shared_long
nn_sr$facet <- factor(paste0(nn_sr$patient, " ", ifelse(nn_sr$region == "Tumor Core", "Tumor Core", "Interface")),
                      levels = c("P13 Tumor Core","P13 Interface","P15 Tumor Core","P15 Interface"))
nn_sr$membership <- ifelse(nn_sr$in_clone1_sender & nn_sr$in_clone2_sender, "两克隆共有",
                    ifelse(nn_sr$in_clone1_sender, "Clone1 独有",
                    ifelse(nn_sr$in_clone2_sender, "Clone2 独有", "均不表达")))
nn_sr$membership <- factor(nn_sr$membership, levels = c("Clone1 独有","Clone2 独有","两克隆共有"))

# 每 region 取两套基因集下活性最高的 10 个配体
nn_sr_top <- do.call(rbind, lapply(split(nn_sr, nn_sr$facet), function(x) {
  ligs <- names(sort(tapply(x$pearson, x$ligand, max), decreasing = TRUE))[1:10]
  x[x$ligand %in% ligs, ]
}))
nn_sr_top$ligand <- factor(nn_sr_top$ligand, levels = unique(nn_sr_top$ligand[order(nn_sr_top$facet, -nn_sr_top$pearson)]))
nn_sr_top$geneset_label <- factor(ifelse(nn_sr_top$geneset_source == "Clone1", "Clone1 基因集", "Clone2 基因集"),
                                  levels = c("Clone2 基因集","Clone1 基因集"))

nn_p4 <- ggplot(nn_sr_top, aes(x = ligand, y = geneset_label, size = pmax(pearson, 0), color = membership)) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = c("Clone1 独有" = nn_pal_qual[1], "Clone2 独有" = nn_pal_qual[2],
                                "两克隆共有" = nn_pal_qual[3]), name = "发送端归属") +
  scale_size_continuous(range = c(0.4, 5), name = "活性 (r)") +
  facet_wrap(~ facet, nrow = 2, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(hjust = 0.5, color = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right") +
  labs(title = "共享受端基因集下的配体活性与发送端归属（敏感性分析）",
       subtitle = "每panel 用同一受端基因集评估两克隆的候选配体；点大小 = 活性，颜色 = 该配体在哪个克隆的发送端表达",
       x = NULL, y = NULL)

ggsave(file.path(nn_fig_dir, "nichenet_shared_readout_membership.pdf"), nn_p4, width = 13, height = 8, bg = "white")
ggsave(file.path(nn_fig_dir, "nichenet_shared_readout_membership.png"), nn_p4, width = 13, height = 8, dpi = 300, bg = "white")
while (dev.cur() > 1) dev.off()

cat("[Figure] nichenet_shared_readout_membership | points=", nrow(nn_sr_top), "| size=13x8\n")
cat("[Check] PNG size:", round(file.info(file.path(nn_fig_dir, "nichenet_shared_readout_membership.png"))$size/1e3, 1), "KB\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits = 3), "=====\n")
