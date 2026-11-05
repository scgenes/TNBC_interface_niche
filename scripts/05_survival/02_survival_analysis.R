# =============================================================================
# Bulk survival analysis of the C1Q+ TAM, SPP1+ TAM and MMP11+ CAF signatures
# (TCGA-BRCA, METABRIC, GSE58812, GSE21653; Figure S6)
#
# Converted from the original analysis notebook (step_survival_v2.ipynb,
# 2026-11-05). Cells are kept in their original order. Cells after the audit
# fix (roughly cell 29 onward) supersede the earlier versions; running the
# script top to bottom leaves the objects/tables in the delivered state.
#
# Paths: set TNBC_SURVIVAL_ROOT to the project directory that contains
# brca_tcga/, brca_metabric/, geo/ and processed/survival_v2/. The bulk data
# are described in data/README.md of this repository.
#
# Environment: R 4.3.3; survival 3.8.3; survminer 0.5.0; maxstat 0.7.26;
# GSVA 1.50.0; limma; data.table; ggplot2; patchwork; ggsci.
# =============================================================================

# 环境自检：确认内核路径、R 版本与关键包可用；创建输出目录结构并清掉临时目录标记文件
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

cat("[Startup Config]\n")
cat("  Language: R\n  Kernel Spec: common_r\n  MCP Server: kernel-dispatcher\n")
cat("  Environment:", R.home(), "\n")
cat("  Notebook: project/project_spatial/TNBC/survival/processed/survival_v2/notebooks/step_survival_v2.ipynb\n")
cat("  Output: TNBC_SURVIVAL_ROOT/processed/survival_v2\n")
cat("  Binding: locked\n")
cat("R version:", R.version.string, "\n")
cat(".libPaths():", paste(.libPaths(), collapse=" | "), "\n\n")

pkgs <- c("survival","survminer","maxstat","GSVA","data.table","ggplot2","patchwork","gridExtra","qs","png")
for (p in pkgs) cat(sprintf("  %-12s %s\n", p, if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "MISSING"))

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT, "processed/survival_v2")
for (d in c(OUT, file.path(OUT,"notebooks"), file.path(OUT,"tables"), file.path(OUT,"figures"),
            file.path(OUT,"figures/km"), file.path(OUT,"figures/forest"), file.path(OUT,"data"),
            file.path(OUT,"data/checkpoints"), file.path(OUT,"review"))) {
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
}
# 删除启动前为创建目录而放置的临时标记文件
marker <- file.path(OUT, "notebooks", ".dir_marker")
if (file.exists(marker)) { file.remove(marker); cat("\n[Cleanup] removed", marker, "\n") }

cat("\n[Output Tree]\n"); cat(paste(list.files(OUT, recursive=FALSE), collapse="\n"), "\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 数据探查（修正连接关闭）：打印各临床表完整列名与关键取值，表达表只打印维度与首列
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
p <- function(x) file.path(ROOT, x)

rd <- function(f, skip=0, n=5, sep="\t") {
  con <- if (grepl("\\.gz$", f)) gzfile(f) else file(f)
  d <- read.delim(con, skip=skip, nrows=n, sep=sep, header=TRUE,
                  na.strings=c("NA","","NaN"), check.names=FALSE, quote="")
  try(close(con), silent=TRUE)
  d
}
show <- function(label, f, skip=0, n=5, sep="\t") {
  d <- rd(f, skip=skip, n=n, sep=sep)
  cat(sprintf("\n=== %s ===\n  ncol=%d  peek_rows=%d\n  columns: %s\n",
              label, ncol(d), nrow(d), paste(names(d), collapse=" | ")))
  invisible(d)
}
expr_dim <- function(label, f, sep="\t") {
  d <- rd(f, n=2, sep=sep)
  cat(sprintf("\n=== %s ===\n  expr peek: ncol=%d (first cols: %s)\n",
              label, ncol(d), paste(head(names(d),4), collapse=" | ")))
  invisible(d)
}

d_tcga_sam <- show("TCGA clinical sample", p("brca_tcga/data_clinical_sample.txt"), skip=4)
expr_dim("TCGA expr", p("brca_tcga/data_mrna_seq_v2_rsem.txt"))
show("METABRIC clinical patient", p("brca_metabric/data_clinical_patient.txt"), skip=4)
show("METABRIC clinical sample",  p("brca_metabric/data_clinical_sample.txt"), skip=4)
expr_dim("METABRIC expr", p("brca_metabric/data_mrna_illumina_microarray.txt"))
show("GSE58812 expr", p("geo/GSE58812/GSE58812_expr_genes.tsv.gz"), n=3)
show("GSE58812 clinical", p("geo/GSE58812/GSE58812_clinical.csv"), n=5, sep=",")
show("GSE21653 expr", p("geo/GSE21653/GSE21653_expr_genes.tsv.gz"), n=3)
show("GSE21653 clinical", p("geo/GSE21653/GSE21653_clinical.csv"), n=5, sep=",")

cat("\n--- TCGA 关键列取值 ---\n")
tp <- rd(p("brca_tcga/data_clinical_patient.txt"), skip=4, n=2000)
for (cn in c("ER_STATUS_BY_IHC","PR_STATUS_BY_IHC","IHC_HER2","HER2_FISH_STATUS",
             "OS_STATUS","DFS_STATUS","AJCC_PATHOLOGIC_TUMOR_STAGE","AGE")) {
  if (cn %in% names(tp)) { cat("\n", cn, ":\n"); print(head(sort(table(tp[[cn]], useNA="ifany"), decreasing=TRUE), 8)) }
}
cat("\n--- TCGA sample 关键列取值 ---\n")
for (cn in c("SAMPLE_TYPE","ONCOTREE_CODE","CANCER_TYPE_DETAILED")) {
  if (cn %in% names(d_tcga_sam)) { cat("\n", cn, ":\n"); print(head(sort(table(d_tcga_sam[[cn]], useNA="ifany"), decreasing=TRUE), 8)) }
}

cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 读取四队列临床表，核对 TNBC 定义所需字段取值、终点字段与样本类型，先不做任何筛选
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
rd <- function(f, skip=0, sep="\t") {
  con <- if (grepl("\\.gz$", f)) gzfile(f) else file(f)
  d <- read.delim(con, skip=skip, sep=sep, header=TRUE, na.strings=c("NA","","NaN"),
                  check.names=FALSE, quote="", stringsAsFactors=FALSE)
  try(close(con), silent=TRUE)
  d
}
tp <- rd(file.path(ROOT,"brca_tcga/data_clinical_patient.txt"), skip=4)
ts <- rd(file.path(ROOT,"brca_tcga/data_clinical_sample.txt"), skip=4)
mp <- rd(file.path(ROOT,"brca_metabric/data_clinical_patient.txt"), skip=4)
ms <- rd(file.path(ROOT,"brca_metabric/data_clinical_sample.txt"), skip=4)
g1c <- rd(file.path(ROOT,"geo/GSE58812/GSE58812_clinical.csv"), sep=",")
g2c <- rd(file.path(ROOT,"geo/GSE21653/GSE21653_clinical.csv"), sep=",")

cat(sprintf("TCGA patient rows=%d, sample rows=%d\n", nrow(tp), nrow(ts)))
cat(sprintf("METABRIC patient rows=%d, sample rows=%d\n", nrow(mp), nrow(ms)))
cat(sprintf("GSE58812 clinical rows=%d ; GSE21653 clinical rows=%d\n\n", nrow(g1c), nrow(g2c)))

cat("--- TCGA sample SAMPLE_TYPE ---\n"); print(table(ts$SAMPLE_TYPE, useNA="ifany"))
cat("\n--- METABRIC sample SAMPLE_TYPE ---\n"); print(table(ms$SAMPLE_TYPE, useNA="ifany"))
cat("\n--- METABRIC ER_STATUS ---\n"); print(table(ms$ER_STATUS, useNA="ifany"))
cat("--- METABRIC PR_STATUS ---\n"); print(table(ms$PR_STATUS, useNA="ifany"))
cat("--- METABRIC HER2_STATUS ---\n"); print(table(ms$HER2_STATUS, useNA="ifany"))
cat("\n--- METABRIC OS_STATUS ---\n"); print(table(mp$OS_STATUS, useNA="ifany"))
cat("--- METABRIC RFS_STATUS ---\n"); print(table(mp$RFS_STATUS, useNA="ifany"))
cat("--- METABRIC TUMOR_STAGE ---\n"); print(table(ms$TUMOR_STAGE, useNA="ifany"))
cat("--- METABRIC GRADE ---\n"); print(table(ms$GRADE, useNA="ifany"))
cat("--- METABRIC AGE_AT_DIAGNOSIS summary ---\n"); print(summary(mp$AGE_AT_DIAGNOSIS))

cat("\n--- GSE58812 clinical 前 12 行 ---\n")
print(head(g1c, 12))
cat("\n--- GSE58812 er/pr/her2 ---\n")
for (cn in c("er-ihc","pr-ihc","her2-ihc","death","meta")) { cat(cn,":\n"); print(table(g1c[[cn]], useNA="ifany")) }

cat("\n--- GSE21653 clinical 前 12 行 ---\n")
print(head(g2c, 12))
cat("\n--- GSE21653 er/pr/erbb2/subtype/tissue ---\n")
for (cn in c("er ihc","pr ihc","erbb2","molecular subtype","tissue","dfs evt","sbr grade")) {
  cat(cn,":\n"); print(table(g2c[[cn]], useNA="ifany"))
}

cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 读取四队列表达矩阵并抽取三套 signature 基因（含 HGNC 新旧符号别名），统计每队列基因覆盖
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(data.table); library(png)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")

mem <- as.numeric(strsplit(grep("^MemAvailable", readLines("/proc/meminfo"), value=TRUE), "\\s+")[[1]][2])/1e6
cat(sprintf("[Memory] MemAvailable = %.1f GB\n", mem))

gs <- list(
  C1Q_TAM  = c("PTGDS","GPX1","SEPP1","CHI3L1","ATP6V0C","C1QB","C1QC","PLD3","APOE","ATP5E",
               "C1QA","FYB","HLA-DRB5","ALDOA","CD68","NR1H3","PSAP","TCEB2","CTSC","LIPA"),
  SPP1_TAM = c("MT1G","MT1H","FABP4","SPP1","MT1X","MMP12","MT2A","MT1E","PLIN2","MMP1",
               "MT1F","MT1M","NUPR1","CSTB","CXCL8","FABP5","CTSL","HK2","CCL7","MARCO"),
  MMP11_CAF= c("COL11A1","FN1","MMP1","SELENOM","ATP5F1E","MMP11","COL12A1","TAGLN","COL10A1","INHBA",
               "RACK1","ACTA2","SULF1","CTHRC1","CST1","CTGF","THBS2","AEBP1","TPM1","ELOB")
)
stopifnot(all(lengths(gs)==20))
cat(sprintf("[Gene sets] sizes: %s ; union unique = %d\n",
            paste(names(gs), lengths(gs), sep="=", collapse=", "), length(unique(unlist(gs)))))

alias <- list(
  SEPP1  = c("SEPP1","SELENOP"),
  ATP5E  = c("ATP5E","ATP5F1E"),
  ATP5F1E= c("ATP5F1E","ATP5E"),
  TCEB2  = c("TCEB2","ELOB"),
  ELOB   = c("ELOB","TCEB2"),
  CTGF   = c("CTGF","CCN2"),
  CXCL8  = c("CXCL8","IL8"),
  FYB    = c("FYB","FYB1"),
  RACK1  = c("RACK1","GNB2L1")
)
canon <- unique(unlist(gs))
alias_of <- function(g) if (!is.null(alias[[g]])) alias[[g]] else g

read_expr_fread <- function(f, sep="\t") {
  fread(f, sep=sep, header=TRUE, showProgress=FALSE, na.strings=c("NA","","NaN"),
        check.names=FALSE, data.table=FALSE)
}
pick_genes <- function(d, cohort) {
  sym <- d[[1]]
  mat_all <- as.matrix(d[,-c(1,2), drop=FALSE])
  storage.mode(mat_all) <- "double"
  rownames(mat_all) <- sym
  res <- list(); cov <- data.frame(gene=character(), matched_symbol=character(), stringsAsFactors=FALSE)
  for (g in canon) {
    cand <- alias_of(g); hit <- cand[cand %in% sym]
    if (length(hit)==0) { cov <- rbind(cov, data.frame(gene=g, matched_symbol=NA_character_)); next }
    h <- hit[1]; idx <- which(sym == h)
    v <- if (length(idx)==1) mat_all[idx, ] else colMeans(mat_all[idx, , drop=FALSE], na.rm=TRUE)
    res[[g]] <- v
    cov <- rbind(cov, data.frame(gene=g, matched_symbol=h, stringsAsFactors=FALSE))
  }
  list(mat = do.call(rbind, res), coverage = cov)
}

cat("\n[Load] TCGA expression ...\n")
d_tcga <- read_expr_fread(file.path(ROOT,"brca_tcga/data_mrna_seq_v2_rsem.txt"))
cat(sprintf("  TCGA raw: %d genes x %d sample columns\n", nrow(d_tcga), ncol(d_tcga)-2))
T_tcga <- pick_genes(d_tcga, "TCGA")
cat(sprintf("  TCGA matched genes: %d/%d\n", sum(!is.na(T_tcga$coverage$matched_symbol)), length(canon)))
rm(d_tcga); gc()

cat("\n[Load] GSE58812 expression ...\n")
d1 <- read.delim(gzfile(file.path(ROOT,"geo/GSE58812/GSE58812_expr_genes.tsv.gz")),
                 header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
cat(sprintf("  GSE58812 raw: %d genes x %d samples\n", nrow(d1), ncol(d1)-1))
names(d1)[1] <- "Hugo_Symbol"; d1$Entrez_Gene_Id <- NA
G1 <- pick_genes(d1[, c("Hugo_Symbol","Entrez_Gene_Id", setdiff(names(d1), c("Hugo_Symbol","Entrez_Gene_Id")))], "GSE58812")
rm(d1); gc()

cat("\n[Load] GSE21653 expression ...\n")
d2 <- read.delim(gzfile(file.path(ROOT,"geo/GSE21653/GSE21653_expr_genes.tsv.gz")),
                 header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
cat(sprintf("  GSE21653 raw: %d genes x %d samples\n", nrow(d2), ncol(d2)-1))
names(d2)[1] <- "Hugo_Symbol"; d2$Entrez_Gene_Id <- NA
G2 <- pick_genes(d2[, c("Hugo_Symbol","Entrez_Gene_Id", setdiff(names(d2), c("Hugo_Symbol","Entrez_Gene_Id")))], "GSE21653")
rm(d2); gc()

cat("\n--- 基因覆盖对照 ---\n")
cov_all <- data.frame(gene = canon,
                      TCGA = T_tcga$coverage$matched_symbol[match(canon, T_tcga$coverage$gene)],
                      GSE58812 = G1$coverage$matched_symbol[match(canon, G1$coverage$gene)],
                      GSE21653 = G2$coverage$matched_symbol[match(canon, G2$coverage$gene)],
                      stringsAsFactors = FALSE)
grp <- sapply(canon, function(g) paste(names(gs)[sapply(gs, function(s) g %in% s)], collapse="+"))
cov_all$set <- grp[cov_all$gene]
print(cov_all, row.names=FALSE)
cat(sprintf("\nCoverage: TCGA %d/%d ; GSE58812 %d/%d ; GSE21653 %d/%d\n",
            sum(!is.na(cov_all$TCGA)), nrow(cov_all),
            sum(!is.na(cov_all$GSE58812)), nrow(cov_all),
            sum(!is.na(cov_all$GSE21653)), nrow(cov_all)))

saveRDS(list(T_tcga=T_tcga, G1=G1, G2=G2, gs=gs, canon=canon, alias=alias, cov_all=cov_all),
        file.path(OUT,"data/checkpoints/expr_signature_genes_partial.rds"))
cat("\n[Save] expr_signature_genes_partial.rds\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 排查 TCGA 矩阵中缺失基因的真实符号写法（大小写/别名/是否被拆到别的行）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
f <- file.path(ROOT,"brca_tcga/data_mrna_seq_v2_rsem.txt")
sym <- fread(f, select=1, header=TRUE, showProgress=FALSE)[[1]]
cat("n symbols:", length(sym), "| n unique:", length(unique(sym)), "\n")
cat("first 5:", paste(head(sym,5), collapse=", "), "\n")
cat("any blank:", sum(is.na(sym) | sym==""), "\n")

targets <- c("CHI3L1","APOE","MT1E","SELENOM","ATP5F1E","COL12A1","TPM1",
             "SELENOP","ATP5E","TCEB2","ELOB","GNB2L1","RACK1","IL8","CXCL8","FYB1","CCN2")
for (t in targets) {
  idx <- grep(paste0("^", t, "$"), sym, ignore.case=TRUE)
  cat(sprintf("%-10s exact=%d  caseInsensitive=%d  %s\n", t, sum(sym %in% t), length(idx),
              if (length(idx)) paste(head(sym[idx],3), collapse=",") else "NOT FOUND"))
}
cat("\nAPOE 子串:", paste(head(grep("APOE", sym, value=TRUE),10), collapse=", "), "\n")
cat("MT1 前缀:", paste(head(grep("^MT1", sym, value=TRUE),25), collapse=", "), "\n")
cat("SELENO 子串:", paste(head(grep("SELENO", sym, value=TRUE),10), collapse=", "), "\n")
cat("ATP5 前缀:", paste(head(grep("^ATP5", sym, value=TRUE),20), collapse=", "), "\n")
cat("TPM 前缀:", paste(head(grep("^TPM", sym, value=TRUE),10), collapse=", "), "\n")
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 打印受限：只统计原始文本中是否出现目标基因，不打印整行（每行含上千样本值）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
f <- file.path(ROOT,"brca_tcga/data_mrna_seq_v2_rsem.txt")
ln <- readLines(f, n=21000, warn=FALSE)
cat("lines read:", length(ln), "| line length range:", range(nchar(ln)), "\n")
for (g in c("APOE","CHI3L1","TPM1","COL12A1","MT1E","SELENOM","SEPP1","SELENOP")) {
  hit <- grep(paste0("^", g, "\t"), ln)
  cat(sprintf("  %-10s leading-tab rows: %d  (first line prefix: %s)\n",
              g, length(hit), if (length(hit)) substr(ln[hit[1]], 1, 40) else "-"))
}
cat("\ncol1 vector check (read first column only via read.delim):\n")
con <- file(f)
x1 <- read.delim(con, header=TRUE, nrows=21000, sep="\t", check.names=FALSE,
                 stringsAsFactors=FALSE, quote="", comment.char="", colClasses="character")
try(close(con), silent=TRUE)
cat("read.delim ncol:", ncol(x1), "\n")
for (g in c("APOE","CHI3L1","TPM1","COL12A1","MT1E")) cat(sprintf("  %-8s in col1: %s\n", g, g %in% x1[[1]]))
cat("blank symbols:", sum(is.na(x1[[1]]) | x1[[1]]==""), "\n")
cat("symbols 1000:1005:", paste(x1[[1]][1000:1005], collapse=", "), "\n")
rm(x1); gc()
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 检查 TCGA 表达矩阵是否真的不含这些基因（R 内固定字符串检索，不打印整行）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
f <- file.path(ROOT,"brca_tcga/data_mrna_seq_v2_rsem.txt")
ln <- readLines(f, warn=FALSE)
cat("total lines:", length(ln), "\n")
probe <- c("APOE","CHI3L1","TPM1","COL12A1","MT1E","SELENOM","GAPDH","ACTB","EPCAM","APOC1","C1QA")
for (g in probe) {
  n_any <- sum(grepl(g, ln, fixed=TRUE))
  n_first <- sum(startsWith(ln, paste0(g, "\t")))
  cat(sprintf("  %-9s substring_rows=%4d  first_field_rows=%4d\n", g, n_any, n_first))
}
# 看看被读成第一个字段的基因总数，并抽 20 个含 "APO" 的行首字段
first_field <- sub("\t.*$", "", ln)
cat("\nfirst_field count:", length(first_field), "| unique:", length(unique(first_field)), "\n")
cat("APO* symbols:", paste(head(grep("^APO", first_field, value=TRUE), 20), collapse=", "), "\n")
cat("C1Q* symbols:", paste(head(grep("^C1Q", first_field, value=TRUE), 20), collapse=", "), "\n")
cat("MT1* symbols:", paste(head(grep("^MT1", first_field, value=TRUE), 20), collapse=", "), "\n")
cat("MT2* symbols:", paste(head(grep("^MT2", first_field, value=TRUE), 20), collapse=", "), "\n")
rm(ln); gc()
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 载入 METABRIC 表达矩阵，并按 "原发肿瘤 + 三阴性" 定义四队列的分析样本（打印筛选级联）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
part <- readRDS(file.path(OUT,"data/checkpoints/expr_signature_genes_partial.rds"))
list2env(part, envir=.GlobalEnv)

rd <- function(f, skip=0, sep="\t") {
  con <- if (grepl("\\.gz$", f)) gzfile(f) else file(f)
  d <- read.delim(con, skip=skip, sep=sep, header=TRUE, na.strings=c("NA","","NaN"),
                  check.names=FALSE, quote="", stringsAsFactors=FALSE)
  try(close(con), silent=TRUE); d
}

cat("\n[Load] METABRIC expression (large file) ...\n")
mb <- fread(file.path(ROOT,"brca_metabric/data_mrna_illumina_microarray.txt"),
            sep="\t", header=TRUE, showProgress=FALSE, na.strings=c("NA","","NaN"),
            check.names=FALSE, data.table=TRUE)
cat(sprintf("  METABRIC raw: %d genes x %d sample columns\n", nrow(mb), ncol(mb)-2))
mb_sym <- mb[[1]]
mb_mat_all <- as.matrix(mb[, -c(1,2), with=FALSE]); storage.mode(mb_mat_all) <- "double"
rownames(mb_mat_all) <- mb_sym
cat(sprintf("  value range (3 random genes): %s\n",
            paste(sprintf("%.2f..%.2f", apply(mb_mat_all[sample(nrow(mb_mat_all),3),], 1, min),
                          apply(mb_mat_all[sample(nrow(mb_mat_all),3),], 1, max)), collapse=" | ")))
rm(mb); gc()

res <- list(); covmb <- data.frame(gene=character(), matched_symbol=character(), stringsAsFactors=FALSE)
for (g in canon) {
  cand <- alias_of(g); hit <- cand[cand %in% mb_sym]
  if (length(hit)==0) { covmb <- rbind(covmb, data.frame(gene=g, matched_symbol=NA_character_)); next }
  h <- hit[1]; idx <- which(mb_sym == h)
  v <- if (length(idx)==1) mb_mat_all[idx,] else colMeans(mb_mat_all[idx,,drop=FALSE], na.rm=TRUE)
  res[[g]] <- v; covmb <- rbind(covmb, data.frame(gene=g, matched_symbol=h, stringsAsFactors=FALSE))
}
MB <- list(mat = do.call(rbind, res), coverage = covmb)
cat(sprintf("  METABRIC matched genes: %d/%d\n", sum(!is.na(covmb$matched_symbol)), length(canon)))
rm(mb_mat_all); gc()

add_alias_rows <- function(m) {
  pairs <- list(c("ATP5F1E","ATP5E"), c("ELOB","TCEB2"))
  for (p in pairs) {
    if (!(p[1] %in% rownames(m)) && (p[2] %in% rownames(m))) {
      newrow <- m[p[2], , drop=FALSE]; rownames(newrow) <- p[1]; m <- rbind(m, newrow)
    }
  }
  m
}
T_tcga$mat <- add_alias_rows(T_tcga$mat)
G1$mat <- add_alias_rows(G1$mat)
G2$mat <- add_alias_rows(G2$mat)
MB$mat <- add_alias_rows(MB$mat)
for (nm in c("T_tcga","G1","G2","MB")) {
  cat(sprintf("  %-8s canonical rows after alias fill: %d/%d\n", nm, nrow(get(nm)$mat), length(canon)))
}

cov_all2 <- data.frame(gene = canon,
                       TCGA = T_tcga$coverage$matched_symbol[match(canon, T_tcga$coverage$gene)],
                       METABRIC = covmb$matched_symbol[match(canon, covmb$gene)],
                       GSE58812 = G1$coverage$matched_symbol[match(canon, G1$coverage$gene)],
                       GSE21653 = G2$coverage$matched_symbol[match(canon, G2$coverage$gene)],
                       stringsAsFactors = FALSE)
cov_all2$set <- sapply(canon, function(g) paste(names(gs)[sapply(gs, function(s) g %in% s)], collapse="+"))
for (nm in c("TCGA","METABRIC","GSE58812","GSE21653")) {
  m <- get(switch(nm, TCGA="T_tcga", METABRIC="MB", GSE58812="G1", GSE21653="G2"))$mat
  miss <- is.na(cov_all2[[nm]]) & canon %in% rownames(m)
  cov_all2[[nm]][miss] <- canon[miss]
}
cat("\n--- 四队列基因覆盖（每群可用基因数）---\n")
for (s in names(gs)) {
  sub <- cov_all2[cov_all2$gene %in% gs[[s]], ]
  cat(sprintf("  %-10s TCGA %2d/20 | METABRIC %2d/20 | GSE58812 %2d/20 | GSE21653 %2d/20\n", s,
              sum(!is.na(sub$TCGA)), sum(!is.na(sub$METABRIC)), sum(!is.na(sub$GSE58812)), sum(!is.na(sub$GSE21653))))
}
cat("\n缺失基因明细（按队列）：\n")
for (nm in c("TCGA","METABRIC","GSE58812","GSE21653")) {
  miss <- cov_all2$gene[is.na(cov_all2[[nm]])]
  cat(sprintf("  %-9s (%d): %s\n", nm, length(miss), if (length(miss)) paste(miss, collapse=", ") else "none"))
}

# ---------------- METABRIC / GEO TNBC 队列定义（TCGA 在下一 cell 依患者级受体状态定义）----------------
cat("\n================ 样本筛选级联（METABRIC / GEO）================\n")
mp <- rd(file.path(ROOT,"brca_metabric/data_clinical_patient.txt"), skip=4)
ms <- rd(file.path(ROOT,"brca_metabric/data_clinical_sample.txt"), skip=4)
expr_mb <- colnames(MB$mat)
ms_p <- ms[ms$SAMPLE_TYPE == "Primary", ]
cat(sprintf("\nMETABRIC: 临床样本 %d -> 原发 %d -> 表达矩阵可匹配 %d\n",
            nrow(ms), nrow(ms_p), sum(ms_p$SAMPLE_ID %in% expr_mb)))
tnbc_mb <- ms_p[!is.na(ms_p$ER_STATUS) & !is.na(ms_p$PR_STATUS) & !is.na(ms_p$HER2_STATUS) &
                ms_p$ER_STATUS == "Negative" & ms_p$PR_STATUS == "Negative" & ms_p$HER2_STATUS == "Negative", ]
cat(sprintf("      ER-/PR-/HER2- 三阴 = %d ; 表达可匹配 = %d ; 唯一患者 = %d\n",
            nrow(tnbc_mb), sum(tnbc_mb$SAMPLE_ID %in% expr_mb), length(unique(tnbc_mb$PATIENT_ID))))
tnbc_mb_expr <- tnbc_mb[tnbc_mb$SAMPLE_ID %in% expr_mb, ]
tnbc_mb_expr <- tnbc_mb_expr[!duplicated(tnbc_mb_expr$PATIENT_ID), ]
cat(sprintf("      去重后患者-level 分析集 = %d\n", nrow(tnbc_mb_expr)))

g1c <- rd(file.path(ROOT,"geo/GSE58812/GSE58812_clinical.csv"), sep=",")
cat(sprintf("\nGSE58812: 临床 %d ; 表达可匹配 %d ; er/pr/her2 全 0 = %d\n",
            nrow(g1c), sum(g1c$sample %in% colnames(G1$mat)),
            sum(g1c[["er-ihc"]]==0 & g1c[["pr-ihc"]]==0 & g1c[["her2-ihc"]]==0)))
g1_keep <- g1c[g1c$sample %in% colnames(G1$mat), ]
cat(sprintf("      分析集 = %d ; OS 事件 = %d ; 随访(天) 中位 = %.0f\n",
            nrow(g1_keep), sum(g1_keep$death==1), median(g1_keep[["os (days)"]], na.rm=TRUE)))

g2c <- rd(file.path(ROOT,"geo/GSE21653/GSE21653_clinical.csv"), sep=",")
g2_keep <- g2c[g2c$sample %in% colnames(G2$mat), ]
g2_tn <- g2_keep[g2_keep[["er ihc"]]==0 & g2_keep[["pr ihc"]]==0 & g2_keep[["erbb2"]]==0, ]
cat(sprintf("\nGSE21653: 临床 %d ; 表达可匹配 %d ; 三阴 = %d ; DFS 事件 = %d ; DFS 时间缺失 = %d\n",
            nrow(g2c), nrow(g2_keep), nrow(g2_tn), sum(g2_tn[["dfs evt"]]==1, na.rm=TRUE),
            sum(is.na(g2_tn[["dfs time (months)"]]))))

saveRDS(list(T_tcga=T_tcga, MB=MB, G1=G1, G2=G2, gs=gs, canon=canon, cov_all=cov_all2,
             tnbc_mb=tnbc_mb_expr, g1_keep=g1_keep, g2_tn=g2_tn,
             mb_patient=mp, mb_sample=ms),
        file.path(OUT,"data/checkpoints/cohorts_partial.rds"))
cat("\n[Save] cohorts_partial.rds\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# TCGA 队列定义：比较几种 HER2 判定口径，确定与预期一致的 TNBC 定义后保存分析集
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
rd <- function(f, skip=0, sep="\t") {
  con <- if (grepl("\\.gz$", f)) gzfile(f) else file(f)
  d <- read.delim(con, skip=skip, sep=sep, header=TRUE, na.strings=c("NA","","NaN"),
                  check.names=FALSE, quote="", stringsAsFactors=FALSE)
  try(close(con), silent=TRUE); d
}
tp <- rd(file.path(ROOT,"brca_tcga/data_clinical_patient.txt"), skip=4)
ts <- rd(file.path(ROOT,"brca_tcga/data_clinical_sample.txt"), skip=4)
ts_p <- ts[ts$SAMPLE_TYPE == "Primary", ]
ts_p <- merge(ts_p, tp[, c("PATIENT_ID","ER_STATUS_BY_IHC","PR_STATUS_BY_IHC","IHC_HER2",
                           "HER2_FISH_STATUS","AGE","AJCC_PATHOLOGIC_TUMOR_STAGE",
                           "OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS")],
              by="PATIENT_ID", all.x=TRUE)
expr_cols <- colnames(T_tcga$mat)
er_pr_neg <- !is.na(ts_p$ER_STATUS_BY_IHC) & !is.na(ts_p$PR_STATUS_BY_IHC) &
             ts_p$ER_STATUS_BY_IHC == "Negative" & ts_p$PR_STATUS_BY_IHC == "Negative"

defs <- list(
  A_IHC_neg_only      = er_pr_neg & ts_p$IHC_HER2 == "Negative",
  B_IHCneg_or_FISHneg = er_pr_neg & ((ts_p$IHC_HER2 == "Negative") | (ts_p$HER2_FISH_STATUS == "Negative")),
  C_IHCneg_or_FISHneg_strict = er_pr_neg & (!(ts_p$IHC_HER2 == "Positive") & !(ts_p$HER2_FISH_STATUS == "Positive")) &
                              ((ts_p$IHC_HER2 == "Negative") | (ts_p$HER2_FISH_STATUS == "Negative")),
  D_IHC_neg_or_equiv   = er_pr_neg & ts_p$IHC_HER2 %in% c("Negative","Equivocal") & !(ts_p$HER2_FISH_STATUS=="Positive")
)
cat("TCGA 原发样本数 =", nrow(ts_p), "\n\n不同 HER2 阴性口径下的三阴计数：\n")
for (nm in names(defs)) {
  sel <- ts_p[defs[[nm]], ]
  sel_expr <- sel[sel$SAMPLE_ID %in% expr_cols, ]
  cat(sprintf("  %-26s 样本 %3d | 唯一患者 %3d | 表达可匹配样本 %3d | 去重患者 %3d\n",
              nm, nrow(sel), length(unique(sel$PATIENT_ID)),
              nrow(sel_expr), length(unique(sel_expr$PATIENT_ID))))
}

# 采用标准口径 B（IHC 阴性或 FISH 阴性且任一阳性即排除）
sel <- ts_p[defs$B_IHCneg_or_FISHneg, ]
sel <- sel[sel$SAMPLE_ID %in% expr_cols, ]
sel <- sel[!duplicated(sel$PATIENT_ID), ]
cat(sprintf("\n[Adopted] 口径 B：患者级分析集 = %d\n", nrow(sel)))
os_ok  <- !is.na(sel$OS_MONTHS)  & sel$OS_STATUS  %in% c("0:LIVING","1:DECEASED")
dfs_ok <- !is.na(sel$DFS_MONTHS) & sel$DFS_STATUS %in% c("0:DiseaseFree","1:Recurred/Progressed")
cat(sprintf("  OS 可分析 %d (事件 %d) ; DFS 可分析 %d (事件 %d)\n",
            sum(os_ok), sum(sel$OS_STATUS[os_ok]=="1:DECEASED"),
            sum(dfs_ok), sum(sel$DFS_STATUS[dfs_ok]=="1:Recurred/Progressed")))

base <- readRDS(file.path(OUT,"data/checkpoints/cohorts_defined.rds"))
base$tnbc_tcga <- sel; base$tcga_patient <- tp; base$tcga_sample <- ts
base$tcga_tnbc_def <- "ER-/PR- (IHC) and HER2 negative (IHC Negative, or FISH Negative without any Positive)"
saveRDS(base, file.path(OUT,"data/checkpoints/cohorts_defined.rds"))
cat("\n[Save] cohorts_defined.rds 完成\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 比较更严格的标准 HER2 判定（IHC 0/1+ 为阴性；IHC 2+ 需 FISH 阴性），核对与预期样本量的一致性
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
tp <- read.delim(file.path(ROOT,"brca_tcga/data_clinical_patient.txt"), skip=4, sep="\t",
                 header=TRUE, na.strings=c("NA","","NaN"), check.names=FALSE, quote="", stringsAsFactors=FALSE)
ts <- read.delim(file.path(ROOT,"brca_tcga/data_clinical_sample.txt"), skip=4, sep="\t",
                 header=TRUE, na.strings=c("NA","","NaN"), check.names=FALSE, quote="", stringsAsFactors=FALSE)
ts_p <- merge(ts[ts$SAMPLE_TYPE=="Primary", ],
              tp[, c("PATIENT_ID","ER_STATUS_BY_IHC","PR_STATUS_BY_IHC","IHC_HER2","HER2_FISH_STATUS")],
              by="PATIENT_ID", all.x=TRUE)
er_pr <- ts_p$ER_STATUS_BY_IHC=="Negative" & ts_p$PR_STATUS_BY_IHC=="Negative"
er_pr[is.na(er_pr)] <- FALSE
expr_cols <- colnames(T_tcga$mat)

her2_E <- er_pr & ((ts_p$IHC_HER2=="Negative") |
                   (ts_p$IHC_HER2=="Equivocal" & ts_p$HER2_FISH_STATUS=="Negative"))
her2_F <- er_pr & ((ts_p$IHC_HER2=="Negative") |
                   (ts_p$IHC_HER2=="Equivocal" & ts_p$HER2_FISH_STATUS %in% c("Negative","[Not Available]")))

for (nm in c("E_standard_IHC0_1_plus_FISH_for_2plus","F_E_plus_equivocal_no_FISH")) {
  sel <- ts_p[get(c("E_standard_IHC0_1_plus_FISH_for_2plus"="her2_E","F_E_plus_equivocal_no_FISH"="her2_F")[[nm]]), ]
  selx <- sel[sel$SAMPLE_ID %in% expr_cols, ]
  cat(sprintf("%-38s 样本 %3d | 患者 %3d | 表达可匹配 %3d | 去重患者 %3d\n",
              nm, nrow(sel), length(unique(sel$PATIENT_ID)), nrow(selx), length(unique(selx$PATIENT_ID))))
}
# 交叉表：IHC x FISH（ER/PR 阴性者）
cat("\nER-/PR- 患者的 HER2 IHC x FISH 交叉表：\n")
print(table(IHC=ts_p$IHC_HER2[er_pr], FISH=ts_p$HER2_FISH_STATUS[er_pr], useNA="ifany"))
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 核查三阴患者中是否有男性乳腺癌（常规乳腺癌生存分析多排除男性），并定位与预期的 1 例差异
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
tp <- read.delim(file.path(ROOT,"brca_tcga/data_clinical_patient.txt"), skip=4, sep="\t",
                 header=TRUE, na.strings=c("NA","","NaN"), check.names=FALSE, quote="", stringsAsFactors=FALSE)
ts <- read.delim(file.path(ROOT,"brca_tcga/data_clinical_sample.txt"), skip=4, sep="\t",
                 header=TRUE, na.strings=c("NA","","NaN"), check.names=FALSE, quote="", stringsAsFactors=FALSE)
ts_p <- merge(ts[ts$SAMPLE_TYPE=="Primary", ],
              tp[, c("PATIENT_ID","SEX","ER_STATUS_BY_IHC","PR_STATUS_BY_IHC","IHC_HER2","HER2_FISH_STATUS")],
              by="PATIENT_ID", all.x=TRUE)
er_pr <- ts_p$ER_STATUS_BY_IHC=="Negative" & ts_p$PR_STATUS_BY_IHC=="Negative"; er_pr[is.na(er_pr)] <- FALSE
her2_neg_E <- (ts_p$IHC_HER2=="Negative") | (ts_p$IHC_HER2=="Equivocal" & ts_p$HER2_FISH_STATUS=="Negative")
sel <- ts_p[er_pr & her2_neg_E, ]
selx <- sel[sel$SAMPLE_ID %in% colnames(T_tcga$mat), ]
cat("三阴患者（口径 E）:", nrow(sel), "| 有表达:", nrow(selx), "\n")
cat("性别分布（全部原发）:\n"); print(table(ts_p$SEX, useNA="ifany"))
cat("三阴患者性别:\n"); print(table(sel$SEX, useNA="ifany"))
cat("三阴且有表达者性别:\n"); print(table(selx$SEX, useNA="ifany"))

# METABRIC 与 GEO 也可顺带检查性别/样本类型是否已为女性原发
cat("\nMETABRIC 性别:\n"); print(table(read.delim(file.path(ROOT,"brca_metabric/data_clinical_patient.txt"), skip=4, sep="\t",
      header=TRUE, na.strings=c("NA","",""), check.names=FALSE, quote="", stringsAsFactors=FALSE)$SEX, useNA="ifany"))
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 核对 TCGA 患者表与样本表的 ID 交集，定位与预期 1097 原发的差异来源
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
tp <- read.delim(file.path(ROOT,"brca_tcga/data_clinical_patient.txt"), skip=4, sep="\t",
                 header=TRUE, na.strings=c("NA","","NaN"), check.names=FALSE, quote="", stringsAsFactors=FALSE)
ts <- read.delim(file.path(ROOT,"brca_tcga/data_clinical_sample.txt"), skip=4, sep="\t",
                 header=TRUE, na.strings=c("NA","","NaN"), check.names=FALSE, quote="", stringsAsFactors=FALSE)
cat("patient rows:", nrow(tp), "| unique PATIENT_ID:", length(unique(tp$PATIENT_ID)), "\n")
cat("sample rows:", nrow(ts), "| unique PATIENT_ID:", length(unique(ts$PATIENT_ID)), "\n")
cat("primary samples:", sum(ts$SAMPLE_TYPE=="Primary"), "| unique patients among primary:",
    length(unique(ts$PATIENT_ID[ts$SAMPLE_TYPE=="Primary"])), "\n")
cat("sample patients not in patient table:", length(setdiff(unique(ts$PATIENT_ID), tp$PATIENT_ID)), "\n")
cat("  examples:", paste(head(setdiff(unique(ts$PATIENT_ID), tp$PATIENT_ID), 6), collapse=", "), "\n")
cat("patient table rows with duplicated PATIENT_ID:", sum(duplicated(tp$PATIENT_ID)), "\n")
cat("patient patients not in sample table:", length(setdiff(tp$PATIENT_ID, unique(ts$PATIENT_ID))), "\n")
cat("  examples:", paste(head(setdiff(tp$PATIENT_ID, unique(ts$PATIENT_ID)), 6), collapse=", "), "\n")

# 以患者表为锚（1097）时的三阴计数
ts_p <- ts[ts$SAMPLE_TYPE=="Primary", ]
m <- merge(ts_p, tp[, c("PATIENT_ID","ER_STATUS_BY_IHC","PR_STATUS_BY_IHC","IHC_HER2","HER2_FISH_STATUS")],
           by="PATIENT_ID", all.x=FALSE)
er_pr <- m$ER_STATUS_BY_IHC=="Negative" & m$PR_STATUS_BY_IHC=="Negative"; er_pr[is.na(er_pr)] <- FALSE
her2_neg <- (m$IHC_HER2=="Negative") | (m$IHC_HER2=="Equivocal" & m$HER2_FISH_STATUS=="Negative")
sel <- m[er_pr & her2_neg, ]
selx <- sel[sel$SAMPLE_ID %in% colnames(T_tcga$mat), ]
cat(sprintf("\n以患者表为锚: 匹配到的原发样本 %d ; 三阴样本 %d ; 三阴患者 %d ; 有表达 %d\n",
            nrow(m), nrow(sel), length(unique(sel$PATIENT_ID)), length(unique(selx$PATIENT_ID))))
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 最终确定四队列分析集与终点变量，统一成 cohort×endpoint 长表（含年龄/分期/分级协变量）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
rd <- function(f, skip=0, sep="\t") {
  con <- if (grepl("\\.gz$", f)) gzfile(f) else file(f)
  d <- read.delim(con, skip=skip, sep=sep, header=TRUE, na.strings=c("NA","","NaN"),
                  check.names=FALSE, quote="", stringsAsFactors=FALSE)
  try(close(con), silent=TRUE); d
}
base <- readRDS(file.path(OUT,"data/checkpoints/cohorts_defined.rds"))

# ---- TCGA：采用标准 HER2 口径（IHC 0/1+ 阴性；IHC 2+ 需 FISH 阴性；IHC 3+ 或 FISH 阳性为阳性）
tp <- base$tcga_patient; ts <- base$tcga_sample
ts_p <- ts[ts$SAMPLE_TYPE=="Primary", ]
tp2 <- tp[, c("PATIENT_ID","SEX","ER_STATUS_BY_IHC","PR_STATUS_BY_IHC","IHC_HER2","HER2_FISH_STATUS",
              "AGE","AJCC_PATHOLOGIC_TUMOR_STAGE","OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS")]
tc <- merge(ts_p[, c("PATIENT_ID","SAMPLE_ID")], tp2, by="PATIENT_ID", all.x=FALSE)
er_pr <- tc$ER_STATUS_BY_IHC=="Negative" & tc$PR_STATUS_BY_IHC=="Negative"; er_pr[is.na(er_pr)] <- FALSE
her2_neg <- (tc$IHC_HER2=="Negative") | (tc$IHC_HER2=="Equivocal" & tc$HER2_FISH_STATUS=="Negative")
tc <- tc[er_pr & her2_neg, ]
tc <- tc[tc$SAMPLE_ID %in% colnames(base$T_tcga$mat), ]
tc <- tc[!duplicated(tc$PATIENT_ID), ]
tc$stage_num <- NA_integer_
st <- tc$AJCC_PATHOLOGIC_TUMOR_STAGE
tc$stage_num[grepl("^Stage I($|A|B|C)", st)]  <- 1
tc$stage_num[grepl("^Stage II($|A|B|C)", st)] <- 2
tc$stage_num[grepl("^Stage III($|A|B|C)", st)]<- 3
tc$stage_num[grepl("^Stage IV", st)]          <- 4
tc$stage_high <- ifelse(is.na(tc$stage_num), NA, as.integer(tc$stage_num >= 3))
cat(sprintf("TCGA 三阴分析集: %d ; 分期可用 %d (III/IV %d) ; OS/DFS 事件 %d/%d\n",
            nrow(tc), sum(!is.na(tc$stage_high)), sum(tc$stage_high==1, na.rm=TRUE),
            sum(tc$OS_STATUS=="1:DECEASED", na.rm=TRUE),
            sum(tc$DFS_STATUS=="1:Recurred/Progressed", na.rm=TRUE)))

# ---- METABRIC
mp <- base$mb_patient; ms <- base$mb_sample
mb <- base$tnbc_mb
mb <- merge(mb[, c("PATIENT_ID","SAMPLE_ID","TUMOR_STAGE","GRADE")],
            mp[, c("PATIENT_ID","AGE_AT_DIAGNOSIS","OS_MONTHS","OS_STATUS","RFS_MONTHS","RFS_STATUS")],
            by="PATIENT_ID", all.x=TRUE)
mb$stage_num <- suppressWarnings(as.integer(mb$TUMOR_STAGE)); mb$stage_num[mb$stage_num==0] <- NA
mb$stage_high <- ifelse(is.na(mb$stage_num), NA, as.integer(mb$stage_num >= 3))
cat(sprintf("METABRIC 三阴分析集: %d ; 分期可用 %d ; OS 事件 %d ; RFS 事件 %d\n",
            nrow(mb), sum(!is.na(mb$stage_high)),
            sum(mb$OS_STATUS=="1:DECEASED", na.rm=TRUE),
            sum(mb$RFS_STATUS=="1:Recurred", na.rm=TRUE)))

# ---- GSE58812
g1 <- base$g1_keep
g1$age <- suppressWarnings(as.numeric(g1[["age at diag"]]))
g1$time_os <- suppressWarnings(as.numeric(g1[["os (days)"]])) / 30.4375
g1$event_os <- as.integer(g1$death)
cat(sprintf("GSE58812 分析集: %d ; OS 事件 %d ; 年龄缺失 %d\n",
            nrow(g1), sum(g1$event_os==1), sum(is.na(g1$age))))

# ---- GSE21653
g2 <- base$g2_tn
g2$age <- suppressWarnings(as.numeric(g2[["age at diagnosis"]]))
g2$grade <- suppressWarnings(as.integer(g2[["sbr grade"]]))
g2$time_dfs <- suppressWarnings(as.numeric(g2[["dfs time (months)"]]))
g2$event_dfs <- as.integer(g2[["dfs evt"]])
cat(sprintf("GSE21653 三阴分析集: %d ; DFS 事件 %d ; 年龄缺失 %d ; 分级缺失 %d\n",
            nrow(g2), sum(g2$event_dfs==1, na.rm=TRUE), sum(is.na(g2$age)), sum(is.na(g2$grade))))

# ---- 统一长表
mk <- function(cohort, endpoint, id, time, event, age, stage_high=NULL, grade=NULL) {
  data.frame(cohort=cohort, endpoint=endpoint, sample=id, time=time, event=event, age=age,
             stage_high=if (is.null(stage_high)) NA_integer_ else stage_high,
             grade=if (is.null(grade)) NA_integer_ else grade, stringsAsFactors=FALSE)
}
long <- rbind(
  mk("TCGA","OS", tc$SAMPLE_ID, suppressWarnings(as.numeric(tc$OS_MONTHS)),
     as.integer(tc$OS_STATUS=="1:DECEASED"), suppressWarnings(as.numeric(tc$AGE)), tc$stage_high),
  mk("TCGA","DFS", tc$SAMPLE_ID, suppressWarnings(as.numeric(tc$DFS_MONTHS)),
     as.integer(tc$DFS_STATUS=="1:Recurred/Progressed"), suppressWarnings(as.numeric(tc$AGE)), tc$stage_high),
  mk("METABRIC","OS", mb$SAMPLE_ID, suppressWarnings(as.numeric(mb$OS_MONTHS)),
     as.integer(mb$OS_STATUS=="1:DECEASED"), suppressWarnings(as.numeric(mb$AGE_AT_DIAGNOSIS)), mb$stage_high),
  mk("METABRIC","RFS", mb$SAMPLE_ID, suppressWarnings(as.numeric(mb$RFS_MONTHS)),
     as.integer(mb$RFS_STATUS=="1:Recurred"), suppressWarnings(as.numeric(mb$AGE_AT_DIAGNOSIS)), mb$stage_high),
  mk("GSE58812","OS", g1$sample, g1$time_os, g1$event_os, g1$age),
  mk("GSE21653","DFS", g2$sample, g2$time_dfs, g2$event_dfs, g2$age, grade=g2$grade)
)
long$event[is.na(long$event)] <- NA
cat("\n--- 各 cohort × endpoint 可分析样本与事件数 ---\n")
agg <- aggregate(cbind(n=rep(1,nrow(long)), events=long$event) ~ cohort + endpoint + (time>0 & !is.na(time) & !is.na(event)), data=long, FUN=sum)
print(agg, row.names=FALSE)
write.csv(long, file.path(OUT,"tables/survival_inputs_long.csv"), row.names=FALSE)

# 保存最终对象
base$tcga_analysis <- tc; base$mb_analysis <- mb; base$g1_analysis <- g1; base$g2_analysis <- g2
saveRDS(base, file.path(OUT,"data/checkpoints/analysis_ready.rds"))
cat("\n[Save] analysis_ready.rds + tables/survival_inputs_long.csv\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 载入四队列全转录组矩阵，仅保留三阴样本列；统一到 log2 尺度并保存为打分输入
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
base <- readRDS(file.path(OUT,"data/checkpoints/analysis_ready.rds"))

load_full_dt <- function(f, keep, sep="\t") {
  d <- fread(f, sep=sep, header=TRUE, showProgress=FALSE, na.strings=c("NA","","NaN"),
             check.names=FALSE, data.table=TRUE)
  sym <- d[[1]]
  cols <- intersect(keep, names(d))
  m <- as.matrix(d[, ..cols]); storage.mode(m) <- "double"
  rownames(m) <- sym
  rm(d); gc()
  list(mat=m, matched=length(cols))
}
load_full_txt <- function(f, keep) {
  d <- read.delim(gzfile(f), header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
  sym <- d[[1]]; names(d)[1] <- "Hugo_Symbol"
  cols <- intersect(keep, names(d))
  m <- as.matrix(d[, cols]); storage.mode(m) <- "double"; rownames(m) <- sym
  rm(d); gc()
  list(mat=m, matched=length(cols))
}

cat("\n[1/4] TCGA ...\n")
x <- load_full_dt(file.path(ROOT,"brca_tcga/data_mrna_seq_v2_rsem.txt"), base$tcga_analysis$SAMPLE_ID)
E_tcga <- x$mat
cat(sprintf("  dim=%d x %d ; range %.2f..%.1f\n", nrow(E_tcga), ncol(E_tcga),
            min(E_tcga,na.rm=TRUE), max(E_tcga,na.rm=TRUE)))
E_tcga <- log2(E_tcga + 1)
cat(sprintf("  after log2(x+1): %.2f..%.2f\n", min(E_tcga,na.rm=TRUE), max(E_tcga,na.rm=TRUE)))

cat("\n[2/4] METABRIC ...\n")
x <- load_full_dt(file.path(ROOT,"brca_metabric/data_mrna_illumina_microarray.txt"), base$mb_analysis$SAMPLE_ID)
E_mb <- x$mat
cat(sprintf("  dim=%d x %d ; range %.2f..%.1f (already log2-scale, no transform)\n",
            nrow(E_mb), ncol(E_mb), min(E_mb,na.rm=TRUE), max(E_mb,na.rm=TRUE)))

cat("\n[3/4] GSE58812 ...\n")
x <- load_full_txt(file.path(ROOT,"geo/GSE58812/GSE58812_expr_genes.tsv.gz"), base$g1_analysis$sample)
E_g1 <- x$mat
cat(sprintf("  dim=%d x %d ; range %.2f..%.1f\n", nrow(E_g1), ncol(E_g1),
            min(E_g1,na.rm=TRUE), max(E_g1,na.rm=TRUE)))
if (max(E_g1, na.rm=TRUE) > 50) { E_g1 <- log2(E_g1 + 1); cat("  -> log2(x+1) applied\n") }

cat("\n[4/4] GSE21653 ...\n")
x <- load_full_txt(file.path(ROOT,"geo/GSE21653/GSE21653_expr_genes.tsv.gz"), base$g2_analysis$sample)
E_g2 <- x$mat
cat(sprintf("  dim=%d x %d ; range %.2f..%.1f\n", nrow(E_g2), ncol(E_g2),
            min(E_g2,na.rm=TRUE), max(E_g2,na.rm=TRUE)))
if (max(E_g2, na.rm=TRUE) > 50) { E_g2 <- log2(E_g2 + 1); cat("  -> log2(x+1) applied\n") }

# 去重基因名（保留表达均值最高的一行）；去掉全 NA / 零方差基因
dedup_genes <- function(m) {
  if (anyDuplicated(rownames(m))) {
    mu <- rowMeans(m, na.rm=TRUE)
    ord <- order(rownames(m), -mu)
    m <- m[ord, ]; m <- m[!duplicated(rownames(m)), ]
  }
  keep <- rowSums(!is.na(m)) >= ncol(m)*0.5
  m <- m[keep, , drop=FALSE]
  sds <- apply(m, 1, sd, na.rm=TRUE)
  m <- m[is.finite(sds) & sds > 0, , drop=FALSE]
  m
}
E_tcga <- dedup_genes(E_tcga); E_mb <- dedup_genes(E_mb)
E_g1 <- dedup_genes(E_g1); E_g2 <- dedup_genes(E_g2)
cat(sprintf("\n清理后维度: TCGA %d x %d | METABRIC %d x %d | GSE58812 %d x %d | GSE21653 %d x %d\n",
            nrow(E_tcga), ncol(E_tcga), nrow(E_mb), ncol(E_mb),
            nrow(E_g1), ncol(E_g1), nrow(E_g2), ncol(E_g2)))

saveRDS(list(TCGA=E_tcga, METABRIC=E_mb, GSE58812=E_g1, GSE21653=E_g2),
        file.path(OUT,"data/checkpoints/expr_tnbc_log2.rds"))
cat("[Save] expr_tnbc_log2.rds\n")
gc()
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 计算 signature 打分：主打分队列内 ssGSEA（GSVA，ssgsea.norm=TRUE），敏感性打分基因 z-score 均值
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(GSVA)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
base <- readRDS(file.path(OUT,"data/checkpoints/analysis_ready.rds"))
E <- readRDS(file.path(OUT,"data/checkpoints/expr_tnbc_log2.rds"))
gs <- base$gs

cat("[GSVA] version:", as.character(packageVersion("GSVA")), "\n")

fix_na <- function(m) {
  n_na <- sum(is.na(m))
  if (n_na > 0) {
    cat(sprintf("  NA cells = %d ; imputed by gene row mean\n", n_na))
    mu <- rowMeans(m, na.rm=TRUE)
    idx <- which(is.na(m), arr.ind=TRUE)
    m[idx] <- mu[idx[,1]]
  }
  m
}
run_ssgsea <- function(m, sets) {
  m <- fix_na(m)
  sets <- lapply(sets, function(s) intersect(s, rownames(m)))
  sets <- sets[lengths(sets) >= 5]
  if (length(sets)==0) return(NULL)
  ok <- tryCatch({
    prm <- ssgseaParam(m, sets, normalize=TRUE, alpha=0.25)
    gsva(prm, verbose=FALSE)
  }, error=function(e) {
    cat("  [fallback] old gsva() API:", conditionMessage(e), "\n")
    gsva(m, sets, method="ssgsea", ssgsea.norm=TRUE, verbose=FALSE)
  })
  list(scores=ok, sizes=lengths(sets))
}
run_zmean <- function(m, sets) {
  m <- fix_na(m)
  out <- list(); sizes <- integer()
  for (nm in names(sets)) {
    s <- intersect(sets[[nm]], rownames(m)); sizes[nm] <- length(s)
    if (length(s) < 5) next
    sub <- m[s, , drop=FALSE]
    z <- t(scale(t(sub)))
    out[[nm]] <- colMeans(z, na.rm=TRUE)
  }
  list(scores=do.call(rbind, out), sizes=sizes)
}

res <- list(); cov_rows <- list()
for (co in names(E)) {
  cat(sprintf("\n[%s] %d genes x %d samples\n", co, nrow(E[[co]]), ncol(E[[co]])))
  ss <- run_ssgsea(E[[co]], gs)
  zz <- run_zmean(E[[co]], gs)
  cat(sprintf("  ssGSEA gene counts: %s\n", paste(names(ss$sizes), ss$sizes, sep="=", collapse=", ")))
  cat(sprintf("  z-mean gene counts: %s\n", paste(names(zz$sizes), zz$sizes, sep="=", collapse=", ")))
  cat(sprintf("  ssGSEA score range: %s\n",
              paste(apply(ss$scores, 1, function(v) sprintf("%.3f..%.3f", min(v), max(v))), collapse=" | ")))
  res[[co]] <- list(ssgsea=ss$scores, zmean=zz$scores)
  cov_rows[[co]] <- data.frame(cohort=co, set=names(gs),
                               n_genes_ssgsea=as.integer(ss$sizes[names(gs)]),
                               n_genes_zmean=as.integer(zz$sizes[names(gs)]), stringsAsFactors=FALSE)
}

# 汇总为长表
score_long <- do.call(rbind, lapply(names(res), function(co) {
  do.call(rbind, lapply(c("ssgsea","zmean"), function(meth) {
    s <- res[[co]][[meth]]
    data.frame(cohort=co, method=meth, sample=colnames(s),
               C1Q_TAM=as.numeric(s["C1Q_TAM",]), SPP1_TAM=as.numeric(s["SPP1_TAM",]),
               MMP11_CAF=as.numeric(s["MMP11_CAF",]), stringsAsFactors=FALSE)
  }))
}))
cat("\n--- 打分长表行数 ---\n"); print(table(score_long$cohort, score_long$method))

cat("\n--- 同一队列内 ssGSEA 与 z-mean 的 Spearman 相关 ---\n")
for (co in names(res)) {
  a <- score_long[score_long$cohort==co & score_long$method=="ssgsea", ]
  b <- score_long[score_long$cohort==co & score_long$method=="zmean", ]
  b <- b[match(a$sample, b$sample), ]
  cat(sprintf("  %-9s C1Q %.2f | SPP1 %.2f | MMP11 %.2f\n", co,
              cor(a$C1Q_TAM, b$C1Q_TAM, method="spearman"),
              cor(a$SPP1_TAM, b$SPP1_TAM, method="spearman"),
              cor(a$MMP11_CAF, b$MMP11_CAF, method="spearman")))
}
cov_tab <- do.call(rbind, cov_rows)
write.csv(cov_tab, file.path(OUT,"tables/signature_gene_coverage.csv"), row.names=FALSE)
write.csv(score_long, file.path(OUT,"tables/signature_scores.csv"), row.names=FALSE)
saveRDS(list(scores=res, long=score_long, coverage=cov_tab),
        file.path(OUT,"data/checkpoints/scores.rds"))
cat("\n[Save] tables/signature_gene_coverage.csv ; tables/signature_scores.csv ; scores.rds\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 验证 maxstat 接口与返回结构，并定义后续统一使用的分层/建模辅助函数
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(survival); library(maxstat)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
sc <- read.csv(file.path(OUT,"tables/signature_scores.csv"), stringsAsFactors=FALSE)

# 取 TCGA-OS / ssGSEA / C1Q 做接口测试
d <- long[long$cohort=="TCGA" & long$endpoint=="OS", ]
s <- sc[sc$cohort=="TCGA" & sc$method=="ssgsea", ]
d <- merge(d, s[, c("sample","C1Q_TAM")], by="sample")
d <- d[!is.na(d$time) & !is.na(d$event) & d$time > 0, ]
cat(sprintf("probe data: n=%d events=%d\n", nrow(d), sum(d$event)))
mt <- maxstat.test(Surv(time, event) ~ C1Q_TAM, data=d, smethod="LogRank", pmethod="Lau92", minprop=0.25)
cat("class:", class(mt), "\n")
cat("names:", paste(names(mt), collapse=", "), "\n")
cat(sprintf("estimate(cutpoint)=%.4f  statistic=%.3f  p.value=%.4g\n",
            as.numeric(mt$estimate), as.numeric(mt$statistic), mt$p.value))
cat("mt$stats head:\n"); print(utils::head(mt$stats, 3))
cat("\nmedian of score:", median(d$C1Q_TAM), "\n")

# 通用：分层（maxstat 或 median）
make_groups <- function(dat, score_col, method, minprop=0.25) {
  x <- dat[[score_col]]
  out <- list(cut=NA_real_, p_select=NA_real_, ok=FALSE, method=method)
  if (method == "maxstat") {
    fit <- tryCatch(maxstat.test(Surv(dat$time, dat$event) ~ x, smethod="LogRank",
                                 pmethod="Lau92", minprop=minprop),
                    error=function(e) NULL)
    if (is.null(fit)) return(out)
    out$cut <- as.numeric(fit$estimate); out$p_select <- fit$p.value; out$ok <- TRUE
  } else {
    out$cut <- median(x, na.rm=TRUE); out$ok <- TRUE
  }
  out
}

cat("\n[Helpers defined] make_groups()\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 主分析循环（修正 maxstat 调用：显式构造含 x 列的 data.frame）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(survival); library(maxstat)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
sc   <- read.csv(file.path(OUT,"tables/signature_scores.csv"), stringsAsFactors=FALSE)

make_groups <- function(dat, score_col, method, minprop=0.25) {
  out <- list(cut=NA_real_, p_select=NA_real_, ok=FALSE, method=method, err="")
  if (method == "maxstat") {
    dd <- data.frame(time=dat$time, event=dat$event, x=dat[[score_col]])
    dd <- dd[stats::complete.cases(dd), ]
    fit <- tryCatch(maxstat.test(Surv(time, event) ~ x, data=dd, smethod="LogRank",
                                 pmethod="Lau92", minprop=minprop),
                    error=function(e) { out$err <<- conditionMessage(e); NULL })
    if (is.null(fit)) return(out)
    out$cut <- as.numeric(fit$estimate); out$p_select <- as.numeric(fit$p.value); out$ok <- TRUE
  } else {
    out$cut <- median(dat[[score_col]], na.rm=TRUE); out$ok <- TRUE
  }
  out
}

covars_of <- function(co) switch(co, TCGA=c("age","stage_high"), METABRIC=c("age","stage_high"),
                                 GSE58812=c("age"), GSE21653=c("age","grade"))

fit_stats <- function(dat, covars) {
  dat <- dat[!is.na(dat$time) & !is.na(dat$event) & dat$time > 0 & !is.na(dat$grp), ]
  n <- nrow(dat); ev <- sum(dat$event==1)
  r <- list(n=n, events=ev, n_high=sum(dat$grp=="High"), n_low=sum(dat$grp=="Low"),
            ev_high=sum(dat$event[dat$grp=="High"]), ev_low=sum(dat$event[dat$grp=="Low"]),
            hr=NA_real_, lo=NA_real_, hi=NA_real_, p_wald=NA_real_, logrank_p=NA_real_,
            hr_adj=NA_real_, lo_adj=NA_real_, hi_adj=NA_real_, p_adj=NA_real_,
            n_adj=NA_integer_, n_excluded_cov=NA_integer_, ph_p=NA_real_, note="")
  if (n < 20 || ev < 5) { r$note <- "insufficient (n<20 or events<5)"; return(r) }
  if (min(r$n_high, r$n_low) < 10) { r$note <- "group n<10"; return(r) }
  sd <- survdiff(Surv(time,event) ~ grp, data=dat)
  r$logrank_p <- 1 - pchisq(sd$chisq, df=length(sd$n)-1)
  fit <- coxph(Surv(time,event) ~ grp, data=dat)
  r$hr <- exp(coef(fit)[["grpHigh"]]); ci <- exp(confint(fit))
  r$lo <- ci["grpHigh",1]; r$hi <- ci["grpHigh",2]
  r$p_wald <- summary(fit)$coefficients["grpHigh","Pr(>|z|)"]
  cv <- covars[covars %in% names(dat)]
  if (length(cv) > 0) {
    dat2 <- dat[complete.cases(dat[, c("time","event","grp",cv)]), ]
    r$n_excluded_cov <- n - nrow(dat2); r$n_adj <- nrow(dat2)
    if (nrow(dat2) >= 20 && sum(dat2$event) >= 5 && min(table(dat2$grp)) >= 10) {
      f2 <- tryCatch(coxph(as.formula(paste("Surv(time,event) ~ grp +", paste(cv, collapse="+"))), data=dat2),
                     error=function(e) NULL)
      if (!is.null(f2)) {
        r$hr_adj <- exp(coef(f2)[["grpHigh"]]); ci2 <- exp(confint(f2))
        r$lo_adj <- ci2["grpHigh",1]; r$hi_adj <- ci2["grpHigh",2]
        r$p_adj <- summary(f2)$coefficients["grpHigh","Pr(>|z|)"]
        zph <- cox.zph(f2)$table
        r$ph_p <- if ("grp" %in% rownames(zph)) as.numeric(zph["grp","p"]) else NA_real_
      }
    } else r$note <- paste(r$note, "adjusted model skipped (n/events)", sep="; ")
  }
  r
}

combos <- data.frame(cohort=c("TCGA","TCGA","METABRIC","METABRIC","GSE58812","GSE21653"),
                     endpoint=c("OS","DFS","OS","RFS","OS","DFS"), stringsAsFactors=FALSE)
res_primary <- list(); res_dual2 <- list(); cuts_tab <- list(); grp_tab <- list()
k <- 0; kd <- 0; kc <- 0; kg <- 0

for (i in seq_len(nrow(combos))) {
  co <- combos$cohort[i]; ep <- combos$endpoint[i]
  d0 <- long[long$cohort==co & long$endpoint==ep, ]
  d0 <- d0[!is.na(d0$time) & !is.na(d0$event) & d0$time > 0, ]
  cvs <- covars_of(co)
  for (meth in c("ssgsea","zmean")) {
    s <- sc[sc$cohort==co & sc$method==meth, ]
    d <- merge(d0, s[, c("sample","C1Q_TAM","SPP1_TAM","MMP11_CAF")], by="sample")
    if (nrow(d) < 20) next
    for (cm in c("maxstat","median")) {
      cuts <- lapply(c("C1Q_TAM","SPP1_TAM","MMP11_CAF"), function(g) make_groups(d, g, cm))
      names(cuts) <- c("C1Q_TAM","SPP1_TAM","MMP11_CAF")
      kc <- kc + 1
      cuts_tab[[kc]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm,
                                   signature=names(cuts),
                                   cutpoint=sapply(cuts, function(z) z$cut),
                                   selection_p=sapply(cuts, function(z) z$p_select),
                                   n=nrow(d), events=sum(d$event), stringsAsFactors=FALSE)
      g1 <- ifelse(d$C1Q_TAM > cuts$C1Q_TAM$cut, "High","Low")
      g2 <- ifelse(d$SPP1_TAM > cuts$SPP1_TAM$cut, "High","Low")
      g3 <- ifelse(d$MMP11_CAF > cuts$MMP11_CAF$cut, "High","Low")
      d$C1Q_grp <- factor(g1, levels=c("Low","High"))
      d$SPP1_grp <- factor(g2, levels=c("Low","High"))
      d$MMP11_grp <- factor(g3, levels=c("Low","High"))
      d$dual_grp <- ifelse(g1=="High" & g3=="High", "DualHigh",
                           ifelse(g1=="Low" & g3=="Low", "DualLow", "Other"))
      kg <- kg + 1
      grp_tab[[kg]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm,
                                  sample=d$sample, C1Q_grp=d$C1Q_grp, SPP1_grp=d$SPP1_grp,
                                  MMP11_grp=d$MMP11_grp, dual_grp=d$dual_grp, stringsAsFactors=FALSE)
      specs <- list(
        list(label="C1Q_TAM",    dat=transform(d, grp=C1Q_grp),  p_sel=cuts$C1Q_TAM$p_select),
        list(label="SPP1_TAM",   dat=transform(d, grp=SPP1_grp), p_sel=cuts$SPP1_TAM$p_select),
        list(label="MMP11_CAF",  dat=transform(d, grp=MMP11_grp),p_sel=cuts$MMP11_CAF$p_select),
        list(label="C1Q_MMP11_DualHigh_vs_Other",
             dat=transform(d, grp=factor(ifelse(dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))),
             p_sel=NA_real_)
      )
      for (sp in specs) {
        st <- fit_stats(sp$dat, cvs)
        k <- k + 1
        adj_used <- cm=="maxstat" && !is.na(sp$p_sel)
        pe <- if (adj_used) sp$p_sel else st$logrank_p
        res_primary[[k]] <- data.frame(
          cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm, analysis_group=sp$label,
          n=st$n, events=st$events, n_high=st$n_high, n_low=st$n_low,
          events_high=st$ev_high, events_low=st$ev_low,
          cutpoint_C1Q=cuts$C1Q_TAM$cut, cutpoint_SPP1=cuts$SPP1_TAM$cut, cutpoint_MMP11=cuts$MMP11_CAF$cut,
          HR=st$hr, HR_lower=st$lo, HR_upper=st$hi, p_wald=st$p_wald, logrank_p=st$logrank_p,
          p_effective=pe,
          p_source=if (adj_used) "maxstat Lau92 (selection-adjusted)" else "log-rank (unadjusted)",
          HR_adjusted=st$hr_adj, HR_adj_lower=st$lo_adj, HR_adj_upper=st$hi_adj, p_adjusted=st$p_adj,
          n_adjusted=st$n_adj, n_excluded_missing_covariate=st$n_excluded_cov,
          covariates=paste(cvs, collapse="+"), ph_test_p=st$ph_p, note=st$note, stringsAsFactors=FALSE)
      }
      dsub <- d[d$dual_grp %in% c("DualHigh","DualLow"), ]
      dsub$grp <- factor(ifelse(dsub$dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))
      st2 <- fit_stats(dsub, cvs)
      kd <- kd + 1
      res_dual2[[kd]] <- data.frame(
        cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm,
        analysis_group="C1Q_MMP11_DualHigh_vs_DualLow",
        n=st2$n, events=st2$events, n_high=st2$n_high, n_low=st2$n_low,
        events_high=st2$ev_high, events_low=st2$ev_low,
        HR=st2$hr, HR_lower=st2$lo, HR_upper=st2$hi, p_wald=st2$p_wald, logrank_p=st2$logrank_p,
        HR_adjusted=st2$hr_adj, HR_adj_lower=st2$lo_adj, HR_adj_upper=st2$hi_adj, p_adjusted=st2$p_adj,
        n_adjusted=st2$n_adj, ph_test_p=st2$ph_p, note=st2$note, stringsAsFactors=FALSE)
    }
  }
}

prim <- do.call(rbind, res_primary)
dual2 <- do.call(rbind, res_dual2)
cuts_df <- do.call(rbind, cuts_tab)
grp_df <- do.call(rbind, grp_tab)

prim$p_fdr <- NA_real_
fam <- paste(prim$score_method, prim$cutoff_method, prim$cohort, prim$endpoint, sep="|")
for (f in unique(fam)) {
  idx <- which(fam == f & !is.na(prim$p_effective))
  if (length(idx) > 0) prim$p_fdr[idx] <- p.adjust(prim$p_effective[idx], method="BH")
}
prim$significant_fdr <- !is.na(prim$p_fdr) & prim$p_fdr < 0.05

write.csv(prim,   file.path(OUT,"tables/survival_results_primary.csv"), row.names=FALSE)
write.csv(dual2,  file.path(OUT,"tables/survival_results_dual_secondary.csv"), row.names=FALSE)
write.csv(cuts_df,file.path(OUT,"tables/stratification_cutpoints.csv"), row.names=FALSE)
write.csv(grp_df, file.path(OUT,"tables/group_assignment_long.csv"), row.names=FALSE)
saveRDS(list(primary=prim, dual_secondary=dual2, cuts=cuts_df, groups=grp_df),
        file.path(OUT,"data/checkpoints/survival_results.rds"))

cat(sprintf("\n[Rows] primary=%d ; dual secondary=%d ; cutpoint rows=%d ; group rows=%d\n",
            nrow(prim), nrow(dual2), nrow(cuts_df), nrow(grp_df)))
cat("\n--- maxstat 切点与选择校正 p（ssGSEA）---\n")
cs <- cuts_df[cuts_df$score_method=="ssgsea" & cuts_df$cutoff_method=="maxstat", ]
print(data.frame(cohort=cs$cohort, endpoint=cs$endpoint, signature=cs$signature,
                 cutpoint=round(cs$cutpoint,3), selection_p=signif(cs$selection_p,3)), row.names=FALSE)
cat("\n--- 未拟合记录 ---\n")
bad <- prim[prim$note != "", c("cohort","endpoint","score_method","cutoff_method","analysis_group","n","events","note")]
if (nrow(bad)==0) cat("  none\n") else { print(unique(bad[, c("cohort","endpoint","note")]), row.names=FALSE) }
cat("\n--- 主分析结果（ssGSEA）---\n")
show <- prim[prim$score_method=="ssgsea", c("cohort","endpoint","cutoff_method","analysis_group","n_high","n_low",
                                            "events_high","events_low","HR","HR_lower","HR_upper","logrank_p","p_effective","p_fdr")]
show[, c("HR","HR_lower","HR_upper")] <- round(show[, c("HR","HR_lower","HR_upper")], 3)
show[, c("logrank_p","p_effective","p_fdr")] <- signif(show[, c("logrank_p","p_effective","p_fdr")], 3)
print(show, row.names=FALSE)
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 绘制 KM 生存曲线（主分析 ssGSEA，两种 cutoff，四个分析组），pdf+png 双格式、白底、ASCII 标注
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(survival); library(survminer); library(ggplot2)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
grp  <- read.csv(file.path(OUT,"tables/group_assignment_long.csv"), stringsAsFactors=FALSE)
prim <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)

NPG_COLORS <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000")
dir.create(file.path(OUT,"figures/km"), recursive=TRUE, showWarnings=FALSE)

save_km <- function(pl, stem) {
  pdf(paste0(stem,".pdf"), width=6, height=5.6, bg="white"); print(pl); dev.off()
  png(paste0(stem,".png"), width=1800, height=1680, res=300, bg="white"); print(pl); dev.off()
}

group_meta <- list(
  C1Q_grp  = list(col="C1Q_grp",  lab="C1Q TAMs"),
  SPP1_grp = list(col="SPP1_grp", lab="SPP1 TAMs"),
  MMP11_grp= list(col="MMP11_grp",lab="MMP11 CAFs"),
  dual_grp = list(col="dual_grp", lab="C1Q+ and MMP11+")
)

n_done <- 0; skipped <- c()
combos <- unique(grp[, c("cohort","endpoint","score_method","cutoff_method")])
for (i in seq_len(nrow(combos))) {
  co <- combos$cohort[i]; ep <- combos$endpoint[i]
  if (combos$score_method[i] != "ssgsea") next
  cm <- combos$cutoff_method[i]
  d0 <- long[long$cohort==co & long$endpoint==ep, ]
  g0 <- grp[grp$cohort==co & grp$endpoint==ep & grp$score_method=="ssgsea" & grp$cutoff_method==cm, ]
  d <- merge(d0[, c("sample","time","event")], g0, by="sample")
  d <- d[!is.na(d$time) & !is.na(d$event) & d$time > 0, ]
  pr <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method=="ssgsea" & prim$cutoff_method==cm, ]
  for (gname in names(group_meta)) {
    gm <- group_meta[[gname]]
    dd <- d
    if (gname == "dual_grp") {
      dd <- dd[dd$dual_grp != "Other", ]
      dd$grp <- factor(ifelse(dd$dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))
    } else {
      dd$grp <- factor(dd[[gm$col]], levels=c("Low","High"))
    }
    dd <- dd[!is.na(dd$grp), ]
    if (nrow(dd) < 10 || sum(dd$event) < 5) { skipped <- c(skipped, paste(co, ep, cm, gname, sep="/")); next }
    lab <- if (gname=="dual_grp") "C1Q+MMP11+ vs others" else gm$lab
    row <- pr[pr$analysis_group == if (gname=="dual_grp") "C1Q_MMP11_DualHigh_vs_Other" else
              switch(gname, C1Q_grp="C1Q_TAM", SPP1_grp="SPP1_TAM", MMP11_grp="MMP11_CAF"), ]
    fit <- survfit(Surv(time,event) ~ grp, data=dd)
    sub <- if (nrow(row)==0 || is.na(row$HR[1])) "model not fitted (small groups)" else
      sprintf("HR=%.2f (95%%CI %.2f-%.2f); log-rank p=%.3g; adj p=%.3g",
              row$HR[1], row$HR_lower[1], row$HR_upper[1], row$logrank_p[1], row$p_effective[1])
    pl <- ggsurvplot(fit, data=dd, palette=NPG_COLORS[1:2], conf.int=FALSE,
                     risk.table=TRUE, risk.table.height=0.24, censor=TRUE,
                     legend.title="", legend.labs=c("Low","High"),
                     xlab="Time (months)", ylab="Survival probability",
                     title=sprintf("%s %s - %s (%s cutoff)", co, ep, lab, cm),
                     subtitle=sub, font.subtitle=c(9,"plain"),
                     ggtheme=theme_bw() + theme(panel.grid=element_blank(),
                                                plot.title=element_text(hjust=0.5, face="bold", size=10),
                                                plot.subtitle=element_text(hjust=0.5),
                                                legend.position="right"))
    stem <- file.path(OUT,"figures/km", sprintf("km_%s_%s_%s_%s", co, ep, cm, gname))
    save_km(pl, stem)
    n_done <- n_done + 1
  }
}
cat(sprintf("[KM] saved %d figures (pdf+png pairs)\n", n_done))
cat("[KM] skipped:", if (length(skipped)) paste(skipped, collapse="; ") else "none", "\n")
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 绘制森林图：全组总森林图（未校正）与"未校正+多因素校正"对照森林图
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(ggplot2)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
prim <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)
dir.create(file.path(OUT,"figures/forest"), recursive=TRUE, showWarnings=FALSE)
NPG_COLORS <- c("#E64B35","#4DBBD5","#00A087","#3C5488")

lab_map <- c(C1Q_TAM="C1Q TAMs", SPP1_TAM="SPP1 TAMs", MMP11_CAF="MMP11 CAFs",
             C1Q_MMP11_DualHigh_vs_Other="C1Q+ & MMP11+ (vs other)")
prim$grp_lab <- lab_map[prim$analysis_group]
prim$star <- ifelse(is.na(prim$p_fdr), "", ifelse(prim$p_fdr < 0.001, "***",
                    ifelse(prim$p_fdr < 0.01, "**", ifelse(prim$p_fdr < 0.05, "*", ""))))

forest_panel <- function(df, xlab_txt, title_txt) {
  df <- df[!is.na(df$HR) & !is.na(df$HR_lower), ]
  if (nrow(df)==0) return(NULL)
  df$row <- factor(df$rowlab, levels=rev(unique(df$rowlab)))
  p <- ggplot(df, aes(x=HR, y=row, color=grp_lab)) +
    geom_vline(xintercept=1, linetype="dashed", color="grey40", linewidth=0.4) +
    geom_errorbarh(aes(xmin=HR_lower, xmax=HR_upper), height=0.18, linewidth=0.5) +
    geom_point(shape=15, size=2.4) +
    geom_text(aes(label=paste0(sprintf("%.2f", HR), star)), vjust=-0.9, size=2.9, show.legend=FALSE) +
    scale_x_log10() +
    scale_color_manual(values=setNames(NPG_COLORS[1:length(unique(df$grp_lab))], unique(df$grp_lab))) +
    labs(x=xlab_txt, y=NULL, title=title_txt, color=NULL) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          panel.border=element_rect(colour="black", fill=NA),
          plot.title=element_text(hjust=0.5, face="bold", size=10),
          legend.position="bottom", legend.text=element_text(size=8),
          axis.text.y=element_text(size=8))
  p
}

n1 <- 0; n2 <- 0
for (co in unique(prim$cohort)) for (ep in unique(prim$endpoint[prim$cohort==co])) for (cm in c("maxstat","median")) {
  sub <- prim[prim$cohort==co & prim$endpoint==ep & prim$cutoff_method==cm & prim$score_method=="ssgsea", ]
  if (nrow(sub)==0) next
  # (1) 未校正总森林图
  sub$rowlab <- paste0(sub$grp_lab, "  [n=", sub$n_high + sub$n_low, "]")
  p1 <- forest_panel(sub, "Hazard ratio (high vs low, log scale)",
                     sprintf("%s %s - %s cutoff", co, ep, cm))
  if (!is.null(p1)) {
    stem <- file.path(OUT,"figures/forest", sprintf("forest_all_%s_%s_%s", co, ep, cm))
    ggsave(paste0(stem,".pdf"), p1, width=5.4, height=2.6, bg="white")
    ggsave(paste0(stem,".png"), p1, width=5.4, height=2.6, dpi=300, bg="white")
    n1 <- n1 + 1
  }
  # (2) 未校正 vs 多因素校正
  d2 <- rbind(
    data.frame(grp_lab=sub$grp_lab, HR=sub$HR, HR_lower=sub$HR_lower, HR_upper=sub$HR_upper,
               star=sub$star, model="unadjusted", stringsAsFactors=FALSE),
    data.frame(grp_lab=sub$grp_lab, HR=sub$HR_adjusted, HR_lower=sub$HR_adj_lower, HR_upper=sub$HR_adj_upper,
               star=ifelse(is.na(sub$p_adjusted),"", ifelse(sub$p_adjusted<0.001,"***",
                    ifelse(sub$p_adjusted<0.01,"**", ifelse(sub$p_adjusted<0.05,"*","")))),
               model="adjusted", stringsAsFactors=FALSE)
  )
  d2 <- d2[!is.na(d2$HR) & !is.na(d2$HR_lower), ]
  if (nrow(d2) > 0) {
    d2$rowlab <- paste0(d2$grp_lab, " - ", d2$model)
    p2 <- forest_panel(d2, "Hazard ratio (high vs low, log scale)",
                       sprintf("%s %s - %s cutoff (adjusted: %s)", co, ep, cm,
                               paste(unique(sub$covariates), collapse=",")))
    if (!is.null(p2)) {
      stem <- file.path(OUT,"figures/forest", sprintf("forest_adjusted_%s_%s_%s", co, ep, cm))
      ggsave(paste0(stem,".pdf"), p2, width=5.6, height=3.6, bg="white")
      ggsave(paste0(stem,".png"), p2, width=5.6, height=3.6, dpi=300, bg="white")
      n2 <- n2 + 1
    }
  }
}
cat(sprintf("[Forest] all-group figures=%d ; adjusted-comparison figures=%d\n", n1, n2))

# --- HR 汇总热图（各 cohort-endpoint-cutoff × 4 组，颜色=log2(HR) 发散）---
suppressPackageStartupMessages(library(data.table))
hm <- prim[prim$score_method=="ssgsea" & !is.na(prim$HR), ]
hm$key <- paste(hm$cohort, hm$endpoint, hm$cutoff_method, sep=" | ")
hm$grp_lab <- factor(hm$grp_lab, levels=unname(lab_map))
ph <- ggplot(hm, aes(x=grp_lab, y=key, fill=log2(HR))) +
  geom_tile(color=NA) +
  geom_text(aes(label=ifelse(is.na(star) | star=="", sprintf("%.2f",HR), paste0(sprintf("%.2f",HR), star))),
            size=2.8) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0) +
  labs(x=NULL, y=NULL, fill="log2(HR)", title="HR (high vs low) by cohort, endpoint and cutoff method") +
  theme_bw() +
  theme(panel.grid=element_blank(), panel.border=element_blank(),
        plot.title=element_text(hjust=0.5, face="bold", size=10),
        axis.text.x=element_text(angle=0, size=8), axis.text.y=element_text(size=7))
ggsave(file.path(OUT,"figures/forest/hr_heatmap.pdf"), ph, width=6.4, height=4.2, bg="white")
ggsave(file.path(OUT,"figures/forest/hr_heatmap.png"), ph, width=6.4, height=4.2, dpi=300, bg="white")
cat("[Forest] hr_heatmap saved\n")
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 生成审查证据表：样本流程、PH 检验、方向校验、两种 cutoff 一致性、两种打分一致性
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
prim <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)
dual2<- read.csv(file.path(OUT,"tables/survival_results_dual_secondary.csv"), stringsAsFactors=FALSE)
sc   <- read.csv(file.path(OUT,"tables/signature_scores.csv"), stringsAsFactors=FALSE)
E    <- readRDS(file.path(OUT,"data/checkpoints/expr_tnbc_log2.rds"))
base <- readRDS(file.path(OUT,"data/checkpoints/analysis_ready.rds"))

cat("[Check] GSE21653 三阴队列:", nrow(base$g2_analysis),
    "| 表达矩阵列:", ncol(E$GSE21653),
    "| 交集:", length(intersect(base$g2_analysis$sample, colnames(E$GSE21653))), "\n")
cat("[Check] GSE21653 缺表达样本:", length(setdiff(base$g2_analysis$sample, colnames(E$GSE21653))), "\n")
cat("[Check] GSE21653 有 DFS 时间与事件:", sum(!is.na(base$g2_analysis$time_dfs) & !is.na(base$g2_analysis$event_dfs)), "\n")

# 1) 样本流程
flow <- data.frame(
  cohort = c("TCGA","METABRIC","GSE58812","GSE21653"),
  raw_samples = c(1108, 2509, 107, 266),
  primary = c(1101, 2509, 107, 266),
  with_expression = c(1093, 1980, 107, 266),
  triple_negative = c(143, 320, 107, 96),
  tnbc_with_expression = c(142, 320, 107, 87),
  endpoint_analysable = c("OS 142 (21 events); DFS 134 (21)","OS 319 (168); RFS 319 (131)","OS 107 (29)","DFS 83 (25)"),
  stringsAsFactors=FALSE)
write.csv(flow, file.path(OUT,"tables/sample_flow.csv"), row.names=FALSE)
print(flow, row.names=FALSE)

# 2) PH 检验表
ph <- prim[, c("cohort","endpoint","score_method","cutoff_method","analysis_group","n_adjusted","ph_test_p")]
ph$ph_violation <- !is.na(ph$ph_test_p) & ph$ph_test_p < 0.05
write.csv(ph, file.path(OUT,"tables/ph_assumption_tests.csv"), row.names=FALSE)
cat(sprintf("\n[PH] 模型数=%d ; PH p<0.05 的模型=%d\n", sum(!is.na(ph$ph_test_p)), sum(ph$ph_violation, na.rm=TRUE)))
if (any(ph$ph_violation)) print(ph[which(ph$ph_violation), ], row.names=FALSE) else cat("  none\n")

# 3) 方向校验：score 与其命名基因表达的相关（应显著为正）
anchor <- list(C1Q_TAM="C1QA", SPP1_TAM="SPP1", MMP11_CAF="MMP11")
dir_rows <- list(); k <- 0
for (co in names(E)) for (scn in names(anchor)) {
  g <- anchor[[scn]]
  if (!(g %in% rownames(E[[co]]))) next
  s <- sc[sc$cohort==co & sc$method=="ssgsea", ]
  s <- s[match(colnames(E[[co]]), s$sample), ]
  rho <- suppressWarnings(cor(s[[scn]], as.numeric(E[[co]][g, ]), method="spearman", use="complete.obs"))
  k <- k + 1
  dir_rows[[k]] <- data.frame(cohort=co, signature=scn, anchor_gene=g, spearman_rho=rho,
                              direction_ok=!is.na(rho) && rho > 0, stringsAsFactors=FALSE)
}
dirchk <- do.call(rbind, dir_rows)
write.csv(dirchk, file.path(OUT,"tables/direction_check.csv"), row.names=FALSE)
cat("\n[Direction] score vs namesake gene:\n"); print(dirchk, row.names=FALSE)

# 4) 两种 cutoff 的一致性（ssGSEA & zmean）
cons <- list(); k <- 0
for (co in unique(prim$cohort)) for (ep in unique(prim$endpoint[prim$cohort==co])) for (meth in c("ssgsea","zmean")) for (ag in unique(prim$analysis_group)) {
  a <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method==meth & prim$cutoff_method=="maxstat" & prim$analysis_group==ag, ]
  b <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method==meth & prim$cutoff_method=="median" & prim$analysis_group==ag, ]
  if (nrow(a)==0 || nrow(b)==0) next
  k <- k + 1
  cons[[k]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, analysis_group=ag,
    HR_maxstat=a$HR, p_maxstat=a$p_effective, q_maxstat=a$p_fdr,
    HR_median=b$HR, p_median=b$logrank_p, q_median=b$p_fdr,
    direction_concordant = ifelse(is.na(a$HR) || is.na(b$HR), NA, sign(log(a$HR))==sign(log(b$HR))),
    sig_maxstat = !is.na(a$p_fdr) & a$p_fdr < 0.05,
    sig_median  = !is.na(b$p_fdr) & b$p_fdr < 0.05,
    stringsAsFactors=FALSE)
}
cons <- do.call(rbind, cons)
cons$significance_concordant <- cons$sig_maxstat == cons$sig_median
write.csv(cons, file.path(OUT,"tables/cutoff_consistency.csv"), row.names=FALSE)
cat(sprintf("\n[Cutoff consistency] 可比较组合=%d ; 方向一致=%d ; 显著性一致=%d\n",
            nrow(cons), sum(cons$direction_concordant, na.rm=TRUE), sum(cons$significance_concordant, na.rm=TRUE)))

# 5) 两种打分方法的一致性
sens <- list(); k <- 0
for (co in unique(prim$cohort)) for (ep in unique(prim$endpoint[prim$cohort==co])) for (cm in c("maxstat","median")) for (ag in unique(prim$analysis_group)) {
  a <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method=="ssgsea" & prim$cutoff_method==cm & prim$analysis_group==ag, ]
  b <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method=="zmean" & prim$cutoff_method==cm & prim$analysis_group==ag, ]
  if (nrow(a)==0 || nrow(b)==0) next
  k <- k + 1
  sens[[k]] <- data.frame(cohort=co, endpoint=ep, cutoff_method=cm, analysis_group=ag,
    HR_ssgsea=a$HR, q_ssgsea=a$p_fdr, HR_zmean=b$HR, q_zmean=b$p_fdr,
    direction_concordant=ifelse(is.na(a$HR)||is.na(b$HR), NA, sign(log(a$HR))==sign(log(b$HR))),
    stringsAsFactors=FALSE)
}
sens <- do.call(rbind, sens)
write.csv(sens, file.path(OUT,"tables/score_method_consistency.csv"), row.names=FALSE)
cat(sprintf("[Score-method consistency] 组合=%d ; 方向一致=%d ; 两者同显著=%d\n",
            nrow(sens), sum(sens$direction_concordant, na.rm=TRUE),
            sum(!is.na(sens$q_ssgsea) & !is.na(sens$q_zmean) & sens$q_ssgsea<0.05 & sens$q_zmean<0.05)))

# 6) 显著结果清单
cat("\n--- FDR<0.05 的主分析结果（两种打分方法）---\n")
sig <- prim[!is.na(prim$p_fdr) & prim$p_fdr < 0.05,
            c("cohort","endpoint","score_method","cutoff_method","analysis_group","HR","HR_lower","HR_upper",
              "p_effective","p_fdr","HR_adjusted","p_adjusted")]
if (nrow(sig)==0) cat("  none\n") else {
  sig[, c("HR","HR_lower","HR_upper","HR_adjusted")] <- round(sig[, c("HR","HR_lower","HR_upper","HR_adjusted")],3)
  sig[, c("p_effective","p_fdr","p_adjusted")] <- signif(sig[, c("p_effective","p_fdr","p_adjusted")],3)
  print(sig, row.names=FALSE)
}
cat("\n--- 双高 vs 双低（次要对比）显著结果 ---\n")
ds <- dual2[!is.na(dual2$logrank_p) & dual2$logrank_p < 0.05, ]
print(ds[, c("cohort","endpoint","score_method","cutoff_method","n_high","n_low","HR","HR_lower","HR_upper","logrank_p")], row.names=FALSE)
if (nrow(ds)==0) cat("  none\n")

cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 检查 GSE21653 临床表是否有重复 sample 行，并核对去重后三阴样本数
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
g2c <- read.delim(file.path(ROOT,"geo/GSE21653/GSE21653_clinical.csv"), sep=",", header=TRUE,
                  na.strings=c("NA","","NaN"), check.names=FALSE, quote="", stringsAsFactors=FALSE)
cat("临床行数:", nrow(g2c), "| 唯一 sample:", length(unique(g2c$sample)),
    "| 重复行:", sum(duplicated(g2c$sample)), "\n")
dup <- g2c$sample[duplicated(g2c$sample)]
cat("重复 sample 示例:", paste(head(unique(dup), 10), collapse=", "), "\n")
if (length(dup)) {
  ex <- g2c[g2c$sample %in% head(unique(dup),2), ]
  cat("\n重复样本的完整记录（比较是否完全一致）:\n")
  print(ex[, c("sample","age at diagnosis","dfs evt","dfs time (months)","molecular subtype")], row.names=FALSE)
  cat("\n完全重复行数（整行相同）:", sum(duplicated(g2c)), "\n")
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 定位 GSE21653 打分输入 87 与三阴 96 的差异来源
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
base <- readRDS(file.path(OUT,"data/checkpoints/analysis_ready.rds"))
g2t  <- base$g2_tn; g2a <- base$g2_analysis
cat("g2_tn rows:", nrow(g2t), "| unique sample:", length(unique(g2t$sample)), "\n")
cat("g2_analysis rows:", nrow(g2a), "| unique sample:", length(unique(g2a$sample)), "\n")
cat("g2_analysis duplicated samples:", sum(duplicated(g2a$sample)), "\n")
cat("g2_analysis first 5 sample:", paste(head(g2a$sample,5), collapse=", "), "\n")
E <- readRDS(file.path(OUT,"data/checkpoints/expr_tnbc_log2.rds"))
cat("E$GSE21653 cols:", ncol(E$GSE21653), "\n")
cat("in expr:", sum(unique(g2a$sample) %in% colnames(E$GSE21653)),
    "| not in expr:", sum(!(unique(g2a$sample) %in% colnames(E$GSE21653))), "\n")
cat("not in expr IDs:", paste(setdiff(unique(g2a$sample), colnames(E$GSE21653)), collapse=", "), "\n")
sc <- read.csv(file.path(OUT,"tables/signature_scores.csv"), stringsAsFactors=FALSE)
cat("score table GSE21653 samples:", length(unique(sc$sample[sc$cohort=="GSE21653"])) , "\n")
# 原始临床三阴计数
g2c <- read.delim(file.path(ROOT,"geo/GSE21653/GSE21653_clinical.csv"), sep=",", header=TRUE,
                  na.strings=c("NA","","NaN"), check.names=FALSE, quote="", stringsAsFactors=FALSE)
tn <- g2c[g2c[["er ihc"]]==0 & g2c[["pr ihc"]]==0 & g2c[["erbb2"]]==0, ]
cat("原始临床三阴:", nrow(tn), "| 其中 sample 在 E 中:", sum(tn$sample %in% colnames(E$GSE21653)), "\n")
cat("前 10 个样本是否在 E:", paste(head(tn$sample,10), tn$sample[1:10] %in% colnames(E$GSE21653), collapse=" | "), "\n")
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 检查 GSE21653 临床表中 sample 缺失的行，确认是否为空白/注释行并从分析集中剔除
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
g2c <- read.delim(file.path(ROOT,"geo/GSE21653/GSE21653_clinical.csv"), sep=",", header=TRUE,
                  na.strings=c("NA","","NaN"), check.names=FALSE, quote="", stringsAsFactors=FALSE)
cat("总行:", nrow(g2c), "| sample 缺失:", sum(is.na(g2c$sample)), "\n")
if (any(is.na(g2c$sample))) {
  idx <- which(is.na(g2c$sample))
  cat("缺失 sample 的行号:", paste(head(idx,20), collapse=", "), "\n")
  print(g2c[idx[1:min(3,length(idx))], 1:8], row.names=FALSE)
  cat("\n这些行 er/pr/erbb2 取值:\n")
  print(table(g2c[["er ihc"]][idx], g2c[["pr ihc"]][idx], useNA="ifany"))
}
# 原始文件的最后几行（不打印长内容）
con <- file(file.path(ROOT,"geo/GSE21653/GSE21653_clinical.csv"))
ln <- readLines(con); close(con)
cat("\n文件总行数（含表头）:", length(ln), "| 最后 3 行:\n")
cat(paste(tail(ln, 3), collapse="\n"), "\n")

# 修正后的 GSE21653 三阴计数
g2_ok <- g2c[!is.na(g2c$sample), ]
tn <- g2_ok[g2_ok[["er ihc"]]==0 & g2_ok[["pr ihc"]]==0 & g2_ok[["erbb2"]]==0, ]
cat(sprintf("\n修正后: 有效临床行 %d ; 三阴 %d ; 唯一 sample %d\n",
            nrow(g2_ok), nrow(tn), length(unique(tn$sample))))
cat("DFS 时间与事件均可得:", sum(!is.na(tn[["dfs time (months)"]]) & !is.na(tn[["dfs evt"]])), "\n")
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 修正 GSE21653 三阴集：剔除因受体缺失而由 NA 条件产生的空行（不影响已建模结果，仅纠正样本计数）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
base <- readRDS(file.path(OUT,"data/checkpoints/analysis_ready.rds"))

cat("修正前 GSE21653 g2_tn 行数:", nrow(base$g2_tn), "| 空 sample 行:", sum(is.na(base$g2_tn$sample)), "\n")
base$g2_tn <- base$g2_tn[!is.na(base$g2_tn$sample), ]
base$g2_analysis <- base$g2_analysis[!is.na(base$g2_analysis$sample), ]
cat("修正后 GSE21653 三阴:", nrow(base$g2_tn), "| 唯一:", length(unique(base$g2_tn$sample)), "\n")
cat("DFS 可分析:", sum(!is.na(base$g2_analysis$time_dfs) & !is.na(base$g2_analysis$event_dfs) & base$g2_analysis$time_dfs > 0),
    "| 事件:", sum(base$g2_analysis$event_dfs==1, na.rm=TRUE), "\n")
saveRDS(base, file.path(OUT,"data/checkpoints/analysis_ready.rds"))

# 清理长表中的空 sample 行
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
cat("\n长表空 sample 行:", sum(is.na(long$sample)), "\n")
long <- long[!is.na(long$sample), ]
write.csv(long, file.path(OUT,"tables/survival_inputs_long.csv"), row.names=FALSE)

# 重建样本流程表（GSE21653 三阴修正为 87）
flow <- data.frame(
  cohort = c("TCGA","METABRIC","GSE58812","GSE21653"),
  raw_samples = c(1108, 2509, 107, 266),
  primary = c(1101, 2509, 107, 266),
  with_expression = c(1093, 1980, 107, 266),
  triple_negative = c(143, 320, 107, 87),
  tnbc_with_expression = c(142, 320, 107, 87),
  endpoint_analysable = c("OS 142 (21 events); DFS 134 (21)","OS 319 (168); RFS 319 (131)","OS 107 (29)","DFS 83 (25)"),
  stringsAsFactors=FALSE)
write.csv(flow, file.path(OUT,"tables/sample_flow.csv"), row.names=FALSE)
print(flow, row.names=FALSE)

# 核对各 cohort × endpoint 的最终 n 与事件数（从长表重算）
cat("\n--- 最终 cohort × endpoint 分析样本 ---\n")
for (co in unique(long$cohort)) for (ep in unique(long$endpoint[long$cohort==co])) {
  d <- long[long$cohort==co & long$endpoint==ep, ]
  ok <- !is.na(d$time) & !is.na(d$event) & d$time > 0
  cat(sprintf("  %-9s %-3s n=%3d events=%3d\n", co, ep, sum(ok), sum(d$event[ok]==1)))
}
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 生成完整多因素 Cox 系数表（主分析 ssGSEA），并汇总关键结果供报告引用
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(survival))
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
grp  <- read.csv(file.path(OUT,"tables/group_assignment_long.csv"), stringsAsFactors=FALSE)
prim <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)
covars_of <- function(co) switch(co, TCGA=c("age","stage_high"), METABRIC=c("age","stage_high"),
                                 GSE58812=c("age"), GSE21653=c("age","grade"))

rows <- list(); k <- 0
for (i in seq_len(nrow(prim))) {
  r <- prim[i, ]
  if (r$score_method != "ssgsea" || is.na(r$HR_adjusted)) next
  d0 <- long[long$cohort==r$cohort & long$endpoint==r$endpoint, ]
  g0 <- grp[grp$cohort==r$cohort & grp$endpoint==r$endpoint & grp$score_method=="ssgsea" & grp$cutoff_method==r$cutoff_method, ]
  d <- merge(d0[, c("sample","time","event","age","stage_high","grade")], g0, by="sample")
  if (r$analysis_group == "C1Q_MMP11_DualHigh_vs_Other") {
    d <- d[d$dual_grp != "Other", ]; d$grp <- factor(ifelse(d$dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))
  } else {
    col <- switch(r$analysis_group, C1Q_TAM="C1Q_grp", SPP1_TAM="SPP1_grp", MMP11_CAF="MMP11_grp")
    d$grp <- factor(d[[col]], levels=c("Low","High"))
  }
  cv <- covars_of(r$cohort); cv <- cv[cv %in% names(d)]
  d <- d[complete.cases(d[, c("time","event","grp",cv)]), ]
  fit <- tryCatch(coxph(as.formula(paste("Surv(time,event) ~ grp +", paste(cv, collapse="+"))), data=d), error=function(e) NULL)
  if (is.null(fit)) next
  s <- summary(fit)
  for (term in rownames(s$coefficients)) {
    k <- k + 1
    rows[[k]] <- data.frame(cohort=r$cohort, endpoint=r$endpoint, cutoff_method=r$cutoff_method,
      analysis_group=r$analysis_group, term=term, n=nrow(d), events=sum(d$event),
      HR=s$conf.int[term,"exp(coef)"], HR_lower=s$conf.int[term,"lower .95"],
      HR_upper=s$conf.int[term,"upper .95"], p=s$coefficients[term,"Pr(>|z|)"], stringsAsFactors=FALSE)
  }
}
coxfull <- do.call(rbind, rows)
coxfull$q_within_model <- ave(coxfull$p, paste(coxfull$cohort, coxfull$endpoint, coxfull$cutoff_method,
                                               coxfull$analysis_group), FUN=function(x) p.adjust(x, "BH"))
write.csv(coxfull, file.path(OUT,"tables/cox_multivariable_full.csv"), row.names=FALSE)
cat(sprintf("[Cox full] rows=%d ; models=%d\n", nrow(coxfull), length(unique(paste(coxfull$cohort,coxfull$endpoint,coxfull$cutoff_method,coxfull$analysis_group)))))

cat("\n--- 主分析 ssGSEA 全部结果（HR / log-rank p / 校正 p / FDR q / 校正HR）---\n")
p <- prim[prim$score_method=="ssgsea", c("cohort","endpoint","cutoff_method","analysis_group",
                                         "n_high","n_low","events_high","events_low",
                                         "HR","HR_lower","HR_upper","logrank_p","p_effective","p_fdr",
                                         "HR_adjusted","p_adjusted")]
p[, c("HR","HR_lower","HR_upper","HR_adjusted")] <- round(p[, c("HR","HR_lower","HR_upper","HR_adjusted")],2)
p[, c("logrank_p","p_effective","p_fdr","p_adjusted")] <- signif(p[, c("logrank_p","p_effective","p_fdr","p_adjusted")],3)
print(p, row.names=FALSE)

cat("\n--- 敏感性打分 zmean：与 ssGSEA 结果对照（仅列 |HR|>1.3 或 q<0.1 的行）---\n")
q <- prim[prim$score_method=="zmean" & (!is.na(prim$p_fdr) & prim$p_fdr < 0.1 | !is.na(prim$HR) & (prim$HR > 1.3 | prim$HR < 0.77)),
          c("cohort","endpoint","cutoff_method","analysis_group","HR","logrank_p","p_effective","p_fdr","HR_adjusted")]
q[, c("HR","HR_adjusted")] <- round(q[, c("HR","HR_adjusted")],2)
q[, c("logrank_p","p_effective","p_fdr")] <- signif(q[, c("logrank_p","p_effective","p_fdr")],3)
print(q, row.names=FALSE)

cat("\n--- 切点（ssGSEA, maxstat & median）---\n")
cuts <- read.csv(file.path(OUT,"tables/stratification_cutpoints.csv"), stringsAsFactors=FALSE)
cs <- cuts[cuts$score_method=="ssgsea", ]
cs$cutpoint <- round(cs$cutpoint,3)
print(reshape(cs[, c("cohort","endpoint","cutoff_method","signature","cutpoint")],
              idvar=c("cohort","endpoint","signature"), timevar="cutoff_method", direction="wide"), row.names=FALSE)
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 图件完整性自检：成对性、文件大小、PNG 非空白（关键图全检 + 其余抽样 20%）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(png))
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
pdfs <- list.files(file.path(OUT,"figures"), pattern="\\.pdf$", recursive=TRUE, full.names=TRUE)
pngs <- list.files(file.path(OUT,"figures"), pattern="\\.png$", recursive=TRUE, full.names=TRUE)
cat(sprintf("[Figures] pdf=%d ; png=%d\n", length(pdfs), length(pngs)))
stem <- function(x) sub("\\.(pdf|png)$", "", x)
cat("PDF 无对应 PNG:", length(setdiff(stem(pdfs), stem(pngs))), "| PNG 无对应 PDF:", length(setdiff(stem(pngs), stem(pdfs))), "\n")
sizes <- file.info(pngs)$size
cat(sprintf("[PNG size] min=%.0f KB ; median=%.0f KB ; <10KB 数量=%d\n",
            min(sizes)/1024, median(sizes)/1024, sum(sizes < 10*1024)))

key <- c("km_METABRIC_OS_maxstat_C1Q_grp.png","km_METABRIC_OS_median_C1Q_grp.png",
         "km_METABRIC_RFS_maxstat_C1Q_grp.png","km_METABRIC_OS_maxstat_dual_grp.png",
         "km_GSE21653_DFS_median_C1Q_grp.png","km_TCGA_OS_median_C1Q_grp.png",
         "km_GSE58812_OS_maxstat_MMP11_grp.png","hr_heatmap.png",
         "forest_all_METABRIC_OS_maxstat.png","forest_adjusted_METABRIC_OS_maxstat.png")
key_paths <- c(file.path(OUT,"figures/km", key[1:7]), file.path(OUT,"figures/forest", key[8:10]))
key_paths <- key_paths[file.exists(key_paths)]
rest <- setdiff(pngs, key_paths)
set.seed(42); samp <- sample(rest, max(3, ceiling(length(rest)*0.2)))
check <- unique(c(key_paths, samp))

measure_one <- function(f) {
  img <- tryCatch(readPNG(f), error=function(e) NULL)
  if (is.null(img)) return(NULL)
  rgb <- if (length(dim(img))==3) img[,,1:3] else array(img, c(dim(img),3))
  v <- as.vector(rgb)
  data.frame(file=basename(f), sd=sd(v), ink_frac=mean(v < 0.98))
}
res <- do.call(rbind, lapply(check, measure_one))
res$pass <- res$sd > 0.02 & res$ink_frac > 0.001
print(res[order(res$sd), ], row.names=FALSE)
cat(sprintf("\n[Pixel check] %d/%d 通过（sd>0.02 且非白通道占比>0.1%%）\n", sum(res$pass), nrow(res)))
write.csv(res, file.path(OUT,"review/figure_pixel_check.csv"), row.names=FALSE)
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 交付清单核对：表格行数、图件数量、检查点与汇总文档是否存在
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
tabs <- c("sample_flow.csv","signature_scores.csv","signature_gene_coverage.csv","survival_inputs_long.csv",
          "stratification_cutpoints.csv","survival_results_primary.csv","survival_results_dual_secondary.csv",
          "cox_multivariable_full.csv","ph_assumption_tests.csv","direction_check.csv",
          "cutoff_consistency.csv","score_method_consistency.csv","group_assignment_long.csv")
cat("--- tables/ ---\n")
for (t in tabs) {
  f <- file.path(OUT,"tables",t)
  if (file.exists(f)) {
    n <- length(readLines(f)) - 1
    cat(sprintf("  %-38s %5d 行  %7.1f KB\n", t, n, file.info(f)$size/1024))
  } else cat(sprintf("  %-38s MISSING\n", t))
}
cat("\n--- data/checkpoints/ ---\n")
for (f in list.files(file.path(OUT,"data/checkpoints"), full.names=TRUE)) {
  cat(sprintf("  %-38s %7.1f MB\n", basename(f), file.info(f)$size/1e6))
}
cat("\n--- 交付文档 ---\n")
doc <- file.path(OUT,"step_survival_v2_汇总.html")
cat(sprintf("  %s : exists=%s, %.1f KB\n", basename(doc), file.exists(doc), file.info(doc)$size/1024))

cat("\n--- 图件 ---\n")
cat("  km:", length(list.files(file.path(OUT,"figures/km"), pattern="\\.pdf$")), "pdf /",
    length(list.files(file.path(OUT,"figures/km"), pattern="\\.png$")), "png\n")
cat("  forest:", length(list.files(file.path(OUT,"figures/forest"), pattern="\\.pdf$")), "pdf /",
    length(list.files(file.path(OUT,"figures/forest"), pattern="\\.png$")), "png\n")

cat("\n--- 表格行数与分组 n 对账（主结果表）---\n")
p <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)
cat("  survival_results_primary.csv 行数:", nrow(p),
    "| 拟合成功(HR 非 NA):", sum(!is.na(p$HR)),
    "| 未拟合:", sum(is.na(p$HR)), "\n")
cat("  未拟合原因:\n"); print(table(p$note[p$note != ""]))
cat("  FDR<0.05 行数:", sum(!is.na(p$p_fdr) & p$p_fdr < 0.05), "\n")

cat("\n[Exit Check]\n")
chk <- data.frame(
  item = c("中间对象(4 个 checkpoint)","核心表格(12 个)","汇总 HTML","KM 图 48 组","森林图 12+12+1"),
  ok = c(length(list.files(file.path(OUT,"data/checkpoints"), pattern="\\.rds$"))==4,
         all(file.exists(file.path(OUT,"tables",tabs))),
         file.exists(doc),
         length(list.files(file.path(OUT,"figures/km"), pattern="\\.pdf$"))==48,
         length(list.files(file.path(OUT,"figures/forest"), pattern="\\.pdf$"))==25)
)
print(chk, row.names=FALSE)
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 修复前的三项核查：旧 signature 分数表结构、METABRIC 唯一 NA 的位置、严格 HER2 口径下的 TCGA 三阴计数
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")

# 1) 旧版分数表
old <- read.csv(file.path(ROOT,"processed/scores.csv"), stringsAsFactors=FALSE)
cat("[old scores] dim:", dim(old), "\n")
cat("columns:", paste(names(old), collapse=" | "), "\n")
cat("head:\n"); print(head(old, 2), row.names=FALSE)

# 2) METABRIC NA 位置
E <- readRDS(file.path(OUT,"data/checkpoints/expr_tnbc_log2.rds"))
mb <- E$METABRIC
na_idx <- which(is.na(mb), arr.ind=TRUE)
cat(sprintf("\n[METABRIC NA] cells=%d\n", nrow(na_idx)))
if (nrow(na_idx)) {
  for (i in seq_len(nrow(na_idx))) {
    cat(sprintf("  gene=%s  sample=%s  value=NA\n", rownames(mb)[na_idx[i,1]], colnames(mb)[na_idx[i,2]]))
  }
  g <- rownames(mb)[na_idx[1,1]]
  cat("  gene in any signature set:", paste(names(which(sapply(E, function(x) TRUE)))[1]), "check:",
      g %in% unlist(readRDS(file.path(OUT,"data/checkpoints/cohorts_defined.rds"))$gs), "\n")
}
cat("other cohorts NA cells:", sapply(E[c("TCGA","GSE58812","GSE21653")], function(m) sum(is.na(m))), "\n")

# 3) 严格 HER2 口径
tp <- read.delim(file.path(ROOT,"brca_tcga/data_clinical_patient.txt"), skip=4, sep="\t", header=TRUE,
                 na.strings=c("NA","","NaN"), check.names=FALSE, quote="", stringsAsFactors=FALSE)
ts <- read.delim(file.path(ROOT,"brca_tcga/data_clinical_sample.txt"), skip=4, sep="\t", header=TRUE,
                 na.strings=c("NA","","NaN"), check.names=FALSE, quote="", stringsAsFactors=FALSE)
ts_p <- merge(ts[ts$SAMPLE_TYPE=="Primary", c("PATIENT_ID","SAMPLE_ID")],
              tp[, c("PATIENT_ID","ER_STATUS_BY_IHC","PR_STATUS_BY_IHC","IHC_HER2","HER2_FISH_STATUS")],
              by="PATIENT_ID", all.x=FALSE)
er_pr <- ts_p$ER_STATUS_BY_IHC=="Negative" & ts_p$PR_STATUS_BY_IHC=="Negative"; er_pr[is.na(er_pr)] <- FALSE
her2_pos <- ts_p$IHC_HER2 %in% c("Positive") | ts_p$HER2_FISH_STATUS %in% c("Positive")
her2_strict <- !her2_pos & (ts_p$IHC_HER2 == "Negative" |
                            (ts_p$IHC_HER2 %in% c("Equivocal","Indeterminate") & ts_p$HER2_FISH_STATUS == "Negative"))
sel <- ts_p[er_pr & her2_strict, ]
expr_cols <- rownames(E$TCGA)  # placeholder replaced below
expr_cols <- readRDS(file.path(OUT,"data/checkpoints/cohorts_defined.rds"))$T_tcga$mat |> colnames()
cat(sprintf("\n[STRICT HER2] primary patients=%d ; triple-negative=%d ; with expression=%d ; unique patient=%d\n",
            nrow(ts_p), nrow(sel), sum(sel$SAMPLE_ID %in% expr_cols), length(unique(sel$SAMPLE_ID[sel$SAMPLE_ID %in% expr_cols]))))
cat("IHC x FISH among ER-/PR-:\n"); print(table(IHC=ts_p$IHC_HER2[er_pr], FISH=ts_p$HER2_FISH_STATUS[er_pr], useNA="ifany"))
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 修复 1/9：按严格 HER2 口径重建 TCGA 三阴队列（IHC 阳性或 FISH 阳性即排除），并重建分析长表
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
base <- readRDS(file.path(OUT,"data/checkpoints/cohorts_defined.rds"))

tp <- base$tcga_patient; ts <- base$tcga_sample
ts_p <- merge(ts[ts$SAMPLE_TYPE=="Primary", c("PATIENT_ID","SAMPLE_ID")],
              tp[, c("PATIENT_ID","SEX","ER_STATUS_BY_IHC","PR_STATUS_BY_IHC","IHC_HER2","HER2_FISH_STATUS",
                     "AGE","AJCC_PATHOLOGIC_TUMOR_STAGE","OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS")],
              by="PATIENT_ID", all.x=FALSE)
er_pr <- ts_p$ER_STATUS_BY_IHC=="Negative" & ts_p$PR_STATUS_BY_IHC=="Negative"; er_pr[is.na(er_pr)] <- FALSE
her2_pos <- ts_p$IHC_HER2 %in% "Positive" | ts_p$HER2_FISH_STATUS %in% "Positive"
her2_neg <- !her2_pos & (ts_p$IHC_HER2 == "Negative" |
                         (ts_p$IHC_HER2 %in% c("Equivocal","Indeterminate") & ts_p$HER2_FISH_STATUS == "Negative"))
tc <- ts_p[er_pr & her2_neg, ]
cat(sprintf("[TCGA strict] primary=%d ; triple-negative=%d ; with expression=%d\n",
            nrow(ts_p), nrow(tc), sum(tc$SAMPLE_ID %in% colnames(base$T_tcga$mat))))
tc <- tc[tc$SAMPLE_ID %in% colnames(base$T_tcga$mat), ]
tc <- tc[!duplicated(tc$PATIENT_ID), ]
tc$stage_num <- NA_integer_
st <- tc$AJCC_PATHOLOGIC_TUMOR_STAGE
tc$stage_num[grepl("^Stage I($|A|B|C)", st)]   <- 1
tc$stage_num[grepl("^Stage II($|A|B|C)", st)]  <- 2
tc$stage_num[grepl("^Stage III($|A|B|C)", st)] <- 3
tc$stage_num[grepl("^Stage IV", st)]           <- 4
tc$stage_high <- ifelse(is.na(tc$stage_num), NA, as.integer(tc$stage_num >= 3))
cat(sprintf("  analysis set=%d ; OS events=%d ; DFS events=%d ; stage available=%d\n",
            nrow(tc), sum(tc$OS_STATUS=="1:DECEASED", na.rm=TRUE),
            sum(tc$DFS_STATUS=="1:Recurred/Progressed", na.rm=TRUE), sum(!is.na(tc$stage_high))))

# METABRIC / GEO 从队列对象重建（保持原定义）
mp <- base$mb_patient
mb <- merge(base$tnbc_mb[, c("PATIENT_ID","SAMPLE_ID","TUMOR_STAGE","GRADE")],
            mp[, c("PATIENT_ID","AGE_AT_DIAGNOSIS","OS_MONTHS","OS_STATUS","RFS_MONTHS","RFS_STATUS")],
            by="PATIENT_ID", all.x=TRUE)
mb$stage_num <- suppressWarnings(as.integer(mb$TUMOR_STAGE)); mb$stage_num[mb$stage_num==0] <- NA
mb$stage_high <- ifelse(is.na(mb$stage_num), NA, as.integer(mb$stage_num >= 3))
g1 <- base$g1_keep
g1$age <- suppressWarnings(as.numeric(g1[["age at diag"]]))
g1$time_os <- suppressWarnings(as.numeric(g1[["os (days)"]])) / 30.4375
g1$event_os <- as.integer(g1$death)
g2 <- base$g2_tn
g2$age <- suppressWarnings(as.numeric(g2[["age at diagnosis"]]))
g2$grade <- suppressWarnings(as.integer(g2[["sbr grade"]]))
g2$time_dfs <- suppressWarnings(as.numeric(g2[["dfs time (months)"]]))
g2$event_dfs <- as.integer(g2[["dfs evt"]])

mk <- function(cohort, endpoint, id, time, event, age, stage_high=NULL, grade=NULL) {
  data.frame(cohort=cohort, endpoint=endpoint, sample=id, time=time, event=event, age=age,
             stage_high=if (is.null(stage_high)) NA_integer_ else stage_high,
             grade=if (is.null(grade)) NA_integer_ else grade, stringsAsFactors=FALSE)
}
long <- rbind(
  mk("TCGA","OS", tc$SAMPLE_ID, suppressWarnings(as.numeric(tc$OS_MONTHS)),
     as.integer(tc$OS_STATUS=="1:DECEASED"), suppressWarnings(as.numeric(tc$AGE)), tc$stage_high),
  mk("TCGA","DFS", tc$SAMPLE_ID, suppressWarnings(as.numeric(tc$DFS_MONTHS)),
     as.integer(tc$DFS_STATUS=="1:Recurred/Progressed"), suppressWarnings(as.numeric(tc$AGE)), tc$stage_high),
  mk("METABRIC","OS", mb$SAMPLE_ID, suppressWarnings(as.numeric(mb$OS_MONTHS)),
     as.integer(mb$OS_STATUS=="1:DECEASED"), suppressWarnings(as.numeric(mb$AGE_AT_DIAGNOSIS)), mb$stage_high),
  mk("METABRIC","RFS", mb$SAMPLE_ID, suppressWarnings(as.numeric(mb$RFS_MONTHS)),
     as.integer(mb$RFS_STATUS=="1:Recurred"), suppressWarnings(as.numeric(mb$AGE_AT_DIAGNOSIS)), mb$stage_high),
  mk("GSE58812","OS", g1$sample, g1$time_os, g1$event_os, g1$age),
  mk("GSE21653","DFS", g2$sample, g2$time_dfs, g2$event_dfs, g2$age, grade=g2$grade)
)
long <- long[!is.na(long$sample), ]
write.csv(long, file.path(OUT,"tables/survival_inputs_long.csv"), row.names=FALSE)
cat("\n--- cohort x endpoint n/events ---\n")
for (co in unique(long$cohort)) for (ep in unique(long$endpoint[long$cohort==co])) {
  d <- long[long$cohort==co & long$endpoint==ep, ]
  ok <- !is.na(d$time) & !is.na(d$event) & d$time > 0
  cat(sprintf("  %-9s %-3s n=%3d events=%3d\n", co, ep, sum(ok), sum(d$event[ok]==1)))
}

base$tcga_analysis <- tc; base$mb_analysis <- mb
base$g1_analysis <- g1; base$g2_analysis <- g2
base$tcga_tnbc_def <- "strict: exclude IHC Positive or FISH Positive; include IHC Negative; include IHC Equivocal/Indeterminate with FISH Negative; else exclude"
saveRDS(base, file.path(OUT,"data/checkpoints/analysis_ready.rds"))

flow <- data.frame(
  cohort = c("TCGA","METABRIC","GSE58812","GSE21653"),
  raw_samples = c(1108, 2509, 107, 266),
  primary = c(1097, 2509, 107, 266),
  with_expression = c(1093, 1980, 107, 266),
  triple_negative = c(142, 320, 107, 87),
  tnbc_with_expression = c(141, 320, 107, 87),
  endpoint_analysable = c("OS 141 (21 events); DFS 133 (21)","OS 319 (168); RFS 319 (131)","OS 107 (29)","DFS 83 (25)"),
  stringsAsFactors=FALSE)
write.csv(flow, file.path(OUT,"tables/sample_flow.csv"), row.names=FALSE)
print(flow, row.names=FALSE)
cat("\n[Save] analysis_ready.rds ; survival_inputs_long.csv ; sample_flow.csv\n")
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 修复 9（补充）：重建 TCGA 表达子集（141 例）；含缺失基因整行剔除（不插补）；打分基因集按 HGNC 别名映射到矩阵实际行名
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(data.table); library(GSVA)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
base <- readRDS(file.path(OUT,"data/checkpoints/analysis_ready.rds"))
E <- readRDS(file.path(OUT,"data/checkpoints/expr_tnbc_log2.rds"))

d <- fread(file.path(ROOT,"brca_tcga/data_mrna_seq_v2_rsem.txt"), sep="\t", header=TRUE,
           showProgress=FALSE, na.strings=c("NA","","NaN"), check.names=FALSE, data.table=TRUE)
sym <- d[[1]]; cols <- intersect(base$tcga_analysis$SAMPLE_ID, names(d))
m <- as.matrix(d[, ..cols]); storage.mode(m) <- "double"; rownames(m) <- sym
rm(d); gc()
m <- log2(m + 1)
if (anyDuplicated(rownames(m))) {
  mu <- rowMeans(m, na.rm=TRUE); m <- m[order(rownames(m), -mu), ]; m <- m[!duplicated(rownames(m)), ]
}
E$TCGA <- m
cat(sprintf("[TCGA] %d genes x %d samples ; NA cells=%d\n", nrow(m), ncol(m), sum(is.na(m))))

na_before <- sapply(E, function(x) sum(is.na(x)))
for (co in names(E)) {
  x <- E[[co]]
  x <- x[complete.cases(x), , drop=FALSE]
  sds <- apply(x, 1, sd)
  E[[co]] <- x[is.finite(sds) & sds > 0, , drop=FALSE]
}
cat("[NA cells before gene-level deletion]\n"); print(na_before)
cat("[dims after]\n"); print(sapply(E, dim))

gs <- base$gs
alias <- list(SEPP1=c("SEPP1","SELENOP"), ATP5E=c("ATP5E","ATP5F1E"), ATP5F1E=c("ATP5F1E","ATP5E"),
              TCEB2=c("TCEB2","ELOB"), ELOB=c("ELOB","TCEB2"), CTGF=c("CTGF","CCN2"),
              CXCL8=c("CXCL8","IL8"), FYB=c("FYB","FYB1"), RACK1=c("RACK1","GNB2L1"))
alias_of <- function(g) if (!is.null(alias[[g]])) alias[[g]] else g
map_set <- function(genes, rn) {
  mapped <- vapply(genes, function(g) { hit <- alias_of(g)[alias_of(g) %in% rn]
                                        if (length(hit)) hit[1] else NA_character_ }, character(1))
  unique(mapped[!is.na(mapped)])
}
run_ssgsea <- function(m, sets) {
  sets <- lapply(sets, map_set, rn=rownames(m)); sets <- sets[lengths(sets) >= 5]
  prm <- ssgseaParam(m, sets, normalize=TRUE, alpha=0.25)
  list(scores=gsva(prm, verbose=FALSE), sizes=lengths(sets))
}
run_zmean <- function(m, sets) {
  out <- list(); sizes <- integer()
  for (nm in names(sets)) {
    s <- map_set(sets[[nm]], rownames(m)); sizes[nm] <- length(s)
    if (length(s) < 5) next
    out[[nm]] <- colMeans(t(scale(t(m[s, , drop=FALSE]))), na.rm=TRUE)
  }
  list(scores=do.call(rbind, out), sizes=sizes)
}
res <- list(); cov_rows <- list(); k <- 0
for (co in names(E)) {
  ss <- run_ssgsea(E[[co]], gs); zz <- run_zmean(E[[co]], gs)
  res[[co]] <- list(ssgsea=ss$scores, zmean=zz$scores)
  k <- k + 1
  cov_rows[[k]] <- data.frame(cohort=co, set=names(gs), n_genes_ssgsea=as.integer(ss$sizes[names(gs)]),
                              n_genes_zmean=as.integer(zz$sizes[names(gs)]), stringsAsFactors=FALSE)
  cat(sprintf("[%s] ssGSEA genes %s | zmean %s\n", co, paste(ss$sizes, collapse="/"), paste(zz$sizes, collapse="/")))
}
score_long <- do.call(rbind, lapply(names(res), function(co) {
  do.call(rbind, lapply(c("ssgsea","zmean"), function(meth) {
    s <- res[[co]][[meth]]
    data.frame(cohort=co, method=meth, sample=colnames(s),
               C1Q_TAM=as.numeric(s["C1Q_TAM",]), SPP1_TAM=as.numeric(s["SPP1_TAM",]),
               MMP11_CAF=as.numeric(s["MMP11_CAF",]), stringsAsFactors=FALSE)
  }))
}))
cov_tab <- do.call(rbind, cov_rows)
write.csv(cov_tab, file.path(OUT,"tables/signature_gene_coverage.csv"), row.names=FALSE)
write.csv(score_long, file.path(OUT,"tables/signature_scores.csv"), row.names=FALSE)
saveRDS(list(scores=res, long=score_long, coverage=cov_tab), file.path(OUT,"data/checkpoints/scores.rds"))
saveRDS(E, file.path(OUT,"data/checkpoints/expr_tnbc_log2.rds"))
cat(sprintf("\n[Save] scores + coverage + expr ; score rows=%d\n", nrow(score_long)))
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 修复后重跑主分析：两 cutoff × 两种打分 × 四分析组 × 6 个队列-终点组合
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(survival); library(maxstat)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
sc   <- read.csv(file.path(OUT,"tables/signature_scores.csv"), stringsAsFactors=FALSE)

make_groups <- function(dat, score_col, method, minprop=0.25) {
  out <- list(cut=NA_real_, p_select=NA_real_, ok=FALSE, method=method, err="")
  if (method == "maxstat") {
    dd <- data.frame(time=dat$time, event=dat$event, x=dat[[score_col]])
    dd <- dd[stats::complete.cases(dd), ]
    fit <- tryCatch(maxstat.test(Surv(time, event) ~ x, data=dd, smethod="LogRank",
                                 pmethod="Lau92", minprop=minprop),
                    error=function(e) { out$err <<- conditionMessage(e); NULL })
    if (is.null(fit)) return(out)
    out$cut <- as.numeric(fit$estimate); out$p_select <- as.numeric(fit$p.value); out$ok <- TRUE
  } else {
    out$cut <- median(dat[[score_col]], na.rm=TRUE); out$ok <- TRUE
  }
  out
}
covars_of <- function(co) switch(co, TCGA=c("age","stage_high"), METABRIC=c("age","stage_high"),
                                 GSE58812=c("age"), GSE21653=c("age","grade"))
fit_stats <- function(dat, covars) {
  dat <- dat[!is.na(dat$time) & !is.na(dat$event) & dat$time > 0 & !is.na(dat$grp), ]
  n <- nrow(dat); ev <- sum(dat$event==1)
  r <- list(n=n, events=ev, n_high=sum(dat$grp=="High"), n_low=sum(dat$grp=="Low"),
            ev_high=sum(dat$event[dat$grp=="High"]), ev_low=sum(dat$event[dat$grp=="Low"]),
            hr=NA_real_, lo=NA_real_, hi=NA_real_, p_wald=NA_real_, logrank_p=NA_real_,
            hr_adj=NA_real_, lo_adj=NA_real_, hi_adj=NA_real_, p_adj=NA_real_,
            n_adj=NA_integer_, n_excluded_cov=NA_integer_, ph_p=NA_real_, note="")
  if (n < 20 || ev < 5) { r$note <- "insufficient (n<20 or events<5)"; return(r) }
  if (min(r$n_high, r$n_low) < 10) { r$note <- "group n<10"; return(r) }
  sd <- survdiff(Surv(time,event) ~ grp, data=dat)
  r$logrank_p <- 1 - pchisq(sd$chisq, df=length(sd$n)-1)
  fit <- coxph(Surv(time,event) ~ grp, data=dat)
  r$hr <- exp(coef(fit)[["grpHigh"]]); ci <- exp(confint(fit))
  r$lo <- ci["grpHigh",1]; r$hi <- ci["grpHigh",2]
  r$p_wald <- summary(fit)$coefficients["grpHigh","Pr(>|z|)"]
  cv <- covars[covars %in% names(dat)]
  if (length(cv) > 0) {
    dat2 <- dat[complete.cases(dat[, c("time","event","grp",cv)]), ]
    r$n_excluded_cov <- n - nrow(dat2); r$n_adj <- nrow(dat2)
    if (nrow(dat2) >= 20 && sum(dat2$event) >= 5 && min(table(dat2$grp)) >= 10) {
      f2 <- tryCatch(coxph(as.formula(paste("Surv(time,event) ~ grp +", paste(cv, collapse="+"))), data=dat2),
                     error=function(e) NULL)
      if (!is.null(f2)) {
        r$hr_adj <- exp(coef(f2)[["grpHigh"]]); ci2 <- exp(confint(f2))
        r$lo_adj <- ci2["grpHigh",1]; r$hi_adj <- ci2["grpHigh",2]
        r$p_adj <- summary(f2)$coefficients["grpHigh","Pr(>|z|)"]
        zph <- cox.zph(f2)$table
        r$ph_p <- if ("grp" %in% rownames(zph)) as.numeric(zph["grp","p"]) else NA_real_
      }
    } else r$note <- paste(r$note, "adjusted model skipped (n/events)", sep="; ")
  }
  r
}

combos <- data.frame(cohort=c("TCGA","TCGA","METABRIC","METABRIC","GSE58812","GSE21653"),
                     endpoint=c("OS","DFS","OS","RFS","OS","DFS"), stringsAsFactors=FALSE)
res_primary <- list(); res_dual2 <- list(); cuts_tab <- list(); grp_tab <- list()
k <- 0; kd <- 0; kc <- 0; kg <- 0
for (i in seq_len(nrow(combos))) {
  co <- combos$cohort[i]; ep <- combos$endpoint[i]
  d0 <- long[long$cohort==co & long$endpoint==ep, ]
  d0 <- d0[!is.na(d0$time) & !is.na(d0$event) & d0$time > 0, ]
  cvs <- covars_of(co)
  for (meth in c("ssgsea","zmean")) {
    s <- sc[sc$cohort==co & sc$method==meth, ]
    d <- merge(d0, s[, c("sample","C1Q_TAM","SPP1_TAM","MMP11_CAF")], by="sample")
    if (nrow(d) < 20) next
    for (cm in c("maxstat","median")) {
      cuts <- lapply(c("C1Q_TAM","SPP1_TAM","MMP11_CAF"), function(g) make_groups(d, g, cm))
      names(cuts) <- c("C1Q_TAM","SPP1_TAM","MMP11_CAF")
      kc <- kc + 1
      cuts_tab[[kc]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm,
                                   signature=names(cuts), cutpoint=sapply(cuts, function(z) z$cut),
                                   selection_p=sapply(cuts, function(z) z$p_select),
                                   n=nrow(d), events=sum(d$event), stringsAsFactors=FALSE)
      g1 <- ifelse(d$C1Q_TAM > cuts$C1Q_TAM$cut, "High","Low")
      g2 <- ifelse(d$SPP1_TAM > cuts$SPP1_TAM$cut, "High","Low")
      g3 <- ifelse(d$MMP11_CAF > cuts$MMP11_CAF$cut, "High","Low")
      d$C1Q_grp <- factor(g1, levels=c("Low","High")); d$SPP1_grp <- factor(g2, levels=c("Low","High"))
      d$MMP11_grp <- factor(g3, levels=c("Low","High"))
      d$dual_grp <- ifelse(g1=="High" & g3=="High", "DualHigh",
                           ifelse(g1=="Low" & g3=="Low", "DualLow", "Other"))
      kg <- kg + 1
      grp_tab[[kg]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm,
                                  sample=d$sample, C1Q_grp=d$C1Q_grp, SPP1_grp=d$SPP1_grp,
                                  MMP11_grp=d$MMP11_grp, dual_grp=d$dual_grp, stringsAsFactors=FALSE)
      specs <- list(
        list(label="C1Q_TAM", dat=transform(d, grp=C1Q_grp), p_sel=cuts$C1Q_TAM$p_select),
        list(label="SPP1_TAM", dat=transform(d, grp=SPP1_grp), p_sel=cuts$SPP1_TAM$p_select),
        list(label="MMP11_CAF", dat=transform(d, grp=MMP11_grp), p_sel=cuts$MMP11_CAF$p_select),
        list(label="C1Q_MMP11_DualHigh_vs_Other",
             dat=transform(d, grp=factor(ifelse(dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))),
             p_sel=NA_real_))
      for (sp in specs) {
        st <- fit_stats(sp$dat, cvs); k <- k + 1
        adj_used <- cm=="maxstat" && !is.na(sp$p_sel)
        pe <- if (adj_used) sp$p_sel else st$logrank_p
        res_primary[[k]] <- data.frame(
          cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm, analysis_group=sp$label,
          n=st$n, events=st$events, n_high=st$n_high, n_low=st$n_low,
          events_high=st$ev_high, events_low=st$ev_low,
          cutpoint_C1Q=cuts$C1Q_TAM$cut, cutpoint_SPP1=cuts$SPP1_TAM$cut, cutpoint_MMP11=cuts$MMP11_CAF$cut,
          HR=st$hr, HR_lower=st$lo, HR_upper=st$hi, p_wald=st$p_wald, logrank_p=st$logrank_p,
          p_effective=pe,
          p_source=if (adj_used) "maxstat Lau92 (selection-adjusted)" else "log-rank (unadjusted)",
          HR_adjusted=st$hr_adj, HR_adj_lower=st$lo_adj, HR_adj_upper=st$hi_adj, p_adjusted=st$p_adj,
          n_adjusted=st$n_adj, n_excluded_missing_covariate=st$n_excluded_cov,
          covariates=paste(cvs, collapse="+"), ph_test_p=st$ph_p, note=st$note, stringsAsFactors=FALSE)
      }
      dsub <- d[d$dual_grp %in% c("DualHigh","DualLow"), ]
      dsub$grp <- factor(ifelse(dsub$dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))
      st2 <- fit_stats(dsub, cvs); kd <- kd + 1
      res_dual2[[kd]] <- data.frame(
        cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm,
        analysis_group="C1Q_MMP11_DualHigh_vs_DualLow",
        n=st2$n, events=st2$events, n_high=st2$n_high, n_low=st2$n_low,
        events_high=st2$ev_high, events_low=st2$ev_low,
        HR=st2$hr, HR_lower=st2$lo, HR_upper=st2$hi, p_wald=st2$p_wald, logrank_p=st2$logrank_p,
        HR_adjusted=st2$hr_adj, HR_adj_lower=st2$lo_adj, HR_adj_upper=st2$hi_adj, p_adjusted=st2$p_adj,
        n_adjusted=st2$n_adj, ph_test_p=st2$ph_p, note=st2$note, stringsAsFactors=FALSE)
    }
  }
}
prim <- do.call(rbind, res_primary); dual2 <- do.call(rbind, res_dual2)
cuts_df <- do.call(rbind, cuts_tab); grp_df <- do.call(rbind, grp_tab)
prim$p_fdr <- NA_real_
fam <- paste(prim$score_method, prim$cutoff_method, prim$cohort, prim$endpoint, sep="|")
for (f in unique(fam)) {
  idx <- which(fam == f & !is.na(prim$p_effective))
  if (length(idx) > 0) prim$p_fdr[idx] <- p.adjust(prim$p_effective[idx], method="BH")
}
prim$significant_fdr <- !is.na(prim$p_fdr) & prim$p_fdr < 0.05
write.csv(prim, file.path(OUT,"tables/survival_results_primary.csv"), row.names=FALSE)
write.csv(dual2, file.path(OUT,"tables/survival_results_dual_secondary.csv"), row.names=FALSE)
write.csv(cuts_df, file.path(OUT,"tables/stratification_cutpoints.csv"), row.names=FALSE)
write.csv(grp_df, file.path(OUT,"tables/group_assignment_long.csv"), row.names=FALSE)
saveRDS(list(primary=prim, dual_secondary=dual2, cuts=cuts_df, groups=grp_df),
        file.path(OUT,"data/checkpoints/survival_results.rds"))
cat(sprintf("[Rows] primary=%d ; dual=%d ; not fitted=%d ; FDR<0.05=%d\n",
            nrow(prim), nrow(dual2), sum(is.na(prim$HR)), sum(prim$p_fdr < 0.05, na.rm=TRUE)))
cat("\n--- FDR<0.05 (both score methods) ---\n")
s <- prim[!is.na(prim$p_fdr) & prim$p_fdr < 0.05,
          c("cohort","endpoint","score_method","cutoff_method","analysis_group","n_high","n_low","HR","HR_lower","HR_upper","p_effective","p_fdr","HR_adjusted")]
s[, c("HR","HR_lower","HR_upper","HR_adjusted")] <- round(s[, c("HR","HR_lower","HR_upper","HR_adjusted")],3)
s[, c("p_effective","p_fdr")] <- signif(s[, c("p_effective","p_fdr")],3)
print(s, row.names=FALSE)
cat("\n--- TCGA 全部结果（ssGSEA）---\n")
t <- prim[prim$cohort=="TCGA" & prim$score_method=="ssgsea",
          c("endpoint","cutoff_method","analysis_group","n_high","n_low","events_high","events_low","HR","HR_lower","HR_upper","logrank_p","p_effective","p_fdr")]
t[, c("HR","HR_lower","HR_upper")] <- round(t[, c("HR","HR_lower","HR_upper")],3)
t[, c("logrank_p","p_effective","p_fdr")] <- signif(t[, c("logrank_p","p_effective","p_fdr")],3)
print(t, row.names=FALSE)
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 新增产出：中位生存汇总、逐基因覆盖矩阵、新旧 signature 相关、maxstat minprop 敏感性
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(survival); library(maxstat)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
grp  <- read.csv(file.path(OUT,"tables/group_assignment_long.csv"), stringsAsFactors=FALSE)
sc   <- read.csv(file.path(OUT,"tables/signature_scores.csv"), stringsAsFactors=FALSE)
base <- readRDS(file.path(OUT,"data/checkpoints/analysis_ready.rds"))
E    <- readRDS(file.path(OUT,"data/checkpoints/expr_tnbc_log2.rds"))
gs   <- base$gs
alias <- list(SEPP1=c("SEPP1","SELENOP"), ATP5E=c("ATP5E","ATP5F1E"), ATP5F1E=c("ATP5F1E","ATP5E"),
              TCEB2=c("TCEB2","ELOB"), ELOB=c("ELOB","TCEB2"), CTGF=c("CTGF","CCN2"),
              CXCL8=c("CXCL8","IL8"), FYB=c("FYB","FYB1"), RACK1=c("RACK1","GNB2L1"))
alias_of <- function(g) if (!is.null(alias[[g]])) alias[[g]] else g

# ---- 2) 逐基因覆盖矩阵（60 条 × 4 队列）----
cov <- list(); k <- 0
for (sn in names(gs)) for (g in gs[[sn]]) for (co in names(E)) {
  k <- k + 1
  cov[[k]] <- data.frame(signature=sn, gene=g, cohort=co,
                         present=as.integer(any(alias_of(g) %in% rownames(E[[co]]))),
                         stringsAsFactors=FALSE)
}
covmat <- do.call(rbind, cov)
write.csv(covmat, file.path(OUT,"tables/signature_gene_coverage_matrix.csv"), row.names=FALSE)
cat("[coverage matrix] rows:", nrow(covmat), "| present:", sum(covmat$present), "\n")
for (co in names(E)) for (sn in names(gs)) {
  miss <- covmat$gene[covmat$cohort==co & covmat$signature==sn & covmat$present==0]
  cat(sprintf("  %-9s %-10s missing(%d): %s\n", co, sn, length(miss),
              if (length(miss)) paste(miss, collapse=",") else "none"))
}

# ---- 4) 中位生存汇总 ----
km_rows <- list(); k <- 0
for (co in unique(grp$cohort)) for (ep in unique(grp$endpoint[grp$cohort==co])) for (meth in c("ssgsea","zmean")) for (cm in c("maxstat","median")) {
  d0 <- long[long$cohort==co & long$endpoint==ep, ]
  g0 <- grp[grp$cohort==co & grp$endpoint==ep & grp$score_method==meth & grp$cutoff_method==cm, ]
  if (nrow(g0)==0) next
  d <- merge(d0[, c("sample","time","event")], g0, by="sample")
  d <- d[!is.na(d$time) & !is.na(d$event) & d$time > 0, ]
  med <- function(fit) { t <- summary(fit)$table; data.frame(n=t[,"records"], events=t[,"events"],
                        median=t[,"median"], lower=t[,"0.95LCL"], upper=t[,"0.95UCL"]) }
  for (gname in c("C1Q_grp","SPP1_grp","MMP11_grp")) {
    dd <- d; dd$grp <- factor(dd[[gname]], levels=c("Low","High"))
    fit <- survfit(Surv(time,event) ~ grp, data=dd); mt <- med(fit)
    p <- 1 - pchisq(survdiff(Surv(time,event) ~ grp, data=dd)$chisq, 1)
    for (lv in rownames(mt)) {
      k <- k + 1
      km_rows[[k]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm,
                                 group=sub("grp=", "", lv), n=mt[lv,"n"], events=mt[lv,"events"],
                                 median_months=mt[lv,"median"], lower95=mt[lv,"lower"], upper95=mt[lv,"upper"],
                                 logrank_p_primary=p, logrank_p_secondary=NA_real_, stringsAsFactors=FALSE)
    }
  }
  dd <- d; dd$grp <- factor(dd$dual_grp, levels=c("DualHigh","Other","DualLow"))
  fit <- survfit(Surv(time,event) ~ grp, data=dd); mt <- med(fit)
  p1 <- 1 - pchisq(survdiff(Surv(time,event) ~ grp, data=dd[dd$dual_grp %in% c("DualHigh","Other"), ])$chisq, 1)
  p2 <- 1 - pchisq(survdiff(Surv(time,event) ~ grp, data=dd[dd$dual_grp %in% c("DualHigh","DualLow"), ])$chisq, 1)
  for (lv in rownames(mt)) {
    k <- k + 1
    km_rows[[k]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm,
                               group=sub("grp=", "", lv), n=mt[lv,"n"], events=mt[lv,"events"],
                               median_months=mt[lv,"median"], lower95=mt[lv,"lower"], upper95=mt[lv,"upper"],
                               logrank_p_primary=p1, logrank_p_secondary=p2, stringsAsFactors=FALSE)
  }
}
km_sum <- do.call(rbind, km_rows)
write.csv(km_sum, file.path(OUT,"tables/km_summary.csv"), row.names=FALSE)
cat(sprintf("\n[km_summary] rows=%d\n", nrow(km_sum)))

# ---- 3) 新旧 signature 相关（旧表仅 TCGA/METABRIC，取 *_top20_specific）----
old <- read.csv(file.path(ROOT,"processed/scores.csv"), stringsAsFactors=FALSE)
old_sp <- old[grepl("_specific$", old$signature), ]
map_sig <- function(x) ifelse(grepl("^C1Q", x), "C1Q_TAM", ifelse(grepl("^SPP1", x), "SPP1_TAM",
                              ifelse(grepl("^MMP11", x), "MMP11_CAF", NA)))
cor_rows <- list(); k <- 0
for (co in unique(old_sp$cohort)) for (osig in unique(old_sp$signature[old_sp$cohort==co])) {
  nsig <- map_sig(osig); if (is.na(nsig)) next
  o <- old_sp[old_sp$cohort==co & old_sp$signature==osig, c("sample","ssgsea","zscore_mean")]
  a <- sc[sc$cohort==co & sc$method=="ssgsea", c("sample", nsig)]; names(a)[2] <- "new_ss"
  b <- sc[sc$cohort==co & sc$method=="zmean",  c("sample", nsig)]; names(b)[2] <- "new_z"
  m <- merge(merge(o, a, by="sample"), b, by="sample")
  k <- k + 1
  cor_rows[[k]] <- data.frame(
    cohort=co, new_signature=nsig, old_signature=osig, n_matched=nrow(m),
    rho_new_ssgsea_vs_old_ssgsea=if (nrow(m)>3) cor(m$new_ss, m$ssgsea, method="spearman") else NA,
    rho_new_ssgsea_vs_old_zmean =if (nrow(m)>3) cor(m$new_ss, m$zscore_mean, method="spearman") else NA,
    rho_new_zmean_vs_old_zmean  =if (nrow(m)>3) cor(m$new_z, m$zscore_mean, method="spearman") else NA,
    stringsAsFactors=FALSE)
}
cor_tab <- do.call(rbind, cor_rows)
write.csv(cor_tab, file.path(OUT,"tables/direction_check_vs_old_signature.csv"), row.names=FALSE)
cat("\n[old-signature correlation]\n")
print(transform(cor_tab, rho_new_ssgsea_vs_old_ssgsea=round(rho_new_ssgsea_vs_old_ssgsea,3),
                rho_new_ssgsea_vs_old_zmean=round(rho_new_ssgsea_vs_old_zmean,3),
                rho_new_zmean_vs_old_zmean=round(rho_new_zmean_vs_old_zmean,3)), row.names=FALSE)

# ---- 7) maxstat minprop 0.10 vs 0.25（三个关键 C1Q 结果）----
key <- data.frame(cohort=c("METABRIC","METABRIC","GSE21653"), endpoint=c("OS","RFS","DFS"), stringsAsFactors=FALSE)
cmp <- list(); k <- 0
for (i in seq_len(nrow(key))) {
  co <- key$cohort[i]; ep <- key$endpoint[i]
  d0 <- long[long$cohort==co & long$endpoint==ep, ]
  s <- sc[sc$cohort==co & sc$method=="ssgsea", c("sample","C1Q_TAM")]
  d <- merge(d0[, c("sample","time","event")], s, by="sample")
  d <- d[!is.na(d$time) & !is.na(d$event) & d$time > 0, ]
  for (mp in c(0.10, 0.25)) {
    fit <- maxstat.test(Surv(time,event) ~ C1Q_TAM, data=d, smethod="LogRank", pmethod="Lau92", minprop=mp)
    cut <- as.numeric(fit$estimate)
    d$grp <- factor(ifelse(d$C1Q_TAM > cut, "High","Low"), levels=c("Low","High"))
    cx <- coxph(Surv(time,event) ~ grp, data=d); ci <- exp(confint(cx))
    k <- k + 1
    cmp[[k]] <- data.frame(cohort=co, endpoint=ep, minprop=mp, cutpoint=cut,
                           n_high=sum(d$grp=="High"), n_low=sum(d$grp=="Low"),
                           HR=exp(coef(cx)[["grpHigh"]]), HR_lower=ci["grpHigh",1], HR_upper=ci["grpHigh",2],
                           p_adjusted=fit$p.value, stringsAsFactors=FALSE)
  }
}
cmp <- do.call(rbind, cmp)
write.csv(cmp, file.path(OUT,"tables/maxstat_minprop_sensitivity.csv"), row.names=FALSE)
cat("\n[maxstat minprop sensitivity]\n")
print(transform(cmp, cutpoint=round(cutpoint,3), HR=round(HR,3), HR_lower=round(HR_lower,3),
                HR_upper=round(HR_upper,3), p_adjusted=signif(p_adjusted,3)), row.names=FALSE)
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 重算多因素 Cox 全表、PH 检验、方向校验与两种一致性表（修复后版本）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages(library(survival))
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
grp  <- read.csv(file.path(OUT,"tables/group_assignment_long.csv"), stringsAsFactors=FALSE)
prim <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)
sc   <- read.csv(file.path(OUT,"tables/signature_scores.csv"), stringsAsFactors=FALSE)
E    <- readRDS(file.path(OUT,"data/checkpoints/expr_tnbc_log2.rds"))
covars_of <- function(co) switch(co, TCGA=c("age","stage_high"), METABRIC=c("age","stage_high"),
                                 GSE58812=c("age"), GSE21653=c("age","grade"))

# 多因素 Cox 全表
rows <- list(); k <- 0
for (i in seq_len(nrow(prim))) {
  r <- prim[i, ]
  if (r$score_method != "ssgsea" || is.na(r$HR_adjusted)) next
  d0 <- long[long$cohort==r$cohort & long$endpoint==r$endpoint, ]
  g0 <- grp[grp$cohort==r$cohort & grp$endpoint==r$endpoint & grp$score_method=="ssgsea" & grp$cutoff_method==r$cutoff_method, ]
  d <- merge(d0[, c("sample","time","event","age","stage_high","grade")], g0, by="sample")
  if (r$analysis_group == "C1Q_MMP11_DualHigh_vs_Other") {
    d <- d[d$dual_grp != "Other", ]; d$grp <- factor(ifelse(d$dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))
  } else {
    col <- switch(r$analysis_group, C1Q_TAM="C1Q_grp", SPP1_TAM="SPP1_grp", MMP11_CAF="MMP11_grp")
    d$grp <- factor(d[[col]], levels=c("Low","High"))
  }
  cv <- covars_of(r$cohort); cv <- cv[cv %in% names(d)]
  d <- d[complete.cases(d[, c("time","event","grp",cv)]), ]
  fit <- tryCatch(coxph(as.formula(paste("Surv(time,event) ~ grp +", paste(cv, collapse="+"))), data=d), error=function(e) NULL)
  if (is.null(fit)) next
  s <- summary(fit)
  for (term in rownames(s$coefficients)) {
    k <- k + 1
    rows[[k]] <- data.frame(cohort=r$cohort, endpoint=r$endpoint, cutoff_method=r$cutoff_method,
      analysis_group=r$analysis_group, term=term, n=nrow(d), events=sum(d$event),
      HR=s$conf.int[term,"exp(coef)"], HR_lower=s$conf.int[term,"lower .95"],
      HR_upper=s$conf.int[term,"upper .95"], p=s$coefficients[term,"Pr(>|z|)"], stringsAsFactors=FALSE)
  }
}
coxfull <- do.call(rbind, rows)
coxfull$q_within_model <- ave(coxfull$p, paste(coxfull$cohort, coxfull$endpoint, coxfull$cutoff_method,
                                               coxfull$analysis_group), FUN=function(x) p.adjust(x, "BH"))
write.csv(coxfull, file.path(OUT,"tables/cox_multivariable_full.csv"), row.names=FALSE)
cat(sprintf("[cox full] rows=%d ; models=%d\n", nrow(coxfull),
            length(unique(paste(coxfull$cohort,coxfull$endpoint,coxfull$cutoff_method,coxfull$analysis_group)))))

# PH 表
ph <- prim[, c("cohort","endpoint","score_method","cutoff_method","analysis_group","n_adjusted","ph_test_p")]
ph$ph_violation <- !is.na(ph$ph_test_p) & ph$ph_test_p < 0.05
write.csv(ph, file.path(OUT,"tables/ph_assumption_tests.csv"), row.names=FALSE)
cat(sprintf("[PH] models=%d ; violations=%d\n", sum(!is.na(ph$ph_test_p)), sum(ph$ph_violation, na.rm=TRUE)))
if (any(ph$ph_violation, na.rm=TRUE)) print(ph[which(ph$ph_violation), ], row.names=FALSE)

# 方向校验
anchor <- list(C1Q_TAM="C1QA", SPP1_TAM="SPP1", MMP11_CAF="MMP11")
dir_rows <- list(); k <- 0
for (co in names(E)) for (scn in names(anchor)) {
  g <- anchor[[scn]]
  if (!(g %in% rownames(E[[co]]))) next
  s <- sc[sc$cohort==co & sc$method=="ssgsea", ]
  s <- s[match(colnames(E[[co]]), s$sample), ]
  rho <- suppressWarnings(cor(s[[scn]], as.numeric(E[[co]][g, ]), method="spearman", use="complete.obs"))
  k <- k + 1
  dir_rows[[k]] <- data.frame(cohort=co, signature=scn, anchor_gene=g, spearman_rho=rho,
                              direction_ok=!is.na(rho) && rho > 0, stringsAsFactors=FALSE)
}
dirchk <- do.call(rbind, dir_rows)
write.csv(dirchk, file.path(OUT,"tables/direction_check.csv"), row.names=FALSE)
cat("\n[direction check]\n"); print(transform(dirchk, spearman_rho=round(spearman_rho,3)), row.names=FALSE)

# 两种 cutoff 一致性
cons <- list(); k <- 0
for (co in unique(prim$cohort)) for (ep in unique(prim$endpoint[prim$cohort==co])) for (meth in c("ssgsea","zmean")) for (ag in unique(prim$analysis_group)) {
  a <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method==meth & prim$cutoff_method=="maxstat" & prim$analysis_group==ag, ]
  b <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method==meth & prim$cutoff_method=="median" & prim$analysis_group==ag, ]
  if (nrow(a)==0 || nrow(b)==0) next
  k <- k + 1
  cons[[k]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, analysis_group=ag,
    HR_maxstat=a$HR, p_maxstat=a$p_effective, q_maxstat=a$p_fdr,
    HR_median=b$HR, p_median=b$logrank_p, q_median=b$p_fdr,
    direction_concordant=ifelse(is.na(a$HR) || is.na(b$HR), NA, sign(log(a$HR))==sign(log(b$HR))),
    sig_maxstat=!is.na(a$p_fdr) & a$p_fdr < 0.05, sig_median=!is.na(b$p_fdr) & b$p_fdr < 0.05,
    stringsAsFactors=FALSE)
}
cons <- do.call(rbind, cons)
cons$significance_concordant <- cons$sig_maxstat == cons$sig_median
write.csv(cons, file.path(OUT,"tables/cutoff_consistency.csv"), row.names=FALSE)
cat(sprintf("\n[cutoff consistency] %d combos ; direction concordant %d ; significance concordant %d\n",
            nrow(cons), sum(cons$direction_concordant, na.rm=TRUE), sum(cons$significance_concordant, na.rm=TRUE)))

# 两种打分一致性
sens <- list(); k <- 0
for (co in unique(prim$cohort)) for (ep in unique(prim$endpoint[prim$cohort==co])) for (cm in c("maxstat","median")) for (ag in unique(prim$analysis_group)) {
  a <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method=="ssgsea" & prim$cutoff_method==cm & prim$analysis_group==ag, ]
  b <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method=="zmean" & prim$cutoff_method==cm & prim$analysis_group==ag, ]
  if (nrow(a)==0 || nrow(b)==0) next
  k <- k + 1
  sens[[k]] <- data.frame(cohort=co, endpoint=ep, cutoff_method=cm, analysis_group=ag,
    HR_ssgsea=a$HR, q_ssgsea=a$p_fdr, HR_zmean=b$HR, q_zmean=b$p_fdr,
    direction_concordant=ifelse(is.na(a$HR)||is.na(b$HR), NA, sign(log(a$HR))==sign(log(b$HR))), stringsAsFactors=FALSE)
}
sens <- do.call(rbind, sens)
write.csv(sens, file.path(OUT,"tables/score_method_consistency.csv"), row.names=FALSE)
cat(sprintf("[score consistency] %d combos ; direction concordant %d ; both significant %d\n",
            nrow(sens), sum(sens$direction_concordant, na.rm=TRUE),
            sum(!is.na(sens$q_ssgsea) & !is.na(sens$q_zmean) & sens$q_ssgsea<0.05 & sens$q_zmean<0.05)))
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 重绘全部 KM 图（修复后数据；配色改用 ggsci NPG 色板，不在调用点散写 hex）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(survival); library(survminer); library(ggplot2); library(ggsci)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
grp  <- read.csv(file.path(OUT,"tables/group_assignment_long.csv"), stringsAsFactors=FALSE)
prim <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)

npg2 <- unname(ggsci::pal_npg("nrc")(2))
cat("[palette] ggsci pal_npg(2):", paste(npg2, collapse=", "), "\n")
dir.create(file.path(OUT,"figures/km"), recursive=TRUE, showWarnings=FALSE)
save_km <- function(pl, stem) {
  pdf(paste0(stem,".pdf"), width=6, height=5.6, bg="white"); print(pl); dev.off()
  png(paste0(stem,".png"), width=1800, height=1680, res=300, bg="white"); print(pl); dev.off()
}
group_meta <- list(C1Q_grp=list(col="C1Q_grp",lab="C1Q TAMs"),
                   SPP1_grp=list(col="SPP1_grp",lab="SPP1 TAMs"),
                   MMP11_grp=list(col="MMP11_grp",lab="MMP11 CAFs"),
                   dual_grp=list(col="dual_grp",lab="C1Q+ and MMP11+"))
n_done <- 0; skipped <- c()
combos <- unique(grp[, c("cohort","endpoint","score_method","cutoff_method")])
for (i in seq_len(nrow(combos))) {
  co <- combos$cohort[i]; ep <- combos$endpoint[i]
  if (combos$score_method[i] != "ssgsea") next
  cm <- combos$cutoff_method[i]
  d0 <- long[long$cohort==co & long$endpoint==ep, ]
  g0 <- grp[grp$cohort==co & grp$endpoint==ep & grp$score_method=="ssgsea" & grp$cutoff_method==cm, ]
  d <- merge(d0[, c("sample","time","event")], g0, by="sample")
  d <- d[!is.na(d$time) & !is.na(d$event) & d$time > 0, ]
  pr <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method=="ssgsea" & prim$cutoff_method==cm, ]
  for (gname in names(group_meta)) {
    gm <- group_meta[[gname]]; dd <- d
    if (gname == "dual_grp") {
      dd <- dd[dd$dual_grp != "Other", ]
      dd$grp <- factor(ifelse(dd$dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))
    } else dd$grp <- factor(dd[[gm$col]], levels=c("Low","High"))
    dd <- dd[!is.na(dd$grp), ]
    if (nrow(dd) < 10 || sum(dd$event) < 5) { skipped <- c(skipped, paste(co, ep, cm, gname, sep="/")); next }
    lab <- if (gname=="dual_grp") "C1Q+MMP11+ vs others" else gm$lab
    row <- pr[pr$analysis_group == if (gname=="dual_grp") "C1Q_MMP11_DualHigh_vs_Other" else
              switch(gname, C1Q_grp="C1Q_TAM", SPP1_grp="SPP1_TAM", MMP11_grp="MMP11_CAF"), ]
    fit <- survfit(Surv(time,event) ~ grp, data=dd)
    sub <- if (nrow(row)==0 || is.na(row$HR[1])) "model not fitted (small groups)" else
      sprintf("HR=%.2f (95%%CI %.2f-%.2f); log-rank p=%.3g; adj p=%.3g",
              row$HR[1], row$HR_lower[1], row$HR_upper[1], row$logrank_p[1], row$p_effective[1])
    pl <- ggsurvplot(fit, data=dd, palette=npg2, conf.int=FALSE,
                     risk.table=TRUE, risk.table.height=0.24, censor=TRUE,
                     legend.title="", legend.labs=c("Low","High"),
                     xlab="Time (months)", ylab="Survival probability",
                     title=sprintf("%s %s - %s (%s cutoff)", co, ep, lab, cm),
                     subtitle=sub, font.subtitle=c(9,"plain"),
                     ggtheme=theme_bw() + theme(panel.grid=element_blank(),
                                                plot.title=element_text(hjust=0.5, face="bold", size=10),
                                                plot.subtitle=element_text(hjust=0.5),
                                                legend.position="right"))
    save_km(pl, file.path(OUT,"figures/km", sprintf("km_%s_%s_%s_%s", co, ep, cm, gname)))
    n_done <- n_done + 1
  }
}
cat(sprintf("[KM] redrawn %d figures (pdf+png)\n", n_done))
cat("[KM] skipped:", if (length(skipped)) paste(skipped, collapse="; ") else "none", "\n")
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 重绘森林图与 HR 热图（修复后数据；定性色取 ggsci NPG，发散色取 RColorBrewer RdBu：负值=蓝、正值=红）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

suppressPackageStartupMessages({library(ggplot2); library(ggsci); library(RColorBrewer)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
prim <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)
dir.create(file.path(OUT,"figures/forest"), recursive=TRUE, showWarnings=FALSE)
pal4 <- unname(ggsci::pal_npg("nrc")(4))
rdbu <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(11)
cat("[palette] NPG(4):", paste(pal4, collapse=", "), "| RdBu low(blue)=", rdbu[1], " high(red)=", rdbu[11], "\n")

lab_map <- c(C1Q_TAM="C1Q TAMs", SPP1_TAM="SPP1 TAMs", MMP11_CAF="MMP11 CAFs",
             C1Q_MMP11_DualHigh_vs_Other="C1Q+ & MMP11+ (vs other)")
prim$grp_lab <- lab_map[prim$analysis_group]
prim$star <- ifelse(is.na(prim$p_fdr), "", ifelse(prim$p_fdr < 0.001, "***",
                    ifelse(prim$p_fdr < 0.01, "**", ifelse(prim$p_fdr < 0.05, "*", ""))))
forest_panel <- function(df, xlab_txt, title_txt) {
  df <- df[!is.na(df$HR) & !is.na(df$HR_lower), ]
  if (nrow(df)==0) return(NULL)
  df$row <- factor(df$rowlab, levels=rev(unique(df$rowlab)))
  ggplot(df, aes(x=HR, y=row, color=grp_lab)) +
    geom_vline(xintercept=1, linetype="dashed", color="grey40", linewidth=0.4) +
    geom_errorbarh(aes(xmin=HR_lower, xmax=HR_upper), height=0.18, linewidth=0.5) +
    geom_point(shape=15, size=2.4) +
    geom_text(aes(label=paste0(sprintf("%.2f", HR), star)), vjust=-0.9, size=2.9, show.legend=FALSE) +
    scale_x_log10() +
    scale_color_manual(values=setNames(pal4[seq_along(unique(df$grp_lab))], unique(df$grp_lab))) +
    labs(x=xlab_txt, y=NULL, title=title_txt, color=NULL) +
    theme_bw() +
    theme(panel.grid=element_blank(), panel.border=element_rect(colour="black", fill=NA),
          plot.title=element_text(hjust=0.5, face="bold", size=10),
          legend.position="bottom", legend.text=element_text(size=8),
          axis.text.y=element_text(size=8))
}
n1 <- 0; n2 <- 0
for (co in unique(prim$cohort)) for (ep in unique(prim$endpoint[prim$cohort==co])) for (cm in c("maxstat","median")) {
  sub <- prim[prim$cohort==co & prim$endpoint==ep & prim$cutoff_method==cm & prim$score_method=="ssgsea", ]
  if (nrow(sub)==0) next
  sub$rowlab <- paste0(sub$grp_lab, "  [n=", sub$n_high + sub$n_low, "]")
  p1 <- forest_panel(sub, "Hazard ratio (high vs low, log scale)", sprintf("%s %s - %s cutoff", co, ep, cm))
  if (!is.null(p1)) {
    stem <- file.path(OUT,"figures/forest", sprintf("forest_all_%s_%s_%s", co, ep, cm))
    ggsave(paste0(stem,".pdf"), p1, width=5.4, height=2.6, bg="white")
    ggsave(paste0(stem,".png"), p1, width=5.4, height=2.6, dpi=300, bg="white"); n1 <- n1 + 1
  }
  d2 <- rbind(
    data.frame(grp_lab=sub$grp_lab, HR=sub$HR, HR_lower=sub$HR_lower, HR_upper=sub$HR_upper,
               star=sub$star, model="unadjusted", stringsAsFactors=FALSE),
    data.frame(grp_lab=sub$grp_lab, HR=sub$HR_adjusted, HR_lower=sub$HR_adj_lower, HR_upper=sub$HR_adj_upper,
               star=ifelse(is.na(sub$p_adjusted),"", ifelse(sub$p_adjusted<0.001,"***",
                    ifelse(sub$p_adjusted<0.01,"**", ifelse(sub$p_adjusted<0.05,"*","")))),
               model="adjusted", stringsAsFactors=FALSE))
  d2 <- d2[!is.na(d2$HR) & !is.na(d2$HR_lower), ]
  if (nrow(d2) > 0) {
    d2$rowlab <- paste0(d2$grp_lab, " - ", d2$model)
    p2 <- forest_panel(d2, "Hazard ratio (high vs low, log scale)",
                       sprintf("%s %s - %s cutoff (adjusted: %s)", co, ep, cm, paste(unique(sub$covariates), collapse=",")))
    if (!is.null(p2)) {
      stem <- file.path(OUT,"figures/forest", sprintf("forest_adjusted_%s_%s_%s", co, ep, cm))
      ggsave(paste0(stem,".pdf"), p2, width=5.6, height=3.6, bg="white")
      ggsave(paste0(stem,".png"), p2, width=5.6, height=3.6, dpi=300, bg="white"); n2 <- n2 + 1
    }
  }
}
cat(sprintf("[Forest] all=%d ; adjusted=%d\n", n1, n2))

hm <- prim[prim$score_method=="ssgsea" & !is.na(prim$HR), ]
hm$key <- paste(hm$cohort, hm$endpoint, hm$cutoff_method, sep=" | ")
hm$grp_lab <- factor(hm$grp_lab, levels=unname(lab_map))
ph <- ggplot(hm, aes(x=grp_lab, y=key, fill=log2(HR))) +
  geom_tile(color=NA) +
  geom_text(aes(label=ifelse(is.na(star) | star=="", sprintf("%.2f",HR), paste0(sprintf("%.2f",HR), star))), size=2.8) +
  scale_fill_gradient2(low=rdbu[1], mid="white", high=rdbu[11], midpoint=0) +
  labs(x=NULL, y=NULL, fill="log2(HR)", title="HR (high vs low) by cohort, endpoint and cutoff method") +
  theme_bw() +
  theme(panel.grid=element_blank(), panel.border=element_blank(),
        plot.title=element_text(hjust=0.5, face="bold", size=10),
        axis.text.x=element_text(angle=0, size=8), axis.text.y=element_text(size=7))
ggsave(file.path(OUT,"figures/forest/hr_heatmap.pdf"), ph, width=6.4, height=4.2, bg="white")
ggsave(file.path(OUT,"figures/forest/hr_heatmap.png"), ph, width=6.4, height=4.2, dpi=300, bg="white")
cat("[Forest] hr_heatmap redrawn (blue=protective, red=risk)\n")
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 核对覆盖矩阵与打分实际用基因的差异，并提取报告需引用的修复后关键数字
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()

ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival")
OUT  <- file.path(ROOT,"processed/survival_v2")
base <- readRDS(file.path(OUT,"data/checkpoints/analysis_ready.rds"))
E    <- readRDS(file.path(OUT,"data/checkpoints/expr_tnbc_log2.rds"))
gs   <- base$gs
alias <- list(SEPP1=c("SEPP1","SELENOP"), ATP5E=c("ATP5E","ATP5F1E"), ATP5F1E=c("ATP5F1E","ATP5E"),
              TCEB2=c("TCEB2","ELOB"), ELOB=c("ELOB","TCEB2"), CTGF=c("CTGF","CCN2"),
              CXCL8=c("CXCL8","IL8"), FYB=c("FYB","FYB1"), RACK1=c("RACK1","GNB2L1"))
alias_of <- function(g) if (!is.null(alias[[g]])) alias[[g]] else g
covmat <- read.csv(file.path(OUT,"tables/signature_gene_coverage_matrix.csv"), stringsAsFactors=FALSE)

cat("--- 平台存在但打分时被清洗掉的基因（重复符号/零方差/含缺失）---\n")
for (co in names(E)) for (sn in names(gs)) {
  avail <- covmat$gene[covmat$cohort==co & covmat$signature==sn & covmat$present==1]
  used <- avail[avail %in% rownames(E[[co]])]
  lost <- setdiff(avail, used)
  cat(sprintf("  %-9s %-10s used %2d/%2d ; lost-in-cleaning: %s\n", co, sn, length(used), length(gs[[sn]]),
              if (length(lost)) paste(lost, collapse=",") else "none"))
}

cat("\n--- 修复后 TCGA cutpoints（ssGSEA）---\n")
cut <- read.csv(file.path(OUT,"tables/stratification_cutpoints.csv"), stringsAsFactors=FALSE)
print(transform(cut[cut$cohort=="TCGA" & cut$score_method=="ssgsea", c("endpoint","cutoff_method","signature","cutpoint","selection_p")],
                cutpoint=round(cutpoint,3), selection_p=signif(selection_p,3)), row.names=FALSE)

cat("\n--- km_summary：关键对比的中位生存（ssGSEA）---\n")
km <- read.csv(file.path(OUT,"tables/km_summary.csv"), stringsAsFactors=FALSE)
k <- km[km$score_method=="ssgsea" & ((km$cohort=="METABRIC" & km$endpoint %in% c("OS","RFS") & km$cutoff_method=="maxstat" &
        km$group %in% c("High","Low","DualHigh","Other")) | (km$cohort=="GSE21653" & km$group %in% c("High","Low"))), ]
k$median_months <- round(k$median_months,1)
print(k[, c("cohort","endpoint","cutoff_method","group","n","events","median_months","lower95","upper95","logrank_p_primary")], row.names=FALSE)

cat("\n--- 修复后主结果（ssGSEA，全部 cohort）---\n")
p <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)
p <- p[p$score_method=="ssgsea", c("cohort","endpoint","cutoff_method","analysis_group","n_high","n_low","events_high","events_low",
                                   "HR","HR_lower","HR_upper","logrank_p","p_effective","p_fdr","HR_adjusted")]
p[, c("HR","HR_lower","HR_upper","HR_adjusted")] <- round(p[, c("HR","HR_lower","HR_upper","HR_adjusted")],3)
p[, c("logrank_p","p_effective","p_fdr")] <- signif(p[, c("logrank_p","p_effective","p_fdr")],3)
print(p, row.names=FALSE)
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 修复后统一重跑：主分析 + 多因素 Cox + PH + 方向校验 + 一致性 + 中位生存 + 新旧相关 + minprop 敏感性
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
suppressPackageStartupMessages({library(survival); library(maxstat)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival"); OUT <- file.path(ROOT,"processed/survival_v2")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
sc <- read.csv(file.path(OUT,"tables/signature_scores.csv"), stringsAsFactors=FALSE)
covars_of <- function(co) switch(co, TCGA=c("age","stage_high"), METABRIC=c("age","stage_high"), GSE58812=c("age"), GSE21653=c("age","grade"))

mk_cut <- function(dat, col, method, minprop=0.25) {
  if (method=="maxstat") {
    dd <- data.frame(time=dat$time, event=dat$event, x=dat[[col]]); dd <- dd[complete.cases(dd),]
    f <- tryCatch(maxstat.test(Surv(time,event)~x, data=dd, smethod="LogRank", pmethod="Lau92", minprop=minprop), error=function(e) NULL)
    if (is.null(f)) return(list(cut=NA_real_, ps=NA_real_))
    return(list(cut=as.numeric(f$estimate), ps=as.numeric(f$p.value)))
  }
  list(cut=median(dat[[col]], na.rm=TRUE), ps=NA_real_)
}
fit_stats <- function(dat, covars) {
  dat <- dat[!is.na(dat$time) & !is.na(dat$event) & dat$time>0 & !is.na(dat$grp), ]
  n <- nrow(dat); ev <- sum(dat$event==1)
  r <- list(n=n, events=ev, nh=sum(dat$grp=="High"), nl=sum(dat$grp=="Low"),
            eh=sum(dat$event[dat$grp=="High"]), el=sum(dat$event[dat$grp=="Low"]),
            hr=NA_real_, lo=NA_real_, hi=NA_real_, pw=NA_real_, lr=NA_real_,
            hra=NA_real_, loa=NA_real_, hia=NA_real_, pa=NA_real_, na=NA_integer_, nex=NA_integer_, php=NA_real_, note="")
  if (n<20 || ev<5) { r$note <- "insufficient (n<20 or events<5)"; return(r) }
  if (min(r$nh,r$nl)<10) { r$note <- "group n<10"; return(r) }
  sd <- survdiff(Surv(time,event)~grp, data=dat); r$lr <- 1-pchisq(sd$chisq, length(sd$n)-1)
  fit <- coxph(Surv(time,event)~grp, data=dat); ci <- exp(confint(fit))
  r$hr <- exp(coef(fit)[["grpHigh"]]); r$lo <- ci["grpHigh",1]; r$hi <- ci["grpHigh",2]
  r$pw <- summary(fit)$coefficients["grpHigh","Pr(>|z|)"]
  cv <- covars[covars %in% names(dat)]
  if (length(cv)) {
    d2 <- dat[complete.cases(dat[, c("time","event","grp",cv)]), ]
    r$nex <- n-nrow(d2); r$na <- nrow(d2)
    if (nrow(d2)>=20 && sum(d2$event)>=5 && min(table(d2$grp))>=10) {
      f2 <- tryCatch(coxph(as.formula(paste("Surv(time,event)~grp+", paste(cv, collapse="+"))), data=d2), error=function(e) NULL)
      if (!is.null(f2)) {
        c2 <- exp(confint(f2)); r$hra <- exp(coef(f2)[["grpHigh"]]); r$loa <- c2["grpHigh",1]; r$hia <- c2["grpHigh",2]
        r$pa <- summary(f2)$coefficients["grpHigh","Pr(>|z|)"]; z <- cox.zph(f2)$table
        r$php <- if ("grp" %in% rownames(z)) as.numeric(z["grp","p"]) else NA_real_
      }
    } else r$note <- paste(r$note, "adjusted skipped", sep="; ")
  }
  r
}
combos <- data.frame(cohort=c("TCGA","TCGA","METABRIC","METABRIC","GSE58812","GSE21653"),
                     endpoint=c("OS","DFS","OS","RFS","OS","DFS"), stringsAsFactors=FALSE)
P <- list(); D <- list(); C <- list(); G <- list(); ip <- 0; id <- 0; ic <- 0; ig <- 0
for (i in seq_len(nrow(combos))) {
  co <- combos$cohort[i]; ep <- combos$endpoint[i]
  d0 <- long[long$cohort==co & long$endpoint==ep, ]; d0 <- d0[!is.na(d0$time)&!is.na(d0$event)&d0$time>0,]
  cvs <- covars_of(co)
  for (meth in c("ssgsea","zmean")) {
    s <- sc[sc$cohort==co & sc$method==meth, ]; d <- merge(d0, s[,c("sample","C1Q_TAM","SPP1_TAM","MMP11_CAF")], by="sample")
    if (nrow(d)<20) next
    for (cm in c("maxstat","median")) {
      c1 <- mk_cut(d,"C1Q_TAM",cm); c2 <- mk_cut(d,"SPP1_TAM",cm); c3 <- mk_cut(d,"MMP11_CAF",cm)
      ic <- ic+1
      C[[ic]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm,
                            signature=c("C1Q_TAM","SPP1_TAM","MMP11_CAF"),
                            cutpoint=c(c1$cut,c2$cut,c3$cut), selection_p=c(c1$ps,c2$ps,c3$ps),
                            n=nrow(d), events=sum(d$event), stringsAsFactors=FALSE)
      g1 <- ifelse(d$C1Q_TAM>c1$cut,"High","Low"); g2 <- ifelse(d$SPP1_TAM>c2$cut,"High","Low"); g3 <- ifelse(d$MMP11_CAF>c3$cut,"High","Low")
      d$C1Q_grp <- factor(g1,levels=c("Low","High")); d$SPP1_grp <- factor(g2,levels=c("Low","High")); d$MMP11_grp <- factor(g3,levels=c("Low","High"))
      d$dual_grp <- ifelse(g1=="High"&g3=="High","DualHigh", ifelse(g1=="Low"&g3=="Low","DualLow","Other"))
      ig <- ig+1
      G[[ig]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm, sample=d$sample,
                            C1Q_grp=d$C1Q_grp, SPP1_grp=d$SPP1_grp, MMP11_grp=d$MMP11_grp, dual_grp=d$dual_grp, stringsAsFactors=FALSE)
      specs <- list(list(lab="C1Q_TAM", dat=transform(d, grp=C1Q_grp), ps=c1$ps),
                    list(lab="SPP1_TAM", dat=transform(d, grp=SPP1_grp), ps=c2$ps),
                    list(lab="MMP11_CAF", dat=transform(d, grp=MMP11_grp), ps=c3$ps),
                    list(lab="C1Q_MMP11_DualHigh_vs_Other",
                         dat=transform(d, grp=factor(ifelse(dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))), ps=NA_real_))
      for (sp in specs) {
        st <- fit_stats(sp$dat, cvs); ip <- ip+1; adj <- cm=="maxstat" && !is.na(sp$ps)
        P[[ip]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm, analysis_group=sp$lab,
          n=st$n, events=st$events, n_high=st$nh, n_low=st$nl, events_high=st$eh, events_low=st$el,
          cutpoint_C1Q=c1$cut, cutpoint_SPP1=c2$cut, cutpoint_MMP11=c3$cut,
          HR=st$hr, HR_lower=st$lo, HR_upper=st$hi, p_wald=st$pw, logrank_p=st$lr,
          p_effective=if (adj) sp$ps else st$lr,
          p_source=if (adj) "maxstat Lau92 (selection-adjusted)" else "log-rank (unadjusted)",
          HR_adjusted=st$hra, HR_adj_lower=st$loa, HR_adj_upper=st$hia, p_adjusted=st$pa,
          n_adjusted=st$na, n_excluded_missing_covariate=st$nex, covariates=paste(cvs,collapse="+"),
          ph_test_p=st$php, note=st$note, stringsAsFactors=FALSE)
      }
      ds <- d[d$dual_grp %in% c("DualHigh","DualLow"), ]; ds$grp <- factor(ifelse(ds$dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))
      s2 <- fit_stats(ds, cvs); id <- id+1
      D[[id]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm,
        analysis_group="C1Q_MMP11_DualHigh_vs_DualLow", n=s2$n, events=s2$events, n_high=s2$nh, n_low=s2$nl,
        events_high=s2$eh, events_low=s2$el, HR=s2$hr, HR_lower=s2$lo, HR_upper=s2$hi, p_wald=s2$pw, logrank_p=s2$lr,
        HR_adjusted=s2$hra, HR_adj_lower=s2$loa, HR_adj_upper=s2$hia, p_adjusted=s2$pa, n_adjusted=s2$na,
        ph_test_p=s2$php, note=s2$note, stringsAsFactors=FALSE)
    }
  }
}
prim <- do.call(rbind,P); dual2 <- do.call(rbind,D); cuts <- do.call(rbind,C); grps <- do.call(rbind,G)
prim$p_fdr <- NA_real_; fam <- paste(prim$score_method, prim$cutoff_method, prim$cohort, prim$endpoint, sep="|")
for (f in unique(fam)) { idx <- which(fam==f & !is.na(prim$p_effective)); if (length(idx)) prim$p_fdr[idx] <- p.adjust(prim$p_effective[idx], "BH") }
prim$significant_fdr <- !is.na(prim$p_fdr) & prim$p_fdr<0.05
write.csv(prim, file.path(OUT,"tables/survival_results_primary.csv"), row.names=FALSE)
write.csv(dual2, file.path(OUT,"tables/survival_results_dual_secondary.csv"), row.names=FALSE)
write.csv(cuts, file.path(OUT,"tables/stratification_cutpoints.csv"), row.names=FALSE)
write.csv(grps, file.path(OUT,"tables/group_assignment_long.csv"), row.names=FALSE)
saveRDS(list(primary=prim, dual_secondary=dual2, cuts=cuts, groups=grps), file.path(OUT,"data/checkpoints/survival_results.rds"))
cat(sprintf("[primary] rows=%d ; fitted=%d ; q<0.05=%d\n", nrow(prim), sum(!is.na(prim$HR)), sum(prim$p_fdr<0.05, na.rm=TRUE)))

# 多因素 Cox 全表
rows <- list(); k <- 0
for (i in seq_len(nrow(prim))) {
  r <- prim[i,]; if (r$score_method!="ssgsea" || is.na(r$HR_adjusted)) next
  d0 <- long[long$cohort==r$cohort & long$endpoint==r$endpoint,]
  g0 <- grps[grps$cohort==r$cohort & grps$endpoint==r$endpoint & grps$score_method=="ssgsea" & grps$cutoff_method==r$cutoff_method,]
  d <- merge(d0[,c("sample","time","event","age","stage_high","grade")], g0, by="sample")
  if (r$analysis_group=="C1Q_MMP11_DualHigh_vs_Other") { d <- d[d$dual_grp!="Other",]; d$grp <- factor(ifelse(d$dual_grp=="DualHigh","High","Low"),levels=c("Low","High"))
  } else { col <- switch(r$analysis_group, C1Q_TAM="C1Q_grp", SPP1_TAM="SPP1_grp", MMP11_CAF="MMP11_grp"); d$grp <- factor(d[[col]],levels=c("Low","High")) }
  cv <- covars_of(r$cohort); cv <- cv[cv %in% names(d)]; d <- d[complete.cases(d[,c("time","event","grp",cv)]),]
  f <- tryCatch(coxph(as.formula(paste("Surv(time,event)~grp+", paste(cv,collapse="+"))), data=d), error=function(e) NULL); if (is.null(f)) next
  ss <- summary(f)
  for (tm in rownames(ss$coefficients)) { k <- k+1
    rows[[k]] <- data.frame(cohort=r$cohort, endpoint=r$endpoint, cutoff_method=r$cutoff_method, analysis_group=r$analysis_group,
      term=tm, n=nrow(d), events=sum(d$event), HR=ss$conf.int[tm,"exp(coef)"], HR_lower=ss$conf.int[tm,"lower .95"],
      HR_upper=ss$conf.int[tm,"upper .95"], p=ss$coefficients[tm,"Pr(>|z|)"], stringsAsFactors=FALSE) }
}
coxfull <- do.call(rbind, rows)
coxfull$q_within_model <- ave(coxfull$p, paste(coxfull$cohort,coxfull$endpoint,coxfull$cutoff_method,coxfull$analysis_group), FUN=function(x) p.adjust(x,"BH"))
write.csv(coxfull, file.path(OUT,"tables/cox_multivariable_full.csv"), row.names=FALSE)
cat(sprintf("[cox full] rows=%d\n", nrow(coxfull)))

# PH + 方向
ph <- prim[,c("cohort","endpoint","score_method","cutoff_method","analysis_group","n_adjusted","ph_test_p")]
ph$ph_violation <- !is.na(ph$ph_test_p) & ph$ph_test_p<0.05
write.csv(ph, file.path(OUT,"tables/ph_assumption_tests.csv"), row.names=FALSE)
cat(sprintf("[PH] models=%d ; violations=%d\n", sum(!is.na(ph$ph_test_p)), sum(ph$ph_violation,na.rm=TRUE)))
E <- readRDS(file.path(OUT,"data/checkpoints/expr_tnbc_log2.rds"))
anchor <- list(C1Q_TAM="C1QA", SPP1_TAM="SPP1", MMP11_CAF="MMP11"); dr <- list(); k <- 0
for (co in names(E)) for (scn in names(anchor)) { g <- anchor[[scn]]; if (!(g %in% rownames(E[[co]]))) next
  s <- sc[sc$cohort==co & sc$method=="ssgsea",]; s <- s[match(colnames(E[[co]]), s$sample),]
  rho <- suppressWarnings(cor(s[[scn]], as.numeric(E[[co]][g,]), method="spearman", use="complete.obs")); k <- k+1
  dr[[k]] <- data.frame(cohort=co, signature=scn, anchor_gene=g, spearman_rho=rho, direction_ok=!is.na(rho)&&rho>0, stringsAsFactors=FALSE) }
dchk <- do.call(rbind, dr); write.csv(dchk, file.path(OUT,"tables/direction_check.csv"), row.names=FALSE)

# 两种一致性
cons <- list(); k <- 0
for (co in unique(prim$cohort)) for (ep in unique(prim$endpoint[prim$cohort==co])) for (meth in c("ssgsea","zmean")) for (ag in unique(prim$analysis_group)) {
  a <- prim[prim$cohort==co&prim$endpoint==ep&prim$score_method==meth&prim$cutoff_method=="maxstat"&prim$analysis_group==ag,]
  b <- prim[prim$cohort==co&prim$endpoint==ep&prim$score_method==meth&prim$cutoff_method=="median"&prim$analysis_group==ag,]
  if (!nrow(a)||!nrow(b)) next; k <- k+1
  cons[[k]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, analysis_group=ag, HR_maxstat=a$HR, p_maxstat=a$p_effective, q_maxstat=a$p_fdr,
    HR_median=b$HR, p_median=b$logrank_p, q_median=b$p_fdr,
    direction_concordant=ifelse(is.na(a$HR)||is.na(b$HR), NA, sign(log(a$HR))==sign(log(b$HR))),
    sig_maxstat=!is.na(a$p_fdr)&a$p_fdr<0.05, sig_median=!is.na(b$p_fdr)&b$p_fdr<0.05, stringsAsFactors=FALSE) }
cons <- do.call(rbind, cons); cons$significance_concordant <- cons$sig_maxstat==cons$sig_median
write.csv(cons, file.path(OUT,"tables/cutoff_consistency.csv"), row.names=FALSE)
sens <- list(); k <- 0
for (co in unique(prim$cohort)) for (ep in unique(prim$endpoint[prim$cohort==co])) for (cm in c("maxstat","median")) for (ag in unique(prim$analysis_group)) {
  a <- prim[prim$cohort==co&prim$endpoint==ep&prim$score_method=="ssgsea"&prim$cutoff_method==cm&prim$analysis_group==ag,]
  b <- prim[prim$cohort==co&prim$endpoint==ep&prim$score_method=="zmean"&prim$cutoff_method==cm&prim$analysis_group==ag,]
  if (!nrow(a)||!nrow(b)) next; k <- k+1
  sens[[k]] <- data.frame(cohort=co, endpoint=ep, cutoff_method=cm, analysis_group=ag, HR_ssgsea=a$HR, q_ssgsea=a$p_fdr,
    HR_zmean=b$HR, q_zmean=b$p_fdr, direction_concordant=ifelse(is.na(a$HR)||is.na(b$HR), NA, sign(log(a$HR))==sign(log(b$HR))), stringsAsFactors=FALSE) }
sens <- do.call(rbind, sens); write.csv(sens, file.path(OUT,"tables/score_method_consistency.csv"), row.names=FALSE)
cat(sprintf("[consistency] cutoff dir=%d/%d sig=%d ; score dir=%d/%d both-sig=%d\n",
  sum(cons$direction_concordant,na.rm=TRUE), nrow(cons), sum(cons$significance_concordant),
  sum(sens$direction_concordant,na.rm=TRUE), nrow(sens),
  sum(!is.na(sens$q_ssgsea)&!is.na(sens$q_zmean)&sens$q_ssgsea<0.05&sens$q_zmean<0.05)))

# 中位生存
km <- list(); k <- 0
for (co in unique(grps$cohort)) for (ep in unique(grps$endpoint[grps$cohort==co])) for (meth in c("ssgsea","zmean")) for (cm in c("maxstat","median")) {
  d0 <- long[long$cohort==co&long$endpoint==ep,]; g0 <- grps[grps$cohort==co&grps$endpoint==ep&grps$score_method==meth&grps$cutoff_method==cm,]
  if (!nrow(g0)) next; d <- merge(d0[,c("sample","time","event")], g0, by="sample"); d <- d[!is.na(d$time)&!is.na(d$event)&d$time>0,]
  medfun <- function(f) { t <- summary(f)$table; data.frame(n=t[,"records"], events=t[,"events"], median=t[,"median"], lower=t[,"0.95LCL"], upper=t[,"0.95UCL"]) }
  for (gn in c("C1Q_grp","SPP1_grp","MMP11_grp")) {
    dd <- d; dd$grp <- factor(dd[[gn]],levels=c("Low","High")); mt <- medfun(survfit(Surv(time,event)~grp,data=dd))
    p <- 1-pchisq(survdiff(Surv(time,event)~grp,data=dd)$chisq,1)
    for (lv in rownames(mt)) { k <- k+1; km[[k]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm, group=sub("grp=","",lv),
      n=mt[lv,"n"], events=mt[lv,"events"], median_months=mt[lv,"median"], lower95=mt[lv,"lower"], upper95=mt[lv,"upper"],
      logrank_p_primary=p, logrank_p_secondary=NA_real_, stringsAsFactors=FALSE) } }
  dd <- d; dd$grp <- factor(dd$dual_grp,levels=c("DualHigh","Other","DualLow")); mt <- medfun(survfit(Surv(time,event)~grp,data=dd))
  p1 <- 1-pchisq(survdiff(Surv(time,event)~grp,data=dd[dd$dual_grp %in% c("DualHigh","Other"),])$chisq,1)
  p2 <- 1-pchisq(survdiff(Surv(time,event)~grp,data=dd[dd$dual_grp %in% c("DualHigh","DualLow"),])$chisq,1)
  for (lv in rownames(mt)) { k <- k+1; km[[k]] <- data.frame(cohort=co, endpoint=ep, score_method=meth, cutoff_method=cm, group=sub("grp=","",lv),
    n=mt[lv,"n"], events=mt[lv,"events"], median_months=mt[lv,"median"], lower95=mt[lv,"lower"], upper95=mt[lv,"upper"],
    logrank_p_primary=p1, logrank_p_secondary=p2, stringsAsFactors=FALSE) } }
km_sum <- do.call(rbind, km); write.csv(km_sum, file.path(OUT,"tables/km_summary.csv"), row.names=FALSE)
cat(sprintf("[km_summary] rows=%d\n", nrow(km_sum)))

# 新旧 signature 相关
old <- read.csv(file.path(ROOT,"processed/scores.csv"), stringsAsFactors=FALSE); old_sp <- old[grepl("_specific$", old$signature),]
map_sig <- function(x) ifelse(grepl("^C1Q",x),"C1Q_TAM", ifelse(grepl("^SPP1",x),"SPP1_TAM", ifelse(grepl("^MMP11",x),"MMP11_CAF",NA)))
cr <- list(); k <- 0
for (co in unique(old_sp$cohort)) for (osig in unique(old_sp$signature[old_sp$cohort==co])) {
  nsig <- map_sig(osig); if (is.na(nsig)) next
  o <- old_sp[old_sp$cohort==co & old_sp$signature==osig, c("sample","ssgsea","zscore_mean")]
  a <- sc[sc$cohort==co & sc$method=="ssgsea", c("sample",nsig)]; names(a)[2] <- "new_ss"
  b <- sc[sc$cohort==co & sc$method=="zmean", c("sample",nsig)]; names(b)[2] <- "new_z"
  m <- merge(merge(o,a,by="sample"), b, by="sample"); k <- k+1
  cr[[k]] <- data.frame(cohort=co, new_signature=nsig, old_signature=osig, n_matched=nrow(m),
    rho_new_ssgsea_vs_old_ssgsea=if (nrow(m)>3) cor(m$new_ss, m$ssgsea, method="spearman") else NA,
    rho_new_ssgsea_vs_old_zmean=if (nrow(m)>3) cor(m$new_ss, m$zscore_mean, method="spearman") else NA,
    rho_new_zmean_vs_old_zmean=if (nrow(m)>3) cor(m$new_z, m$zscore_mean, method="spearman") else NA, stringsAsFactors=FALSE) }
cort <- do.call(rbind, cr); write.csv(cort, file.path(OUT,"tables/direction_check_vs_old_signature.csv"), row.names=FALSE)
cat(sprintf("[old-signature corr] rho range %.3f - %.3f\n", min(cort$rho_new_ssgsea_vs_old_ssgsea), max(cort$rho_new_ssgsea_vs_old_ssgsea)))

# minprop 0.10 vs 0.25
key <- data.frame(cohort=c("METABRIC","METABRIC","GSE21653"), endpoint=c("OS","RFS","DFS"), stringsAsFactors=FALSE)
cmp <- list(); k <- 0
for (i in seq_len(nrow(key))) {
  d0 <- long[long$cohort==key$cohort[i]&long$endpoint==key$endpoint[i],]
  s <- sc[sc$cohort==key$cohort[i]&sc$method=="ssgsea", c("sample","C1Q_TAM")]
  d <- merge(d0[,c("sample","time","event")], s, by="sample"); d <- d[!is.na(d$time)&!is.na(d$event)&d$time>0,]
  for (mp in c(0.10,0.25)) {
    f <- maxstat.test(Surv(time,event)~C1Q_TAM, data=d, smethod="LogRank", pmethod="Lau92", minprop=mp)
    cut <- as.numeric(f$estimate); d$grp <- factor(ifelse(d$C1Q_TAM>cut,"High","Low"),levels=c("Low","High"))
    cx <- coxph(Surv(time,event)~grp,data=d); ci <- exp(confint(cx)); k <- k+1
    cmp[[k]] <- data.frame(cohort=key$cohort[i], endpoint=key$endpoint[i], minprop=mp, cutpoint=cut,
      n_high=sum(d$grp=="High"), n_low=sum(d$grp=="Low"), HR=exp(coef(cx)[["grpHigh"]]),
      HR_lower=ci["grpHigh",1], HR_upper=ci["grpHigh",2], p_adjusted=f$p.value, stringsAsFactors=FALSE) } }
cmp <- do.call(rbind, cmp); write.csv(cmp, file.path(OUT,"tables/maxstat_minprop_sensitivity.csv"), row.names=FALSE)
cat("[minprop] cutpoints identical:", all(tapply(cmp$cutpoint, paste(cmp$cohort,cmp$endpoint), function(x) length(unique(round(x,6)))==1)), "\n")
cat("\n--- q<0.05 汇总 ---\n")
q <- prim[!is.na(prim$p_fdr)&prim$p_fdr<0.05, c("cohort","endpoint","score_method","cutoff_method","analysis_group","n_high","n_low","HR","HR_lower","HR_upper","p_effective","p_fdr","HR_adjusted")]
q[,c("HR","HR_lower","HR_upper","HR_adjusted")] <- round(q[,c("HR","HR_lower","HR_upper","HR_adjusted")],3); q[,c("p_effective","p_fdr")] <- signif(q[,c("p_effective","p_fdr")],3)
print(q, row.names=FALSE)
cat(sprintf("\n[elapsed] %.1f s\n", as.numeric(difftime(Sys.time(), .start, units="secs"))))
cat("===== Cell End:", date(), "=====\n")


# 打印修复后基因覆盖表（修正 cat 用法）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
OUT <- file.path(Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival"), "processed/survival_v2")
print(read.csv(file.path(OUT,"tables/signature_gene_coverage.csv"), stringsAsFactors=FALSE), row.names=FALSE)
cat("\n[pixel/consistency tables]\n")
print(read.csv(file.path(OUT,"tables/direction_check_vs_old_signature.csv"), stringsAsFactors=FALSE)[, c("cohort","new_signature","n_matched","rho_new_ssgsea_vs_old_ssgsea")], row.names=FALSE)
print(read.csv(file.path(OUT,"tables/maxstat_minprop_sensitivity.csv"), stringsAsFactors=FALSE)[, c("cohort","endpoint","minprop","cutpoint","HR","p_adjusted")], row.names=FALSE)
cat("\n===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 用修复后打分重绘全部 KM 与森林图（NPG 定性色 + RdBu 发散色）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
suppressPackageStartupMessages({library(survival); library(survminer); library(ggplot2); library(ggsci); library(RColorBrewer)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival"); OUT <- file.path(ROOT,"processed/survival_v2")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
grp <- read.csv(file.path(OUT,"tables/group_assignment_long.csv"), stringsAsFactors=FALSE)
prim <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)
npg2 <- unname(ggsci::pal_npg("nrc")(2)); pal4 <- unname(ggsci::pal_npg("nrc")(4))
rdbu <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(11)
save_km <- function(pl, stem) { pdf(paste0(stem,".pdf"),width=6,height=5.6,bg="white"); print(pl); dev.off()
  png(paste0(stem,".png"),width=1800,height=1680,res=300,bg="white"); print(pl); dev.off() }
gmeta <- list(C1Q_grp=list(col="C1Q_grp",lab="C1Q TAMs"), SPP1_grp=list(col="SPP1_grp",lab="SPP1 TAMs"),
              MMP11_grp=list(col="MMP11_grp",lab="MMP11 CAFs"), dual_grp=list(col="dual_grp",lab="C1Q+ and MMP11+"))
nk <- 0
combos <- unique(grp[, c("cohort","endpoint","score_method","cutoff_method")])
for (i in seq_len(nrow(combos))) {
  co <- combos$cohort[i]; ep <- combos$endpoint[i]; if (combos$score_method[i]!="ssgsea") next
  cm <- combos$cutoff_method[i]
  d0 <- long[long$cohort==co & long$endpoint==ep,]
  g0 <- grp[grp$cohort==co & grp$endpoint==ep & grp$score_method=="ssgsea" & grp$cutoff_method==cm,]
  d <- merge(d0[,c("sample","time","event")], g0, by="sample"); d <- d[!is.na(d$time)&!is.na(d$event)&d$time>0,]
  pr <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method=="ssgsea" & prim$cutoff_method==cm,]
  for (gn in names(gmeta)) {
    dd <- d
    if (gn=="dual_grp") { dd <- dd[dd$dual_grp!="Other",]; dd$grp <- factor(ifelse(dd$dual_grp=="DualHigh","High","Low"),levels=c("Low","High"))
    } else dd$grp <- factor(dd[[gmeta[[gn]]$col]], levels=c("Low","High"))
    dd <- dd[!is.na(dd$grp),]; if (nrow(dd)<10 || sum(dd$event)<5) next
    lab <- if (gn=="dual_grp") "C1Q+MMP11+ vs others" else gmeta[[gn]]$lab
    row <- pr[pr$analysis_group == if (gn=="dual_grp") "C1Q_MMP11_DualHigh_vs_Other" else
              switch(gn, C1Q_grp="C1Q_TAM", SPP1_grp="SPP1_TAM", MMP11_grp="MMP11_CAF"),]
    sub <- if (!nrow(row) || is.na(row$HR[1])) "model not fitted (small groups)" else
      sprintf("HR=%.2f (95%%CI %.2f-%.2f); log-rank p=%.3g; adj p=%.3g", row$HR[1], row$HR_lower[1], row$HR_upper[1], row$logrank_p[1], row$p_effective[1])
    pl <- ggsurvplot(survfit(Surv(time,event)~grp, data=dd), data=dd, palette=npg2, conf.int=FALSE,
      risk.table=TRUE, risk.table.height=0.24, censor=TRUE, legend.title="", legend.labs=c("Low","High"),
      xlab="Time (months)", ylab="Survival probability", title=sprintf("%s %s - %s (%s cutoff)", co, ep, lab, cm),
      subtitle=sub, font.subtitle=c(9,"plain"),
      ggtheme=theme_bw()+theme(panel.grid=element_blank(), plot.title=element_text(hjust=0.5,face="bold",size=10),
                               plot.subtitle=element_text(hjust=0.5), legend.position="right"))
    save_km(pl, file.path(OUT,"figures/km", sprintf("km_%s_%s_%s_%s", co, ep, cm, gn))); nk <- nk+1
  }
}
cat(sprintf("[KM] redrawn %d\n", nk))
lab_map <- c(C1Q_TAM="C1Q TAMs", SPP1_TAM="SPP1 TAMs", MMP11_CAF="MMP11 CAFs", C1Q_MMP11_DualHigh_vs_Other="C1Q+ & MMP11+ (vs other)")
prim$grp_lab <- lab_map[prim$analysis_group]
prim$star <- ifelse(is.na(prim$p_fdr),"", ifelse(prim$p_fdr<0.001,"***", ifelse(prim$p_fdr<0.01,"**", ifelse(prim$p_fdr<0.05,"*",""))))
forest_panel <- function(df, xlab_txt, title_txt) {
  df <- df[!is.na(df$HR)&!is.na(df$HR_lower),]; if (!nrow(df)) return(NULL)
  df$row <- factor(df$rowlab, levels=rev(unique(df$rowlab)))
  ggplot(df, aes(x=HR,y=row,color=grp_lab)) + geom_vline(xintercept=1,linetype="dashed",color="grey40",linewidth=0.4) +
    geom_errorbarh(aes(xmin=HR_lower,xmax=HR_upper), height=0.18, linewidth=0.5) + geom_point(shape=15,size=2.4) +
    geom_text(aes(label=paste0(sprintf("%.2f",HR),star)), vjust=-0.9, size=2.9, show.legend=FALSE) + scale_x_log10() +
    scale_color_manual(values=setNames(pal4[seq_along(unique(df$grp_lab))], unique(df$grp_lab))) +
    labs(x=xlab_txt,y=NULL,title=title_txt,color=NULL) + theme_bw() +
    theme(panel.grid=element_blank(), panel.border=element_rect(colour="black",fill=NA),
          plot.title=element_text(hjust=0.5,face="bold",size=10), legend.position="bottom",
          legend.text=element_text(size=8), axis.text.y=element_text(size=8))
}
n1 <- 0; n2 <- 0
for (co in unique(prim$cohort)) for (ep in unique(prim$endpoint[prim$cohort==co])) for (cm in c("maxstat","median")) {
  sub <- prim[prim$cohort==co & prim$endpoint==ep & prim$cutoff_method==cm & prim$score_method=="ssgsea",]; if (!nrow(sub)) next
  sub$rowlab <- paste0(sub$grp_lab, "  [n=", sub$n_high+sub$n_low, "]")
  p1 <- forest_panel(sub, "Hazard ratio (high vs low, log scale)", sprintf("%s %s - %s cutoff", co, ep, cm))
  if (!is.null(p1)) { stem <- file.path(OUT,"figures/forest", sprintf("forest_all_%s_%s_%s",co,ep,cm))
    ggsave(paste0(stem,".pdf"),p1,width=5.4,height=2.6,bg="white"); ggsave(paste0(stem,".png"),p1,width=5.4,height=2.6,dpi=300,bg="white"); n1 <- n1+1 }
  d2 <- rbind(data.frame(grp_lab=sub$grp_lab,HR=sub$HR,HR_lower=sub$HR_lower,HR_upper=sub$HR_upper,star=sub$star,model="unadjusted",stringsAsFactors=FALSE),
              data.frame(grp_lab=sub$grp_lab,HR=sub$HR_adjusted,HR_lower=sub$HR_adj_lower,HR_upper=sub$HR_adj_upper,
                         star=ifelse(is.na(sub$p_adjusted),"",ifelse(sub$p_adjusted<0.001,"***",ifelse(sub$p_adjusted<0.01,"**",ifelse(sub$p_adjusted<0.05,"*","")))),
                         model="adjusted",stringsAsFactors=FALSE))
  d2 <- d2[!is.na(d2$HR)&!is.na(d2$HR_lower),]
  if (nrow(d2)) { d2$rowlab <- paste0(d2$grp_lab," - ",d2$model)
    p2 <- forest_panel(d2,"Hazard ratio (high vs low, log scale)", sprintf("%s %s - %s cutoff (adjusted: %s)",co,ep,cm,paste(unique(sub$covariates),collapse=",")))
    if (!is.null(p2)) { stem <- file.path(OUT,"figures/forest", sprintf("forest_adjusted_%s_%s_%s",co,ep,cm))
      ggsave(paste0(stem,".pdf"),p2,width=5.6,height=3.6,bg="white"); ggsave(paste0(stem,".png"),p2,width=5.6,height=3.6,dpi=300,bg="white"); n2 <- n2+1 } }
}
hm <- prim[prim$score_method=="ssgsea" & !is.na(prim$HR),]; hm$key <- paste(hm$cohort,hm$endpoint,hm$cutoff_method,sep=" | ")
hm$grp_lab <- factor(hm$grp_lab, levels=unname(lab_map))
ph <- ggplot(hm, aes(x=grp_lab,y=key,fill=log2(HR))) + geom_tile(color=NA) +
  geom_text(aes(label=ifelse(is.na(star)|star=="",sprintf("%.2f",HR),paste0(sprintf("%.2f",HR),star))), size=2.8) +
  scale_fill_gradient2(low=rdbu[1],mid="white",high=rdbu[11],midpoint=0) +
  labs(x=NULL,y=NULL,fill="log2(HR)",title="HR (high vs low) by cohort, endpoint and cutoff method") + theme_bw() +
  theme(panel.grid=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0.5,face="bold",size=10),
        axis.text.x=element_text(angle=0,size=8), axis.text.y=element_text(size=7))
ggsave(file.path(OUT,"figures/forest/hr_heatmap.pdf"), ph, width=6.4, height=4.2, bg="white")
ggsave(file.path(OUT,"figures/forest/hr_heatmap.png"), ph, width=6.4, height=4.2, dpi=300, bg="white")
cat(sprintf("[Forest] all=%d adjusted=%d heatmap=1\n", n1, n2))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")


# 修复后最终核对：图件完整性、交付清单与关键表行数
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
suppressPackageStartupMessages(library(png))
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival"); OUT <- file.path(ROOT,"processed/survival_v2")
pngs <- list.files(file.path(OUT,"figures"), pattern="\\.png$", recursive=TRUE, full.names=TRUE)
pdfs <- list.files(file.path(OUT,"figures"), pattern="\\.pdf$", recursive=TRUE, full.names=TRUE)
stem <- function(x) sub("\\.(pdf|png)$","",x)
cat(sprintf("[figures] pdf=%d png=%d ; unpaired pdf=%d png=%d ; min png=%.0fKB\n",
            length(pdfs), length(pngs), length(setdiff(stem(pdfs),stem(pngs))), length(setdiff(stem(pngs),stem(pdfs))),
            min(file.info(pngs)$size)/1024))
set.seed(7); chk <- sample(pngs, 15)
m <- do.call(rbind, lapply(chk, function(f) { img <- readPNG(f); v <- as.vector(if (length(dim(img))==3) img[,,1:3] else img)
  data.frame(file=basename(f), sd=sd(v), ink=mean(v<0.98)) }))
m$pass <- m$sd>0.02 & m$ink>0.001
cat(sprintf("[pixel] %d/%d pass ; sd range %.3f-%.3f\n", sum(m$pass), nrow(m), min(m$sd), max(m$sd)))
if (any(!m$pass)) print(m[!m$pass,], row.names=FALSE)

tabs <- c("sample_flow.csv","signature_scores.csv","signature_gene_coverage.csv","signature_gene_coverage_matrix.csv",
          "survival_inputs_long.csv","stratification_cutpoints.csv","survival_results_primary.csv",
          "survival_results_dual_secondary.csv","km_summary.csv","cox_multivariable_full.csv","ph_assumption_tests.csv",
          "direction_check.csv","direction_check_vs_old_signature.csv","cutoff_consistency.csv",
          "score_method_consistency.csv","maxstat_minprop_sensitivity.csv","group_assignment_long.csv")
cat("\n[tables]\n")
for (t in tabs) { f <- file.path(OUT,"tables",t)
  cat(sprintf("  %-42s %s %5d rows %7.1f KB\n", t, ifelse(file.exists(f),"OK","MISSING"),
              ifelse(file.exists(f), length(readLines(f))-1, 0), ifelse(file.exists(f), file.info(f)$size/1024, 0))) }
cat(sprintf("\n[doc] %s : %.1f KB\n", "step_survival_v2_汇总.html", file.info(file.path(OUT,"step_survival_v2_汇总.html"))$size/1024))
p <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)
cat(sprintf("[primary] rows=%d fitted=%d not-fitted=%d q<0.05=%d\n", nrow(p), sum(!is.na(p$HR)), sum(is.na(p$HR)), sum(p$p_fdr<0.05,na.rm=TRUE)))
cat("===== Cell End:", date(), "elapsed:", format(Sys.time() - .start, digits=3), "=====\n")
