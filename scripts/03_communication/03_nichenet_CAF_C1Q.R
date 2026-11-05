# =============================================================================
# TNBC interface-niche code repository
# Original source: FigureAI/scripts/03_nichenet_C1Q_MMP11CAF.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

#!/usr/bin/env Rscript
## ============================================================
## NicheNet — C1Q TAM  <->  MMP11 CAF  (spatial, per zone)
##
## Question
##   Direction A: C1Q TAM as LIGAND SENDER  ->  MMP11 CAF as receiver
##   Direction B: MMP11 CAF as LIGAND SENDER ->  C1Q TAM as receiver
##   For each of Tumor Core / Interface / Stroma, and compare zones.
##
## Data (read-only)
##   SCT log-normalised spatial matrix : 23,947 genes x 16,607 spots
##   RCTD subcluster proportions       : 35 subclusters per spot
##   Zone assignment                   : Tumor Core / Interface / Stroma
##
## Receiver geneset (user choice): genes whose expression, among receiver-
##   positive spots of that zone, correlates positively with the SENDER
##   proportion. This asks "which genes in the receiver track sender abundance".
##
## Sender / receiver spots
##   C1Q TAM-positive spot     : C1Q TAMs proportion >= 0.02
##   MMP11 CAF-positive spot   : MMP11 CAFs proportion >= 0.05
##
## Reproduce
##   /path/to/envs/common/bin/Rscript scripts/03_communication/03_nichenet_CAF_C1Q.R
## ============================================================

suppressPackageStartupMessages({
  library(Matrix); library(nichenetr); library(ggplot2)
  library(dplyr); library(tidyr); library(tibble)
})

set.seed(42)

## ---------------- paths ----------------
ROOT <- "/path/to/project/project_spatial/FigureAI"
SRC  <- "/path/to/project/project_singlecell/TNBC_result"
SCT  <- file.path(SRC, "final_figures/V2/spatial/spatial_sct_matrix.rds")
PROP <- file.path(SRC, "figure18/intermediates/subcluster_prop.csv.gz")
REG  <- file.path(SRC, "_archive_legacy/figure2/region_assignment_extended.csv")
LT   <- "/path/to/data/other_ans_project/SeekSpace_code/ligand_target_matrix.rds"
LRN  <- "/path/to/data/other_ans_project/SeekSpace_code/lr_network.rds"

OUT  <- file.path(ROOT, "results/nichenet_C1Q_MMP11CAF")
TDIR <- file.path(OUT, "tables"); FDIR <- file.path(OUT, "figures")
dir.create(TDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FDIR, showWarnings = FALSE, recursive = TRUE)

## ---------------- style ----------------
FONT <- "sans"
theme_set(theme_classic(base_size = 7, base_family = FONT) +
            theme(axis.text = element_text(size = 6.4, colour = "black"),
                  axis.title = element_text(size = 7.5, colour = "black"),
                  plot.title = element_text(size = 8.4, face = "bold", hjust = 0),
                  plot.subtitle = element_text(size = 6.2, colour = "grey30", hjust = 0),
                  strip.text = element_text(size = 6.4, face = "bold"),
                  legend.title = element_text(size = 6.4),
                  legend.text = element_text(size = 6),
                  legend.key.size = unit(2.6, "mm"),
                  panel.grid = element_blank()))

## ---------------- data ----------------
cat("[1] loading\n")
sct  <- readRDS(SCT)
prop <- read.csv(gzfile(PROP), check.names = FALSE)
reg  <- read.csv(REG, check.names = FALSE)
stopifnot(all(colnames(sct) %in% prop$barcode), all(colnames(sct) %in% reg$barcode))
prop <- prop[match(colnames(sct), prop$barcode), ]
reg  <- reg[match(colnames(sct), reg$barcode), ]
stopifnot(identical(prop$barcode, colnames(sct)), identical(reg$barcode, colnames(sct)))

lt  <- readRDS(LT)
lrn <- readRDS(LRN)

C1Q  <- as.numeric(prop[["C1Q TAMs"]])
MMP  <- as.numeric(prop[["MMP11 CAFs"]])
zone <- as.character(reg$region)
zone[zone == "Tumor Core"] <- "Tumor Core"        # keep original labels
samples <- reg$sample

TH_C1Q <- 0.02
TH_MMP <- 0.05
ZONES  <- c("Tumor Core", "Interface", "Stroma")

spot_counts <- data.frame(
  zone = ZONES,
  n_spots = sapply(ZONES, function(z) sum(zone == z)),
  n_C1Q_pos = sapply(ZONES, function(z) sum(zone == z & C1Q >= TH_C1Q)),
  n_MMP11_pos = sapply(ZONES, function(z) sum(zone == z & MMP >= TH_MMP)))
write.csv(spot_counts, file.path(TDIR, "00_spot_counts.csv"), row.names = FALSE)
print(spot_counts)

## ---------------- helpers ----------------
## Pearson r of every gene against a numeric vector, computed from sparse sums
cor_with_vector <- function(X, y) {
  n  <- ncol(X)
  my <- mean(y); sy2 <- sum((y - my)^2)
  sx  <- Matrix::rowSums(X)
  sx2 <- Matrix::rowSums(X * X)
  sxy <- as.numeric(X %*% y)
  num <- sxy - n * (sx / n) * my
  den <- sqrt(pmax(sx2 - n * (sx / n)^2, 0)) * sqrt(sy2)
  r <- ifelse(den > 0, num / den, NA_real_)
  p <- ifelse(is.na(r), NA_real_,
              2 * pt(-abs(r * sqrt((n - 2) / pmax(1 - r^2, 1e-12))), df = n - 2))
  data.frame(gene = rownames(X), r = r, p = p,
             frac = Matrix::rowMeans(X > 0), stringsAsFactors = FALSE)
}

run_one <- function(zone_name, sender_prop, receiver_prop,
                    sender_lab, receiver_lab, sender_th, receiver_th, tag) {
  idx_z  <- zone == zone_name
  idx_s  <- which(idx_z & sender_prop   >= sender_th)
  idx_r  <- which(idx_z & receiver_prop >= receiver_th)
  n_s <- length(idx_s); n_r <- length(idx_r)
  cat(sprintf("  [%s] %s->%s | sender n=%d | receiver n=%d\n",
              zone_name, sender_lab, receiver_lab, n_s, n_r))
  if (n_s < 30 || n_r < 30) {
    cat("    [skip] fewer than 30 sender or receiver spots\n")
    return(NULL)
  }

  Xr <- sct[, idx_r, drop = FALSE]
  ## receiver enrichment: keep genes higher in receiver-positive than -negative
  ## spots of the same zone, so the geneset reflects the receiver rather than
  ## transcripts of the sender that merely co-localise (spot-mixture control).
  idx_rn <- which(idx_z & receiver_prop < receiver_th)
  if (length(idx_rn) >= 20) {
    mu_pos <- Matrix::rowMeans(Xr)
    mu_neg <- Matrix::rowMeans(sct[, idx_rn, drop = FALSE])
    rec_lfc <- mu_pos - mu_neg
  } else {
    rec_lfc <- setNames(rep(NA_real_, nrow(sct)), rownames(sct))
  }
  ## receiver geneset: correlation with sender abundance inside receiver spots
  cr <- cor_with_vector(Xr, sender_prop[idx_r])
  cr$padj <- p.adjust(cr$p, method = "BH")
  cr$receiver_lfc <- rec_lfc[cr$gene]
  ## sender enrichment, used only as a specificity filter below
  idx_sn <- which(idx_z & sender_prop < sender_th)
  if (length(idx_sn) >= 20) {
    sn_lfc <- Matrix::rowMeans(sct[, idx_s, drop = FALSE]) -
      Matrix::rowMeans(sct[, idx_sn, drop = FALSE])
  } else {
    sn_lfc <- setNames(rep(NA_real_, nrow(sct)), rownames(sct))
  }
  cr$sender_lfc <- sn_lfc[cr$gene]
  cr <- cr[!is.na(cr$r) & cr$frac >= 0.10, ]
  ## receiver-enriched AND more receiver-specific than sender-specific:
  ## removes genes that merely co-localise with the sender in the same spot
  cr <- cr[!is.na(cr$receiver_lfc) & cr$receiver_lfc >= 0.25, ]
  ## NOTE: sender_lfc is kept as a diagnostic column, not a hard filter. A strict
  ## "receiver_lfc - sender_lfc >= 0.1" filter collapses the geneset to 0-88 genes
  ## because the two populations are strongly co-localised at spot resolution
  ## (see summary); the filter is therefore reported as a sensitivity check.
  cr <- cr[order(-cr$r), ]
  write.csv(cr, file.path(TDIR, sprintf("geneset_%s_%s.csv", tag, gsub(" ", "", zone_name))),
            row.names = FALSE)

  geneset <- cr$gene[cr$r >= 0.30 & cr$padj < 0.05]
  geneset <- head(intersect(geneset, rownames(lt)), 200)
  if (length(geneset) < 20) {
    cat(sprintf("    [skip] geneset too small (%d genes)\n", length(geneset)))
    return(NULL)
  }

  ## ligands expressed in sender spots
  Xs <- sct[, idx_s, drop = FALSE]
  lig_frac <- Matrix::rowMeans(Xs > 0)
  lig_mean <- Matrix::rowMeans(Xs)
  potential <- intersect(unique(lrn$from), rownames(sct))
  potential <- potential[lig_frac[potential] >= 0.10 & lig_mean[potential] > 0.05]
  potential <- intersect(potential, colnames(lt))

  ## background = every gene expressed in >=10% of receiver spots (independent
  ## of the geneset filter, so AUPR stays comparable across zones)
  background <- rownames(Xr)[Matrix::rowMeans(Xr > 0) >= 0.10]
  background <- intersect(background, rownames(lt))

  cat(sprintf("    geneset=%d | potential ligands=%d | background=%d\n",
              length(geneset), length(potential), length(background)))

  act <- predict_ligand_activities(
    geneset                  = geneset,
    background_expressed_genes = background,
    ligand_target_matrix     = lt,
    potential_ligands        = potential) %>%
    as_tibble() %>%
    arrange(desc(aupr), desc(auroc))

  act$zone <- zone_name; act$direction <- tag
  act$n_sender <- n_s; act$n_receiver <- n_r
  act$receiver_geneset_n <- length(geneset)

  ## top ligand -> predicted targets in the receiver
  top_lig <- act$test_ligand[1]
  tg <- data.frame(gene = geneset,
                   weight = as.numeric(lt[geneset, top_lig]),
                   r_with_sender = cr$r[match(geneset, cr$gene)]) %>%
    arrange(desc(weight))
  write.csv(tg, file.path(TDIR, sprintf("targets_%s_%s.csv", tag, gsub(" ", "", zone_name))),
            row.names = FALSE)

  list(act = act, geneset = geneset, top_ligand = top_lig,
       n_s = n_s, n_r = n_r, targets = tg)
}

## ---------------- run both directions ----------------
cat("[2] direction A: C1Q TAM (ligand) -> MMP11 CAF (receiver)\n")
resA <- lapply(ZONES, function(z)
  run_one(z, C1Q, MMP, "C1Q TAMs", "MMP11 CAFs", TH_C1Q, TH_MMP, "A_C1Q_to_MMP11"))
names(resA) <- ZONES

cat("[3] direction B: MMP11 CAF (ligand) -> C1Q TAM (receiver)\n")
resB <- lapply(ZONES, function(z)
  run_one(z, MMP, C1Q, "MMP11 CAFs", "C1Q TAMs", TH_MMP, TH_C1Q, "B_MMP11_to_C1Q"))
names(resB) <- ZONES

## ---------------- collect ----------------
collect <- function(res, tag) {
  bind_rows(lapply(names(res), function(z) if (!is.null(res[[z]])) res[[z]]$act))
}
actA <- collect(resA, "A_C1Q_to_MMP11")
actB <- collect(resB, "B_MMP11_to_C1Q")
write.csv(actA, file.path(TDIR, "ligand_activity_A_C1Q_to_MMP11.csv"), row.names = FALSE)
write.csv(actB, file.path(TDIR, "ligand_activity_B_MMP11_to_C1Q.csv"), row.names = FALSE)

## ---------------- zone comparison ----------------
zone_cmp <- function(act, top_n = 20) {
  zs <- unique(act$zone)
  top_by_zone <- lapply(zs, function(z)
    act %>% filter(zone == z) %>% arrange(desc(aupr)) %>% head(top_n) %>% pull(test_ligand))
  names(top_by_zone) <- zs
  shared <- Reduce(intersect, top_by_zone)
  ## rank correlation on ligands measured in every zone
  wide <- act %>% select(test_ligand, zone, aupr) %>%
    pivot_wider(names_from = zone, values_from = aupr) %>%
    filter(complete.cases(.))
  cors <- if (nrow(wide) >= 5 && ncol(wide) >= 3) {
    m <- as.matrix(wide[, zs]); cor(m, method = "spearman")
  } else NULL
  list(top_by_zone = top_by_zone, shared = shared, wide = wide, cors = cors)
}
cmpA <- zone_cmp(actA); cmpB <- zone_cmp(actB)

saveRDS(list(A = resA, B = resB, cmpA = cmpA, cmpB = cmpB,
             spot_counts = spot_counts),
        file.path(OUT, "nichenet_results.rds"))

cat(sprintf("[4] shared top-20 ligands across zones: A=%d, B=%d\n",
            length(cmpA$shared), length(cmpB$shared)))
if (!is.null(cmpA$cors)) { cat("A rank cor:\n"); print(round(cmpA$cors, 2)) }
if (!is.null(cmpB$cors)) { cat("B rank cor:\n"); print(round(cmpB$cors, 2)) }

## ---------------- figures ----------------
plot_activity <- function(act, ttl, sub, n_show = 18) {
  if (nrow(act) == 0) return(NULL)
  top <- act %>% group_by(zone) %>% arrange(desc(aupr)) %>% slice_head(n = n_show) %>%
    ungroup() %>% distinct(test_ligand) %>% pull(test_ligand)
  d <- act %>% filter(test_ligand %in% top) %>%
    mutate(test_ligand = factor(test_ligand, levels = rev(top)),
           zone = factor(zone, levels = ZONES))
  ggplot(d, aes(zone, test_ligand, fill = aupr)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradientn(colours = c("#F7F7F7", "#FDDBC7", "#E8789A", "#C51B7D"),
                         name = "AUPR", na.value = "grey90") +
    labs(title = ttl, subtitle = sub, x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          axis.text.y = element_text(size = 5.8))
}
pA <- plot_activity(actA, "C1Q TAM -> MMP11 CAF", "spatial NicheNet, per zone")
pB <- plot_activity(actB, "MMP11 CAF -> C1Q TAM", "spatial NicheNet, per zone")

if (!is.null(pA)) ggsave(file.path(FDIR, "fig_nichenet_activity_A.pdf"), pA,
                         width = 3.0, height = 4.4, device = cairo_pdf, family = FONT, bg = "white")
if (!is.null(pA)) ggsave(file.path(FDIR, "fig_nichenet_activity_A.png"), pA,
                         width = 3.0, height = 4.4, dpi = 600, device = ragg::agg_png, bg = "white")
if (!is.null(pB)) ggsave(file.path(FDIR, "fig_nichenet_activity_B.pdf"), pB,
                         width = 3.0, height = 4.4, device = cairo_pdf, family = FONT, bg = "white")
if (!is.null(pB)) ggsave(file.path(FDIR, "fig_nichenet_activity_B.png"), pB,
                         width = 3.0, height = 4.4, dpi = 600, device = ragg::agg_png, bg = "white")

save_cor <- function(cors, tag) {
  if (is.null(cors)) return(NULL)
  d <- as.data.frame(as.table(cors))
  names(d) <- c("zone1", "zone2", "rho")
  p <- ggplot(d, aes(zone1, zone2, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 2.4) +
    scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#C51B7D",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman") +
    labs(title = sprintf("Zone rank agreement (%s)", tag), x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(FDIR, sprintf("fig_zone_agreement_%s.pdf", tag)), p,
         width = 2.8, height = 2.6, device = cairo_pdf, family = FONT, bg = "white")
  ggsave(file.path(FDIR, sprintf("fig_zone_agreement_%s.png", tag)), p,
         width = 2.8, height = 2.6, dpi = 600, device = ragg::agg_png, bg = "white")
}
save_cor(cmpA$cors, "A"); save_cor(cmpB$cors, "B")

## ---------------- summary table ----------------
summ <- bind_rows(
  actA %>% group_by(zone) %>% slice_max(aupr, n = 10) %>%
    summarise(zone = first(zone), direction = "C1Q TAM -> MMP11 CAF",
              top_ligands = paste(test_ligand, collapse = ", "),
              best = test_ligand[which.max(aupr)],
              best_aupr = round(max(aupr), 3),
              n_sender = first(n_sender), n_receiver = first(n_receiver),
              geneset_n = first(receiver_geneset_n), .groups = "drop"),
  actB %>% group_by(zone) %>% slice_max(aupr, n = 10) %>%
    summarise(zone = first(zone), direction = "MMP11 CAF -> C1Q TAM",
              top_ligands = paste(test_ligand, collapse = ", "),
              best = test_ligand[which.max(aupr)],
              best_aupr = round(max(aupr), 3),
              n_sender = first(n_sender), n_receiver = first(n_receiver),
              geneset_n = first(receiver_geneset_n), .groups = "drop"))
write.csv(summ, file.path(TDIR, "summary_top_ligands_by_zone.csv"), row.names = FALSE)
print(as.data.frame(summ))

cat("\n[Done] NicheNet C1Q TAM <-> MMP11 CAF\n")
cat("  outputs:", OUT, "\n")
