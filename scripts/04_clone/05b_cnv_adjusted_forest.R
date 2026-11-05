# =============================================================================
# TNBC interface-niche code repository
# Original source: project_spatial/TNBC/patient_ALL/CNV_confounding_preview/plot_adjusted_model_forest.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

#!/usr/bin/env Rscript
## Adjusted-model forest plot: clone effect on TAM-high fraction, P13 / P15
## Models (all include slice as a fixed effect; only 2 slices per patient):
##   Base (+slice), +CNV, +Region, +CNV+Region
## Outcome: binary TAM-high (top 25% per slice); linear probability model.
## Coefficient shown: Clone2 - Clone1 (difference in TAM-high fraction), 95% CI.

suppressPackageStartupMessages({
  library(qs); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

ROOT <- "/path/to/project/project_spatial"
OUT  <- file.path(ROOT, "TNBC/patient_ALL/CNV_confounding_preview")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

g <- qread(file.path(ROOT, "TNBC/TNBC_analysis/SpatialInferCNV/results/step3_spatial_tme_groups/data/spatial_tme_groups.qs"))
g$bc <- sub("_[0-9]+$", "", g$barcode)

fit_one <- function(d, y) {
  d$clone <- factor(d$clone, levels = c("Clone1", "Clone2"))
  d$region2 <- factor(d$region2, levels = c("Core", "Interface"))
  d$slice <- factor(d$slice)
  forms <- list(
    "Base (+slice)"   = as.formula(paste(y, "~ clone + slice")),
    "+CNV"            = as.formula(paste(y, "~ clone + slice + cnv_burden")),
    "+Region"         = as.formula(paste(y, "~ clone + slice + region2")),
    "+CNV +Region"    = as.formula(paste(y, "~ clone + slice + cnv_burden + region2"))
  )
  do.call(rbind, lapply(names(forms), function(nm) {
    fit <- lm(forms[[nm]], data = d)
    cf <- summary(fit)$coefficients
    est <- cf["cloneClone2", "Estimate"]; se <- cf["cloneClone2", "Std. Error"]
    data.frame(model = nm, est = est, lo = est - 1.96 * se, hi = est + 1.96 * se,
               p = cf["cloneClone2", "Pr(>|t|)"], n = nrow(d))
  }))
}

plot_forest <- function(d, xlab) {
  d$model <- factor(d$model, levels = rev(c("Base (+slice)", "+CNV", "+Region", "+CNV +Region")))
  d$col <- ifelse(d$est > 0, "#E64B35", "#3C5488")
  ggplot(d, aes(est, model, colour = col)) +
    geom_vline(xintercept = 0, linetype = "22", colour = "grey40", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.12, linewidth = 0.4) +
    geom_point(size = 2.0) +
    geom_text(aes(label = sprintf("%+.3f", est)), vjust = -0.9, size = 2.4, colour = "black") +
    scale_colour_identity() +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 7.2),
          plot.margin = margin(2, 8, 2, 2))
}

for (p in c("P13", "P15")) {
  cl <- read.csv(file.path(ROOT, sprintf("TNBC/patient_ALL/clone/%s/clone_membership.csv", p)),
                 stringsAsFactors = FALSE)
  cl$bc <- sub("^[^_]+_", "", cl$cell)
  m <- merge(g, cl[, c("slice", "bc", "clone", "cnv_burden")], by = c("slice", "bc")) %>%
    filter(region %in% c("Tumor Core", "Interface")) %>%
    mutate(region2 = ifelse(region == "Tumor Core", "Core", "Interface")) %>%
    group_by(slice) %>%
    mutate(c1q_hi = as.numeric(C1QC_TAM >= quantile(C1QC_TAM, 0.75)),
           spp1_hi = as.numeric(SPP1_TAM >= quantile(SPP1_TAM, 0.75))) %>%
    ungroup()

  f1 <- fit_one(m, "c1q_hi"); f2 <- fit_one(m, "spp1_hi")
  pa <- plot_forest(f1, "Clone2 - Clone1 (C1Q TAM-high fraction)")
  pb <- plot_forest(f2, "Clone2 - Clone1 (SPP1 TAM-high fraction)")

  fig <- (pa | pb) + plot_annotation(title = p, tag_levels = "a") &
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.tag = element_text(size = 9, face = "bold"))

  stem <- file.path(OUT, sprintf("%s_adjusted_model_forest", p))
  ggsave(paste0(stem, ".pdf"), fig, width = 7.2, height = 3.2, device = cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".png"), fig, width = 7.2, height = 3.2, dpi = 600, bg = "white")

  bind_rows(mutate(f1, outcome = "C1Q_TAM"), mutate(f2, outcome = "SPP1_TAM")) %>%
    mutate(patient = p) %>% write.csv(file.path(OUT, sprintf("%s_adjusted_model_values.csv", p)), row.names = FALSE)
  cat("[Done]", p, "\n")
}
