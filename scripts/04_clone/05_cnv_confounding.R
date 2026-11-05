# =============================================================================
# TNBC interface-niche code repository
# Original source: project_spatial/TNBC/patient_ALL/CNV_confounding_preview/plot_cnv_confounding.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

#!/usr/bin/env Rscript
## CNV confounding preview: P13 / P15
## a/b: clone comparison within matched CNV bins (Interface, faceted by slice)
## c/d: within-clone Core vs Interface comparison (paired by slice x clone)
## Read-only inputs; output to this folder.

suppressPackageStartupMessages({
  library(qs); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

ROOT <- "/path/to/project/project_spatial"
OUT  <- file.path(ROOT, "TNBC/patient_ALL/CNV_confounding_preview")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

g <- qread(file.path(ROOT, "TNBC/TNBC_analysis/SpatialInferCNV/results/step3_spatial_tme_groups/data/spatial_tme_groups.qs"))
g$bc <- sub("_[0-9]+$", "", g$barcode)

wilson <- function(x, n, z = 1.96) {
  if (n == 0) return(c(NA, NA))
  p <- x / n; den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(max(0, ctr - half), min(1, ctr + half))
}

summ <- function(d, grp, value) {
  d %>% group_by(across(all_of(grp))) %>%
    summarise(n = n(), x = sum(.data[[value]]), frac = mean(.data[[value]]), .groups = "drop") %>%
    rowwise() %>%
    mutate(lo = wilson(x, n)[1], hi = wilson(x, n)[2]) %>%
    ungroup()
}

plot_bin <- function(d, value, ylab) {
  s <- summ(d, c("slice", "bin_mid", "clone"), value) %>% filter(n >= 10)
  ggplot(s, aes(bin_mid, frac, colour = clone, fill = clone)) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.6) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.004, linewidth = 0.3) +
    facet_wrap(~ slice, nrow = 1) +
    scale_colour_manual(values = c(Clone1 = "#E64B35", Clone2 = "#3C5488"), name = NULL) +
    scale_fill_manual(values = c(Clone1 = "#E64B35", Clone2 = "#3C5488"), guide = "none") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(x = "CNV score (0.02 bins)", y = ylab) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 7.5),
          legend.position = "none")
}

plot_region <- function(d, value, ylab) {
  s <- summ(d, c("slice", "clone", "region2"), value)
  w <- s %>% select(slice, clone, region2, n) %>%
    pivot_wider(names_from = region2, values_from = n)
  if (!all(c("Core", "Interface") %in% names(w))) return(NULL)
  keep <- w %>% filter(!is.na(Core), !is.na(Interface), Core >= 10, Interface >= 10)
  s <- s %>% semi_join(keep, by = c("slice", "clone")) %>%
    mutate(grp = interaction(slice, clone, sep = " | "))
  ggplot(s, aes(region2, frac, colour = clone, linetype = slice, shape = slice, group = grp)) +
    geom_line(linewidth = 0.45, alpha = 0.9) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = c(Clone1 = "#E64B35", Clone2 = "#3C5488"), name = NULL) +
    scale_linetype_discrete(name = "Slice") +
    scale_shape_discrete(name = "Slice") +
    scale_x_discrete(limits = c("Core", "Interface")) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(),
          legend.position = "bottom", legend.key.size = unit(2.5, "mm"))
}

for (p in c("P13", "P15")) {
  cl <- read.csv(file.path(ROOT, sprintf("TNBC/patient_ALL/clone/%s/clone_membership.csv", p)),
                 stringsAsFactors = FALSE)
  cl$bc <- sub("^[^_]+_", "", cl$cell)
  m <- merge(g, cl[, c("slice", "bc", "clone", "cnv_burden")], by = c("slice", "bc"))

  ## TAM-high threshold per slice (all spots), as in the clone-niche analysis
  m <- m %>% group_by(slice) %>%
    mutate(c1q_hi = C1QC_TAM >= quantile(C1QC_TAM, 0.75),
           spp1_hi = SPP1_TAM >= quantile(SPP1_TAM, 0.75)) %>% ungroup()

  m_if <- m %>% filter(region == "Interface") %>%
    mutate(bin_mid = floor(cnv_burden / 0.02) * 0.02 + 0.01)
  m_rg <- m %>% filter(region %in% c("Tumor Core", "Interface")) %>%
    mutate(region2 = ifelse(region == "Tumor Core", "Core", "Interface"))

  pa <- plot_bin(m_if, "c1q_hi", "C1Q TAM-high fraction (top 25%)")
  pb <- plot_bin(m_if, "spp1_hi", "SPP1 TAM-high fraction (top 25%)")
  pc <- plot_region(m_rg, "c1q_hi", "C1Q TAM-high fraction (top 25%)")
  pd <- plot_region(m_rg, "spp1_hi", "SPP1 TAM-high fraction (top 25%)")

  fig <- (pa | pb) / (pc | pd) +
    plot_annotation(title = p, tag_levels = "a") &
    theme(legend.position = "bottom",
          plot.title = element_text(size = 9, face = "bold"),
          plot.tag = element_text(size = 9, face = "bold"))
  stem <- file.path(OUT, sprintf("%s_CNV_matched_clone_TAM", p))
  ggsave(paste0(stem, ".pdf"), fig, width = 7.2, height = 6.4, device = cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".png"), fig, width = 7.2, height = 6.4, dpi = 600, bg = "white")

  ## values
  v1 <- summ(m_if, c("slice", "bin_mid", "clone"), "c1q_hi") %>% mutate(set = "CNV_bin_C1Q")
  v2 <- summ(m_if, c("slice", "bin_mid", "clone"), "spp1_hi") %>% mutate(set = "CNV_bin_SPP1")
  v3 <- summ(m_rg, c("slice", "clone", "region2"), "c1q_hi") %>% mutate(set = "Region_C1Q")
  v4 <- summ(m_rg, c("slice", "clone", "region2"), "spp1_hi") %>% mutate(set = "Region_SPP1")
  bind_rows(v1, v2, v3, v4) %>% mutate(patient = p) %>%
    write.csv(file.path(OUT, sprintf("%s_values.csv", p)), row.names = FALSE)
  cat("[Done]", p, "\n")
}
