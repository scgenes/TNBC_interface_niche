# =============================================================================
# TNBC interface-niche code repository
# Original source: FigureAI/Figure3/scripts/29_Activation_Interface_CD69_GZMB_PRF1.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

## Interface-stratified CD8 activation score: CD69 + GZMB + PRF1.
## Same CD8-spot definition as the original analysis (top 25% CD8_sum per slice),
## then restricted to the Interface region. Each gene is z-scored across the
## selected Interface spots within slice; the score is the mean of the three z-scores.
## Outputs:
##   panels/fig3i_Activation_CD69_GZMB_PRF1_Interface_box.{png,pdf}
##   panels/fig3i_Activation_CD69_GZMB_PRF1_Interface_heatmap.{png,pdf}
##   tables/fig3i_Activation_Interface_per_slice.csv
##   tables/fig3i_Activation_Interface_per_gene.csv
##   tables/fig3i_Activation_Interface_tests.csv
##   tables/fig3i_Activation_Interface_per_spot.csv
suppressPackageStartupMessages({library(qs); library(Matrix); library(dplyr); library(tidyr); library(ggplot2)})

OUT <- "/path/to/project/project_spatial/FigureAI/Figure3"
d <- qread(file.path(OUT,"data/spot_table.qs"))
m <- readRDS("/path/to/project/project_singlecell/TNBC_result/final_figures/V2/spatial/spatial_sct_matrix.rds")

d$cat <- with(d, ifelse(C1QC_pos & MMP11_pos, "Double+",
                 ifelse(C1QC_pos, "C1QC only",
                 ifelse(MMP11_pos, "MMP11 only", "Other"))))

## original CD8-spot definition: top 25% CD8_sum per slice
d$CD8pos <- FALSE
for (s in unique(d$slice)) {
  i <- which(d$slice == s)
  d$CD8pos[i] <- d$CD8_sum[i] >= quantile(d$CD8_sum[i], 0.75)
}

genes <- c("CD69","GZMB","PRF1")
sc <- list()
for (s in unique(d$slice)) {
  ds <- d[d$slice == s & d$CD8pos & d$region == "Interface", ]
  if (nrow(ds) < 10) next
  X <- as.matrix(m[genes, ds$barcode, drop = FALSE])
  Z <- t(scale(t(X))); Z[is.na(Z)] <- 0
  sc[[length(sc)+1]] <- data.frame(
    barcode = ds$barcode, slice = s, cat = ds$cat,
    CD69 = Z["CD69",], GZMB = Z["GZMB",], PRF1 = Z["PRF1",])
}
sc <- bind_rows(sc)
sc$score <- (sc$CD69 + sc$GZMB + sc$PRF1) / 3
write.csv(sc, file.path(OUT,"tables/fig3i_Activation_Interface_per_spot.csv"), row.names = FALSE)

## per-slice means; keep slices with >=5 spots in all three niches
ps <- sc %>% group_by(slice, cat) %>%
  summarise(n = n(), CD69 = mean(CD69), GZMB = mean(GZMB),
            PRF1 = mean(PRF1), score = mean(score), .groups = "drop")
write.csv(ps, file.path(OUT,"tables/fig3i_Activation_Interface_per_slice.csv"), row.names = FALSE)

w   <- ps %>% select(slice, cat, score) %>% pivot_wider(names_from = cat, values_from = score)
cnt <- ps %>% select(slice, cat, n)     %>% pivot_wider(names_from = cat, values_from = n)
niches <- c("Double+","C1QC only","MMP11 only")
keep <- complete.cases(w[, niches]) &
        apply(cnt[, niches], 1, min) >= 5
w <- w[keep, ]; ps_keep <- ps[ps$slice %in% w$slice, ]

## per-gene summary (equal slice weighting)
pg <- ps_keep %>% group_by(cat) %>%
  summarise(n_slices = n(), CD69 = mean(CD69), GZMB = mean(GZMB),
            PRF1 = mean(PRF1), score = mean(score), .groups = "drop")
write.csv(pg, file.path(OUT,"tables/fig3i_Activation_Interface_per_gene.csv"), row.names = FALSE)
cat("\n[per-gene / composite, Interface]\n"); print(as.data.frame(pg), digits = 3)

## paired tests
prs <- list(c("Double+","C1QC only"), c("Double+","MMP11 only"), c("C1QC only","MMP11 only"))
tt <- bind_rows(lapply(prs, function(p) {
  a <- w[[p[1]]]; b <- w[[p[2]]]
  wt <- wilcox.test(a, b, paired = TRUE, exact = FALSE)
  data.frame(comparison = paste(p[1], "vs", p[2]),
             n = nrow(w), mean_diff = mean(a - b),
             n_pos = sum(a > b), p = wt$p.value)
})) %>% mutate(padj = p.adjust(p, "BH"))
write.csv(tt, file.path(OUT,"tables/fig3i_Activation_Interface_tests.csv"), row.names = FALSE)
cat("\n[paired tests]\n"); print(as.data.frame(tt), digits = 4)

## ---------------- box plot ----------------
long <- w %>% pivot_longer(all_of(niches), names_to = "cat", values_to = "score") %>%
  mutate(cat = factor(cat, levels = niches))
fills <- c("Double+" = "#B2182B", "C1QC only" = "#E6AB02", "MMP11 only" = "#1B9E77")

lab <- function(i) sprintf("p = %.3f", tt$p[i])
ymax <- max(long$score); ymin <- min(long$score); gap <- ymax - ymin

bracket <- function(x1, x2, y, label, tip = gap * 0.016) {
  list(annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.42),
       annotate("segment", x = x1, xend = x1, y = y - tip, yend = y, linewidth = 0.42),
       annotate("segment", x = x2, xend = x2, y = y - tip, yend = y, linewidth = 0.42),
       annotate("text", x = (x1 + x2) / 2, y = y + gap * 0.012, label = label,
                size = 2.85, vjust = 0))
}

p_box <- ggplot(long, aes(cat, score)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.35) +
  geom_boxplot(aes(fill = cat), width = 0.55, outlier.shape = NA,
               alpha = 0.32, color = "grey25", linewidth = 0.45) +
  geom_line(aes(group = slice), color = "grey55", linewidth = 0.42, alpha = 0.85) +
  geom_point(size = 1.7, color = "grey15", stroke = 0) +
  scale_fill_manual(values = fills, guide = "none") +
  scale_x_discrete(limits = niches, labels = c("Double+", "C1Q only", "MMP11 only")) +
  scale_y_continuous(limits = c(ymin - gap * 0.05, ymax + gap * 0.36),
                     breaks = pretty(c(ymin, ymax), n = 5)) +
  bracket(1, 2, ymax + gap * 0.08, lab(1)) +
  bracket(2, 3, ymax + gap * 0.17, lab(3)) +
  bracket(1, 3, ymax + gap * 0.26, lab(2)) +
  labs(x = NULL, y = "Activation score (mean z)",
       title = "CD8 T-cell activation in Interface") +
  theme_classic(base_size = 10) +
  theme(axis.ticks = element_line(colour = "grey20"),
        axis.text.x = element_text(face = "bold", size = 10.5, color = "black"),
        axis.text.y = element_text(size = 9.5, color = "black"),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 10.8))

ggsave(file.path(OUT,"panels/fig3i_Activation_CD69_GZMB_PRF1_Interface_box.png"), p_box,
       width = 5.6, height = 4.6, dpi = 600, bg = "white")
ggsave(file.path(OUT,"panels/fig3i_Activation_CD69_GZMB_PRF1_Interface_box.pdf"), p_box,
       width = 5.6, height = 4.6, bg = "white")

## ---------------- per-gene heatmap ----------------
hm <- pg %>% filter(cat %in% niches) %>%
  select(cat, CD69, GZMB, PRF1, score) %>%
  pivot_longer(-cat, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("CD69","GZMB","PRF1","score")),
         cat = factor(cat, levels = niches))
hm$metric_lab <- factor(hm$metric,
  levels = c("CD69","GZMB","PRF1","score"),
  labels = c("CD69","GZMB","PRF1","Composite"))

p_hm <- ggplot(hm, aes(cat, metric_lab, fill = value)) +
  geom_tile(color = "grey88", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-0.10, 0.10),
                       oob = scales::squish, name = "Mean z") +
  scale_x_discrete(limits = niches, labels = c("Double+", "C1Q only", "MMP11 only")) +
  scale_y_discrete(limits = rev(c("CD69","GZMB","PRF1","Composite"))) +
  labs(title = "Interface: activation genes by niche", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10))

ggsave(file.path(OUT,"panels/fig3i_Activation_CD69_GZMB_PRF1_Interface_heatmap.png"), p_hm,
       width = 4.8, height = 3.4, dpi = 600, bg = "white")
ggsave(file.path(OUT,"panels/fig3i_Activation_CD69_GZMB_PRF1_Interface_heatmap.pdf"), p_hm,
       width = 4.8, height = 3.4, bg = "white")

cat("\n[Done] Interface CD69+GZMB+PRF1 activation\n")
cat("Slices in paired plot:", nrow(w), "\n")
print(as.data.frame(w), digits = 3)
