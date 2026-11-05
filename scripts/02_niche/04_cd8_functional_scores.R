# =============================================================================
# TNBC interface-niche code repository
# Original source: FigureAI/Figure3/scripts/24_Tcell_function_CD8.R
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

## T-cell gene-set scores restricted to CD8 T-cell spots, by C1QC/MMP11 niche.
suppressPackageStartupMessages({library(qs); library(Matrix); library(dplyr); library(tidyr); library(ggplot2)})
OUT <- "/path/to/project/project_spatial/FigureAI/Figure3"
d <- qread(file.path(OUT,"data/spot_table.qs"))
m <- readRDS("/path/to/project/project_singlecell/TNBC_result/final_figures/V2/spatial/spatial_sct_matrix.rds")
sets <- list(Cytotoxicity=c("GZMB","PRF1","NKG7","GNLY"),
             Exhaustion=c("PDCD1","LAG3","TIGIT","HAVCR2","TOX","CXCL13"),
             Treg=c("FOXP3","IL2RA","CTLA4","IKZF2"),
             TCR_activation=c("IL7R","CD69","NR4A1"))

d$cat <- with(d, ifelse(C1QC_pos & MMP11_pos, "Double+",
                 ifelse(C1QC_pos, "C1QC only",
                 ifelse(MMP11_pos, "MMP11 only", "Other"))))
d$CD8pos <- FALSE
for (s in unique(d$slice)) { i<-which(d$slice==s); d$CD8pos[i] <- d$CD8_sum[i] >= quantile(d$CD8_sum[i],0.75) }
cat("[CD8 spots by category]\n"); print(table(d$cat[d$CD8pos]))
cat("[CD8 spots by slice x category]\n"); print(table(d$slice[d$CD8pos], d$cat[d$CD8pos]))

sc <- list()
for (s in unique(d$slice)) {
  ds <- d[d$slice==s & d$CD8pos, ]
  X <- as.matrix(m[, ds$barcode, drop=FALSE]); Z <- t(scale(t(X))); Z[is.na(Z)] <- 0
  for (nm in names(sets)) {
    g <- intersect(sets[[nm]], rownames(Z))
    sc[[length(sc)+1]] <- data.frame(barcode=ds$barcode, slice=s, cat=ds$cat,
                                     set=nm, score=colMeans(Z[g,,drop=FALSE]))
  }
}
sc <- bind_rows(sc)
write.csv(sc, file.path(OUT,"tables/fig3i_CD8_scores_per_spot.csv"), row.names=FALSE)
mean_df <- sc %>% group_by(slice,cat,set) %>% summarise(mean_score=mean(score), n=n(), .groups="drop")
write.csv(mean_df, file.path(OUT,"tables/fig3i_CD8_scores_per_slice.csv"), row.names=FALSE)

comp <- list(Double_vs_Other=c("Double+","Other"),
             C1QC_vs_Other=c("C1QC only","Other"),
             MMP11_vs_Other=c("MMP11 only","Other"))
tests <- list()
for (nm in names(sets)) for (cn in names(comp)) {
  w <- mean_df %>% filter(set==nm, cat %in% comp[[cn]]) %>%
    select(slice,cat,mean_score) %>% pivot_wider(names_from=cat, values_from=mean_score)
  a <- w[[comp[[cn]][1]]]; b <- w[[comp[[cn]][2]]]; ok <- complete.cases(a,b)
  if (sum(ok) < 4) next
  wt <- wilcox.test(a[ok], b[ok], paired=TRUE, exact=FALSE)
  tests[[length(tests)+1]] <- data.frame(set=nm, comparison=cn, n=sum(ok),
    mean_diff=mean(a[ok]-b[ok]), n_pos=sum(a[ok]>b[ok]), p=wt$p.value)
}
tests <- bind_rows(tests) %>% mutate(padj=p.adjust(p,"BH"))
write.csv(tests, file.path(OUT,"tables/fig3i_CD8_score_tests.csv"), row.names=FALSE)
cat("\n[mean score by category]\n")
print(as.data.frame(mean_df %>% group_by(set,cat) %>% summarise(mean=mean(mean_score), .groups="drop") %>%
  pivot_wider(names_from=cat, values_from=mean)), digits=3)
cat("\n[tests]\n"); print(as.data.frame(tests), digits=3)

## heatmap
agg <- mean_df %>% group_by(set,cat) %>% summarise(mean_score=mean(mean_score), .groups="drop")
setlev <- c("Cytotoxicity","Exhaustion","TCR_activation","Treg")
catlev <- c("Double+","C1QC only","MMP11 only","Other")
map <- c(Double_vs_Other="Double+", C1QC_vs_Other="C1QC only", MMP11_vs_Other="MMP11 only")
ann <- tests %>% mutate(cat=map[comparison],
  star=ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,"*",ifelse(p<0.05,".",""))))) %>%
  select(set,cat,star)
agg <- agg %>% left_join(ann, by=c("set","cat")); agg$star[is.na(agg$star)] <- ""
agg$txt <- sprintf("%.2f%s", agg$mean_score, agg$star)
p <- ggplot(agg, aes(cat, set, fill=mean_score)) +
  geom_tile(color="white", linewidth=0.8) + geom_text(aes(label=txt), size=3.4) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0,
                       limits=c(-0.30,0.30), oob=scales::squish, name="Mean score") +
  scale_x_discrete(limits=catlev) + scale_y_discrete(limits=rev(setlev)) +
  labs(title="CD8 T-cell functional scores by niche",
       subtitle="Mean z-score in CD8 T-cell spots; * BH<0.05; . nominal p<0.05",
       x=NULL, y=NULL) +
  theme_minimal(base_size=10) + theme(panel.grid=element_blank(),
    axis.text.x=element_text(size=10, face="bold"), axis.text.y=element_text(size=10))
ggsave(file.path(OUT,"panels/fig3i_CD8_function_heatmap.png"), p, width=6.6, height=3.8, dpi=600, bg="white")
ggsave(file.path(OUT,"panels/fig3i_CD8_function_heatmap.pdf"), p, width=6.6, height=3.8, bg="white")
cat("[Done] CD8-restricted scores\n")
