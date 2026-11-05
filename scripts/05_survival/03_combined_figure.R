# =============================================================================
# Combined discovery/validation KM figure (Figure S6) and submission export.
# Runs after 02_survival_analysis.R (reads tables from processed/survival_v2).
#
# Outputs:
#   processed/survival_v2/figures/combined/fig_km_discovery_validation_2cohorts_4groups.*
#   ${PAPER_FIGURE_DIR}/FigureS6_survival_180mm.{pdf,png,tif}
#
# Set TNBC_SURVIVAL_ROOT and PAPER_FIGURE_DIR (see 02_survival_analysis.R).
# =============================================================================

# 新增 2×4 组合图：Discovery GSE21653 DFS / Validation GSE58812 OS（不覆盖旧图）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
suppressPackageStartupMessages({library(survival); library(ggplot2); library(patchwork); library(ggsci)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival"); OUT <- file.path(ROOT,"processed/survival_v2")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
grp  <- read.csv(file.path(OUT,"tables/group_assignment_long.csv"), stringsAsFactors=FALSE)
prim <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)
pal2 <- unname(ggsci::pal_lancet()(2)); names(pal2) <- c("Low","High")
combos <- data.frame(cohort=c("GSE21653","GSE58812"), endpoint=c("DFS","OS"),
                     score_method=c("ssgsea","ssgsea"), cutoff_method=c("median","median"),
                     stringsAsFactors=FALSE)
groups <- c("C1Q_TAM","C1Q_MMP11_DualHigh_vs_Other","SPP1_TAM","MMP11_CAF")
gcol <- c(C1Q_TAM="C1Q_grp", SPP1_TAM="SPP1_grp", MMP11_CAF="MMP11_grp")
src <- list(); k <- 0
for (i in seq_len(nrow(combos))) for (g in groups) {
  r <- prim[prim$cohort==combos$cohort[i] & prim$endpoint==combos$endpoint[i] &
            prim$score_method==combos$score_method[i] & prim$cutoff_method==combos$cutoff_method[i] &
            prim$analysis_group==g, ]
  stopifnot(nrow(r)==1); k <- k+1
  src[[k]] <- data.frame(cohort=r$cohort, endpoint=r$endpoint, score_method=r$score_method,
    cutoff_method=r$cutoff_method, analysis_group=r$analysis_group, n_high=r$n_high, n_low=r$n_low,
    events_high=r$events_high, events_low=r$events_low, HR=r$HR, HR_lower=r$HR_lower, HR_upper=r$HR_upper,
    p_effective=r$p_effective, p_fdr=r$p_fdr, stringsAsFactors=FALSE)
}
src <- do.call(rbind, src); write.csv(src, file.path(OUT,"tables/combined_panel_source_2cohorts.csv"), row.names=FALSE)
cat("[source] rows:", nrow(src), "| pairs:", paste(unique(paste(src$cohort, src$endpoint)), collapse=" -> "), "\n")

panel <- function(co, ep, meth, cm, g, bottom=FALSE) {
  d0 <- long[long$cohort==co & long$endpoint==ep, ]
  g0 <- grp[grp$cohort==co & grp$endpoint==ep & grp$score_method==meth & grp$cutoff_method==cm, ]
  d <- merge(d0[,c("sample","time","event")], g0, by="sample")
  d <- d[!is.na(d$time) & !is.na(d$event) & d$time > 0, ]
  if (g == "C1Q_MMP11_DualHigh_vs_Other") {
    dd <- d; dd$grp <- factor(ifelse(dd$dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))
  } else { dd <- d; dd$grp <- factor(dd[[gcol[[g]]]], levels=c("Low","High")) }
  dd <- dd[!is.na(dd$grp), ]
  r <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method==meth & prim$cutoff_method==cm & prim$analysis_group==g, ]
  lab <- if (g=="C1Q_MMP11_DualHigh_vs_Other") "C1Q+MMP11+ vs Other" else
         switch(g, C1Q_TAM="C1Q TAMs", SPP1_TAM="SPP1 TAMs", MMP11_CAF="MMP11 CAFs")
  qstar <- if (nrow(r) && !is.na(r$p_fdr[1]) && r$p_fdr[1] < 0.05) "*" else ""
  sub <- if (!nrow(r) || is.na(r$HR[1])) "not fitted" else
    sprintf("HR=%.2f (95%%CI %.2f-%.2f)\np_eff=%.3g; q=%.3g%s", r$HR[1], r$HR_lower[1], r$HR_upper[1], r$p_effective[1], r$p_fdr[1], qstar)
  fit <- survfit(Surv(time,event) ~ grp, data=dd)
  sname <- rep(sub("grp=","",names(fit$strata)), fit$strata)
  km <- data.frame(time=fit$time, surv=fit$surv, grp=sname, cens=fit$n.censor)
  km <- rbind(data.frame(time=0, surv=1, grp=c("Low","High"), cens=0), km)
  km <- km[order(km$grp, km$time), ]
  cen <- km[km$cens > 0, ]
  tmax <- max(dd$time, na.rm=TRUE)
  brk <- pretty(c(0, tmax), n=4); brk <- sort(unique(brk[brk >= 0 & brk <= tmax]))
  if (length(brk) < 4) brk <- sort(unique(c(brk, tmax)))
  tp <- brk
  rk <- expand.grid(time=tp, grp=c("Low","High"), KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  rk$n <- mapply(function(t, g2) sum(dd$time >= t & as.character(dd$grp) == g2), rk$time, rk$grp)
  rk$grp <- factor(rk$grp, levels=c("Low","High"))
  xpad <- expansion(mult=c(0.05,0.10))
  base_th <- theme_bw() + theme(panel.grid=element_blank(), panel.border=element_rect(colour="black", fill=NA),
    plot.title=element_text(hjust=0.5, face="bold", size=8.5, lineheight=0.95),
    plot.subtitle=element_text(hjust=0.5, size=7, lineheight=0.95),
    axis.title=element_text(size=8), axis.text=element_text(size=7), legend.position="none",
    plot.margin=margin(3,4,1,3))
  pkm <- ggplot(km, aes(x=time, y=surv, color=grp)) +
    geom_step(linewidth=0.9, direction="hv") +
    geom_point(data=cen, aes(x=time, y=surv, color=grp), shape=3, size=0.9, stroke=0.35, show.legend=FALSE) +
    scale_color_manual(values=pal2, limits=c("Low","High")) +
    scale_x_continuous(limits=c(0, tmax), breaks=tp, expand=xpad) +
    scale_y_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1)) +
    labs(title=sprintf("%s %s\n%s", co, ep, lab), subtitle=sub, x=NULL, y="Survival probability") +
    base_th + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  prt <- ggplot(rk, aes(x=time, y=grp, label=format(n, big.mark=","))) +
    geom_text(size=6.5/.pt) +
    scale_x_continuous(limits=c(0, tmax), breaks=tp, expand=xpad) +
    labs(x=if (bottom) "Time (months)" else NULL, y=NULL) +
    theme_bw() + theme(panel.grid=element_blank(), panel.border=element_rect(colour="black", fill=NA),
      axis.text.y=element_text(size=6.5), axis.text.x=element_text(size=7), axis.title=element_text(size=8),
      axis.ticks.y=element_blank(), plot.margin=margin(0,4,3,3),
      axis.title.x=if (bottom) element_text(size=8) else element_blank())
  pkm / prt + plot_layout(heights=c(3.1,1))
}
P <- list(); q <- 0
for (i in seq_len(nrow(combos))) for (g in groups) { q <- q+1
  P[[q]] <- panel(combos$cohort[i], combos$endpoint[i], combos$score_method[i], combos$cutoff_method[i], g, bottom=(i==2)) }
strip <- function(txt) ggplot(data.frame(x=1,y=1,lab=txt), aes(x,y,label=lab)) +
  geom_text(angle=90, size=9/.pt) + coord_cartesian(clip="off") + theme_void() + theme(plot.margin=margin(2,0,2,0))
sA <- strip("Discovery (GSE21653)"); sB <- strip("Validation (GSE58812)")
pleg <- ggplot(data.frame(x=c(0,1,0,1), y=c(0,1,0,1),
                          grp=factor(rep(c("Low","High"), each=2), levels=c("Low","High"))), aes(x,y,color=grp)) +
  geom_line(linewidth=0.9, alpha=0) + scale_color_manual(values=pal2, limits=c("Low","High")) +
  theme_void() + theme(legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=8)) +
  guides(color=guide_legend(nrow=1, override.aes=list(linewidth=0.9, alpha=1), title=NULL))
row1 <- P[[1]] | P[[2]] | P[[3]] | P[[4]]
row2 <- P[[5]] | P[[6]] | P[[7]] | P[[8]]
right <- row1 / row2
left <- sA / sB + plot_layout(heights=c(1,1))
fig <- ((left | right) + plot_layout(widths=c(0.045,0.955))) / pleg + plot_layout(heights=c(1,0.07))
stem <- file.path(OUT,"figures/combined","fig_km_discovery_validation_2cohorts_4groups")
cairo_pdf(paste0(stem,".pdf"), width=13, height=7.6, bg="white"); print(fig); dev.off()
png(paste0(stem,".png"), width=13, height=7.6, units="in", res=600, bg="white", type="cairo"); print(fig); dev.off()
cat("[saved]", basename(stem), "| pages:", qpdf::pdf_length(paste0(stem,".pdf")),
    "| png KB:", round(file.info(paste0(stem,".png"))$size/1024), "\n")
cat("===== Cell End:", date(), "=====\n")


# 180×135 mm 投稿图（终版：短组名防裁切、面板字母上移至标题带左侧）
cat("===== Cell Start:", date(), "=====\n")
.start <- Sys.time()
suppressPackageStartupMessages({library(survival); library(ggplot2); library(patchwork); library(ggsci); library(magick)})
ROOT <- Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival"); OUT <- file.path(ROOT,"processed/survival_v2")
SUB <- Sys.getenv("PAPER_FIGURE_DIR", "/path/to/project/project_spatial/paper/Figure/submission")
long <- read.csv(file.path(OUT,"tables/survival_inputs_long.csv"), stringsAsFactors=FALSE)
grp  <- read.csv(file.path(OUT,"tables/group_assignment_long.csv"), stringsAsFactors=FALSE)
prim <- read.csv(file.path(OUT,"tables/survival_results_primary.csv"), stringsAsFactors=FALSE)
pal2 <- unname(ggsci::pal_lancet()(2)); names(pal2) <- c("Low","High")
combos <- data.frame(cohort=c("GSE21653","GSE58812"), endpoint=c("DFS","OS"), score_method=c("ssgsea","ssgsea"),
                     cutoff_method=c("median","median"), stringsAsFactors=FALSE)
groups <- c("C1Q_TAM","C1Q_MMP11_DualHigh_vs_Other","SPP1_TAM","MMP11_CAF")
gcol <- c(C1Q_TAM="C1Q_grp", SPP1_TAM="SPP1_grp", MMP11_CAF="MMP11_grp")
panel <- function(co, ep, meth, cm, g, tag, bottom=FALSE) {
  d0 <- long[long$cohort==co & long$endpoint==ep, ]; g0 <- grp[grp$cohort==co & grp$endpoint==ep & grp$score_method==meth & grp$cutoff_method==cm, ]
  d <- merge(d0[,c("sample","time","event")], g0, by="sample"); d <- d[!is.na(d$time)&!is.na(d$event)&d$time>0,]
  if (g == "C1Q_MMP11_DualHigh_vs_Other") { dd <- d; dd$grp <- factor(ifelse(dd$dual_grp=="DualHigh","High","Low"), levels=c("Low","High"))
  } else { dd <- d; dd$grp <- factor(dd[[gcol[[g]]]], levels=c("Low","High")) }
  dd <- dd[!is.na(dd$grp), ]
  r <- prim[prim$cohort==co & prim$endpoint==ep & prim$score_method==meth & prim$cutoff_method==cm & prim$analysis_group==g, ]
  lab <- if (g=="C1Q_MMP11_DualHigh_vs_Other") "Dual C1Q+MMP11+" else switch(g, C1Q_TAM="C1Q TAMs", SPP1_TAM="SPP1 TAMs", MMP11_CAF="MMP11 CAFs")
  qstar <- if (nrow(r) && !is.na(r$p_fdr[1]) && r$p_fdr[1] < 0.05) "*" else ""
  sub <- sprintf("HR=%.2f (95%%CI %.2f-%.2f)\nq=%.3g%s", r$HR[1], r$HR_lower[1], r$HR_upper[1], r$p_fdr[1], qstar)
  fit <- survfit(Surv(time,event) ~ grp, data=dd); sn <- rep(sub("grp=","",names(fit$strata)), fit$strata)
  km <- data.frame(time=fit$time, surv=fit$surv, grp=sn, cens=fit$n.censor)
  km <- rbind(data.frame(time=0, surv=1, grp=c("Low","High"), cens=0), km); km <- km[order(km$grp, km$time),]
  cen <- km[km$cens>0,]; tmax <- max(dd$time, na.rm=TRUE)
  brk <- pretty(c(0,tmax), n=4); brk <- sort(unique(brk[brk>=0 & brk<=tmax])); if (length(brk)<4) brk <- sort(unique(c(brk,tmax)))
  rk <- expand.grid(time=brk, grp=c("Low","High"), KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  rk$n <- mapply(function(t,g2) sum(dd$time>=t & as.character(dd$grp)==g2), rk$time, rk$grp)
  rk$grp <- factor(rk$grp, levels=c("Low","High")); xpad <- expansion(mult=c(0.05,0.10))
  th <- theme_bw() + theme(panel.grid=element_blank(), panel.border=element_rect(colour="black", fill=NA),
    plot.title=element_text(hjust=0.5, face="bold", size=7.5, lineheight=0.92),
    plot.subtitle=element_text(hjust=0.5, size=6.8, lineheight=0.92),
    axis.title=element_text(size=7), axis.text=element_text(size=6.8), legend.position="none", plot.margin=margin(2,0,1,0))
  pkm <- ggplot(km, aes(x=time, y=surv, color=grp)) + geom_step(linewidth=0.8, direction="hv") +
    geom_point(data=cen, aes(x=time,y=surv,color=grp), shape=3, size=0.7, stroke=0.3, show.legend=FALSE) +
    scale_color_manual(values=pal2, limits=c("Low","High")) +
    scale_x_continuous(limits=c(0,tmax), breaks=brk, expand=xpad) +
    scale_y_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1)) +
    labs(title=sprintf("%s %s\n%s", co, ep, lab), subtitle=sub, x=NULL, y="Survival probability") +
    th + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    annotate("text", x=-Inf, y=Inf, label=tag, hjust=-0.25, vjust=4.2, size=9/.pt, fontface="bold") + coord_cartesian(clip="off")
  prt <- ggplot(rk, aes(x=time, y=grp, label=format(n, big.mark=","))) + geom_text(size=6.8/.pt) +
    scale_x_continuous(limits=c(0,tmax), breaks=brk, expand=xpad) + labs(x=if (bottom) "Time (months)" else NULL, y=NULL) +
    theme_bw() + theme(panel.grid=element_blank(), panel.border=element_rect(colour="black", fill=NA),
      axis.text.y=element_text(size=6.8), axis.text.x=element_text(size=6.8), axis.title=element_text(size=7),
      axis.ticks.y=element_blank(), plot.margin=margin(0,0,2,0), axis.title.x=if (bottom) element_text(size=7) else element_blank())
  pkm / prt + plot_layout(heights=c(3.2,1))
}
P <- list(); q <- 0; letters8 <- letters[1:8]
for (i in 1:2) for (j in 1:4) { q <- q+1
  P[[q]] <- panel(combos$cohort[i], combos$endpoint[i], combos$score_method[i], combos$cutoff_method[i], groups[j], letters8[q], bottom=(i==2)) }
strip <- function(txt) ggplot(data.frame(x=1,y=1,lab=txt), aes(x,y,label=lab)) + geom_text(angle=90, size=7/.pt) +
  coord_cartesian(clip="off") + theme_void() + theme(plot.margin=margin(1,0,1,0))
sA <- strip("Discovery (GSE21653)"); sB <- strip("Validation (GSE58812)")
pleg <- ggplot(data.frame(x=c(0,1,0,1), y=c(0,1,0,1), grp=factor(rep(c("Low","High"), each=2), levels=c("Low","High"))), aes(x,y,color=grp)) +
  geom_line(linewidth=0.8, alpha=0) + scale_color_manual(values=pal2, limits=c("Low","High")) + theme_void() +
  theme(legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=7.5)) +
  guides(color=guide_legend(nrow=1, override.aes=list(linewidth=0.8, alpha=1), title=NULL))
right <- (P[[1]]|P[[2]]|P[[3]]|P[[4]]) / (P[[5]]|P[[6]]|P[[7]]|P[[8]]); left <- sA / sB + plot_layout(heights=c(1,1))
fig <- ((left | right) + plot_layout(widths=c(0.05,0.95))) / pleg + plot_layout(heights=c(1,0.075))
W <- 180/25.4; H <- 135/25.4; stem1 <- file.path(SUB,"FigureS6_survival_180mm")
cairo_pdf(paste0(stem1,".pdf"), width=W, height=H, bg="white"); print(fig); dev.off()
png(paste0(stem1,".png"), width=(4252+0.2)/600, height=(3189+0.2)/600, units="in", res=600, bg="white", type="cairo"); print(fig); dev.off()
py <- sprintf("import fitz, os\np=r'%s.pdf'\nd=fitz.open(p); pg=d[0]; pg.set_mediabox(fitz.Rect(0,0,510.2362,382.6772))\nt=p+'.tmp'; d.save(t, garbage=4, deflate=True); d.close(); os.replace(t,p)\nprint('mediabox set')", stem1)
invisible(suppressWarnings(system2("python3", c("-c", shQuote(py)), stdout=TRUE, stderr=TRUE)))
stem2 <- file.path(OUT,"figures/combined","fig_km_discovery_validation_2cohorts_4groups_180mm")
file.copy(paste0(stem1,".pdf"), paste0(stem2,".pdf"), overwrite=TRUE); file.copy(paste0(stem1,".png"), paste0(stem2,".png"), overwrite=TRUE)
cat("[saved] pages:", qpdf::pdf_length(paste0(stem1,".pdf")), "| png:", paste(magick::image_info(magick::image_read(paste0(stem1,".png")))[,c("width","height")], collapse="x"), "\n")
cat("===== Cell End:", date(), "=====\n")


# 生成 LZW TIFF 并做全项校验
cat("===== Cell Start:", date(), "=====\n")
suppressPackageStartupMessages({library(magick)})
OUT <- file.path(Sys.getenv("TNBC_SURVIVAL_ROOT", "/path/to/project/project_spatial/TNBC/survival"), "processed/survival_v2")
SUB <- Sys.getenv("PAPER_FIGURE_DIR", "/path/to/project/project_spatial/paper/Figure/submission")
pdf <- file.path(SUB,"FigureS6_survival_180mm.pdf"); png <- file.path(SUB,"FigureS6_survival_180mm.png"); tif <- file.path(SUB,"FigureS6_survival_180mm.tif")
info <- magick::image_info(magick::image_read(png))
cat(sprintf("[png] %dx%d ; %.0f KB\n", info$width, info$height, file.info(png)$size/1024))
py <- sprintf("from PIL import Image\nim=Image.open(r'%s').convert('RGB')\nim.save(r'%s', format='TIFF', compression='tiff_lzw', dpi=(600,600))\nprint('saved', im.size, im.mode)", png, tif)
cat(paste(suppressWarnings(system2("python3", c("-c", shQuote(py)), stdout=TRUE, stderr=TRUE)), collapse="\n"), "\n")
py2 <- sprintf("from PIL import Image\nim=Image.open(r'%s')\nprint('TIFF', im.size, im.mode, im.info.get('compression'), im.info.get('dpi'))", tif)
cat(paste(suppressWarnings(system2("python3", c("-c", shQuote(py2)), stdout=TRUE, stderr=TRUE)), collapse="\n"), "\n")
py3 <- sprintf("import fitz\nd=fitz.open(r'%s')\np=d[0]\nprint('PAGE', round(p.rect.width,2), round(p.rect.height,2), 'pages', d.page_count)\nfonts=p.get_fonts(full=True)\nprint('FONTS', len(fonts), [(f[3], f[1]) for f in fonts])\nt=p.get_text()\nprint('NONASCII', sum(1 for ch in t if ord(ch)>127))\nsizes=[]\nfor b in p.get_text('dict')['blocks']:\n    for l in b.get('lines',[]):\n        for s in l.get('spans',[]): sizes.append(round(s['size'],2))\nprint('MINSIZE', min(sizes), 'MAXSIZE', max(sizes))\nprint('SIZES', sorted(set(sizes)))", pdf)
cat(paste(suppressWarnings(system2("python3", c("-c", shQuote(py3)), stdout=TRUE, stderr=TRUE)), collapse="\n"), "\n")
magick::image_write(magick::image_resize(magick::image_read(png), "1700x1276"), file.path(OUT,"review","FigureS6_180mm_preview.png"))
cat("[preview] saved\n")
cat("===== Cell End:", date(), "=====\n")
