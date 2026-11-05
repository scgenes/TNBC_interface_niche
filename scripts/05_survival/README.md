# 05 · Survival analysis (Figure S6)

Bulk survival validation of the C1Q⁺ TAM, SPP1⁺ TAM and MMP11⁺ CAF signatures in
triple-negative breast cancer.

## Inputs

| Input | Description | Location |
|---|---|---|
| Mac_sub / Fib_sub DE tables | Single-cell differential-expression tables (`genes, log2FC, cluster, pct, avgExpr`) used to build the 20-gene signatures | Not included; export from the single-cell project |
| `data/signatures/gene_list_top20_*.{csv,tsv}` | Delivered 20-gene lists (3 signatures, 59 unique genes) | Included |
| TCGA-BRCA | cBioPortal files: `data_clinical_patient.txt`, `data_clinical_sample.txt`, `data_mrna_seq_v2_rsem.txt` | `${TNBC_SURVIVAL_ROOT}/brca_tcga/` |
| METABRIC | cBioPortal files: `data_clinical_patient.txt`, `data_clinical_sample.txt`, `data_mrna_illumina_microarray.txt` | `${TNBC_SURVIVAL_ROOT}/brca_metabric/` |
| GSE58812 | `GSE58812_expr_genes.tsv.gz`, `GSE58812_clinical.csv` | `${TNBC_SURVIVAL_ROOT}/geo/GSE58812/` |
| GSE21653 | `GSE21653_expr_genes.tsv.gz`, `GSE21653_clinical.csv` | `${TNBC_SURVIVAL_ROOT}/geo/GSE21653/` |

Only **primary tumours** and **triple-negative** samples (ER⁻/PR⁻/HER2⁻; HER2 by
IHC/FISH) are retained. No subsampling, no imputation.

> **Scope of the manuscript figure.** Figure S6 presents only GSE21653 (DFS,
> discovery) and GSE58812 (OS, validation), scored with ssGSEA and split at the
> median. The TCGA-BRCA and METABRIC code paths are retained from the
> exploratory analysis but are not shown in the paper.

## Scripts and run order

| Script | Purpose |
|---|---|
| `01_build_gene_lists.py` | Rebuild the 20-gene lists from the two DE tables (top20 by log2FC after removing technical/sex/mitochondrial and IGKV/IGLV genes) |
| `02_survival_analysis.R` | Sample definition, ssGSEA scoring, median/maxstat cutpoints, KM/log-rank, uni-/multivariable Cox, PH checks, BH correction, tables and per-cohort figures |
| `03_combined_figure.R` | 2 × 4 discovery/validation KM figure (GSE21653 DFS discovery, GSE58812 OS validation) and 180 mm submission export (PDF/PNG/TIFF) |

```bash
# 1) gene lists (optional; the delivered lists are already in data/signatures/)
python 01_build_gene_lists.py \
  --mac /path/to/Mac_sub_差异列表.tsv \
  --fib /path/to/Fib_sub_差异列表.tsv \
  --out ../../data/signatures

# 2) analysis (run from the repository root)
export TNBC_SURVIVAL_ROOT=/path/to/project/project_spatial/TNBC/survival
Rscript scripts/05_survival/02_survival_analysis.R

# 3) combined figure + submission export
export PAPER_FIGURE_DIR=/path/to/project/project_spatial/paper/Figure/submission
Rscript scripts/05_survival/03_combined_figure.R
```

## Key parameters

- **Signatures**: 20 genes per cell state, ranked by `log2FC` from the Mac_sub /
  Fib_sub DE tables; excluded `^(MT-|RPL|RPS|HB[AB]|MALAT1$|XIST$|RPS4Y1$|DDX3Y$)`
  and `^(IGKV|IGLV)`.
- **Scoring**: ssGSEA within each TNBC subset (`GSVA::gsva`, `ssgsea.norm=TRUE`,
  α=0.25); gene-wise z-score mean as sensitivity analysis.
- **Cutoffs**: median split (`High = score > median`) and maxstat LogRank optimal
  cutpoint (`maxstat.test`, `smethod="LogRank"`, `pmethod="Lau92"`,
  `minprop=0.25`, `maxprop=0.9`).
- **Dual-high**: C1Q and MMP11 both above their own cutpoint (same cutoff method);
  all other samples form the comparison group.
- **Models**: KM + log-rank; Cox univariable and multivariable (TCGA/METABRIC:
  age + stage; GSE58812: age; GSE21653: age + grade); `cox.zph`; BH correction
  within score-method × cutoff × cohort × endpoint (4 groups per family).
- **Software**: R 4.3.3, survival 3.8.3, survminer 0.5.0, maxstat 0.7.26,
  GSVA 1.50.0, limma, data.table, ggplot2, patchwork, ggsci.

## Outputs

- `processed/survival_v2/tables/survival_results_primary.csv` — primary results
  (96 rows: cohort × endpoint × score × cutoff × 4 groups).
- `processed/survival_v2/tables/stratification_cutpoints.csv` — cutpoints and
  group sizes.
- `processed/survival_v2/tables/combined_panel_source_2cohorts.csv` — the 8 rows
  used in Figure S6.
- `processed/survival_v2/figures/km/` — per-cohort KM panels.
- `processed/survival_v2/figures/forest/` — forest plots.
- `processed/survival_v2/figures/combined/` — the 2 × 4 discovery/validation
  figure and its 180 mm version.
- Submission files: `${PAPER_FIGURE_DIR}/FigureS6_survival_180mm.{pdf,png,tif}`.

Expected reference outputs from the original run are in
`reference_outputs/survival_v2/` (three tables) for comparison.

## Notes

- `02_survival_analysis.R` is a converted notebook script. It retains the
  in-notebook QC/verification cells from the original analysis; cells after the
  audit fix (roughly cell 29 onward) supersede earlier versions. Running the
  script top to bottom leaves the objects and tables in the delivered state.
- The full write-up of the survival analysis is in the manuscript Methods
  (§2.11) and Results (§3.8), and in the figure legend for Figure S6.
- `03_combined_figure.R` reproduces the manuscript figure: discovery =
  GSE21653 DFS (ssGSEA, median), validation = GSE58812 OS (ssGSEA, median),
  four groups (C1Q, dual-high, SPP1, MMP11).
