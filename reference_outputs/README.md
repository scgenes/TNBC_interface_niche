# Reference outputs

Expected outputs from the original analysis run, provided for comparison and
regression checking. They are **not** inputs to the pipeline.

| File | Produced by | Content |
|---|---|---|
| `survival_v2/survival_results_primary.csv` | `scripts/05_survival/02_survival_analysis.R` | 96 rows: cohort × endpoint × scoring method × cutoff method × 4 analysis groups (HR, CI, p_eff, p_adj, q, PH) |
| `survival_v2/stratification_cutpoints.csv` | `scripts/05_survival/02_survival_analysis.R` | median/maxstat cutpoints, selection p, group sizes |
| `survival_v2/combined_panel_source_2cohorts.csv` | `scripts/05_survival/03_combined_figure.R` | 8 rows used in Figure S6 (GSE21653 DFS discovery, GSE58812 OS validation) |

To regenerate:

```bash
export TNBC_SURVIVAL_ROOT=/path/to/project/project_spatial/TNBC/survival
Rscript scripts/05_survival/02_survival_analysis.R
Rscript scripts/05_survival/03_combined_figure.R
```
