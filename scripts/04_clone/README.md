# 04_clone

Clone inference and clone–TAM analyses.

| Script | Purpose |
|---|---|
| `00_spatial_infercnv_patient.R` | Patient-level spatial inferCNV (converted from the original R notebook): CNV matrix, 6-state HMM, CNV-based zone summaries and QC. |
| `01_clone_inference_per_sample.R` | Per-sample clone inference: arm-level CNV matrix, silhouette-based k selection (k = 2–10), cluster stability, arm consensus, Newick tree. |
| `02_clone_tam_colocalisation.R` | Clone–TAM spatial proximity with two complementary statistics: (A) radius-neighbourhood observed edges with 1,000 permutations, and (B) kNN Ro/e with 1,000 clone-label permutations. Stratified by core/interface. Reads `clone_membership.csv`, `region_assignment_extended.csv` and `spatial_tme_groups.qs`. |
| `03_clone_differential_expression.R` | Clone1 vs Clone2 differential expression in P13 and P15 (`FindMarkers`, Wilcoxon, `logfc.threshold = 0`, BH `padj < 0.05`). Reads `clone_seurat.rds`. |
| `04_clone_tam_cellchat.R` | Clone-level CellChat in the interface of P13 and P15: Clone1/Clone2 → C1Q_TAM / SPP1_TAM, with TAM-enriched spots defined by the within-section q75 of the corresponding proportion. |
| `05_cnv_confounding.R` | CNV confounding diagnostics and linear probability models (TAM-high ~ clone + slice, adding CNV score and region term by term). |
| `05b_cnv_adjusted_forest.R` | Forest plot of CNV/region-adjusted clone effects. |
| `06_clone_tam_figure.py` | Main and supplementary clone–TAM figure panels from the tables in `TAM_clone/tables/`. |

## Run order

1. `00_spatial_infercnv_patient.R` (or use the archived per-patient clone outputs)
2. `01_clone_inference_per_sample.R`
3. `02_clone_tam_colocalisation.R`
4. `03_clone_differential_expression.R`
5. `04_clone_tam_cellchat.R`
6. `05_cnv_confounding.R` / `05b_cnv_adjusted_forest.R`
7. `06_clone_tam_figure.py`

## Notes

- `00_`, `02_`, `03_`, `04_` and `02_nichenet_clone_LR.R` in `03_communication/` were converted from R notebooks with `jupyter nbconvert --to script`. They retain the original in-notebook QC/audit cells; these are not required for the scientific workflow.
- The per-patient tree-node clone script referenced in Methods §3.2 (`run_clone_patient.R`) was not present in the archived project. The per-patient clone outputs (`clone_membership.csv`, `clone_k_selection.csv`, `clone_arm_consensus.csv`, etc.) are inputs to the downstream scripts here.
- Clone labels are patient-specific and must not be compared across patients.

## Paths

`/path/to/project` and `/path/to/data` in the scripts are placeholders. Replace them with the local roots before running.
