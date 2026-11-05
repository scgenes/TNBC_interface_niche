# 03_communication

Ligand–receptor communication analyses for the C1Q⁺ TAM ↔ MMP11⁺ CAF axis.

| Script | Purpose |
|---|---|
| `01_cellchat_spatial.R` | CellChat on the spatial/single-cell Seurat object, at two annotation levels (`Celltype_V2` major classes and `Celltype_sub` subtypes). Runs `createCellChat` → `computeCommunProb` (`triMean`, `population.size = FALSE`, `nboot = 100`) → `filterCommunication(min.cells = 10)` → `computeCommunProbPathway` → `aggregateNet`. Reads `sc_lean.rds`. |
| `02_nichenet_clone_LR.R` | Clone-level NicheNet implementation (converted from the original R notebook): receiver gene sets, candidate ligands and `nichenetr::predict_ligand_activities`. |
| `03_nichenet_CAF_C1Q.R` | **Main NicheNet analysis** for C1Q TAM ↔ MMP11 CAF in each zone. Direction A: C1Q TAM → MMP11 CAF; direction B: MMP11 CAF → C1Q TAM. Receiver gene set: genes correlated with sender abundance among receiver-positive spots (Pearson r ≥ 0.30, BH padj < 0.05, detection ≥ 10%, enrichment ≥ 0.25), top 200. Outputs ligand activity, receiver gene sets, predicted targets and `nichenet_results.rds`. |
| `04_cellchat_CAF_C1Q_by_region.R` | **Main CellChat analysis** for C1Q TAM ↔ MMP11 CAF, run per sample × zone with `distance.use = TRUE`, `interaction.range = 500`, `triMean`, `nboot = 100`, `min.cells = 10`. Spot assignment: C1Q TAM ≥ 2%, MMP11 CAF ≥ 5%; double-positive spots go to the higher proportion; others are `Other`. Outputs edge, LR and pathway tables by region. |
| `05_caf_c1q_panels.R` | Figure 5 / Figure S5 panels a–c: NicheNet activity and target heatmaps, and the CellChat spatial LR dot plot. Reads the result tables produced by `03_` and `04_`. |
| `06_nichenet_B_standalone.R` | NicheNet direction B standalone target heatmap, with collagen/stromal genes removed from the target panel (the relative-weight baseline is kept at the original 25-column scale). |

## Notes

- The **main CAF → C1Q TAM NicheNet** analysis is `03_nichenet_CAF_C1Q.R`; `02_nichenet_clone_LR.R` is the clone-level implementation.
- The single-cell CellChat panel (Figure S5d) is assembled from the precomputed `patient_ALL/cellchat_C1Q_MMP11` tables; the upstream CellChat computation is in `01_cellchat_spatial.R` (spatial) and the clone-level `../04_clone/04_clone_tam_cellchat.R`.
- Clone-level CellChat (Clone1/Clone2 → C1Q_TAM / SPP1_TAM) is in `../04_clone/04_clone_tam_cellchat.R`.

## Run order

1. `03_nichenet_CAF_C1Q.R`
2. `04_cellchat_CAF_C1Q_by_region.R`
3. `01_cellchat_spatial.R` (if the spatial/single-cell networks are needed)
4. `05_caf_c1q_panels.R` / `06_nichenet_B_standalone.R` (figure assembly)
5. `02_nichenet_clone_LR.R` (clone-level; requires clone inputs)

## Paths

`/path/to/project` and `/path/to/data` in the scripts are placeholders. Replace them with the local roots before running.
