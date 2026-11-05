# 02_niche

Spatial neighbourhood and CD8⁺ T-cell functional analyses.

| Script | Purpose |
|---|---|
| `01_spatial_lag_niche.R` | kNN (K = 6) spatial-lag framework: CAF-subtype specificity of the C1Q⁺ TAM niche, TAM → T-cell coupling, and spot-level export for representative-section maps. Uses within-sample permutation and Fisher-z cross-sample combination. Reads `figure20/subproportions.rds` and `figure20/spatial_meta.rds` via `lib/common.R`. |
| `02_spatial_colocalisation.py` | Continuous-proportion spatial cross-correlation on the Squidpy graph, with 1,000 permutations. Reads `bridge/subproportions.csv`, `bridge/barcodes.tsv` and `region_assignment_extended.csv`. |
| `03_neighbourhood_enrichment.py` | Independent re-implementation of the neighbourhood claim with `squidpy.gr.nhood_enrichment`, using two graph constructions (Visium hex grid and kNN = 6). Reads `panels/figure3/fig3_spots.csv`. |
| `04_cd8_functional_scores.R` | CD8⁺ T-cell functional scores by C1QC/MMP11 niche. CD8⁺ spots are the top 25% of RCTD CD8 proportion per slice; each gene is z-scored across CD8⁺ spots within slice, then averaged per gene set (cytotoxicity, exhaustion, TCR activation, Treg). Reads `Figure3/data/spot_table.qs` and `spatial_sct_matrix.rds`. Outputs per-spot and per-slice score tables plus paired Wilcoxon tests. |
| `05_cd8_activation_interface.R` | Interface-restricted CD8⁺ T-cell activation score using CD69 + GZMB + PRF1, with the same CD8⁺-spot definition. Outputs per-spot, per-slice, per-gene and test tables. |

## Run order

1. `01_spatial_lag_niche.R` (produces the spot-level table used below)
2. `02_spatial_colocalisation.py`
3. `03_neighbourhood_enrichment.py`
4. `04_cd8_functional_scores.R`
5. `05_cd8_activation_interface.R`

## Outputs

- `panels/figure3/*.csv` — spatial-lag summaries, per-sample results, and spot-level export.
- `squidpy_coloc.csv`, `squidpy_coloc_matrix.csv` — pairwise colocalisation.
- `nhood_z_*.csv`, `nhood_count_*.csv` — Squidpy neighbourhood enrichment.
- `fig3i_CD8_scores_per_spot.csv`, `fig3i_CD8_scores_per_slice.csv`, `fig3i_CD8_score_tests.csv` — CD8⁺ functional scores.
- `fig3i_Activation_Interface_*.csv` — interface activation scores and tests.

## Paths

`/path/to/project` and `/path/to/data` in the scripts are placeholders. Replace them with the local roots before running.
