# Inputs and outputs

The table lists the intermediate objects expected by the scripts and the main outputs produced. File names follow the original analysis projects.

| File | Produced by | Used by | Content / key columns |
|---|---|---|---|
| `integrated_seurat_12.rds` | upstream single-cell/spatial integration | `01_zoning/02_banksy_domains.R` | Seurat object with SCT assay and spot metadata |
| `region_assignment.rds` | upstream RCTD + zone definition | `01_zoning/01_zone_and_tam_composition.R` | spot × `sample`, `x`, `y`, `tumor_prop`, `region` (`Tumor Core` / `Interface` / `Stroma`) |
| `region_assignment_extended.csv` | upstream RCTD + zone definition | `01_zoning/03_morans_i.py`, `02_niche/02_spatial_colocalisation.py`, `04_clone/02_clone_tam_colocalisation.R` | barcode, sample, region, tumor_prop, x, y |
| `subproportions.rds` | upstream RCTD (35 subtypes) | `02_niche/01_spatial_lag_niche.R`, `01_zoning/01_zone_and_tam_composition.R` | spots × 35 subtype proportions |
| `spatial_meta.rds` | upstream spatial metadata | `02_niche/01_spatial_lag_niche.R` | barcode, sample, x, y, region |
| `spatial_tme_groups.qs` | upstream RCTD (9 major classes + subtypes) | `04_clone/02_clone_tam_colocalisation.R`, `04_clone/04_clone_tam_cellchat.R` | spots × cell-type proportions, including `C1QC_TAM`, `SPP1_TAM`, `MMP11_CAF` |
| `spatial_sct_matrix.rds` | upstream | `02_niche/04_cd8_functional_scores.R`, `02_niche/05_cd8_activation_interface.R`, `03_communication/03_nichenet_CAF_C1Q.R`, `04_clone/04_clone_tam_cellchat.R` | SCT log-normalized matrix, genes × spots |
| `sc_lean.rds` | upstream | `03_communication/01_cellchat_spatial.R` | Seurat object with `Celltype_V2` and `Celltype_sub` |
| `spot_table.qs` | upstream | `02_niche/04_cd8_functional_scores.R`, `05_cd8_activation_interface.R` | barcode, slice, region, `CD8_sum`, `CD4_sum`, `C1QC_pos`, `MMP11_pos` |
| `subcluster_prop.csv.gz` | upstream RCTD (35 subclusters) | `03_communication/03_nichenet_CAF_C1Q.R`, `04_cellchat_CAF_C1Q_by_region.R` | spots × 35 subcluster proportions |
| `region_assignment.csv` | upstream | `03_communication/04_cellchat_CAF_C1Q_by_region.R` | barcode, region, x, y |
| `ligand_target_matrix.rds`, `lr_network.rds` | NicheNet v2 priors | `03_communication/03_nichenet_CAF_C1Q.R` | ligand–target matrix and ligand–receptor network |
| `clone_membership.csv` | `04_clone/01_clone_inference_per_sample.R` / upstream per-patient clone pipeline | `04_clone/02_clone_tam_colocalisation.R`, `04_clone/03_clone_differential_expression.R`, `04_clone/04_clone_tam_cellchat.R` | `cell` = `<slice>_<barcode>`, clone label |
| `clone_seurat.rds` | upstream clone pipeline | `04_clone/03_clone_differential_expression.R` | Seurat object with clone labels (P13, P15) |
| `figure20/subproportions.rds` | upstream RCTD | `02_niche/01_spatial_lag_niche.R` | spots × subtype proportions (input to `load_data()` in `lib/common.R`) |
| `figure20/spatial_meta.rds` | upstream | `02_niche/01_spatial_lag_niche.R` | barcode, sample, x, y, region |

## Main outputs

| Output | Script | Content |
|---|---|---|
| `squidpy_svg.csv` | `01_zoning/03_morans_i.py` | genome-wide Moran's I, p-values, FDR |
| `squidpy_svg_keygenes.csv` | `01_zoning/03_morans_i.py` | key-gene permutation p-values and zone means |
| `squidpy_coloc.csv` / `_matrix.csv` | `02_niche/02_spatial_colocalisation.py` | pairwise spatial colocalisation correlations and permutation p-values |
| `nhood_z_*.csv` / `nhood_count_*.csv` | `02_niche/03_neighbourhood_enrichment.py` | Squidpy neighbourhood enrichment z-scores and counts |
| `fig3i_CD8_scores_per_spot.csv` / `_per_slice.csv` / `_tests.csv` | `02_niche/04_cd8_functional_scores.R` | CD8⁺ T-cell functional scores and paired tests |
| `fig3i_Activation_Interface_*.csv` | `02_niche/05_cd8_activation_interface.R` | Interface activation scores, per-gene values and tests |
| `nichenet_C1Q_MMP11CAF/tables/*` | `03_communication/03_nichenet_CAF_C1Q.R` | ligand activity, receiver gene sets and predicted targets |
| `cellchat_by_region/tables/*` | `03_communication/04_cellchat_CAF_C1Q_by_region.R` | CellChat edges, LR and pathway summaries by region |
| `panels/figure3/*.csv` | `02_niche/01_spatial_lag_niche.R` | spatial-lag summaries, per-sample results, spot-level export |
| `cellchat_*.rds` | `03_communication/01_cellchat_spatial.R` | CellChat objects |
| clone–TAM tables | `04_clone/02_clone_tam_colocalisation.R` | radius and kNN Ro/e, per-slice and region-stratified results |
| DE tables | `04_clone/03_clone_differential_expression.R` | Clone1 vs Clone2 Wilcoxon results |
| model tables / forest plot | `04_clone/05_cnv_confounding.R`, `05b_cnv_adjusted_forest.R` | CNV-adjusted linear probability models |
| figure panels | `04_clone/06_clone_tam_figure.py` | main and supplementary clone–TAM figures |
| `survival_results_primary.csv` | `05_survival/02_survival_analysis.R` | Cox/KM results for 4 cohorts × 2 endpoints × 2 scoring methods × 2 cutoff methods × 4 analysis groups |
| `stratification_cutpoints.csv` | `05_survival/02_survival_analysis.R` | median/maxstat cutpoints, selection p, group sizes |
| `combined_panel_source_2cohorts.csv` | `05_survival/03_combined_figure.R` | 8 rows used in Figure S6 (GSE21653 DFS discovery, GSE58812 OS validation) |
| `FigureS6_survival_180mm.{pdf,png,tif}` | `05_survival/03_combined_figure.R` | submission figure (180 mm, 600 dpi LZW TIFF + PDF) |
