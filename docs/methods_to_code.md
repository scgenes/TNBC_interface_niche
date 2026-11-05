# Methods-to-code map

This table maps each analysis in the manuscript to the corresponding script in this repository. Paths refer to `scripts/`.

| Methods | Analysis | Script | Main inputs | Main outputs |
|---|---|---|---|---|
| §3.2 | Spatial inferCNV, patient-level CNV matrix, HMM states | `04_clone/00_spatial_infercnv_patient.R` | `run.final.infercnv_obj`, HMM gene-level predictions, gene order, spatial coordinates | CNV matrix, HMM state table, zone-by-CNV summaries |
| §3.2 | Clone inference (k selection, tree, arm consensus) | `04_clone/01_clone_inference_per_sample.R` | `run.final.infercnv_obj`, HMM predictions, gene order | `clone_membership.csv`, `clone_k_selection.csv`, `clone_arm_consensus.csv`, Newick trees |
| §3.3, §7 | Three-zone summaries, TAM composition, distance-to-core gradient | `01_zoning/01_zone_and_tam_composition.R` | `region_assignment.rds`, RCTD subtype proportions | Zone maps, composition tables, gradient summaries, figure panels |
| §4 | RCTD deconvolution | upstream (objects are inputs) | scRNA reference, Visium counts | `subproportions.rds`, `spatial_meta.rds`, `spatial_tme_groups.qs` |
| §5 | BANKSY spatial domains and ARI vs zones | `01_zoning/02_banksy_domains.R` | `integrated_seurat_12.rds`, `region_assignment.csv` | `banksy_ari.csv`, domain tables, ARI/kmeans figures |
| §6 | Moran's I / spatially variable genes | `01_zoning/03_morans_i.py` | `bridge/counts.mtx`, `features.tsv`, `barcodes.tsv`, `region_assignment_extended.csv` | `squidpy_svg.csv`, `squidpy_svg_keygenes.csv`, figure |
| §8 | Spatial lag, CAF/TAM coupling, per-sample permutation | `02_niche/01_spatial_lag_niche.R` | `figure20/subproportions.rds`, `figure20/spatial_meta.rds` | Fig. 3 data tables in `panels/figure3/` |
| §8 | Independent Squidpy neighbourhood enrichment check | `02_niche/03_neighbourhood_enrichment.py` | `panels/figure3/fig3_spots.csv` | `nhood_z_*.csv`, `nhood_count_*.csv`, label counts |
| §9 | Region-stratified spatial colocalisation | `02_niche/02_spatial_colocalisation.py` | `bridge/subproportions.csv`, `barcodes.tsv`, `region_assignment_extended.csv` | `squidpy_coloc.csv`, `squidpy_coloc_matrix.csv`, figures |
| §10 | CD8⁺ T-cell cytotoxicity / activation scores | `02_niche/04_cd8_functional_scores.R`, `02_niche/05_cd8_activation_interface.R` | `Figure3/data/spot_table.qs`, `spatial_sct_matrix.rds` | `fig3i_CD8_scores_per_spot.csv`, `fig3i_CD8_scores_per_slice.csv`, `fig3i_CD8_score_tests.csv`, `fig3i_Activation_Interface_*` |
| §11.1 | CellChat, spatial and single-cell | `03_communication/01_cellchat_spatial.R` | `sc_lean.rds` (Seurat object with `Celltype_V2`, `Celltype_sub`) | `cellchat_*.rds`, communication tables, logs |
| §11.1 | CellChat, C1Q TAM ↔ MMP11 CAF by region | `03_communication/04_cellchat_CAF_C1Q_by_region.R` | `integrated_seurat_12.rds`, `subcluster_prop.csv.gz`, `region_assignment.csv` | `cellchat_by_region/tables/*`, figures |
| §11.1 | CellChat, clone level | `04_clone/04_clone_tam_cellchat.R` | spatial SCT matrix, TME proportions, `clone_membership.csv` | CellChat objects, LR/pathway tables, figures |
| §11.2 | NicheNet, C1Q TAM ↔ MMP11 CAF (main) | `03_communication/03_nichenet_CAF_C1Q.R` | spatial SCT matrix, `subcluster_prop.csv.gz`, `region_assignment_extended.csv`, NicheNet priors | `nichenet_C1Q_MMP11CAF/tables/*`, `nichenet_results.rds` |
| §11.2 | NicheNet, clone-level implementation | `03_communication/02_nichenet_clone_LR.R` | clone/spatial expression, receiver gene sets | Ligand activity tables, target weight tables |
| §12.1 | Clone–TAM colocalisation (radius + kNN Ro/e) | `04_clone/02_clone_tam_colocalisation.R` | `clone_membership.csv`, `region_assignment_extended.csv`, `spatial_tme_groups.qs` | `clone_tam_radius_perm_by_slice.csv`, `clone_tam_radius_perm_region_stratified.csv`, `clone_tam_coloc_samespot.csv`, `clone_tam_coloc_neighbor.csv` |
| §12.2 | Clone1 vs Clone2 differential expression | `04_clone/03_clone_differential_expression.R` | `clone_seurat.rds` (P13, P15) | DE tables, volcano plots |
| §12 | CNV confounding / adjusted models | `04_clone/05_cnv_confounding.R`, `04_clone/05b_cnv_adjusted_forest.R` | Clone membership, CNV scores, TAM niche scores | Model tables, forest plot |
| Figure 5 | CAF–C1Q TAM communication panels | `03_communication/05_caf_c1q_panels.R`, `06_nichenet_B_standalone.R` | NicheNet/CellChat result tables | Figure 5 / Figure S5 panels |
| Figure | Clone–TAM figure assembly | `04_clone/06_clone_tam_figure.py` | Tables in `TAM_clone/tables/` | Main and supplementary figure panels |
| §2.11 | Signature gene lists from the Mac_sub / Fib_sub DE tables | `05_survival/01_build_gene_lists.py` | `Mac_sub_差异列表.tsv`, `Fib_sub_差异列表.tsv` | `gene_list_top20_long.tsv`, `gene_list_top20_wide.csv`, `gene_list_unique.txt` |
| §2.11, §3.8 | Bulk survival analysis in TCGA/METABRIC/GSE58812/GSE21653 (Figure S6) | `05_survival/02_survival_analysis.R`, `05_survival/03_combined_figure.R` | cBioPortal clinical/expression files, GEO matrices, `data/signatures/gene_list_top20_wide.csv` | `survival_results_primary.csv`, `stratification_cutpoints.csv`, Figure S6 PDF/TIFF |

## Notes

- **Converted notebooks.** `00_spatial_infercnv_patient.R`, `02_clone_tam_colocalisation.R`, `03_clone_differential_expression.R`, `04_clone_tam_cellchat.R` and `02_nichenet_clone_LR.R` were converted from R notebooks with `jupyter nbconvert --to script`. They retain the in-notebook QC/audit cells from the original analysis environment; those cells are not part of the scientific pipeline. The CD8⁺ scoring and CAF–C1Q communication scripts from FigureAI are original `.R` scripts.
- **RCTD upstream.** The RCTD deconvolution and the final spot-level zone table were produced before the `final_figures` stage. The zone table (`region_assignment.rds` / `region_assignment_extended.csv`) and the RCTD proportion objects are inputs to this repository.
- **Missing script.** The per-patient tree-node clone script (`run_clone_patient.R`) referenced in the Methods was not present in the archived project; per-patient clone outputs are inputs here.
