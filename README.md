# TNBC interface-niche spatial analysis

Code accompanying the manuscript:

> **A spatially coupled C1Q⁺ TAM–MMP11⁺ CAF niche at the tumour–immune interface of triple-negative breast cancer**

This repository contains the analysis code used to define the tumour core / interface / stroma zones, quantify TAM and CAF spatial distributions, characterise the C1Q⁺ TAM–MMP11⁺ CAF neighbourhood, and relate tumour clones to TAM niches in triple-negative breast cancer (TNBC). It also contains the bulk survival validation shown in Figure S6, using GSE21653 (DFS, discovery) and GSE58812 (OS, validation); the survival script also retains TCGA-BRCA and METABRIC code paths from the exploratory analysis, which are not shown in the paper.

The code was assembled from three analysis projects:

- `project_singlecell/TNBC_result/final_figures` — zoning, TAM/CAF distribution, spatial niche and communication analyses.
- `project_spatial/TNBC/patient_ALL` — clone-level analyses (clone–TAM colocalisation, clone differential expression, clone-level CellChat/NicheNet, CNV confounding).
- `project_spatial/FigureAI` — CD8⁺ T-cell functional scoring and the C1Q TAM ↔ MMP11 CAF NicheNet/CellChat analyses (final figure-stage scripts).

Raw data and large intermediate objects are **not** included; see [`data/README.md`](data/README.md).

---

## Repository layout

```
.
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── config/
│   └── paths.example.yml
├── data/
│   └── README.md
├── docs/
│   ├── methods_to_code.md
│   ├── inputs_outputs.md
│   └── methods_CN.md
├── environment/
│   ├── environment.yml
│   ├── requirements.txt
│   └── sessionInfo.txt
└── scripts/
    ├── lib/
    │   ├── common.R
    │   └── shared_style.R
    ├── 01_zoning/
    │   ├── 01_zone_and_tam_composition.R
    │   ├── 02_banksy_domains.R
    │   ├── 03_morans_i.py
    │   └── README.md
    ├── 02_niche/
    │   ├── 01_spatial_lag_niche.R
    │   ├── 02_spatial_colocalisation.py
    │   ├── 03_neighbourhood_enrichment.py
    │   ├── 04_cd8_functional_scores.R
    │   ├── 05_cd8_activation_interface.R
    │   └── README.md
    ├── 03_communication/
    │   ├── 01_cellchat_spatial.R
    │   ├── 02_nichenet_clone_LR.R
    │   ├── 03_nichenet_CAF_C1Q.R
    │   ├── 04_cellchat_CAF_C1Q_by_region.R
    │   ├── 05_caf_c1q_panels.R
    │   ├── 06_nichenet_B_standalone.R
    │   └── README.md
    ├── 04_clone/
    │   ├── 00_spatial_infercnv_patient.R
    │   ├── 01_clone_inference_per_sample.R
    │   ├── 02_clone_tam_colocalisation.R
    │   ├── 03_clone_differential_expression.R
    │   ├── 04_clone_tam_cellchat.R
    │   ├── 05_cnv_confounding.R
    │   ├── 05b_cnv_adjusted_forest.R
    │   ├── 06_clone_tam_figure.py
    │   └── README.md
    └── 05_survival/
        ├── 01_build_gene_lists.py
        ├── 02_survival_analysis.R
        ├── 03_combined_figure.R
        └── README.md
```

All scripts carry a header with their original source file. Absolute paths in the original analysis environment have been replaced by `/path/to/project` and `/path/to/data`; update these to your local roots before running.

---

## Data

| Data | Source | Description |
|---|---|---|
| Single-cell RNA-seq | GEO **GSE176078** | 12 TNBC samples, 38,947 cells; Wu et al., *Nat Genet* 2021 |
| Spatial transcriptomics | GEO **GSE210616** | 12 Visium sections from 6 TNBC patients, 16,607 spots; Bassiouni et al., *Cancer Res* 2023 |

Expected intermediate objects are listed in [`data/README.md`](data/README.md) and [`docs/inputs_outputs.md`](docs/inputs_outputs.md). They include the integrated Seurat object, RCTD proportion matrices, the spot-level zone table and per-patient clone membership tables. These objects are available from the authors on reasonable request.

---

## Environment

The analyses used R 4.3.3 with Seurat 4.4.0, spacexr 2.2.1, inferCNV 1.22.0, BANKSY 0.99.13, CellChat 2.1.2 and nichenetr 2.2.1, and Python 3.10 with Squidpy 1.2.2 (1.6.1 in some individual scripts; see script headers). Exact versions are listed in [`environment/sessionInfo.txt`](environment/sessionInfo.txt) and [`environment/environment.yml`](environment/environment.yml).

```bash
# Create the conda environment (R + Python)
conda env create -f environment/environment.yml
conda activate tnbc-interface-niche

# Or, for a pip-only Python environment:
pip install -r environment/requirements.txt
```

Squidpy analyses were run in a Python environment with `squidpy==1.2.2`, `scanpy`, `anndata`, `pandas`, `numpy`, `matplotlib` and `seaborn`.

---

## Quick start

1. Clone the repository.
2. Set the path placeholders in the scripts (or copy `config/paths.example.yml` and use it as a checklist):
   - `/path/to/project` → directory containing `project_singlecell/TNBC_result` and `project_spatial/TNBC`.
   - `/path/to/data` → directory containing the Visium/CNV inputs.
3. Place the intermediate objects listed in [`data/README.md`](data/README.md).
4. Run the modules in order **from the repository root** (scripts source helpers from `scripts/lib/`):

```text
01_zoning        →  02_niche        →  03_communication
                                      ↘
                                       04_clone        →  05_survival
```

Each module has its own `README.md` with inputs, outputs and run order.

---

## Methods-to-code map

| Methods section | Analysis | Main script(s) |
|---|---|---|
| §3.2 | Spatial inferCNV and clone assignment | `04_clone/00_spatial_infercnv_patient.R`, `04_clone/01_clone_inference_per_sample.R` |
| §3.3 | Three-zone definition (RCTD `tumor_prop`) | `01_zoning/01_zone_and_tam_composition.R` (uses the upstream `region_assignment` table) |
| §4 | RCTD deconvolution | upstream; spot-level proportion matrices are inputs |
| §5 | BANKSY domain concordance | `01_zoning/02_banksy_domains.R` |
| §6 | Moran's I / spatially variable genes | `01_zoning/03_morans_i.py` |
| §7 | TAM composition and distance-to-core gradient | `01_zoning/01_zone_and_tam_composition.R` |
| §8 | Neighbourhood enrichment (Ro/e) | `02_niche/01_spatial_lag_niche.R`, `02_niche/03_neighbourhood_enrichment.py` |
| §9 | Region-stratified spatial colocalisation | `02_niche/02_spatial_colocalisation.py` |
| §10 | CD8⁺ T-cell functional scoring | `02_niche/04_cd8_functional_scores.R`, `02_niche/05_cd8_activation_interface.R` |
| §11.1 | CellChat (spatial) | `03_communication/01_cellchat_spatial.R` |
| §11.1 | CellChat (C1Q TAM ↔ MMP11 CAF, by region) | `03_communication/04_cellchat_CAF_C1Q_by_region.R` |
| §11.1 | CellChat (clone level) | `04_clone/04_clone_tam_cellchat.R` |
| §11.2 | NicheNet (C1Q TAM ↔ MMP11 CAF, main) | `03_communication/03_nichenet_CAF_C1Q.R` |
| §11.2 | NicheNet (clone level) | `03_communication/02_nichenet_clone_LR.R` |
| §12.1 | Clone–TAM colocalisation | `04_clone/02_clone_tam_colocalisation.R` |
| §12.2 | Clone differential expression | `04_clone/03_clone_differential_expression.R` |
| §12 | CNV confounding | `04_clone/05_cnv_confounding.R`, `04_clone/05b_cnv_adjusted_forest.R` |
| Figure 5 | CAF–C1Q TAM communication panels | `03_communication/05_caf_c1q_panels.R`, `03_communication/06_nichenet_B_standalone.R` |
| Figure assembly | Clone–TAM figures | `04_clone/06_clone_tam_figure.py` |
| §2.11, §3.8 | Bulk survival analysis of the three signatures (Figure S6) | `05_survival/02_survival_analysis.R`, `05_survival/03_combined_figure.R` |
| §2.11 | Signature gene lists from the single-cell DE tables | `05_survival/01_build_gene_lists.py` |

A more detailed table, including inputs and outputs, is in [`docs/methods_to_code.md`](docs/methods_to_code.md).

---

## Notes on completeness

- **Converted notebooks.** Several scripts were converted from the original R notebooks with `jupyter nbconvert --to script` (the `04_clone` scripts converted from notebooks, `03_communication/02_nichenet_clone_LR.R`, and `05_survival/02_survival_analysis.R`). They retain the original in-notebook QC/audit cells from the analysis environment. These cells are not required for the scientific workflow and can be skipped or deleted. The CD8⁺ scoring and CAF–C1Q communication scripts from `FigureAI` are original `.R` scripts.
- **Upstream objects.** Single-cell preprocessing (Seurat/Harmony), RCTD deconvolution and the final three-zone table were produced upstream of `final_figures`; the corresponding objects are inputs here. The script that plots and summarises the zones (`01_zone_and_tam_composition.R`) is included.
- **Per-patient clone inference.** The included clone script is a per-sample implementation. The per-patient tree-node implementation referenced in the Methods (`run_clone_patient.R`) was not present in the archived project; the per-patient clone outputs (`clone_membership.csv`, `clone_k_selection.csv`, etc.) are inputs to the downstream scripts.

---

## Citation

If you use this code, please cite the manuscript above. A `CITATION.cff` template is provided; update author and DOI fields when the paper is published.

## License

MIT License. See [`LICENSE`](LICENSE). Change the license if your institution requires a different one.
