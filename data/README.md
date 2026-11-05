# Data

Raw and processed data are **not** included in this repository. This file describes what to obtain and where to place it.

## Raw data

| Dataset | Accession | Content |
|---|---|---|
| Single-cell RNA-seq | [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) | 12 TNBC samples, 38,947 cells (Wu et al., *Nat Genet* 2021) |
| Spatial transcriptomics | [GSE210616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616) | 12 10x Visium sections from 6 TNBC patients (P9, P11–P15), 16,607 spots (Bassiouni et al., *Cancer Res* 2023) |

The 12 Visium sections are 092A, 092B, 093A–093D, 094A–094D, 095A and 095B.

Download the processed matrices from the GEO supplementary files (or from the authors' processed objects). Raw FASTQ files are not required for the analyses in this repository.

## Bulk survival cohorts (Figure S6)

| Dataset | Source | Content |
|---|---|---|
| TCGA-BRCA | cBioPortal `brca_tcga` | clinical patient/sample tables, `data_mrna_seq_v2_rsem.txt` |
| METABRIC | cBioPortal `brca_metabric` | clinical patient/sample tables, `data_mrna_illumina_microarray.txt` |
| GSE58812 | GEO | expression matrix, clinical table (OS) |
| GSE21653 | GEO | expression matrix, clinical table (DFS) |

The 20-gene signature lists used in the survival analysis are included in
[`signatures/`](signatures/): `gene_list_top20_wide.csv`,
`gene_list_top20_long.tsv`, `gene_list_unique.txt`, and `gene_list_stats.csv`.
They were derived from the Mac_sub / Fib_sub differential-expression tables with
`scripts/05_survival/01_build_gene_lists.py`.

## Expected directory layout

```
data/
├── scRNA/                         # GSE176078 processed matrices / Seurat object
├── visium/                        # GSE210616 per-section filtered_feature_bc_matrix + spatial/
├── signatures/                    # 20-gene bulk signatures (included)
└── intermediate/
    ├── integrated_seurat_12.rds
    ├── region_assignment.rds
    ├── region_assignment_extended.csv
    ├── subproportions.rds
    ├── spatial_meta.rds
    ├── spatial_tme_groups.qs
    ├── spatial_sct_matrix.rds
    ├── sc_lean.rds
    ├── spot_table.qs
    ├── subcluster_prop.csv.gz
    ├── ligand_target_matrix.rds
    ├── lr_network.rds
    └── clones/
        ├── P9/clone_membership.csv
        ├── P12/clone_membership.csv
        ├── P13/clone_membership.csv
        ├── P14/clone_membership.csv
        └── P15/clone_membership.csv
```

See [`../docs/inputs_outputs.md`](../docs/inputs_outputs.md) for the contents and key columns of each intermediate file.

## Availability

The intermediate objects (RCTD proportions, spot-level zone table, clone membership tables) are available from the authors on reasonable request. They are omitted here because of size and because they can be regenerated from the public raw data with the upstream pipeline.
