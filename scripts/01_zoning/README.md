# 01_zoning

Zone definition, TAM/CAF composition and spatial-domain analyses.

| Script | Purpose |
|---|---|
| `01_zone_and_tam_composition.R` | Figure 2 analyses: spatial maps of RCTD tumour proportion, three-zone assignment summaries, TAM composition across zones, and the distance-to-core gradient. Reads `region_assignment.rds` and RCTD subtype proportions. |
| `02_banksy_domains.R` | BANKSY spatial-domain clustering and comparison with the RCTD three-zone labels (ARI), plus a non-spatial k-means control. Reads `integrated_seurat_12.rds` and `region_assignment.csv`. |
| `03_morans_i.py` | Squidpy Moran's I for genome-wide spatially variable genes and the six key genes; zone-wise mean expression. Reads `bridge/counts.mtx`, `bridge/features.tsv`, `bridge/barcodes.tsv` and `region_assignment_extended.csv`. |

## Run order

1. `01_zone_and_tam_composition.R`
2. `02_banksy_domains.R`
3. `03_morans_i.py`

## Inputs

The spot-level three-zone table (`region_assignment.rds` / `region_assignment_extended.csv`) and the RCTD proportion objects are produced upstream. See `../../data/README.md` and `../../docs/inputs_outputs.md`.

## Paths

`/path/to/project` and `/path/to/data` in the scripts are placeholders. Replace them with the local roots before running.
