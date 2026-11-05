# lib

Shared helper functions used by the analysis scripts.

| File | Purpose |
|---|---|
| `common.R` | `load_data()` (RCTD subtype proportions + spatial metadata), `build_W()` (per-sample kNN spatial graph), `dist_to_core_um()`, spatial-lag helpers and colour palettes. Sourced by `02_niche/01_spatial_lag_niche.R`. |
| `shared_style.R` | Shared colours and ggplot theme for the figure panels. Sourced by `01_zoning/01_zone_and_tam_composition.R`. |

The `DATA_DIR` in `common.R` is a placeholder (`/path/to/project/project_singlecell/TNBC_result`); update it to the local project root.
