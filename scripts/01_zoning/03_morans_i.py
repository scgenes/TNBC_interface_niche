# =============================================================================
# TNBC interface-niche code repository
# Original source: project_singlecell/TNBC_result/final_figures/V2/spatial/scripts/03_squidpy_svg.py
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

#!/usr/bin/env python3
"""
03_squidpy_svg.py — Squidpy spatial variable genes (Moran's I) for TNBC V2 spatial project.

Scientific question:
    Which genes show significant spatial autocorrelation (i.e. spatial variability)?
    Are the Fig2/Fig4 key genes (SPP1 / C1Q / MMP11 / YBX1) enriched at the
    tumor-immune boundary (Interface region)?

Input:
    - bridge/counts.mtx   : SCT log-normalized counts, genes x spots (23947 x 16607)
                            -> read then transposed to spots x genes
    - bridge/features.tsv : gene id / symbol
    - bridge/barcodes.tsv : spot barcodes
    - _archive_legacy/figure2/region_assignment_extended.csv
                            : barcode / sample / region / tumor_prop / x / y

Output (to V2/spatial/):
    - squidpy_svg.csv           : genome-wide Moran's I table (I, pval_norm, FDR), key-gene flag
    - squidpy_svg_keygenes.csv  : key-gene permutation p-values + region (Core/Interface/Stroma) means
    - fig_squidpy_svg.png/.pdf  : summary figure

Method:
    sq.gr.spatial_neighbors   : per-sample (12 samples, separate coordinate systems),
                                coord_type='generic', kNN k=7 minus self => 6 spatial neighbors
                                (block-diagonal connectivity, no cross-sample edges)
    sq.gr.spatial_autocorr    : mode='moran'; all genes with analytic p-value + FDR;
                                key genes additionally with n_perms=1000 permutation p-values

Environment:
    /path/to/data/project_singlecell/apps/miniconda3/envs/Squidpy/bin/python3
    (squidpy 1.6.1, scanpy 1.9.8, anndata 0.10.5, numpy 1.25.2 — Moran's I is self-contained,
     no esda/libpysal dependency)
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import io, stats
import scanpy as sc
import squidpy as sq
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------- paths
BRIDGE = "/path/to/project/project_singlecell/TNBC_result/final_figures/V2/spatial/bridge"
RA = "/path/to/project/project_singlecell/TNBC_result/_archive_legacy/figure2/region_assignment_extended.csv"
OUT = "/path/to/project/project_singlecell/TNBC_result/final_figures/V2/spatial"
os.makedirs(OUT, exist_ok=True)

KEY_GENES = ["SPP1", "C1QA", "C1QB", "C1QC", "MMP11", "YBX1"]
SEED = 42

# ----------------------------------------------------------------------------- load
print("[SVG] reading counts.mtx (genes x spots -> transpose)")
X = io.mmread(f"{BRIDGE}/counts.mtx").tocsr().T  # genes x spots -> spots x genes
var = pd.read_csv(f"{BRIDGE}/features.tsv", sep="\t", header=None)
var.columns = ["gene_ids", "gene_symbols"]
obs = pd.read_csv(f"{BRIDGE}/barcodes.tsv", sep="\t", header=None)
obs.columns = ["barcode"]

adata = sc.AnnData(X, obs=obs, var=var)
adata.obs_names = obs["barcode"].astype(str).values
adata.var_names = var["gene_symbols"].astype(str).values
print(f"[SVG] loaded {adata.shape[0]} spots x {adata.shape[1]} genes, nnz={adata.X.nnz}")

# metadata (sample / region / tumor_prop) + spatial coords (x/y)
meta = pd.read_csv(RA).set_index("barcode").loc[adata.obs_names]
adata.obs["sample"] = meta["sample"].astype(str).values
adata.obs["region"] = meta["region"].astype(str).values
adata.obs["tumor_prop"] = meta["tumor_prop"].astype(float).values
adata.obsm["spatial"] = meta[["x", "y"]].astype(float).values

print(f"[SVG] n_samples = {adata.obs['sample'].nunique()}")
print("[SVG] region counts:\n" + adata.obs["region"].value_counts().to_string())

# ----------------------------------------------------------------------------- per-sample spatial graph
# 12 samples have independent x/y coordinate systems -> build the graph PER SAMPLE
# (block-diagonal) so no cross-sample edges are created.
from scipy.sparse import lil_matrix

n = adata.n_obs
conn = lil_matrix((n, n), dtype=np.float32)
for s in sorted(adata.obs["sample"].unique()):
    idx = np.where(adata.obs["sample"].values == s)[0]
    sub = adata[idx].copy()
    sq.gr.spatial_neighbors(sub, coord_type="generic", n_neighs=7, key_added="spatial")
    c = sub.obsp["spatial_connectivities"].tolil()
    c.setdiag(0)  # remove self-loop -> 6 spatial neighbors (Visium hexagonal ring)
    c = c.tocsr()
    c.eliminate_zeros()
    conn[np.ix_(idx, idx)] = c
    print(f"[SVG] sample {s}: {len(idx)} spots, edges={c.nnz}")

adata.obsp["spatial_connectivities"] = conn.tocsr()
print(f"[SVG] total spatial edges (no self-loops): {adata.obsp['spatial_connectivities'].nnz}")

# ----------------------------------------------------------------------------- genome-wide Moran's I (analytic)
print("[SVG] Moran's I on all genes (analytic p-value + FDR), this may take a few minutes")
df = sq.gr.spatial_autocorr(
    adata,
    genes=None,          # all genes (no 'highly_variable' column present)
    mode="moran",
    n_perms=None,        # analytic p-value only (fast)
    corr_method="fdr_bh",
    n_jobs=1,
    seed=SEED,
    copy=True,
    show_progress_bar=True,
)
df = df.reset_index().rename(columns={"index": "gene"})
df["key_gene"] = df["gene"].isin(KEY_GENES)
df.to_csv(f"{OUT}/squidpy_svg.csv", index=False)
print(f"[SVG] genome-wide table: {df.shape[0]} genes -> squidpy_svg.csv")
print("[SVG] top 15 genes by Moran's I:")
print(df.head(15).to_string(index=False))

# ----------------------------------------------------------------------------- key genes: permutation p-values
print("[SVG] key-gene permutation test (n_perms=1000)")
kdf = sq.gr.spatial_autocorr(
    adata,
    genes=KEY_GENES,
    mode="moran",
    n_perms=1000,
    corr_method="fdr_bh",
    n_jobs=1,
    seed=SEED,
    copy=True,
    show_progress_bar=False,
)
kdf = kdf.reset_index().rename(columns={"index": "gene"})

# region means for key genes (boundary enrichment readout)
# expression matrix is SCT log-normalized; use sparse column means per region
Xk = adata[:, KEY_GENES].X
region_order = ["Tumor Core", "Interface", "Stroma"]
for r in region_order:
    m = (adata.obs["region"].values == r)
    kdf[f"mean_{r.replace(' ', '_')}"] = np.asarray(Xk[m].mean(axis=0)).ravel()
# log2FC Interface vs Tumor Core (boundary vs core)
kdf["log2FC_Interface_vs_Core"] = np.log2(
    (kdf["mean_Interface"] + 1e-6) / (kdf["mean_Tumor_Core"] + 1e-6)
)
# Mann-Whitney U, Interface vs Core (descriptive; Moran's I is the primary spatial test)
mw_p, mw_u = [], []
for g in KEY_GENES:
    gv = np.asarray(adata[:, g].X.toarray()).ravel()
    u, p = stats.mannwhitneyu(
        gv[adata.obs["region"].values == "Interface"],
        gv[adata.obs["region"].values == "Tumor Core"],
        alternative="two-sided",
    )
    mw_u.append(u)
    mw_p.append(p)
kdf["mwu_Interface_vs_Core_p"] = mw_p

kdf.to_csv(f"{OUT}/squidpy_svg_keygenes.csv", index=False)
print("[SVG] key-gene table:")
print(kdf.to_string(index=False))

# ----------------------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# Panel 1: volcano of Moran's I vs -log10 FDR p
ax = axes[0]
fdr = df["pval_norm_fdr_bh"].clip(lower=1e-300)
ax.scatter(df["I"], -np.log10(fdr), s=4, c="0.85", edgecolors="none", label="all genes")
kk = df[df["key_gene"]]
ax.scatter(kk["I"], -np.log10(kk["pval_norm_fdr_bh"].clip(lower=1e-300)),
           s=70, c="tab:red", edgecolors="black", linewidths=0.6, zorder=5, label="key genes")
for _, r in kk.iterrows():
    ax.annotate(r["gene"], (r["I"], -np.log10(max(r["pval_norm_fdr_bh"], 1e-300))),
                textcoords="offset points", xytext=(6, 4), fontsize=8)
ax.set_xlabel("Moran's I")
ax.set_ylabel("-log10(FDR p-value)")
ax.set_title("Spatial autocorrelation (Moran's I), all genes")
ax.legend(loc="upper left", fontsize=8)

# Panel 2: key gene Moran's I bar
ax = axes[1]
kdf_sorted = kdf.sort_values("I")
cols = ["tab:red" if v > 0 else "tab:blue" for v in kdf_sorted["I"]]
ax.barh(kdf_sorted["gene"], kdf_sorted["I"], color=cols)
ax.axvline(0, color="black", lw=0.8)
ax.set_xlabel("Moran's I")
ax.set_title("Key genes: spatial autocorrelation")
for i, (_, r) in enumerate(kdf_sorted.iterrows()):
    sig = "***" if r["pval_sim_fdr_bh"] < 1e-3 else ("**" if r["pval_sim_fdr_bh"] < 0.01 else ("*" if r["pval_sim_fdr_bh"] < 0.05 else "ns"))
    ax.text(r["I"] + 0.002, i, sig, va="center", fontsize=9)

# Panel 3: key gene mean expression by region (boundary enrichment)
ax = axes[2]
regions = ["Tumor_Core", "Interface", "Stroma"]
xpos = np.arange(len(KEY_GENES))
w = 0.26
for j, r in enumerate(regions):
    ax.bar(xpos + (j - 1) * w, kdf[f"mean_{r}"].values, w, label=r.replace("_", " "))
ax.set_xticks(xpos)
ax.set_xticklabels(KEY_GENES, rotation=45, ha="right")
ax.set_ylabel("mean SCT log-norm expression")
ax.set_title("Key genes: mean expression by region (Interface = tumor boundary)")
ax.legend(fontsize=8)

fig.suptitle("Squidpy spatial variable genes — TNBC V2 (tumor-immune boundary)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{OUT}/fig_squidpy_svg.png", dpi=150)
fig.savefig(f"{OUT}/fig_squidpy_svg.pdf")
print("[SVG] figure saved -> fig_squidpy_svg.png/.pdf")
print("[SVG] DONE")
