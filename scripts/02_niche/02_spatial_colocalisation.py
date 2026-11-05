# =============================================================================
# TNBC interface-niche code repository
# Original source: project_singlecell/TNBC_result/final_figures/V2/spatial/scripts/04_squidpy_coloc.py
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

#!/usr/bin/env python3
"""
04_squidpy_coloc.py — Squidpy-graph spatial co-localization of deconvolved cell types,
TNBC V2 spatial project.

Scientific question:
    From 35 subpopulation deconvolution proportions, which cell-type pairs co-localize
    (spatially co-vary) or mutually exclude? Focus: do T-cell subtypes spatially exclude
    SPP1 TAMs, and co-localize with C1Q TAMs (Fig3 immune-exclusion at the neighborhood level)?

Method:
    Continuous-proportion spatial cross-correlation on the squidpy spatial graph.

    WHY NOT argmax + sq.gr.nhood_enrichment:
        The dominant cell type per spot (argmax over 35 proportions) is degenerate here —
        "Tumor Epithelial cells" wins 15428/16607 spots (92.9%), so immune/stromal types are
        never dominant and nhood_enrichment collapses to ~3 tumor-dominated categories.
        Co-localization of minority cell types MUST be measured on the continuous proportions.

    Statistic (basic statistic — correlation on the spatial-lag):
        P  : (n_spots x 35) deconvolution proportions
        W  : row-normalized per-sample spatial adjacency (squidpy sq.gr.spatial_neighbors,
             kNN k=7 minus self -> 6 neighbors, block-diagonal across 12 samples)
        L  = W @ P                    : neighborhood-mean proportion per cell type
        C[A,B] = Pearson(P[:,A], L[:,B])   : A's own abundance vs B's neighborhood abundance
        C_sym = (C + C.T) / 2         : symmetrized pairwise co-localization matrix
        C_sym > 0 -> co-localize ; C_sym < 0 -> mutual exclusion ; diagonal = self-aggregation
    Significance: permutation test (n_perms=1000), shuffling spot-to-proportion assignment
        (rows of P) to break spatial structure; two-tailed p = (#(|null|>=|obs|)+1)/(n_perms+1).

Input:
    - bridge/subproportions.csv : barcode + 35 subpopulation proportions
    - bridge/barcodes.tsv       : spot barcodes
    - _archive_legacy/figure2/region_assignment_extended.csv : sample / x / y

Output (to V2/spatial/):
    - squidpy_coloc.csv        : long-format pairwise table (A, B, corr, p_perm, significant)
    - squidpy_coloc_matrix.csv : full symmetrized cross-correlation matrix (35 x 35)
    - fig_squidpy_coloc.png/.pdf : heatmaps (full + T-cell x TAM focused)

Environment:
    /path/to/data/project_singlecell/apps/miniconda3/envs/Squidpy/bin/python3
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse import lil_matrix
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

SEED = 42
N_PERMS = 1000

# ----------------------------------------------------------------------------- load
print("[COLOC] loading subproportions + barcodes + coords")
subp = pd.read_csv(f"{BRIDGE}/subproportions.csv")
obs = pd.read_csv(f"{BRIDGE}/barcodes.tsv", sep="\t", header=None)
obs.columns = ["barcode"]

subp = subp.set_index("barcode").loc[obs["barcode"].astype(str).values]
sub_cols = [c for c in subp.columns if c != "barcode"]
P = subp[sub_cols].values.astype(np.float64)  # (n_spots, 35)
print(f"[COLOC] {P.shape[0]} spots x {P.shape[1]} subpopulations")

meta = pd.read_csv(RA).set_index("barcode").loc[obs["barcode"].astype(str).values]
sample = meta["sample"].astype(str).values
coords = meta[["x", "y"]].astype(float).values

# ----------------------------------------------------------------------------- per-sample spatial graph
# build block-diagonal row-normalized adjacency W (n x n); X content is unused by
# spatial_neighbors (it only reads obsm['spatial']), so a minimal AnnData suffices.
n = P.shape[0]
conn = lil_matrix((n, n), dtype=np.float32)
for s in sorted(np.unique(sample)):
    idx = np.where(sample == s)[0]
    sub = sc.AnnData(X=np.zeros((len(idx), 1)), obsm={"spatial": coords[idx]})
    sq.gr.spatial_neighbors(sub, coord_type="generic", n_neighs=7, key_added="spatial")
    c = sub.obsp["spatial_connectivities"].tolil()
    c.setdiag(0)
    c = c.tocsr()
    c.eliminate_zeros()
    conn[np.ix_(idx, idx)] = c
W = conn.tocsr().astype(np.float64)
# row-normalize (each row sums to 1 -> L = W @ P is the neighborhood-mean proportion)
rs = np.asarray(W.sum(axis=1)).ravel()
rs[rs == 0] = 1.0
W = sp.diags(1.0 / rs) @ W
print(f"[COLOC] spatial edges (no self-loops): {W.nnz}")

# ----------------------------------------------------------------------------- statistic
def cross_corr(Pmat, Lmat):
    """(n,k) x (n,k) -> (k,k) Pearson correlation matrix (columns standardized)."""
    Pc = Pmat - Pmat.mean(axis=0, keepdims=True)
    Lc = Lmat - Lmat.mean(axis=0, keepdims=True)
    Pn = Pc / np.sqrt((Pc ** 2).sum(axis=0, keepdims=True) + 1e-300)
    Ln = Lc / np.sqrt((Lc ** 2).sum(axis=0, keepdims=True) + 1e-300)
    return Pn.T @ Ln

k = P.shape[1]
L = W @ P
C_obs = cross_corr(P, L)                      # asymmetric (A own vs B neighborhood)
C_sym_obs = (C_obs + C_obs.T) / 2.0           # symmetrized pairwise matrix

print(f"[COLOC] permutation test n_perms={N_PERMS}")
rng = np.random.default_rng(SEED)
C_sym_null = np.empty((N_PERMS, k, k))
for i in range(N_PERMS):
    perm = rng.permutation(n)
    P_perm = P[perm]
    L_perm = W @ P_perm
    C = cross_corr(P_perm, L_perm)
    C_sym_null[i] = (C + C.T) / 2.0
    if (i + 1) % 250 == 0:
        print(f"  ... {i+1}/{N_PERMS} perms done")

# two-tailed permutation p-value (with +1 adjustment)
p_perm = (np.abs(C_sym_null) >= np.abs(C_sym_obs)).sum(axis=0) + 1
p_perm = p_perm / (N_PERMS + 1)

C_sym = pd.DataFrame(C_sym_obs, index=sub_cols, columns=sub_cols)
Pval = pd.DataFrame(p_perm, index=sub_cols, columns=sub_cols)

# ----------------------------------------------------------------------------- save
rows = []
for a in sub_cols:
    for b in sub_cols:
        rows.append({
            "cell_type_A": a, "cell_type_B": b,
            "spatial_corr": C_sym.loc[a, b],
            "p_perm": Pval.loc[a, b],
            "significant": Pval.loc[a, b] < 0.05,
        })
pair_df = pd.DataFrame(rows)
pair_df.to_csv(f"{OUT}/squidpy_coloc.csv", index=False)
C_sym.to_csv(f"{OUT}/squidpy_coloc_matrix.csv")
print(f"[COLOC] pairwise table -> squidpy_coloc.csv ({len(pair_df)} pairs)")

# ----------------------------------------------------------------------------- readout: T cells vs SPP1/C1Q TAMs
tcell_cols = [c for c in sub_cols if "T cells" in c]
tam_cols = [c for c in sub_cols if "TAMs" in c]
cd8_cols = [c for c in tcell_cols if c.startswith("CD8")]
cd4_cols = [c for c in tcell_cols if c.startswith("CD4")]
print(f"[COLOC] T-cell subtypes: {len(tcell_cols)} (CD4={len(cd4_cols)}, CD8={len(cd8_cols)}, other={len(tcell_cols)-len(cd4_cols)-len(cd8_cols)})")
print(f"[COLOC] TAM subtypes: {len(tam_cols)}")

print("\n[COLOC] T-cell x SPP1/C1Q TAM readout (spatial_corr; + co-localize / - exclude):")
for tam in ["SPP1 TAMs", "C1Q TAMs"]:
    print(f"  --- {tam} ---")
    for t in sorted(tcell_cols):
        v = C_sym.loc[t, tam]
        p = Pval.loc[t, tam]
        sig = "***" if p < 1e-3 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
        print(f"    {t:24s} corr={v:+.3f}  p_perm={p:.3f} ({sig})")

# aggregate readouts
def agg(types, tam):
    vals = [C_sym.loc[t, tam] for t in types]
    ps = [Pval.loc[t, tam] for t in types]
    n_sig_pos = sum(1 for v, p in zip(vals, ps) if v > 0 and p < 0.05)
    n_sig_neg = sum(1 for v, p in zip(vals, ps) if v < 0 and p < 0.05)
    return np.mean(vals), n_sig_pos, n_sig_neg

print("\n[COLOC] summary (mean corr, #sig-pos, #sig-neg):")
for grp, cols in [("CD8", cd8_cols), ("CD4", cd4_cols), ("all T", tcell_cols)]:
    for tam in ["SPP1 TAMs", "C1Q TAMs"]:
        m, npos, nneg = agg(cols, tam)
        print(f"    {grp:7s} x {tam:12s}: mean_corr={m:+.3f}, sig_pos={npos}, sig_neg={nneg}")

# ----------------------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 2, figsize=(20, 9), gridspec_kw={"width_ratios": [1.3, 1]})

# Panel 1: full cross-correlation heatmap
ax = axes[0]
vmax = np.nanmax(np.abs(C_sym.values))
im = ax.imshow(C_sym.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(k)); ax.set_yticks(range(k))
ax.set_xticklabels(sub_cols, rotation=90, fontsize=6)
ax.set_yticklabels(sub_cols, fontsize=6)
ax.set_title("Spatial cross-correlation (all 35 subpopulations)", fontsize=11)
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="spatial cross-correlation (co-loc + / excl -)")

# Panel 2: focused T-cell x TAM heatmap
ax = axes[1]
submat = C_sym.loc[sorted(tcell_cols), tam_cols]
subp_ = Pval.loc[sorted(tcell_cols), tam_cols]
vmax2 = np.nanmax(np.abs(submat.values)) if submat.size else 1.0
im2 = ax.imshow(submat.values, cmap="RdBu_r", vmin=-vmax2, vmax=vmax2, aspect="auto")
ax.set_xticks(range(len(tam_cols))); ax.set_yticks(range(len(submat.index)))
ax.set_xticklabels(tam_cols, rotation=45, ha="right", fontsize=7)
ax.set_yticklabels(submat.index, fontsize=7)
for i in range(len(submat.index)):
    for j in range(len(tam_cols)):
        v = submat.values[i, j]
        p = subp_.values[i, j]
        star = "*" if p < 0.05 else ""
        ax.text(j, i, f"{v:.2f}{star}", ha="center", va="center", fontsize=6)
ax.set_title("T-cell subtypes x TAM subtypes", fontsize=11)
fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04, label="spatial cross-correlation")

fig.suptitle("Squidpy spatial co-localization — TNBC V2 (Fig3 immune-exclusion, neighborhood level)",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{OUT}/fig_squidpy_coloc.png", dpi=150)
fig.savefig(f"{OUT}/fig_squidpy_coloc.pdf")
print("[COLOC] figure saved -> fig_squidpy_coloc.png/.pdf")
print("[COLOC] DONE")
