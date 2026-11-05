# =============================================================================
# TNBC interface-niche code repository
# Original source: project_singlecell/TNBC_result/final_figures/V3/scripts/figure3/squidpy_check/run_nhood.py
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

#!/usr/bin/env python3
"""Independent re-implementation of Figure 3's neighbourhood claim with Squidpy.

The figure's panel e is a hand-written permutation test: among T-cell-enriched
spots, what fraction of the six nearest neighbours is enriched for X, against
1000 within-section label shuffles.  This script asks the same question with
`sq.gr.nhood_enrichment` on the same Q75 labels, so that agreement (or not)
cannot be attributed to our own code.

Two label columns, each mutually exclusive by construction (a spot belongs to
exactly one class), which is what nhood_enrichment requires:

    coloc   neither / C1Q only / MMP11 only / both      -> C1Q-TAM vs MMP11-CAF
    tspp1   neither / T only / SPP1 only / both         -> SPP1-TAM vs T cells

`sample` is passed as library_key, so both the graph and the label shuffles are
built within each section — the same protection our own test uses.  (Squidpy
1.6.1 shuffles per library via _shuffle_group; without library_key it would
shuffle globally, which is the wrong null for 12 sections with different
composition.)

Two graph constructions are run, because Visium has two defensible ones:

    grid     squidpy's Visium hex grid, n_rings=1  (canonical for Visium)
    generic  kNN, n_neighs=6                       (what our W uses)

Outputs (V3/panels/figure3/squidpy_check/):
    nhood_z_<graph>_<labels>.csv      z-scores, rows = source, cols = neighbour
    nhood_count_<graph>_<labels>.csv  observed edge counts
    nhood_label_counts.csv            spots per class

    /path/to/data/project_singlecell/apps/miniconda3/envs/Squidpy/bin/python run_nhood.py
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq

HERE = os.path.dirname(os.path.abspath(__file__))
V3 = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
SPOTS = os.path.join(V3, 'panels', 'figure3', 'fig3_spots.csv')
OUT = os.path.join(V3, 'panels', 'figure3', 'squidpy_check')
N_PERMS = 1000
SEED = 42

CLASSES = {
    'coloc': ('c1q', 'mmp11', ['C1Q only', 'MMP11 only']),
    'tspp1': ('tcell', 'spp1', ['T only', 'SPP1 only']),
}


def build():
    spots = pd.read_csv(SPOTS)
    # .values matters: a Series keeps its RangeIndex and pandas would align
    # it against the barcode index, leaving every sample NaN
    obs = pd.DataFrame({'sample': pd.Categorical(spots['sample'].values)},
                       index=spots['barcode'].values)
    q75 = {}
    for name, (k1, k2, _) in CLASSES.items():
        t1, t2 = np.quantile(spots[k1], 0.75), np.quantile(spots[k2], 0.75)
        q75[name] = (t1, t2)
        a = spots[k1].values >= t1
        b = spots[k2].values >= t2
        lab = np.full(len(spots), 'neither', dtype=object)
        lab[a & ~b] = CLASSES[name][2][0]
        lab[b & ~a] = CLASSES[name][2][1]
        lab[a & b] = 'both'
        obs[name] = pd.Categorical(lab)
        print(f'[labels] {name}: thr {k1}={t1:.4f} {k2}={t2:.4f} ->',
              dict(pd.Series(lab).value_counts()))

    A = ad.AnnData(X=np.zeros((len(spots), 1), dtype=np.float32), obs=obs)
    A.obsm['spatial'] = spots[['x', 'y']].values.astype(np.float64)
    return A, q75


def main():
    os.makedirs(OUT, exist_ok=True)
    A, _ = build()
    A.obs['sample'] = A.obs['sample'].cat.remove_unused_categories()

    labels = pd.concat([A.obs['sample'], A.obs['coloc'], A.obs['tspp1']], axis=1)
    labels.to_csv(f'{OUT}/nhood_label_counts.csv')
    print('[spots]', A.n_obs, '| sections', A.obs['sample'].nunique())

    for graph, kw in [('grid', dict(coord_type='grid', n_rings=1)),
                      ('generic', dict(coord_type='generic', n_neighs=6))]:
        # rebuild the graph for every variant: spatial_neighbors overwrites the
        # key, and the two variants must not leak into each other
        sq.gr.spatial_neighbors(A, library_key='sample', **kw)
        C = A.obsp['spatial_connectivities']
        print(f'[graph] {graph}: {C.shape}, nnz={C.nnz}, '
              f'mean degree={C.getnnz(axis=1).mean():.2f}')

        for name in CLASSES:
            sq.gr.nhood_enrichment(A, cluster_key=name, library_key='sample',
                                   n_perms=N_PERMS, seed=SEED, n_jobs=4,
                                   show_progress_bar=False)
            z = pd.DataFrame(A.uns[f'{name}_nhood_enrichment']['zscore'],
                             index=A.obs[name].cat.categories,
                             columns=A.obs[name].cat.categories)
            cnt = pd.DataFrame(A.uns[f'{name}_nhood_enrichment']['count'],
                               index=z.index, columns=z.columns)
            z.to_csv(f'{OUT}/nhood_z_{graph}_{name}.csv')
            cnt.to_csv(f'{OUT}/nhood_count_{graph}_{name}.csv')
            print(f'[z] {graph}/{name}')
            print(z.round(2).to_string())


if __name__ == '__main__':
    main()
