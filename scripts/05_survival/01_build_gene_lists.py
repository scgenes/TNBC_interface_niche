#!/usr/bin/env python3
"""Build the 20-gene bulk signatures used in the survival analysis (Figure S6).

Inputs are the two differential-expression tables exported from the single-cell
analysis (Mac_sub and Fib_sub). Each table has columns:
    genes, log2FC, cluster, pct, avgExpr

For each target cluster we take the top 20 genes by log2FC after removing
technical / sex / mitochondrial genes and immunoglobulin variable-region genes.

Original source: project_spatial/TNBC/survival/processed/signature_from_DE
(ad-hoc extraction); this script reproduces the delivered gene lists.

Usage
-----
python 01_build_gene_lists.py \
    --mac  /path/to/Mac_sub_差异列表.tsv \
    --fib  /path/to/Fib_sub_差异列表.tsv \
    --out  ../../data/signatures
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

# Genes removed before ranking (same rule as the delivered signatures).
TECH_RE = re.compile(r"^(MT-|RPL|RPS|HB[AB]|MALAT1$|XIST$|RPS4Y1$|DDX3Y$)")
IG_RE = re.compile(r"^(IGKV|IGLV)")

TARGETS = [
    ("Mac_sub_差异列表.tsv", ["C1Q TAMs", "SPP1 TAMs"]),
    ("Fib_sub_差异列表.tsv", ["MMP11 CAFs"]),
]


def build_one(df: pd.DataFrame, cluster: str, source: str, top_n: int = 20) -> tuple[pd.DataFrame, dict]:
    x = df[df["cluster"] == cluster].copy()
    x["excluded_tech"] = x["genes"].astype(str).str.match(TECH_RE)
    x["excluded_ig"] = x["genes"].astype(str).str.match(IG_RE)
    keep = x[~(x["excluded_tech"] | x["excluded_ig"])]
    kept = (
        keep.sort_values(["log2FC", "pct"], ascending=[False, False])
        .head(top_n)
        .reset_index(drop=True)
    )
    kept["rank"] = range(1, len(kept) + 1)
    kept["cell_type"] = cluster
    kept["source_file"] = source
    stats = {
        "cell_type": cluster,
        "n_tested": len(x),
        "n_excluded_tech": int(x["excluded_tech"].sum()),
        "n_excluded_ig": int(x["excluded_ig"].sum()),
        "n_after_excl": len(keep),
        "n_top20_pct_ge_025": int((kept["pct"] >= 0.25).sum()),
        "top20_min_pct": round(float(kept["pct"].min()), 3),
    }
    return kept[["cell_type", "rank", "genes", "log2FC", "pct", "avgExpr", "source_file"]], stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mac", required=True, help="Mac_sub differential-expression table (TSV)")
    ap.add_argument("--fib", required=True, help="Fib_sub differential-expression table (TSV)")
    ap.add_argument("--out", required=True, help="Output directory (data/signatures)")
    ap.add_argument("--top-n", type=int, default=20)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    files = {"Mac_sub_差异列表.tsv": Path(args.mac), "Fib_sub_差异列表.tsv": Path(args.fib)}

    rows, stats = [], []
    for source, clusters in TARGETS:
        df = pd.read_csv(files[source], sep="\t")
        for cluster in clusters:
            kept, st = build_one(df, cluster, source, args.top_n)
            rows.append(kept)
            stats.append(st)

    long = pd.concat(rows, ignore_index=True)
    long.to_csv(out / "gene_list_top20_long.tsv", sep="\t", index=False)
    long.pivot(index="rank", columns="cell_type", values="genes").to_csv(
        out / "gene_list_top20_wide.csv"
    )
    unique = sorted(long["genes"].unique())
    (out / "gene_list_unique.txt").write_text("\n".join(unique) + "\n")
    pd.DataFrame(stats).to_csv(out / "gene_list_stats.csv", index=False)

    print(f"[done] {len(long)} genes across {long['cell_type'].nunique()} signatures")
    print(f"[done] unique genes: {len(unique)}")
    print(f"[out ] {out}")


if __name__ == "__main__":
    main()
