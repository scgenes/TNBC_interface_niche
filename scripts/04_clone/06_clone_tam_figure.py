# =============================================================================
# TNBC interface-niche code repository
# Original source: project_spatial/TNBC/patient_ALL/TAM_clone/scripts/fig_mainA_C1Q.py
# Inputs/outputs and run order: see scripts/<subdir>/README.md
# Paths: /path/to/project and /path/to/data are placeholders; update
#        them to your local roots before running.
# =============================================================================

import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

STEP = Path("/path/to/project/project_spatial/TNBC/patient_ALL/TAM_clone")
TB = STEP / "tables"
FIG = STEP / "figures"
FIG.mkdir(exist_ok=True)

# ---- nature palette (ggsci NPG) ----
RED, BLUE = "#E64B35", "#4DBBD5"

def num(s):
    try:
        v = float(s)
        return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan

recs = []
def push(r, zone, method, roe_key, pe_key, pd_key):
    recs.append(dict(patient=r["patient"], slice=r["slice"], clone=r["clone_id"], zone=zone,
                     method=method, roe=num(r[roe_key]), pE=num(r[pe_key]), pD=num(r[pd_key]),
                     est=str(r["estimable"]).strip().upper() == "TRUE"))

for r in pd.read_csv(TB / "clone_tam_radius_perm_by_slice.csv", keep_default_na=False).to_dict("records"):
    if r["tam"] == "C1Q TAMs":
        push(r, "All", "A radius", "enrichment", "p_perm", "p_deplete")
for r in pd.read_csv(TB / "clone_tam_radius_perm_region_stratified.csv", keep_default_na=False).to_dict("records"):
    if r["tam"] == "C1Q TAMs":
        push(r, r["region"], "A radius", "enrichment", "p_perm", "p_deplete")
for r in pd.read_csv(TB / "clone_tam_coloc_samespot.csv", keep_default_na=False).to_dict("records"):
    if r["tam"] == "C1Q TAMs":
        push(r, r["zone"], "B same spot", "roe", "p_enrich", "p_deplete")
for r in pd.read_csv(TB / "clone_tam_coloc_neighbor.csv", keep_default_na=False).to_dict("records"):
    if r["tam"] == "C1Q TAMs" and num(r["k"]) == 6:
        push(r, r["zone"], "B k=6", "roe", "p_enrich", "p_deplete")

df = pd.DataFrame(recs)
df = df[df.patient.isin(["P13", "P15"])].copy()
df["clone_short"] = df["clone"].str.split("_").str[1]

def direction(row):
    if not row["est"]:
        return "not estimable"
    if row["roe"] > 1 and row["pE"] < 0.05:
        return "enriched"
    if row["roe"] < 1 and row["pD"] < 0.05:
        return "depleted"
    return "not significant"

df["direction"] = df.apply(direction, axis=1)
df["significant"] = df["direction"].isin(["enriched", "depleted"])

out = df[["patient", "slice", "clone", "zone", "method", "roe", "pE", "pD", "est", "direction", "significant"]].copy()
out.columns = ["patient", "slice", "clone", "zone", "method", "Ro_e", "p_enrich", "p_deplete",
               "estimable", "direction", "significant"]
out.to_csv(TB / "mainA_C1Q_plot_values.csv", index=False)

# ---- plotting ----
PATIENTS = [("P13", ["094A", "094B"]), ("P15", ["095A", "095B"])]
ZONES = ["All", "Tumor Core", "Interface"]
METHODS = ["A radius", "B same spot", "B k=6"]
MK = {"A radius": "o", "B same spot": "s", "B k=6": "^"}
MOFF = {"A radius": -0.055, "B same spot": 0.0, "B k=6": 0.055}
SLICE_OFF = {"0": -0.16, "1": 0.16}
CLONES = ["Clone1", "Clone2"]
CLONE_LAB = {"Clone1": "Clone1\n(CNV-low)", "Clone2": "Clone2\n(CNV-high)"}
XL = {"Clone1": 0.0, "Clone2": 1.0}
NA_Y = 0.052

fig, axes = plt.subplots(3, 2, figsize=(12, 12), sharey=True,
                         gridspec_kw=dict(left=0.075, right=0.985, top=0.935, bottom=0.115,
                                          wspace=0.10, hspace=0.16))
fig.patch.set_facecolor("white")

for ci, (pat, slices) in enumerate(PATIENTS):
    scol = {s: (RED if i == 0 else BLUE) for i, s in enumerate(slices)}
    for ri, zone in enumerate(ZONES):
        ax = axes[ri, ci]
        ax.set_facecolor("white")
        ax.axhline(1.0, color="#555555", lw=0.9, ls=(0, (4, 3)), zorder=1)

        for clone in CLONES:
            for si, sl in enumerate(slices):
                xc = XL[clone] + SLICE_OFF[str(si)]
                sub = df[(df.patient == pat) & (df.slice == sl) & (df.clone_short == clone) & (df.zone == zone)]
                if sub.empty:
                    continue
                plotted = []
                for m in METHODS:
                    r = sub[sub.method == m]
                    if r.empty:
                        continue
                    r = r.iloc[0]
                    x = xc + MOFF[m]
                    if r["est"] and np.isfinite(r["roe"]):
                        plotted.append((x, r["roe"], r["significant"], m))
                    else:
                        ax.plot(x, NA_Y, marker=MK[m], ms=4.2, mfc="none", mec="#9A9A9A",
                                mew=0.9, ls="none", zorder=3)
                ys = [p[1] for p in plotted]
                if len(ys) > 1:
                    ax.plot([np.mean([p[0] for p in plotted])] * 2, [min(ys), max(ys)],
                            color=scol[sl], lw=0.8, alpha=0.30, zorder=2, solid_capstyle="round")
                for x, y, sg, m in plotted:
                    ax.plot(x, y, marker=MK[m], ms=7.0, ls="none", zorder=4,
                            mfc=scol[sl] if sg else "white", mec=scol[sl], mew=1.4)

        ax.set_yscale("log")
        ax.set_ylim(0.042, 6.5)
        ax.set_yticks([0.05, 0.1, 0.25, 0.5, 1, 2, 4])
        ax.set_yticklabels(["0.05", "0.1", "0.25", "0.5", "1", "2", "4"], fontsize=9)
        ax.minorticks_off()
        ax.set_xlim(-0.42, 1.42)
        ax.set_xticks([XL[c] for c in CLONES])
        if ri == 2:
            ax.set_xticklabels([CLONE_LAB[c] for c in CLONES], fontsize=9.5)
        else:
            ax.set_xticklabels([])
        for s in ax.spines.values():
            s.set_linewidth(0.9); s.set_color("#222222")
        ax.tick_params(axis="both", width=0.8, length=3.5, colors="#222222")

        if ri == 0:
            ax.set_title(f"{pat}   ({' / '.join(slices)})", fontsize=13, pad=9, color="#111111")
        if ci == 0:
            ax.set_ylabel(f"{zone}\nRo/e  (C1Q+ TAM)", fontsize=10.5)
        else:
            ax.text(0.985, 0.965, zone, transform=ax.transAxes, ha="right", va="top",
                    fontsize=8.5, color="#666666")

        if ri == 0:
            handles = [Line2D([], [], marker="o", ls="none", ms=6.5, mfc=scol[sl], mec=scol[sl],
                              label=f"{sl}") for sl in slices]
            ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 1.0), frameon=False,
                      fontsize=8.5, handletextpad=0.35, borderpad=0.1, labelspacing=0.25)

fig.text(0.5, 0.075,
         "Ro/e = observed / expected overlap with C1Q+ TAMs.   "
         "Filled = p < 0.05 (within-slice permutation);  open = not significant;  "
         "grey marker near axis = not estimable.   Each slice tested independently.",
         ha="center", va="center", fontsize=8.6, color="#333333")

legend_elems = [
    Line2D([], [], marker="o", ls="none", ms=7, mfc="#444444", mec="#444444", label="Method A: radius-neighbourhood"),
    Line2D([], [], marker="s", ls="none", ms=7, mfc="#444444", mec="#444444", label="Method B: same spot"),
    Line2D([], [], marker="^", ls="none", ms=7, mfc="#444444", mec="#444444", label="Method B: k = 6 neighbours"),
    Line2D([], [], marker="o", ls="none", ms=7, mfc="#444444", mec="#444444", label="p < 0.05"),
    Line2D([], [], marker="o", ls="none", ms=7, mfc="white", mec="#444444", mew=1.4, label="not significant"),
    Line2D([], [], marker="o", ls="none", ms=5, mfc="none", mec="#9A9A9A", label="not estimable"),
]
fig.legend(handles=legend_elems, loc="lower center", bbox_to_anchor=(0.5, 0.018),
           ncol=6, frameon=False, fontsize=8.4, handletextpad=0.3, columnspacing=1.5)

stem = FIG / "mainA_C1Q_clone_niche"
fig.savefig(str(stem) + ".png", dpi=300, facecolor="white", transparent=False, bbox_inches="tight")
fig.savefig(str(stem) + ".pdf", facecolor="white", transparent=False, bbox_inches="tight")
print("saved:", stem)
print(df.groupby(["patient", "zone"]).size().to_string())
