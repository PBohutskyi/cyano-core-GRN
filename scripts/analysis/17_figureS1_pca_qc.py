#!/usr/bin/env python3
"""
17_figureS1_pca_qc.py

Global quality check of the core expression compendium (SynCOREexpress):
principal component analysis, association of each PC with species, BioProject
and condition, and an assessment of BioProject (batch) against biological
structure. Supports Supplementary Figure S1 and Methods 2.2.

Method: log2 TPM, genes z-scored, PCA with 10 components (scikit-learn), no
batch correction. eta2 = share of a PC's variance explained by a one-way
grouping; permuted null = mean eta2 over 200 label permutations. Distances are
Euclidean in PC1 to PC10.

Inputs : expression/expression_core.csv                   (genes x samples)
         datasets/Dataset_S2_PostQC_RNAseq_samples.xlsx   (SRA, BioProject, Condition)
Outputs: tables/figS1_pca_scores.csv          PC1-PC10 per sample, with species,
                                              BioProject and condition
         tables/figS1_pc_associations.csv    variance explained and eta2 per PC
         tables/figS1_batch_assessment.csv   nearest-neighbour, silhouette and
                                              distance metrics
         figures/figS1_pca_pairwise.pdf / .png
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

EXPR = os.path.join("expression", "expression_core.csv")
META = os.path.join("datasets", "Dataset_S2_PostQC_RNAseq_samples.xlsx")
META_SHEETS = {"PCC 7942": "S2-1. PCC 7942",
               "PCC 7002": "S2-2. PCC 7002",
               "PCC 6803": "S2-3. PCC 6803"}
TAB = os.path.join("outputs", "tables")
FIG = os.path.join("outputs", "figures")

SPECIES_ORDER = ["PCC 7942", "PCC 6803", "PCC 7002"]
# same strain colours as Figures 1, 2 and S2 (02_figure1A.R, 04_figure2.R, 13_figureS2)
SPECIES_COLORS = {"PCC 7942": "#D55E00", "PCC 6803": "#009E73", "PCC 7002": "#0072B2"}
N_PC = 10
N_PERM = 200
K_NN = 15
SEED = 0


def load():
    for f in (EXPR, META):
        if not os.path.exists(f):
            sys.exit("Input not found: %s" % f)
    expr = pd.read_csv(EXPR, index_col=0)
    # two split-run samples are written SRRxxx.1 in the matrix and SRRxxx_1 in Dataset S2
    expr.columns = [c.replace(".", "_") if c.startswith(("SRR12129543", "SRR12129544")) else c
                    for c in expr.columns]
    frames = []
    for sp, sheet in META_SHEETS.items():
        d = pd.read_excel(META, sheet_name=sheet)[["SRA", "BioProject", "Condition"]].copy()
        d["species"] = sp
        frames.append(d)
    meta = pd.concat(frames, ignore_index=True).set_index("SRA")
    missing = set(expr.columns) ^ set(meta.index)
    if missing:
        sys.exit("STOP: samples differ between matrix and Dataset S2: %s" % sorted(missing)[:10])
    meta = meta.loc[expr.columns]
    if expr.isna().any().any():
        sys.exit("STOP: missing values in the expression matrix")
    print("core genes x samples: %d x %d" % expr.shape)
    print(meta.species.value_counts().to_string())
    return expr, meta


def eta2(y, g):
    y = np.asarray(y, dtype=float)
    gm = pd.Series(y).groupby(np.asarray(g, dtype=object)).transform("mean").values
    return ((gm - y.mean()) ** 2).sum() / ((y - y.mean()) ** 2).sum()


def associations(scores, ve, meta, expr):
    rng = np.random.default_rng(SEED)
    mean_expr = expr.mean(axis=0).values
    rows = []
    for i in range(5):
        y = scores[:, i]
        r = {"PC": "PC%d" % (i + 1),
             "Variance_explained_percent": round(ve[i], 2),
             "Cumulative_percent": round(ve[: i + 1].sum(), 2)}
        for col, lab in (("species", "Species"), ("BioProject", "BioProject"),
                         ("Condition", "Condition")):
            g = meta[col].values
            r[lab + "_eta2"] = round(eta2(y, g), 3)
            r[lab + "_eta2_permuted_null"] = round(
                np.mean([eta2(y, rng.permutation(g)) for _ in range(N_PERM)]), 3)
        w = 0.0
        for _, idx in meta.groupby("species").indices.items():
            w += len(idx) * eta2(y[idx], meta.Condition.values[idx])
        r["Within_species_Condition_eta2"] = round(w / len(y), 3)
        r["r_with_mean_log2TPM"] = round(float(np.corrcoef(y, mean_expr)[0, 1]), 3)
        rows.append(r)
    return pd.DataFrame(rows)


def batch_assessment(scores, meta):
    P = scores[:, :N_PC]
    sp = np.asarray(meta.species.values, dtype=object)
    bp = np.asarray(meta.BioProject.values, dtype=object)
    cond = np.asarray(meta.Condition.values, dtype=object)

    _, ind = NearestNeighbors(n_neighbors=K_NN + 1).fit(P).kneighbors(P)
    ind = ind[:, 1:]
    same_sp = sp[ind] == sp[:, None]
    diff_bp = same_sp & (bp[ind] != bp[:, None])

    D = squareform(pdist(P))
    upper = np.triu(np.ones_like(D, dtype=bool), 1)
    SS = sp[:, None] == sp[None, :]
    SC = cond[:, None] == cond[None, :]
    SB = bp[:, None] == bp[None, :]
    multi = {(s, c) for (s, c), g in meta.groupby(["species", "Condition"])
             if g.BioProject.nunique() >= 2}
    in_multi = np.array([(s, c) in multi for s, c in zip(sp, cond)])
    M = in_multi[:, None]

    d_within = D[upper & SS & SC & SB & M].mean()
    d_across = D[upper & SS & SC & ~SB & M].mean()
    d_random = D[upper & SS].mean()
    rows = [
        ("Fraction of nearest neighbours (k=%d) from the same species" % K_NN, same_sp.mean()),
        ("Fraction of same-species nearest neighbours from a different BioProject",
         diff_bp.sum() / same_sp.sum()),
        ("Silhouette width by species (PC1-PC%d)" % N_PC, silhouette_score(P, sp)),
        ("Mean PC distance: same condition, same BioProject", d_within),
        ("Mean PC distance: same condition, different BioProject", d_across),
        ("Ratio, cross-BioProject / within-BioProject (same condition)", d_across / d_within),
        ("Mean PC distance: same-species pairs", d_random),
        ("Cross-BioProject distance as a fraction of same-species distance", d_across / d_random),
        ("Condition labels represented in two or more BioProjects within a species", len(multi)),
    ]
    return pd.DataFrame(rows, columns=["Metric", "Value"]).round({"Value": 3})


def figure(pc, ve, meta):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "axes.linewidth": 0.6, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 300})
    bps = sorted(meta.BioProject.unique())
    cmap = plt.cm.gist_ncar(np.linspace(0, 1, len(bps)))
    bp_col = [dict(zip(bps, cmap))[b] for b in meta.BioProject]
    pairs = [(1, 2), (1, 3), (2, 3), (3, 4), (3, 5), (4, 5)]
    fig, axes = plt.subplots(3, 4, figsize=(7.5, 6.6))
    li = 0
    for k, (a, b) in enumerate(pairs):
        r, c = divmod(k, 2)
        ax_sp, ax_bp = axes[r, 2 * c], axes[r, 2 * c + 1]
        for s in SPECIES_ORDER:
            m = (meta.species == s).values
            ax_sp.scatter(pc[m, a - 1], pc[m, b - 1], s=2.0, c=SPECIES_COLORS[s],
                          alpha=0.55, linewidths=0, rasterized=True)
        ax_bp.scatter(pc[:, a - 1], pc[:, b - 1], s=2.0, c=bp_col, alpha=0.55,
                      linewidths=0, rasterized=True)
        for ax, lab in ((ax_sp, "species"), (ax_bp, "BioProject")):
            ax.set_xlabel("PC%d (%.1f%%)" % (a, ve[a - 1]), fontsize=6.5, labelpad=1.5)
            ax.set_ylabel("PC%d (%.1f%%)" % (b, ve[b - 1]), fontsize=6.5, labelpad=1.5)
            ax.set_title(lab, fontsize=6.8, loc="left", pad=2.5)
            ax.tick_params(labelsize=5.8, length=2, pad=1)
            ax.text(-0.30, 1.16, "ABCDEFGHIJKL"[li], transform=ax.transAxes,
                    fontsize=10, fontweight="bold", va="top")
            li += 1
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=SPECIES_COLORS[s],
                      markersize=4.5, label="%s (n=%d)" % (s, (meta.species == s).sum()))
               for s in SPECIES_ORDER]
    handles.append(Line2D([0], [0], marker="o", color="none", markerfacecolor="0.55",
                          markersize=4.5,
                          label="BioProject panels: %d projects, one color each" % len(bps)))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=7,
               handletextpad=0.3, columnspacing=1.4, bbox_to_anchor=(0.5, -0.012))
    fig.tight_layout(w_pad=1.4, h_pad=2.0, rect=[0, 0.035, 1, 1])
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, "figS1_pca_pairwise." + ext), bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    expr, meta = load()
    X = StandardScaler().fit_transform(expr.T.values)
    pca = PCA(n_components=N_PC, random_state=SEED)
    pc = pca.fit_transform(X)
    ve = pca.explained_variance_ratio_ * 100
    print("variance explained PC1-PC5 (%%): %s" % ", ".join("%.2f" % v for v in ve[:5]))

    scores = pd.DataFrame(pc, index=meta.index, columns=["PC%d" % (i + 1) for i in range(N_PC)])
    scores = meta[["species", "BioProject", "Condition"]].join(scores)
    scores.index.name = "SRA"
    scores.round(4).to_csv(os.path.join(TAB, "figS1_pca_scores.csv"))

    assoc = associations(pc, ve, meta, expr)
    assoc.to_csv(os.path.join(TAB, "figS1_pc_associations.csv"), index=False)
    print(assoc.to_string(index=False))

    batch = batch_assessment(pc, meta)
    batch.to_csv(os.path.join(TAB, "figS1_batch_assessment.csv"), index=False)
    print(batch.to_string(index=False))
    print("group counts: species %d, BioProject %d, Condition %d"
          % tuple(meta[c].nunique() for c in ("species", "BioProject", "Condition")))

    figure(pc, ve, meta)
    print("wrote tables/figS1_*.csv and figures/figS1_pca_pairwise.pdf/.png")


if __name__ == "__main__":
    main()
