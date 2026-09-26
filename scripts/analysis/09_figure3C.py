#!/usr/bin/env python3
"""
09_figure3C.py

Figure 3C: dendrogram of the core GRN regulators by Jaccard distance between
their target sets (average linkage), with their centrality values below it.
Regulators sharing no target with any other regulator have distance 1 to all
others and are left out (7 of 38); the script prints which.

Inputs : networks/core/network_of_record.graphml
         outputs/tables/centralities_core.csv   (05_centralities.R)
         curation/core_regulator_curation.csv
Outputs: outputs/figures/fig3C_assembled.{png,pdf}
         outputs/tables/fig3C_values.csv       leaf order, colour group, values
         outputs/tables/fig3C_leaf_order.txt
"""

import os
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage, set_link_color_palette

# ---- options ---------------------------------------------------------------
ALL_REGULATORS = False        # True: all 38. False: regulators sharing a target.

# ---- style ----------------------------------------------------------------
FONT       = "Times New Roman"
PALETTE    = ["c", "m", "y", "k", "g", "b", "purple"]
ABOVE_COL  = "red"
CUT_H      = 0.99
FIGSIZE    = (6, 8)
BAR_ALPHA  = 0.8
BAR_EDGE   = "black"
BAR_LW     = 1
BAR_WIDTH  = 0.8
YLAB_SIZE  = 12
XLAB_SIZE  = 10
LOCUS_STRIP = "SYNPCC7942_"   # dropped from the tick labels
DPI        = 600

# column, y-axis label
ROWS = [("degree",      "Degree"),
        ("betweenness", "BC"),
        ("eigenvector", "EC"),
        ("k_core",      "K-Core"),
        ("stress",      "Stress"),
        ("IC",          "IV")]

GRAPHML = os.path.join("networks", "core", "network_of_record.graphml")
CENT    = os.path.join("outputs", "tables", "centralities_core.csv")
CURCORE = os.path.join("curation", "core_regulator_curation.csv")
FIGDIR  = os.path.join("outputs", "figures")
TABDIR  = os.path.join("outputs", "tables")


def jaccard(a, b):
    return 1 - len(a & b) / len(a | b)


def main():
    for path in (GRAPHML, CENT, CURCORE):
        if not os.path.exists(path):
            sys.exit("Input not found: %s" % path)
    os.makedirs(FIGDIR, exist_ok=True)
    os.makedirs(TABDIR, exist_ok=True)

    G = nx.read_graphml(GRAPHML)
    tfs = [n for n, d in G.nodes(data=True) if str(d.get("TF")) == "1"]
    print("read %s: %d nodes, %d edges, %d regulators"
          % (os.path.basename(GRAPHML), G.number_of_nodes(),
             G.number_of_edges(), len(tfs)))

    cen = pd.read_csv(CENT).rename(columns={"kcore": "k_core"})
    cur = pd.read_csv(CURCORE)[["locus_tag", "TF_name"]]
    s6 = cen.drop(columns=["TF_name"], errors="ignore").merge(cur, on="locus_tag")
    if len(s6) != len(cen):
        sys.exit("centralities_core.csv and the core curation list different regulators")
    name = dict(zip(s6.locus_tag, s6.TF_name))

    targets = {t: set(G.successors(t)) for t in tfs}
    if ALL_REGULATORS:
        keep = list(tfs)
        dropped = []
    else:
        non_tf = {t: {u for u in v if u not in set(tfs)}
                  for t, v in targets.items()}
        counts = Counter(u for v in non_tf.values() for u in v)
        shared = {u for u, c in counts.items() if c > 1}
        keep = [t for t in tfs if non_tf[t] & shared]
        dropped = [t for t in tfs if t not in keep]
    print("%d regulators clustered, %d dropped for having no shared target: %s"
          % (len(keep), len(dropped),
             ", ".join(sorted(name.get(t, t) for t in dropped)) or "none"))

    labels = [name.get(t, t) for t in keep]
    sets = {name.get(t, t): targets[t] for t in keep}

    n = len(labels)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = jaccard(sets[labels[i]], sets[labels[j]])

    linked = linkage(squareform(D), method="average", optimal_ordering=True)

    matplotlib.rcParams["font.family"] = FONT
    if FONT not in {f.name for f in matplotlib.font_manager.fontManager.ttflist}:
        print("NOTE: '%s' not found, matplotlib will substitute a serif font."
              % FONT)

    fig = plt.figure(figsize=FIGSIZE)
    gs = matplotlib.gridspec.GridSpec(nrows=len(ROWS) + 1, ncols=1, wspace=-1)

    ax_dendro = fig.add_subplot(gs[0])
    set_link_color_palette(PALETTE)
    dend = dendrogram(linked, labels=labels, orientation="top",
                      distance_sort="ascending", ax=ax_dendro,
                      color_threshold=CUT_H, above_threshold_color=ABOVE_COL)
    ax_dendro.axis("off")

    order = dend["ivl"]
    colours = dend["leaves_color_list"]
    print("leaf order: " + ", ".join(order))

    df = s6.set_index("TF_name").loc[order].reset_index()
    x = np.arange(len(order))
    ticklab = ["%s (%s)" % (r.TF_name, str(r.locus_tag).replace(LOCUS_STRIP, ""))
               for r in df.itertuples()]

    for k, (col, ylab) in enumerate(ROWS):
        ax = fig.add_subplot(gs[k + 1])
        ax.bar(x, df[col].astype(float).values, width=BAR_WIDTH,
               color=colours, alpha=BAR_ALPHA, edgecolor=BAR_EDGE,
               linewidth=BAR_LW)
        ax.set_ylabel(ylab, fontsize=YLAB_SIZE, weight="bold", labelpad=14)
        ax.set_xlim(-0.6, len(order) - 0.4)
        ax.set_xlabel("")
        if k < len(ROWS) - 1:
            ax.set_xticks([])
        else:
            ax.set_xticks(x)
            ax.set_xticklabels(ticklab, rotation=90, fontsize=XLAB_SIZE,
                               weight="bold")

    plt.tight_layout()
    stem = os.path.join(FIGDIR, "fig3C_assembled")
    fig.savefig(stem + ".png", format="png", bbox_inches="tight", dpi=DPI)
    fig.savefig(stem + ".pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)

    out = df[["TF_name", "locus_tag"] + [c for c, _ in ROWS]].copy()
    out.insert(0, "order", np.arange(1, len(order) + 1))
    out["colour_group"] = colours
    out.to_csv(os.path.join(TABDIR, "fig3C_values.csv"), index=False)
    with open(os.path.join(TABDIR, "fig3C_leaf_order.txt"), "w") as fh:
        fh.write(", ".join(order) + "\n")

    print("wrote %s.png and .pdf" % stem)
    print("wrote %s" % os.path.join(TABDIR, "fig3C_values.csv"))
    print("wrote %s" % os.path.join(TABDIR, "fig3C_leaf_order.txt"))


if __name__ == "__main__":
    main()
