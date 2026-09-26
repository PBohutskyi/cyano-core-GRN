#!/usr/bin/env python3
"""
08_figure3B_table.py

Figure 3B table. Rows are Louvain clusters, ordered by the IC of the
highest-ranked regulator in each; the cluster with no module (0, shown as n/a)
is the last row. A cluster holding several regulators is written as
"Lead (other, other)" with the IC values in the same order. COG letters are
not in the inputs and are set below, keyed by one regulator in each cluster.

Inputs : outputs/tables/centralities_core.csv   (05_centralities.R)
         curation/core_regulator_curation.csv
         curation/S6_source_annotation.xlsx     Louvain cluster of each regulator
Output : outputs/tables/figure3B_table.xlsx
"""

import os
import sys

import pandas as pd

# COG letters, keyed by a regulator that sits in the cluster.
COG = {"CyAbrB2": "C", "DnaA": "F", "XylR": "G", "NtcA": "M",
       "RpaB": "C, F", "BolA": "J", "TetR": "C, F", "NusG": "J",
       "SigA1": "J", "NusA": "J", "NtcB": "K"}

CENT    = os.path.join("outputs", "tables", "centralities_core.csv")
CURCORE = os.path.join("curation", "core_regulator_curation.csv")
S6SRC   = os.path.join("curation", "S6_source_annotation.xlsx")
TAB   = os.path.join("outputs", "tables")
SHEET = "S6-2. Transcription regulators"


def main():
    for p in (CENT, CURCORE, S6SRC):
        if not os.path.exists(p):
            sys.exit("Input not found: %s" % p)
    os.makedirs(TAB, exist_ok=True)

    cen = pd.read_csv(CENT)[["locus_tag", "IC", "IC_rank"]]
    cur = pd.read_csv(CURCORE)
    cur["Stress_coupled"] = cur.Stress_category.fillna("").str.strip().ne("").map({True: "Y", False: "N"})
    lou = (pd.read_excel(S6SRC, sheet_name=SHEET)[["locus_tag", "Louvain_cluster"]])
    d = (cur[["locus_tag", "TF_name", "Stress_coupled"]]
           .merge(cen, on="locus_tag").merge(lou, on="locus_tag"))
    if len(d) != 38:
        sys.exit("expected 38 core regulators after joining inputs, got %d" % len(d))
    need = {"TF_name", "locus_tag", "Louvain_cluster", "IC", "IC_rank"}
    missing = need - set(d.columns)
    if missing:
        sys.exit("Sheet %s is missing columns: %s" % (SHEET, ", ".join(missing)))
    print("read centralities, curation and Louvain clusters: %d regulators" % len(d))

    cog_by_cluster = {}
    for tf, letter in COG.items():
        hit = d.loc[d.TF_name == tf, "Louvain_cluster"]
        if len(hit):
            cog_by_cluster[int(hit.iloc[0])] = letter
        else:
            print("NOTE: %s not found, its COG letter is unused" % tf)

    rows = []
    for cl, grp in d.groupby("Louvain_cluster"):
        grp = grp.sort_values("IC", ascending=False)
        lead, others = grp.iloc[0], grp.iloc[1:]
        name = lead.TF_name if others.empty else "%s (%s)" % (
            lead.TF_name, ", ".join(others.TF_name))
        ic = "%.2f" % lead.IC if others.empty else "%.2f (%s)" % (
            lead.IC, ", ".join("%.2f" % v for v in others.IC))
        rows.append({"Cluster": "n/a" if cl == 0 else int(cl),
                     "COG": cog_by_cluster.get(int(cl), "-"),
                     "Regulator": name,
                     "Integrated_Centrality": ic,
                     "_lead_ic": float(lead.IC),
                     "_is_na": cl == 0})

    tab = pd.DataFrame(rows)
    tab = tab.sort_values(["_is_na", "_lead_ic"], ascending=[True, False])
    tab = tab.drop(columns=["_lead_ic", "_is_na"])
    tab.insert(0, "Row", range(1, len(tab) + 1))

    out = os.path.join(TAB, "figure3B_table.xlsx")
    with pd.ExcelWriter(out) as writer:
        tab.to_excel(writer, sheet_name="Figure 3B table", index=False)
        (d[["TF_name", "locus_tag", "Louvain_cluster", "Stress_coupled",
            "IC", "IC_rank"]].sort_values("IC_rank")
           .to_excel(writer, sheet_name="Per regulator (S6-2)", index=False))

    print(tab.to_string(index=False))
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
