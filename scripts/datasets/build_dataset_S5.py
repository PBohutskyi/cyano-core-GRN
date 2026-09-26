#!/usr/bin/env python3
"""
build_dataset_S5.py

Builds Dataset S5 from the curated regulator workbook. Every sheet is copied;
in S5-1 to S5-3 two columns are added, name_source and functional_category,
filled for the regulators shown in Figure 6 and blank for all others. TF_name
comes from the curated workbook; if curation/figure6_categories.csv gives a
different name for the same locus, the run stops.

Inputs : curation/S5_curated_regulators.xlsx
         curation/figure6_categories.csv
         outputs/tables/figure6_symbols.csv   (14_figures5_6.py)
Output : datasets/Dataset_S5_Curated_TFs.xlsx
"""

import os
import sys

import pandas as pd

SRC = os.path.join("curation", "S5_curated_regulators.xlsx")
CUR = os.path.join("curation", "figure6_categories.csv")
SYM = os.path.join("outputs", "tables", "figure6_symbols.csv")
DST = os.path.join("datasets", "Dataset_S5_Curated_TFs.xlsx")

SHEETS = ["S5-1. Curated TFs PCC 7942",
          "S5-2. Curated TFs PCC 7002",
          "S5-3. Curated TFs PCC 6803"]


def blank(v):
    return not isinstance(v, str) or not v.strip() or v.strip().lower() == "nan"


def main():
    for p in (SRC, CUR, SYM):
        if not os.path.exists(p):
            sys.exit("Input not found: %s" % p)

    cur = pd.read_csv(CUR)
    in_fig6 = set(pd.read_csv(SYM).Species_locus)
    if set(cur.Species_locus) != in_fig6:
        sys.exit("STOP: figure6_categories.csv and figure6_symbols.csv list different regulators")
    name = dict(zip(cur.Species_locus, cur.Regulator))
    src  = dict(zip(cur.Species_locus, cur.name_source))
    cat  = dict(zip(cur.Species_locus, cur.category))

    book = pd.ExcelFile(SRC)
    out, clashes, filled = {}, [], 0
    for sheet in book.sheet_names:
        if sheet not in SHEETS and not sheet.startswith("S5-4"):
            out[sheet] = book.parse(sheet, header=None)     # legend: copied as is
            continue
        d = book.parse(sheet)
        if sheet in SHEETS:
            d["name_source"] = ""
            d["functional_category"] = ""
            for i, row in d.iterrows():
                loc = row["locus_tag"]
                if loc not in name:
                    continue
                old, new = row.get("TF_name"), name[loc]
                if not blank(new) and not blank(old) and str(old).strip() != str(new).strip():
                    clashes.append((loc, old, new))
                d.at[i, "name_source"] = "" if blank(src.get(loc)) else src[loc]
                d.at[i, "functional_category"] = "" if blank(cat.get(loc)) else cat[loc]
                filled += 1
        out[sheet] = d

    if clashes:
        print("STOP: TF_name differs between the two curation files, nothing written:")
        for loc, old, new in clashes:
            print("   %-20s S5 curation %-12s figure6_categories %s" % (loc, old, new))
        sys.exit(1)
    missing = set(name) - {l for s in SHEETS for l in out[s].locus_tag}
    if missing:
        sys.exit("STOP: Figure 6 regulators missing from S5 lists: %s" % ", ".join(sorted(missing)))

    with pd.ExcelWriter(DST) as writer:
        for sheet, d in out.items():
            d.to_excel(writer, sheet_name=sheet, index=False,
                       header=not sheet.endswith("Legend"))
    print("rows given name_source and functional_category: %d" % filled)
    for s in SHEETS:
        print("  %s: %d regulators" % (s, len(out[s])))


if __name__ == "__main__":
    main()
