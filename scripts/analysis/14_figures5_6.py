#!/usr/bin/env python3
"""
14_figures5_6.py

Symbol tables behind Figure 5 (cross-validated regulators) and Figure 6
(high-influence species regulators that are not cross-validated), with the
counts given in the captions. Cutoffs and measures are those of
13_crossvalidation.R: core top 12 and species top 10 per measure, ties broken
by IC.

A placement is one regulator in one cell (species x measure group). Figure 5
categories come from the core curation; Figure 6 names and categories from
curation/figure6_categories.csv. Every tri-homolog in Figure 6 must carry the
same category as in the core curation, and every Figure 6 regulator must be
curated; otherwise the run stops.

Shapes are assigned from the product name by shape_of(). Symbol and font sizes
are spread by rank within each figure, ties broken by the regulator's IC in its
own species network.

Inputs : outputs/tables/centralities_{core,s7942,s6803,s7002}.csv  (05_centralities.R)
         outputs/tables/ortholog_groups.csv                        (01_regulator_counts.R)
         curation/S5_curated_regulators.xlsx
         curation/core_regulator_curation.csv
         curation/figure6_categories.csv
Outputs: outputs/tables/figure5_symbols.csv
         outputs/tables/figure6_symbols.csv
         outputs/tables/figure6_regulators.csv   one row per Figure 6 regulator
         outputs/tables/figure56_counts.txt      numbers for the captions
"""

import os
import re
import sys

import numpy as np
import pandas as pd

CORE_K, SPECIES_N = 12, 10          # as in 13_crossvalidation.R
FONT_MIN, FONT_MAX = 11.0, 18.5      # font size of the regulator name
FONT_STEP = 0.5
SYM_MIN            = 1.2             # symbol size
SYM_MAX_FIG5       = 2.2
SYM_MAX_FIG6       = 2.0

# measure column in the centrality files -> normalized column, group label
MEASURES = {
    "degree":      ("degree_norm",      "local (degree/k-core)"),
    "kcore":       ("kcore_norm",       "local (degree/k-core)"),
    "betweenness": ("betweenness_norm", "global (betweenness/stress)"),
    "stress":      ("stress_norm",      "global (betweenness/stress)"),
    "eigenvector": ("eigenvector_norm", "community (eigenvector)"),
    "IC":          ("IC",               "integrated (IC)"),
}
GROUP_ORDER = ["local (degree/k-core)", "global (betweenness/stress)",
               "community (eigenvector)", "integrated (IC)"]

SPECIES = {                          # label: (centrality file, S5-4 column, S5 sheet)
    "S. elongatus PCC 7942":
        ("s7942", "PCC7942_locus_tag", "S5-1. Curated TFs PCC 7942"),
    "Synechocystis sp. PCC 6803":
        ("s6803", "PCC6803_locus_tag", "S5-3. Curated TFs PCC 6803"),
    "Picosynechococcus sp. PCC 7002":
        ("s7002", "PCC7002_locus_tag", "S5-2. Curated TFs PCC 7002"),
}

TAB     = os.path.join("outputs", "tables")
S5      = os.path.join("curation", "S5_curated_regulators.xlsx")
CURCORE = os.path.join("curation", "core_regulator_curation.csv")
CURFIG6 = os.path.join("curation", "figure6_categories.csv")


SHAPES_IN_KEY = {"oval", "sigma", "hexagon", "starburst", "diamond"}


def shape_of(product):
    """Regulator type, from the product name. Checked by eye before drawing.

    Order matters. "sigma factor" is tested before the starburst rule, since
    sigma factor product names also contain "RNA polymerase". Sigma-54
    activators (NtrC) are DNA-binding enhancer-binding proteins, not sigma
    factors, and fall through to oval.
    """
    p = str(product).lower()
    if "circadian clock protein kaic" in p:
        return "diamond"
    if "response regulator" in p:
        return "hexagon"
    if "sigma factor" in p:
        return "sigma"
    if ("antitermination" in p or "transcription termination" in p
            or "rna-binding" in p or "rna binding" in p):
        return "starburst"
    return "oval"


STRESS_CATEGORIES = {"multi-stress", "nutrient stress", "metal stress",
                     "uncharacterized TCS stress response",
                     "DNA-damage and replication stress"}


def core_category(stress_cat, non_stress_cat):
    """S6-2 classification -> figure category. Returns (category, note)."""
    s = "" if pd.isna(stress_cat) else str(stress_cat).strip().lower()
    n = "" if pd.isna(non_stress_cat) else str(non_stress_cat).strip().lower()
    if s:
        for key, cat in (("dna-damage", "DNA-damage and replication stress"),
                         ("metal", "metal stress"),
                         ("nutrient", "nutrient stress"),
                         ("two-component", "uncharacterized TCS stress response"),
                         ("multi", "multi-stress")):
            if key in s:
                return cat, ""
        return "", "unmapped S6-2 stress category: %s" % stress_cat
    if n:
        first = n.split(";")[0]
        for key, cat in (("circadian", "circadian clock"),
                         ("carbon", "carbon metabolism"),
                         ("housekeeping", "housekeeping and transcription"),
                         ("morphology", "cell morphology and motility")):
            if key in first:
                return cat, ""
        return "", "unmapped S6-2 category: %s" % non_stress_cat
    return "", "no S6-2 category"


def top_by(df, col, n):
    x = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
    # ties broken by IC, as in 13_crossvalidation.R
    ic = pd.to_numeric(df["IC"], errors="coerce").fillna(0.0)
    order = (pd.DataFrame({"v": x, "ic": ic})
               .sort_values(["v", "ic"], ascending=False, kind="stable").index)
    return df.loc[order[:n], "locus_tag"].tolist()


def main():
    for p in (S5, CURCORE, CURFIG6):
        if not os.path.exists(p):
            sys.exit("Input not found: %s" % p)

    core = pd.read_csv(os.path.join(TAB, "centralities_core.csv"))
    cen = {}
    for label, (stem, _, _) in SPECIES.items():
        f = os.path.join(TAB, "centralities_%s.csv" % stem)
        if not os.path.exists(f):
            sys.exit("Input not found: %s" % f)
        cen[label] = pd.read_csv(f)

    s54 = pd.read_excel(S5, sheet_name="S5-4. Conserved regulators")
    s6  = pd.read_csv(CURCORE)
    s6["Stress_coupled"] = (s6.Stress_category.fillna("").astype(str).str.strip()
                              .ne("").map({True: "Y", False: "N"}))
    core_name   = dict(zip(s6.locus_tag, s6.TF_name))
    core_stress = dict(zip(s6.locus_tag, s6.Stress_coupled))
    core_cat, core_note = {}, {}
    for _, r in s6.iterrows():
        core_cat[r.locus_tag], core_note[r.locus_tag] = core_category(
            r.get("Stress_category"), r.get("Non_stress_category"))

    fig5, fig6 = [], []
    for label, (stem, s5col, s5sheet) in SPECIES.items():
        d = cen[label]
        s5 = pd.read_excel(S5, sheet_name=s5sheet)
        sp_name = dict(zip(s5.locus_tag, s5.TF_name))
        sp_prod = dict(zip(s5.locus_tag, s5.product_name))
        # species locus -> core locus, for the tri-homologs
        to_core = dict(zip(s54[s5col], s54.PCC7942_locus_tag))

        for measure, (norm_col, group) in MEASURES.items():
            core_top = set(top_by(core, measure, CORE_K))
            sp_top   = top_by(d, measure, SPECIES_N)
            ic_max   = pd.to_numeric(d["IC"], errors="coerce").max()

            for loc in sp_top:
                row = d.loc[d.locus_tag == loc].iloc[0]
                score = float(row[norm_col])
                if measure == "IC":
                    score = score / ic_max
                core_loc = to_core.get(loc)
                crossval = core_loc in core_top if core_loc else False
                name = (core_name.get(core_loc) if crossval or core_loc
                        else None) or sp_name.get(loc)
                rec = dict(Species=label, Group=group, Measure=measure,
                           Species_locus=loc, Core_locus=core_loc,
                           Regulator=name if pd.notna(name) else "",
                           product_name=sp_prod.get(loc, ""),
                           Stress_coupled=core_stress.get(core_loc, ""),
                           score=round(score, 3),
                           ic_tiebreak=round(
                               float(row["IC"]) / ic_max if ic_max else 0.0, 4))
                (fig5 if crossval else fig6).append(rec)

    def rank_sizes(g, sym_max=SYM_MAX_FIG5):
        """Spread sizes by rank within the figure, top entry at the ceiling."""
        g = g.sort_values(["score", "ic_tiebreak"], ascending=False).copy()
        n = len(g)
        if n == 0:
            return g
        frac = 1.0 if n == 1 else 1.0 - np.arange(n) / (n - 1.0)
        g["font_size"] = [round((FONT_MIN + (FONT_MAX - FONT_MIN) * f)
                                / FONT_STEP) * FONT_STEP for f in frac]
        g["symbol_size"] = [round(SYM_MIN + (sym_max - SYM_MIN) * f, 2)
                            for f in frac]
        return g

    def collapse(recs, sym_max=SYM_MAX_FIG5):
        if not recs:
            return pd.DataFrame()
        df = pd.DataFrame(recs)
        g = (df.groupby(["Species", "Group", "Species_locus", "Core_locus",
                         "Regulator", "product_name", "Stress_coupled"],
                        dropna=False, as_index=False)
               .agg(measures=("Measure", lambda s: "+".join(sorted(set(s)))),
                    score=("score", "max"),
                    ic_tiebreak=("ic_tiebreak", "max")))
        g = rank_sizes(g, sym_max)
        g["Group"] = pd.Categorical(g.Group, GROUP_ORDER, ordered=True)
        return g.sort_values(["Species", "Group", "score"],
                             ascending=[True, True, False])

    f5 = collapse(fig5, SYM_MAX_FIG5)
    f6 = collapse(fig6, SYM_MAX_FIG6)
    f6["needs_name"] = np.where(f6.Regulator.fillna("").astype(str).str.strip() == "",
                                "Y", "")

    # curated names and categories for Figure 6
    warnings = []
    cur = pd.read_csv(CURFIG6)
    name_by = dict(zip(cur.Species_locus, cur.Regulator))
    cat_by  = dict(zip(cur.Species_locus, cur.category))
    for g in (f5, f6):
        g["Regulator"] = [name_by.get(l) if pd.notna(name_by.get(l, None))
                          and str(name_by.get(l)).strip() else r
                          for l, r in zip(g.Species_locus, g.Regulator)]
    f6["category"] = [cat_by.get(l, "") for l in f6.Species_locus]
    missing = sorted(set(f6.Species_locus) - set(cur.Species_locus))
    extra   = sorted(set(cur.Species_locus) - set(f6.Species_locus))
    if missing or extra:
        print("STOP: curation/figure6_categories.csv does not match Figure 6.")
        if missing: print("  not curated: " + ", ".join(missing))
        if extra:   print("  curated but not in Figure 6: " + ", ".join(extra))
        sys.exit(1)

    # Figure 5 takes its category from the core curation (Table 2)
    f5["category"] = [core_cat.get(c, "") for c in f5.Core_locus]
    for c in sorted(set(f5.Core_locus)):
        if core_note.get(c):
            warnings.append("Figure 5, %s (%s): %s"
                            % (core_name.get(c), c, core_note[c]))

    # Figure 6 tri-homologs must agree with the core curation
    seen = set()
    for _, r in f6.iterrows():
        c = r.Core_locus
        if not isinstance(c, str) or not c.strip() or r.Species_locus in seen:
            continue
        seen.add(r.Species_locus)
        ref = core_cat.get(c, "")
        got = "" if pd.isna(r.category) else str(r.category).strip()
        if ref and got != ref:
            warnings.append("Figure 6, %s (%s): curated category '%s' but "
                            "S6-2 gives '%s'" % (r.Regulator, r.Species_locus,
                                                 got, ref))

    # A core regulator must carry one classification everywhere. Disagreement
    # between the core curation and figure6_categories.csv, or between a
    # Figure 5 category and its stress flag, stops the run before any output
    # is written, so no figure or count is built on conflicting curation.
    fatal = [w for w in warnings if w.startswith("Figure 6, ") and "S6-2 gives" in w]
    st5_chk = f5.category.isin(STRESS_CATEGORIES)
    bad5 = f5.loc[st5_chk != (f5.Stress_coupled == "Y"), "Regulator"]
    if len(bad5):
        fatal.append("Figure 5: stress category and Stress_coupled disagree for %s"
                     % ", ".join(sorted(set(bad5))))
    if fatal:
        print("STOP: curation conflict, nothing written. Fix the curation files:")
        for w in fatal:
            print("  " + w)
        sys.exit(1)

    for g in (f5, f6):
        g["shape"] = g.product_name.map(shape_of)
        for loc, name, shp in zip(g.Species_locus, g.Regulator, g["shape"]):
            if shp not in SHAPES_IN_KEY:
                warnings.append("%s (%s): shape '%s' is not in the figure key"
                                % (name, loc, shp))
        g["label"] = [r.Regulator if str(r.Regulator).strip()
                      else re.sub(r"^(SYNPCC7942_|SYNPCC7002_|SGL_)", "",
                                  str(r.Species_locus))
                      for _, r in g.iterrows()]
    warnings = list(dict.fromkeys(warnings))

    f5.to_csv(os.path.join(TAB, "figure5_symbols.csv"), index=False)
    f6.to_csv(os.path.join(TAB, "figure6_symbols.csv"), index=False)

    # homology class of every regulator, from the ortholog groups
    og = pd.read_csv(os.path.join(TAB, "ortholog_groups.csv"))
    COLBY = {"S. elongatus PCC 7942": "PCC7942",
             "Synechocystis sp. PCC 6803": "PCC6803",
             "Picosynechococcus sp. PCC 7002": "PCC7002"}

    def homology(species, locus):
        col = COLBY[species]
        hit = og[og[col] == locus]
        if hit.empty:
            return "species-unique", ""
        r = hit.iloc[0]
        present = [c for c in ("PCC7942", "PCC7002", "PCC6803")
                   if pd.notna(r[c]) and str(r[c]).strip()]
        partners = [str(r[c]) for c in present if c != col]
        if len(present) == 3:
            return "tri-homolog", "; ".join(partners)
        if len(present) == 2:
            other = [c for c in present if c != col]
            return ("bi-homolog with %s" % other[0].replace("PCC", "PCC "),
                    "; ".join(partners))
        return "species-unique", ""

    # protein accession, from the per-species curated TF sheets
    prot = {}
    for _, (_, _, sheet) in SPECIES.items():
        s5s = pd.read_excel(S5, sheet_name=sheet)
        prot.update(dict(zip(s5s.locus_tag, s5s.protein_id)))

    # one row per Figure 6 regulator
    cat = (f6.drop_duplicates("Species_locus")
             [["Species", "Species_locus", "Regulator", "product_name",
               "shape"]].copy())
    cat["protein_id"] = cat.Species_locus.map(prot)
    hc = [homology(r.Species, r.Species_locus) for _, r in cat.iterrows()]
    cat["homology"] = [h[0] for h in hc]
    cat["homolog_partners"] = [h[1] for h in hc]
    cat["category"] = cat.Species_locus.map(cat_by)
    cat["appears_in_cells"] = [
        int((f6.Species_locus == loc).sum()) for loc in cat.Species_locus]
    cat = cat.sort_values(["Species", "appears_in_cells", "Species_locus"],
                          ascending=[True, False, True])
    cat.to_csv(os.path.join(TAB, "figure6_regulators.csv"), index=False)

    lines = []
    lines.append("Figure 5, cross-validated regulators")
    lines.append("  placements (species x measure group): %d" % len(f5))
    lines.append("  matches (regulator x species x measure): %d"
                 % sum(len(str(m).split("+")) for m in f5.measures))
    lines.append("  distinct core regulators: %d" % f5.Core_locus.nunique())
    lines.append("  stress-coupled placements: %d (%.0f%%)"
                 % ((f5.Stress_coupled == "Y").sum(),
                    100 * (f5.Stress_coupled == "Y").mean()))
    st5 = f5.category.isin(STRESS_CATEGORIES)
    if (st5 != (f5.Stress_coupled == "Y")).any():
        warnings.append("Figure 5: stress category and S6-2 Stress_coupled "
                        "disagree for %s" % ", ".join(sorted(set(
                            f5.loc[st5 != (f5.Stress_coupled == "Y"),
                                   "Regulator"]))))
    lines.append("")
    lines.append("Figure 6, species-specific regulators")
    lines.append("  placements: %d" % len(f6))
    lines.append("  distinct regulators: %d" % f6.Species_locus.nunique())
    unnamed = f6.Regulator.fillna("").astype(str).str.strip().isin(["", "nan"])
    lines.append("  placements without a regulator name: %d (%d distinct regulators)"
                 % (unnamed.sum(), f6.loc[unnamed, "Species_locus"].nunique()))
    if "category" in f6.columns and (f6.category.fillna("").astype(str).str.strip() != "").any():
        STRESS = STRESS_CATEGORIES
        st = f6.category.isin(STRESS)
        blank = f6.category.fillna("").astype(str).str.strip().isin(["", "nan"])
        lines.append("  stress-coupled placements: %d of %d (%.0f%%)"
                     % (st.sum(), len(f6), 100 * st.mean()))
        sd = f6.loc[st, "Species_locus"].nunique()
        lines.append("  stress-coupled distinct regulators: %d of %d"
                     % (sd, f6.Species_locus.nunique()))
        if blank.any():
            lines.append("  placements with no category: %d" % blank.sum())
        lines.append("    %-38s %8s %8s" % ("category", "distinct", "placements"))
        ent = f6.category.value_counts()
        for c, n in f6.drop_duplicates("Species_locus").category.value_counts().items():
            lines.append("    %-38s %8d %8d" % (c, n, ent.get(c, 0)))
    else:
        lines.append("  no curated categories found")
    if warnings:
        lines.append("")
        lines.append("WARNINGS (%d), resolve before drawing:" % len(warnings))
        lines.extend("  " + w for w in warnings)
    txt = "\n".join(lines)
    with open(os.path.join(TAB, "figure56_counts.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
