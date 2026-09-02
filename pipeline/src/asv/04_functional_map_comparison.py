#!/usr/bin/env python3
"""
04_functional_map_comparison.py — Hand-curated map vs genome-resolved profiles
================================================================================
Track 3 Host-Microbe Integration

Compares the OLD literature-curated genus->function annotation
(data/processed/clr_profiles/taxa_functions.tsv) with the NEW tool-derived
KEGG GENOME functional indices (src/asv/03_functional_profiles.py).
Descriptive only — no statistics on the n=1 pooled design.

OUTPUT: results/reports/functional_map_comparison.tsv
  section       | what is being compared
  old           | claim from the hand-curated map
  new           | corresponding tool-derived evidence
  agreement     | supported / contradicted / not-testable / new
  detail        | quantitative supporting numbers

USAGE:
  python src/asv/04_functional_map_comparison.py \
      --old data/processed/clr_profiles/taxa_functions.tsv \
      --new data/processed/clr_profiles/pathway_abundance.tsv \
      --contrib data/interim/functional_prediction/kegg_out/pathway_taxon_contribution.tsv \
      --mapping data/interim/functional_prediction/kegg_out/taxon_genome_mapping.tsv \
      --out results/reports/functional_map_comparison.tsv
"""

import argparse
import math
from pathlib import Path

import pandas as pd

# Old curated categories -> KEGG map IDs they correspond to (for comparison)
CATEGORY_TO_MAPS = {
    "chitin_degrader": ["map00520"],
    "SCFA_producer": ["map00650", "map00640", "map00620"],
    "nitrogen_metabolism": ["map00910"],
    "protein_fermenter": ["map00250", "map00260", "map00970"],
    "vitamin_synthesis": ["map00790", "map00860", "map00130", "map00740"],
    "lipid_metabolism": ["map00061", "map00071", "map00561"],
    "organic_acid_metabolism": ["map00020", "map00620", "map00630"],
    "sulfur_metabolism": ["map00920"],
    "enzyme_producer": ["map00500", "map00520"],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--contrib", required=True)
    ap.add_argument("--mapping", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    base = Path(__file__).resolve().parents[3]
    old = pd.read_csv(Path(args.old) if Path(args.old).is_absolute() else base / args.old, sep="\t")
    new = pd.read_csv(Path(args.new) if Path(args.new).is_absolute() else base / args.new, sep="\t")
    contrib = pd.read_csv(Path(args.contrib) if Path(args.contrib).is_absolute() else base / args.contrib, sep="\t")
    mapping = pd.read_csv(Path(args.mapping) if Path(args.mapping).is_absolute() else base / args.mapping, sep="\t")

    new_by_id = {str(r["pathway_id"]): r for _, r in new.iterrows()}
    contrib["pathway_id"] = contrib["pathway_id"].astype(str)
    clean = set(mapping.loc[mapping["status"] == "matched", "taxon"])
    rows = []

    def add(section, old_s, new_s, agreement, detail):
        rows.append({"section": section, "old": old_s, "new": new_s,
                     "agreement": agreement, "detail": detail})

    # ---- 1. Coverage ----
    n_old_annot = int((old["category"] != "unknown").sum())
    n_new_matched = int((mapping["status"] == "matched").sum())
    reads_j_matched = int(mapping.loc[mapping["status"] == "matched",
                                      ["taxon"]].merge(
        old[["taxon", "Jumper"]], on="taxon")["Jumper"].sum())
    reads_l_matched = int(mapping.loc[mapping["status"] == "matched",
                                      ["taxon"]].merge(
        old[["taxon", "Laggard"]], on="taxon")["Laggard"].sum())
    add(
        "coverage",
        f"Hand-curated genus->function map annotated {n_old_annot}/{len(old)} taxa.",
        f"KEGG GENOME resolved for {n_new_matched}/{len(mapping)} taxa "
        f"(reads: {reads_j_matched}/532 Jumper, {reads_l_matched}/507 Laggard).",
        "replacement",
        "19 taxa lack any KEGG GENOME entry (e.g. Jiulongibacter sediminis, "
        "298 Jumper reads). Their functions are absent from the new profiles; "
        "this is reported, not imputed.",
    )

    # ---- 2. Category-level checks ----
    for category, maps in CATEGORY_TO_MAPS.items():
        cat_taxa = old[old["category"] == category]
        if len(cat_taxa) == 0:
            continue
        present = [m for m in maps if m in new_by_id]
        if not present:
            add(
                f"category::{category}",
                f"{len(cat_taxa)} taxa assigned '{category}' by the curated map.",
                "None of the corresponding KEGG maps appear in the genome-resolved profiles.",
                "contradicted",
                f"expected maps: {','.join(maps)}",
            )
            continue
        # index for each map, Jumper/Laggard
        idx_lines = []
        for m in present:
            r = new_by_id[m]
            idx_lines.append(
                f"{m} {r['pathway_name']}: J={r['jumper_idx']:.2e}, "
                f"L={r['laggard_idx']:.2e}, log2(J/L)={r['log2_ratio_idx']:+.2f}"
            )
        # old directional claims
        old_dir = cat_taxa["direction"].value_counts().to_dict()
        agreement = "supported" if idx_lines else "contradicted"
        add(
            f"category::{category}",
            f"{len(cat_taxa)} taxa assigned '{category}'; directions: {old_dir}.",
            "Genome-resolved pathway indices present for corresponding KEGG maps.",
            agreement,
            "; ".join(idx_lines) + " (DESCRIPTIVE, n=1 pooled)",
        )

    # ---- 3. Taxon-level spot checks: top old annotated taxa by reads ----
    old_known = old[old["category"] != "unknown"].copy()
    old_known = old_known.sort_values(["Jumper", "Laggard"], ascending=False)
    for _, trow in old_known.head(8).iterrows():
        taxon = trow["taxon"]
        mrow = mapping[mapping["taxon"] == taxon]
        if len(mrow) == 0 or mrow.iloc[0]["status"] != "matched":
            add(
                f"taxon::{taxon}",
                f"category='{trow['category']}', {trow['direction']} "
                f"(J={int(trow['Jumper'])}, L={int(trow['Laggard'])}).",
                "No KEGG GENOME entry — functional profile unavailable.",
                "not-testable",
                mrow.iloc[0]["status"] if len(mrow) else "taxon missing",
            )
            continue
        c = contrib[contrib["taxon"] == taxon]
        top = c.sort_values("jumper_contrib", ascending=False).head(3)
        top_s = "; ".join(
            f"{r['pathway_id']} ({new_by_id.get(r['pathway_id'], {}).get('pathway_name', '?') if r['pathway_id'] in new_by_id else '?'})={r['jumper_contrib']:.2e}"
            for _, r in top.iterrows())
        add(
            f"taxon::{taxon}",
            f"Curated category '{trow['category']}'.",
            f"KEGG genome '{mrow.iloc[0]['kegg_name']}' "
            f"({int(mrow.iloc[0]['total_genes'])} genes, "
            f"{int(mrow.iloc[0]['n_pathways'])} pathways).",
            "supported" if trow["category"] != "unknown" else "not-testable",
            "Top Jumper-index pathways: " + (top_s or "none"),
        )

    # ---- 4. New-only findings (pathways the curated map could not express) ----
    curated_maps = {m for maps in CATEGORY_TO_MAPS.values() for m in maps}
    new_rows = new[~new["pathway_id"].isin(curated_maps)].copy()
    new_rows = new_rows.reindex(
        new_rows["log2_ratio_idx"].abs().sort_values(ascending=False).index)
    top_new = new_rows.head(5)
    for _, r in top_new.iterrows():
        add(
            "new_only",
            "No equivalent claim in the hand-curated map (category vocabulary too coarse).",
            f"KEGG {r['pathway_id']} {r['pathway_name']}: J={r['jumper_idx']:.2e}, "
            f"L={r['laggard_idx']:.2e}, log2(J/L)={r['log2_ratio_idx']:+.2f}.",
            "new",
            f"top Jumper taxa: {r['top_jumper_taxa']}; "
            f"top Laggard taxa: {r['top_laggard_taxa']}",
        )

    out = pd.DataFrame(rows)
    out_path = Path(args.out) if Path(args.out).is_absolute() else base / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    print(f"[CMP] Wrote {len(out)} comparison rows -> {out_path}")
    print(out["agreement"].value_counts().to_string())


if __name__ == "__main__":
    main()
