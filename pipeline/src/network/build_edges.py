#!/usr/bin/env python3
"""
build_edges.py — Build Integration Network Edges
==================================================
Track 3 Host-Microbe Integration

INPUTS:
  --taxa:  taxa_functions.tsv (from 02_picrust_inference.py)
  --genes: hub_genes.tsv (from WGCNA — may be empty if host step not yet run)

OUTPUT:  network_edges.csv with columns:
  edge_id, source_type, source_name, target_type, target_name,
  edge_basis, evidence_summary, direction_concordance, weight_meaning

CRITICAL RULES:
  - NO cross-layer covariance or correlation coefficients
  - Edges are typed: "phenotype_concordance", "predicted_function_overlap",
    or "curated_literature_interaction"
  - Every edge declares exactly what its weight means
  - Evidence must be traceable to source data OR literature

USAGE:  python src/network/build_edges.py \
            --taxa data/processed/clr_profiles/taxa_functions.tsv \
            --genes data/processed/wgcna/hub_genes.tsv \
            --out results/shortlist/network_edges.csv
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


# ================================================================
# MICROBIAL FUNCTION ↔ HOST PATHWAY BRIDGE
# ================================================================
# The microbial pathway assignments come from KEGG GENOME (tool-derived,
# src/asv/03_functional_profiles.py) — NOT from a hand-curated genus
# dictionary. The bridge from a microbial KEGG pathway to a host KEGG
# pathway is literature-curated (biology), anchored to concrete KEGG map
# IDs in data/interim/literature_tables/kegg_pathway_bridge.tsv.
# ================================================================

def load_bridge(bridge_path: Path) -> pd.DataFrame:
    """Load the curated microbial->host KEGG pathway bridge table."""
    return pd.read_csv(bridge_path, sep="\t")


def build_function_edges_from_profiles(pathway_df: pd.DataFrame,
                                       contrib_df: pd.DataFrame,
                                       taxa_df: pd.DataFrame,
                                       bridge_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build function-overlap edges from genome-resolved pathway indices.

    For each bridge row whose microbial KEGG pathway is present in the
    profiles:
      - source = microbial KEGG pathway (map ID), with pool indices
      - target = host KEGG pathway (name + map ID)
      - edge_confidence = share of the enriched pool's pathway index
        contributed by LOW-contamination-risk taxa (honest weighting,
        no stats on the n=1 pools).
    """
    edges = []
    edge_id = 0

    # clean-taxon lookup
    taxa_df = taxa_df.copy()
    taxa_df["clean"] = (~taxa_df.get("contaminant_flag", pd.Series([False] * len(taxa_df))).fillna(False)) & \
                       (taxa_df.get("contamination_score", pd.Series([0.0] * len(taxa_df))).fillna(0.0) < 0.85)
    clean_map = dict(zip(taxa_df["taxon"], taxa_df["clean"]))

    by_id = {str(r["pathway_id"]): r for _, r in pathway_df.iterrows()}
    contrib_by_id = {}
    for pid, grp in contrib_df.groupby("pathway_id"):
        contrib_by_id[str(pid)] = grp

    for _, bridge in bridge_df.iterrows():
        mid = str(bridge["microbial_map"])
        row = by_id.get(mid)
        if row is None:
            continue
        j_idx = float(row["jumper_idx"])
        l_idx = float(row["laggard_idx"])
        if j_idx <= 0 and l_idx <= 0:
            continue
        ratio = float(row["log2_ratio_idx"])
        if ratio > 0:
            pool = "Jumper"
            concord = "Jumper-enriched"
        elif ratio < 0:
            pool = "Laggard"
            concord = "Laggard-enriched"
        else:
            pool = "Jumper"
            concord = "balanced"

        # contamination-aware share for the enriched pool
        contrib = contrib_by_id.get(mid)
        clean_share = 1.0
        top_taxa = ""
        if contrib is not None:
            col = "jumper_contrib" if pool == "Jumper" else "laggard_contrib"
            c = contrib[contrib[col] > 0].copy()
            total = c[col].sum()
            if total > 0:
                c["clean"] = c["taxon"].map(clean_map).fillna(True)
                clean_share = float((c.loc[c["clean"], col].sum()) / total)
                top = c.sort_values(col, ascending=False).head(8)
                top_taxa = "; ".join(f"{t} ({v:.2e})"
                                     for t, v in zip(top["taxon"], top[col]))

        risk = "LOW" if clean_share >= 0.8 else ("MEDIUM" if clean_share >= 0.5 else "HIGH")
        edge_id += 1
        edges.append({
            "edge_id": f"E{edge_id:04d}",
            "source_type": "microbial_pathway",
            "source_name": f"{mid} ({bridge['microbial_name']})",
            "target_type": "host_pathway",
            "target_name": bridge["host_pathway"],
            "host_map": bridge["host_map"],
            "edge_basis": "predicted_function_overlap",
            "evidence_summary": (
                f"Tool-derived KEGG GENOME functional index (not a literature "
                f"taxon map): microbial pathway {mid} '{bridge['microbial_name']}' "
                f"has pool indices Jumper={j_idx:.5f}, Laggard={l_idx:.5f} "
                f"(log2 ratio {ratio:+.2f}, DESCRIPTIVE ONLY, n=1 pooled per "
                f"group). {bridge['host_effect']}"
            ),
            "direction_concordance": concord,
            "is_contaminant_source": False,
            "contamination_score": round(1.0 - clean_share, 3),
            "contamination_risk": risk,
            "edge_confidence": round(clean_share, 3),
            "weight_meaning": (
                f"Confidence = share of the {pool}-pool pathway index "
                f"contributed by LOW-contamination-risk taxa "
                f"({clean_share:.2f}). Pathway assignment from KEGG GENOME; "
                f"host-pathway bridge literature-curated. "
                f"NO cross-layer correlation claimed."
            ),
            "literature_refs": bridge["refs"],
            "source_taxa": top_taxa,
        })

    return pd.DataFrame(edges)


def build_phenotype_concordance_edges(taxa_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build edges based on phenotype concordance: taxa enriched in Jumpers
    are concordant with the fast-growth phenotype.
    
    v2: Edges include contamination_score and edge_confidence.
    """
    edges = []
    edge_id = 10000  # start at a different range

    jumper_taxa = taxa_df[taxa_df["direction"].str.contains("Jumper", na=False)]

    for _, row in jumper_taxa.iterrows():
        taxon = row["taxon"]
        contam_score = row.get("contamination_score", 0.0)
        contam_risk = row.get("contamination_risk", "LOW")
        edge_id += 1
        edge_weight = round(1.0 - float(contam_score), 3)
        edges.append({
            "edge_id": f"E{edge_id:04d}",
            "source_type": "microbial_taxon",
            "source_name": taxon,
            "target_type": "host_phenotype",
            "target_name": "fast_growth_Jumper",
            "edge_basis": "phenotype_concordance",
            "evidence_summary": f"Taxon enriched in Jumper (fast-growth) group. "
                                f"CLR fold-diff: {row['fold_diff']:.2f}. "
                                f"No statistical inference possible (n=1 pooled).",
            "direction_concordance": row["direction"],
            "is_contaminant_source": row.get("contaminant_flag", False),
            "contamination_score": contam_score,
            "contamination_risk": contam_risk,
            "edge_confidence": edge_weight,
            "weight_meaning": f"Confidence weight = 1 − contamination_score "
                              f"({edge_weight:.3f}). Directional association only. "
                              f"NOT a correlation.",
            "literature_refs": "",
        })

    return pd.DataFrame(edges)


def main():
    parser = argparse.ArgumentParser(description="Build integration network edges")
    parser.add_argument("--taxa", required=True, help="Taxa functions TSV")
    parser.add_argument("--genes", required=True, help="Hub genes TSV (from WGCNA)")
    parser.add_argument("--out", required=True, help="Output edge table CSV")
    parser.add_argument("--pathways",
                        default="data/processed/clr_profiles/pathway_abundance.tsv",
                        help="Genome-resolved pathway index table (03_functional_profiles.py)")
    parser.add_argument("--contrib",
                        default="data/interim/functional_prediction/kegg_out/pathway_taxon_contribution.tsv",
                        help="Pathway x taxon contribution long table")
    parser.add_argument("--bridge",
                        default="data/interim/literature_tables/kegg_pathway_bridge.tsv",
                        help="Curated microbial->host KEGG pathway bridge")
    parser.add_argument("--exclude-contaminants", action="store_true",
                        help="Exclude contaminant-flagged taxa from phenotype-edge building")
    args = parser.parse_args()

    # ---- Load taxa ----
    print(f"[EDGES] Loading taxa functions: {args.taxa}")
    taxa_df = pd.read_csv(args.taxa, sep="\t")
    print(f"[EDGES] {len(taxa_df)} taxa loaded")

    # ---- Filter contaminants if requested (v2: supports score-based filtering) ----
    n_contam_flag = taxa_df["contaminant_flag"].sum() if "contaminant_flag" in taxa_df.columns else 0
    n_high_risk = (taxa_df.get("contamination_risk", pd.Series(["LOW"]*len(taxa_df))) == "HIGH").sum()
    if args.exclude_contaminants and n_contam_flag > 0:
        # Exclude taxonomy-flagged AND high-risk scored taxa
        mask = ~taxa_df["contaminant_flag"]
        if "contamination_score" in taxa_df.columns:
            mask = mask & (taxa_df["contamination_score"] < 0.85)
        clean = taxa_df[mask].copy()
        n_excluded = len(taxa_df) - len(clean)
        print(f"[EDGES] Excluding {n_excluded} contaminant taxa "
              f"(flagged={n_contam_flag}, high-risk-scored={n_high_risk}) "
              f"→ {len(clean)} clean taxa for edge building")
        edge_taxa = clean
    else:
        if n_contam_flag > 0 or n_high_risk > 0:
            print(f"[EDGES] INFO: {n_contam_flag} contaminant-flagged, "
                  f"{n_high_risk} high-risk-scored taxa present. "
                  f"Use --exclude-contaminants to filter. "
                  f"Edges from high-score taxa carry low edge_confidence.")
        edge_taxa = taxa_df

    # ---- Build function-overlap edges from genome-resolved profiles ----
    print("[EDGES] Building predicted-function-overlap edges (KEGG GENOME profiles)...")
    pathway_path = Path(args.pathways)
    contrib_path = Path(args.contrib)
    bridge_path = Path(args.bridge)
    if pathway_path.exists() and bridge_path.exists():
        pathway_df = pd.read_csv(pathway_path, sep="\t")
        bridge_df = load_bridge(bridge_path)
        contrib_df = (pd.read_csv(contrib_path, sep="\t")
                      if contrib_path.exists()
                      else pd.DataFrame(columns=["pathway_id", "taxon",
                                                 "jumper_contrib", "laggard_contrib"]))
        func_edges = build_function_edges_from_profiles(
            pathway_df, contrib_df, taxa_df, bridge_df)
    else:
        print(f"[EDGES] WARNING: pathway profile ({pathway_path}) or bridge "
              f"({bridge_path}) missing. Run src/asv/03_functional_profiles.py "
              f"first. Building phenotype edges only.")
        func_edges = pd.DataFrame()
    print(f"[EDGES] {len(func_edges)} function-overlap edges")

    # ---- Build phenotype concordance edges ----
    print("[EDGES] Building phenotype-concordance edges...")
    pheno_edges = build_phenotype_concordance_edges(edge_taxa)
    print(f"[EDGES] {len(pheno_edges)} phenotype-concordance edges")

    # ---- Combine ----
    all_edges = pd.concat([func_edges, pheno_edges], ignore_index=True)
    # normalise dtypes for concat with empty frames
    for col in ["contamination_score", "contamination_risk", "edge_confidence"]:
        if col in all_edges.columns:
            all_edges[col] = all_edges[col].astype(object)

    if len(all_edges) == 0:
        print("[EDGES] WARNING: No edges built. Check functional annotations.")

    # ---- Try loading hub genes (may not exist yet) ----
    genes_path = Path(args.genes)
    if genes_path.exists():
        print(f"[EDGES] Loading hub genes: {args.genes}")
        hub_df = pd.read_csv(args.genes, sep="\t")
        print(f"[EDGES] {len(hub_df)} hub genes — host-layer edges deferred to WGCNA output")
        # Note: host-layer edges from WGCNA module-trait correlations
        # are built separately. This script bridges microbial→host.
    else:
        print(f"[EDGES] Hub genes file not found: {args.genes}")
        print("[EDGES] Host-layer edges will be added after WGCNA completes.")
        print("[EDGES] This is expected if host RNA-seq has not yet been processed.")

    # ---- Write output ----
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    all_edges.to_csv(args.out, index=False)
    print(f"[EDGES] Output: {args.out} ({len(all_edges)} edges)")
    print("[EDGES] Edge type counts:")
    print(all_edges["edge_basis"].value_counts().to_string())
    print("[EDGES] Done.")


if __name__ == "__main__":
    main()
