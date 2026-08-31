#!/usr/bin/env python3
"""
02_picrust_inference.py — Taxonomic Functional Prediction
===========================================================
Track 3 Host-Microbe Integration (Approach 2)

INPUT:  taxa_direction.tsv (from 01_clr_transform.py)

OUTPUT: taxa_functions.tsv with columns:
        taxon, direction, contaminant_flag,
        predicted_function, evidence_basis, literature_refs

METHOD:
  - Uses a CURATED LOOKUP TABLE of genus-level functional traits
    derived from published literature on gut microbiome functions.
  - This is NOT a PICRUSt2 run (which would require a reference tree
    and 16S sequences). Instead it's a literature-curated mapping
    since only genus-level ASV names are available, not sequences.
  - Each functional annotation is tagged with its evidence basis:
    "literature-curated" or "taxonomy-inferred"

CRITICAL RULES:
  - All function claims MUST be backed by the internal literature map
  - NO invented functions — if a genus isn't in the map, it gets "unknown"
  - Evidence basis is ALWAYS declared

USAGE:  python src/asv/02_picrust_inference.py \
            --in data/processed/clr_profiles/taxa_direction.tsv \
            --out data/processed/clr_profiles/taxa_functions.tsv
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


# ================================================================
# CURATED LITERATURE-BASED FUNCTIONAL MAP
# ================================================================
# Each entry maps a GENUS (lowercase) to known/predicted functions.
# Sources are noted inline. This is a STARTING POINT and should be
# expanded with systematic literature search.
# ================================================================

GENUS_FUNCTION_MAP = {
    # --- Short-chain fatty acid (SCFA) producers ---
    "acetobacter": {
        "functions": ["Acetate production", "Carbohydrate fermentation"],
        "pathways": ["Pyruvate metabolism", "Glycolysis / Gluconeogenesis"],
        "refs": ["Komagata et al. 2014, The Genus Acetobacter"],
        "category": "SCFA_producer",
    },
    "lactobacillus": {
        "functions": ["Lactate production", "Bacteriocin synthesis", "Vitamin B12 synthesis"],
        "pathways": ["Lactic acid fermentation", "Folate biosynthesis"],
        "refs": ["Wang et al. 2021, Lactobacillus in aquaculture"],
        "category": "SCFA_producer",
    },
    "bifidobacterium": {
        "functions": ["Acetate and lactate production", "Folate biosynthesis"],
        "pathways": ["Bifid shunt", "Folate biosynthesis"],
        "refs": ["O'Callaghan & van Sinderen 2016, Bifidobacteria"],
        "category": "SCFA_producer",
    },
    "clostridium": {
        "functions": ["Butyrate production", "Fiber degradation", "SCFA synthesis"],
        "pathways": ["Butyrate fermentation", "Cellulose degradation"],
        "refs": ["Louis & Flint 2017, Clostridium clusters in gut"],
        "category": "SCFA_producer",
    },
    "faecalibacterium": {
        "functions": ["Butyrate production", "Anti-inflammatory metabolite synthesis"],
        "pathways": ["Butyrate kinase pathway"],
        "refs": ["Sokol et al. 2008, Faecalibacterium prausnitzii"],
        "category": "SCFA_producer",
    },

    # --- Nitrogen metabolism ---
    "rhizobium": {
        "functions": ["Nitrogen fixation", "Ammonia assimilation"],
        "pathways": ["Nitrogen fixation (nif genes)", "Glutamine synthetase pathway"],
        "refs": ["Poole et al. 2018, Rhizobia"],
        "category": "nitrogen_metabolism",
    },
    "bradyrhizobium": {
        "functions": ["Nitrogen fixation", "Denitrification"],
        "pathways": ["Nitrogen fixation", "Denitrification pathway"],
        "refs": ["Schulte et al. 2021, Bradyrhizobium"],
        "category": "nitrogen_metabolism",
    },
    "azotobacter": {
        "functions": ["Nitrogen fixation (aerobic)", "PHA accumulation"],
        "pathways": ["Nitrogen fixation (nif)", "PHA biosynthesis"],
        "refs": ["Jiménez et al. 2019, Azotobacter as PGPR"],
        "category": "nitrogen_metabolism",
    },
    "endothiovibrio": {
        "functions": ["Sulfur oxidation", "Nitrogen fixation (putative)"],
        "pathways": ["Sulfur oxidation (sox)", "Nitrogen fixation (nif)"],
        "refs": ["Bazylinski et al. 2016, Endothiovibrio diazotrophicus"],
        "category": "sulfur_metabolism",
    },

    # --- Chitin degradation (KEY for crustacean gut) ---
    "vibrio": {
        "functions": ["Chitin degradation", "Chitinase secretion", "Protease secretion"],
        "pathways": ["Chitin degradation to GlcNAc", "Amino sugar metabolism"],
        "refs": ["Meibom et al. 2004, Vibrio chitin utilization"],
        "category": "chitin_degrader",
    },
    "aeromonas": {
        "functions": ["Chitin degradation", "Chitinase production", "Amylase secretion"],
        "pathways": ["Chitin degradation", "Starch and sucrose metabolism"],
        "refs": ["Janda & Abbott 2010, The Genus Aeromonas"],
        "category": "chitin_degrader",
    },
    "pseudoalteromonas": {
        "functions": ["Chitin degradation", "Extracellular enzyme secretion", "Bioactive compound synthesis"],
        "pathways": ["Chitin degradation", "Secondary metabolite biosynthesis"],
        "refs": ["Holmström & Kjelleberg 1999, Marine Pseudoalteromonas"],
        "category": "chitin_degrader",
    },

    # --- Protein/amino acid metabolism ---
    "pseudomonas": {
        "functions": ["Protein degradation", "Amino acid catabolism", "Biofilm formation"],
        "pathways": ["Branched-chain amino acid degradation", "TCA cycle"],
        "refs": ["Silby et al. 2011, Pseudomonas"],
        "category": "protein_fermenter",
    },
    "bacillus": {
        "functions": ["Protease secretion", "Amylase secretion", "Antimicrobial peptide synthesis"],
        "pathways": ["Extracellular protease", "Starch degradation"],
        "refs": ["Cutting 2011, Bacillus probiotics in aquaculture"],
        "category": "enzyme_producer",
    },
    "citrobacter": {
        "functions": ["Citrate utilization", "Mixed acid fermentation", "Hydrogen sulfide production"],
        "pathways": ["Citrate cycle", "Mixed acid fermentation"],
        "refs": ["Borenshtein & Schauer 2006, The Genus Citrobacter"],
        "category": "organic_acid_metabolism",
    },
    "klebsiella": {
        "functions": ["Mixed acid fermentation", "Nitrogen fixation", "Siderophore production"],
        "pathways": ["Mixed acid fermentation", "Nitrogen fixation (nif)"],
        "refs": ["Podschun & Ullmann 1998, Klebsiella spp."],
        "category": "organic_acid_metabolism",
    },

    # --- Vitamin / cofactor synthesis ---
    "flavobacterium": {
        "functions": ["Vitamin B12 synthesis (some spp.)", "Organic matter degradation"],
        "pathways": ["Cobalamin biosynthesis", "Extracellular enzyme secretion"],
        "refs": ["McBride 2014, The Family Flavobacteriaceae"],
        "category": "vitamin_synthesis",
    },

    # --- Sulfur metabolism ---
    "thiohalospira": {
        "functions": ["Sulfur oxidation", "Halophilic metabolism"],
        "pathways": ["Sulfur oxidation (sox)"],
        "refs": ["Sorokin et al. 2008, Thiohalospira"],
        "category": "sulfur_metabolism",
    },

    # --- Known contaminants (still noted but functionally annotated) ---
    "cutibacterium": {
        "functions": ["Lipid metabolism", "Propionate fermentation"],
        "pathways": ["Propionate fermentation", "Lipase activity"],
        "refs": ["Brüggemann et al. 2021, Cutibacterium acnes in skin"],
        "category": "lipid_metabolism",
        "contaminant_note": "Skin commensal; likely contaminant in gut samples",
    },
    "staphylococcus": {
        "functions": ["Facultative anaerobic metabolism", "Biofilm formation"],
        "pathways": ["Glycolysis", "Biofilm polysaccharide synthesis"],
        "refs": ["Otto 2018, Staphylococcal biofilms"],
        "category": "biofilm_former",
        "contaminant_note": "Skin/respiratory; possible contaminant",
    },
    "stenotrophomonas": {
        "functions": ["Multidrug resistance (intrinsic)", "Biofilm formation"],
        "pathways": ["Efflux pump systems", "Biofilm matrix"],
        "refs": ["Brooke 2012, Stenotrophomonas maltophilia"],
        "category": "biofilm_former",
        "contaminant_note": "Common reagent contaminant",
    },

    # --- Sediment/environmental bacteria ---
    "jiulongibacter": {
        "functions": ["Anaerobic metabolism", "Fermentation (putative)"],
        "pathways": ["Anaerobic fermentation"],
        "refs": ["Liu et al. 2016, Jiulongibacter sediminis gen. nov."],
        "category": "anaerobic_fermenter",
        "contaminant_note": "Sediment bacterium; environmental origin likely",
    },
}


def lookup_genus(genus: str) -> dict | None:
    """Look up functional annotations for a genus. Returns dict or None."""
    return GENUS_FUNCTION_MAP.get(genus.lower())


def annotate_taxon(taxon: str) -> dict:
    """
    Annotate a single taxon with functional predictions.
    Returns a dict with keys: predicted_function, evidence_basis,
    pathways, literature_refs, category, contaminant_note.
    """
    genus = taxon.split()[0] if taxon else ""
    entry = lookup_genus(genus)

    if entry:
        return {
            "predicted_function": "; ".join(entry.get("functions", [])),
            "evidence_basis": "literature-curated",
            "pathways": "; ".join(entry.get("pathways", [])),
            "literature_refs": "; ".join(entry.get("refs", [])),
            "category": entry.get("category", "unknown"),
            "contaminant_note": entry.get("contaminant_note", ""),
        }
    else:
        # No literature entry — mark as unknown, not invented
        return {
            "predicted_function": "unknown (no curated entry)",
            "evidence_basis": "none",
            "pathways": "",
            "literature_refs": "",
            "category": "unknown",
            "contaminant_note": "",
        }


def main():
    parser = argparse.ArgumentParser(description="Literature-based functional annotation of ASV taxa")
    parser.add_argument("--in", dest="input_file", required=True, help="CLR direction TSV")
    parser.add_argument("--out", dest="output_file", required=True, help="Output annotated TSV")
    args = parser.parse_args()

    # ---- Load ----
    print(f"[FUNC] Loading direction table: {args.input_file}")
    df = pd.read_csv(args.input_file, sep="\t")
    print(f"[FUNC] {len(df)} taxa loaded")

    # ---- Annotate each taxon ----
    print("[FUNC] Annotating taxa from literature map...")
    annotations = []
    for taxon in df["taxon"]:
        ann = annotate_taxon(taxon)
        annotations.append(ann)

    ann_df = pd.DataFrame(annotations)
    df_out = pd.concat([df, ann_df], axis=1)

    # ---- Quality check ----
    n_known = (df_out["evidence_basis"] == "literature-curated").sum()
    n_unknown = (df_out["evidence_basis"] == "none").sum()
    print(f"[FUNC] Annotated: {n_known} from literature, {n_unknown} unknown")

    # ---- Select output columns ----
    out_cols = [
        "taxon", "Jumper", "Laggard",
        "rel_abund_jumper_pct", "rel_abund_laggard_pct",
        "fold_diff", "direction", "exclusive_to",
        "predicted_function", "pathways", "category",
        "evidence_basis", "literature_refs",
        "contaminant_flag", "contaminant_reason", "contaminant_note",
    ]
    # Keep only columns that exist
    out_cols = [c for c in out_cols if c in df_out.columns]
    df_out = df_out[out_cols]

    # ---- Write ----
    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(args.output_file, sep="\t", index=False, na_rep="")
    print(f"[FUNC] Output written: {args.output_file}")
    print(f"[FUNC] Literature map covers {len(GENUS_FUNCTION_MAP)} genera")
    print("[FUNC] Done.")


if __name__ == "__main__":
    main()
