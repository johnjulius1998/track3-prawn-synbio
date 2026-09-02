#!/usr/bin/env python3
"""
m6_gene_annotation.py — M6 Module Gene Functional Annotation & GSEA
=====================================================================
Track 3 Host-Microbe Integration (v3.2)

PURPOSE:
  Annotate M6 hub genes (72 genes, the strongest growth-associated module)
  and test whether they are enriched for specific biological functions
  that would explain the negative correlation with weight gain.

  M6 partial_r(WG|sex) = −0.590: higher M6 expression → lower WG.
  Hypotheses: protein degradation, immune response, stress response,
  or metabolic cost pathways.

METHOD:
  1. Query NCBI eutils for XM_ accession protein products
  2. Extract GO terms and protein descriptions
  3. Manual curation of functional categories
  4. If sufficient annotations found, run clusterProfiler GSEA in R

OUTPUTS:
  results/reports/
    m6_gene_annotations.tsv     — Per-gene NCBI annotations
    m6_functional_summary.tsv   — Functional category counts

USAGE:
  python src/network/m6_gene_annotation.py \
      --hub-genes data/processed/wgcna/hub_genes.tsv \
      --modules results/tables/wgcna_modules.csv \
      --out-dir results/reports/
"""

import argparse
import sys
import json
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict
from urllib.request import urlopen, Request
from urllib.error import HTTPError

import numpy as np
import pandas as pd


# ================================================================
# NCBI E-utilities query
# ================================================================

def fetch_ncbi_gene_info(accession, retries=3, delay=0.5):
    """
    Fetch gene/protein info for an XM_ accession from NCBI.
    Returns dict with title, organism, and any available annotation.
    """
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=nuccore&id={accession}&retmode=json"
    )
    
    for attempt in range(retries):
        try:
            req = Request(url, headers={"User-Agent": "Track3-HostMicrobe/1.0"})
            with urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read())
            
            result = data.get("result", {})
            uids = result.get("uids", [])
            if not uids:
                return {"accession": accession, "status": "not_found"}
            
            uid = uids[0]
            record = result.get(uid, {})
            
            return {
                "accession": accession,
                "status": "found",
                "title": record.get("title", ""),
                "organism": record.get("organism", ""),
                "taxid": record.get("taxid", ""),
                "length": record.get("slen", 0),
                "gb_division": record.get("gbdiv", ""),
                "update_date": record.get("update date", ""),
            }
        except HTTPError as e:
            if e.code == 429:
                time.sleep(2 * (attempt + 1))
            elif e.code == 404:
                return {"accession": accession, "status": "not_found"}
            else:
                time.sleep(delay * (attempt + 1))
        except Exception as e:
            time.sleep(delay * (attempt + 1))
    
    return {"accession": accession, "status": "error"}


def fetch_protein_for_mrna(accession, retries=3, delay=0.5):
    """
    Try to find the protein product (XP_ accession) linked to an XM_ mRNA.
    Uses NCBI elink from nuccore to protein database.
    """
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
        f"?dbfrom=nuccore&db=protein&id={accession}&retmode=json"
    )
    
    for attempt in range(retries):
        try:
            req = Request(url, headers={"User-Agent": "Track3-HostMicrobe/1.0"})
            with urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read())
            
            linksets = data.get("linksets", [])
            if not linksets:
                return None
            
            links = linksets[0].get("linksetdbs", [])
            for link_db in links:
                if link_db.get("linkname") == "nuccore_protein":
                    protein_ids = link_db.get("links", [])
                    if protein_ids:
                        return protein_ids[0]  # Return first protein ID
            return None
        except HTTPError:
            time.sleep(delay * (attempt + 1))
        except Exception:
            time.sleep(delay * (attempt + 1))
    
    return None


def fetch_protein_title(protein_uid, retries=3, delay=0.5):
    """Fetch protein title from NCBI."""
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=protein&id={protein_uid}&retmode=json"
    )
    
    for attempt in range(retries):
        try:
            req = Request(url, headers={"User-Agent": "Track3-HostMicrobe/1.0"})
            with urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read())
            
            result = data.get("result", {})
            uids = result.get("uids", [])
            if uids:
                record = result.get(uids[0], {})
                return record.get("title", "")
            return ""
        except Exception:
            time.sleep(delay * (attempt + 1))
    
    return ""


# ================================================================
# Functional keyword classification
# ================================================================

# Curated keyword → functional category mapping for crustacean biology
FUNCTIONAL_KEYWORDS = {
    # Protein degradation / turnover
    "proteasome": "protein_degradation",
    "ubiquitin": "protein_degradation",
    "protease": "protein_degradation",
    "peptidase": "protein_degradation",
    "cathepsin": "protein_degradation",
    "lysosomal": "protein_degradation",
    "autophagy": "protein_degradation",
    "e3 ubiquitin": "protein_degradation",
    "f-box": "protein_degradation",
    
    # Immune response
    "immune": "immune_response",
    "immunoglobulin": "immune_response",
    "lectin": "immune_response",
    "antimicrobial": "immune_response",
    "defensin": "immune_response",
    "toll": "immune_response",
    "lysozyme": "immune_response",
    "complement": "immune_response",
    "prophenoloxidase": "immune_response",
    "hemocyanin": "immune_response",
    
    # Stress response
    "heat shock": "stress_response",
    "hsp": "stress_response",
    "chaperone": "stress_response",
    "stress": "stress_response",
    "antioxidant": "stress_response",
    "glutathione": "stress_response",
    "peroxidase": "stress_response",
    "catalase": "stress_response",
    "detoxif": "stress_response",
    "cytochrome p450": "stress_response",
    
    # Energy metabolism
    "atp": "energy_metabolism",
    "mitochondri": "energy_metabolism",
    "oxidase": "energy_metabolism",
    "dehydrogenase": "energy_metabolism",
    "cytochrome c": "energy_metabolism",
    "nadh": "energy_metabolism",
    "citrate": "energy_metabolism",
    
    # Growth / molting
    "ecdys": "growth_molting",
    "molt": "growth_molting",
    "chitin": "growth_molting",
    "cuticle": "growth_molting",
    "insulin": "growth_molting",
    "growth factor": "growth_molting",
    
    # Transcription / signaling
    "transcription factor": "transcription_signaling",
    "kinase": "transcription_signaling",
    "phosphatase": "transcription_signaling",
    "gtpase": "transcription_signaling",
    "receptor": "transcription_signaling",
    "signaling": "transcription_signaling",
    
    # Translation / ribosome
    "ribosom": "translation_ribosome",
    "translation": "translation_ribosome",
    "trna": "translation_ribosome",
    "elongation factor": "translation_ribosome",
}


def classify_function(title):
    """Classify a protein/gene title into functional categories."""
    if not title:
        return ["unknown"]
    
    title_lower = title.lower()
    categories = set()
    
    for keyword, category in FUNCTIONAL_KEYWORDS.items():
        if keyword in title_lower:
            categories.add(category)
    
    if not categories:
        categories.add("other_uncharacterized")
    
    return sorted(categories)


# ================================================================
# Main
# ================================================================

def main():
    parser = argparse.ArgumentParser(description="M6 gene annotation and GSEA")
    parser.add_argument("--hub-genes", required=True)
    parser.add_argument("--modules", required=True)
    parser.add_argument("--out-dir", default="results/reports")
    parser.add_argument("--max-queries", type=int, default=72,
                        help="Max NCBI queries (rate-limited)")
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("M6 MODULE GENE FUNCTIONAL ANNOTATION (v3.2)")
    print("=" * 70)
    
    # Load M6 hub genes
    hub = pd.read_csv(args.hub_genes, sep="\t")
    modules = pd.read_csv(args.modules)
    
    m6_genes = hub[hub["module"] == "M6"]["gene"].values
    print(f"\n[1/3] M6 module: {len(m6_genes)} hub genes")
    print(f"  Partial r(WG|sex) = −0.590 (higher expression → lower WG)")
    
    # Query NCBI (with rate limiting)
    print(f"\n[2/3] Querying NCBI for {min(len(m6_genes), args.max_queries)} genes...")
    print(f"  (Rate-limited: ~0.5-1 sec per query)")
    
    annotations = []
    protein_hits = 0
    
    for i, gene in enumerate(m6_genes[:args.max_queries]):
        sys.stdout.write(f"\r  {i+1}/{min(len(m6_genes), args.max_queries)}: {gene[:30]}...")
        sys.stdout.flush()
        
        # Fetch gene info
        info = fetch_ncbi_gene_info(gene)
        
        # Try to get protein product
        protein_uid = None
        protein_title = ""
        if info.get("status") == "found":
            protein_uid = fetch_protein_for_mrna(gene)
            if protein_uid:
                protein_title = fetch_protein_title(protein_uid)
                protein_hits += 1
        
        # Classify function
        categories = classify_function(info.get("title", ""))
        if protein_title:
            prot_categories = classify_function(protein_title)
            categories = sorted(set(categories + prot_categories))
        
        annotations.append({
            "gene": gene,
            "ncbi_status": info.get("status", "error"),
            "ncbi_title": info.get("title", "")[:200],
            "organism": info.get("organism", ""),
            "seq_length": info.get("length", 0),
            "protein_accession": f"protein:{protein_uid}" if protein_uid else "",
            "protein_title": protein_title[:200],
            "functional_categories": "; ".join(categories),
            "primary_category": categories[0] if categories else "unknown",
        })
        
        time.sleep(0.3)  # Rate limit for NCBI
    
    print(f"\n  Queried {len(annotations)} genes")
    print(f"  Found on NCBI: {sum(1 for a in annotations if a['ncbi_status'] == 'found')}")
    print(f"  With protein products: {protein_hits}")
    
    # Summarize functional categories
    print(f"\n[3/3] Functional category summary:")
    ann_df = pd.DataFrame(annotations)
    
    cat_counts = defaultdict(int)
    cat_genes = defaultdict(list)
    for _, row in ann_df.iterrows():
        for cat in row["functional_categories"].split("; "):
            cat = cat.strip()
            if cat:
                cat_counts[cat] += 1
                cat_genes[cat].append(row["gene"])
    
    print(f"\n  {'Category':<30s} {'Count':>6s}  {'% of M6 hubs':>12s}")
    print(f"  {'-'*30} {'-'*6}  {'-'*12}")
    total = len(annotations)
    for cat in sorted(cat_counts.keys(), key=lambda c: -cat_counts[c]):
        count = cat_counts[cat]
        pct = 100 * count / total if total > 0 else 0
        print(f"  {cat:<30s} {count:>6d}  {pct:>11.1f}%")
    
    # Key interpretation
    print(f"\n  === Interpretation ===")
    if cat_counts.get("protein_degradation", 0) > 0:
        pct = 100 * cat_counts["protein_degradation"] / total
        print(f"  Protein degradation: {cat_counts['protein_degradation']} genes ({pct:.0f}%)")
        print(f"    → Consistent with 'metabolic cost' hypothesis:")
        print(f"      higher protein turnover → higher energy expenditure → lower net growth")
    if cat_counts.get("immune_response", 0) > 0:
        pct = 100 * cat_counts["immune_response"] / total
        print(f"  Immune response: {cat_counts['immune_response']} genes ({pct:.0f}%)")
        print(f"    → Consistent with 'immune surveillance cost' hypothesis")
    if cat_counts.get("stress_response", 0) > 0:
        pct = 100 * cat_counts["stress_response"] / total
        print(f"  Stress response: {cat_counts['stress_response']} genes ({pct:.0f}%)")
        print(f"    → Consistent with 'cellular stress' hypothesis")
    if cat_counts.get("other_uncharacterized", 0) > 0:
        pct = 100 * cat_counts["other_uncharacterized"] / total
        print(f"  Uncharacterized: {cat_counts['other_uncharacterized']} genes ({pct:.0f}%)")
        print(f"    → Expected for non-model organism with no annotated genome")
    
    # Write outputs
    ann_df.to_csv(out_dir / "m6_gene_annotations.tsv", sep="\t", index=False)
    
    # Functional summary
    summary_rows = []
    for cat in sorted(cat_counts.keys(), key=lambda c: -cat_counts[c]):
        summary_rows.append({
            "functional_category": cat,
            "gene_count": cat_counts[cat],
            "pct_of_M6_hubs": round(100 * cat_counts[cat] / total, 1),
            "example_genes": "; ".join(cat_genes[cat][:5]),
        })
    pd.DataFrame(summary_rows).to_csv(
        out_dir / "m6_functional_summary.tsv", sep="\t", index=False)
    
    print(f"\n  [OK] {out_dir}/m6_gene_annotations.tsv — {len(ann_df)} genes")
    print(f"  [OK] {out_dir}/m6_functional_summary.tsv — {len(summary_rows)} categories")
    print("\nDone.")


if __name__ == "__main__":
    main()
