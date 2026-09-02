#!/usr/bin/env python3
"""
03_functional_profiles.py — Genome-resolved functional profiling of pooled ASV taxa
====================================================================================
Track 3 Host-Microbe Integration

WHY THIS EXISTS
---------------
The supplied ASV table (ASV_table_Jumpers_Laggards.GKAQUA.csv) carries
species-level names only: no ASV IDs, no 16S sequences, and the Track 3
Data Note states raw FASTQs are not being released. A literal PICRUSt2 /
Tax4Fun2 run is therefore impossible (both require representative
sequences). The Data Note sanctions the alternative explicitly:
    "Functional inference from taxonomy does not require replication —
     a predicted functional profile can be derived for each pool and the
     two compared descriptively."

This script replaces the hand-curated genus→function map
(02_picrust_inference.py) with TOOL-DERIVED functional profiles:

  species name
    -> NCBI Taxonomy ID            (NCBI E-utilities esearch)
    -> KEGG GENOME entry           (KEGG REST, binomial name match)
    -> gene->KO links              (KEGG REST link/ko/{org})
    -> per-genome KO capacity      (copies normalised by total genes)
    -> pool-level functional indices (weighted by pool read counts)
    -> reference-pathway capacity  (KEGG REST link/pathway/ko: KO x map
                                    membership product, comparable across taxa)
    -> EC numbers                  (KEGG link/reaction/ko + link/ec/rn chain)

MetaCyc mapping is intentionally skipped: BioCyc web-service function
queries and EC frames are hCaptcha-gated (documented in
metacyc_out/README.txt). KEGG pathway + EC tables are the canonical outputs.

CRITICAL RULES (enforced)
-------------------------
  - n = 1 pooled sample per group: ALL outputs are DESCRIPTIVE indices
    and log2 ratios of indices. NO p-values, NO FDR, NO enrichment tests.
  - Every table documents which taxa contributed and how many lacked a
    KEGG genome (coverage is reported, not hidden).
  - All network responses are cached so the run is reproducible offline.

USAGE
-----
  python src/asv/03_functional_profiles.py \
      --asv data/raw/supplied/ASV_table_Jumpers_Laggards.GKAQUA.csv \
      --out-dir data/interim/functional_prediction \
      --final-dir data/processed/clr_profiles \
      [--skip-metacyc] \
      [--cache-dir data/interim/literature_tables/api_cache]
"""

import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import pandas as pd

USER_AGENT = "Mozilla/5.0 (research; track3_host_microbe)"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
KEGG = "https://rest.kegg.jp"
BIOCYC = "https://websvc.biocyc.org"

SLEEP_EUTILS = 0.35   # NCBI asks <=3 requests/sec without an API key
SLEEP_KEGG = 0.20
SLEEP_BIOCYC = 0.30


# ----------------------------------------------------------------------
# HTTP helpers (with on-disk cache for reproducibility)
# ----------------------------------------------------------------------

class ApiCache:
    """Persistent GET cache keyed by sha256 of URL."""

    def __init__(self, cache_dir: Path):
        self.dir = Path(cache_dir)
        self.dir.mkdir(parents=True, exist_ok=True)

    def _path(self, url: str) -> Path:
        return self.dir / (hashlib.sha256(url.encode()).hexdigest() + ".txt")

    def get(self, url: str) -> str | None:
        p = self._path(url)
        if p.exists():
            return p.read_text(encoding="utf-8", errors="replace")
        return None

    def put(self, url: str, text: str) -> None:
        self._path(url).write_text(text, encoding="utf-8")


def http_get(url: str, cache: ApiCache, sleep_s: float = 0.2,
             retries: int = 4, timeout: int = 60) -> str:
    """GET with cache, retry + backoff, and polite throttling."""
    cached = cache.get(url)
    if cached is not None:
        return cached
    import urllib.request
    time.sleep(sleep_s)
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                text = resp.read().decode("utf-8", errors="replace")
            cache.put(url, text)
            return text
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"GET failed after {retries} tries: {url} ({last_err})")


# ----------------------------------------------------------------------
# Step 1: species name -> NCBI Taxonomy ID
# ----------------------------------------------------------------------

def resolve_taxid(name: str, cache: ApiCache) -> int | None:
    term = name.replace(" ", "+")
    url = (f"{EUTILS}/esearch.fcgi?db=taxonomy&term={term}%5BScientific%20Name%5D"
           f"&retmode=json&retmax=3")
    text = http_get(url, cache, sleep_s=SLEEP_EUTILS)
    try:
        payload = json.loads(text)
        ids = payload["esearchresult"].get("idlist", [])
    except Exception:
        # sometimes esearch returns XML-ish errors; parse fallback
        ids = re.findall(r"<Id>(\d+)</Id>", text)
    if ids:
        return int(ids[0])
    # fallback: search without field qualifier
    term2 = name.replace(" ", "+")
    url2 = (f"{EUTILS}/esearch.fcgi?db=taxonomy&term={term2}&retmode=json&retmax=3")
    text2 = http_get(url2, cache, sleep_s=SLEEP_EUTILS)
    try:
        ids = json.loads(text2)["esearchresult"].get("idlist", [])
    except Exception:
        ids = re.findall(r"<Id>(\d+)</Id>", text2)
    return int(ids[0]) if ids else None


# ----------------------------------------------------------------------
# Step 2: KEGG GENOME list -> binomial-name index
# ----------------------------------------------------------------------

def load_kegg_genomes(cache: ApiCache) -> dict[str, list[tuple[str, str]]]:
    """
    Download `list genome` once and index by lowercase binomial
    ('genus species') -> list of (org_code, full_name).
    """
    text = http_get(f"{KEGG}/list/genome", cache, sleep_s=SLEEP_KEGG)
    index: dict[str, list[tuple[str, str]]] = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        rest = parts[1].strip()
        if ";" not in rest:
            continue
        org, name = rest.split(";", 1)
        org, name = org.strip(), name.strip()
        if not org or not name:
            continue
        tokens = name.split()
        if len(tokens) >= 2:
            binomial = f"{tokens[0].lower()} {tokens[1].lower()}"
            index.setdefault(binomial, []).append((org, name))
        elif tokens:
            index.setdefault(tokens[0].lower(), []).append((org, name))
    return index


def match_kegg_org(taxon: str, genome_index: dict) -> tuple[str | None, str | None]:
    """Match 'Genus species [strain]' to a KEGG organism code."""
    tokens = taxon.strip().split()
    if len(tokens) >= 2:
        key = f"{tokens[0].lower()} {tokens[1].lower()}"
    elif tokens:
        key = tokens[0].lower()
    else:
        return None, None
    if key in genome_index:
        return genome_index[key][0]
    # fuzzy: genus match only, first candidate
    if len(tokens) >= 2:
        genus = tokens[0].lower()
        for k, v in genome_index.items():
            if k.startswith(genus + " "):
                return v[0]
    return None, None


# ----------------------------------------------------------------------
# Step 3: KEGG link operations
# ----------------------------------------------------------------------

def kegg_link(org: str, target: str, cache: ApiCache) -> list[tuple[str, str]]:
    """rest.kegg.jp/link/{target}/{org} -> list of (gene_id, target_id)."""
    url = f"{KEGG}/link/{target}/{org}"
    text = http_get(url, cache, sleep_s=SLEEP_KEGG)
    pairs = []
    for line in text.splitlines():
        if "\t" not in line:
            continue
        g, t = line.split("\t", 1)
        pairs.append((g.strip(), t.strip()))
    return pairs


def kegg_list(db: str, cache: ApiCache) -> list[tuple[str, str]]:
    """rest.kegg.jp/list/{db} -> list of (id, name)."""
    text = http_get(f"{KEGG}/list/{db}", cache, sleep_s=SLEEP_KEGG)
    out = []
    for line in text.splitlines():
        if "\t" not in line:
            continue
        i, n = line.split("\t", 1)
        out.append((i.strip(), n.strip()))
    return out


# ----------------------------------------------------------------------
# Step 4 (optional): MetaCyc pathway mapping via BioCyc web service
# ----------------------------------------------------------------------

def ec_to_metacyc_pathways(ec: str, cache: ApiCache) -> list[str]:
    """Best-effort: EC number -> MetaCyc pathway frame IDs (BioCyc API)."""
    if not re.fullmatch(r"\d+\.\d+\.\d+\.\d+", ec):
        return []
    results = []
    for frame in (f"META%3A%3AEC-{ec}", f"META%3A%3A{ec}"):
        url = f"{BIOCYC}/apixml?fn=pathways-of-reaction&id={frame}&detail=none"
        try:
            text = http_get(url, cache, sleep_s=SLEEP_BIOCYC, timeout=40)
        except Exception:
            continue
        if "<Error>" in text or "pathways-of-reaction" not in text:
            continue
        hits = re.findall(r"frameid=\"([^\"]+)\"", text)
        if hits:
            results = hits
            break
    return results


def metacyc_pathway_name(pw_id: str, cache: ApiCache) -> str:
    url = (f"{BIOCYC}/getxml?META%3A%3A{pw_id}"
           f"&detail=low")
    try:
        text = http_get(url, cache, sleep_s=SLEEP_BIOCYC, timeout=40)
    except Exception:
        return ""
    m = re.search(r"<common-name>(.*?)</common-name>", text, re.S)
    return m.group(1).strip() if m else ""


# ----------------------------------------------------------------------
# Aggregation
# ----------------------------------------------------------------------

def aggregate_pool_index(taxa_rows: list[dict], capacity: dict,
                         pool: str) -> float:
    """Sum over taxa: rel_abundance(t,pool) * capacity(t)."""
    total = 0.0
    for row in taxa_rows:
        rel = row["reads"][pool] / row["total"][pool]
        cap = capacity.get(row["taxon"], 0.0)
        total += rel * cap
    return total


def main():
    ap = argparse.ArgumentParser(description="Genome-resolved functional profiling (KEGG GENOME)")
    ap.add_argument("--asv", required=True, help="Supplied pooled ASV table (CSV)")
    ap.add_argument("--out-dir", default="data/interim/functional_prediction")
    ap.add_argument("--final-dir", default="data/processed/clr_profiles")
    ap.add_argument("--skip-metacyc", action="store_true")
    ap.add_argument("--cache-dir",
                    default="data/interim/literature_tables/api_cache")
    args = ap.parse_args()

    base = Path(__file__).resolve().parents[3]
    asv_path = Path(args.asv)
    out_dir = Path(args.out_dir)
    final_dir = Path(args.final_dir)
    if not out_dir.is_absolute():
        out_dir = base / out_dir
    if not final_dir.is_absolute():
        final_dir = base / final_dir
    kegg_out = out_dir / "kegg_out"
    metacyc_out = out_dir / "metacyc_out"
    kegg_out.mkdir(parents=True, exist_ok=True)
    metacyc_out.mkdir(parents=True, exist_ok=True)
    final_dir.mkdir(parents=True, exist_ok=True)

    cache = ApiCache(Path(args.cache_dir) if Path(args.cache_dir).is_absolute()
                     else base / args.cache_dir)

    print("[FUNC] Loading ASV table:", asv_path)
    asv = pd.read_csv(asv_path)
    taxa = [(str(t).strip(), int(j), int(l))
            for t, j, l in zip(asv["ASV"], asv["Jumper"], asv["Laggard"])]
    taxa = [t for t in taxa if t[0].lower() not in ("unknown", "nan", "")]
    total_j = sum(t[1] for t in taxa)
    total_l = sum(t[2] for t in taxa)
    print(f"[FUNC] {len(taxa)} named taxa | Jumper total={total_j} | "
          f"Laggard total={total_l}")

    # ---- resolve taxids + KEGG genomes (network) ----
    print("[FUNC] Resolving species -> KEGG GENOME ...")
    genome_index = load_kegg_genomes(cache)
    print(f"[FUNC] KEGG genome index: {len(genome_index)} binomial keys")

    rows = []           # per-taxon genome summary
    ko_capacity = {}    # taxon -> {ko: copies}
    ko_genes = {}       # taxon -> total genes
    path_capacity = {}  # taxon -> {map: gene count}

    for i, (taxon, jr, lr) in enumerate(taxa, 1):
        print(f"[FUNC] ({i}/{len(taxa)}) {taxon}", flush=True)
        taxid = resolve_taxid(taxon, cache)
        org, kegg_name = match_kegg_org(taxon, genome_index)
        rec = {
            "taxon": taxon,
            "ncbi_taxid": taxid or "",
            "kegg_org": org or "",
            "kegg_name": kegg_name or "",
            "total_genes": 0,
            "n_ko": 0,
            "n_pathways": 0,
            "status": "",
        }
        if org is None:
            rec["status"] = "no_kegg_genome"
            rows.append(rec)
            ko_capacity[taxon] = {}
            path_capacity[taxon] = {}
            ko_genes[taxon] = 0
            continue
        try:
            ko_links = kegg_link(org, "ko", cache)
        except Exception as exc:  # noqa: BLE001
            print(f"[FUNC]   KEGG link failed for {org}: {exc}")
            rec["status"] = "kegg_link_failed"
            rows.append(rec)
            ko_capacity[taxon] = {}
            path_capacity[taxon] = {}
            ko_genes[taxon] = 0
            continue

        loci = set()
        ko_counts: dict[str, int] = {}
        for g, ko in ko_links:
            if ko.startswith("ko:"):
                loci.add(g)
                ko_counts[ko] = ko_counts.get(ko, 0) + 1

        rec.update({
            "total_genes": len(loci),
            "n_ko": len(ko_counts),
            "n_pathways": 0,  # filled after global KO->map membership load
            "status": "matched",
        })
        rows.append(rec)
        ko_capacity[taxon] = ko_counts
        path_capacity[taxon] = {}
        ko_genes[taxon] = len(loci)

    # Global KO -> reference pathway membership (one request, cached).
    # /link/pathway/{org} returns ORGANISM-specific path ids, which are not
    # comparable across taxa; /link/pathway/ko returns KO -> path:mapXXXXX
    # membership, which is. Pathway capacity per taxon is therefore the
    # KO-profile x KO->map membership product.
    print("[FUNC] Loading KO -> KEGG reference pathway membership ...")
    ko_to_maps: dict[str, set[str]] = {}
    text = http_get(f"{KEGG}/link/pathway/ko", cache, sleep_s=SLEEP_KEGG)
    for line in text.splitlines():
        if "\t" in line:
            ko, pth = line.split("\t", 1)
            if ko.startswith("ko:") and pth.startswith("path:map"):
                ko_to_maps.setdefault(ko, set()).add(pth)
    print(f"[FUNC] {len(ko_to_maps)} KOs mapped to reference pathways")

    for rec, (taxon, _jr, _lr) in zip(rows, taxa):
        if rec["status"] != "matched":
            continue
        pc_ap: dict[str, int] = {}
        for ko, copies in ko_capacity[taxon].items():
            for pth in ko_to_maps.get(ko, ()):
                pc_ap[pth] = pc_ap.get(pth, 0) + copies
        path_capacity[taxon] = pc_ap
        rec["n_pathways"] = len(pc_ap)

    mapping_df = pd.DataFrame(rows)
    mapping_df.to_csv(kegg_out / "taxon_genome_mapping.tsv", sep="\t", index=False)
    n_matched = (mapping_df["status"] == "matched").sum()
    n_reads_matched_j = sum(t[1] for t, r in zip(taxa, rows) if r["status"] == "matched")
    n_reads_matched_l = sum(t[2] for t, r in zip(taxa, rows) if r["status"] == "matched")
    print(f"[FUNC] KEGG GENOME matched: {n_matched}/{len(taxa)} taxa "
          f"(reads covered: {n_reads_matched_j}/{total_j} Jumper, "
          f"{n_reads_matched_l}/{total_l} Laggard)")

    # ---- descriptive indices (NO statistics) ----
    taxa_rows = [{"taxon": t[0],
                  "reads": {"Jumper": t[1], "Laggard": t[2]},
                  "total": {"Jumper": total_j, "Laggard": total_l}}
                 for t in taxa]

    def pool_capacity(taxon, cap, ngenes):
        return {ko: cnt / ngenes for ko, cnt in cap.items()} if ngenes else {}

    print("[FUNC] Aggregating pool-level functional indices ...")
    ko_names = dict(kegg_list("ko", cache))
    path_names = {f"path:{k}": v for k, v in kegg_list("pathway", cache)}

    # KO index per pool
    ko_rows = []
    all_kos = sorted({ko for cap in ko_capacity.values() for ko in cap})
    for ko in all_kos:
        j_idx = 0.0
        l_idx = 0.0
        for tr in taxa_rows:
            cap = pool_capacity(tr["taxon"], ko_capacity.get(tr["taxon"], {}),
                                ko_genes.get(tr["taxon"], 0))
            rel_j = tr["reads"]["Jumper"] / tr["total"]["Jumper"]
            rel_l = tr["reads"]["Laggard"] / tr["total"]["Laggard"]
            j_idx += rel_j * cap.get(ko, 0.0)
            l_idx += rel_l * cap.get(ko, 0.0)
        ko_rows.append({
            "ko": ko,
            "ko_name": ko_names.get(ko, ""),
            "jumper_idx": j_idx,
            "laggard_idx": l_idx,
        })
    ko_df = pd.DataFrame(ko_rows)
    pc = 1e-6
    if len(ko_df):
        ko_df["log2_ratio_idx"] = (ko_df["jumper_idx"] + pc) / (ko_df["laggard_idx"] + pc)
        ko_df["log2_ratio_idx"] = ko_df["log2_ratio_idx"].map(lambda x: __import__("math").log2(x))
    else:
        ko_df = pd.DataFrame(columns=["ko", "ko_name", "jumper_idx",
                                      "laggard_idx", "log2_ratio_idx"])
    ko_df.to_csv(kegg_out / "ko_abundance.tsv", sep="\t", index=False)
    print(f"[FUNC] {len(ko_df)} KOs profiled")

    # Pathway index per pool
    path_rows = []
    contrib_rows = []
    all_paths = sorted({p for cap in path_capacity.values() for p in cap})
    for pth in all_paths:
        j_idx = 0.0
        l_idx = 0.0
        j_taxa = []
        l_taxa = []
        for tr in taxa_rows:
            cap = pool_capacity(tr["taxon"], path_capacity.get(tr["taxon"], {}),
                                ko_genes.get(tr["taxon"], 0))
            rel_j = tr["reads"]["Jumper"] / tr["total"]["Jumper"]
            rel_l = tr["reads"]["Laggard"] / tr["total"]["Laggard"]
            cj = rel_j * cap.get(pth, 0.0)
            cl = rel_l * cap.get(pth, 0.0)
            if cj > 0:
                j_taxa.append((tr["taxon"], cj))
            if cl > 0:
                l_taxa.append((tr["taxon"], cl))
            if cj > 0 or cl > 0:
                contrib_rows.append({
                    "pathway_id": pth.replace("path:", ""),
                    "taxon": tr["taxon"],
                    "jumper_contrib": cj,
                    "laggard_contrib": cl,
                })
            j_idx += cj
            l_idx += cl
        j_taxa.sort(key=lambda x: -x[1])
        l_taxa.sort(key=lambda x: -x[1])
        path_rows.append({
            "pathway_id": pth.replace("path:", ""),
            "pathway_name": path_names.get(pth, ""),
            "jumper_idx": j_idx,
            "laggard_idx": l_idx,
            "n_contributing_taxa_jumper": len(j_taxa),
            "n_contributing_taxa_laggard": len(l_taxa),
            "top_jumper_taxa": "; ".join(f"{t} ({v:.2e})" for t, v in j_taxa[:5]),
            "top_laggard_taxa": "; ".join(f"{t} ({v:.2e})" for t, v in l_taxa[:5]),
        })
    path_df = pd.DataFrame(path_rows)
    if len(path_df):
        path_df["log2_ratio_idx"] = ((path_df["jumper_idx"] + pc)
                                     / (path_df["laggard_idx"] + pc))
        path_df["log2_ratio_idx"] = path_df["log2_ratio_idx"].map(
            lambda x: __import__("math").log2(x))
    else:
        path_df = pd.DataFrame(columns=["pathway_id", "pathway_name",
                                        "jumper_idx", "laggard_idx",
                                        "log2_ratio_idx",
                                        "n_contributing_taxa_jumper",
                                        "n_contributing_taxa_laggard",
                                        "top_jumper_taxa", "top_laggard_taxa"])
    path_df.to_csv(kegg_out / "pathway_abundance.tsv", sep="\t", index=False)
    contrib_df = pd.DataFrame(contrib_rows)
    if len(contrib_df):
        contrib_df.to_csv(kegg_out / "pathway_taxon_contribution.tsv",
                          sep="\t", index=False)
    print(f"[FUNC] {len(path_df)} KEGG pathways profiled "
          f"({len(contrib_df)} pathway x taxon contributions)")

    # EC index via the KO -> reaction -> EC chain (global tables, cached)
    print("[FUNC] Mapping KO -> EC (via reactions) ...")
    ko_to_ecs: dict[str, set[str]] = {}
    try:
        ko_to_rns: dict[str, set[str]] = {}
        text = http_get(f"{KEGG}/link/reaction/ko", cache, sleep_s=SLEEP_KEGG)
        for line in text.splitlines():
            if "\t" in line:
                ko, rn = line.split("\t", 1)
                if ko.startswith("ko:") and rn.startswith("rn:"):
                    ko_to_rns.setdefault(ko, set()).add(rn)
        rn_to_ecs: dict[str, set[str]] = {}
        text2 = http_get(f"{KEGG}/link/ec/rn", cache, sleep_s=SLEEP_KEGG)
        for line in text2.splitlines():
            if "\t" in line:
                rn, ec = line.split("\t", 1)
                if rn.startswith("rn:") and ec.startswith("ec:"):
                    rn_to_ecs.setdefault(rn, set()).add(ec)
        for ko, rns in ko_to_rns.items():
            for rn in rns:
                for ec in rn_to_ecs.get(rn, ()):
                    ko_to_ecs.setdefault(ko, set()).add(ec)
        print(f"[FUNC] {len(ko_to_ecs)} KOs have EC numbers")
    except Exception as exc:  # noqa: BLE001
        print(f"[FUNC] EC mapping failed ({exc}) — EC table skipped")

    ec_rows = {}
    for ko in all_kos:
        for ec in ko_to_ecs.get(ko, []):
            j_idx = 0.0
            l_idx = 0.0
            for tr in taxa_rows:
                cap = pool_capacity(tr["taxon"], ko_capacity.get(tr["taxon"], {}),
                                    ko_genes.get(tr["taxon"], 0))
                rel_j = tr["reads"]["Jumper"] / tr["total"]["Jumper"]
                rel_l = tr["reads"]["Laggard"] / tr["total"]["Laggard"]
                j_idx += rel_j * cap.get(ko, 0.0)
                l_idx += rel_l * cap.get(ko, 0.0)
            if ec not in ec_rows:
                ec_rows[ec] = {"ec": ec, "jumper_idx": 0.0, "laggard_idx": 0.0}
            ec_rows[ec]["jumper_idx"] += j_idx
            ec_rows[ec]["laggard_idx"] += l_idx
    ec_df = pd.DataFrame(list(ec_rows.values()))
    if len(ec_df):
        ec_df["log2_ratio_idx"] = ((ec_df["jumper_idx"] + pc)
                                   / (ec_df["laggard_idx"] + pc))
        ec_df["log2_ratio_idx"] = ec_df["log2_ratio_idx"].map(
            lambda x: __import__("math").log2(x))
    ec_df.to_csv(kegg_out / "ec_abundance.tsv", sep="\t", index=False)
    print(f"[FUNC] {len(ec_df)} EC numbers profiled")

    # ---- MetaCyc ----
    # BioCyc web-service function queries (apixml fn=...) and EC-number frames
    # are hCaptcha-gated (verified 2026-08-16: getxml for META:PWY-* works but
    # EC/reaction frames and pathways-of-reaction return an HTML challenge).
    # MetaCyc mapping is therefore skipped and documented here; the KEGG
    # pathway + EC tables are the canonical outputs.
    metacyc_done = False
    note = (
        "MetaCyc pathway mapping not available for this run.\n"
        "Reason: the BioCyc web service gates function queries and EC-number\n"
        "frames behind an hCaptcha challenge (verified 2026-08-16). KEGG\n"
        "pathway (map IDs) and EC profiles are provided instead via KEGG\n"
        "GENOME -> KO -> reference-pathway/EC mapping.\n"
    )
    (metacyc_out / "README.txt").write_text(note)

    # ---- final deliverables in processed/ ----
    final_path = path_df.copy()
    final_path["source"] = "KEGG_GENOME"
    final_path["method_note"] = ("Genome-resolved functional index: "
                                 "sum over taxa of (pool rel. abundance x "
                                 "pathway genes / genome total genes). "
                                 "DESCRIPTIVE ONLY (n=1 pooled per group; no "
                                 "statistical inference).")
    cols = ["pathway_id", "pathway_name", "source", "jumper_idx", "laggard_idx",
            "log2_ratio_idx", "n_contributing_taxa_jumper",
            "n_contributing_taxa_laggard", "top_jumper_taxa",
            "top_laggard_taxa", "method_note"]
    final_path = final_path[[c for c in cols if c in final_path.columns]]
    final_path.to_csv(final_dir / "pathway_abundance.tsv", sep="\t", index=False)

    final_ko = ko_df.copy()
    final_ko["source"] = "KEGG_GENOME"
    final_ko.to_csv(final_dir / "ko_abundance.tsv", sep="\t", index=False)

    print(f"[FUNC] Final outputs:")
    print(f"[FUNC]   {final_dir}/pathway_abundance.tsv ({len(final_path)} rows)")
    print(f"[FUNC]   {final_dir}/ko_abundance.tsv ({len(final_ko)} rows)")
    print(f"[FUNC]   {kegg_out}/taxon_genome_mapping.tsv "
          f"({n_matched}/{len(taxa)} taxa matched)")
    print("[FUNC] MetaCyc: skipped (BioCyc API gated; see metacyc_out/README.txt)")
    print("[FUNC] DONE — no p-values were computed (n=1 design).")


if __name__ == "__main__":
    main()
