#!/usr/bin/env python3
"""
map_uniprot_sets.py — Swiss-Prot best hits -> KEGG / GO gene sets
==================================================================
Track 3 Host-Microbe Integration (M6 GSEA support)

Takes DIAMOND blastp best hits of prawn proteins against Swiss-Prot and
builds organism-agnostic gene sets for GSEA:

  prawn gene -> best Swiss-Prot hit
             -> UniProt xref_kegg (hsa:xxxx)  -> KEGG pathway (KEGG REST)
             -> UniProt go_p (GO BP terms)    -> GO term sets

OUTPUTS (data/interim/literature_tables/):
  prawn_gene_uniprot_map.tsv   gene, uniprot_acc, kegg_ids, pident, evalue, bitscore
  uniprot_gene_sets_kegg.tsv   gene, pathway_id, pathway_name
  uniprot_gene_sets_go.tsv     gene, go_id

USAGE:
  python src/network/map_uniprot_sets.py \
      --diamond data/interim/literature_tables/diamond_swissprot.tsv \
      --out-dir data/interim/literature_tables \
      [--cache-dir data/interim/literature_tables/api_cache]
"""

import argparse
import hashlib
import re
import time
from pathlib import Path

import pandas as pd

USER_AGENT = "Mozilla/5.0 (research; track3_host_microbe)"
UNIPROT = "https://rest.uniprot.org"
KEGG = "https://rest.kegg.jp"
SLEEP_S = 0.25

COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


class FileCache:
    def __init__(self, cache_dir: Path):
        self.dir = Path(cache_dir)
        self.dir.mkdir(parents=True, exist_ok=True)

    def _path(self, url: str) -> Path:
        return self.dir / (hashlib.sha256(url.encode()).hexdigest() + ".txt")

    def get(self, url: str) -> str | None:
        p = self._path(url)
        return p.read_text(encoding="utf-8", errors="replace") if p.exists() else None

    def put(self, url: str, text: str) -> None:
        self._path(url).write_text(text, encoding="utf-8")


def http_get(url: str, cache: FileCache, retries: int = 4, timeout: int = 120) -> str:
    cached = cache.get(url)
    if cached is not None:
        return cached
    import urllib.request
    time.sleep(SLEEP_S)
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
    raise RuntimeError(f"GET failed: {url} ({last_err})")


def fetch_uniprot_fields(accessions: list[str], cache: FileCache) -> pd.DataFrame:
    """Batch-fetch accession,xref_kegg,go_p,ec for Swiss-Prot accessions."""
    frames = []
    for i in range(0, len(accessions), 100):
        batch = accessions[i:i + 100]
        url = (f"{UNIPROT}/uniprotkb/accessions?accessions={','.join(batch)}"
               f"&fields=accession,xref_kegg,go_p&format=tsv")
        text = http_get(url, cache)
        lines = text.splitlines()
        if len(lines) < 2:
            continue
        header = lines[0].split("\t")
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) >= len(header):
                frames.append(dict(zip(header, parts)))
    up = pd.DataFrame(frames)
    # UniProt TSV uses display names for the fields
    up = up.rename(columns={
        "Entry": "accession",
        "KEGG": "xref_kegg",
        "Gene Ontology (biological process)": "go_p",
    })
    for col in ("accession", "xref_kegg", "go_p"):
        if col not in up.columns:
            up[col] = None
    return up


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diamond", required=True,
                    help="DIAMOND blastp outfmt6 TSV")
    ap.add_argument("--out-dir", default="data/interim/literature_tables")
    ap.add_argument("--cache-dir",
                    default="data/interim/literature_tables/api_cache")
    args = ap.parse_args()

    base = Path(__file__).resolve().parents[2]
    diamond_path = Path(args.diamond)
    out_dir = Path(args.out_dir)
    if not diamond_path.is_absolute():
        diamond_path = base / diamond_path
    if not out_dir.is_absolute():
        out_dir = base / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = FileCache(Path(args.cache_dir) if Path(args.cache_dir).is_absolute()
                      else base / args.cache_dir)

    print("[SETS] Loading DIAMOND output:", diamond_path)
    blast = pd.read_csv(diamond_path, sep="\t", header=None, names=COLS)
    print(f"[SETS] {len(blast)} alignments for {blast['qseqid'].nunique()} queries")

    # best hit per query by bitscore
    blast = blast.sort_values("bitscore", ascending=False)
    best = blast.drop_duplicates("qseqid", keep="first").copy()
    best["uniprot_acc"] = best["sseqid"].str.split("|").str[1]
    print(f"[SETS] {len(best)} unique best hits")

    accs = sorted(best["uniprot_acc"].dropna().unique().tolist())
    print(f"[SETS] Fetching UniProt fields for {len(accs)} accessions ...")
    up = fetch_uniprot_fields(accs, cache)
    print(f"[SETS] UniProt returned {len(up)} records")
    if len(up) == 0:
        print("[SETS] ERROR: no UniProt annotations returned — check network")
        return 1

    best = best.merge(up.rename(columns={"accession": "uniprot_acc"}),
                      on="uniprot_acc", how="left")
    print(f"[SETS] {best['xref_kegg'].notna().sum()} hits carry KEGG cross-refs")

    # ---- KEGG pathway sets ----
    # collect hsa gene ids
    def parse_hsa(cell):
        if pd.isna(cell):
            return []
        return re.findall(r"hsa:(\d{4,6})", str(cell))

    gene_to_hsa: dict[str, list[str]] = {}
    for gene, cell in zip(best["qseqid"], best["xref_kegg"]):
        ids = parse_hsa(cell)
        if ids:
            gene_to_hsa[gene] = ids

    hsa_ids = sorted({i for ids in gene_to_hsa.values() for i in ids})
    print(f"[SETS] {len(hsa_ids)} human KEGG gene ids from {len(gene_to_hsa)} genes")

    # /link/pathway/hsa:{id} returns organism-specific path ids; instead use
    # the global tables: hsa gene -> KO (one request) then KO -> map (one
    # request, same table as src/asv/03_functional_profiles.py).
    hsa_to_kos: dict[str, set[str]] = {}
    ko_to_maps: dict[str, set[str]] = {}
    if hsa_ids:
        text = http_get(f"{KEGG}/link/ko/hsa", cache)
        for line in text.splitlines():
            if "\t" in line:
                g, ko = line.split("\t", 1)
                if g.startswith("hsa:") and ko.startswith("ko:"):
                    hsa_to_kos.setdefault(g.strip().split(":")[-1], set()).add(ko)
        text2 = http_get(f"{KEGG}/link/pathway/ko", cache)
        for line in text2.splitlines():
            if "\t" in line:
                ko, pth = line.split("\t", 1)
                if ko.startswith("ko:") and pth.startswith("path:map"):
                    ko_to_maps.setdefault(ko, set()).add(pth)

    hsa_to_path: dict[str, set[str]] = {}
    for h in hsa_ids:
        maps = set()
        for ko in hsa_to_kos.get(h, ()):
            maps |= ko_to_maps.get(ko, set())
        if maps:
            hsa_to_path[h] = maps
    print(f"[SETS] {len(hsa_to_path)} hsa genes mapped to KEGG map pathways")

    # pathway names
    path_names: dict[str, str] = {}
    if hsa_to_path:
        text = http_get(f"{KEGG}/list/pathway", cache)
        for line in text.splitlines():
            if "\t" in line:
                pid, name = line.split("\t", 1)
                if pid.startswith("path:map"):
                    path_names[pid] = name.strip()

    kegg_rows = []
    for gene, ids in gene_to_hsa.items():
        seen = set()
        for h in ids:
            for pid in hsa_to_path.get(h, ()):
                if pid not in seen:
                    seen.add(pid)
                    kegg_rows.append({
                        "gene": gene,
                        "pathway_id": pid.replace("path:", ""),
                        "pathway_name": path_names.get(pid, ""),
                    })
    kegg_df = pd.DataFrame(kegg_rows)
    if len(kegg_df):
        kegg_df.to_csv(out_dir / "uniprot_gene_sets_kegg.tsv", sep="\t", index=False)
        print(f"[SETS] KEGG gene sets: {len(kegg_df)} rows, "
              f"{kegg_df['gene'].nunique()} genes, "
              f"{kegg_df['pathway_id'].nunique()} pathways")
    else:
        kegg_df = pd.DataFrame(columns=["gene", "pathway_id", "pathway_name"])
        kegg_df.to_csv(out_dir / "uniprot_gene_sets_kegg.tsv", sep="\t", index=False)
        print("[SETS] KEGG gene sets: 0 rows (no hsa->map mapping)")

    # ---- GO BP sets ----
    def parse_go(cell):
        if pd.isna(cell):
            return []
        return re.findall(r"GO:\d{7}", str(cell))

    go_rows = []
    for gene, cell in zip(best["qseqid"], best["go_p"]):
        for go in set(parse_go(cell)):
            go_rows.append({"gene": gene, "go_id": go})
    go_df = pd.DataFrame(go_rows)
    if len(go_df):
        go_df.to_csv(out_dir / "uniprot_gene_sets_go.tsv", sep="\t", index=False)
        print(f"[SETS] GO BP gene sets: {len(go_df)} rows, "
              f"{go_df['gene'].nunique()} genes, {go_df['go_id'].nunique()} terms")
    else:
        go_df = pd.DataFrame(columns=["gene", "go_id"])
        go_df.to_csv(out_dir / "uniprot_gene_sets_go.tsv", sep="\t", index=False)
        print("[SETS] GO BP gene sets: 0 rows")

    # ---- mapping table ----
    map_df = best[["qseqid", "uniprot_acc", "pident", "evalue", "bitscore"]].copy()
    map_df.columns = ["gene", "uniprot_acc", "pident", "evalue", "bitscore"]
    map_df["kegg_ids"] = [
        ",".join(gene_to_hsa.get(g, [])) for g in map_df["gene"]]
    map_df.to_csv(out_dir / "prawn_gene_uniprot_map.tsv", sep="\t", index=False)
    print(f"[SETS] Mapping written: {out_dir}/prawn_gene_uniprot_map.tsv")
    print("[SETS] DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
