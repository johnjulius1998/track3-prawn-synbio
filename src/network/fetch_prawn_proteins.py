#!/usr/bin/env python3
"""
fetch_prawn_proteins.py — Build the prawn protein set for orthology mapping
============================================================================
Track 3 Host-Microbe Integration (M6 GSEA support)

WHY: the M. rosenbergii reference transcriptome is unannotated and hub genes
are NCBI nucleotide accessions only. GSEA needs proteins to map against a
reference proteome (Swiss-Prot). This script produces, for every gene in the
merged counts matrix:

  1. RefSeq mRNAs (XM_/XR_): the annotated CDS protein via NCBI efetch
     `rettype=fasta_cds_na` (headers pair XM_ -> XP_), translated locally.
  2. EST/TSA/GenBank accessions: the deposited nucleotide sequence,
     translated in-silico using the longest open reading frame.

OUTPUTS:
  data/interim/literature_tables/prawn_proteins.faa   (>{gene_id})
  data/interim/literature_tables/prawn_gene_protein_map.tsv
      (gene_id, protein_acc, source, n_aa, orf_frame, status)

USAGE:
  python src/network/fetch_prawn_proteins.py \
      --counts data/processed/gene_expression/merged_counts.tsv \
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
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SLEEP_S = 0.35  # NCBI asks <=3 req/sec without API key
BATCH = 200

REFSEQ_RNA_PREFIXES = ("XM_", "XR_")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(seq: str) -> str:
    seq = seq.upper()
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X")
                   for i in range(0, len(seq) - 2, 3))


def longest_orf(seq: str) -> tuple[str, int]:
    """Longest ATG...stop ORF across the three forward frames."""
    best, best_frame = "", -1
    for frame in range(3):
        prot = translate(seq[frame:])
        # split at stop codons; keep ORFs starting with M
        orfs = prot.split("*")
        # consider ORFs that begin with M and, if none, any ORF
        for start_ok in (True, False):
            for orf in orfs:
                if start_ok and not orf.startswith("M"):
                    continue
                if len(orf) > len(best):
                    best, best_frame = orf, frame
            if best:
                break
        if best:
            break
    if not best:
        # no stop codon: translate whole frame
        prot = translate(seq[0:])
        best = prot.rstrip("*") or prot
        best_frame = 0
    return best, best_frame


class FileCache:
    def __init__(self, cache_dir: Path):
        self.dir = Path(cache_dir)
        self.dir.mkdir(parents=True, exist_ok=True)

    def path(self, url: str) -> Path:
        return self.dir / (hashlib.sha256(url.encode()).hexdigest() + ".txt")

    def get(self, url: str) -> str | None:
        p = self.path(url)
        return p.read_text(encoding="utf-8", errors="replace") if p.exists() else None

    def put(self, url: str, text: str) -> None:
        self.path(url).write_text(text, encoding="utf-8")


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


def parse_fasta(text: str) -> list[tuple[str, str]]:
    recs, header, seq = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                recs.append((header, "".join(seq)))
            header, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if header is not None:
        recs.append((header, "".join(seq)))
    return recs


def fetch_refseq_cds(ids: list[str], cache: FileCache) -> dict[str, tuple[str, str]]:
    """XM_ -> (XP_ accession, protein seq) via fasta_cds_na."""
    out: dict[str, tuple[str, str]] = {}
    url = (f"{EUTILS}/efetch.fcgi?db=nuccore&id={','.join(ids)}"
           f"&rettype=fasta_cds_na&retmode=text")
    text = http_get(url, cache)
    for header, seq in parse_fasta(text):
        xm = re.search(r"(XM_\d+\.\d+|XR_\d+\.\d+)", header)
        xp = re.search(r"(XP_\d+\.\d+|YP_\d+\.\d+)", header)
        if xm:
            out[xm.group(1)] = (xp.group(1) if xp else "", seq)
    return out


def fetch_nucleotide_fasta(ids: list[str], cache: FileCache) -> dict[str, str]:
    """Any accession -> deposited nucleotide sequence."""
    out: dict[str, str] = {}
    url = (f"{EUTILS}/efetch.fcgi?db=nuccore&id={','.join(ids)}"
           f"&rettype=fasta&retmode=text")
    text = http_get(url, cache)
    for header, seq in parse_fasta(text):
        acc = header.split()[0]
        out[acc] = seq
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--out-dir", default="data/interim/literature_tables")
    ap.add_argument("--cache-dir",
                    default="data/interim/literature_tables/api_cache")
    args = ap.parse_args()

    base = Path(__file__).resolve().parents[2]
    counts_path = Path(args.counts)
    out_dir = Path(args.out_dir)
    if not counts_path.is_absolute():
        counts_path = base / counts_path
    if not out_dir.is_absolute():
        out_dir = base / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = FileCache(Path(args.cache_dir) if Path(args.cache_dir).is_absolute()
                      else base / args.cache_dir)

    print("[PROT] Loading counts:", counts_path)
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    genes = list(counts.index.astype(str))
    print(f"[PROT] {len(genes)} genes")

    refseq = [g for g in genes if g.startswith(REFSEQ_RNA_PREFIXES)]
    other = [g for g in genes if not g.startswith(REFSEQ_RNA_PREFIXES)]
    print(f"[PROT] RefSeq mRNAs: {len(refseq)} | other accessions: {len(other)}")

    results: dict[str, dict] = {}  # gene -> {protein_acc, seq, source, orf_frame, status}

    def batched(ids, n=BATCH):
        return [ids[i:i + n] for i in range(0, len(ids), n)]

    print("[PROT] Fetching RefSeq CDS proteins ...")
    for b in batched(refseq):
        try:
            got = fetch_refseq_cds(b, cache)
        except Exception as exc:  # noqa: BLE001
            print(f"[PROT]   batch failed: {exc}")
            got = {}
        for g in b:
            if g in got and got[g][1]:
                prot, frame = longest_orf(got[g][1])  # translate CDS nuc -> protein
                results[g] = {
                    "protein_acc": got[g][0], "seq": prot,
                    "source": "refseq_cds", "orf_frame": frame,
                    "status": "ok" if len(prot) >= 20 else "short_orf",
                }
            else:
                results[g] = {
                    "protein_acc": "", "seq": "",
                    "source": "refseq_cds", "orf_frame": -1,
                    "status": "no_cds_retrieved",
                }
        print(f"[PROT]   {min(b[-1], refseq[-1]) if b else ''} "
              f"... {len(results)} processed", flush=True)

    print("[PROT] Fetching other accessions (nucleotide) ...")
    for b in batched(other):
        try:
            got = fetch_nucleotide_fasta(b, cache)
        except Exception as exc:  # noqa: BLE001
            print(f"[PROT]   batch failed: {exc}")
            got = {}
        for g in b:
            seq = got.get(g, "")
            if seq:
                prot, frame = longest_orf(seq)
                results[g] = {
                    "protein_acc": "", "seq": prot,
                    "source": "longest_orf", "orf_frame": frame,
                    "status": "ok" if len(prot) >= 20 else "short_orf",
                }
            else:
                results[g] = {
                    "protein_acc": "", "seq": "", "source": "longest_orf",
                    "orf_frame": -1, "status": "no_sequence",
                }
        print(f"[PROT]   {len(results)} processed", flush=True)

    # fill any refseq with no CDS by falling back to nucleotide+ORF
    for g in refseq:
        if results.get(g, {}).get("status") == "no_cds_retrieved":
            try:
                got = fetch_nucleotide_fasta([g], cache)
            except Exception:  # noqa: BLE001
                got = {}
            seq = got.get(g, "")
            if seq:
                prot, frame = longest_orf(seq)
                results[g] = {
                    "protein_acc": "", "seq": prot, "source": "longest_orf",
                    "orf_frame": frame, "status": "ok" if len(prot) >= 20 else "short_orf",
                }
    print("[PROT] RefSeq fallback ORF pass done")

    # ---- write outputs ----
    faa = out_dir / "prawn_proteins.faa"
    n_written = 0
    with open(faa, "w") as f:
        for g in genes:
            seq = results[g]["seq"]
            if seq:
                f.write(f">{g}\n")
                for i in range(0, len(seq), 70):
                    f.write(seq[i:i + 70] + "\n")
                n_written += 1
    print(f"[PROT] Wrote {n_written}/{len(genes)} proteins -> {faa}")

    rows = [{
        "gene_id": g,
        "protein_acc": results[g].get("protein_acc", ""),
        "source": results[g].get("source", ""),
        "orf_frame": results[g].get("orf_frame", ""),
        "n_aa": len(results[g].get("seq", "")),
        "status": results[g].get("status", "missing"),
    } for g in genes]
    map_df = pd.DataFrame(rows)
    map_df.to_csv(out_dir / "prawn_gene_protein_map.tsv", sep="\t", index=False)
    print(f"[PROT] Map written: {out_dir}/prawn_gene_protein_map.tsv")
    print("[PROT] Status counts:")
    print(map_df["status"].value_counts().to_string())
    print("[PROT] DONE")


if __name__ == "__main__":
    main()
