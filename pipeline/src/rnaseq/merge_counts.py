#!/usr/bin/env python3
"""
merge_counts.py — Merge Salmon quant.sf Files into Counts Matrix
==================================================================
Track 3 Host-Microbe Integration

Reads all quant.sf files from data/interim/host_counts/*/quant.sf,
extracts TPM and NumReads, and produces:
  - merged_tpm.tsv     (genes × samples, TPM values)
  - merged_counts.tsv  (genes × samples, estimated read counts)

Also filters: removes genes with zero expression in ALL samples,
and genes with mean TPM < min_tpm.

USAGE: python src/rnaseq/merge_counts.py \
           --quant-dir data/interim/host_counts \
           --out-dir data/processed/gene_expression \
           --metadata data/raw/sra/PRJNA875278/metadata.tsv \
           --min-tpm 1.0
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np


def load_quant(quant_file: Path) -> pd.Series:
    """Load a single quant.sf file, return TPM and NumReads as Series."""
    df = pd.read_csv(quant_file, sep="\t", usecols=["Name", "TPM", "NumReads"])
    df = df.set_index("Name")
    return df["TPM"], df["NumReads"]


def main():
    parser = argparse.ArgumentParser(description="Merge Salmon quant outputs")
    parser.add_argument("--quant-dir", default="data/interim/host_counts",
                        help="Directory containing per-sample quant folders")
    parser.add_argument("--out-dir", default="data/processed/gene_expression",
                        help="Output directory for merged matrices")
    parser.add_argument("--metadata", default="data/raw/sra/PRJNA875278/metadata.tsv",
                        help="Sample metadata TSV")
    parser.add_argument("--min-tpm", type=float, default=1.0,
                        help="Minimum mean TPM to retain a gene (default: 1.0)")
    parser.add_argument("--min-samples", type=int, default=3,
                        help="Minimum number of samples a gene must be expressed in")
    args = parser.parse_args()

    quant_dir = Path(args.quant_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Discover samples ----
    sample_dirs = sorted(quant_dir.glob("*"))
    samples = []
    for d in sample_dirs:
        qf = d / "quant.sf"
        if qf.exists():
            samples.append((d.name, qf))

    if not samples:
        sys.exit("ERROR: No quant.sf files found. Run Salmon quant first.")

    print(f"[MERGE] Found {len(samples)} quantified samples:")
    for name, _ in samples:
        print(f"  {name}")

    # ---- Load all samples ----
    tpm_data = {}
    count_data = {}
    genes_all = set()

    for name, qf in samples:
        tpm, counts = load_quant(qf)
        tpm_data[name] = tpm
        count_data[name] = counts
        genes_all.update(tpm.index.tolist())

    genes_all = sorted(genes_all)
    print(f"[MERGE] Total unique genes: {len(genes_all):,}")

    # ---- Build matrices ----
    tpm_matrix = pd.DataFrame(index=genes_all, columns=[s[0] for s in samples], dtype=float)
    count_matrix = pd.DataFrame(index=genes_all, columns=[s[0] for s in samples], dtype=float)

    for name, _ in samples:
        tpm_matrix[name] = tpm_data[name].reindex(genes_all).fillna(0.0)
        count_matrix[name] = count_data[name].reindex(genes_all).fillna(0.0)

    # ---- Filter low-expression genes ----
    # Keep genes with mean TPM >= min_tpm AND expressed in >= min_samples
    mean_tpm = tpm_matrix.mean(axis=1)
    n_expressed = (tpm_matrix > 0).sum(axis=1)

    keep = (mean_tpm >= args.min_tpm) & (n_expressed >= args.min_samples)

    tpm_filtered = tpm_matrix.loc[keep]
    count_filtered = count_matrix.loc[keep]

    print(f"[MERGE] Genes retained after filtering: {len(tpm_filtered):,} / {len(genes_all):,} "
          f"({100*len(tpm_filtered)/len(genes_all):.1f}%)")
    print(f"[MERGE]   Filter: mean TPM >= {args.min_tpm}, expressed in >= {args.min_samples} samples")

    # ---- Write ----
    tpm_out = out_dir / "merged_tpm.tsv"
    count_out = out_dir / "merged_counts.tsv"

    tpm_filtered.to_csv(tpm_out, sep="\t", float_format="%.4f")
    count_filtered.to_csv(count_out, sep="\t", float_format="%.2f")

    print(f"[MERGE] TPM matrix: {tpm_out} ({tpm_filtered.shape[0]} genes × {tpm_filtered.shape[1]} samples)")
    print(f"[MERGE] Counts matrix: {count_out} ({count_filtered.shape[0]} genes × {count_filtered.shape[1]} samples)")

    # ---- Quick QC ----
    print(f"\n[MERGE] === Quick QC ===")
    print(f"  TPM range: {tpm_filtered.min().min():.2f} – {tpm_filtered.max().max():.2f}")
    print(f"  TPM median (per gene): {tpm_filtered.median(axis=1).median():.2f}")
    print(f"  Total counts (millions): {count_filtered.sum(axis=0).div(1e6).round(1).to_dict()}")

    # ---- Check metadata alignment ----
    meta_path = Path(args.metadata)
    if meta_path.exists():
        meta = pd.read_csv(meta_path, sep="\t", index_col=0)
        common = sorted(set(tpm_filtered.columns) & set(meta.index))
        missing_in_meta = sorted(set(tpm_filtered.columns) - set(meta.index))
        missing_in_counts = sorted(set(meta.index) - set(tpm_filtered.columns))

        print(f"\n[MERGE] === Metadata alignment ===")
        print(f"  Samples in counts: {len(tpm_filtered.columns)}")
        print(f"  Samples in metadata: {len(meta)}")
        print(f"  Common: {len(common)}")
        if missing_in_meta:
            print(f"  WARNING: In counts but NOT in metadata: {missing_in_meta}")
        if missing_in_counts:
            print(f"  WARNING: In metadata but NOT in counts: {missing_in_counts}")

        # Save aligned metadata
        meta_aligned = meta.loc[common] if common else meta
        meta_out = out_dir / "metadata_aligned.tsv"
        meta_aligned.to_csv(meta_out, sep="\t")
        print(f"  Aligned metadata: {meta_out}")

    print("[MERGE] Done.")


if __name__ == "__main__":
    main()
