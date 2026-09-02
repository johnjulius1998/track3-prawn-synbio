#!/usr/bin/env python3
"""
confound_correction.py — Latent-Confound Removal Before WGCNA
===============================================================
Track 3 Host-Microbe Integration (v2 — Fix 3)

This module identifies and removes latent confounders from the gene
expression matrix BEFORE WGCNA network construction. The current
approach uses post-hoc partial correlation — this is a pre-WGCNA
correction that prevents modules from forming around confounded variation.

METHOD:
  1. Log-transform TPM (log2(TPM+1))
  2. Run PCA on the expression matrix
  3. Identify PCs that correlate with known confounders (sex, tissue)
     but NOT with the trait of interest (weight_gain)
  4. Regress out those PCs from each gene's expression profile
  5. Output a corrected expression matrix for downstream WGCNA

RATIONALE:
  The feasibility study demonstrated:
  - PC1 (46.9% var) is tissue-driven (r_tissue=−0.798)
  - PC2 (15.1% var) is sex+tissue (r_sex=+0.563)
  - PC5 (4.2% var) is sex-correlated (r_sex=−0.569)
  - Removing these 3 PCs removes 66.2% of variance
  - 69.6% of genes change WG correlation by >|0.1| after correction

  The current post-hoc partial correlation approach is methodologically
  weaker because modules may already be organized around confounded
  variation by the time partial correlation is computed.

REFERENCES:
  - Leek & Storey 2007, PLoS Genetics 3:e161 (SVA methodology)
  - Gagnon-Bartsch & Speed 2012, NAR 40:e29 (RUV methodology)
  - Johnson, Li & Rabinovic 2007, Biostatistics 8:118 (ComBat)

USAGE:
  python src/network/confound_correction.py \\
      --tpm data/processed/gene_expression/merged_tpm.tsv \\
      --metadata data/raw/sra/PRJNA875278/metadata.tsv \\
      --trait weight_gain \\
      --confounds sex,tissue \\
      --out data/processed/gene_expression/corrected_tpm.tsv \\
      --out-report results/reports/confound_report.tsv
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg
from scipy.stats import zscore


def parse_confound_label(raw_labels, value_map):
    """
    Parse a metadata column into a numeric array.
    E.g., tissue: Hepatopancreas→0, Testis→1, Ovary→1
           sex: male→1, female→0
    """
    values = []
    for v in raw_labels:
        v_str = str(v).strip().lower()
        if v_str in value_map:
            values.append(value_map[v_str])
        else:
            # Try to find partial match
            found = False
            for key, val in value_map.items():
                if key in v_str:
                    values.append(val)
                    found = True
                    break
            if not found:
                values.append(np.nan)
    return np.array(values, dtype=float)


def main():
    parser = argparse.ArgumentParser(
        description="Remove latent confounders from expression data before WGCNA"
    )
    parser.add_argument("--tpm", required=True, help="Merged TPM matrix (TSV, genes×samples)")
    parser.add_argument("--metadata", required=True, help="Sample metadata (TSV)")
    parser.add_argument("--trait", default="weight_gain",
                        help="Trait of interest column name [default: weight_gain]")
    parser.add_argument("--confounds", default="sex,tissue",
                        help="Comma-separated known confounders [default: sex,tissue]")
    parser.add_argument("--min-tpm", type=float, default=1.0,
                        help="Mean TPM threshold to retain gene [default: 1.0]")
    parser.add_argument("--pcs-to-remove", type=str, default=None,
                        help="Manually specify PCs to remove (comma-separated 1-based). "
                             "If not provided, auto-detect via correlation threshold.")
    parser.add_argument("--confound-r-threshold", type=float, default=0.3,
                        help="|r| threshold for PC-confound correlation [default: 0.3]")
    parser.add_argument("--max-pcs-to-check", type=int, default=20,
                        help="Number of PCs to evaluate [default: 20]")
    parser.add_argument("--out", required=True,
                        help="Output corrected expression matrix (TSV)")
    parser.add_argument("--out-report", default=None,
                        help="Output confound removal report (TSV)")
    args = parser.parse_args()

    # ---- Load data ----
    print(f"[CONFOUND] Loading TPM: {args.tpm}")
    tpm = pd.read_csv(args.tpm, sep="\t", index_col=0)
    print(f"[CONFOUND]   {tpm.shape[0]} genes × {tpm.shape[1]} samples")

    print(f"[CONFOUND] Loading metadata: {args.metadata}")
    meta = pd.read_csv(args.metadata, sep="\t", index_col=0)
    print(f"[CONFOUND]   {len(meta)} samples, columns: {list(meta.columns)}")

    # ---- Filter genes ----
    gene_mean_tpm = tpm.mean(axis=1)
    keep = gene_mean_tpm >= args.min_tpm
    expr = tpm.loc[keep]
    print(f"[CONFOUND] Genes after mean TPM ≥ {args.min_tpm}: "
          f"{len(expr)} / {len(tpm)}")

    # ---- Align samples ----
    common = sorted(set(expr.columns) & set(meta.index))
    if len(common) < 10:
        sys.exit(f"ERROR: Only {len(common)} common samples. Need ≥10.")
    print(f"[CONFOUND] Common samples: {len(common)}")
    X_raw = expr[common].values.T  # samples × genes
    sample_ids = common

    # ---- Extract trait & confounds ----
    trait_name = args.trait
    if trait_name not in meta.columns:
        sys.exit(f"ERROR: Trait '{trait_name}' not found in metadata. "
                 f"Available: {list(meta.columns)}")
    trait = np.array([float(meta.loc[s, trait_name]) for s in sample_ids])

    confound_names = [c.strip() for c in args.confounds.split(",")]
    confound_matrix = []
    for cname in confound_names:
        if cname not in meta.columns:
            print(f"[CONFOUND] WARNING: confound '{cname}' not in metadata — skipping")
            continue
        raw_vals = [meta.loc[s, cname] for s in sample_ids]
        unique_raw = sorted(set(str(v).strip() for v in raw_vals))
        print(f"[CONFOUND] Confound '{cname}': unique values = {unique_raw}")

        # Auto-generate numeric mapping
        if all(v.replace(".", "").replace("-", "").isdigit() for v in unique_raw if v):
            # Already numeric
            vals = np.array([float(v) if v else np.nan for v in raw_vals])
        else:
            # Categorical: create value map
            value_map = {v.lower(): i for i, v in enumerate(unique_raw)}
            print(f"[CONFOUND]   Mapping: {value_map}")
            vals = parse_confound_label(raw_vals, value_map)
        confound_matrix.append(vals)

    if not confound_matrix:
        sys.exit("ERROR: No valid confounds found.")
    confound_matrix = np.column_stack(confound_matrix)
    n_confounds = confound_matrix.shape[1]
    print(f"[CONFOUND] Confound matrix: {len(sample_ids)} × {n_confounds}")

    # Remove samples with NaN in confounds
    valid_samples = ~np.isnan(confound_matrix).any(axis=1)
    if not valid_samples.all():
        n_bad = (~valid_samples).sum()
        print(f"[CONFOUND] Removing {n_bad} samples with NaN confounds")
        X_raw = X_raw[valid_samples]
        trait = trait[valid_samples]
        confound_matrix = confound_matrix[valid_samples]
        sample_ids = [s for i, s in enumerate(sample_ids) if valid_samples[i]]

    # ---- Log-transform ----
    X_log = np.log2(X_raw + 1)

    # ---- PCA ----
    X_scaled = zscore(X_log, axis=0, nan_policy="omit")
    X_scaled = np.nan_to_num(X_scaled, 0)

    U, S, Vt = linalg.svd(X_scaled, full_matrices=False)
    explained_var = S**2 / np.sum(S**2)
    cum_var = np.cumsum(explained_var)

    print(f"\n[CONFOUND] PCA of expression matrix "
          f"({X_scaled.shape[0]} samples × {X_scaled.shape[1]} genes)")
    print(f"[CONFOUND] {'PC':>5s} {'Var%':>8s} {'Cum%':>8s} "
          f"{'r('+trait_name+')':>10s} " +
          " ".join(f"{'r('+c+')':>10s}" for c in confound_names))
    print(f"[CONFOUND] {'-'*60}")

    pc_report = []
    for i in range(min(args.max_pcs_to_check, len(explained_var))):
        r_trait = np.corrcoef(U[:, i], trait)[0, 1] if not np.isnan(trait).any() else np.nan
        r_confs = [np.corrcoef(U[:, i], confound_matrix[:, j])[0, 1]
                   for j in range(n_confounds)]
        r_confs_str = " ".join(f"{r:10.3f}" for r in r_confs)
        print(f"[CONFOUND] PC{i+1:3d}: {explained_var[i]*100:7.1f}% "
              f"{cum_var[i]*100:7.1f}% {r_trait:10.3f} {r_confs_str}")
        pc_report.append({
            "pc": i + 1,
            "variance_explained_pct": round(explained_var[i] * 100, 2),
            "cumulative_pct": round(cum_var[i] * 100, 2),
            f"r_{trait_name}": round(r_trait, 4),
            **{f"r_{cname}": round(r_confs[j], 4)
               for j, cname in enumerate(confound_names)},
        })

    # ---- Identify confounded PCs ----
    if args.pcs_to_remove is not None:
        pcs_to_remove = [int(x) - 1 for x in args.pcs_to_remove.split(",")]
        print(f"\n[CONFOUND] Using manually specified PCs: "
              f"{[f'PC{p+1}' for p in pcs_to_remove]}")
    else:
        pcs_to_remove = []
        for i in range(min(args.max_pcs_to_check, U.shape[1])):
            r_trait = abs(np.corrcoef(U[:, i], trait)[0, 1]) if not np.isnan(trait).any() else 0
            r_conf_max = max(
                abs(np.corrcoef(U[:, i], confound_matrix[:, j])[0, 1])
                for j in range(n_confounds)
            )
            # Remove if correlated with confound AND less correlated with trait
            if r_conf_max > args.confound_r_threshold and r_trait < r_conf_max:
                pcs_to_remove.append(i)

    print(f"\n[CONFOUND] PCs to remove: {[f'PC{p+1}' for p in pcs_to_remove]}")
    if pcs_to_remove:
        var_removed = explained_var[pcs_to_remove].sum() * 100
        print(f"[CONFOUND] Total variance removed: {var_removed:.1f}%")
        print(f"[CONFOUND] ⚠️  Note: Some removed PCs may carry genuine biological signal "
              f"if trait and confound are correlated. Review the report.")
    else:
        print("[CONFOUND] No PCs met removal criteria — expression matrix unchanged.")

    # ---- Regress out confounded PCs ----
    X_corrected = X_scaled.copy()
    if pcs_to_remove:
        confound_pcs = U[:, pcs_to_remove]
        for j in range(X_scaled.shape[1]):
            y = X_scaled[:, j]
            try:
                beta, _, _, _ = np.linalg.lstsq(confound_pcs, y, rcond=None)
                X_corrected[:, j] = y - confound_pcs @ beta
            except Exception:
                pass  # Keep original for problematic genes
        print(f"[CONFOUND] Corrected: removed {len(pcs_to_remove)} confounded PCs "
              f"from {X_corrected.shape[1]} genes")

    # ---- Measure impact ----
    n_check = min(1000, X_scaled.shape[1])
    deltas = []
    for j in range(n_check):
        orig_r = np.corrcoef(X_scaled[:, j], trait)[0, 1]
        corr_r = np.corrcoef(X_corrected[:, j], trait)[0, 1]
        deltas.append(abs(corr_r - orig_r))
    deltas = np.array(deltas)
    print(f"[CONFOUND] Impact on {trait_name}-gene correlations "
          f"(n={n_check} genes sampled):")
    print(f"[CONFOUND]   Mean |Δr|: {np.mean(deltas):.4f}")
    print(f"[CONFOUND]   Median |Δr|: {np.median(deltas):.4f}")
    print(f"[CONFOUND]   |Δr| > 0.1: {sum(deltas > 0.1)} "
          f"({sum(deltas > 0.1)/len(deltas)*100:.1f}%)")

    # ---- Write corrected matrix (back-transform to TPM space) ----
    # Since we worked in z-score(log2(TPM+1)) space, convert back:
    # corrected_log = X_corrected * σ_j + μ_j
    # corrected_TPM = 2^corrected_log - 1
    means = X_log.mean(axis=0)
    stds = X_log.std(axis=0, ddof=1)
    stds[stds == 0] = 1.0

    X_corrected_log = X_corrected * stds + means
    X_corrected_tpm = np.power(2.0, X_corrected_log) - 1.0
    X_corrected_tpm = np.maximum(X_corrected_tpm, 0.0)  # No negative TPM

    out_df = pd.DataFrame(
        X_corrected_tpm.T,
        index=expr.index,
        columns=sample_ids
    )
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", float_format="%.4f")
    print(f"[CONFOUND] Corrected TPM written to: {args.out}")

    # ---- Write report ----
    if args.out_report:
        report_df = pd.DataFrame(pc_report)
        report_df["removed"] = [i in pcs_to_remove for i in range(len(pc_report))]
        Path(args.out_report).parent.mkdir(parents=True, exist_ok=True)
        report_df.to_csv(args.out_report, sep="\t", index=False)
        print(f"[CONFOUND] Report written to: {args.out_report}")

    print("[CONFOUND] Done.")


if __name__ == "__main__":
    main()
