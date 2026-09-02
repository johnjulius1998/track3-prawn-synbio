#!/usr/bin/env python3
"""
01_clr_transform.py — CLR-Normalize Pooled ASV Table
=====================================================
Track 3 Host-Microbe Integration (Approach 2)

INPUT:  ASV_table_Jumpers_Laggards.GKAQUA.csv
        (191 taxa, raw counts; n=1 pooled per group)

OUTPUT: taxa_direction.tsv with columns:
        taxon, jumper_raw, laggard_raw, clr_jumper, clr_laggard,
        fold_diff, direction, exclusive_to, contaminant_flag, contaminant_reason

CRITICAL RULES (enforced):
  - NO p-values, NO FDR, NO statistical tests (n=1 prohibits inference)
  - Direction is descriptive ONLY: "Jumper-enriched" / "Laggard-enriched"
  - Fold-difference is computed as CLR difference, not raw ratio
  - Contaminant screening uses pipeline/config/contaminant_genera.txt
  - All data originates ONLY from the supplied ASV table

USAGE:  python src/asv/01_clr_transform.py \
            --in data/raw/supplied/ASV_table_Jumpers_Laggards.GKAQUA.csv \
            --out data/processed/clr_profiles/taxa_direction.tsv \
            --min-reads 100 \
            --contaminants pipeline/config/contaminant_genera.txt
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def geometric_mean(series: pd.Series) -> float:
    """Geometric mean of non-zero values in a series."""
    nonzero = series[series > 0]
    if len(nonzero) == 0:
        return 0.0
    return np.exp(np.mean(np.log(nonzero)))


def clr_transform(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """
    Centered log-ratio (CLR) transformation.
    CLR(x_i) = ln(x_i / g(x)) where g(x) is the geometric mean of the sample.

    With n=1 per group, each column is treated as a single composition.
    """
    result = df.copy()
    for col in columns:
        gmean = geometric_mean(df[col])
        if gmean > 0:
            result[f"clr_{col}"] = np.log(df[col] / gmean)
        else:
            result[f"clr_{col}"] = np.nan
        # Handle zeros: CLR is undefined for zero. Use a small pseudocount
        # based on the minimum non-zero value in the table, or set to NaN.
        # Here we replace 0 with NaN in CLR space and flag them.
        result.loc[df[col] == 0, f"clr_{col}"] = np.nan
    return result


def load_contaminant_list(path: str) -> set[str]:
    """Load contaminant genera from config file. Returns set of lowercase genus names."""
    contaminants = set()
    if not Path(path).exists():
        print(f"[WARN] Contaminant list not found: {path}", file=sys.stderr)
        return contaminants
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Take first word (genus name), strip trailing comment
            genus = line.split("#")[0].strip().lower()
            if genus:
                contaminants.add(genus)
    return contaminants


def extract_genus(taxon: str) -> str:
    """Extract genus from 'Genus species' formatted taxon name."""
    return taxon.split()[0].lower() if taxon else ""


def classify_direction(
    df: pd.DataFrame, total_jumper: int, total_laggard: int
) -> pd.DataFrame:
    """
    Classify each taxon's direction using CLR fold-difference.
    Rules:
      - Direction = sign of (clr_jumper - clr_laggard)
      - If exclusive to one group: use "Jumper-exclusive" / "Laggard-exclusive"
      - If both zero: "undetermined"
      - fold_diff = clr_jumper - clr_laggard (NaN-safe)
    """
    # Compute CLR fold difference
    df["fold_diff"] = df["clr_Jumper"] - df["clr_Laggard"]

    conditions = []

    for idx, row in df.iterrows():
        j_raw = row["Jumper"]
        l_raw = row["Laggard"]
        fd = row["fold_diff"]

        if j_raw > 0 and l_raw == 0:
            conditions.append("Jumper-exclusive")
        elif l_raw > 0 and j_raw == 0:
            conditions.append("Laggard-exclusive")
        elif j_raw == 0 and l_raw == 0:
            conditions.append("absent")
        elif pd.notna(fd):
            if fd > 0:
                conditions.append("Jumper-enriched")
            elif fd < 0:
                conditions.append("Laggard-enriched")
            else:
                conditions.append("balanced")
        else:
            conditions.append("undetermined")

    df["direction"] = conditions
    df["exclusive_to"] = df["direction"].apply(
        lambda x: (
            "Jumper"
            if "Jumper" in str(x) and "exclusive" in str(x)
            else ("Laggard" if "Laggard" in str(x) and "exclusive" in str(x) else "none")
        )
    )
    return df


def flag_contaminants(df: pd.DataFrame, contaminant_genera: set[str]) -> pd.DataFrame:
    """Flag taxa whose genus appears in the contaminant list."""
    flags = []
    reasons = []

    for taxon in df["taxon"]:
        genus = extract_genus(taxon)
        if genus in contaminant_genera:
            flags.append(True)
            reasons.append(f"Genus '{genus}' is a known reagent/skin contaminant")
        else:
            flags.append(False)
            reasons.append("")

    df["contaminant_flag"] = flags
    df["contaminant_reason"] = reasons
    return df


def compute_contamination_score(df: pd.DataFrame, contaminant_genera: set[str]) -> pd.DataFrame:
    """
    Compute a continuous contamination score (0=clean, 1=highly suspect).

    Method: Two-component mixture model on log10(total_reads) combined with
    abundance percentile. The low-abundance component is treated as the
    contaminant source distribution.

    The score supplements, rather than replaces, the binary genus-level flag.
    Taxa with high scores AND in the contaminant list are the strongest suspects.
    Taxa with high scores NOT in the list may be stochastic low-abundance
    contaminants missed by the static list.

    References:
      - Davis et al. 2018, Microbiome 6:226 (decontam R package)
      - Salter et al. 2014, BMC Biology 12:87
      - Eisenhofer et al. 2019, Trends in Microbiology 27:117
    """
    from scipy.stats import norm as scipy_norm

    log_counts = np.log10(df["total_reads"].values + 1)
    median_log = np.median(log_counts)
    low_comp = log_counts <= median_log
    high_comp = log_counts > median_log

    mu_low = log_counts[low_comp].mean() if low_comp.any() else 0.0
    sigma_low = max(float(log_counts[low_comp].std()), 0.01) if low_comp.sum() > 1 else 0.1
    mu_high = log_counts[high_comp].mean() if high_comp.any() else 1.0
    sigma_high = max(float(log_counts[high_comp].std()), 0.01) if high_comp.sum() > 1 else 0.1

    p_low = max(float(low_comp.mean()), 0.01)
    loglik_low = scipy_norm.logpdf(np.clip(log_counts, -5, 5), mu_low, sigma_low)
    loglik_high = scipy_norm.logpdf(np.clip(log_counts, -5, 5), mu_high, sigma_high)
    log_joint_low = np.log(p_low) + loglik_low
    log_joint_high = np.log(1.0 - p_low) + loglik_high
    log_sum = np.logaddexp(log_joint_low, log_joint_high)
    score_mixture = np.exp(np.clip(log_joint_low - log_sum, -50, 0))

    # Combine mixture posterior with abundance percentile
    abund_percentile = df["total_reads"].rank(pct=True)
    score = 0.5 * score_mixture + 0.5 * (1.0 - abund_percentile)

    df["contamination_score"] = score.round(4)

    # Risk category based on score
    df["contamination_risk"] = df["contamination_score"].apply(
        lambda s: "HIGH" if s >= 0.85 else ("MEDIUM" if s >= 0.50 else "LOW")
    )

    # Detailed rationale string
    rationales = []
    for _, row in df.iterrows():
        parts = []
        genus_in_list = extract_genus(row["taxon"]) in contaminant_genera
        if genus_in_list:
            parts.append(f"Genus '{extract_genus(row['taxon'])}' in contaminant database")
        if row["contamination_score"] >= 0.85:
            parts.append(f"Low abundance (score={row['contamination_score']:.2f})")
        if row["total_reads"] <= 2:
            parts.append("Singleton/doubleton — possible stochastic contaminant")
        if not parts:
            parts.append("No contaminant indicators")
        rationales.append("; ".join(parts))

    df["contamination_rationale"] = rationales
    return df


def compute_relative_abundance(df: pd.DataFrame, total_j: int, total_l: int) -> pd.DataFrame:
    """Add relative abundance columns (percentage of total reads)."""
    df["rel_abund_jumper_pct"] = (df["Jumper"] / total_j) * 100
    df["rel_abund_laggard_pct"] = (df["Laggard"] / total_l) * 100
    return df


def main():
    parser = argparse.ArgumentParser(description="CLR-transform pooled ASV table")
    parser.add_argument("--in", dest="input_file", required=True, help="Input ASV CSV file")
    parser.add_argument("--out", dest="output_file", required=True, help="Output TSV file")
    parser.add_argument("--min-reads", type=int, default=100,
                        help="Minimum total reads to retain a taxon (default: 100)")
    parser.add_argument("--contaminants", default="pipeline/config/contaminant_genera.txt",
                        help="Path to contaminant genera list")
    args = parser.parse_args()

    # ---- Load data ----
    print(f"[CLR] Loading ASV table: {args.input_file}")
    raw = pd.read_csv(args.input_file)

    # Validate columns
    required_cols = {"ASV", "Jumper", "Laggard"}
    if not required_cols.issubset(raw.columns):
        missing = required_cols - set(raw.columns)
        sys.exit(f"ERROR: Missing columns: {missing}")

    # Rename ASV -> taxon for clarity
    raw = raw.rename(columns={"ASV": "taxon"})

    # Ensure numeric
    raw["Jumper"] = pd.to_numeric(raw["Jumper"], errors="coerce").fillna(0).astype(int)
    raw["Laggard"] = pd.to_numeric(raw["Laggard"], errors="coerce").fillna(0).astype(int)

    n_total = len(raw)
    total_jumper = raw["Jumper"].sum()
    total_laggard = raw["Laggard"].sum()
    print(f"[CLR] Loaded {n_total} taxa")
    print(f"[CLR] Total reads — Jumper: {total_jumper}, Laggard: {total_laggard}")

    # ---- Filter low-read taxa ----
    raw["total_reads"] = raw["Jumper"] + raw["Laggard"]
    before = len(raw)
    df = raw[raw["total_reads"] >= args.min_reads].copy()
    after = len(df)
    if before > after:
        print(f"[CLR] Filtered {before - after} taxa below {args.min_reads} total reads")
    else:
        print(f"[CLR] All taxa pass min-reads filter (>= {args.min_reads})")

    if len(df) == 0:
        sys.exit("ERROR: No taxa remain after filtering.")

    # ---- CLR Transform ----
    print("[CLR] Computing CLR transformation...")
    df = clr_transform(df, ["Jumper", "Laggard"])

    # ---- Classify direction ----
    print("[CLR] Classifying enrichment direction...")
    df = classify_direction(df, total_jumper, total_laggard)

    # ---- Relative abundance ----
    df = compute_relative_abundance(df, total_jumper, total_laggard)

    # ---- Contaminant screening (v2: binary + continuous score) ----
    contaminant_genera = load_contaminant_list(args.contaminants)
    print(f"[CLR] Loaded {len(contaminant_genera)} contaminant genera")
    df = flag_contaminants(df, contaminant_genera)
    df = compute_contamination_score(df, contaminant_genera)
    n_high = (df["contamination_risk"] == "HIGH").sum()
    n_med = (df["contamination_risk"] == "MEDIUM").sum()
    n_low = (df["contamination_risk"] == "LOW").sum()
    print(f"[CLR] Contamination risk — HIGH: {n_high}, MEDIUM: {n_med}, LOW: {n_low}")

    # ---- Sort by absolute fold-difference (Jumper-enriched first) ----
    df["abs_fold_diff"] = df["fold_diff"].abs()
    df = df.sort_values(["direction", "abs_fold_diff"], ascending=[True, False])

    # ---- Select & order output columns (v2: includes contamination_score) ----
    out_cols = [
        "taxon",
        "Jumper", "Laggard",
        "rel_abund_jumper_pct", "rel_abund_laggard_pct",
        "clr_Jumper", "clr_Laggard",
        "fold_diff", "direction", "exclusive_to",
        "contaminant_flag", "contamination_score", "contamination_risk",
        "contaminant_reason", "contamination_rationale",
    ]
    df_out = df[out_cols].copy()

    # Round for readability
    float_cols = ["rel_abund_jumper_pct", "rel_abund_laggard_pct",
                  "clr_Jumper", "clr_Laggard", "fold_diff"]
    for col in float_cols:
        if col in df_out.columns:
            df_out[col] = df_out[col].round(4)

    # ---- Write output ----
    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(args.output_file, sep="\t", index=False, na_rep="NA")
    print(f"[CLR] Output written: {args.output_file}")

    # ---- Summary ----
    n_jumper = (df["direction"].str.contains("Jumper")).sum()
    n_laggard = (df["direction"].str.contains("Laggard")).sum()
    n_contam = df["contaminant_flag"].sum()
    print(f"[CLR] Summary — Jumper-associated: {n_jumper}, "
          f"Laggard-associated: {n_laggard}, Flagged contaminants: {n_contam}")
    print("[CLR] Done.")


if __name__ == "__main__":
    main()
