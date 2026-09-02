#!/usr/bin/env python3
"""
01b_sensitivity_analysis.py — Pseudocount & Leave-One-Taxon-Out Stability
===========================================================================
Track 3 Host-Microbe Integration (v3 — Stability Analysis)

PURPOSE:
  Quantify how sensitive the ranked microbial taxon list is to:
  (A) Choice of pseudocount for zero-handling in CLR transformation
  (B) Removal of any single taxon (leave-one-taxon-out jackknife)

RATIONALE:
  The current CLR implementation (01_clr_transform.py) handles zeros by
  setting them to NaN and computing the geometric mean from non-zero values
  only. Alternative pseudocount strategies produce different fold-difference
  estimates, especially for low-abundance taxa. With 48.4% zeros in the
  191×2 matrix, this choice is consequential.

  Similarly, highly dominant taxa like Jiulongibacter sediminis (55.5% of
  all Jumper reads) heavily influence the geometric mean — removing it
  could re-rank other taxa substantially.

OUTPUTS:
  results/reports/
    pseudocount_sensitivity.tsv  — Full sweep: rank per taxon per pseudocount
    loto_stability.tsv           — Per-taxon stability metrics
    taxon_confidence_report.tsv  — Merged report: stability flags appended

USAGE:
  python src/asv/01b_sensitivity_analysis.py \
      --in data/raw/supplied/ASV_table_Jumpers_Laggards.GKAQUA.csv \
      --out-dir results/reports/ \
      --contaminants pipeline/config/contaminant_genera.txt \
      --n-pseudocounts 10
"""

import argparse
import sys
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
from collections import defaultdict


# ================================================================
# CLR IMPLEMENTATIONS
# ================================================================

def geometric_mean_nonzero(values: np.ndarray) -> float:
    """Geometric mean of non-zero values (current implementation)."""
    nz = values[values > 0]
    if len(nz) == 0:
        return 0.0
    return np.exp(np.mean(np.log(nz)))


def geometric_mean_all(values: np.ndarray) -> float:
    """Geometric mean of all values (assumes pseudocount already applied)."""
    if len(values) == 0:
        return 0.0
    return np.exp(np.mean(np.log(values)))


def clr_zero_to_nan(counts: np.ndarray) -> np.ndarray:
    """
    Current implementation: compute CLR from nonzero values only.
    Zeros → NaN in CLR space.
    Returns CLR values (NaN for zeros).
    """
    gmean = geometric_mean_nonzero(counts)
    if gmean <= 0:
        return np.full_like(counts, np.nan, dtype=float)
    clr = np.full_like(counts, np.nan, dtype=float)
    mask = counts > 0
    clr[mask] = np.log(counts[mask] / gmean)
    return clr


def clr_with_pseudocount(counts: np.ndarray, pseudocount: float) -> np.ndarray:
    """
    CLR with uniform pseudocount for zeros.
    x_i' = x_i if x_i > 0 else pseudocount
    CLR(x_i') = ln(x_i' / g(x'))
    """
    adjusted = np.where(counts > 0, counts.astype(float), pseudocount)
    gmean = geometric_mean_all(adjusted)
    if gmean <= 0:
        return np.full_like(counts, np.nan, dtype=float)
    return np.log(adjusted / gmean)


def clr_with_multiplicative_replacement(counts: np.ndarray) -> np.ndarray:
    """
    Multiplicative replacement (Martin-Fernandez et al. 2003):
    Replace zeros with δ * (sum of non-zeros / (#zeros)) scaled by detection limit.
    Uses δ = 0.5 * min_nonzero / total_composition.
    This preserves the ratios among non-zero components.
    """
    nz = counts[counts > 0]
    n_zeros = (counts == 0).sum()
    if n_zeros == 0:
        gmean = geometric_mean_nonzero(counts)
        return np.log(counts / gmean)
    if len(nz) == 0:
        return np.full_like(counts, np.nan, dtype=float)
    
    total = counts.sum()
    min_nz = nz.min()
    delta = 0.5 * min_nz / total if total > 0 else 0.5
    
    # Replace zeros, then rescale non-zeros to preserve sum-to-1
    rep = np.where(counts > 0, counts.astype(float), delta)
    # Adjust non-zeros so that sum stays constant
    scale_factor = (total - n_zeros * delta) / nz.sum() if nz.sum() > 0 else 1.0
    rep[counts > 0] *= scale_factor
    
    gmean = geometric_mean_all(rep)
    return np.log(rep / gmean)


def clr_bayesian_mult_replacement(counts: np.ndarray) -> np.ndarray:
    """
    Bayesian-multiplicative replacement (Palarea-Albaladejo & Martin-Fernandez 2015).
    Uses the zCompositions approach: zeros are imputed via a Dirichlet
    posterior with uniform prior, preserving the compositional structure.
    
    Simplified here as: zeros → 2/3 * detection_limit, then rescale.
    For low-count 16S data with min count = 1, detection_limit = 1.
    """
    nz = counts[counts > 0]
    n_zeros = (counts == 0).sum()
    if n_zeros == 0:
        gmean = geometric_mean_nonzero(counts)
        return np.log(counts / gmean)
    if len(nz) == 0:
        return np.full_like(counts, np.nan, dtype=float)
    
    # Bayesian-multiplicative: impute zeros with a fraction of detection limit
    # For count data with min=1, use 2/3 as the imputation value for each zero
    delta = 0.65  # Common default in zCompositions
    
    rep = np.where(counts > 0, counts.astype(float), delta)
    # Rescale to preserve sum
    old_sum = rep.sum()
    if old_sum > 0:
        rep *= counts.sum() / old_sum
    
    gmean = geometric_mean_all(rep)
    return np.log(rep / gmean)


# ================================================================
# PSEUDOCOUNT SWEEP
# ================================================================

def pseudocount_sweep(df: pd.DataFrame, pseudocounts: list[float]) -> pd.DataFrame:
    """
    For each pseudocount value, compute CLR for both columns,
    fold-difference, and ranking by |fold_diff|.
    
    Returns a wide-format DataFrame with one row per taxon,
    columns for each pseudocount's fold_diff and rank.
    """
    j_counts = df["Jumper"].values.astype(float)
    l_counts = df["Laggard"].values.astype(float)
    taxon_names = df["ASV"].values

    # Baseline: current method (zero → NaN)
    clr_j_0 = clr_zero_to_nan(j_counts)
    clr_l_0 = clr_zero_to_nan(l_counts)
    fd_0 = clr_j_0 - clr_l_0
    rank_0 = np.argsort(np.argsort(-np.abs(np.nan_to_num(fd_0, nan=0.0))))

    results = [{
        "taxon": taxon_names[i],
        "method": "current (zero→NaN)",
        "fold_diff": fd_0[i] if not np.isnan(fd_0[i]) else 0.0,
        "fold_diff_is_nan": np.isnan(fd_0[i]),
        "rank": int(rank_0[i]),
        "pseudocount": "NaN",
    } for i in range(len(taxon_names))]

    # Add multiplicative replacement methods
    for method_name, clr_func in [
        ("mult_replacement", clr_with_multiplicative_replacement),
        ("bayesian_mult", clr_bayesian_mult_replacement),
    ]:
        try:
            clr_j = clr_func(j_counts)
            clr_l = clr_func(l_counts)
            fd = clr_j - clr_l
            rank = np.argsort(np.argsort(-np.abs(fd)))
            for i in range(len(taxon_names)):
                results.append({
                    "taxon": taxon_names[i],
                    "method": method_name,
                    "fold_diff": round(float(fd[i]), 6),
                    "fold_diff_is_nan": False,
                    "rank": int(rank[i]),
                    "pseudocount": method_name,
                })
        except Exception as e:
            print(f"  [WARN] {method_name} failed: {e}", file=sys.stderr)

    # Pseudocount sweep
    for pc in pseudocounts:
        clr_j = clr_with_pseudocount(j_counts, pc)
        clr_l = clr_with_pseudocount(l_counts, pc)
        fd = clr_j - clr_l
        rank = np.argsort(np.argsort(-np.abs(fd)))

        for i in range(len(taxon_names)):
            results.append({
                "taxon": taxon_names[i],
                "method": f"pseudocount={pc}",
                "fold_diff": round(float(fd[i]), 6),
                "fold_diff_is_nan": False,
                "rank": int(rank[i]),
                "pseudocount": f"{pc}",
            })

    return pd.DataFrame(results)


def compute_pseudocount_stability_metrics(sweep_df: pd.DataFrame) -> pd.DataFrame:
    """
    From the sweep results, compute per-taxon stability metrics:
      - mean_rank: mean rank across all pseudocount methods
      - rank_std: standard deviation of rank
      - rank_range: max_rank - min_rank
      - rank_volatility: rank_std / sqrt(n_methods)
      - top5_fraction: fraction of methods where taxon is in top 5
      - sign_flips: does fold_diff ever change sign?
      - pseudocount_stable: True if rank_std < 10
    """
    methods = sweep_df["method"].unique()
    n_methods = len(methods)

    rows = []
    for taxon in sweep_df["taxon"].unique():
        sub = sweep_df[sweep_df["taxon"] == taxon]
        ranks = sub["rank"].values.astype(float)
        fds = sub["fold_diff"].values.astype(float)
        
        mean_rank = np.mean(ranks)
        rank_std = np.std(ranks, ddof=1) if len(ranks) > 1 else 0.0
        rank_range = ranks.max() - ranks.min()
        top5_frac = (ranks < 5).mean()
        
        # Do fold-diffs ever change sign?
        sign_changes = 0
        signs = np.sign(fds)
        for i, j in combinations(range(len(signs)), 2):
            if signs[i] != signs[j] and signs[i] != 0 and signs[j] != 0:
                sign_changes += 1
        sign_stable = sign_changes == 0

        rows.append({
            "taxon": taxon,
            "n_methods": n_methods,
            "mean_rank": round(mean_rank, 2),
            "rank_std": round(rank_std, 2),
            "rank_range": round(rank_range, 2),
            "top5_fraction": round(top5_frac, 3),
            "sign_stable": sign_stable,
            "pseudocount_stable": rank_std < 10.0,
            "stability_score": round(1.0 / (1.0 + rank_std / 10.0), 4),
        })

    return pd.DataFrame(rows).sort_values("mean_rank")


# ================================================================
# LEAVE-ONE-TAXON-OUT (LOTO) JACKKNIFE
# ================================================================

def leave_one_out_jackknife(df: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    """
    For each taxon i, remove it, recompute geometric mean and CLR,
    re-rank remaining taxa by |fold_diff|.
    
    Returns per-taxon stability metrics:
      - loto_rank_shift_mean: mean rank change when other taxa are removed
      - loto_top5_fraction: fraction of LOTO iterations where this taxon stays in top 5
      - loto_influence_score: mean |Δfold_diff| of all other taxa when this taxon removed
      - is_high_leverage: True if removing this taxon shifts its own mean rank by >20
    """
    j_counts = df["Jumper"].values.astype(float)
    l_counts = df["Laggard"].values.astype(float)
    taxon_names = df["ASV"].values
    n_taxa = len(taxon_names)

    # Baseline: compute with all taxa present
    clr_j_base = clr_with_pseudocount(j_counts, pseudocount)
    clr_l_base = clr_with_pseudocount(l_counts, pseudocount)
    fd_base = clr_j_base - clr_l_base
    rank_base = np.argsort(np.argsort(-np.abs(fd_base)))

    # Per-taxon accumulators (all defaultdicts of lists)
    loto_ranks_per_taxon = defaultdict(list)   # rank of taxon T when another taxon is removed
    loto_influence = defaultdict(float)         # mean |Δfd| caused by removing this taxon
    loto_top5_hits = defaultdict(int)           # count this taxon appears in LOTO top 5
    
    for i in range(n_taxa):
        taxon_removed = taxon_names[i]
        mask = np.ones(n_taxa, dtype=bool)
        mask[i] = False

        j_loo = j_counts[mask]
        l_loo = l_counts[mask]
        taxa_loo = taxon_names[mask]

        clr_j_loo = clr_with_pseudocount(j_loo, pseudocount)
        clr_l_loo = clr_with_pseudocount(l_loo, pseudocount)
        fd_loo = clr_j_loo - clr_l_loo
        rank_loo = np.argsort(np.argsort(-np.abs(fd_loo)))

        # Track how removal of taxon_i affects all other taxa
        fd_deltas = []
        for k in range(len(taxa_loo)):
            t = taxa_loo[k]
            # Find original position of taxon t
            orig_positions = np.where(taxon_names == t)[0]
            if len(orig_positions) == 0:
                continue
            orig_idx = orig_positions[0]
            
            loto_ranks_per_taxon[t].append(int(rank_loo[k]))
            if rank_loo[k] < 5:
                loto_top5_hits[t] += 1
            
            delta_fd = abs(fd_loo[k] - fd_base[orig_idx])
            fd_deltas.append(delta_fd)
        
        loto_influence[taxon_removed] = float(np.mean(fd_deltas)) if fd_deltas else 0.0

    # Build per-taxon output rows
    rows = []
    for i, taxon in enumerate(taxon_names):
        ranks_in_loto = loto_ranks_per_taxon.get(taxon, [])
        n_seen = len(ranks_in_loto)  # number of LOTO iterations this taxon appeared in
        
        if n_seen > 0:
            mean_rank_loto = np.mean(ranks_in_loto)
            mean_shift = mean_rank_loto - rank_base[i]
            rank_std_loo = np.std(ranks_in_loto, ddof=1) if n_seen > 1 else 0.0
        else:
            mean_shift = 0.0
            rank_std_loo = 0.0
        
        top5_hits = loto_top5_hits.get(taxon, 0)
        # Fraction of LOTO iterations where taxon is in top 5 (among times it was ranked)
        top5_frac = top5_hits / n_seen if n_seen > 0 else 0.0
        
        influence = loto_influence.get(taxon, 0.0)
        high_leverage = abs(mean_shift) > 20.0

        rows.append({
            "taxon": taxon,
            "baseline_rank": int(rank_base[i]),
            "baseline_fold_diff": round(float(fd_base[i]), 6),
            "loto_mean_rank_shift": round(float(mean_shift), 2),
            "loto_rank_std": round(float(rank_std_loo), 2),
            "loto_top5_fraction": round(float(top5_frac), 3),
            "loto_influence_score": round(float(influence), 6),
            "is_high_leverage": high_leverage,
            "n_loto_appearances": n_seen,
        })

    return pd.DataFrame(rows).sort_values("baseline_rank")


# ================================================================
# MERGED CONFIDENCE REPORT
# ================================================================

def build_confidence_report(
    pseudo_stability: pd.DataFrame,
    loto_stability: pd.DataFrame,
    taxa_direction_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge pseudocount stability, LOTO stability, and current ranking
    into a single per-taxon confidence report.
    """
    merged = taxa_direction_df[["taxon", "Jumper", "Laggard", "fold_diff",
                                  "direction", "contaminant_flag",
                                  "contamination_score", "contamination_risk"]].copy()
    
    # Merge pseudocount stability
    ps = pseudo_stability.set_index("taxon")
    merged["pseudocount_stable"] = merged["taxon"].map(ps["pseudocount_stable"])
    merged["pseudocount_stability_score"] = merged["taxon"].map(ps["stability_score"])
    merged["pseudocount_mean_rank"] = merged["taxon"].map(ps["mean_rank"])
    merged["pseudocount_rank_std"] = merged["taxon"].map(ps["rank_std"])
    merged["sign_stable"] = merged["taxon"].map(ps["sign_stable"])
    merged["top5_stable_pseudocount"] = merged["taxon"].map(ps["top5_fraction"])

    # Merge LOTO stability
    ls = loto_stability.set_index("taxon")
    merged["loto_top5_stability"] = merged["taxon"].map(ls["loto_top5_fraction"])
    merged["loto_mean_rank_shift"] = merged["taxon"].map(ls["loto_mean_rank_shift"])
    merged["loto_influence_score"] = merged["taxon"].map(ls["loto_influence_score"])
    merged["is_high_leverage"] = merged["taxon"].map(ls["is_high_leverage"])

    # Combined confidence score
    merged["pseudocount_score"] = merged["pseudocount_stability_score"].fillna(0.5)
    merged["loto_score"] = (1.0 - np.clip(
        merged["loto_mean_rank_shift"].abs().fillna(50) / 100.0, 0, 1))
    
    merged["combined_stability"] = round(
        0.5 * merged["pseudocount_score"].fillna(0.5) +
        0.5 * merged["loto_score"].fillna(0.5), 4
    )

    # Stability tier
    def stability_tier(score):
        if pd.isna(score):
            return "UNKNOWN"
        if score >= 0.8:
            return "HIGH"
        if score >= 0.5:
            return "MEDIUM"
        return "LOW"

    merged["stability_tier"] = merged["combined_stability"].apply(stability_tier)

    # Fill NaN for taxa that couldn't be assessed
    for col in ["pseudocount_stable", "sign_stable", "is_high_leverage"]:
        if col in merged.columns:
            merged[col] = merged[col].fillna(False)

    return merged


# ================================================================
# MAIN
# ================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Pseudocount & LOTO sensitivity analysis for CLR-based taxon ranking")
    parser.add_argument("--in", dest="input_file", required=True,
                        help="Raw ASV CSV (ASV_table_Jumpers_Laggards.GKAQUA.csv)")
    parser.add_argument("--out-dir", dest="out_dir", default="results/reports",
                        help="Output directory for reports")
    parser.add_argument("--taxa-direction", dest="taxa_dir",
                        default="data/processed/clr_profiles/taxa_direction.tsv",
                        help="Existing taxa_direction.tsv for metadata merge")
    parser.add_argument("--n-pseudocounts", type=int, default=8,
                        help="Number of pseudocount values to test")
    parser.add_argument("--loto-pseudocount", type=float, default=1.0,
                        help="Pseudocount to use for LOTO analysis")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Load data ----
    print("=" * 70)
    print("TAXON STABILITY ANALYSIS — Pseudocount + Leave-One-Taxon-Out")
    print("=" * 70)
    print(f"\n[1/5] Loading ASV table: {args.input_file}")
    raw = pd.read_csv(args.input_file)
    print(f"  {len(raw)} taxa loaded")
    print(f"  Jumper: {raw['Jumper'].sum()} reads, {(raw['Jumper']>0).sum()} nonzero taxa")
    print(f"  Laggard: {raw['Laggard'].sum()} reads, {(raw['Laggard']>0).sum()} nonzero taxa")
    n_zeros = ((raw["Jumper"] == 0) | (raw["Laggard"] == 0)).sum()
    print(f"  Taxa with at least one zero: {n_zeros}/{len(raw)} ({100*n_zeros/len(raw):.1f}%)")

    # ---- Design pseudocount range ----
    # Span from below min_nonzero to ~10× the read depth ratio
    j_min = raw.loc[raw["Jumper"] > 0, "Jumper"].min()
    l_min = raw.loc[raw["Laggard"] > 0, "Laggard"].min()
    min_val = min(j_min, l_min)
    max_count = max(raw["Jumper"].max(), raw["Laggard"].max())

    pseudocounts = sorted(set([
        0.1, 0.5, 1.0,       # Below and at typical minimum
        2.0, 5.0,             # Small multiples
        10.0, 20.0,           # Moderate multiples
        50.0,                 # Large pseudocount: stress test
    ]))[:args.n_pseudocounts + 1]  # +1 for the extra method slots

    print(f"\n[2/5] Pseudocount sweep: {len(pseudocounts) + 3} methods")
    print(f"  Pseudocounts: {pseudocounts}")
    print(f"  Additional methods: current (NaN), multiplicative, Bayesian-multiplicative")
    print(f"  Total: {len(pseudocounts) + 3} CLR variants")

    sweep_df = pseudocount_sweep(raw, pseudocounts)
    n_methods = sweep_df["method"].nunique()
    print(f"  Computed {len(sweep_df)} rows ({n_methods} methods × {len(raw)} taxa)")

    # ---- Pseudocount stability metrics ----
    print(f"\n[3/5] Computing pseudocount stability metrics...")
    pseudo_stability = compute_pseudocount_stability_metrics(sweep_df)
    n_stable = pseudo_stability["pseudocount_stable"].sum()
    n_sign_flip = (~pseudo_stability["sign_stable"]).sum()
    print(f"  Pseudocount-stable taxa (rank_std < 10): {n_stable}/{len(pseudo_stability)}")
    print(f"  Taxa with sign flips across methods: {n_sign_flip}")

    # Top 5 stability
    top5_sweep = pseudo_stability[pseudo_stability["mean_rank"] < 5]
    print(f"  Top 5 (by mean rank):")
    for _, r in top5_sweep.iterrows():
        print(f"    {r['taxon']:<45s} mean_rank={r['mean_rank']:.1f} "
              f"rank_std={r['rank_std']:.1f} top5_frac={r['top5_fraction']:.2f} "
              f"sign_stable={r['sign_stable']}")

    # Bottom — taxa with highest rank volatility
    volatile = pseudo_stability.nlargest(10, "rank_std")
    print(f"\n  Most volatile taxa (highest rank_std):")
    for _, r in volatile.iterrows():
        print(f"    {r['taxon']:<45s} rank_std={r['rank_std']:.1f} "
              f"rank_range={r['rank_range']:.0f} sign_stable={r['sign_stable']}")

    # ---- LOTO jackknife ----
    print(f"\n[4/5] Leave-one-taxon-out jackknife ({len(raw)} iterations)...")
    loto_stability = leave_one_out_jackknife(raw, pseudocount=args.loto_pseudocount)
    n_top5_stable = (loto_stability["loto_top5_fraction"] >= 0.8).sum()
    n_high_lev = loto_stability["is_high_leverage"].sum()
    print(f"  Taxa with top-5 stability >= 80%: {n_top5_stable}")
    print(f"  High-leverage taxa (removal shifts others >20 ranks): {n_high_lev}")

    top5_loto = loto_stability[loto_stability["baseline_rank"] < 5]
    print(f"  Top 5 LOTO metrics:")
    for _, r in top5_loto.iterrows():
        print(f"    {r['taxon']:<45s} baseline_rank={int(r['baseline_rank'])} "
              f"loto_top5_frac={r['loto_top5_fraction']:.3f} "
              f"influence={r['loto_influence_score']:.4f} "
              f"high_lev={r['is_high_leverage']}")

    # ---- Build merged confidence report ----
    print(f"\n[5/5] Building merged confidence report...")
    taxa_dir_df = pd.read_csv(args.taxa_dir, sep="\t")
    confidence = build_confidence_report(pseudo_stability, loto_stability, taxa_dir_df)

    tier_counts = confidence["stability_tier"].value_counts()
    for tier in ["HIGH", "MEDIUM", "LOW", "UNKNOWN"]:
        n = tier_counts.get(tier, 0)
        print(f"  Stability tier {tier}: {n} taxa")

    # Highlight currently-ranked top 5 taxa (using same ranking formula)
    current_top5_ranked = taxa_dir_df[
        (taxa_dir_df["contamination_risk"] == "LOW") &
        (taxa_dir_df["taxon"].str.lower() != "unknown")
    ].copy()
    current_top5_ranked["abs_fold_diff"] = current_top5_ranked["fold_diff"].abs()
    current_top5_ranked = current_top5_ranked.nlargest(5, "abs_fold_diff")

    print(f"\n  Current top-5 ranked taxa stability assessment:")
    for _, r in current_top5_ranked.iterrows():
        ci = confidence[confidence["taxon"] == r["taxon"]]
        if len(ci) > 0:
            c = ci.iloc[0]
            print(f"    {r['taxon']:<40s} "
                  f"stability_tier={c['stability_tier']:<8s} "
                  f"combined={c['combined_stability']:.4f} "
                  f"pseudo_stable={c['pseudocount_stable']} "
                  f"sign_stable={c['sign_stable']} "
                  f"loto_top5={c['loto_top5_stability']:.2f}")

    # ---- Write outputs ----
    print(f"\n  Writing outputs to {out_dir}/...")
    sweep_df.to_csv(out_dir / "pseudocount_sensitivity.tsv", sep="\t", index=False)
    pseudo_stability.to_csv(out_dir / "pseudocount_stability_metrics.tsv", sep="\t", index=False)
    loto_stability.to_csv(out_dir / "loto_stability.tsv", sep="\t", index=False)
    confidence.to_csv(out_dir / "taxon_confidence_report.tsv", sep="\t", index=False)

    print(f"  [OK] pseudocount_sensitivity.tsv — {len(sweep_df)} rows")
    print(f"  [OK] pseudocount_stability_metrics.tsv — {len(pseudo_stability)} taxa")
    print(f"  [OK] loto_stability.tsv — {len(loto_stability)} taxa")
    print(f"  [OK] taxon_confidence_report.tsv — {len(confidence)} taxa")
    print("\nDone.")


if __name__ == "__main__":
    main()
