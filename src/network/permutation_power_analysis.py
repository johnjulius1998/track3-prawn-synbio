#!/usr/bin/env python3
"""
permutation_power_analysis.py — Permutation Null + Power Analysis
===================================================================
Track 3 Host-Microbe Integration (v3.3)

PURPOSE:
  1. Permutation null: Randomly shuffle weight_gain 1000×, recompute
     module eigengene partial correlations. Test whether observed
     partial_r values are unusual under the null.

  2. Power analysis: Given observed effect sizes, compute the sample
     size needed for 80% power at α=0.05.

METHOD:
  Permutation:
    - Shuffle WG labels across the 20 samples (preserving sex/tissue pairing)
    - Recompute partial_r(WG|sex) for each module eigengene
    - Build null distribution: 95th/99th percentile of |partial_r|
    - Flag modules exceeding the null threshold as "non-null"

  Power:
    - Use observed partial_r as effect size estimate
    - Compute required n for 80% power with sex as covariate
    - Also compute: "what n would we need for the bootstrap CI to exclude zero?"

OUTPUTS:
  results/reports/
    permutation_null.tsv    — Per-module null comparison
    power_analysis.tsv      — Per-module power estimates

USAGE:
  python src/network/permutation_power_analysis.py \
      --me data/processed/wgcna/me_matrix.tsv \
      --metadata data/raw/sra/PRJNA875278/metadata.tsv \
      --out-dir results/reports/ \
      --n-perms 1000
"""

import argparse, sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")


def partial_corr(x, y, z):
    """Partial correlation r(x,y | z)."""
    r_xy = np.corrcoef(x, y)[0, 1]
    r_xz = np.corrcoef(x, z)[0, 1]
    r_yz = np.corrcoef(y, z)[0, 1]
    denom = np.sqrt((1 - r_xz**2) * (1 - r_yz**2))
    if denom == 0 or np.isnan(denom):
        return 0.0
    return np.clip((r_xy - r_xz * r_yz) / denom, -1.0, 1.0)


def run_permutation_null(me_df, meta_df, n_perms=1000, seed=42):
    """
    Permutation test: shuffle WG, recompute partial_r for each module.
    Returns observed stats + null distribution thresholds.
    """
    print(f"\n{'='*60}")
    print(f"PERMUTATION NULL TEST ({n_perms} permutations)")
    print(f"{'='*60}")
    
    rng = np.random.default_rng(seed)
    common = sorted(set(me_df.index) & set(meta_df.index))
    me = me_df.loc[common]
    meta = meta_df.loc[common]
    
    wg = meta["weight_gain"].values.astype(float)
    sex = np.array([1.0 if str(meta.loc[s, "sex"]).lower() == "male" else 0.0
                     for s in common])
    
    module_cols = [c for c in me.columns if c.startswith("ME")]
    n_modules = len(module_cols)
    
    # Observed
    obs_stats = {}
    for col in module_cols:
        mv = me[col].values.astype(float)
        obs_stats[col] = partial_corr(mv, wg, sex)
    
    # Permutation null
    null_dist = {col: [] for col in module_cols}
    
    for i in range(n_perms):
        wg_perm = rng.permutation(wg)
        for col in module_cols:
            mv = me[col].values.astype(float)
            rp_perm = partial_corr(mv, wg_perm, sex)
            null_dist[col].append(rp_perm)
    
    # Build results
    rows = []
    for col in module_cols:
        mod = col.replace("ME", "M")
        obs = obs_stats[col]
        nulls = np.array(null_dist[col])
        
        # Null thresholds
        abs_nulls = np.abs(nulls)
        p95 = np.percentile(abs_nulls, 95)
        p99 = np.percentile(abs_nulls, 99)
        p999 = np.percentile(abs_nulls, 99.9)
        
        # Empirical p-value
        emp_p = (np.abs(nulls) >= np.abs(obs)).mean()
        
        # Is observed beyond null?
        sig_95 = np.abs(obs) > p95
        sig_99 = np.abs(obs) > p99
        
        rows.append({
            "module": mod,
            "partial_r_obs": round(obs, 4),
            "null_95th_pct": round(p95, 4),
            "null_99th_pct": round(p99, 4),
            "null_99_9th_pct": round(p999, 4),
            "null_mean": round(np.mean(nulls), 4),
            "null_std": round(np.std(nulls, ddof=1), 4),
            "empirical_p": round(emp_p, 4),
            "exceeds_95th": sig_95,
            "exceeds_99th": sig_99,
            "n_perms": n_perms,
        })
    
    result = pd.DataFrame(rows).sort_values("partial_r_obs", key=abs, ascending=False)
    
    # Summary
    n_sig_95 = result["exceeds_95th"].sum()
    n_sig_99 = result["exceeds_99th"].sum()
    print(f"\n  Modules exceeding null 95th percentile: {n_sig_95}/{n_modules}")
    print(f"  Modules exceeding null 99th percentile: {n_sig_99}/{n_modules}")
    
    print(f"\n  {'Module':<6s} {'partial_r':>9s} {'95th_null':>9s} {'99th_null':>9s} {'emp_p':>7s} {'Sig?':>5s}")
    print(f"  {'-'*6} {'-'*9} {'-'*9} {'-'*9} {'-'*7} {'-'*5}")
    for _, r in result.iterrows():
        marker = "**" if r["exceeds_99th"] else ("*" if r["exceeds_95th"] else "")
        print(f"  {r['module']:<6s} {r['partial_r_obs']:>+9.4f} "
              f"{r['null_95th_pct']:>9.4f} {r['null_99th_pct']:>9.4f} "
              f"{r['empirical_p']:>7.4f} {marker:>5s}")
    
    print(f"\n  *  = exceeds permutation 95th percentile")
    print(f"  ** = exceeds permutation 99th percentile")
    
    return result


def run_power_analysis(perm_result, n_current=20, alpha=0.05, target_power=0.80):
    """
    Power analysis: given observed partial_r, what n is needed?
    
    For partial correlation with 1 confounder (sex), df = n - 3.
    Fisher z-transformation used for sample size calculation.
    """
    print(f"\n{'='*60}")
    print(f"POWER ANALYSIS (α={alpha}, target power={target_power})")
    print(f"{'='*60}")
    
    rows = []
    for _, r in perm_result.iterrows():
        obs_r = r["partial_r_obs"]
        mod = r["module"]
        
        if abs(obs_r) < 0.01:
            # Near-zero effect — would need infinite n
            n_needed = float("inf")
            power_at_current = 0.05  # roughly alpha
        else:
            # Fisher z-transform
            z_obs = np.arctanh(abs(obs_r))
            
            # For partial correlation with k=1 confounder, effective df = n - k - 2 = n - 3
            # Standard error of z: 1/sqrt(n - 3 - 1) = 1/sqrt(n - 4) for simple
            # More precisely: SE(z) = 1/sqrt(n - 3) for partial correlation
            
            # Required n: n = ((z_α + z_β) / z_obs)^2 + 3
            z_alpha = stats.norm.ppf(1 - alpha / 2)  # two-sided
            z_beta = stats.norm.ppf(target_power)
            
            n_needed_raw = ((z_alpha + z_beta) / z_obs) ** 2 + 3
            n_needed = int(np.ceil(n_needed_raw))
            
            # Power at current n
            se_current = 1.0 / np.sqrt(n_current - 3) if n_current > 3 else 1.0
            z_current = z_obs / se_current
            power_current = stats.norm.cdf(z_current - z_alpha) + stats.norm.cdf(-z_current - z_alpha)
            # Simplified: power = P(|Z| > z_alpha) where Z ~ N(z_obs * sqrt(n-3), 1)
            power_at_current = 1 - stats.norm.cdf(z_alpha - z_obs * np.sqrt(n_current - 3))
        
        # Bootstrap CI: when would CI exclude zero?
        # CI excludes zero when |r| > z_α / sqrt(n - 3)
        # So n_needed_for_ci = (z_α / r)^2 + 3
        if abs(obs_r) > 0.01:
            n_for_ci = int(np.ceil((z_alpha / abs(obs_r))**2 + 3))
        else:
            n_for_ci = float("inf")
        
        rows.append({
            "module": mod,
            "partial_r_obs": obs_r,
            "n_current": n_current,
            "power_at_n20": round(min(power_at_current, 1.0), 3),
            "n_for_80pct_power": n_needed,
            "n_for_ci_exclude_zero": n_for_ci,
            "alpha": alpha,
            "target_power": target_power,
        })
    
    result = pd.DataFrame(rows).sort_values("partial_r_obs", key=abs, ascending=False)
    
    print(f"\n  {'Module':<6s} {'|r|':>7s} {'Power@n=20':>10s} {'n for 80%':>9s} {'n for CI≠0':>10s}")
    print(f"  {'-'*6} {'-'*7} {'-'*10} {'-'*9} {'-'*10}")
    for _, r in result.iterrows():
        n80_str = f"{r['n_for_80pct_power']}" if r['n_for_80pct_power'] != float('inf') else "∞"
        nci_str = f"{r['n_for_ci_exclude_zero']}" if r['n_for_ci_exclude_zero'] != float('inf') else "∞"
        print(f"  {r['module']:<6s} {abs(r['partial_r_obs']):>7.4f} "
              f"{r['power_at_n20']:>10.3f} {n80_str:>9s} {nci_str:>10s}")
    
    # Key recommendation
    top_mod = result.iloc[0]
    print(f"\n  Key recommendation:")
    print(f"    To confirm {top_mod['module']} (r={top_mod['partial_r_obs']:+.3f})"
          f" with 80% power: n ≥ {top_mod['n_for_80pct_power']} samples")
    
    # Also report: how many modules are adequately powered at n=20?
    powered = (result["power_at_n20"] >= 0.80).sum()
    print(f"    Modules with ≥80% power at n=20: {powered}/{len(result)}")
    if powered == 0:
        print(f"    → No module is adequately powered at current sample size.")
        print(f"    → This is the core statistical limitation of this study.")
    
    return result


def main():
    parser = argparse.ArgumentParser(description="Permutation null + power analysis")
    parser.add_argument("--me", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out-dir", default="results/reports")
    parser.add_argument("--n-perms", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("PERMUTATION NULL + POWER ANALYSIS (v3.3)")
    print("=" * 60)
    
    # Load
    print("\n[1/3] Loading data...")
    me = pd.read_csv(args.me, sep="\t", index_col=0)
    meta = pd.read_csv(args.metadata, sep="\t", index_col=0)
    print(f"  ME matrix: {me.shape[0]} samples × {me.shape[1]} columns")
    
    # Permutation null
    print(f"\n[2/3] Running {args.n_perms} permutations...")
    perm_result = run_permutation_null(me, meta, n_perms=args.n_perms, seed=args.seed)
    
    # Power analysis
    print(f"\n[3/3] Computing power...")
    power_result = run_power_analysis(perm_result)
    
    # Write outputs
    perm_result.to_csv(out_dir / "permutation_null.tsv", sep="\t", index=False)
    power_result.to_csv(out_dir / "power_analysis.tsv", sep="\t", index=False)
    print(f"\n  [OK] {out_dir}/permutation_null.tsv")
    print(f"  [OK] {out_dir}/power_analysis.tsv")
    print("\nDone.")


if __name__ == "__main__":
    main()
