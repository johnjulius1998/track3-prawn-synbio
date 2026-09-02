#!/usr/bin/env python3
"""
wgcna_bootstrap_stability.py — Bootstrap Module-Trait Stability
=================================================================
Track 3 Host-Microbe Integration (v3 — Internal Validation)

PURPOSE:
  Quantify the stability of module-trait associations via bootstrap
  resampling of the 20 samples. Does NOT re-run WGCNA — uses the
  existing module eigengenes and resamples the sample labels to
  compute confidence intervals on partial correlations.

METHOD:
  1. Load existing ME matrix and metadata (n=20)
  2. Bootstrap resample (with replacement) n=20 samples, N=1000 iterations
  3. For each bootstrap, recompute partial_r(ME, weight_gain | sex)
  4. Report: 95% CI, fraction of bootstraps where module is in top 5,
     and bootstrap p-value (fraction where partial_r sign flips)

OUTPUTS:
  results/reports/
    wgcna_bootstrap_stability.tsv  — Per-module bootstrap statistics
    wgcna_bootstrap_top5.tsv       — Top-5 co-occurrence matrix

USAGE:
  python src/network/wgcna_bootstrap_stability.py \
      --me data/processed/wgcna/me_matrix.tsv \
      --metadata data/raw/sra/PRJNA875278/metadata.tsv \
      --out-dir results/reports/ \
      --n-bootstrap 1000
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from collections import defaultdict


def partial_corr(x, y, z):
    """
    Compute partial correlation r(x,y | z).
    r_xy_z = (r_xy - r_xz * r_yz) / sqrt((1 - r_xz^2) * (1 - r_yz^2))
    """
    r_xy = np.corrcoef(x, y)[0, 1]
    r_xz = np.corrcoef(x, z)[0, 1]
    r_yz = np.corrcoef(y, z)[0, 1]
    
    denom = np.sqrt((1 - r_xz**2) * (1 - r_yz**2))
    if denom == 0 or np.isnan(denom):
        return 0.0
    result = (r_xy - r_xz * r_yz) / denom
    return np.clip(result, -1.0, 1.0)


def bootstrap_module_stability(
    me_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    n_bootstrap: int = 1000,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Bootstrap module-trait partial correlations.
    
    Returns per-module statistics:
      - partial_r_obs: observed partial r (full dataset)
      - partial_r_mean: bootstrap mean
      - partial_r_ci_lower / ci_upper: 95% percentile CI
      - top5_fraction: fraction of bootstraps module is in top 5 by |partial_r|
      - sign_flip_rate: fraction of bootstraps where sign(partial_r) ≠ sign(observed)
      - bootstrap_p_value: fraction where |partial_r_boot| >= |partial_r_obs|
    """
    rng = np.random.default_rng(seed)
    
    # Align samples
    common = sorted(set(me_df.index) & set(meta_df.index))
    me = me_df.loc[common]
    meta = meta_df.loc[common]
    
    wg = meta["weight_gain"].values.astype(float)
    sex = np.array([1.0 if str(meta.loc[s, "sex"]).lower() == "male" else 0.0
                     for s in common])
    
    n_samples = len(common)
    module_cols = [c for c in me.columns if c.startswith("ME")]
    
    # Observed partial correlations
    obs_stats = {}
    for col in module_cols:
        mv = me[col].values.astype(float)
        rp = partial_corr(mv, wg, sex)
        obs_stats[col] = rp
    
    # Bootstrap
    bootstrap_results = defaultdict(list)
    n_in_top5 = defaultdict(int)
    n_sign_flip = defaultdict(int)
    
    for b in range(n_bootstrap):
        idx = rng.choice(n_samples, size=n_samples, replace=True)
        wg_boot = wg[idx]
        sex_boot = sex[idx]
        
        boot_partial_rs = {}
        for col in module_cols:
            mv_boot = me[col].values.astype(float)[idx]
            rp_boot = partial_corr(mv_boot, wg_boot, sex_boot)
            bootstrap_results[col].append(rp_boot)
            boot_partial_rs[col] = rp_boot
            
            # Sign flip?
            if np.sign(rp_boot) != np.sign(obs_stats[col]) and abs(rp_boot) > 0.01:
                n_sign_flip[col] += 1
        
        # Which modules are in top 5 by |partial_r|?
        top5_boot = sorted(boot_partial_rs.keys(), 
                          key=lambda k: abs(boot_partial_rs[k]), 
                          reverse=True)[:5]
        for col in top5_boot:
            n_in_top5[col] += 1
    
    # Build output
    rows = []
    for col in module_cols:
        rp_obs = obs_stats[col]
        boots = np.array(bootstrap_results[col])
        
        ci_lower = np.percentile(boots, 2.5)
        ci_upper = np.percentile(boots, 97.5)
        rp_mean = np.mean(boots)
        rp_std = np.std(boots, ddof=1)
        
        top5_frac = n_in_top5.get(col, 0) / n_bootstrap
        sign_flip = n_sign_flip.get(col, 0) / n_bootstrap
        bootstrap_p = (np.abs(boots) >= np.abs(rp_obs)).mean()
        
        # Stability tier
        if top5_frac >= 0.8 and sign_flip < 0.05:
            tier = "HIGH"
        elif top5_frac >= 0.5 and sign_flip < 0.2:
            tier = "MEDIUM"
        else:
            tier = "LOW"
        
        rows.append({
            "module": col.replace("ME", "M"),
            "n_samples": n_samples,
            "partial_r_obs": round(rp_obs, 4),
            "partial_r_bootstrap_mean": round(rp_mean, 4),
            "partial_r_bootstrap_std": round(rp_std, 4),
            "partial_r_ci95_lower": round(ci_lower, 4),
            "partial_r_ci95_upper": round(ci_upper, 4),
            "ci95_excludes_zero": (ci_lower > 0) or (ci_upper < 0),
            "top5_fraction": round(top5_frac, 3),
            "sign_flip_rate": round(sign_flip, 3),
            "bootstrap_p_value": round(bootstrap_p, 3),
            "stability_tier": tier,
            "n_bootstraps": n_bootstrap,
        })
    
    return pd.DataFrame(rows).sort_values("top5_fraction", ascending=False)


def build_top5_heatmap_data(boot_results, n_bootstrap=1000):
    """
    Build top-5 co-occurrence matrix: for each pair of modules,
    what fraction of bootstraps do they co-occur in the top 5?
    """
    # This requires the raw bootstrap data, not the summary
    # We'll return a simpler pairwise co-rank matrix from the summary data
    pass


def main():
    parser = argparse.ArgumentParser(
        description="Bootstrap stability of WGCNA module-trait associations")
    parser.add_argument("--me", required=True,
                        help="Module eigengene matrix (me_matrix.tsv)")
    parser.add_argument("--metadata", required=True,
                        help="Sample metadata with weight_gain and sex")
    parser.add_argument("--out-dir", default="results/reports",
                        help="Output directory")
    parser.add_argument("--n-bootstrap", type=int, default=1000,
                        help="Number of bootstrap iterations")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility")
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("WGCNA MODULE-TRAIT BOOTSTRAP STABILITY")
    print("=" * 70)
    
    # Load data
    print(f"\n[1/3] Loading data...")
    me = pd.read_csv(args.me, sep="\t", index_col=0)
    meta = pd.read_csv(args.metadata, sep="\t", index_col=0)
    n_modules = len([c for c in me.columns if c.startswith("ME")])
    print(f"  ME matrix: {me.shape[0]} samples × {n_modules} modules")
    print(f"  Metadata: {meta.shape[0]} samples")
    print(f"  Weight gain range: {meta['weight_gain'].min():.2f}–{meta['weight_gain'].max():.2f}")
    
    # Bootstrap
    print(f"\n[2/3] Running {args.n_bootstrap} bootstrap iterations...")
    results = bootstrap_module_stability(
        me, meta, n_bootstrap=args.n_bootstrap, seed=args.seed)
    
    # Summary
    print(f"\n[3/3] Results:")
    high = results[results["stability_tier"] == "HIGH"]
    med = results[results["stability_tier"] == "MEDIUM"]
    low = results[results["stability_tier"] == "LOW"]
    print(f"  HIGH stability: {len(high)} modules")
    print(f"  MEDIUM stability: {len(med)} modules")
    print(f"  LOW stability: {len(low)} modules")
    
    print(f"\n  Module stability ranking (by top5_fraction):")
    for _, r in results.head(10).iterrows():
        ci_str = f"[{r['partial_r_ci95_lower']:+.3f}, {r['partial_r_ci95_upper']:+.3f}]"
        marker = "✓" if r['ci95_excludes_zero'] else " "
        print(f"    {r['module']:<6s} partial_r={r['partial_r_obs']:+.4f} "
              f"95%CI={ci_str} {marker} "
              f"top5_rate={r['top5_fraction']:.3f} "
              f"sign_flip={r['sign_flip_rate']:.3f} "
              f"tier={r['stability_tier']}")
    
    # Check current top modules (M6, M8, M13, M7, M15 from the report)
    current_top = ["M6", "M8", "M13", "M7", "M15"]
    print(f"\n  Current top-5 growth modules stability:")
    for m in current_top:
        row = results[results["module"] == m]
        if len(row) > 0:
            r = row.iloc[0]
            ci_str = f"[{r['partial_r_ci95_lower']:+.3f}, {r['partial_r_ci95_upper']:+.3f}]"
            print(f"    {m}: partial_r={r['partial_r_obs']:+.4f} 95%CI={ci_str} "
                  f"top5_rate={r['top5_fraction']:.3f} "
                  f"sign_flip={r['sign_flip_rate']:.3f} "
                  f"tier={r['stability_tier']}")
    
    # Write output
    out_path = out_dir / "wgcna_bootstrap_stability.tsv"
    results.to_csv(out_path, sep="\t", index=False)
    print(f"\n  [OK] {out_path} — {len(results)} modules")
    
    # Also save a compact version for the host gene ranking
    compact = results[["module", "partial_r_obs", "top5_fraction", 
                        "sign_flip_rate", "stability_tier"]].copy()
    compact.columns = ["module", "partial_r_wg_given_sex", 
                        "module_top5_stability", "module_sign_flip_rate",
                        "module_stability_tier"]
    compact_path = out_dir / "wgcna_module_stability_summary.tsv"
    compact.to_csv(compact_path, sep="\t", index=False)
    print(f"  [OK] {compact_path} — compact summary for downstream ranking")
    
    print("\nDone.")


if __name__ == "__main__":
    main()
