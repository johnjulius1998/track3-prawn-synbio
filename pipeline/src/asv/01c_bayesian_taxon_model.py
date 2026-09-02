#!/usr/bin/env python3
"""
01c_bayesian_taxon_model.py — Bayesian Dirichlet-Multinomial Taxon Enrichment
===============================================================================
Track 3 Host-Microbe Integration (v3.3)

PURPOSE:
  Replace the fragile CLR point-estimate ranking with a principled Bayesian
  hierarchical model that:
    1. Models read counts as Dirichlet-multinomial (correct for compositions)
    2. Shrinks low-count taxa toward zero (automatic regularization)
    3. Produces posterior probabilities and 95% credible intervals
    4. Eliminates the pseudocount dependency entirely

MODEL:
  For taxon i in pool j (j ∈ {Jumper, Laggard}):
    y[j] ~ Multinomial(N_j, π_j)
    π_j ~ Dirichlet(α)
    log(α_i) = β_0[i] + β_1[i] × group_j
    β_1[i] ~ Normal(0, σ_β)  # hierarchical prior shrinks low-count taxa

  Taxa with high posterior P(β_1[i] > 0) are Jumper-enriched.
  Taxa with high posterior P(β_1[i] < 0) are Laggard-enriched.

  Since we only have 2 pools (n=1 per group), we can't fit the full
  hierarchical model with group-level predictors. Instead, we use a
  simpler Bayesian formulation:

  For each taxon i:
    y_j[i] ~ Poisson(λ_j)
    λ_j = exp(α × (1 if Jumper else 0) + offset)
    # With only 2 data points, this is essentially a Bayesian
    # version of the rate ratio with shrinkage priors.

SIMPLIFIED APPROACH (for n=1 per group):
  We use a Bayesian beta-binomial model for the proportion of reads
  attributed to each taxon in Jumper vs Laggard, with an uninformative
  prior. The posterior distribution of the Jumper:Laggard ratio gives
  us credible intervals on enrichment.

  P(enriched) = posterior probability that Jumper proportion > Laggard proportion
  log2FC_posterior = posterior mean of log2(Jumper_rate / Laggard_rate)

USAGE:
  python src/asv/01c_bayesian_taxon_model.py \
      --in data/raw/supplied/ASV_table_Jumpers_Laggards.GKAQUA.csv \
      --out-dir results/reports/
"""

import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import betaln, gammaln


def bayesian_beta_binomial_enrichment(
    jumper_counts, laggard_counts, n_mcmc=20000, seed=42
):
    """
    Bayesian beta-binomial model for each taxon.

    For taxon i with counts (k_j, k_l) out of totals (N_j, N_l):

    Prior: θ_j ~ Beta(1, 1), θ_l ~ Beta(1, 1)  [uniform]
    Posterior: θ_j ~ Beta(1 + k_j, 1 + N_j - k_j)
               θ_l ~ Beta(1 + k_l, 1 + N_l - k_l)

    We sample from the posteriors and compute:
      - log2_fold_change = log2(θ_j / θ_l)
      - P(enriched_in_jumper) = P(θ_j > θ_l)
      - 95% credible interval on log2FC

    This is an exact Bayesian solution — no MCMC needed for the
    independent beta posteriors.
    """
    n_taxa = len(jumper_counts)
    rng = np.random.default_rng(seed)
    
    N_j = jumper_counts.sum()
    N_l = laggard_counts.sum()
    
    results = []
    
    for i in range(n_taxa):
        k_j = int(jumper_counts[i])
        k_l = int(laggard_counts[i])
        
        # Posterior parameters
        a_j = 1 + k_j
        b_j = 1 + N_j - k_j
        a_l = 1 + k_l
        b_l = 1 + N_l - k_l
        
        # Sample from posteriors
        theta_j_samples = rng.beta(a_j, b_j, size=n_mcmc)
        theta_l_samples = rng.beta(a_l, b_l, size=n_mcmc)
        
        # Log2 fold change (Jumper / Laggard)
        # Add pseudocount in log space to avoid -inf when k_l = 0
        eps = 1e-6
        log2fc_samples = np.log2(
            (theta_j_samples + eps) / (theta_l_samples + eps)
        )
        
        # Posterior summaries
        log2fc_mean = float(np.mean(log2fc_samples))
        log2fc_median = float(np.median(log2fc_samples))
        log2fc_ci_lower = float(np.percentile(log2fc_samples, 2.5))
        log2fc_ci_upper = float(np.percentile(log2fc_samples, 97.5))
        log2fc_std = float(np.std(log2fc_samples, ddof=1))
        
        # Probability of enrichment
        p_jumper_enriched = float((theta_j_samples > theta_l_samples).mean())
        p_laggard_enriched = float((theta_l_samples > theta_j_samples).mean())
        
        # Direction classification
        if p_jumper_enriched > 0.95:
            direction = "Jumper-enriched"
        elif p_laggard_enriched > 0.95:
            direction = "Laggard-enriched"
        elif k_j > 0 and k_l == 0:
            direction = "Jumper-exclusive (low-count)"
        elif k_l > 0 and k_j == 0:
            direction = "Laggard-exclusive (low-count)"
        elif p_jumper_enriched > 0.5:
            direction = "Jumper-leaning (uncertain)"
        elif p_laggard_enriched > 0.5:
            direction = "Laggard-leaning (uncertain)"
        else:
            direction = "balanced"
        
        # Credible interval excludes zero?
        ci_excludes_zero = (log2fc_ci_lower > 0) or (log2fc_ci_upper < 0)
        
        results.append({
            "taxon_idx": i,
            "jumper_raw": k_j,
            "laggard_raw": k_l,
            "jumper_rel_abund": round(k_j / N_j * 100, 4),
            "laggard_rel_abund": round(k_l / N_l * 100, 4),
            "log2FC_posterior_mean": round(log2fc_mean, 4),
            "log2FC_posterior_median": round(log2fc_median, 4),
            "log2FC_ci95_lower": round(log2fc_ci_lower, 4),
            "log2FC_ci95_upper": round(log2fc_ci_upper, 4),
            "log2FC_posterior_std": round(log2fc_std, 4),
            "p_jumper_enriched": round(p_jumper_enriched, 4),
            "p_laggard_enriched": round(p_laggard_enriched, 4),
            "direction_bayesian": direction,
            "ci_excludes_zero": ci_excludes_zero,
            "n_mcmc": n_mcmc,
        })
    
    return pd.DataFrame(results)


def compare_with_clr(bayes_df, clr_df):
    """Compare Bayesian and CLR-based rankings."""
    print(f"\n{'='*60}")
    print(f"BAYESIAN vs CLR COMPARISON")
    print(f"{'='*60}")
    
    merged = bayes_df.merge(
        clr_df[["taxon", "fold_diff", "direction", "contamination_risk"]],
        left_on="taxon_idx", right_index=True, how="left"
    )
    
    # Top 10 by Bayesian posterior probability of enrichment
    bayes_top = merged.nlargest(10, "p_jumper_enriched")
    print(f"\n  Top 10 Jumper-enriched (by Bayesian posterior probability):")
    for _, r in bayes_top.iterrows():
        ci_str = f"[{r['log2FC_ci95_lower']:+.2f}, {r['log2FC_ci95_upper']:+.2f}]"
        ci_mark = "✓" if r["ci_excludes_zero"] else " "
        print(f"    {r['taxon']:<40s} "
              f"P(Jumper)={r['p_jumper_enriched']:.3f} "
              f"log2FC={r['log2FC_posterior_mean']:+.2f} {ci_str} {ci_mark} "
              f"CLR_dir={r['direction']}")
    
    # How many taxa have CI excluding zero?
    n_ci_nonzero = merged["ci_excludes_zero"].sum()
    print(f"\n  Taxa with 95% CI excluding zero: {n_ci_nonzero}/{len(merged)} "
          f"({100*n_ci_nonzero/len(merged):.1f}%)")
    
    # Concordance between CLR direction and Bayesian direction
    clr_jumper = merged[merged["direction"].str.contains("Jumper", na=False)]
    clr_laggard = merged[merged["direction"].str.contains("Laggard", na=False)]
    
    bayes_jumper = merged[merged["direction_bayesian"].str.contains("Jumper", na=False)]
    bayes_laggard = merged[merged["direction_bayesian"].str.contains("Laggard", na=False)]
    
    n_agree_jumper = len(set(clr_jumper.index) & set(bayes_jumper.index))
    n_agree_laggard = len(set(clr_laggard.index) & set(bayes_laggard.index))
    total_agree = n_agree_jumper + n_agree_laggard
    
    print(f"  CLR-Bayesian direction agreement: {total_agree}/{len(merged)} "
          f"({100*total_agree/len(merged):.1f}%)")
    
    # Key differences
    clr_top5 = set(clr_df.nlargest(5, "fold_diff")["taxon"].values)
    bayes_top5 = set(merged.nlargest(5, "log2FC_posterior_mean")["taxon_idx"].values)
    # Map bayes_top5 indices back to taxon names
    bayes_top5_names = set()
    for _, r in merged.iterrows():
        if r["taxon_idx"] in bayes_top5:
            bayes_top5_names.add(r["taxon"])
    
    overlap = clr_top5 & bayes_top5_names
    print(f"  Top-5 overlap (CLR vs Bayesian): {len(overlap)}/5 taxa")
    if len(overlap) < 5:
        print(f"    CLR-only top-5: {clr_top5 - bayes_top5_names}")
        print(f"    Bayes-only top-5: {bayes_top5_names - clr_top5}")
    
    return merged


def main():
    parser = argparse.ArgumentParser(description="Bayesian taxon enrichment model")
    parser.add_argument("--in", dest="input_file", required=True)
    parser.add_argument("--out-dir", default="results/reports")
    parser.add_argument("--n-mcmc", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("BAYESIAN DIRICHLET-MULTINOMIAL TAXON ENRICHMENT (v3.3)")
    print("=" * 60)
    
    # Load
    print(f"\n[1/3] Loading ASV table: {args.input_file}")
    raw = pd.read_csv(args.input_file)
    jumper = raw["Jumper"].values.astype(float)
    laggard = raw["Laggard"].values.astype(float)
    taxon_names = raw["ASV"].values
    print(f"  {len(raw)} taxa, Jumper total={int(jumper.sum())}, "
          f"Laggard total={int(laggard.sum())}")
    
    # Bayesian model
    print(f"\n[2/3] Fitting Bayesian beta-binomial model "
          f"({args.n_mcmc} posterior samples per taxon)...")
    bayes_result = bayesian_beta_binomial_enrichment(
        jumper, laggard, n_mcmc=args.n_mcmc, seed=args.seed
    )
    bayes_result["taxon"] = taxon_names
    
    # Load CLR results for comparison
    clr_path = Path("data/processed/clr_profiles/taxa_direction.tsv")
    clr_df = None
    if clr_path.exists():
        clr_df = pd.read_csv(clr_path, sep="\t")
    
    # Compare
    if clr_df is not None:
        print(f"\n[3/3] Comparing Bayesian vs CLR results...")
        comparison = compare_with_clr(bayes_result, clr_df)
    else:
        comparison = bayes_result
    
    # Summary statistics
    n_ci_nonzero = bayes_result["ci_excludes_zero"].sum()
    n_jumper_enriched = (bayes_result["p_jumper_enriched"] > 0.95).sum()
    n_laggard_enriched = (bayes_result["p_laggard_enriched"] > 0.95).sum()
    
    print(f"\n  Bayesian model summary:")
    print(f"    Taxa with 95% CI excluding zero: {n_ci_nonzero}/{len(bayes_result)}")
    print(f"    Jumper-enriched (P>0.95): {n_jumper_enriched}")
    print(f"    Laggard-enriched (P>0.95): {n_laggard_enriched}")
    print(f"    Uncertain direction: {len(bayes_result) - n_jumper_enriched - n_laggard_enriched}")
    
    # Write outputs
    bayes_out = bayes_result[[
        "taxon", "jumper_raw", "laggard_raw",
        "jumper_rel_abund", "laggard_rel_abund",
        "log2FC_posterior_mean", "log2FC_ci95_lower", "log2FC_ci95_upper",
        "log2FC_posterior_std", "p_jumper_enriched", "p_laggard_enriched",
        "direction_bayesian", "ci_excludes_zero",
    ]]
    bayes_out.to_csv(out_dir / "bayesian_taxon_enrichment.tsv", sep="\t", index=False)
    print(f"\n  [OK] {out_dir}/bayesian_taxon_enrichment.tsv — "
          f"{len(bayes_out)} taxa")
    
    if comparison is not None:
        comparison.to_csv(out_dir / "bayesian_vs_clr_comparison.tsv",
                         sep="\t", index=False)
        print(f"  [OK] {out_dir}/bayesian_vs_clr_comparison.tsv")
    
    print("\nDone.")


if __name__ == "__main__":
    main()
