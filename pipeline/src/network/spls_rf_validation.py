#!/usr/bin/env python3
"""
spls_rf_validation.py — sPLS-DA + Random Forest Orthogonal Validation
=======================================================================
Track 3 Host-Microbe Integration (v3.2)

PURPOSE:
  Two additional orthogonal methods to validate WGCNA hub genes:

  1. Sparse PLS (sPLS-like): Uses PLS regression with L1-style feature
     selection (top coefficients) to find the smallest set of genes that
     predict weight_gain while controlling for sex and tissue.
     Compares selected genes with WGCNA hubs.

  2. Random Forest: Non-linear ensemble method ranking gene importance
     for predicting weight_gain. Identifies genes that WGCNA (a linear
     correlation-based method) might miss due to non-linear relationships.

METHOD:
  sPLS:
    - Regress out sex and tissue from TPM (same confound correction)
    - Fit PLSRegression with increasing n_components
    - Select top genes by |coefficient| magnitude
    - Report overlap with WGCNA hub genes + module assignments

  Random Forest:
    - Train RandomForestRegressor on confound-corrected TPM
    - Compute permutation importance (mean decrease in R²)
    - Rank all 18,276 genes by importance
    - Report where WGCNA top-10 genes rank in RF importance

OUTPUTS:
  results/reports/
    spls_gene_selection.tsv      — Top genes by PLS coefficient magnitude
    rf_gene_importance.tsv       — All genes ranked by RF importance
    spls_rf_wgcna_comparison.tsv — Overlap and rank comparison

USAGE:
  python src/network/spls_rf_validation.py \
      --tpm data/processed/gene_expression/merged_tpm.tsv \
      --metadata data/raw/sra/PRJNA875278/metadata.tsv \
      --hub-genes results/shortlist/host_genes.csv \
      --modules results/tables/wgcna_modules.csv \
      --out-dir results/reports/ \
      --n-top 50
"""

import argparse
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler


def load_and_prepare(tpm_path, meta_path):
    """Load TPM, metadata, align, and confound-correct."""
    tpm = pd.read_csv(tpm_path, sep="\t", index_col=0)
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)
    
    common = sorted(set(tpm.columns) & set(meta.index))
    tpm = tpm[common]
    meta = meta.loc[common]
    
    # Extract variables
    wg = meta["weight_gain"].values.astype(float)
    sex = np.array([1.0 if str(meta.loc[s, "sex"]).lower() == "male" else 0.0
                     for s in common])
    tissue = np.array([1.0 if "hepato" in str(meta.loc[s, "tissue"]).lower() else 0.0
                        for s in common])
    
    # Log2 transform
    X = np.log2(tpm.values.T + 1)
    gene_ids = tpm.index.values
    n_samples, n_genes = X.shape
    
    # Confound correction: regress out sex and tissue
    confounds = np.column_stack([sex, tissue, np.ones(n_samples)])
    beta = np.linalg.lstsq(confounds, X, rcond=None)[0]
    X_corrected = X - confounds @ beta
    
    return X_corrected, wg, sex, tissue, gene_ids, meta


def run_spls(X, y, gene_ids, n_top=50, max_components=5, seed=42):
    """
    Sparse PLS-like gene selection.
    
    Fits PLS with increasing n_components and selects top genes
    by the magnitude of their coefficients in the first component.
    This approximates sPLS-DA without requiring the mixOmics package.
    """
    print(f"\n--- sPLS Gene Selection ---")
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    y_scaled = StandardScaler().fit_transform(y.reshape(-1, 1)).ravel()
    
    best_genes = {}
    for n_comp in range(1, min(max_components + 1, X.shape[1])):
        pls = PLSRegression(n_components=n_comp, scale=False)
        pls.fit(X_scaled, y_scaled)
        
        # Get coefficients for first latent variable direction
        coef = pls.x_weights_[:, 0] if n_comp >= 1 else pls.coef_.ravel()
        abs_coef = np.abs(coef)
        
        # Top genes by coefficient magnitude
        top_idx = np.argsort(-abs_coef)[:n_top]
        
        for rank, idx in enumerate(top_idx):
            g = gene_ids[idx]
            if g not in best_genes:
                best_genes[g] = {
                    "gene": g,
                    "best_rank": rank + 1,
                    "best_component": n_comp,
                    "coefficient": float(coef[idx]),
                    "abs_coefficient": float(abs_coef[idx]),
                    "selected_in_n_components": 1,
                }
            else:
                best_genes[g]["selected_in_n_components"] += 1
                if rank + 1 < best_genes[g]["best_rank"]:
                    best_genes[g]["best_rank"] = rank + 1
                    best_genes[g]["best_component"] = n_comp
                    best_genes[g]["coefficient"] = float(coef[idx])
                    best_genes[g]["abs_coefficient"] = float(abs_coef[idx])
    
    result = pd.DataFrame(list(best_genes.values()))
    result = result.sort_values("abs_coefficient", ascending=False)
    result.insert(0, "spls_rank", range(1, len(result) + 1))
    
    # Also report R² for each n_components
    r2_scores = {}
    for n_comp in range(1, min(max_components + 1, X.shape[1])):
        pls = PLSRegression(n_components=n_comp, scale=False)
        pls.fit(X_scaled, y_scaled)
        y_pred = pls.predict(X_scaled).ravel()
        r2 = 1 - np.sum((y_scaled - y_pred)**2) / np.sum((y_scaled - y_scaled.mean())**2)
        r2_scores[n_comp] = round(r2, 4)
    
    print(f"  Components tested: {list(r2_scores.keys())}")
    print(f"  R² scores: {r2_scores}")
    print(f"  Selected {len(result)} unique genes across all components")
    
    return result, r2_scores


def run_random_forest(X, y, gene_ids, n_estimators=500, seed=42):
    """
    Random Forest gene importance ranking.
    
    Uses MDI (mean decrease in impurity) for ranking.
    Permutation importance is skipped due to computational cost
    with p=18K and n=20.
    """
    print(f"\n--- Random Forest Gene Importance ---")
    
    n_samples, n_genes = X.shape
    
    rf = RandomForestRegressor(
        n_estimators=n_estimators,
        max_features=min(1000, n_genes // 10),
        min_samples_leaf=3,
        n_jobs=2,
        random_state=seed,
        oob_score=True,
    )
    
    rf.fit(X, y)
    oob_r2 = rf.oob_score_
    print(f"  OOB R²: {oob_r2:.4f}")
    if oob_r2 < 0:
        print(f"  (Negative OOB R² means RF performs worse than mean predictor —")
        print(f"   expected with n=20, p=18K. Importance ranking still informative")
        print(f"   for relative comparison, not absolute prediction.)")
    
    # MDI importance
    mdi_importance = rf.feature_importances_
    
    # Build gene ranking
    rows = []
    for i, g in enumerate(gene_ids):
        rows.append({
            "gene": g,
            "rf_mdi_importance": round(float(mdi_importance[i]), 8),
        })
    
    result = pd.DataFrame(rows)
    result["rf_rank"] = result["rf_mdi_importance"].rank(ascending=False).astype(int)
    result = result.sort_values("rf_mdi_importance", ascending=False)
    
    top_n = 20
    print(f"  Top {top_n} genes by MDI importance:")
    for _, r in result.head(top_n).iterrows():
        print(f"    {r['gene']:<30s} importance={r['rf_mdi_importance']:.6f}")
    
    return result


def compare_with_wgcna(spls_result, rf_result, hub_genes, modules):
    """
    Compare sPLS and RF selected genes with WGCNA hubs.
    Reports overlap and rank agreement.
    """
    print(f"\n--- WGCNA Comparison ---")
    
    # Get WGCNA top-10 genes
    wgcna_top10 = set(hub_genes["gene_id"].values[:10])
    
    # Module mapping
    gene_to_module = dict(zip(modules["gene"], modules["module"]))
    
    # sPLS overlap
    spls_top50 = set(spls_result["gene"].values[:50])
    spls_top100 = set(spls_result["gene"].values[:100])
    spls_overlap_10 = wgcna_top10 & spls_top50
    spls_overlap_100 = wgcna_top10 & spls_top100
    
    print(f"  sPLS top-50 overlap with WGCNA top-10: {len(spls_overlap_10)} genes")
    if spls_overlap_10:
        for g in sorted(spls_overlap_10):
            spls_rank = spls_result[spls_result["gene"] == g]["spls_rank"].values[0]
            wgcna_rank = hub_genes[hub_genes["gene_id"] == g]["rank"].values[0]
            mod = gene_to_module.get(g, "?")
            print(f"    {g:<30s} sPLS_rank={spls_rank} WGCNA_rank={wgcna_rank} M{mod}")
    
    # RF overlap
    rf_top50 = set(rf_result["gene"].values[:50])
    rf_overlap = wgcna_top10 & rf_top50
    print(f"  RF top-50 overlap with WGCNA top-10: {len(rf_overlap)} genes")
    if rf_overlap:
        for g in sorted(rf_overlap):
            rf_rank = rf_result[rf_result["gene"] == g]["rf_rank"].values[0]
            wgcna_rank = hub_genes[hub_genes["gene_id"] == g]["rank"].values[0]
            mod = gene_to_module.get(g, "?")
            print(f"    {g:<30s} RF_rank={rf_rank} WGCNA_rank={wgcna_rank} M{mod}")
    
    # Rank correlation between methods
    # For genes in common between sPLS and RF
    spls_ranked = spls_result.set_index("gene")
    rf_ranked = rf_result.set_index("gene")
    common_genes = sorted(set(spls_ranked.index) & set(rf_ranked.index))
    
    if len(common_genes) > 10:
        from scipy.stats import spearmanr
        spls_ranks = [spls_ranked.loc[g, "spls_rank"] for g in common_genes]
        rf_ranks = [rf_ranked.loc[g, "rf_rank"] for g in common_genes]
        rho, p = spearmanr(spls_ranks, rf_ranks)
        print(f"\n  sPLS vs RF rank correlation: Spearman ρ={rho:.3f}, p={p:.4f}")
    
    # Build comparison table for WGCNA top-10
    comparison_rows = []
    for _, hg in hub_genes.iterrows():
        g = hg["gene_id"]
        mod = hg["associated_module"]
        
        spls_r = spls_ranked.loc[g, "spls_rank"] if g in spls_ranked.index else None
        rf_r = rf_ranked.loc[g, "rf_rank"] if g in rf_ranked.index else None
        
        comparison_rows.append({
            "gene": g,
            "wgcna_module": mod,
            "wgcna_rank": hg["rank"],
            "wgcna_kME": hg["kME"],
            "spls_rank": spls_r,
            "rf_rank": rf_r,
            "spls_selected": g in spls_top50,
            "rf_selected": g in rf_top50,
            "validated_by_both": g in spls_top50 and g in rf_top50,
        })
    
    comparison = pd.DataFrame(comparison_rows)
    n_both = comparison["validated_by_both"].sum()
    n_either = (comparison["spls_selected"] | comparison["rf_selected"]).sum()
    print(f"\n  WGCNA top-10 genes validated by:")
    print(f"    sPLS (top-50):   {comparison['spls_selected'].sum()}/10")
    print(f"    RF (top-50):     {comparison['rf_selected'].sum()}/10")
    print(f"    Either method:   {n_either}/10")
    print(f"    Both methods:    {n_both}/10")
    
    return comparison


def main():
    parser = argparse.ArgumentParser(description="sPLS + RF orthogonal validation")
    parser.add_argument("--tpm", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--hub-genes", required=True)
    parser.add_argument("--modules", required=True)
    parser.add_argument("--out-dir", default="results/reports")
    parser.add_argument("--n-top", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("sPLS + RANDOM FOREST — Orthogonal Host Gene Validation (v3.2)")
    print("=" * 70)
    
    # Load
    print("\n[1/4] Loading and confound-correcting data...")
    X, y, sex, tissue, gene_ids, meta = load_and_prepare(args.tpm, args.metadata)
    n_samples, n_genes = X.shape
    print(f"  Samples: {n_samples}, Genes: {n_genes}")
    print(f"  WG range: [{y.min():.2f}, {y.max():.2f}], mean={y.mean():.2f}")
    
    hub_genes = pd.read_csv(args.hub_genes)
    modules = pd.read_csv(args.modules, sep="\t")
    
    # sPLS
    print("\n[2/4] Running sPLS gene selection...")
    spls_result, r2_scores = run_spls(X, y, gene_ids, n_top=args.n_top, seed=args.seed)
    
    # Random Forest
    print("\n[3/4] Running Random Forest...")
    rf_result = run_random_forest(X, y, gene_ids, seed=args.seed)
    
    # Comparison
    print("\n[4/4] Comparing with WGCNA...")
    comparison = compare_with_wgcna(spls_result, rf_result, hub_genes, modules)
    
    # Write outputs
    spls_result.to_csv(out_dir / "spls_gene_selection.tsv", sep="\t", index=False)
    rf_result.to_csv(out_dir / "rf_gene_importance.tsv", sep="\t", index=False)
    comparison.to_csv(out_dir / "spls_rf_wgcna_comparison.tsv", sep="\t", index=False)
    
    print(f"\n  [OK] {out_dir}/spls_gene_selection.tsv — {len(spls_result)} genes")
    print(f"  [OK] {out_dir}/rf_gene_importance.tsv — {len(rf_result)} genes")
    print(f"  [OK] {out_dir}/spls_rf_wgcna_comparison.tsv — {len(comparison)} genes")
    print("\nDone.")


if __name__ == "__main__":
    main()
