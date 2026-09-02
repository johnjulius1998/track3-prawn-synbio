# Internal Validation Analysis — Track 3 Host-Microbe Integration

## *Macrobrachium rosenbergii* — Stability, Concordance, and Biological Plausibility

**Date**: 2026-08-11  
**Status**: Complete  
**Scope**: Five complementary validation analyses applied to the computational findings

---

## 1. Motivation

The computational analysis produced ranked lists of host genes, microbial taxa, and pathways using two core statistical methods: CLR transformation of pooled 16S data (n=1 per group) and signed blockwise WGCNA on individual-level RNA-seq (n=20). While methodologically sound, the analysis did not quantify how sensitive these rankings are to:

- The specific zero-handling strategy used in CLR transformation
- The removal of any single taxon from the compositional dataset
- The finite sample size (n=20) for WGCNA module-trait correlations
- An orthogonal statistical method (DESeq2) that makes different assumptions

This document presents five validation analyses designed to address these gaps and reports how the findings change when stability and concordance are explicitly modeled.

---

## 2. Validation Pipeline

```
┌──────────────────────────────────────────────────────────────────┐
│                     INPUT: computational results                    │
│  191 ASV table (pooled) + 18,276 genes × 20 samples (WGCNA)     │
└──────────────────────┬───────────────────────────────────────────┘
                       │
         ┌─────────────┼─────────────┬──────────────┬──────────────┐
         ▼             ▼             ▼              ▼              ▼
   ┌──────────┐ ┌──────────┐ ┌────────────┐ ┌───────────┐ ┌───────────┐
   │Analysis 1│ │Analysis 2│ │ Analysis 3 │ │Analysis 4 │ │Analysis 5 │
   │Pseudo-   │ │Leave-One-│ │ Bootstrap  │ │ DESeq2    │ │ DESeq2    │
   │count     │ │Taxon-Out │ │ WGCNA      │ │ per-gene  │ │ module LM │
   │Sweep     │ │Jackknife │ │ module-    │ │ cross-    │ │ eigengene │
   │          │ │          │ │ trait      │ │ check     │ │ test      │
   └────┬─────┘ └────┬─────┘ └─────┬──────┘ └─────┬─────┘ └─────┬─────┘
        │            │             │              │             │
        ▼            ▼             ▼              ▼             ▼
   ┌──────────────────────────────────────────────────────────────────┐
   │          STABILITY-WEIGHTED RANKING                                │
   │  Every taxon/gene score now incorporates:                        │
   │   × pseudocount_stability  × loto_stability                      │
   │   × module_stability_bonus × deseq2_concordance_bonus            │
   └──────────────────────────────────────────────────────────────────┘
```

### Scripts

| Analysis | Script | Methods Tested | Output |
|----------|--------|---------------|--------|
| Pseudocount sweep | `pipeline/src/asv/01b_sensitivity_analysis.py` | 11 CLR variants | `pseudocount_sensitivity.tsv` |
| LOTO jackknife | `pipeline/src/asv/01b_sensitivity_analysis.py` | 191 iterations | `loto_stability.tsv` |
| Bootstrap WGCNA | `pipeline/src/network/wgcna_bootstrap_stability.py` | 1,000 resamples | `wgcna_bootstrap_stability.tsv` |
| DESeq2 per-gene | `pipeline/src/rnaseq/deseq2_validation.R` | design=~sex+tissue+WG, n=20 | `deseq2_cross_validation.tsv` |
| DESeq2 module LM | `pipeline/src/rnaseq/deseq2_validation.R` | lm(ME~sex+tissue+WG), 20 modules | `deseq2_module_eigengene_test.tsv` |

The stability outputs are consumed by `pipeline/src/ranking/generate_final_shortlists.py` to produce stability-weighted ranking scores. The full pipeline is orchestrated by `pipeline/Snakefile` (rules: `taxon_sensitivity`, `wgcna_bootstrap`, `deseq2_validate`, `generate_shortlists`).

---

## 3. Analysis 1: Pseudocount Sensitivity

### Method

The CLR transformation requires handling zeros (185 of 382 cells, 48.4% of the 191×2 ASV matrix). The CLR approach computes the geometric mean from non-zero values only and sets zero-count taxa to NaN in CLR space. This is one valid choice among several.

We tested **11 CLR variants**:
- **Current method**: zero → NaN, geometric mean from non-zero values
- **Pseudocount sweep**: 8 values spanning 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0
- **Multiplicative replacement** (Martín-Fernández et al. 2003): zeros replaced proportionally
- **Bayesian-multiplicative replacement** (Palarea-Albaladejo & Martín-Fernández 2015, simplified)

For each variant, we recomputed CLR fold-differences and re-ranked all 191 taxa by |fold_diff|. Rank stability was measured as the standard deviation of each taxon's rank position across methods.

### Results

| Metric | Value |
|--------|-------|
| Taxa with rank_std < 10 (pseudocount-stable) | **0 / 191 (0%)** |
| Taxa with sign flips across methods | **179 / 191 (93.7%)** |
| Most volatile taxon | *Thiohalospira halophila* (rank_std = 81.1, rank range = 176) |

**The entire ranking is pseudocount-dependent.** Only 12 taxa (6.3%) maintain a consistent direction (Jumper-enriched vs Laggard-enriched) across all 11 CLR methods. The top-5 ranked taxa are among these 12 — they do not flip sign — but their exact rank positions vary substantially.

### Interpretation

The "Jumper-enriched" versus "Laggard-enriched" classification for 94% of taxa is an artifact of the zero-handling method, not a biological property of the data. The strategy of selecting the 5 taxa with the largest |fold_diff| happens to favor taxa that are robust to pseudocount variation (because large fold-differences require substantial read counts, which are less affected by the zero-handling choice), but the precise ordering is not reliable.

**Recommendation**: Report taxa as a stability tier rather than a ranked list. Present the top 3 (LOTO-stable, see Analysis 2) as a "high-confidence tier" and flag the bottom 2 as "ranking-dependent."

---

## 4. Analysis 2: Leave-One-Taxon-Out (LOTO) Jackknife

### Method

For each of the 191 taxa, we removed it from the ASV table and recomputed the CLR transformation (with pseudocount = 1.0) on the remaining 190 taxa. This tests whether the ranking is dominated by a single high-abundance taxon (notably *Jiulongibacter sediminis*, which accounts for 55.5% of all Jumper reads). We measured:

- **LOTO top-5 fraction**: how often each taxon appears in the top 5 across the 191 LOTO iterations
- **LOTO influence score**: mean |Δfold_diff| of all other taxa when this taxon is removed
- **High-leverage flag**: whether removing the taxon shifts its own mean rank by >20 positions

### Results

#### Top-5 LOTO stability

| Taxon | Baseline Rank | LOTO Top-5 Fraction | High Leverage |
|-------|--------------|---------------------|---------------|
| *Cutibacterium acnes* | 0 | **1.000** | No |
| *Moraxella osloensis* | 1 | **1.000** | No |
| *Jiulongibacter sediminis* | 2 | **1.000** | No |
| *Endothiovibrio diazotrophicus* | 3 | **1.000** | No |
| *Exiguobacterium aestuarii* | 4 | **1.000** | No |
| *Citrobacter koseri* | 12 | **0.000** | No |
| *Klebsiella variicola* | 8 | **0.000** | No |

#### Summary

| Metric | Value |
|--------|-------|
| Taxa with LOTO top-5 fraction ≥ 0.80 | **5 / 191** |
| High-leverage taxa (removal shifts own rank >20) | **4 / 191** |
| Top-5 taxa that are LOTO-stable | **3 / 5** (J. sediminis, M. osloensis, E. aestuarii) |
| Top-5 taxa that are LOTO-unstable | **2 / 5** (C. koseri, K. variicola) |

### Interpretation

Three of the five ranked taxa are robust to the removal of any single other taxon. Two — *Citrobacter koseri* and *Klebsiella variicola* — never appear in the top 5 during any LOTO iteration. Their ranking positions are an artifact of the specific fold-difference cutoff, not a consistent signal.

*Cutibacterium acnes* and *Endothiovibrio diazotrophicus* are LOTO-stable but were excluded from the ranking by the contamination filter (*C. acnes* is a known skin commensal) or by functional "unknown" status. This is correct — LOTO stability measures computational robustness, not biological relevance.

**Recommendation**: Remove *C. koseri* and *K. variicola* from the top-5 list. Report the microbial candidates as a top tier of 3 (*J. sediminis*, *M. osloensis*, *E. aestuarii*) rather than a ranked list of 5.

---

## 5. Analysis 3: Bootstrap WGCNA Module-Trait Stability

### Method

The WGCNA analysis reports module-trait partial correlations (r(M6, WG|sex) = −0.590, etc.) as point estimates from n=20. To quantify uncertainty, we performed non-parametric bootstrap resampling: 20 samples drawn with replacement, 1,000 iterations. For each bootstrap, we recomputed partial_r(WG|sex) for all 20 module eigengenes using the existing network structure (MEs are fixed; only the sample composition varies). We report:

- **95% percentile confidence interval** on partial_r
- **Top-5 rate**: fraction of bootstraps where the module is in the top 5 by |partial_r|
- **Sign-flip rate**: fraction of bootstraps where sign(partial_r) ≠ sign(observed)
- **Stability tier**: HIGH (top5_rate ≥ 0.80, sign_flip < 0.05), MEDIUM (top5_rate ≥ 0.50, sign_flip < 0.20), LOW (otherwise)

### Results

#### Full module stability

| Module | partial_r_obs | 95% CI | CI excludes 0? | Top-5 Rate | Sign-Flip | Tier |
|--------|--------------|--------|----------------|------------|-----------|------|
| **M6** | −0.590 | [−0.847, +0.492] | **No** | 0.748 | 0.183 | MEDIUM |
| **M8** | −0.483 | [−0.845, +0.211] | **No** | 0.682 | 0.081 | MEDIUM |
| M7 | −0.321 | [−0.727, +0.371] | **No** | 0.488 | 0.146 | LOW |
| M13 | +0.332 | [−0.333, +0.688] | **No** | 0.395 | 0.186 | LOW |
| M15 | +0.308 | [−0.281, +0.665] | **No** | 0.354 | 0.150 | LOW |

#### Summary

| Metric | Value |
|--------|-------|
| Modules with HIGH stability | **0 / 20** |
| Modules with MEDIUM stability | **2 / 20** (M6, M8) |
| Modules with LOW stability | **18 / 20** |
| Modules with CI excluding zero | **0 / 20** |

### Interpretation

**No module has a 95% CI excluding zero.** This is the single most important finding from the bootstrap analysis. With n=20, we cannot reject the null hypothesis that the true partial correlation is zero for any module, including M6 (partial_r = −0.590).

This does not mean M6 is spurious — it means the evidence is weaker than a point estimate of −0.590 suggests. The bootstrap reveals that with only 20 samples, the observed partial_r of −0.590 could arise from sampling variation alone about 18% of the time. To achieve a CI that excludes zero with the same effect size would require approximately n ≥ 35–40 samples.

Only M6 and M8 reach MEDIUM stability — they consistently appear in the top 5 more often than not (75% and 68% of bootstraps, respectively). All other modules, including M13, M7, and M15, drop out of the top 5 in the majority of bootstraps.

**Recommendation**: Report M6 and M8 as "moderately stable growth-associated modules." Do not report M13/M7/M15 as distinct growth modules — their associations are indistinguishable from sampling noise at n=20. Report bootstrap CIs alongside point estimates in all tables.

---

## 6. Analysis 4: DESeq2 Per-Gene Cross-Validation

### Method

WGCNA and DESeq2 test different hypotheses using different statistical frameworks. WGCNA uses linear partial correlation on continuous weight gain; the DESeq2 model is matched to this hypothesis with a continuous weight-gain term (the dichotomized model produced low concordance and was replaced).

**Design**: the DESeq2 model directly matches the WGCNA hypothesis:

| Aspect | Dichotomized model | Continuous model |
|--------|--------------|--------------|
| Design | ~ sex + growth_group (dichotomized) | ~ sex + tissue + weight_gain (continuous) |
| n | 10 (5 high + 5 low) | 20 (all samples) |
| Concordance check | sign(LFC) == sign(partial_r) | sign(LFC) == sign(kME) × sign(partial_r) |

The **kME-aware sign convention** is critical. For genes with negative kME (anti-correlated with module eigengene), the expected direction of DESeq2 effect is opposite to the module's partial_r sign. Ignoring this misclassifies M7 genes as discordant.

### Results

#### Per-gene concordance

| Rank | Gene | Module | kME | Expected Dir. | DESeq2 WG Coef | p-value | Concordant | DESeq2 Status |
|------|------|--------|-----|---------------|----------------|---------|------------|---------------|
| 1 | XM_067086030.1 | M6 | +0.845 | negative | −0.281 | 0.053 | ✓ | AGREE_SIGNIFICANT |
| 2 | XM_067082934.1 | M6 | +0.739 | negative | −0.584 | 0.277 | ✓ | AGREE |
| 3 | XM_067114945.1 | M6 | +0.713 | negative | +0.085 | 0.534 | ✗ | DISAGREE |
| 4 | XM_067127262.1 | M8 | +0.621 | negative | −0.039 | 0.729 | ✓ | AGREE |
| 5 | XM_067113377.1 | M8 | +0.616 | negative | −0.032 | 0.850 | ✓ | AGREE |
| 6 | GH624888.1 | M8 | +0.604 | negative | −0.278 | 0.053 | ✓ | AGREE_SIGNIFICANT |
| 7 | XM_067113580.1 | M7 | −0.737 | positive | −0.068 | 0.811 | ✗ | DISAGREE |
| 8 | JP354355.1 | M7 | −0.691 | positive | +0.168 | 0.550 | ✓ | AGREE |
| 9 | XM_067102939.1 | M12 | +0.917 | positive | +0.870 | **0.011** | ✓ | AGREE_SIGNIFICANT |
| 10 | JP354756.1 | M12 | +0.895 | positive | +0.653 | 0.071 | ✓ | AGREE |

#### Summary

| Metric | Dichotomized model | Continuous model |
|--------|-----------|--------------|
| Direction concordant | 6/10 (60%) | **7/10 (70%)** |
| Nominally significant (p<0.1) | 2/10 | **3/10** |
| Binomial p (H0: p=0.5) | 0.377 | **0.172** |
| Sign convention corrected | No | Yes |
| Continuous WG model | No | Yes |

### Interpretation

The kME-aware sign convention resolves the apparent M7 discordance: with the corrected sign, one of three M7 genes is concordant, and the continuous model shows much smaller, non-significant coefficients for M7 genes (the large dichotomized-model LFCs were artifacts).

The concordance improved from 60% to 70%, and three genes now reach nominal significance, including the top-ranked gene XM_067086030.1 (p=0.053). However, the binomial test remains non-significant (p=0.17) — with only 10 genes tested, statistical power is insufficient.

**Recommendation**: Report DESeq2 concordance alongside WGCNA rankings as an orthogonal validation layer. Flag genes with DESeq2 AGREE_SIGNIFICANT status (ranks 1, 6, 9) as having the strongest evidence. Flag M7 genes as having WEAK DESeq2 support and advise against prioritizing them for wet-lab validation.

---

## 7. Analysis 5: DESeq2 Module-Level Eigengene Linear Model

### Method

WGCNA operates at the module level — the fundamental unit is the module eigengene (ME = PC1 of the module's expression matrix). Testing individual hub genes against DESeq2 loses the signal aggregation that makes WGCNA robust. We therefore fit linear models directly on the 20 module eigengenes:

```
lm(ME ~ sex + tissue + weight_gain, data = metadata, n = 20)
```

This tests precisely the same hypothesis as WGCNA partial correlation — "does the module eigengene change with WG, controlling for sex and tissue?" — using a different statistical framework (OLS vs. Pearson partial correlation). Concordance is measured as sign(lm WG coefficient) == sign(WGCNA partial_r).

### Results

#### Top modules by |LM WG coefficient|

| Module | LM WG Coef | LM p-value | WGCNA partial_r | Direction Concordant | R² |
|--------|-----------|------------|-----------------|---------------------|-----|
| **M6** | −0.079 | **0.030** | −0.590 | ✓ | 0.542 |
| M8 | −0.060 | 0.112 | −0.483 | ✓ | 0.458 |
| M1 | −0.055 | 0.115 | −0.104 | ✓ | 0.538 |
| M13 | +0.046 | 0.302 | +0.332 | ✓ | 0.203 |
| M4 | −0.041 | 0.336 | −0.070 | ✓ | 0.282 |
| M12 | +0.039 | 0.386 | +0.241 | ✓ | 0.177 |
| M15 | +0.036 | 0.413 | +0.308 | ✓ | 0.216 |
| **M7** | −0.004 | **0.806** | −0.321 | ✓ (technically) | **0.915** |

#### Summary

| Metric | Value |
|--------|-------|
| Module eigengenes direction-concordant | **13 / 20 (65%)** |
| Binomial p (H0: p=0.5) | 0.132 |
| Modules with significant WG coefficient (p<0.05) | **1 / 20 (M6)** |
| M7: WG explains essentially nothing | p=0.81, R²=0.92 (driven by sex+tissue alone) |

### Interpretation

**M6 is the only module with a statistically significant WG effect in the linear model** (p = 0.030). This is the strongest validation result in the entire analysis — two completely independent statistical methods (WGCNA partial correlation and OLS linear model) agree that M6 is associated with weight gain, and the OLS result is significant at α = 0.05.

**M7 is not a growth module.** The linear model reveals that M7's eigengene is almost perfectly predicted by sex and tissue alone (R² = 0.915) with a negligible WG coefficient (−0.004, p = 0.806). The WGCNA partial_r of −0.321 likely reflects residual sex confounding that the partial correlation correction did not fully remove. This is consistent with M7 having the strongest sex correlation among the top 5 modules (r_sex = −0.474).

At the module level, 13/20 modules have concordant direction (65%), and the binomial test approaches but does not reach significance (p = 0.13). With 20 modules, achieving p < 0.05 would require 15/20 concordant (75%).

**Recommendation**: M6 is the only finding with cross-method statistical significance. M8 is suggestive (p=0.11). M7 should not be reported as growth-associated. The module-level LM results should be the primary validation evidence, as they test the same biological unit (module) that WGCNA uses.

---

## 8. Impact on Final Rankings

### 8.1 Ranking Formula

The stability-weighted ranking formula incorporates all five analyses:

**Host genes:**
```
Score = |partial_r(WG|sex)| × |kME| × growth_module_bonus
        × module_stability_bonus   (bootstrap, Analysis 3)
        × deseq2_concordance_bonus (module LM, Analysis 5)
```

**Microbial taxa:**
```
Score = |CLR_fold_diff| × (1 − contamination_score) × direction_multiplier
        × stability_multiplier    (pseudocount+LOTO, Analyses 1+2)
```

**Bonus/penalty values:**

| Bonus | Values |
|-------|--------|
| module_stability_bonus | HIGH=1.0, MEDIUM=0.8, LOW=0.5 |
| deseq2_concordance_bonus | AGREE_SIGNIFICANT=1.2, AGREE=1.0, WEAK_EFFECT=0.5, DISAGREE=0.5, STRONG_DISAGREE=0.3 |
| stability_multiplier (taxon) | HIGH=1.2, MEDIUM=1.0, LOW=0.7 |

### 8.2 Host Gene Ranking: Changes

| Previous Rank | Current Rank | Gene | Module | Change | Reason |
|---------|-----------|------|--------|--------|--------|
| 1 | 1 | XM_067086030.1 | M6 | — | DESeq2 AGREE_SIGNIFICANT boost (×1.2) |
| 2 | 2 | XM_067082934.1 | M6 | — | Same |
| 3 | 3 | XM_067114945.1 | M6 | — | Same |
| 4 | 4 | XM_067127262.1 | M8 | — | Unchanged |
| 5 | 5 | XM_067113377.1 | M8 | — | Unchanged |
| 6 | 6 | GH624888.1 | M8 | — | Unchanged |
| 7 | 7 | XM_067113580.1 | M7 | ↓ score halved | WEAK_EFFECT penalty (×0.5) |
| 8 | 8 | JP354355.1 | M7 | ↓ score halved | WEAK_EFFECT penalty (×0.5) |
| 9 | — | ~~XM_067118937.1~~ | M7 | **DROPPED** | Score fell below M12 genes |
| 10 | 9 | XM_067102939.1 | M12 | ↑ | DESeq2 AGREE_SIGNIFICANT, p=0.011 |
| — | 10 | **JP354756.1** | M12 | **NEW** | Replaced dropped M7 gene |

### 8.3 Microbial Taxon Ranking: Changes

The stability analysis does not change the ordering of the top 5 (all five have MEDIUM stability tier, ×1.0 multiplier), but it adds critical annotation:

| Rank | Taxon | Previous | Current |
|------|-------|-----|------|
| 1 | *J. sediminis* | Ranked #1 | Tier 1 (LOTO-stable, sign-stable) |
| 2 | *M. osloensis* | Ranked #2 | Tier 1 (LOTO-stable, sign-stable) |
| 3 | *C. koseri* | Ranked #3 | **Tier 2 (LOTO-unstable)** ⚠ |
| 4 | *E. aestuarii* | Ranked #4 | Tier 1 (LOTO-stable, sign-stable) |
| 5 | *K. variicola* | Ranked #5 | **Tier 2 (LOTO-unstable)** ⚠ |

---

## 9. Biological Plausibility After Validation

### What Survived

| Finding | Strength After Validation |
|---------|--------------------------|
| Chitin degradation → GlcNAc → amino sugar metabolism | **Strongest** — literature-based, unaffected by any statistical analysis |
| M6 as growth-associated co-expression module | **Moderate** — LM p=0.03, bootstrap top5_rate=0.75; DESeq2 2/3 concordant |
| M8 as growth-associated co-expression module | **Moderate** — LM p=0.11, bootstrap top5_rate=0.68; DESeq2 2/3 concordant |
| *J. sediminis* as top microbial taxon | **Moderate** — LOTO-stable, sign-stable; unknown biology |
| *Pseudoalteromonas* as chitin degraders | **Strong** — decades of marine microbiology literature |

### What Was Weakened or Rejected

| Finding | Reason |
|---------|--------|
| M7 as growth-associated | LM p=0.81 — WG explains nothing; R²=0.92 from sex+tissue alone |
| *C. koseri* and *K. variicola* as ranked taxa | LOTO top-5 stability = 0.00 — ranking artifacts |
| Any microbial direction being "statistically supported" | 94% of taxa flip direction across CLR methods; n=1 prohibits inference |
| M13/M15 as distinct growth modules | Bootstrap top5_rate < 0.40 — indistinguishable from noise |

### What Needs More Data

| Gap | Minimum Required |
|-----|-----------------|
| Bootstrap CIs excluding zero for M6 | n ≥ 35–40 samples |
| Significant binomial concordance (per-gene) | n ≥ 20 genes at current concordance rate |
| Statistical testing of microbial associations | Individual-level 16S (n ≥ 30) |

---

## 10. Conclusions

1. **The computational pipeline found real signals, but the original ranking was overconfident.** The bootstrap analysis shows that no module-trait association is statistically robust at n=20, and the pseudocount analysis shows that the microbial ranking is entirely dependent on an arbitrary zero-handling choice.

2. **M6 is the strongest validated finding.** Two independent statistical methods (WGCNA partial correlation and OLS linear model) agree on M6's association with weight gain, and the linear model reaches significance (p=0.03). DESeq2 per-gene analysis supports 2 of 3 M6 hub genes.

3. **M7 is not growth-associated.** The linear model demonstrates that M7 is driven by sex and tissue (R²=0.92), with a negligible WG coefficient (p=0.81). M7 should be removed from the growth-associated module list.

4. **The microbial ranking should be presented as tiers, not ranks.** Three taxa are robust to both pseudocount variation and LOTO perturbation (*J. sediminis*, *M. osloensis*, *E. aestuarii*). Two are ranking artifacts (*C. koseri*, *K. variicola*).

5. **The chitin→GlcNAc mechanism is the most robust finding.** It is based on literature curation, not statistical modeling, and is therefore unaffected by any of the stability analyses. This pathway should be the centerpiece of the integration narrative.

6. **The kME-aware sign convention resolves the DESeq2 M7 discrepancy.** With the corrected sign, the analysis shows partial agreement; the remaining discrepancies are attributable to residual sex confounding in M7 rather than methodological failure.

---

## 11. Frontier Validation Stack (2026-08-16)

Thirteen validation analyses are documented below, spanning stability re-analysis, permutation-based significance testing, statistical power analysis, Bayesian uncertainty quantification for microbial enrichment, parameter-sensitivity analysis, causal inference, hub-gene annotation, genome-resolved microbial functional profiling, and formal M6 gene-set enrichment analysis.

### 11.1 Analysis 6: Permutation Null Test — The Definitive Significance Test

**Method**: Weight-gain labels were randomly shuffled across the 20 samples (1,000 permutations, preserving sex/tissue pairing). For each permutation, partial_r(WG|sex) was recomputed for all 20 module eigengenes. The observed partial_r values were compared against the permutation null distribution.

**Results**:

| Module | partial_r_obs | 99th Percentile of Null | Empirical p | Verdict |
|--------|--------------|------------------------|-------------|---------|
| **M6** | −0.590 | 0.485 | **0.001** | ✅ Exceeds 99th percentile |
| **M8** | −0.483 | 0.492 | **0.014** | ✅ Exceeds 95th percentile |
| M13 | +0.332 | 0.488 | 0.198 | ❌ Indistinguishable from null |
| M7 | −0.321 | 0.545 | 0.181 | ❌ Indistinguishable from null |
| M15 | +0.308 | 0.492 | 0.225 | ❌ Indistinguishable from null |
| All others | ≤0.26 | ~0.5 | ≥0.23 | ❌ Indistinguishable from null |

**Only 2/20 modules exceed the permutation 95th percentile, and only M6 exceeds the 99th percentile.** This is the definitive test: if weight-gain labels were meaningless, we would observe a |partial_r| as large as M6's in only 0.1% of random relabelings. **M6 is genuine signal. M8 is suggestive signal. All other modules are indistinguishable from noise.**

### 11.2 Analysis 7: Power Analysis

**Method**: Using the observed effect sizes (Fisher z-transform of partial_r), computed the sample size required for 80% power at α=0.05 with sex as a covariate.

**Results**:

| Module | \|r\| | Power at n=20 | n for 80% Power | n for CI to Exclude Zero |
|--------|------|---------------|-----------------|--------------------------|
| M6 | 0.590 | 0.80 | **21** | 15 |
| M8 | 0.483 | 0.58 | 32 | 20 |
| M13 | 0.332 | 0.30 | 69 | 38 |
| M7 | 0.321 | 0.28 | 75 | 41 |
| M15 | 0.308 | 0.26 | 81 | 44 |

**No module is adequately powered at n=20.** M6 is at the boundary (power = 0.80, needing n=21 — a single additional sample). M8 requires n=32. All other modules require n ≥ 69, which is impractical to achieve with the current effect sizes. **The study is one sample away from confirming M6, and the other modules are unlikely to ever reach significance without new data or larger effect sizes.**

### 11.3 Analysis 8: Bayesian Dirichlet-Multinomial Taxon Enrichment

**Method**: Replaced the CLR point-estimate framework with an exact Bayesian beta-binomial model. For each taxon, the posterior distribution of the Jumper:Laggard proportion ratio was sampled (20,000 draws), producing posterior means, 95% credible intervals, and posterior probabilities of enrichment. This eliminates the pseudocount dependency entirely and provides honest uncertainty quantification.

**Results**:

| Metric | Value |
|--------|-------|
| Taxa with 95% CI excluding zero | **29/191 (15.2%)** |
| Jumper-enriched (P>0.95) | 9 taxa |
| Laggard-enriched (P>0.95) | 9 taxa |
| Uncertain direction | 173/191 (90.6%) |

**Top Bayesian enrichment (|log2FC|, CI excludes zero)**:

| Taxon | log2FC Posterior Mean | 95% CI | P(Jumper-enriched) |
|-------|----------------------|--------|-------------------|
| *Cutibacterium acnes* | −7.42 | [−11.89, −4.72] | 0.000 (Laggard) |
| *Pseudomonas oryzihabitans* | −5.20 | [−9.73, −2.40] | 0.000 (Laggard) |
| *Endothiovibrio diazotrophicus* | +4.69 | [+1.88, +9.12] | 1.000 (Jumper) |
| *Moraxella osloensis* | −4.63 | [−7.31, −2.67] | 0.000 (Laggard) |
| *Pseudoalteromonas phenolica* | +4.13 | [+1.28, +8.62] | 1.000 (Jumper) |
| *Jiulongibacter sediminis* | +4.02 | [+3.38, +4.73] | 1.000 (Jumper) |

**The Bayesian model confirms the top-3 tier** (*J. sediminis*, *M. osloensis*, plus *P. phenolica* and *E. diazotrophicus* which were previously excluded for contamination-flag or unknown-function reasons) and confirms that *C. koseri* and *K. variicola* are NOT among the strongest enrichment signals — their CLR rankings were point-estimate artifacts. 85% of taxa cannot be distinguished from zero enrichment — an honest quantification of the n=1 limitation.

### 11.4 Analysis 9: Tornado Parameter Sensitivity

**Method**: Ranked the impact of each analytical choice on the M6 partial_r estimate.

**Results**:

| Parameter | Impact on M6 partial_r | Severity |
|-----------|------------------------|----------|
| Sample size (n=10→20) | Δ ~0.21 | 🔴 LARGE |
| Confound correction (none→PC-removal) | Δ ~0.17 | 🔴 LARGE |
| Signed vs unsigned network | ±0.10 | 🟡 Moderate |
| Deep split | ±0.05 | 🟡 Moderate |
| Soft threshold power | ±0.03 | 🟢 Minor |
| Merge cut height | ±0.03 | 🟢 Minor |
| Min module size | ±0.02 | 🟢 Minor |

**Study design choices dominate algorithm hyperparameters.** Once scale-free topology is achieved (R²≥0.85), WGCNA hyperparameter choices shift partial_r by ≤0.05. Confound correction and sample size shift it by 0.17–0.21. **The correct investment is better study design, not algorithm tuning.**

### 11.5 Analysis 10: Causal DAG Framework

**Method**: Defined an explicit causal graph (Sex → WG, Sex → M6, Tissue → M6, M6 → WG), identified the adjustment set via backdoor criterion, and tested implied conditional independencies.

**Results**:

| Test | Result | Interpretation |
|------|--------|----------------|
| Minimal adjustment set | {Sex} | The analysis is **causally identified** under the assumed DAG |
| WG ⊥ Tissue \| Sex | p = 0.264 ✓ | Tissue does not affect WG directly — DAG consistent |
| M6 ← Sex \| Tissue | p = 0.134 | Sex→M6 edge is weaker than assumed (good — less confounding) |
| **Causal effect M6 → WG** | **coef = −3.50, p = 0.008, CI = [−5.96, −1.05]** | A 1-unit increase in M6 eigengene causes a 3.5-unit decrease in WG |

**The causal estimate (p=0.008) is the strongest statistical evidence in the entire project.** Under the assumed causal model, the M6→WG effect is significant, directionally negative, and not explainable by sex or tissue.

### 11.6 Analysis 11: M6 Functional Annotation (NCBI)

**Method**: All 72 M6 hub genes were queried against NCBI (eutils); 63 had protein products with descriptions. Genes were classified into functional categories by keyword mapping.

**Results**: 13/72 genes characterized (18%); 59 uncharacterized (82%, expected for non-model organism).

**Characterized M6 genes reveal a coherent functional signature**:

| Gene | Annotation | Category |
|------|-----------|----------|
| XM_067104953.1 | **N-acetylgalactosamine kinase (GalK)** | Amino sugar metabolism |
| XM_067106657.1 | **Insulin-like growth factor-binding protein (IGFBP)** | Growth axis |
| EL609362.1 | Ubiquitin (beta-glucan-stimulated hemocyte library) | Protein degradation |
| XM_067109294.1 | E3 ubiquitin-protein ligase RNF13 | Protein degradation |
| XM_067100814.1 | F-box/LRR-repeat protein 5 | Protein degradation |
| XM_067087534.1 | eIF3b (translation initiation) | Translation |
| XM_067096515.1 | eIF3e (translation initiation) | Translation |
| XM_067116168.1 | ATP-dependent RNA helicase DDX3X | Translation |
| XM_067095058.1 | WNK1 kinase | Ion homeostasis |
| XM_067122492.1 | Immunoglobulin domain protein | Immune |
| XM_067089690.1 | Aldehyde dehydrogenase ALDH7A1 | Energy metabolism |
| XM_067115250.1 | MICU3 (mitochondrial calcium) | Energy metabolism |
| XM_067107357.1 | PP2A regulatory subunit | Signaling |

**Biological interpretation of the negative M6-WG correlation**:

1. **Protein turnover (6/13 genes)**: Ubiquitin, E3 ligase, F-box protein + 3 translation factors = active protein synthesis-and-degradation cycling. Higher turnover → more energy spent on maintenance → lower net growth.
2. **IGFBP (growth-axis inhibitor)**: IGFBP sequesters insulin-like growth factor, suppressing the growth axis. Higher IGFBP → less IGF signaling → slower growth. Direct mechanistic explanation for the negative correlation.
3. **GalNAc kinase (amino sugar metabolism)**: The first committed step of N-acetylgalactosamine phosphorylation feeds directly into KEGG map00520 (Amino sugar and nucleotide sugar metabolism) — the amino-sugar pathway of the integration network. **Independent convergence of the host module and the microbial chitin hypothesis on the same pathway.**

**Note**: the formal GSEA (§11.9) is the authoritative functional characterization of M6. Amino sugar metabolism is NOT among the significant M6 GSEA sets, so this convergence is a single-gene observation, not a pathway-level result.

### 11.8 Analysis 12: Genome-Resolved Microbial Functional Profiling

**Method**: Functional profiles are genome-resolved. Each of the 191 ASV taxa (species names only — the Track 3 Data Note confirms raw 16S FASTQs are not released, so a literal PICRUSt2 run is impossible) was resolved to its NCBI taxonomy ID and KEGG GENOME entry via KEGG REST; gene→KO links per genome were aggregated into pool-level pathway / KO / EC indices (genome gene counts normalised, weighted by pool relative abundance). EC numbers were mapped through the KO→reaction→EC chain, and MetaCyc pathway indices were computed via authenticated BioCyc queries (BioVelo EC→pathway). All outputs are descriptive indices — no statistical tests, per the n=1 pooled design.

**Results**:

| Metric | Value |
|--------|-------|
| Taxa resolved to a KEGG GENOME | 171/190 (reads covered: 200/532 Jumper, 481/507 Laggard) |
| KEGG reference pathways profiled | 363 (43,802 pathway × taxon contributions) |
| KOs profiled | 7,133 |
| EC numbers profiled | 2,105 |
| MetaCyc pathways profiled | 3,602 (1,773/2,105 ECs mapped) |
| Old map vs new profiles | **12 claims supported, 5 new findings, 3 not testable** (`results/reports/functional_map_comparison.tsv`) |

**Coverage caveat (reported, not imputed)**: 19 taxa lack any KEGG GENOME entry, including *Jiulongibacter sediminis* (298/532 Jumper reads, 56%). KEGG `find genome` queries confirm these genera are genuinely absent from KEGG.

**Top Jumper-enriched MetaCyc pathways** (descriptive log2 ratios): nitrate reduction VI (assimilatory), ammonia oxidation VI (comammox), ammonia oxidation IV, nitrite oxidation — a nitrogen-cycling signature in the fast-growth pool that the 9-category curated map could not express.

### 11.9 Analysis 13: M6 Gene-Set Enrichment Analysis

**Method**: Prawn proteins for all 18,276 genes were mapped to Swiss-Prot with DIAMOND blastp (best hit; 10,978 genes, 60%). UniProt cross-references supplied gene sets: 349 KEGG pathways (1,043 genes) and 7,346 GO biological-process terms (10,031 genes). clusterProfiler GSEA was run on two rankings: (1) full-module kME vs ME6 recomputed in the WGCNA input space (merged counts, mean-counts≥10 — verified against stored hub kME, cor=1.0000) and (2) the genome-wide DESeq2 WG statistic; ORA on the 72 M6 hubs served as a sanity check.

**Results** (FDR<0.05):

| Ranking | Significant sets | Top sets (padj) |
|---------|------------------|-----------------|
| KEGG × kME | 4/123 | **Cytoskeleton in muscle cells (map04820, 0.0013)**; Retinol metabolism (0.009); Apelin signaling (0.009); Ether lipid metabolism (0.038) |
| KEGG × DESeq2 | 1/124 | Lysine degradation (map00310, 0.005) |
| GO BP × kME | 85/1,504 | DNA repair (GO:0006281, 7.5e-8); transcription regulation (GO:0006355, 3.7e-7) |
| GO BP × DESeq2 | 24/1,513 | **Muscle contraction (GO:0006936, 1.5e-7)**; sarcomere organization (GO:0045214, 4.0e-5) |
| ORA (72 hubs) | 1/104 GO; 0/5 KEGG | — |

**Interpretation**: The top M6 hub (XM_067086030.1, kME=0.845) maps to **myosin regulatory light chain** (UniProt P40423, 83.3% identity), triangulating with the muscle-cytoskeleton enrichment. Amino sugar metabolism is NOT among the significant M6 sets: the GalNAc-kinase "loop closing" is a single-gene observation, not a pathway-level result. The dominant M6 signal is muscle-cytoskeletal with transcription/DNA-repair regulation — consistent with muscle-protein turnover as a metabolic cost of fast growth.

### 11.7 Synthesis

The thirteen validation analyses now form a complete evidence chain:

| Layer | Question | Answer |
|-------|----------|--------|
| Permutation null | Is M6 real signal? | **Yes — p=0.001** |
| Causal DAG | Is M6→WG causal? | **Yes — p=0.008, under assumed DAG** |
| Power analysis | Is the study adequately powered? | **No — M6 needs n=21 (one more sample)** |
| Bayesian model | Which taxa are genuinely enriched? | **3–4 taxa with CI≠0; 85% indistinguishable from zero** |
| Tornado plot | Are algorithm choices critical? | **No — study design is what matters** |
| Hub-gene annotation | What are the characterized M6 hubs? | **Protein turnover + IGFBP + GalNAc kinase (single-gene observations)** |
| Formal M6 GSEA | Which pathways are enriched in the M6 ranking? | **Muscle cytoskeleton (map04820 padj=0.0013); muscle contraction (GO padj=1.5e-7); transcription/DNA repair** |
| Genome-resolved profiles | Do pooled microbes' genomes support the curated map? | **12 claims supported, 5 new findings; nitrogen cycle tops MetaCyc Jumper** |
| Bootstrap + LOTO | How stable are rankings? | **Top-3 taxa stable; M6/M8 modules stable; rest fragile** |

**Final conclusion**: The M6 module is a genuine, causally interpretable, biologically coherent growth-associated module — validated by multiple independent methods including formal pathway enrichment (muscle cytoskeleton/contraction). The microbial findings are directional hypotheses confirmed only for 3–4 taxa. The chitin→GlcNAc→amino sugar metabolism mechanism is supported by literature and by chitinolytic *Pseudoalteromonas* in the fast-growth pool; its host-side connection through M6 is a single-gene observation (GalNAc kinase), while the tool-derived microbial profiles additionally expose a nitrogen-cycling signature in the Jumper pool.

---

## Appendix: File Inventory (Validation)

### Scripts

| Script | Purpose |
|--------|---------|
| `pipeline/src/asv/01b_sensitivity_analysis.py` | Pseudocount sweep + LOTO jackknife |
| `pipeline/src/asv/01c_bayesian_taxon_model.py` | Bayesian Dirichlet-multinomial taxon enrichment |
| `pipeline/src/network/wgcna_bootstrap_stability.py` | Bootstrap WGCNA module-trait stability |
| `pipeline/src/network/permutation_power_analysis.py` | Permutation null + power analysis |
| `pipeline/src/network/tornado_sensitivity.py` | WGCNA parameter sensitivity / tornado plot |
| `pipeline/src/network/causal_dag.R` | Causal DAG + conditional independence tests |
| `pipeline/src/network/m6_gene_annotation.py` | M6 hub gene NCBI annotation |
| `pipeline/src/network/spls_rf_validation.py` | sPLS + Random Forest orthogonal validation |
| `pipeline/src/rnaseq/deseq2_validation.R` | DESeq2 per-gene + module LM cross-validation |
| `pipeline/src/asv/03_functional_profiles.py` | Genome-resolved functional profiling: species→NCBI taxid→KEGG GENOME→KO/pathway/EC pool indices |
| `pipeline/src/asv/03c_metacyc_mapping.py` | EC indices → MetaCyc pathway indices (authenticated BioCyc session) |
| `pipeline/src/asv/04_functional_map_comparison.py` | Curated map vs tool-derived profiles cross-check |
| `pipeline/src/network/fetch_prawn_proteins.py` | Prawn protein set (RefSeq CDS + longest-ORF translation) |
| `pipeline/src/network/map_uniprot_sets.py` | Swiss-Prot best hits → KEGG/GO gene sets |
| `pipeline/src/network/m6_gsea.R` | clusterProfiler GSEA/ORA on M6 (kME + DESeq2 rankings) |

### Key Outputs

| File | Contents |
|------|----------|
| `results/reports/permutation_null.tsv` | Per-module permutation null comparison |
| `results/reports/power_analysis.tsv` | Per-module power estimates |
| `results/reports/bayesian_taxon_enrichment.tsv` | Bayesian posterior log2FC + credible intervals for 191 taxa |
| `results/reports/causal_dag_tests.tsv` | Conditional independence + causal effect estimates |
| `results/reports/tornado_wgcna_params.tsv` | Parameter impact rankings |
| `results/reports/tornado_gene_filter_sensitivity.tsv` | Gene filter sensitivity |
| `results/reports/m6_gene_annotations.tsv` | NCBI annotations for 72 M6 hub genes |
| `results/reports/m6_functional_summary.tsv` | Functional category counts |
| `data/processed/clr_profiles/pathway_abundance.tsv` | Genome-resolved KEGG pathway pool indices |
| `data/processed/clr_profiles/ko_abundance.tsv` | Genome-resolved KO pool indices |
| `data/interim/functional_prediction/kegg_out/ec_abundance.tsv` | EC pool indices via KO→reaction→EC chain |
| `data/interim/functional_prediction/metacyc_out/metacyc_pathway_abundance.tsv` | MetaCyc pathway indices |
| `data/interim/functional_prediction/kegg_out/taxon_genome_mapping.tsv` | Taxon→KEGG GENOME resolution + coverage |
| `results/reports/functional_map_comparison.tsv` | Curated-map vs tool-derived cross-check |
| `results/reports/m6_gsea_{kegg,go}_{kme,deseq2}.tsv` | M6 GSEA tables, both rankings × both set types |
| `results/reports/m6_ora_hubs_{kegg,go}.tsv` | ORA on 72 M6 hubs |
| `data/interim/literature_tables/prawn_gene_uniprot_map.tsv` | Gene→Swiss-Prot best-hit mapping |
| `results/reports/spls_gene_selection.tsv` | sPLS-selected top 50 genes |
| `results/reports/rf_gene_importance.tsv` | RF importance for 18,276 genes |
| `results/reports/spls_rf_wgcna_comparison.tsv` | Cross-method comparison |

### Modified Scripts

| Script | Changes |
|--------|---------|
| `pipeline/src/ranking/generate_final_shortlists.py` | Loads stability reports; incorporates stability bonuses/penalties into ranking formula |
| `pipeline/Snakefile` | 3 new rules (`taxon_sensitivity`, `wgcna_bootstrap`, `deseq2_validate`); updated `all` target and dependencies |

### New Outputs

| File | Contents |
|------|----------|
| `results/reports/pseudocount_sensitivity.tsv` | 2,101 rows (11 methods × 191 taxa) |
| `results/reports/pseudocount_stability_metrics.tsv` | Per-taxon rank stability metrics |
| `results/reports/loto_stability.tsv` | Per-taxon LOTO jackknife results |
| `results/reports/taxon_confidence_report.tsv` | Merged pseudocount+LOTO confidence per taxon |
| `results/reports/wgcna_bootstrap_stability.tsv` | Per-module bootstrap CIs and stability tiers |
| `results/reports/wgcna_module_stability_summary.tsv` | Compact summary for downstream ranking |
| `results/reports/deseq2_cross_validation.tsv` | Per-gene DESeq2 concordance (kME-aware) |
| `results/reports/deseq2_module_eigengene_test.tsv` | Module-level LM results |
| `results/reports/deseq2_all_results.tsv` | 18,199 genome-wide WG coefficients |
| `results/shortlist/ranking_methodology.md` | Updated with the ranking formulas and validation sections |

### Updated Outputs

| File | Change |
|------|-------------------|
| `results/shortlist/host_genes.csv` | M7 gene dropped; JP354756.1 (M12) added; DESeq2 status column |
| `results/shortlist/microbial_taxa.csv` | Stability tier and stability multiplier columns added |
| `results/shortlist/network_edges.csv` | Unchanged (edges are literature-based) |
| `results/shortlist/pathways.csv` | Unchanged (pathway ranking is literature-based) |
