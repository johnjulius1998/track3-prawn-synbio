# Ranking Methodology -- Track 3 Host-Microbe Integration (v3.1, n=20)

> **v5 correction (2026-08-17)**: the "Selected microbial taxa" and "Selected host genes" lists
> below are the v3.1 outputs and are **superseded** by `results/shortlist/*.csv` /
> `track3_shortlists.xlsx`. Two validation-driven gates, run after this document was written,
> remove some v3.1 selections:
>
> 1. **Bayesian Dirichlet-multinomial taxon enrichment** (`docs/VALIDATION_ANALYSIS.md` §11.3,
>    `src/asv/01c_bayesian_taxon_model.py`) replaces the pseudocount-dependent CLR ranking.
>    *Citrobacter koseri* and *Klebsiella variicola* (v3.1 ranks #3 and #5 below) are dropped:
>    both are LOTO-unstable (jackknife stability 0.000) and rank well outside the top tier by
>    Bayesian effect size. *Pseudoalteromonas phenolica* and *Endothiovibrio diazotrophicus*
>    replace them — both have a 95% credible interval excluding zero and a larger effect size,
>    but were absent from this v3.1 ranking because the CLR/LOTO framework flagged them as
>    unstable for the opposite reason (low raw counts, not point-estimate sensitivity).
>    **Final microbial shortlist: 4 taxa** (J. sediminis, M. osloensis, P. phenolica,
>    E. diazotrophicus) — one short of the 5-taxon ceiling, deliberately.
> 2. **Permutation null testing** (`results/reports/permutation_null.tsv`,
>    `src/network/permutation_power_analysis.py`) shows only M6 (empirical p=0.001) and M8
>    (p=0.014) exceed chance across 1,000 label shuffles; M7 (p=0.181) and M12 (p=0.335), whose
>    hub genes are v3.1 ranks #7-#10 below, are statistically indistinguishable from noise and
>    are dropped. **Final host-gene shortlist: 6 genes** (M6 x3, M8 x3) — four short of the
>    10-gene ceiling, deliberately (brief section 3a: "fewer, better-supported candidates is
>    not penalised").
>
> The pathway ranking (§3 below) is unaffected by either gate and stands as originally computed.
> See `results/shortlist/ranking_weights.csv` for the disclosed v5 ranking function, including
> the two gates above.

## 1. Host Gene Ranking

**Formula**: `Score = |partial_r(WG|sex)| × |kME| × growth_module_bonus × module_stability_bonus × deseq2_concordance_bonus`

| Component | Source | n | v3.1 Bonus |
|-----------|--------|---|------------|
| partial_r(WG\|sex) | Partial correlation of module eigengene with weight_gain, controlling for sex | 20 RNA-seq libraries | — |
| kME | Module membership: Pearson r(gene, module eigengene) | 20 samples | — |
| growth_module_bonus | 1.0 if module in top 5 by \|partial_r\|, else 0.5 | — | — |
| module_stability_bonus | Bootstrap stability tier (HIGH=1.0, MEDIUM=0.8, LOW=0.5) | 1000 bootstraps | v3 |
| deseq2_concordance_bonus | DESeq2 module eigengene LM test (AGREE_SIGNIFICANT=1.2, AGREE=1.0, WEAK_EFFECT=0.5, DISAGREE=0.5, STRONG_DISAGREE=0.3) | 20 samples, design=~sex+tissue+WG | v3.1 |

**Constraint**: Max 3 genes per module for biological diversity.

**Top 5 growth modules (n=20)**:
  - **M6**: partial r(WG|sex)=-0.5898, r(WG)=-0.1084, r(sex)=-0.3311
  - **M8**: partial r(WG|sex)=-0.4828, r(WG)=-0.1741, r(sex)=-0.1758
  - **M13**: partial r(WG|sex)=+0.3324, r(WG)=+0.3468, r(sex)=-0.1770
  - **M7**: partial r(WG|sex)=-0.3205, r(WG)=+0.1783, r(sex)=-0.4742
  - **M15**: partial r(WG|sex)=+0.3079, r(WG)=+0.3392, r(sex)=-0.1881

**Selected host genes**:
  - **#1 XM_067086030.1** (MM6): kME=+0.8454, score=0.478672, DESeq2=AGREE_SIGNIFICANT
  - **#2 XM_067082934.1** (MM6): kME=+0.7388, score=0.418314, DESeq2=AGREE_SIGNIFICANT
  - **#3 XM_067114945.1** (MM6): kME=+0.7127, score=0.403536, DESeq2=AGREE_SIGNIFICANT
  - **#4 XM_067127262.1** (MM8): kME=+0.6207, score=0.239739, DESeq2=AGREE
  - **#5 XM_067113377.1** (MM8): kME=+0.6161, score=0.237962, DESeq2=AGREE
  - **#6 GH624888.1** (MM8): kME=+0.6044, score=0.233443, DESeq2=AGREE
  - **#7 XM_067113580.1** (MM7): kME=-0.7367, score=0.059028, DESeq2=WEAK_EFFECT
  - **#8 JP354355.1** (MM7): kME=-0.6908, score=0.055350, DESeq2=WEAK_EFFECT
  - **#9 XM_067102939.1** (MM12): kME=+0.9169, score=0.055220, DESeq2=AGREE
  - **#10 JP354756.1** (MM12): kME=+0.8952, score=0.053913, DESeq2=AGREE

---

## 2. Microbial Taxon Ranking

**Formula**: `Score = |CLR_fold_diff| × (1 − contamination_score) × direction_multiplier × stability_multiplier`

| Component | Source | n | v3.1 |
|-----------|--------|---|------|
| CLR_fold_diff | CLR(Jumper) − CLR(Laggard) | 2 pooled samples | — |
| contamination_score | Log-normal mixture model on read counts (Fix 1) | 191 taxa | — |
| direction_multiplier | 1.0 (Jumper-associated), 0.7 (Laggard-associated) | — | — |
| stability_multiplier | Pseudocount+LOTO combined stability tier (HIGH=1.2, MEDIUM=1.0, LOW=0.7) | 11 CLR variants × 191 LOTO iterations | v3 |

**Exclusions**: HIGH contamination risk (score ≥0.85), "Unknown" taxon.

**CRITICAL**: All microbial associations are DIRECTIONAL HYPOTHESES ONLY.
No p-values, FDR, or statistical tests. n=1 pooled per group.

**Selected microbial taxa**:
  - **#1 Jiulongibacter sediminis**: CLR fold-diff=+3.1066, Jumper-enriched, contam_score=0.000, stability=MEDIUM, score=3.106600
  - **#2 Moraxella osloensis**: CLR fold-diff=-3.3126, Laggard-enriched, contam_score=0.008, stability=MEDIUM, score=2.300501
  - **#3 Citrobacter koseri**: CLR fold-diff=+2.1146, Jumper-enriched, contam_score=0.022, stability=MEDIUM, score=2.067444
  - **#4 Exiguobacterium aestuarii**: CLR fold-diff=-2.7530, Laggard-enriched, contam_score=0.013, stability=MEDIUM, score=1.901855
  - **#5 Klebsiella variicola**: CLR fold-diff=+1.5719, Jumper-enriched, contam_score=0.005, stability=MEDIUM, score=1.563726

---

## 3. Pathway Ranking

**Formula**: `Score = edge_count × |module_partial_r| × biological_relevance_factor`

| Component | Source | n |
|-----------|--------|---|
| edge_count | Number of microbial KEGG pathways with function-overlap edges to this pathway | KEGG GENOME profiles |
| module_partial_r | \|Partial r(WG\|sex)\| of strongest growth module | 20 samples |
| biological_relevance | 1.5 (chitin/amino-sugar), 1.0 (TCA/energy), 0.7 (other) | Expert curation |

**Selected pathways**:
  - **#1 PPAR signaling pathway**: 3 microbial edges, bio_relevance=1.0, score=1.769400
  - **#2 Amino sugar and nucleotide sugar metabolism**: 1 microbial edges, bio_relevance=1.5, score=0.884700
  - **#3 Lysosome**: 1 microbial edges, bio_relevance=1.5, score=0.884700

---

## 4. Validation Layers (v3.1)

### 4.1 Microbial Taxon Stability

| Analysis | Method | Key Finding |
|----------|--------|-------------|
| Pseudocount sensitivity | 11 CLR variants (pseudocount 0.1–50.0 + multiplicative + Bayesian) | 0/191 taxa pseudocount-stable; 179/191 flip direction |
| Leave-one-taxon-out (LOTO) | 191 jackknife iterations | Top 3 taxa LOTO-stable (1.000); C. koseri & K. variicola LOTO-unstable (0.000) |

### 4.2 Host Module Stability

| Analysis | Method | Key Finding |
|----------|--------|-------------|
| Bootstrap module-trait | 1000 resamples of n=20 with replacement | No module has 95% CI excluding zero; M6 top5_rate=0.748 |
| DESeq2 module eigengene LM | lm(ME ~ sex + tissue + WG), all 20 samples | M6: p=0.03 (significant); M7: p=0.81 (WG explains nothing) |
| DESeq2 per-gene cross-check | design=~sex+tissue+WG, continuous n=20 | 7/10 hub genes direction-concordant (kME-aware); 3/10 nominally significant |

### 4.3 Concordance Summary

| Layer | Concordance | Binomial p | Assessment |
|-------|-------------|------------|------------|
| Hub gene direction (v3.1 fixed) | 7/10 (70%) | p=0.17 | Not significant at α=0.05 with n=10 genes |
| Module eigengene direction | 13/20 (65%) | p=0.13 | Not significant at α=0.05 with n=20 modules |
| M6 (strongest module) | ✓ Concordant | LM p=0.03 | **Significant** — strongest validated finding |

---

## 5. What Each Edge Means

| Edge Type | Meaning |
|-----------|---------|
| `predicted_function_overlap` | KEGG GENOME-derived microbial pathway index → curated host pathway bridge. **NO cross-layer correlation.** |
| `phenotype_concordance` | Taxon more abundant in Jumper pool. **Directional hypothesis only.** |

## 6. Confound Handling

- **Sex-WG confound**: r(WG, sex) = −0.762 (females larger).
- **Post-hoc**: Partial correlation r(WG\|sex) on module eigengenes.
- **Pre-WGCNA (Fix 3)**: PCA removal of sex/tissue PCs → 66.2% variance removed.
- **v3.1**: DESeq2 design includes sex + tissue as covariates; module LM confirms sex/tissue separation.
