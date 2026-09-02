# Wet-Lab Validation Plan

## Track 3 Host-Microbe Integration — *Macrobrachium rosenbergii*

**Status**: Proposed experimental pipeline to validate the computational findings of `FINAL_REPORT.md` (Section 3.6) and `MANUSCRIPT.md` (Section 6.4). Designed to be executable at GK Aqua R&D — the assays below assume standard aquaculture-R&D-scale equipment (qPCR, a shared/outsourced GC-MS and LC-MS/MS service, standard microbiology culture, indoor juvenile-prawn rearing tanks) rather than specialised core-facility instrumentation, and are costed and scheduled accordingly (see Timeline).

---

## Overview

This tiered plan is designed to test the key computational findings of the study:

1. **Host hub genes** — M6/M8 co-expression module hub genes associated with weight gain.
2. **Ranked microbial taxa** — *Jiulongibacter sediminis*, *Moraxella osloensis*, *Pseudoalteromonas phenolica*, *Endothiovibrio diazotrophicus* (the 4 taxa whose 95% Bayesian credible interval excludes zero; see `results/shortlist/microbial_taxa.csv`. *Citrobacter koseri* and *Klebsiella variicola* were dropped from the ranked shortlist — both are LOTO-unstable point-estimate artifacts, see "Withdrawn from the plan" below).
3. **Core mechanism** — chitin degradation → GlcNAc → host amino sugar metabolism → growth.

The tiers are deliberately ordered to de-risk the strongest claims first (cheap, uses existing samples) before committing to animal trials and large cohorts.

---

## Tier 1 — Molecular Validation (qPCR)

Confirms the host hub genes and microbial taxa in independent measurements.

| Target | Method | Sample Requirements | Success Criteria |
|--------|--------|---------------------|------------------|
| Top 5 M6 hub genes (e.g., `XM_067086030.1`) | RT-qPCR (SYBR Green) on hepatopancreas cDNA | Existing RNA aliquots from all 20 individuals | \|ΔΔCt\| direction concordant with WGCNA sign; Spearman ρ ≥ 0.6 between qPCR fold-change and RNA-seq TPM |
| *Jiulongibacter sediminis* | Genus-specific 16S qPCR on gut DNA | Newly collected individual gut contents (n≥30) | Detection in ≥50% of samples; correlation with weight gain |
| *Pseudoalteromonas phenolica* | Species-specific 16S qPCR | Same gut DNA as above | Detection in ≥50% of samples; correlation with host chitinase gene expression |
| *Moraxella osloensis* | Species-specific 16S qPCR | Same gut DNA as above | Confirmation of presence; abundance vs weight gain (Laggard-side confirmation arm) |
| *Endothiovibrio diazotrophicus* | Species-specific 16S qPCR | Same gut DNA as above | Confirmation of presence; abundance vs weight gain. Largest Bayesian effect size of the 4 confirmed taxa (log2FC +4.69) but functionally uncharacterised (no KEGG GENOME match) — this assay is confirmatory only, not mechanism-testing |

---

## Tier 2 — Functional Validation (Metabolomics)

Tests the biological mechanism (chitin → GlcNAc → SCFA) at the metabolite level.

| Target | Method | Sample Requirements | Success Criteria |
|--------|--------|---------------------|------------------|
| SCFAs (acetate, butyrate, propionate) | GC-MS on hepatopancreas homogenate | 20 individual hepatopancreas samples (existing tissue) | Significantly higher SCFA levels in high-WG individuals |
| N-acetylglucosamine (GlcNAc) | LC-MS/MS targeted assay | Same samples | Higher GlcNAc in individuals with detectable *Pseudoalteromonas* |
| Chitin degradation products | Chitinase activity assay + GlcNAc quantification | Gut content + hepatopancreas | Positive correlation between chitinolytic bacteria and GlcNAc |

---

## Tier 3 — Causal Validation (In Vivo / In Vitro)

Establishes causality, not just association.

| Experiment | Design | Expected Outcome |
|-----------|--------|------------------|
| **Probiotic supplementation trial** | n=40 juvenile prawns, randomised into 4 groups of 10 (basal feed control; *Pseudoalteromonas phenolica*-supplemented feed, low dose; high dose; high dose + chitin-enriched feed); 8-week growth trial with weekly weighing | Improved weight gain and FCR in supplemented groups relative to control, dose-dependent; increased hepatopancreatic GlcNAc in the chitin-enriched arm. (*Citrobacter*/*Klebsiella* supplement arms from the earlier v3.1 design are dropped — both genera were withdrawn from the ranked shortlist after validation, see "Withdrawn from the plan" below) |
| **Chitin digestibility assay** | In vitro: incubate sterile chitin with *Pseudoalteromonas phenolica* isolates (+ uninoculated chitin as negative control); measure GlcNAc release over time | Dose-dependent GlcNAc release confirming chitinolytic capability, absent in the uninoculated control |
| **SCFA production assay** | Anaerobic culture of *Jiulongibacter sediminis* and *P. phenolica* isolates with prawn feed substrate (+ uninoculated substrate as negative control); measure SCFA by GC-MS | Detectable acetate, propionate, butyrate production above the uninoculated baseline |
| **Individual-level 16S sequencing** | Collect gut contents from n≥30 individual prawns with measured WG; sequence V3–V4 16S | Replace n=1 pooled design; enable Spearman correlation between taxon abundance and WG |

---

## Tier 4 — Independent Replication Cohort

A fully independent cohort of **n≥60** *M. rosenbergii* with:

- Individual 16S gut microbiome (not pooled)
- Individual RNA-seq (hepatopancreas + gut epithelium)
- Individual weight gain measurements
- Both sexes represented
- Multiple time points through the molt cycle

This would enable:

- (a) statistical testing of all microbial associations,
- (b) direct cross-layer correlation (taxon abundance ↔ host gene expression),
- (c) sex-stratified networks,
- (d) longitudinal tracking of microbiome–transcriptome dynamics through molting.

---

## Key Hypotheses Being Tested

1. **M6/M8 hub genes** are genuinely growth-associated (negative partial correlation with weight gain), not artifacts of sex confounding or small-sample WGCNA.
2. **Jumper-associated taxa** (especially *Jiulongibacter sediminis*, *Pseudoalteromonas* spp.) are real gut colonizers, not environmental/reagent contaminants.
3. **Chitin degradation → GlcNAc → amino sugar metabolism** is a functional pathway linking gut microbes to host growth.

---

## Validation-Informed Priorities (2026-08-16)

The internal validation stack (permutation null, Bayesian modeling, causal DAG, functional annotation, **formal M6 GSEA**, and **genome-resolved microbial functional profiling** — see `docs/VALIDATION_ANALYSIS.md` §11) refines what the wet-lab should test first:

### Priority 1 — Directly test the now-validated M6 biology

NCBI annotation of M6 hub genes plus formal KEGG/GO enrichment identified specific, testable molecular targets:

| Target | Gene | qPCR Assay | Hypothesis |
|--------|------|-----------|------------|
| **IGFBP** | XM_067106657.1 | RT-qPCR on hepatopancreas cDNA (existing n=20) | Higher IGFBP → lower IGF signaling → lower WG. Confirm expression inversely correlates with WG |
| **GalNAc kinase** | XM_067104953.1 | RT-qPCR + enzyme activity assay | M6's amino-sugar enzyme links the module to the chitin→GlcNAc pathway. Confirm activity correlates with WG |
| Ubiquitin / E3 ligase / F-box | EL609362.1, XM_067109294.1, XM_067100814.1 | RT-qPCR panel | Confirm the "protein turnover" signature of M6 — higher expression in low-WG animals |
| eIF3b/e | XM_067087534.1, XM_067096515.1 | RT-qPCR | Translation machinery confirming turnover hypothesis |
| **Myosin regulatory light chain (top M6 hub)** | XM_067086030.1 | RT-qPCR + Western blot | **GSEA-derived**: muscle cytoskeleton (map04820, padj=0.0013) and muscle contraction (GO:0006936, padj=1.5e-7) are the module's top enriched sets. Confirm protein abundance inversely correlates with WG |
| Muscle/cytoskeletal panel | top hub + sarcomere-set members | RT-qPCR panel | Formal enrichment identifies muscle/cytoskeletal biology — assay the sarcomere organization genes |

**GSEA caveat**: amino sugar metabolism is NOT among the significant M6 GSEA sets. The GalNAc-kinase target above is therefore a single-gene (exploratory) assay, not a pathway-level prediction.

### Priority 2 — One additional RNA-seq sample to confirm M6

Power analysis shows M6 needs n=21 for 80% power (currently n=20). **A single additional hepatopancreas RNA-seq sample with measured weight gain would push M6 to the confirmation threshold.** This is the highest-value, lowest-cost experiment available.

### Priority 3 — Individual-level 16S for the validated top-tier taxa

The Bayesian model confirmed 4 taxa with 95% credible intervals excluding zero (29/191 taxa overall have any nonzero CI; these 4 are the largest-effect, LOW-contamination-risk subset — see `results/shortlist/microbial_taxa.csv`). The individual-level 16S cohort (Tier 3) should prioritize, in order of how directly each links to a testable mechanism:
- *Pseudoalteromonas phenolica* (log2FC +4.13, CI [+1.28, +8.62]) — carries the chitin-degradation mechanism (Priority 1 experimentally); genome-resolved.
- *Jiulongibacter sediminis* (log2FC +4.02, CI [+3.38, +4.73]) — the single dominant Jumper taxon (55.5% of Jumper reads); tightest CI of the four; unresolved to a KEGG GENOME entry, so confirmatory only.
- *Moraxella osloensis* (log2FC −4.63, CI [−7.31, −2.67]) — the reciprocal Laggard-side signal; genome-resolved, contributes to the same PPAR/amino-sugar pool indices from the opposite direction.
- *Endothiovibrio diazotrophicus* (log2FC +4.69, CI [+1.88, +9.12]) — largest point estimate of the four but the widest CI (Jumper-exclusive, 0 Laggard reads) and unresolved to a KEGG GENOME entry; confirmatory only, no functional assay currently attaches to it.

*Citrobacter koseri* and *Klebsiella variicola* are **deprioritized** — both failed LOTO stability (0.00) and, while their Bayesian CIs also exclude zero, their effect sizes (log2FC +2.35 and +1.78 respectively) rank well outside this top tier.

### Priority 4 — Functional profiling-derived hypotheses

The tool-derived microbial functional profiles (KEGG GENOME + MetaCyc; descriptive, n=1 pools) generate two testable functional hypotheses:

| Hypothesis | Evidence | Assay |
|------------|----------|-------|
| **SCFA → PPAR signaling bridge (rank #1)** | Microbial butanoate/propanoate metabolism and fatty-acid degradation index → host PPAR | GC-MS SCFA quantification (acetate/butyrate/propionate) in gut content + PPAR target-gene qPCR in hepatopancreas |
| **Nitrogen-cycle signature in the Jumper pool** | Top MetaCyc Jumper-enriched pathways: nitrate reduction VI (assimilatory), ammonia oxidation (comammox/IV), nitrite oxidation | Nitrate/ammonium/nitrite assays on gut content; functional-gene qPCR (amoA, nirK, narG) on pooled DNA |

Both are reported descriptively (n=1 pools); the assays above would test them at individual level.

### Withdrawn from the plan

| Item | Reason |
|------|--------|
| M7 hub genes (XM_067113580.1, JP354355.1, XM_067118937.1) | DESeq2 LM shows WG explains nothing (p=0.81); permutation null p=0.18 |
| *Citrobacter koseri*, *Klebsiella variicola* 16S qPCR | LOTO-unstable (CLR ranking was a pseudocount artifact); both do have a nonzero Bayesian CI but rank well below the top-4 tier by effect size (log2FC +2.35 and +1.78 vs +4.0-4.7 for the confirmed tier) |
| Probiotic trial arm with *Citrobacter* supplement | Genus no longer ranked after validation |

### Retained unchanged

The chitin metabolomics (GlcNAc LC-MS/MS, chitinase activity) and the *Pseudoalteromonas* probiotic/chitin-digestibility experiments remain unchanged — the chitin mechanism is now *stronger* than before, reinforced by the GalNAc kinase finding inside M6.

---

## Timeline

| Tier | Scope | Estimated Duration | Sample Requirements |
|------|-------|--------------------|---------------------|
| 1 | Molecular (qPCR) | 3–6 months | Uses existing samples |
| 2 | Functional metabolomics | 3–6 months | Uses existing samples |
| 3 | Causal (animal trial) | 6–12 months | Requires animal trial |
| 4 | Replication cohort | 12–18 months | New cohort + sequencing budget |

> **Caveat**: All microbial associations in the computational analysis are directional hypotheses only (n=1 pooled per growth group). The wet-lab plan, especially Tier 3 (individual-level 16S) and Tier 4 (replication), is what converts these hypotheses into statistically testable associations. The validation stack has narrowed the focus to M6 (host), the top-3 microbial tier, and the chitin→amino sugar axis.
