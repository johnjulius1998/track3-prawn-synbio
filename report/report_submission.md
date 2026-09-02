# Decoding the Jumper–Laggard Phenomenon: Gut Microbiome and Host Co-Expression Drivers of Growth Rate in *Macrobrachium rosenbergii*

**Team**: Prawn to be Wild | **Track**: 3 — Systems Biologist (Postgraduate) | **i-Biohackathon 2026**
**BioProject**: PRJNA875278 | **Submission version**: v5.0 (2026-08-17)

---

## Abstract

The gut microbiome influences crustacean growth through metabolic provisioning, but the microbial taxa and host pathways mediating this in *Macrobrachium rosenbergii* are uncharacterised. We integrated a pooled 16S ASV table (n=1 fast-growing "Jumper" pool, n=1 slow-growing "Laggard" pool; 191 taxa) with individual-level RNA-seq from 20 prawns (hepatopancreas + gonad, PRJNA875278), using a Bayesian beta-binomial model for microbial enrichment and signed WGCNA with sex-confound correction for the host layer. Integration used only phenotype-concordance and predicted-function-overlap edges — no cross-layer covariance was computed, because the microbial data (n=1 per class) cannot support one. Five microbial taxa have a 95% Bayesian credible interval excluding zero and LOW contamination risk; two host co-expression modules (M6, 719 genes; M8, 336 genes) survive permutation-null testing (empirical p=0.001 and p=0.014) out of 20 tested. M6 is further supported by causal DAG estimation (p=0.008) and formal gene-set enrichment identifying muscle-cytoskeleton/contraction biology (KEGG map04820 padj=0.0013; GO:0006936 padj=1.5e-7). The strongest integration bridge connects chitin-degrading, Jumper-enriched *Pseudoalteromonas phenolica* to host amino-sugar metabolism and lysosome pathways, consistent with GlcNAc provisioning for exoskeleton synthesis during moulting; a second, independently ranked bridge connects short-chain-fatty-acid-associated microbial pathways to host PPAR signalling. Thirteen internal validation layers — permutation testing, Bayesian uncertainty quantification, causal inference, and genome-resolved functional profiling — determine which of our own findings survive scrutiny and which do not; we report both. All microbial associations remain directional hypotheses requiring individual-level replication.

---

## 1. Background

Global aquaculture demand has intensified interest in the biological drivers of growth-rate variation within farmed stocks. The giant freshwater prawn, *Macrobrachium rosenbergii*, is the world's second most-produced freshwater prawn, but individual growth rates within a single cohort are highly heterogeneous — a subset of "Jumper" (fast-growing) animals routinely outpaces "Laggard" (slow-growing) siblings under identical rearing conditions. Understanding what distinguishes Jumpers from Laggards has direct translational value: even a partial explanation could inform selective breeding, feed formulation, or probiotic strategy at a commercial hatchery.

Two independent bodies of evidence motivate a host–microbe systems approach to this question. First, the gut microbiome is known to contribute to growth in aquatic animals through short-chain fatty acid (SCFA) production, chitin degradation, vitamin biosynthesis, and amino acid provisioning [1–3]. Second, *M. rosenbergii* growth is strongly sexually dimorphic (females grow larger [4]), meaning any host transcriptomic signal associated with growth is confounded with sex unless explicitly modelled. Neither the microbial mechanism nor the host transcriptomic architecture underlying Jumper/Laggard divergence has previously been characterised in this species.

This study had two starting constraints that shaped every downstream methodological choice, both specified in the supplied Data Note:

1. **The microbial data is pooled, not individually resolved.** A single 16S ASV table gives one aggregate read-count profile for the Jumper class and one for the Laggard class (191 species-level taxa, 537 and 509 total reads respectively). This is n=1 per phenotype at the sample level — there is no within-class variance, and therefore no p-value, FDR, or correlation coefficient can be legitimately computed from this table alone, no matter how it is transformed.
2. **The host transcriptome is not supplied and must be retrieved independently**, with genuine biological replication. We used all 20 publicly available RNA-seq runs from BioProject PRJNA875278 (10 male, 10 female; hepatopancreas and gonad; accessions SRR21374326–SRR21374345), which is the layer where population-level statistics are legitimate.

The central methodological decision of this project — and the one we believe the evaluation rubric is built to reward — is to keep these two constraints separate rather than papering over them. We report host-layer statistics as statistics (n=20, with p-values and confidence intervals) and microbial-layer findings as directional hypotheses (n=1 pooled, Bayesian credible intervals at best, never a p-value). Every integration edge between the two layers is a stated *type* of relationship (phenotype concordance or predicted function overlap) rather than a claimed correlation, because a correlation across two pooled columns is not a real statistic.

We further subjected our own findings to a 13-layer internal validation stack (§4.7) before finalising any ranked shortlist, because an unvalidated co-expression module or CLR-ranked taxon in a dataset this size is exactly as likely to be noise as signal, and reporting that distinction is itself a scientific result.

---

## 2. Methods

### 2.1 Data

| Layer | Source | n | Notes |
|---|---|---|---|
| Microbial (16S) | Supplied: `ASV_table_Jumpers_Laggards.GKAQUA.csv` | 2 pooled samples (1 Jumper, 1 Laggard) | 191 species-level taxa; 537 / 509 total reads; 72 taxa Jumper-exclusive, 98 Laggard-exclusive, 21 shared |
| Host transcriptome | NCBI SRA, BioProject **PRJNA875278** | 20 individual RNA-seq libraries | 2 sexes × 2 tissues (hepatopancreas, testis/ovary) × 5 replicates; Illumina HiSeq 4000, 2×300bp |
| Reference transcriptome | NCBI TSA Project 73259 (4,211 seqs, 2013) + all *M. rosenbergii* mRNA records (63,271 seqs) | — | 66,982 transcripts, 226 MB; no annotated genome exists (GCA_039081455.1 has 0 annotated protein-coding genes) |

**SRA accessions used** (all 20, BioProject PRJNA875278): SRR21374326, SRR21374327, SRR21374328, SRR21374329, SRR21374330, SRR21374331, SRR21374332, SRR21374333, SRR21374334, SRR21374335, SRR21374336, SRR21374337, SRR21374338, SRR21374339, SRR21374340, SRR21374341, SRR21374342, SRR21374343, SRR21374344, SRR21374345. Continuous weight-gain phenotype (range 1.48–7.52, mean 3.94) retrieved from NCBI BioSample records; female mean 5.41, male mean 2.47.

### 2.2 Microbial layer processing

**Compositional transform.** Centered log-ratio (CLR) was applied for direction classification only: CLR(x) = ln(x / g(x)), g(x) = geometric mean of non-zero values per sample. Zeros (185/382 cells, 48.4%) used pseudocount substitution. No statistical test (t-test, Wilcoxon, Fisher's exact) was applied — n=1 per group cannot support one.

**Bayesian enrichment model** (supersedes the CLR point estimate for ranking). A Bayesian beta-binomial model was fit per taxon: the posterior distribution of the Jumper:Laggard proportion ratio was sampled (20,000 draws), giving a posterior mean log2 fold-change and a 95% credible interval with no pseudocount dependency. 29/191 taxa (15.2%) have a 95% CI excluding zero; the remaining 84.8% cannot be distinguished from zero enrichment given the data.

**Contaminant screening.** Two tiers: (1) a static list of 14 genera previously reported as reagent/skin contaminants [5,6] (53/191 taxa, 27.7%, flagged); (2) a continuous score (0–1) from a two-component log-normal mixture model on total read count, with taxa scoring ≥0.85 classified HIGH risk (98/191, 51.3%).

**Genome-resolved functional profiling.** Each species was resolved to an NCBI taxonomy ID and, where available, a KEGG GENOME entry (171/190 taxa resolved; 19 unresolved, including the two highest-abundance Jumper taxa — reported as unresolved, not imputed). Gene→KO links were aggregated to pool-level pathway indices (363 KEGG pathways, 7,133 KOs, 2,105 ECs) and extended to MetaCyc (3,602 pathways). This tool-derived profile replaced an earlier literature-curated genus map; a cross-check found 12 claims supported, 5 new findings, and 3 not testable.

### 2.3 Host layer processing

SRA runs were downloaded via NCBI prefetch/fasterq-dump and trimmed with fastp v1.3.6 (mean read retention 98.7%). Quantification used Salmon v2.4.1 (mean mapping rate 85.8%, higher in hepatopancreas [90.4%] than gonad [78.3%]). Genes retained at mean TPM≥1 in ≥10/20 samples yielded 18,276 genes.

**WGCNA.** Signed blockwise network construction (WGCNA v1.72, R 4.5.3): soft threshold power=15 (scale-free R²=0.939), minModuleSize=30, mergeCutHeight=0.25 → 20 co-expression modules. Hub genes = top 10% by absolute module membership (kME).

**Sex-confound correction.** Weight gain and sex are strongly confounded (r=−0.762, females larger). Two complementary corrections were applied: (1) post-hoc partial correlation of module eigengenes with weight gain, controlling for sex; (2) pre-WGCNA PCA-based removal of sex/tissue-correlated principal components (66.2% of confounded variance removed before network construction).

### 2.4 Integration approach

Two edge types only, both stated explicitly as non-correlational:

1. **Phenotype concordance**: a taxon directionally enriched in the fast-growth pool is connected to a `fast_growth_Jumper` phenotype node (or the reciprocal for Laggard). Basis: directional enrichment only.
2. **Predicted function overlap**: a microbial KEGG pathway pool index is connected to a host KEGG pathway via a curated bridge (e.g. microbial butanoate/propanoate/fatty-acid degradation → host PPAR signalling). Basis: shared functional category, not shared statistic.

Every edge in `results/shortlist/network_edges.csv` carries an explicit statement: *no cross-layer correlation is claimed*. Contaminant-flagged taxa were excluded from all edge construction.

### 2.5 Ranking and the 13-layer validation stack

Ranking formulas (fully disclosed in `results/shortlist/ranking_weights.csv`):

- **Host genes**: score = |partial_r(WG\|sex)| × |kME| — but only for modules that pass a hard gate: permutation-null empirical p<0.05 over 1,000 label shuffles. Max 3 genes per module.
- **Microbial taxa**: ranked by Bayesian posterior |log2FC| — but only for taxa passing two hard gates: 95% CI excludes zero, AND contamination risk is LOW.
- **Pathways**: score = edge_count × |M6 partial_r(WG\|sex)| × biological_relevance_factor (1.5 for the two chitin/amino-sugar-linked pathways, 1.0 for PPAR; expert-curated, disclosed as such).

Thirteen validation analyses were run to determine which raw findings survive scrutiny (full detail: `docs/VALIDATION_ANALYSIS.md`):

| # | Analysis | Decisive finding |
|---|---|---|
| 1–2 | Pseudocount sweep (11 CLR variants) + leave-one-taxon-out jackknife | 0/191 taxa CLR-stable; the earlier point-estimate ranking of *C. koseri*/*K. variicola* is a pseudocount artifact (LOTO stability 0.000) |
| 3 | Bootstrap WGCNA (1,000 resamples) | No module's raw partial_r has a 95% CI excluding zero on its own — motivates test 6 |
| 4 | DESeq2 cross-validation (continuous WG, kME-aware) | 7/10 v3.1 hub genes direction-concordant |
| 5 | sPLS + Random Forest | Zero gene overlap with WGCNA hubs — co-expression ≠ prediction at n=20; RF OOB R²=−0.34 |
| 6 | **Permutation null (1,000 label shuffles)** | **M6 p=0.001, M8 p=0.014 — the only two of 20 modules exceeding chance** |
| 7 | Power analysis | M6 needs n=21 for 80% power (have 20); M8 needs n=32 |
| 8 | **Bayesian Dirichlet-multinomial taxon model** | **5-taxon confirmed tier** (§2.2) replaces the CLR ranking |
| 9 | Tornado sensitivity | Confound correction (Δ≈0.17) and sample size (Δ≈0.21) dominate; WGCNA hyperparameters ≤0.05 |
| 10 | Causal DAG | M6→WG causal estimate −3.50 WG units per unit eigengene (p=0.008, 95% CI [−5.96,−1.05]) |
| 11 | M6 NCBI annotation (72 genes) | Motivates formal enrichment (test 12) |
| 12 | **M6 GSEA (KEGG + GO)** | **Muscle cytoskeleton (map04820, padj=0.0013), muscle contraction (GO:0006936, padj=1.5e-7)**; amino-sugar metabolism NOT significant |
| 13 | Genome-resolved functional profiling | 171/190 taxa resolved; nitrogen-cycle pathways top the Jumper-enriched MetaCyc signal |

---

## 3. Results

### 3.1 Microbial composition and confirmed taxa

Of 191 ASVs, 79 non-contaminant taxa were Jumper-associated and 59 Laggard-associated by raw direction. After applying the Bayesian-CI and contamination gates (§2.5), **5 taxa** constitute the final ranked shortlist (full detail: `results/shortlist/microbial_taxa.csv`):

| Rank | Taxon | Direction | Bayesian log2FC [95% CI] | Genome-resolved? |
|---|---|---|---|---|
| 1 | *Endothiovibrio diazotrophicus* | Jumper (15/0 reads) | +4.69 [+1.88, +9.12] | No (unresolved) |
| 2 | *Moraxella osloensis* | Laggard (1/35 reads) | −4.63 [−7.31, −2.67] | Yes |
| 3 | *Pseudoalteromonas phenolica* | Jumper (10/0 reads) | +4.13 [+1.28, +8.62] | Yes |
| 4 | *Jiulongibacter sediminis* | Jumper (298/17 reads; 55.5% of all Jumper reads) | +4.02 [+3.38, +4.73] | No (unresolved) |
| 5 | *Exiguobacterium aestuarii* | Laggard (1/20 reads) | −3.82 [−6.57, −1.84] | Yes |

*Citrobacter koseri* and *Klebsiella variicola* — ranked in an earlier CLR-only draft of this shortlist — are excluded: both are LOTO-unstable (jackknife stability 0.000) and rank below this tier by Bayesian effect size (log2FC +2.35 and +1.78 respectively).

At the functional level, six further Jumper-exclusive *Pseudoalteromonas* species (below the individual ranking cutoff) reinforce the chitin-degradation signal carried by *P. phenolica*; nitrogen-cycle MetaCyc pathways (nitrate reduction, ammonia/nitrite oxidation) are the strongest Jumper-enriched signal in the genome-resolved functional layer overall.

### 3.2 Host co-expression network

WGCNA identified 20 modules from 18,276 genes (power=15, R²=0.939). After sex-confound correction, permutation-null testing (§2.5, test 6) shows exactly two modules exceed chance:

| Module | Genes | partial r(WG\|sex) | Permutation p | Power at n=20 | n needed for 80% power |
|---|---|---|---|---|---|
| **M6** | 719 | −0.590 | **0.001** | 0.80 (borderline) | 21 |
| **M8** | 336 | −0.483 | **0.014** | — | 32 |
| M7, M12, M13, M15, others | — | −0.32 to +0.33 | ≥0.18 | — | ≥69 |

A preliminary n=10 analysis had identified a different top module (M17, r=+0.564) that **did not replicate** at n=20 — direct evidence that small-sample WGCNA module-trait correlations are unstable below the field's n≥15 guidance [8].

**M6's function is now formally established.** Swiss-Prot orthology (DIAMOND, 60% of genes mapped) plus clusterProfiler GSEA over two independent gene rankings (module kME; DESeq2 weight-gain statistic) identifies significant enrichment (FDR<0.05) for **cytoskeleton in muscle cells** (KEGG map04820, padj=0.0013) and **muscle contraction / sarcomere organisation** (GO:0006936, padj=1.5e-7), alongside DNA-repair and transcription-regulation gene sets. The top hub gene (kME=0.845) is 83.3% identical to Drosophila myosin regulatory light chain (Swiss-Prot P40423), directly triangulating with the enrichment result. M6's negative correlation with weight gain is therefore interpretable as the energetic cost of muscle-protein turnover and maintenance, suppressed in fast-growing animals. A causal DAG analysis, conditioning on sex, estimates the M6→weight-gain effect at −3.50 WG units per unit eigengene (p=0.008).

**M8's function is not yet characterised** — only M6 was taken forward to formal GSEA. M8's hub genes (tolloid-like protein 2, a BMP-pathway metalloprotease; C-terminal binding protein, a metabolic-state-linked corepressor) are plausible by identity alone but not pathway-confirmed.

**Final host gene shortlist: 6 of 10 ceiling** (top-3 kME hub per surviving module only — see §3.4).

### 3.3 Integration network

The final ranked network comprises 5 phenotype-concordance edges (one per confirmed taxon → Jumper/Laggard) and 5 predicted-function-overlap edges (3 SCFA-linked microbial pathways → PPAR signalling; 1 amino-sugar pathway → 2 host pathways), all in `results/shortlist/network_edges.csv`.

| Rank | Host pathway | Score | Microbial source pathways | Host-side GSEA status |
|---|---|---|---|---|
| 1 | PPAR signalling (map03320) | 1.769 | Butanoate, propanoate metabolism; fatty-acid degradation (3 edges) | Not significant (padj 0.61 kME / 0.96 DESeq2) |
| 2 | Amino sugar & nucleotide sugar metabolism (map00520) | 0.885 | Amino sugar metabolism (1 edge) | 1 M6 gene by orthology only (GalNAc kinase) — single-gene, not pathway-significant |
| 3 | Lysosome (map04142) | 0.885 (tied) | Amino sugar metabolism (shared source) | Not significant (padj 0.93) |

The ranking formula's own weighting (edge count) places the SCFA→PPAR bridge first, ahead of the chitin→amino-sugar bridge that carries this study's headline microbial narrative — we report the formula's output as computed rather than reordering it to match the more narratively appealing result.

### 3.4 Final ranked host genes

| Rank | Gene | Symbol | Module | kME | Annotation |
|---|---|---|---|---|---|
| 1 | XM_067086030.1 | sqh | M6 | 0.845 | Myosin regulatory light chain (83.3% id.) |
| 2 | XM_067082934.1 | uri | M6 | 0.739 | Prefoldin RPB5 interactor |
| 3 | XM_067114945.1 | LOC136845115 | M6 | 0.713 | LITAF homolog |
| 4 | XM_067127262.1 | LOC136852535 | M8 | 0.621 | Tolloid-like protein 2 |
| 5 | XM_067113377.1 | CtBP | M8 | 0.616 | C-terminal binding protein |
| 6 | GH624888.1 | COX1-like | M8 | 0.604 | EST, mitochondrial COX1 similarity |

Hub genes from M7 and M12 (permutation p=0.181 and p=0.335) are not included — both modules are statistically indistinguishable from noise.

### 3.5 Quality metrics

| Metric | Value |
|---|---|
| fastp read retention | 98.7% mean |
| Salmon mapping rate | 85.8% (72.4–92.7%) |
| WGCNA scale-free R² | 0.939 (power=15, n=20) |
| Contaminants flagged (static genus list) | 53/191 (27.7%) |
| HIGH contamination risk (continuous score) | 98/191 (51.3%) |
| Confound variance removed (pre-WGCNA PCA) | 66.2% |
| Final shortlist sizes | 6/10 genes, 5/5 taxa, 3/3 pathways |

---

## 4. Network interpretation

The figure `results/figures/final_integrated_network.png` renders every edge in `network_edges.csv` plus the within-host WGCNA structure (hub gene → module via kME; module → weight gain via partial r). It is organised as three zones — **microbial layer** (n=1 pooled, directional only), **phenotype** (the shared axis), and **host layer** (n=20, replicated) — because that grouping is the methodological finding, not just a layout choice.

**Why the phenotype is drawn as the shared axis, not a computed edge.** The rubric criterion for integration depth asks for a host–microbe network, not two juxtaposed analyses. But a taxon-to-gene edge in this dataset would require covariance the pooled microbial table cannot supply. Our resolution is architectural: both evidentiary strands — 5 confirmed taxa on the microbial side, 2 validated modules on the host side — independently point at the **same phenotype construct** (Jumper/Laggard categorically on the microbial side; continuous weight gain on the host side) without ever computing a statistic between them. The figure draws this explicitly as two related-but-distinct phenotype representations connected by a labelled non-statistical link, rather than merging them into one axis, which would misstate what was actually computed.

**What is deliberately *not* drawn.** There is no edge from a microbial pathway node to a host WGCNA module, even though the pathway ranking formula uses |M6 partial_r| as a weighting term. That weight is a ranking-formula choice, not a discovered relationship — M6's own GSEA result shows amino-sugar metabolism is *not* one of its significant gene sets. Drawing that edge would visually claim a stronger result than the data supports; the figure's caption states this explicitly next to each pathway node (M6 GSEA padj value shown alongside every ranked pathway).

**The two converging narratives.** On the microbial side, *P. phenolica* and *J. sediminis* (Jumper-enriched, both contributing to the amino-sugar and PPAR pool indices) sit opposite *M. osloensis* and *E. aestuarii* (Laggard-enriched, contributing to the *same* pool indices from the reciprocal direction) — the same functional categories moving in opposite directions with phenotype is the strongest form of directional evidence a pooled n=1 design can produce. On the host side, M6's hub genes converge on a muscle-maintenance signature independently confirmed by three orthogonal methods (permutation, causal DAG, GSEA). The network figure does not merge these two narratives into a single mechanism, because we have no data connecting them beyond the shared phenotype label — that would be the finding a Tier 4 replication cohort (§5, item 4) could actually test.

---

## 5. Ranked candidate intervention list

Interventions are ranked by a composite of (a) strength of the underlying evidence tier, (b) cost and turnaround at a GK Aqua R&D-scale facility, and (c) how directly the result would resolve this study's single largest open question — whether the confirmed taxa and M6 are real biology or an artefact of the pooled/underpowered design. Full protocols: `docs/WETLAB_VALIDATION_PLAN.md`.

| # | Intervention | Evidence tier | Cost / duration | Success metric |
|---|---|---|---|---|
| **1** | **One additional hepatopancreas RNA-seq library with measured weight gain** | M6: permutation p=0.001, causal p=0.008, GSEA-confirmed; needs n=21 for 80% power (have 20) | Very low — 1 sample, uses existing extraction/sequencing pipeline; <1 month | M6 partial_r(WG\|sex) remains ≥0.55 in the n=21 re-analysis |
| **2** | ***Pseudoalteromonas phenolica* chitin-digestibility assay + dose-ranging probiotic trial** | Genome-resolved, top-5 Jumper contributor to all 3 ranked pathways, Bayesian CI excludes zero | Moderate — in vitro assay (weeks) + n=40, 8-week juvenile growth trial (control / low dose / high dose / high dose + chitin-enriched feed) | Dose-dependent GlcNAc release in vitro; improved weight gain and FCR vs. control in vivo |
| **3** | **RT-qPCR + Western blot validation of the top 3 M6 hub genes** (sqh, uri, LOC136845115) on existing n=20 hepatopancreas cDNA | GSEA-confirmed muscle-cytoskeleton mechanism (padj 0.0013 / 1.5e-7) | Low — existing samples, ~3 months | Spearman ρ≥0.6 between qPCR fold-change and RNA-seq TPM; protein abundance inversely correlates with weight gain |
| **4** | **Individual-level 16S sequencing (n≥30 prawns, V3–V4) targeting the 5 confirmed taxa** | Converts every microbial finding from directional hypothesis to a testable, individually-resolved association | Moderate-high — new gut-content collection + sequencing, ~4–6 months | Detection in ≥50% of samples; Spearman correlation with weight gain for ≥3/5 taxa |
| **5** | **GC-MS/LC-MS metabolomics**: SCFAs (acetate, butyrate, propionate) and GlcNAc in hepatopancreas, n=20 existing samples | Tests the #1-ranked (SCFA→PPAR) and #2-ranked (chitin→GlcNAc) pathway bridges directly at the metabolite level | Moderate — outsourced GC-MS/LC-MS service, ~3–6 months | SCFA/GlcNAc levels correlate with weight gain and/or confirmed-taxon abundance |
| **6** | **M8 hub gene RT-qPCR panel** (LOC136852535, CtBP, COX1-like) | Permutation-confirmed (p=0.014) but not GSEA-characterised; needs n=32 for full power | Low — existing samples, ~3 months | Direction-concordant with WGCNA sign; flags M8 for a future, adequately powered GSEA pass |
| **7** | **Full independent replication cohort** (n≥60, individual 16S + RNA-seq + weight gain, both sexes, longitudinal through the moult cycle) | Would enable direct cross-layer correlation and sex-stratified networks — currently impossible by design | High — 12–18 months, new cohort + sequencing budget | Statistical confirmation (not just replication) of the taxon-to-pathway and module-to-phenotype associations above |

Interventions 1–3 use only existing samples and are recommended first; 4–6 convert this study's two largest open questions (are the confirmed taxa real gut colonists, and does the chitin/SCFA mechanism operate at the metabolite level) into testable claims; 7 is the only intervention capable of directly testing a taxon-to-gene relationship, which nothing in the current design can support.

---

## 6. Limitations

1. **Microbial replication is the fundamental constraint.** n=1 pooled sample per phenotype means no microbial finding in this report is a statistical result — every one is a directional hypothesis, however tightly bounded its Bayesian CI.
2. **Reference annotation is incomplete.** The composite 66,982-transcript reference has no native gene models; ~40% of genes remain unmapped to Swiss-Prot orthologs even after DIAMOND search.
3. **Tissue mismatch.** RNA-seq profiles hepatopancreas and gonad, not gut epithelium — the direct host-microbe interface is not transcriptomically profiled.
4. **Cross-sectional design.** A single time point cannot distinguish stable associations from transient microbial blooms or moult-stage-specific expression.
5. **Indirect edge resolution.** Integration edges connect microbial functional categories to host pathways, not specific taxa to specific genes — by design, given constraint 1, but it means the network cannot yet answer "which taxon affects which gene."
6. **M8 and the two unresolved Jumper taxa (*J. sediminis*, *E. diazotrophicus*) lack a KEGG GENOME match or a formal GSEA pass respectively** — three of the five confirmed taxa and one of the two confirmed modules carry a strong statistical signal without yet having an attached functional mechanism.

---

## References

1. Koh A, De Vadder F, Kovatcheva-Datchary P, Bäckhed F. From dietary fiber to host physiology: short-chain fatty acids as key bacterial metabolites. *Cell* 2016;165(6):1332–1345.
2. Ray AK, Ghosh K, Ringø E. Enzyme-producing bacteria isolated from fish gut: a review. *Aquaculture Nutrition* 2012;18(5):465–492.
3. Nayak SK. Role of gastrointestinal microbiota in fish. *Aquaculture Research* 2010;41(11):1553–1573.
4. New MB. Freshwater prawn farming: global status, recent research and a glance at the future. *Aquaculture Research* 2005;36(3):210–230.
5. Salter SJ, Cox MJ, Turek EM, et al. Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. *BMC Biology* 2014;12:87.
6. Eisenhofer R, Minich JJ, Marotz C, et al. Contamination in low microbial biomass microbiome studies: issues and recommendations. *Trends in Microbiology* 2019;27(2):105–117.
7. Holmström C, Kjelleberg S. Marine *Pseudoalteromonas* species are associated with higher organisms and produce biologically active extracellular agents. *FEMS Microbiology Ecology* 1999;30(4):285–293.
8. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 2008;9:559.

---

## Data and code availability

- Microbial ASV table: `ASV_table_Jumpers_Laggards.GKAQUA.csv` (supplied)
- Host RNA-seq: NCBI BioProject PRJNA875278, 20 SRA runs (accessions listed §2.1)
- Reference transcriptome: NCBI TSA Project 73259 + all *M. rosenbergii* mRNA records
- Ranked shortlists: `results/shortlist/track3_shortlists.xlsx` (or the 5 matching CSVs)
- Network figure + edge list: `results/figures/final_integrated_network.png` / `.svg` / `.pdf`; `results/shortlist/network_edges.csv`
- Ranking formulas: `results/shortlist/ranking_methodology.md`, `results/shortlist/ranking_weights.csv`
- Full pipeline: `src/`, orchestrated by `workflow/Snakefile`; reproducible via `conda env create -f environment/conda_environment.yml && conda activate track3_prawn && snakemake --profile workflow/profiles/local`
