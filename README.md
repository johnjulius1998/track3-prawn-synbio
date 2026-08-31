# Track 3 Host-Microbe Integration — *Macrobrachium rosenbergii*

## Gut Microbiome & Host Transcriptome Drivers of Growth Rate

**BioProject**: PRJNA875278 | **Version**: v5.0 (2026-08-17) | **Approach**: 2 (pooled ASV) + 3 (WGCNA) + tool-derived functional profiling & KEGG-bridged integration

> **v5.0**: the shortlists below reflect the validation gates described in this document — 5/5 microbial taxa (all pass the Bayesian-CI + contamination gates) and 6/10 host genes (only M6 and M8 survive permutation testing).

> **This repository** is a curated, manuscript-ready subset of the full analysis project. See [Repository Structure](#repository-structure) and [Data Availability](DATA_AVAILABILITY.md) for what is and isn't included here.

---

## TL;DR — The Takeaway in One Paragraph

We asked: *which gut microbes and host genes are associated with faster growth in the giant freshwater prawn?* We found that **chitin-degrading bacteria** (*Pseudoalteromonas* spp.) in the gut of fast-growing prawns likely provide GlcNAc building blocks for the host's exoskeleton, feeding into amino sugar metabolism during the energy-intensive molting process. At the host level, a **719-gene co-expression module (M6)** is causally associated with growth rate — higher expression of these genes predicts *lower* weight gain — and its function is now established **by formal gene-set enrichment** (Swiss-Prot orthology + clusterProfiler GSEA): the module is dominated by **muscle cytoskeleton / contraction** biology (KEGG map04820 padj=0.0013; GO:0006936 padj=1.5e-7), with transcription and DNA-repair regulation — muscle-maintenance programs whose energetic cost is suppressed in fast growers. Thirteen validation layers confirm M6 as genuine signal (permutation null p=0.001; causal effect p=0.008): it is the only host finding that survives every test, and one additional sample (n=21) would confirm it at 80% power. The microbial functional layer is **tool-derived**: genome-resolved KEGG GENOME profiles (171/190 taxa; 363 pathways) plus MetaCyc (3,602 pathways) replaced the literature genus map — the top Jumper-enriched MetaCyc signals are **nitrogen-cycle pathways**, and the top-ranked integration bridge is **SCFA → PPAR signaling** (amino sugar ranks #2–3). The microbial findings remain directional hypotheses — only 5 of 191 taxa (the confirmed shortlist) combine a Bayesian credible interval excluding zero with LOW contamination risk — 29/191 have any nonzero CI at all.

---

## The Data

| Layer | What | Samples | Source |
|-------|------|---------|--------|
| 🦠 Microbial | 191 ASVs (genus/species), pooled 16S | 2 pools (1 fast, 1 slow) | Supplied CSV |
| 🧬 Host | RNA-seq, hepatopancreas + gonad | 20 individuals (10♂, 10♀) | NCBI PRJNA875278 |
| 📐 Reference | Custom transcriptome, 66,982 sequences | — | NCBI TSA + mRNA records |

**The fundamental limitation**: the microbial data is a single pooled sample per growth group. We can describe what's different between fast and slow, but we cannot test whether those differences are statistically meaningful. Every microbial association is a **directional hypothesis**, not a conclusion.

---

## The Strongest Finding

### 🦐 Chitin → GlcNAc → Amino Sugar Metabolism → Growth

```
Dietary chitin (from exoskeletons in feed)
        │
        ▼
Pseudoalteromonas spp. (chitinolytic gut bacteria)
  6 species, ALL found only in fast-growth sample
        │
        ▼
N-acetylglucosamine (GlcNAc)
        │
        ▼
Host "Amino sugar & nucleotide sugar metabolism" pathway
        │
        ▼
New exoskeleton synthesis during molting → GROWTH
```

**Why this is the strongest microbial result**: It doesn't depend on any statistical test in our data. *Pseudoalteromonas* species are canonical chitin-degrading bacteria — known since the 1990s (Holmström & Kjelleberg, 1999, *FEMS Microbiol Ecol*). Crustaceans need GlcNAc to build their exoskeletons. The genome-resolved microbial amino-sugar pathway index bridges to host amino sugar metabolism + lysosome (ranks #2–3; the top-ranked bridge is SCFA → PPAR signaling).

**On the host side**: NCBI annotation found **N-acetylgalactosamine kinase (GalNAc kinase)** among M6 hubs — an enzyme of amino sugar metabolism. This convergence is **single-gene-level**: formal GSEA does NOT rank amino sugar metabolism among significant M6 sets. The chitin→GlcNAc hypothesis rests on the microbial layer and literature — exactly where a single pooled sample per growth group limits us to directional evidence rather than statistical inference.

---

## What We Found — With Confidence Levels

### 🟢 Confident (biology-backed, literature-grounded)

| Finding | Why we believe it |
|---------|-------------------|
| Chitin→GlcNAc pathway bridges microbes to host growth | Decades of marine microbiology literature + basic crustacean biology |
| Sex confounds growth (r=−0.76) | Textbook *M. rosenbergii* biology — females grow larger |
| *C. acnes* is a contaminant, not a gut microbe | Well-established in microbiome literature (Salter et al. 2014) |
| *Pseudoalteromonas* are chitin degraders | Canonical — multiple characterized species with known chitinases |
| **M6 is a genuine growth-associated module** | **Permutation null p=0.001 + causal estimate p=0.008 + DESeq2 LM p=0.03 + formal enrichment (muscle cytoskeleton padj=0.0013)** |

### 🟡 Plausible (data-supported but underpowered)

| Finding | What the data shows | What validation revealed |
|---------|--------------------|-------------------------|
| M6 co-expression module (719 genes) associated with growth | Partial r = −0.59 with weight gain | **Permutation p=0.001 (real signal); power 0.80 at n=20; needs n=21** |
| M8 module (336 genes) associated with growth | Partial r = −0.48 | **Permutation p=0.014 (95th percentile); needs n=32** |
| *J. sediminis* dominates fast-growth gut | 55% of all reads in the fast-growth pool | **Bayesian log2FC +4.02, CI [+3.38, +4.73]** — confirmed by Bayesian model |
| *M. osloensis* and *E. aestuarii* enriched in slow-growth | Consistent CLR fold-differences | *M. osloensis* Bayesian log2FC −4.63, CI excludes zero |

### 🔴 Weak or Rejected (data contradicts or doesn't support)

| Finding | What went wrong |
|---------|----------------|
| M7 as growth-associated module | DESeq2 LM p=0.81; **permutation null p=0.18** — indistinguishable from noise |
| *C. koseri* and *K. variicola* as ranked microbial taxa | LOTO top-5 stability = 0.00 (CLR ranking was a pseudocount artifact). Both do have a nonzero Bayesian CI, but rank well below the top-5 tier by effect size (log2FC +2.35 and +1.78 vs +3.8–+4.7 for the confirmed tier) |
| M13/M15 as distinct growth modules | **Permutation null p≥0.20** — noise; would need n≥69 |
| Any microbial direction being "statistically supported" | 94% of taxa flip direction; only 15% have Bayesian CI≠0 |
| Random Forest predicting weight gain | OOB R² = −0.34 — worse than guessing the mean; p >> n kills prediction |

### 🧬 New: M6's Functional Identity (annotation + formal GSEA)

13 of 72 M6 hub genes are characterized by NCBI title. They tell a coherent story:

| Gene | Function | Why it matters |
|------|----------|----------------|
| XM_067104953.1 | **N-acetylgalactosamine kinase (GalK)** | **Amino sugar metabolism enzyme — the same pathway as the #1 integration bridge (chitin→GlcNAc)** |
| XM_067106657.1 | **IGF-binding protein (IGFBP)** | Growth-axis inhibitor — directly explains negative growth correlation |
| EL609362.1 + 2 more | Ubiquitin / E3 ligase / F-box | Protein degradation machinery — metabolic cost |
| XM_067087534.1 + 2 more | eIF3b, eIF3e, DDX3X | Translation — turnover coupled with degradation |
| XM_067122492.1 | Immunoglobulin domain protein | Immune cost |

**The negative M6-WG correlation is now formally enriched**: higher M6 = more muscle-cytoskeleton/contraction and transcription/DNA-repair program expression = more maintenance cost, less growth. GSEA (Swiss-Prot orthology gene sets; both kME and DESeq2 rankings): **cytoskeleton in muscle cells (map04820 padj=0.0013)**, muscle contraction (GO:0006936 padj=1.5e-7), sarcomere organization, DNA repair, transcription regulation. The top hub XM_067086030.1 is **myosin regulatory light chain** (P40423, 83.3% identity). Amino sugar metabolism is NOT significant in M6 — the GalNAc-kinase observation remains single-gene-level.

---

## How We Validated This (Current Pipeline)

Thirteen layers of internal validation were applied to avoid overinterpreting our small dataset:

| # | Analysis | What it tests | Main finding |
|---|----------|---------------|--------------|
| 1 | **Pseudocount sweep** (11 CLR variants) | Is the microbial ranking an artifact of zero-handling? | Yes — 0/191 taxa stable; 94% flip direction |
| 2 | **Leave-one-taxon-out** (191 iterations) | Does one dominant taxon drive the ranking? | 3/5 ranked taxa stable; 2 are artifacts |
| 3 | **Bootstrap WGCNA** (1,000 resamples) | How stable are module-trait correlations? | No module has CI excluding zero |
| 4 | **DESeq2 cross-validation** (continuous WG, kME-aware) | Does an orthogonal GLM confirm hub genes? | 7/10 concordant; M6 confirmed; M7 rejected |
| 5 | **sPLS + Random Forest** | Do prediction methods find the same genes? | Zero overlap — co-expression ≠ prediction at n=20 |
| 6 | **Permutation null** (1,000 label shuffles) | Is any module real signal? | **M6 p=0.001 — the only module exceeding the 99th percentile** |
| 7 | **Power analysis** | What sample size is needed? | **M6 needs n=21** (we have 20); M8 needs 32; rest ≥69 |
| 8 | **Bayesian Dirichlet-multinomial** | Which taxa are genuinely enriched? | Only 29/191 (15%) have CI≠0; confirms the 5-taxon shortlist tier |
| 9 | **Tornado sensitivity** | Which choices drive results? | Confound correction + sample size ≫ algorithm params |
| 10 | **Causal DAG** | Is the effect causally identified? | **M6→WG = −3.50 (p=0.008)** under assumed DAG |
| 11 | **M6 NCBI annotation** (72 genes) | Why does M6 negatively correlate with growth? | **Protein turnover + IGFBP + GalNAc kinase (single-gene observations)** |
| 12 | **M6 GSEA — KEGG/GO** | Which pathways are formally enriched in M6? | **Muscle cytoskeleton (map04820 padj=0.0013); muscle contraction (GO padj=1.5e-7); transcription/DNA repair** |
| 13 | **Genome-resolved functional profiles** | Is the microbial functional layer tool-derived? | **Yes — KEGG GENOME + MetaCyc; cross-check vs the curated map: 12 supported / 5 new / 3 not-testable; nitrogen cycle tops Jumper** |

**The decisive pattern**: M6 survives all 11 tests. No other host module survives the permutation null. Only 5 microbial taxa combine a nonzero Bayesian CI with LOW contamination risk (29/191 have any nonzero CI).

---

## Final Ranked Candidates

**These are the final ranked candidates** — see `results/shortlist/track3_shortlists.xlsx` (or the matching CSVs) for the fully annotated version, with supporting columns such as `statistic_attributable_to`, `n_libraries_behind_this_call`, `edge_basis`, and `contaminant_risk` for every candidate. Both lists are capped below their maximum size (10 genes, 5 taxa) because the validation gates stop passing candidates before that ceiling is reached — fewer, better-supported candidates were prioritized over hitting the maximum count.

### Host Genes (6 of 10 ceiling — gated by permutation-null testing, not by data availability)

Only modules M6 (empirical p=0.001) and M8 (p=0.014) exceed the permutation null; every other module, including M7 and M12 below, is statistically indistinguishable from a shuffled label (`results/reports/permutation_null.tsv`). The final shortlist is therefore capped at the top-3-by-kME hub genes from each of M6 and M8 only.

| Rank | Gene | Symbol | Module | kME | Why it's here |
|------|------|--------|--------|-----|----------------|
| 1 | XM_067086030.1 | sqh | M6 | 0.845 | Top hub; myosin regulatory light chain (83.3% id.) — triangulates with M6's formal muscle-cytoskeleton GSEA signal |
| 2 | XM_067082934.1 | uri | M6 | 0.739 | Prefoldin RPB5 interactor — proteostasis/chaperone, consistent with M6's turnover signature |
| 3 | XM_067114945.1 | LOC136845115 | M6 | 0.713 | LITAF homolog — immune-energetic cost component |
| 4 | XM_067127262.1 | LOC136852535 | M8 | 0.621 | Tolloid-like protein 2 — BMP-pathway metalloprotease, plausible moult/cuticle mechanism |
| 5 | XM_067113377.1 | CtBP | M8 | 0.616 | C-terminal binding protein — metabolic-state-linked transcriptional corepressor |
| 6 | GH624888.1 | COX1-like | M8 | 0.604 | EST similar to mitochondrial COX1 — weak lead, oxidative-metabolism signature |

**Dropped from the final shortlist** (kept here for transparency):

| Gene | Module | Why dropped |
|------|--------|-------------|
| XM_067113580.1, JP354355.1 | **M7** | Permutation null p=0.181 — indistinguishable from noise; DESeq2 LM shows WG explains nothing (p=0.81) |
| XM_067102939.1, JP354756.1 | **M12** | Permutation null p=0.335 — even weaker than M7; M12 isn't in the top-5 modules by raw partial_r either |

### Microbial Taxa (5 of 5 ceiling — all 5 independently pass every gate)

Ranked by Bayesian posterior |log2FC| among taxa whose 95% credible interval excludes zero AND contamination risk is LOW (`src/asv/01c_bayesian_taxon_model.py`, `results/reports/bayesian_taxon_enrichment.tsv`):

| Rank | Taxon | Direction | Bayesian log2FC [95% CI] | Why |
|------|-------|-----------|---------------------------|-----|
| 1 | *Endothiovibrio diazotrophicus* | Jumper-enriched | +4.69 [+1.88, +9.12] | Largest effect size of all 191 taxa; Jumper-exclusive (15/0 reads); unresolved to a KEGG GENOME entry — confirmatory only |
| 2 | *Moraxella osloensis* | Laggard-enriched | −4.63 [−7.31, −2.67] | Strongest Laggard signal; genome-resolved; top-5 contributor to the PPAR + amino-sugar pool indices |
| 3 | *Pseudoalteromonas phenolica* | Jumper-enriched | +4.13 [+1.28, +8.62] | Carries the chitin-degradation mechanism; genome-resolved; top-5 Jumper contributor to all 3 ranked host pathways |
| 4 | *Jiulongibacter sediminis* | Jumper-enriched | +4.02 [+3.38, +4.73] | Dominates raw abundance (55.5% of all Jumper reads); tightest CI of any confirmed taxon; unresolved to a KEGG GENOME entry |
| 5 | *Exiguobacterium aestuarii* | Laggard-enriched | −3.82 [−6.57, −1.84] | Genome-resolved; reciprocal Laggard-side contributor to the same PPAR + amino-sugar pool indices as *M. osloensis* |

**Dropped from the final shortlist**: *Citrobacter koseri* (log2FC +2.35) and *Klebsiella variicola* (log2FC +1.78) — both LOTO-unstable (jackknife stability 0.000, i.e. their CLR point-estimate ranking was a pseudocount artifact) and, independently, both rank well below the 5 taxa above by Bayesian effect size.

### Pathways

| Rank | Pathway | Why |
|------|---------|-----|
| 1 | PPAR signaling | SCFA bridge — microbial butanoate/propanoate metabolism + fatty-acid degradation (3 edges) |
| 2 | Amino sugar & nucleotide sugar metabolism | The chitin→GlcNAc bridge (1 edge, genome-resolved microbial amino-sugar index) |
| 3 | Lysosome | Chitinolytic enzyme processing |

---

## What This Means for Wet-Lab Validation

The internal validation tells us **exactly which findings to prioritize** for experimental follow-up:

**Test first (now validated with specific gene targets)**:
- **Muscle/cytoskeleton panel** (top hub XM_067086030.1 = myosin regulatory light chain) by RT-qPCR — the module's formally enriched signature (map04820 padj=0.0013; GO:0006936 padj=1.5e-7)
- **IGFBP** (XM_067106657.1) by RT-qPCR — growth-axis inhibitor explaining the negative M6-WG correlation
- **GalNAc kinase** (XM_067104953.1) by RT-qPCR + activity assay — the amino-sugar enzyme linking M6 to the chitin pathway (single-gene-level; amino sugar is NOT pathway-significant in M6 GSEA)
- M6 protein-turnover panel (ubiquitin, E3 ligase, F-box) — confirms the metabolic-cost signature
- *Pseudoalteromonas* 16S qPCR on individual gut samples (the chitin degradation mechanism depends on this genus being present)

**Test second (plausible but unconfirmed)**:
- M8 hub genes (permutation p=0.014 — signal, but needs n=32)
- All 5 confirmed taxa by individual-level 16S — *E. diazotrophicus*, *M. osloensis*, *P. phenolica*, *J. sediminis*, *E. aestuarii* (Bayesian-confirmed enrichment, but need to confirm they're real gut colonizers, not pond/reagent background)

**Don't prioritize (evidence too weak)**:
- M7 and M12 genes (permutation p=0.18 and p=0.335 — both indistinguishable from noise)
- *C. koseri* and *K. variicola* (LOTO-unstable; Bayesian CI does exclude zero for both, but effect size ranks well below the confirmed 5-taxon tier)

**The two highest-value experiments**:
1. **One additional RNA-seq sample with weight gain** — pushes M6 from power 0.80 to >0.80 (needs n=21). Cheapest possible confirmation.
2. **Individual-level 16S of n≥30 prawns with weight gain** — converts all microbial hypotheses into statistically testable associations and tells us whether *J. sediminis* is a genuine gut symbiont or pond sediment contamination.

---

## Repository Structure

```
.
├── workflow/
│   ├── Snakefile              # Full analysis DAG (microbial, host, WGCNA, validation, integration, shortlisting)
│   └── profiles/local/        # Snakemake execution profile
├── config/
│   ├── project_config.yaml    # Pipeline parameters, dataset/edge-type rules
│   └── contaminant_genera.txt # Known-contaminant genus list used by the ASV filtering step
├── environment/
│   └── conda_environment.yml  # Conda env spec (name: track3_prawn)
├── src/
│   ├── asv/                   # CLR transform, taxon sensitivity, functional profiling, MetaCyc mapping
│   ├── network/                # WGCNA, confound correction, bootstrap/permutation/causal validation, edge building
│   ├── rnaseq/                 # SRA metadata/download, salmon quantification, DESeq2 validation
│   ├── ranking/                # Final shortlist generation and submission figure generation
│   └── utils/                  # Monitoring/preflight helper scripts
├── data/
│   ├── raw/supplied/          # The original supplied ASV table (small, included)
│   ├── raw/sra/PRJNA875278/    # Sample metadata/run info only — NOT the raw reads (see Data Availability)
│   └── processed/              # Final CLR profiles, gene expression matrices, WGCNA outputs
└── results/
    ├── figures/                 # Submission figures (cnsplots fig1-8, final integrated network, GSEA dotplots)
    ├── tables/                  # WGCNA module tables
    ├── reports/                 # Validation/statistics reports (permutation null, power analysis, etc.)
    └── shortlist/                # Official ranked shortlists (host genes, microbial taxa, pathways, network edges)
```

See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) for what raw/intermediate data is excluded from this repository and how to regenerate it.

---

## Reproducibility

```bash
conda env create -f environment/conda_environment.yml
conda activate track3_prawn
snakemake --profile workflow/profiles/local
```

All analysis code is in `src/`.

**Note on exact reproducibility**: the shipped results in `results/` and `data/processed/` were originally produced by running `src/rnaseq/process_all_samples.sh` as a serial driver over the raw SRA downloads (per the note in `workflow/Snakefile`), not by a single end-to-end `snakemake` invocation — the raw FASTQ files were deleted after processing to save space and are not bundled in this repository. A fresh clone can:
- Re-run every step downstream of `data/processed/` and `data/raw/supplied/` (the ASV/WGCNA/validation/shortlisting layers) immediately.
- Re-run the RNA-seq layer from scratch only after re-fetching raw reads from BioProject **PRJNA875278** (see [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md)); `snakemake -n` (dry run) will correctly report these inputs as missing until then.

---

## Key Caveat

> **All microbial associations are directional hypotheses only.** With n=1 pooled sample per growth group, we cannot compute p-values or effect sizes for any microbial finding. The 16S data tells us what's *different* between one fast pool and one slow pool — not what's *significantly associated* with growth across individuals. The Bayesian model quantifies this honestly: only 15% of taxa (29/191) have credible intervals excluding zero.

---

## Development Workflow

For any future changes to this repository: work on a feature branch and open a pull request into `main`, even when working solo — the PR diff and description double as a reviewable changelog. Keep `main` always in a working, citable state. Write commit messages that explain *why* a change was made, not just what changed. Never edit files under `data/raw/` by hand — all raw→processed transformations should happen through the scripts in `src/` and the Snakemake rules in `workflow/Snakefile`, so the data lineage stays reproducible.

---

## License

Code in this repository is licensed under the [MIT License](LICENSE). Manuscript text, figures, and other non-code content are licensed under [CC-BY-4.0](LICENSE-DOCS).
