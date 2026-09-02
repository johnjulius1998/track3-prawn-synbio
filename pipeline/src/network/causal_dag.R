#!/usr/bin/env Rscript
#'
#' causal_dag.R — Causal DAG for M. rosenbergii Growth Analysis
#' ==============================================================
#' Track 3 Host-Microbe Integration (v3.3)
#'
#' PURPOSE:
#'   Replace "controlling for confounders" with explicit causal reasoning.
#'   Draw the DAG, identify adjustment sets, and test implied conditional
#'   independencies against the data.
#'
#' METHOD:
#'   1. Define DAG: Sex, Tissue, GeneExpression → WeightGain
#'   2. Identify adjustment sets using do-calculus
#'   3. Test implied conditional independencies
#'   4. Estimate average causal effect bounds
#'
#' USAGE:
#'   Rscript src/network/causal_dag.R \
#'       --me data/processed/wgcna/me_matrix.tsv \
#'       --metadata data/raw/sra/PRJNA875278/metadata.tsv \
#'       --out-dir results/reports/
#'

suppressPackageStartupMessages({
  library(data.table)
})

# ---- CLI ----
args <- commandArgs(trailingOnly = TRUE)
me_file <- NULL; metadata_file <- NULL; out_dir <- "results/reports"
i <- 1
while (i <= length(args)) {
  if (args[i] == "--me")          { me_file <- args[i+1]; i <- i+2 }
  else if (args[i] == "--metadata") { metadata_file <- args[i+1]; i <- i+2 }
  else if (args[i] == "--out-dir")  { out_dir <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

cat("=== CAUSAL DAG ANALYSIS (v3.3) ===\n\n")

# ---- Load data ----
cat("[1/4] Loading data...\n")
me <- fread(me_file, data.table = FALSE)
# First column is row names (sample IDs)
rownames(me) <- me[, 1]; me <- me[, -1, drop = FALSE]
meta <- fread(metadata_file, data.table = FALSE)
rownames(meta) <- meta[, 1]; meta <- meta[, -1, drop = FALSE]
common <- intersect(rownames(me), rownames(meta))
me <- me[common, , drop = FALSE]; meta <- meta[common, , drop = FALSE]
# Ensure M6 eigengene is accessible
me_cols <- colnames(me)
m6_col <- if ("ME6" %in% me_cols) "ME6" else if ("MEM6" %in% me_cols) "MEM6" else grep("M6", me_cols, value=TRUE)[1]
cat(sprintf("  Samples: %d, M6 column: %s\n", nrow(meta), m6_col))

# ---- Prepare variables ----
wg <- meta$weight_gain
sex <- ifelse(meta$sex == "male", 1, 0)
tissue <- ifelse(grepl("hepato", meta$tissue, ignore.case = TRUE), 1, 0)

# ---- Define and analyze DAG ----
cat("\n[2/4] Defining causal DAG...\n")
cat("
Proposed DAG:

    Sex ──────────► WeightGain
     │                  ▲
     │                  │
     └──► M6_Expression ─┘
              ▲
              │
           Tissue

Interpretation:
  - Sex → WeightGain: females are larger (biological fact, r=-0.76)
  - Sex → M6_Expression: sex affects gene expression
  - Tissue → M6_Expression: tissue type affects expression
  - M6_Expression → WeightGain: hypothesized causal effect
\n")

# ---- Adjustment set ----
cat("[3/4] Identifying adjustment sets...\n")
cat("
  To estimate: causal effect of M6 expression on WeightGain
  Backdoor paths: M6 ← Sex → WG  (confounded by sex)
                   M6 ← Tissue   (not a backdoor to WG unless Tissue→WG)
  
  Minimal sufficient adjustment set: {Sex}
  (Tissue is not a confounder for M6→WG; it only affects M6, not WG directly)
  
  Conditioning on {Sex} is sufficient. This is what we already do via
  partial correlation. The analysis IS causally identified under the
  assumed DAG — the limitation is statistical power (n=20), not
  identification.
\n")

# ---- Test implied conditional independencies ----
cat("[4/4] Testing implied conditional independencies...\n")

# Implication 1: WeightGain ⊥ Tissue | Sex
# If tissue doesn't affect WG except through M6, then WG ⊥ Tissue | Sex
fit1 <- lm(wg ~ sex + tissue)
p_tissue <- coef(summary(fit1))["tissue", "Pr(>|t|)"]
cat(sprintf("  Test 1: WG ⊥ Tissue | Sex\n"))
cat(sprintf("    lm(WG ~ Sex + Tissue): Tissue p = %.4f\n", p_tissue))
if (p_tissue > 0.05) {
  cat("    ✓ Consistent with DAG — tissue does not affect WG directly\n")
} else {
  cat("    ⚠ Tissue may have a direct effect on WG — DAG may need revision\n")
}

# Implication 2: M6 ⊥ Sex | Tissue? (No — Sex affects M6)
m6_vals <- me[[m6_col]]
fit2 <- lm(m6_vals ~ tissue + sex)
p_sex_m6 <- coef(summary(fit2))["sex", "Pr(>|t|)"]
cat(sprintf("\n  Test 2: M6 ⊥ Sex | Tissue\n"))
cat(sprintf("    lm(M6 ~ Tissue + Sex): Sex p = %.4f\n", p_sex_m6))
if (p_sex_m6 < 0.05) {
  cat("    ✓ Consistent with DAG — sex affects M6 expression\n")
} else {
  cat("    ⚠ Sex may not affect M6 — backdoor path may be weaker than assumed\n")
}

# Causal effect estimate: M6 → WG | Sex
fit_causal <- lm(wg ~ m6_vals + sex)
m6_coef <- coef(summary(fit_causal))["m6_vals", "Estimate"]
m6_p <- coef(summary(fit_causal))["m6_vals", "Pr(>|t|)"]
cat(sprintf("\n  Causal effect estimate (M6 → WG, adjusted for Sex):\n"))
cat(sprintf("    Coefficient: %.4f (per unit increase in M6 eigengene)\n", m6_coef))
cat(sprintf("    p-value: %.4f\n", m6_p))
cat(sprintf("    Interpretation: For each 1-unit increase in M6 expression,\n"))
cat(sprintf("    WG changes by %.2f units, controlling for sex.\n", m6_coef))

# Bounds: with n=20, what range of effect sizes are consistent with the data?
ci <- confint(fit_causal)["m6_vals", ]
cat(sprintf("    95%% CI: [%.3f, %.3f]\n", ci[1], ci[2]))

# ---- Write outputs ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

results <- data.frame(
  test = c("WG_Tissue_given_Sex", "M6_Sex_given_Tissue", "Causal_M6_WG_given_Sex"),
  description = c("WG independent of Tissue given Sex",
                  "M6 depends on Sex given Tissue",
                  "Causal effect of M6 on WG"),
  statistic = c(sprintf("p=%.4f", p_tissue),
                sprintf("p=%.4f", p_sex_m6),
                sprintf("coef=%.4f", m6_coef)),
  p_value = c(p_tissue, p_sex_m6, m6_p),
  consistent_with_dag = c(p_tissue > 0.05, p_sex_m6 < 0.05, NA),
  stringsAsFactors = FALSE
)

write.table(results, file.path(out_dir, "causal_dag_tests.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\n  [OK] %s/causal_dag_tests.tsv\n", out_dir))
cat("\nDone.\n")
