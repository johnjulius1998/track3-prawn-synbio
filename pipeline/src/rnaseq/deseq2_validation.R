#!/usr/bin/env Rscript
#'
#' deseq2_cross_validate.R — DESeq2 Orthogonal Validation of WGCNA Hub Genes
#' ==========================================================================
#' Track 3 Host-Microbe Integration (v3.1 — Fixed kME sign convention +
#' continuous WG model + module-level eigengene test + tissue covariate)
#'
#' PURPOSE:
#'   Cross-validate WGCNA hub genes using DESeq2 with a model that directly
#'   matches the WGCNA hypothesis: continuous weight_gain as predictor,
#'   controlling for sex and tissue.  Also tests module eigengenes via LM.
#'
#' v3.1 FIXES:
#'   1. Concordance now accounts for kME sign: expected_direction = sign(kME) × sign(partial_r)
#'      (the v3.0 version ignored kME sign, misclassifying M7 genes as discordant)
#'   2. Continuous WG model (n=20) replaces dichotomized high-vs-low (n=10)
#'   3. Module-level eigengene LM test added as orthogonal validation
#'   4. Tissue added as covariate in DESeq2 design
#'
#' METHOD:
#'   1. Load merged counts + metadata + hub genes + module eigengenes
#'   2. DESeq2 with design = ~ sex + tissue + weight_gain (continuous, n=20)
#'   3. Extract WG coefficient (log2FC per unit WG) for top 10 hub genes
#'   4. Check direction concordance: sign(DESeq2_WG_coef) == sign(kME) × sign(partial_r)
#'   5. Module-level: lm(ME ~ sex + tissue + WG) for all 20 modules
#'   6. Report per-gene and per-module concordance
#'
#' OUTPUTS:
#'   results/reports/
#'     deseq2_cross_validation.tsv       — Per-gene results for top 10 hub genes
#'     deseq2_module_eigengene_test.tsv  — Per-module LM results
#'     deseq2_all_results.tsv            — Genome-wide WG coefficient results
#'
#' USAGE:
#'   Rscript src/rnaseq/deseq2_validation.R \
#'       --counts data/processed/gene_expression/merged_counts.tsv \
#'       --metadata data/raw/sra/PRJNA875278/metadata.tsv \
#'       --hub-genes results/shortlist/host_genes.csv \
#'       --me data/processed/wgcna/me_matrix.tsv \
#'       --out-dir results/reports/
#'

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

# ---- CLI ----
args <- commandArgs(trailingOnly = TRUE)
counts_file <- NULL; metadata_file <- NULL; hub_file <- NULL
me_file <- NULL; out_dir <- "results/reports"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--counts")           { counts_file <- args[i+1]; i <- i+2 }
  else if (args[i] == "--metadata")    { metadata_file <- args[i+1]; i <- i+2 }
  else if (args[i] == "--hub-genes")   { hub_file <- args[i+1]; i <- i+2 }
  else if (args[i] == "--me")          { me_file <- args[i+1]; i <- i+2 }
  else if (args[i] == "--out-dir")     { out_dir <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

# ---- Load data ----
cat("=== DESeq2 Orthogonal Validation (v3.1 — Fixed + Continuous WG + Module LM) ===\n\n")

cat("[1/6] Loading data...\n")
counts <- fread(counts_file, data.table = FALSE)
rownames(counts) <- counts[, 1]; counts <- counts[, -1, drop = FALSE]
counts <- round(as.matrix(counts))

metadata <- fread(metadata_file, data.table = FALSE)
rownames(metadata) <- metadata[, 1]; metadata <- metadata[, -1, drop = FALSE]

hub_genes <- read.csv(hub_file, stringsAsFactors = FALSE)

cat(sprintf("  Counts: %d genes x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("  Metadata: %d samples\n", nrow(metadata)))
cat(sprintf("  Hub genes to validate: %d\n", nrow(hub_genes)))

# Align samples
common_samples <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]
cat(sprintf("  Common samples: %d\n", length(common_samples)))

# ---- Prepare continuous model covariates ----
cat("\n[2/6] Preparing continuous WG model (n=20, design = ~ sex + tissue + weight_gain)...\n")

# Sex as factor
metadata$sex <- factor(metadata$sex)

# Tissue as factor
metadata$tissue <- factor(metadata$tissue)

# WG as continuous (keep raw values)
wg_values <- metadata$weight_gain
cat(sprintf("  WG: mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
            mean(wg_values), sd(wg_values), min(wg_values), max(wg_values)))

# ---- Filter low-count genes ----
cat("\n[3/6] Filtering low-count genes...\n")
keep_genes <- rowSums(counts >= 10) >= 3
counts_filt <- counts[keep_genes, , drop = FALSE]
cat(sprintf("  Kept %d / %d genes\n", nrow(counts_filt), sum(keep_genes)))

# ---- DESeq2: continuous WG (n=20) + tissue covariate ----
cat("\n[4/6] Running DESeq2 with design = ~ sex + tissue + weight_gain (n=20)...\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData = metadata,
  design = ~ sex + tissue + weight_gain
)
dds <- DESeq(dds, quiet = TRUE)

# Extract WG coefficient (log2 fold change per unit increase in weight_gain)
res <- results(dds, name = "weight_gain", alpha = 0.1)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
rownames(res_df) <- NULL

# Also run the dichotomized model for comparison
cat("  Also running dichotomized model for comparison...\n")
wg <- metadata$weight_gain
wg_q <- quantile(wg, probs = c(0.25, 0.75), na.rm = TRUE)
metadata$growth_group <- ifelse(wg >= wg_q[2], "high",
                         ifelse(wg <= wg_q[1], "low", "mid"))
metadata$growth_group <- factor(metadata$growth_group, levels = c("low", "mid", "high"))
keep_dich <- metadata$growth_group %in% c("low", "high")
dds_dich <- DESeqDataSetFromMatrix(
  countData = counts_filt[, keep_dich, drop = FALSE],
  colData = metadata[keep_dich, , drop = FALSE],
  design = ~ sex + tissue + growth_group
)
dds_dich <- DESeq(dds_dich, quiet = TRUE)
res_dich <- results(dds_dich, contrast = c("growth_group", "high", "low"), alpha = 0.1)

# ---- Cross-validate hub genes (fixed sign convention) ----
cat("\n[5/6] Cross-validating top 10 hub genes (kME-aware sign convention)...\n")

hub_validated <- data.frame(stringsAsFactors = FALSE)

for (i in seq_len(nrow(hub_genes))) {
  g <- hub_genes$gene_id[i]
  mod <- hub_genes$associated_module[i]
  kME <- hub_genes$kME[i]
  mod_pr <- hub_genes$module_partial_r_wg_given_sex[i]
  
  # FIXED: expected direction = sign(kME) × sign(partial_r)
  # If kME > 0 and partial_r < 0: gene→ME positive, ME→WG negative → gene→WG negative
  # If kME < 0 and partial_r < 0: gene→ME negative, ME→WG negative → gene→WG positive
  expected_dir <- sign(kME) * sign(mod_pr)
  expected_dir_str <- ifelse(expected_dir > 0, "positive",
                      ifelse(expected_dir < 0, "negative", "zero"))
  
  # Also old WGCNA direction for comparison
  wgcna_dir_old <- ifelse(mod_pr < 0, "negative", "positive")
  
  # Find in continuous DESeq2 results
  match_row <- which(res_df$gene_id == g)
  
  if (length(match_row) == 1) {
    lfc_cont <- res_df$log2FoldChange[match_row]   # per-unit-WG coefficient
    pval_cont <- res_df$pvalue[match_row]
    padj_cont <- res_df$padj[match_row]
    basemean <- res_df$baseMean[match_row]
    
    # Find in dichotomized results for comparison
    match_dich <- which(rownames(res_dich) == g)
    lfc_dich <- if (length(match_dich) == 1) res_dich$log2FoldChange[match_dich] else NA_real_
    pval_dich <- if (length(match_dich) == 1) res_dich$pvalue[match_dich] else NA_real_
    
    # Concordance: sign(DESeq2 WG coefficient) == expected_dir
    if (!is.na(lfc_cont) && expected_dir != 0) {
      concordant <- (sign(lfc_cont) == expected_dir)
    } else {
      concordant <- NA
    }
    
    # Also check old logic for comparison
    concordant_old <- NA
    if (!is.na(lfc_dich)) {
      if (wgcna_dir_old == "negative") {
        concordant_old <- (lfc_dich < 0)
      } else {
        concordant_old <- (lfc_dich > 0)
      }
    }
    
    nom_sig <- !is.na(pval_cont) && pval_cont < 0.1
    
    hub_validated <- rbind(hub_validated, data.frame(
      rank = hub_genes$rank[i],
      gene_id = g,
      wgcna_module = mod,
      wgcna_kME = kME,
      wgcna_module_partial_r = mod_pr,
      expected_direction = expected_dir_str,
      wgcna_direction_old = wgcna_dir_old,
      deseq2_WG_coef = round(lfc_cont, 4),
      deseq2_pvalue_cont = round(pval_cont, 6),
      deseq2_padj_cont = round(padj_cont, 6),
      deseq2_lfc_dichotomized = round(lfc_dich, 4),
      deseq2_pvalue_dich = round(pval_dich, 6),
      deseq2_baseMean = round(basemean, 1),
      direction_concordant = concordant,
      concordant_old_logic = concordant_old,
      nominally_significant = nom_sig,
      stringsAsFactors = FALSE
    ))
  } else {
    hub_validated <- rbind(hub_validated, data.frame(
      rank = hub_genes$rank[i], gene_id = g,
      wgcna_module = mod, wgcna_kME = kME,
      wgcna_module_partial_r = mod_pr,
      expected_direction = expected_dir_str,
      wgcna_direction_old = wgcna_dir_old,
      deseq2_WG_coef = NA_real_, deseq2_pvalue_cont = NA_real_,
      deseq2_padj_cont = NA_real_, deseq2_lfc_dichotomized = NA_real_,
      deseq2_pvalue_dich = NA_real_, deseq2_baseMean = NA_real_,
      direction_concordant = NA, concordant_old_logic = NA,
      nominally_significant = NA, stringsAsFactors = FALSE
    ))
  }
}

# ---- Module-level eigengene LM test ----
cat("\n[6/6] Module-level eigengene linear model test...\n")

module_lm_results <- data.frame(stringsAsFactors = FALSE)

if (!is.null(me_file) && file.exists(me_file)) {
  me <- fread(me_file, data.table = FALSE)
  rownames(me) <- me[, 1]; me <- me[, -1, drop = FALSE]
  
  me_cols <- grep("^ME", colnames(me), value = TRUE)
  common_me <- intersect(rownames(me), rownames(metadata))
  me <- me[common_me, , drop = FALSE]
  meta_me <- metadata[common_me, , drop = FALSE]
  
  for (col in me_cols) {
    me_vals <- me[[col]]
    fit <- lm(me_vals ~ sex + tissue + weight_gain, data = meta_me)
    s <- summary(fit)$coefficients
    wg_row <- s["weight_gain", ]
    
    mod_name <- sub("^ME", "M", col)
    module_lm_results <- rbind(module_lm_results, data.frame(
      module = mod_name,
      wg_coef = round(wg_row["Estimate"], 6),
      wg_se = round(wg_row["Std. Error"], 6),
      wg_t = round(wg_row["t value"], 4),
      wg_pvalue = round(wg_row["Pr(>|t|)"], 6),
      r_squared = round(summary(fit)$r.squared, 4),
      stringsAsFactors = FALSE
    ))
  }
  
  # Add WGCNA partial_r for comparison
  # (computed from me and metadata, matching generate_final_shortlists logic)
  wg_vec <- meta_me$weight_gain
  sex_vec <- ifelse(meta_me$sex == "male", 1, 0)
  
  for (i in seq_len(nrow(module_lm_results))) {
    mn <- module_lm_results$module[i]
    me_col <- paste0("ME", substring(mn, 2))
    if (me_col %in% colnames(me)) {
      mv <- me[[me_col]]
      r_wg <- cor(mv, wg_vec)
      r_sex <- cor(mv, sex_vec)
      r_wg_sex <- cor(wg_vec, sex_vec)
      num <- r_wg - r_sex * r_wg_sex
      denom <- sqrt((1 - r_sex^2) * (1 - r_wg_sex^2))
      partial_r <- if (denom != 0) num / denom else 0
      module_lm_results$wgcna_partial_r[i] <- round(partial_r, 4)
      # Direction concordance: sign(lm coefficient) == sign(partial_r)?
      module_lm_results$direction_concordant[i] <- 
        sign(module_lm_results$wg_coef[i]) == sign(partial_r)
    }
  }
}

# ---- Summary ----
cat("\n=== PER-GENE RESULTS (fixed kME-aware sign convention) ===\n")
n_concordant <- sum(hub_validated$direction_concordant, na.rm = TRUE)
n_tested <- sum(!is.na(hub_validated$deseq2_WG_coef))
n_sig <- sum(hub_validated$nominally_significant, na.rm = TRUE)
n_old_conc <- sum(hub_validated$concordant_old_logic, na.rm = TRUE)

cat(sprintf("  Hub genes tested: %d / %d\n", n_tested, nrow(hub_validated)))
cat(sprintf("  Direction concordant (fixed):    %d / %d (%.1f%%)\n",
            n_concordant, n_tested, 100 * n_concordant / max(n_tested, 1)))
cat(sprintf("  Direction concordant (old logic): %d / %d (%.1f%%)\n",
            n_old_conc, n_tested, 100 * n_old_conc / max(n_tested, 1)))
cat(sprintf("  Nominally significant (p<0.1): %d / %d\n", n_sig, n_tested))

cat("\n  Per-gene validation:\n")
for (i in seq_len(nrow(hub_validated))) {
  r <- hub_validated[i, ]
  if (is.na(r$deseq2_WG_coef)) {
    cat(sprintf("    #%d %-25s NOT FOUND\n", r$rank, r$gene_id))
  } else {
    conc_str <- ifelse(isTRUE(r$direction_concordant), "✓ concordant", "✗ discordant")
    old_str  <- ifelse(isTRUE(r$concordant_old_logic), "(old:✓)", "(old:✗)")
    sig_str <- ifelse(r$nominally_significant, " (p<0.1)", "")
    cat(sprintf("    #%d %-25s M%-5s kME=%+.3f exp_dir=%s DESeq2_coef=%+.4f p=%.4f %s %s%s\n",
                r$rank, r$gene_id, r$wgcna_module, r$wgcna_kME,
                r$expected_direction, r$deseq2_WG_coef, r$deseq2_pvalue_cont,
                conc_str, old_str, sig_str))
  }
}

# Binomial test (fixed logic)
if (n_tested > 0) {
  binom_p <- binom.test(n_concordant, n_tested, p = 0.5, alternative = "greater")
  cat(sprintf("\n  Binomial test (fixed logic, H0: p=0.5): p = %.4f\n", binom_p$p.value))
  if (binom_p$p.value < 0.05) {
    cat("  → Hub genes ARE significantly concordant between WGCNA and DESeq2\n")
  } else {
    cat("  → Hub genes NOT significantly concordant (but may improve with n>20)\n")
  }
}

# ---- Module-level summary ----
if (nrow(module_lm_results) > 0) {
  n_mod_conc <- sum(module_lm_results$direction_concordant, na.rm = TRUE)
  n_mod <- nrow(module_lm_results)
  cat(sprintf("\n=== MODULE-LEVEL EIGENGENE LM RESULTS ===\n"))
  cat(sprintf("  Module eigengenes direction-concordant: %d / %d (%.1f%%)\n",
              n_mod_conc, n_mod, 100 * n_mod_conc / max(n_mod, 1)))
  
  mod_binom <- binom.test(n_mod_conc, n_mod, p = 0.5, alternative = "greater")
  cat(sprintf("  Binomial test: p = %.4f\n", mod_binom$p.value))
  
  cat("\n  Top growth-associated modules (by |LM WG coefficient|):\n")
  mod_sorted <- module_lm_results[order(-abs(module_lm_results$wg_coef)), ]
  for (i in seq_len(min(8, nrow(mod_sorted)))) {
    r <- mod_sorted[i, ]
    conc_str <- ifelse(isTRUE(r$direction_concordant), "✓", "✗")
    cat(sprintf("    %-6s LM_coef=%+.4f p=%.4f WGCNA_pr=%+.4f %s R²=%.3f\n",
                r$module, r$wg_coef, r$wg_pvalue, r$wgcna_partial_r, conc_str, r$r_squared))
  }
}

# ---- Write outputs ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.table(hub_validated, file.path(out_dir, "deseq2_cross_validation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\n  [OK] %s/deseq2_cross_validation.tsv\n", out_dir))

if (nrow(module_lm_results) > 0) {
  write.table(module_lm_results, file.path(out_dir, "deseq2_module_eigengene_test.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  [OK] %s/deseq2_module_eigengene_test.tsv\n", out_dir))
}

# Full genome-wide results
res_sorted <- res_df[order(-abs(res_df$log2FoldChange)), ]
write.table(res_sorted, file.path(out_dir, "deseq2_all_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  [OK] %s/deseq2_all_results.tsv (%d genes)\n", out_dir, nrow(res_sorted)))

cat("\nDone.\n")
