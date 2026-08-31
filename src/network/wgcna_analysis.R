#!/usr/bin/env Rscript
# ================================================================
# MEMORY-SAFE WGCNA — Track 3 Host-Microbe Integration
# ================================================================
# Designed for WSL/Ubuntu with 7.7 GB RAM.
# Uses blockwise network construction to stay within memory budget.
#
# Usage:
#   # Standard: runs on merged counts
#   Rscript src/network/wgcna_analysis.R \
#     --counts merged_counts.tsv \
#     --metadata metadata.tsv \
#     --out-hub hub_genes.tsv \
#     --out-modules modules.csv \
#     --out-me me_matrix.tsv \
#     --max-block-size 8000 \
#     --min-counts 10 \
#     --threads 2
#
#   # With confound correction (v2): runs on corrected counts
#   Rscript src/network/wgcna_analysis.R \
#     --counts corrected_counts.tsv \
#     --metadata metadata.tsv \
#     --out-hub hub_genes_corrected.tsv \
#     --out-modules modules_corrected.csv \
#     --out-me me_matrix_corrected.tsv \
#     --corrected TRUE
# ================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(optparse)
})

# --- Parse arguments ---
option_list <- list(
  make_option("--counts", type="character", help="Path to merged counts TSV"),
  make_option("--metadata", type="character", help="Path to metadata TSV"),
  make_option("--out-hub", type="character", help="Output hub genes"),
  make_option("--out-modules", type="character", help="Output module assignments"),
  make_option("--out-me", type="character", help="Output module eigengenes"),
  make_option("--max-block-size", type="integer", default=8000,
              help="Max genes per block [default: 8000]"),
  make_option("--min-counts", type="integer", default=10,
              help="Min mean counts to retain a gene [default: 10]"),
  make_option("--threads", type="integer", default=2,
              help="Number of threads [default: 2]"),
  make_option("--corrected", type="character", default="FALSE",
              help="Whether input counts have been confound-corrected [default: FALSE]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# Log whether using confound-corrected counts
is_corrected <- toupper(opts$corrected) %in% c("TRUE", "YES", "1")
cat(sprintf("[WGCNA] Input type: %s\n",
            ifelse(is_corrected, "CONFOUND-CORRECTED", "STANDARD")))
if (is_corrected) {
  cat("[WGCNA] NOTE: Running on confound-corrected expression matrix.\n")
  cat("[WGCNA]       Modules are expected to be less driven by sex/tissue confounds.\n")
}

# --- Guardrail: check available memory ---
meminfo <- system("grep MemAvailable /proc/meminfo | awk '{print $2}'", intern=TRUE)
avail_mb <- as.numeric(meminfo) / 1024
cat(sprintf("[WGCNA] Available RAM: %.0f MB\n", avail_mb))
if (avail_mb < 2000) {
  stop(sprintf("Insufficient RAM (%.0f MB < 2000 MB). Free memory and retry.", avail_mb))
}

# Enable multi-threading
allowWGCNAThreads(nThreads = opts$threads)

# --- Load counts ---
cat("[WGCNA] Loading counts from:", opts$counts, "\n")
counts <- read.table(opts$counts, header=TRUE, row.names=1, check.names=FALSE)

# --- Filter low-count genes (saves RAM) ---
min_counts_val <- opts$`min-counts`
mean_counts <- rowMeans(counts)
keep <- mean_counts >= min_counts_val
cat(sprintf("[WGCNA] Filtering genes: %d / %d pass min-counts >= %d\n",
            sum(keep), nrow(counts), min_counts_val))
counts <- counts[keep, ]

if (nrow(counts) < 100) {
  stop("Too few genes after filtering. Lower --min-counts or check input.")
}

# --- Transpose: WGCNA expects genes as columns, samples as rows ---
datExpr <- t(counts)

# --- Check for missing values ---
if (any(is.na(datExpr))) {
  cat("[WGCNA] Warning: NA values found. Imputing with column medians.\n")
  datExpr <- apply(datExpr, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm=TRUE); x
  })
}

# --- Soft threshold power ---
cat("[WGCNA] Picking soft threshold...\n")
powers <- c(1:20, seq(22, 30, 2))
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=2, networkType="signed")

# Pick power at R^2 >= 0.8, or default to 6
if (!is.na(sft$powerEstimate)) {
  softPower <- sft$powerEstimate
} else {
  softPower <- 6
  cat(sprintf("[WGCNA] No power reached R^2=0.8, using default power=%d\n", softPower))
}
cat(sprintf("[WGCNA] Soft threshold power: %d\n", softPower))

# --- Blockwise network (MEMORY-SAFE) ---
cat(sprintf("[WGCNA] Building blockwise network (maxBlockSize=%d)...\n", opts$`max-block-size`))
cat("[WGCNA] This is the most memory-intensive step. May take 15-60 min.\n")

# Memory guard: check RAM before blockwise
gc(reset=TRUE)
meminfo2 <- system("grep MemAvailable /proc/meminfo | awk '{print $2}'", intern=TRUE)
avail_mb2 <- as.numeric(meminfo2) / 1024
cat(sprintf("[WGCNA] RAM before blockwise: %.0f MB\n", avail_mb2))

net <- blockwiseModules(
  datExpr,
  power = softPower,
  maxBlockSize = opts$`max-block-size`,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,       # Don't save TOM to disk (saves I/O)
  saveTOMFileBase = NULL,
  verbose = 3,
  nThreads = opts$threads,
  deepSplit = 2,
  detectCutHeight = 0.995
)

cat(sprintf("[WGCNA] Detected %d modules\n", length(unique(net$colors))))

# --- Module Eigengenes ---
MEs <- net$MEs
colnames(MEs) <- paste0("ME", seq_len(ncol(MEs)))

# --- Load trait data ---
cat("[WGCNA] Loading metadata:", opts$metadata, "\n")
meta <- read.table(opts$metadata, header=TRUE, row.names=1, check.names=FALSE, sep="\t")

# Align samples
common_samples <- intersect(rownames(datExpr), rownames(meta))
cat(sprintf("[WGCNA] Common samples: %d\n", length(common_samples)))
if (length(common_samples) < 10) {
  stop("Too few common samples. Check sample IDs match between counts and metadata.")
}
MEs <- MEs[common_samples, , drop=FALSE]
meta <- meta[common_samples, , drop=FALSE]

# --- Module-Trait Correlations ---
cat("[WGCNA] Computing module-trait correlations...\n")
module_trait_cor <- cor(MEs, meta, use="pairwise.complete.obs")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(MEs))

# --- Identify hub genes (highest module membership per module) ---
cat("[WGCNA] Identifying hub genes...\n")
gene_module_colors <- net$colors
names(gene_module_colors) <- colnames(datExpr)

hub_genes <- data.frame(
  gene = character(),
  module = character(),
  kME = numeric(),
  trait_cor = numeric(),
  stringsAsFactors = FALSE
)

for (mod in setdiff(unique(gene_module_colors), 0)) {  # 0 = unassigned
  mod_genes <- names(gene_module_colors)[gene_module_colors == mod]
  if (length(mod_genes) < 3) next

  # Module membership
  mod_me <- MEs[, paste0("ME", mod)]
  kME <- cor(datExpr[, mod_genes, drop=FALSE], mod_me, use="pairwise.complete.obs")
  kME <- kME[, 1]

  # Top 10% as hubs
  top_n <- max(3, ceiling(length(mod_genes) * 0.1))
  top_idx <- order(abs(kME), decreasing=TRUE)[1:min(top_n, length(kME))]

  for (i in top_idx) {
    hub_genes <- rbind(hub_genes, data.frame(
      gene = mod_genes[i],
      module = paste0("M", mod),
      kME = round(kME[i], 4),
      trait_cor = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
}

cat(sprintf("[WGCNA] Identified %d hub genes across %d modules\n",
            nrow(hub_genes), length(unique(hub_genes$module))))

# --- Write outputs ---
write.table(hub_genes, opts$`out-hub`, sep="\t", row.names=FALSE, quote=FALSE)
cat("[WGCNA] Hub genes written to:", opts$`out-hub`, "\n")

module_df <- data.frame(
  gene = names(gene_module_colors),
  module = paste0("M", gene_module_colors),
  stringsAsFactors = FALSE
)
write.table(module_df, opts$`out-modules`, sep="\t", row.names=FALSE, quote=FALSE)
cat("[WGCNA] Module assignments written to:", opts$`out-modules`, "\n")

write.table(MEs, opts$`out-me`, sep="\t", row.names=TRUE, quote=FALSE)
cat("[WGCNA] Eigengenes written to:", opts$`out-me`, "\n")

cat("[WGCNA] Done.\n")
