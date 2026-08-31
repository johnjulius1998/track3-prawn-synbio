#!/usr/bin/env Rscript
###############################################################################
# m6_gsea.R — Formal GSEA / ORA on the M6 module
# ==============================================
# Track 3 Host-Microbe Integration
#
# Replaces the keyword-guess annotation (m6_gene_annotation.py) with a formal
# gene-set enrichment analysis:
#
#   ranked lists:
#     (1) kME  — module membership of ALL expressed genes vs the M6 eigengene,
#               recomputed in the WGCNA INPUT SPACE (merged_counts.tsv,
#               mean-counts >= 10; verified to reproduce the stored M6 hub
#               kME with cor = 1.0000)
#     (2) DESeq2 — signed genome-wide weight-gain statistic
#               (results/reports/deseq2_all_results.tsv)
#
#   gene sets (organism-agnostic, from Swiss-Prot orthology mapping):
#     - KEGG pathways  data/interim/literature_tables/uniprot_gene_sets_kegg.tsv
#     - GO BP terms    data/interim/literature_tables/uniprot_gene_sets_go.tsv
#
#   tests:
#     - clusterProfiler::GSEA on both rankings x both set types
#     - clusterProfiler::enricher ORA on the 72 M6 hub genes (sanity check)
#
# OUTPUTS (results/reports/, results/figures/):
#   m6_gsea_kegg_kme.tsv      m6_gsea_kegg_deseq2.tsv
#   m6_gsea_go_kme.tsv        m6_gsea_go_deseq2.tsv
#   m6_ora_hubs_kegg.tsv      m6_ora_hubs_go.tsv
#   figures/m6_gsea_*.png
#
# USAGE: Rscript src/network/m6_gsea.R [BASE_DIR]
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[[1]] else "."

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(GO.db)
})

msg <- function(...) cat(sprintf("[GSEA] %s\n", sprintf(...)))

out_reports <- file.path(BASE, "results", "reports")
out_figs <- file.path(BASE, "results", "figures")
dir.create(out_reports, showWarnings = FALSE, recursive = TRUE)
dir.create(out_figs, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Expression data + eigengenes
# ---------------------------------------------------------------------------
# The M6 eigengene (ME6) comes from the canonical WGCNA run, which was built
# on merged_counts.tsv filtered to mean counts >= 10 (see wgcna_analysis.R).
# kME is therefore recomputed in THAT space — verified to reproduce the
# stored hub kME exactly (cor = 1.0000). corrected_tpm.tsv belongs to a
# different (confound-corrected) pipeline and must NOT be mixed in here.
counts_path <- file.path(BASE, "data", "processed", "gene_expression", "merged_counts.tsv")
me_path <- file.path(BASE, "data", "processed", "wgcna", "me_matrix.tsv")
hub_path <- file.path(BASE, "data", "processed", "wgcna", "hub_genes.tsv")

msg("loading WGCNA input counts: %s", counts_path)
counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
keep <- rowMeans(counts) >= 10
msg("genes passing mean-counts >= 10 filter: %d / %d", sum(keep), nrow(counts))
counts <- counts[keep, , drop = FALSE]

msg("loading eigengene matrix: %s", me_path)
me <- read.delim(me_path, row.names = 1, check.names = FALSE)
stopifnot("ME6" %in% colnames(me))

common <- intersect(rownames(me), colnames(counts))
msg("%d samples shared between counts columns and eigengene rows", length(common))
me <- me[common, , drop = FALSE]
counts <- counts[, common, drop = FALSE]
me6 <- as.numeric(me[, "ME6"])
names(me6) <- rownames(me)

# ---------------------------------------------------------------------------
# 2. Ranked list 1: full-module kME vs ME6 (all genes, not just hubs)
# ---------------------------------------------------------------------------
msg("recomputing kME for all %d genes vs ME6 (WGCNA input space)", nrow(counts))
cnt_mat <- as.matrix(counts)
kme_all <- suppressWarnings(
  apply(cnt_mat, 1, function(r) cor(r, me6, use = "complete.obs")))
kme_all <- kme_all[is.finite(kme_all)]
kme_rank <- sort(kme_all, decreasing = TRUE)
msg("kME ranking: %d genes (range %.3f .. %.3f)", length(kme_rank),
    tail(kme_rank, 1), head(kme_rank, 1))

# cross-check against the stored hub kME (should be ~1.0 for M6)
if (file.exists(hub_path)) {
  hubs0 <- read.delim(hub_path)
  m6h <- hubs0[hubs0$module == "M6", ]
  cc <- cor(m6h$kME, kme_all[m6h$gene], use = "complete.obs")
  msg("sanity check: cor(stored M6 hub kME, recomputed) = %.4f", cc)
}

# ---------------------------------------------------------------------------
# 3. Ranked list 2: DESeq2 signed statistic (weight-gain effect)
# ---------------------------------------------------------------------------
deseq_path <- file.path(BASE, "results", "reports", "deseq2_all_results.tsv")
if (file.exists(deseq_path)) {
  deseq <- read.delim(deseq_path)
  if (!"gene_id" %in% colnames(deseq) && "gene" %in% colnames(deseq)) {
    colnames(deseq)[colnames(deseq) == "gene"] <- "gene_id"
  }
  deseq_rank <- setNames(as.numeric(deseq$stat), deseq$gene_id)
  deseq_rank <- sort(deseq_rank[is.finite(deseq_rank)], decreasing = TRUE)
  msg("DESeq2 ranking: %d genes", length(deseq_rank))
} else {
  stop("DESeq2 results not found: ", deseq_path)
}

# ---------------------------------------------------------------------------
# 4. Gene sets (Swiss-Prot orthology mapping)
# ---------------------------------------------------------------------------
load_sets <- function(kind) {
  p <- file.path(BASE, "data", "interim", "literature_tables",
                 sprintf("uniprot_gene_sets_%s.tsv", kind))
  if (!file.exists(p)) return(NULL)
  d <- read.delim(p, stringsAsFactors = FALSE)
  if (!nrow(d)) return(NULL)
  list(
    t2g = data.frame(term = d[[2]], gene = d[[1]], stringsAsFactors = FALSE),
    t2n = if (ncol(d) >= 3) unique(data.frame(term = d[[2]], name = d[[3]],
                                               stringsAsFactors = FALSE)) else NULL
  )
}

kegg_sets <- load_sets("kegg")
go_sets_raw <- load_sets("go")
go_sets <- NULL
if (!is.null(go_sets_raw)) {
  go_names <- tryCatch({
    ids <- unique(go_sets_raw$t2g$term)
    nms <- suppressMessages(AnnotationDbi::Term(GO.db::GOTERM[ids]))
    nms[is.na(nms)] <- ids[is.na(nms)]
    data.frame(term = ids, name = as.character(nms), stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  go_sets <- list(t2g = go_sets_raw$t2g, t2n = go_names)
}

# ---------------------------------------------------------------------------
# 5. GSEA runs
# ---------------------------------------------------------------------------
run_gsea <- function(ranked, sets, label, min_gs = 10, max_gs = 500) {
  out_file <- file.path(out_reports, paste0("m6_gsea_", label, ".tsv"))
  tryCatch({
    g <- clusterProfiler::GSEA(
      geneList = ranked,
      TERM2GENE = sets$t2g,
      TERM2NAME = sets$t2n,
      minGSSize = min_gs,
      maxGSSize = max_gs,
      pvalueCutoff = 1.0,
      pAdjustMethod = "BH",
      seed = TRUE,
      verbose = FALSE
    )
    res <- as.data.frame(g)
    if (nrow(res)) res <- res[order(res$pvalue), ]
    write.table(res, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
    msg("%s: %d enriched sets -> %s", label, nrow(res), out_file)
    if (nrow(res)) {
      tryCatch({
        p <- suppressMessages(enrichplot::dotplot(g, showCategory = min(20, nrow(res))))
        png(file.path(out_figs, paste0("m6_gsea_dotplot_", label, ".png")),
            width = 1400, height = 1100, res = 150)
        print(p)
        dev.off()
      }, error = function(e) msg("%s dotplot failed: %s", label, conditionMessage(e)))
    }
    invisible(res)
  }, error = function(e) {
    msg("%s FAILED: %s (writing empty table)", label, conditionMessage(e))
    write.table(data.frame(), out_file, sep = "\t", row.names = FALSE)
    invisible(NULL)
  })
}

if (!is.null(kegg_sets)) {
  run_gsea(kme_rank, kegg_sets, "kegg_kme")
  run_gsea(deseq_rank, kegg_sets, "kegg_deseq2")
} else {
  msg("KEGG gene sets missing — skipping KEGG GSEA")
}
if (!is.null(go_sets)) {
  run_gsea(kme_rank, go_sets, "go_kme")
  run_gsea(deseq_rank, go_sets, "go_deseq2")
} else {
  msg("GO gene sets missing — skipping GO GSEA")
}

# ---------------------------------------------------------------------------
# 6. ORA sanity check on the 72 M6 hub genes
# ---------------------------------------------------------------------------
msg("ORA on M6 hub genes")
if (file.exists(hub_path)) {
  hubs <- read.delim(hub_path)
  if ("module" %in% colnames(hubs)) {
    hub_genes <- hubs$gene[hubs$module == "M6"]
  } else {
    hub_genes <- hubs$gene
  }
  msg("M6 hub genes: %d", length(hub_genes))

  run_ora <- function(genes, universe, sets, label) {
    out_file <- file.path(out_reports, paste0("m6_ora_hubs_", label, ".tsv"))
    tryCatch({
      en <- clusterProfiler::enricher(
        gene = genes,
        universe = universe,
        TERM2GENE = sets$t2g,
        TERM2NAME = sets$t2n,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1.0,
        pAdjustMethod = "BH"
      )
      res <- as.data.frame(en)
      if (nrow(res)) res <- res[order(res$pvalue), ]
      write.table(res, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
      msg("ORA %s: %d sets -> %s", label, nrow(res), out_file)
      invisible(res)
    }, error = function(e) {
      msg("ORA %s FAILED: %s", label, conditionMessage(e))
      write.table(data.frame(), out_file, sep = "\t", row.names = FALSE)
      invisible(NULL)
    })
  }

  universe <- intersect(names(kme_rank), rownames(counts))
  if (!is.null(kegg_sets)) {
    kegg_universe <- intersect(universe, unique(kegg_sets$t2g$gene))
    run_ora(intersect(hub_genes, kegg_universe), kegg_universe, kegg_sets, "kegg")
  }
  if (!is.null(go_sets)) {
    go_universe <- intersect(universe, unique(go_sets$t2g$gene))
    run_ora(intersect(hub_genes, go_universe), go_universe, go_sets, "go")
  }
} else {
  msg("hub_genes.tsv missing — skipping ORA")
}

msg("DONE")
