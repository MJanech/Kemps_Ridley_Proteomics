# ==============================================================================
# GSEA_collapse_redundant_pathways.R
#
# Purpose
#   Collapse redundant significant Reactome pathways from the completed
#   tie-robust nonparametric preranked GSEA using fgsea::collapsePathways().
#
# Relationship to the primary protein-level analysis
#   - limma-trend remains the PRIMARY protein-level differential-abundance test.
#   - This script operates only on completed GSEA pathway results.
#   - It does not rerun protein-level statistics and does not use limma.
#
# Method
#   fgsea::collapsePathways() identifies a smaller set of representative
#   "main pathways" by conditional enrichment/dependence testing.
#   This replaces the earlier Jaccard-similarity pathway reduction.
#
# Required existing folder in the working directory:
#   NP_ranked_fgsea_tie_robust
#
# Required R packages:
#   fgsea, msigdbr, data.table
# ==============================================================================

required <- c("fgsea", "msigdbr", "data.table")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Please install required package(s): ", paste(missing_pkgs, collapse=", "))
}

library(fgsea)
library(msigdbr)
library(data.table)

indir <- "NP_ranked_fgsea_tie_robust"
outdir <- "NP_ranked_fgsea_collapsed_main_pathways"
if (!dir.exists(indir)) {
  stop("Cannot find folder '", indir, "' in working directory: ", getwd())
}
dir.create(outdir, showWarnings=FALSE)

# Same Reactome collection used in the previous fgsea analysis.
cat("Loading MSigDB Reactome gene sets...\n")
reactome_df <- tryCatch(
  msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME"),
  error=function(e) msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
)
reactome_df <- reactome_df[!is.na(reactome_df$gene_symbol) & reactome_df$gene_symbol != "", ]
pathways <- split(reactome_df$gene_symbol, reactome_df$gs_name)
pathways <- lapply(pathways, unique)
cat("Reactome gene sets loaded:", length(pathways), "\n")

pretty_pathway <- function(x) {
  x <- sub("^REACTOME_", "", x)
  gsub("_", " ", x)
}

direction_notes <- c(
  S4="Positive NES = higher in OC4_0; negative NES = higher in OC2_0",
  S5="Positive NES = higher in OC4_1; negative NES = higher in OC2_1",
  S6="Positive NES = higher at OC2_1 release; negative NES = higher at OC2_0 admission",
  S7="Positive NES = higher at OC4_1 release; negative NES = higher at OC4_0 admission"
)

summary_parts <- list()
combined_main <- list()
combined_mapping <- list()

set.seed(20260822)

for (nm in c("S4","S5","S6","S7")) {
  rank_file <- file.path(indir, paste0(nm, "_tie_robust_ranked_genes.csv"))
  res_file  <- file.path(indir, paste0(nm, "_Reactome_tie_robust_fgsea_ALL.csv"))
  if (!file.exists(rank_file)) stop("Missing: ", rank_file)
  if (!file.exists(res_file)) stop("Missing: ", res_file)

  ranks_df <- read.csv(rank_file, stringsAsFactors=FALSE, check.names=FALSE)
  res <- read.csv(res_file, stringsAsFactors=FALSE, check.names=FALSE)

  if (!all(c("Gene", "NP_tie_robust_score") %in% names(ranks_df))) {
    stop("Expected Gene and NP_tie_robust_score columns in ", rank_file)
  }
  if (!all(c("pathway","pval","padj","ES","NES","size") %in% names(res))) {
    stop("Expected fgsea result columns missing in ", res_file)
  }

  ranks <- as.numeric(ranks_df$NP_tie_robust_score)
  names(ranks) <- ranks_df$Gene
  ok <- is.finite(ranks) & !is.na(names(ranks)) & names(ranks) != ""
  ranks <- sort(ranks[ok], decreasing=TRUE)

  sig <- res[!is.na(res$padj) & res$padj < 0.05, , drop=FALSE]
  sig <- sig[order(sig$pval, sig$padj), , drop=FALSE]

  cat("\n", nm, ": ", nrow(sig), " significant Reactome pathways at BH FDR < 0.05", sep="")

  if (nrow(sig) == 0) {
    cat("; nothing to collapse.\n")
    empty <- sig
    if (!"Pathway_label" %in% names(empty)) empty$Pathway_label <- character(0)
    if (!"Collapse_status" %in% names(empty)) empty$Collapse_status <- character(0)
    if (!"Parent_pathway" %in% names(empty)) empty$Parent_pathway <- character(0)
    write.csv(empty, file.path(outdir, paste0(nm, "_collapsed_main_pathways.csv")), row.names=FALSE)
    write.csv(empty, file.path(outdir, paste0(nm, "_collapse_mapping_ALL_significant.csv")), row.names=FALSE)

    summary_parts[[nm]] <- data.frame(
      Comparison=nm,
      Significant_before_collapse=0,
      Main_pathways_after_collapse=0,
      Redundant_pathways_collapsed=0,
      Positive_main=0,
      Negative_main=0,
      Direction_definition=direction_notes[[nm]],
      stringsAsFactors=FALSE
    )
    next
  }

  # collapsePathways expects an fgsea-style table. The complete significant
  # result table is passed; only same-direction pathways are reduced together.
  sig_dt <- as.data.table(sig)

  # pval.threshold=0.05 is a conditional enrichment/dependence criterion,
  # NOT a Jaccard overlap cutoff. nperm is increased above the package default
  # for a more stable collapse decision.
  collapsed <- fgsea::collapsePathways(
    fgseaRes=sig_dt,
    pathways=pathways,
    stats=ranks,
    pval.threshold=0.05,
    nperm=2000,
    gseaParam=1
  )

  main_ids <- collapsed$mainPathways
  parent_vec <- collapsed$parentPathways

  mapping <- sig
  mapping$Comparison <- nm
  mapping$Pathway_label <- pretty_pathway(mapping$pathway)
  mapping$Collapse_status <- ifelse(mapping$pathway %in% main_ids, "Main", "Collapsed_redundant")
  mapping$Parent_pathway <- NA_character_
  if (!is.null(parent_vec) && length(parent_vec) > 0) {
    pv <- as.character(parent_vec)
    names(pv) <- names(parent_vec)
    idx <- match(mapping$pathway, names(pv))
    mapping$Parent_pathway <- pv[idx]
  }
  mapping$Parent_pathway_label <- ifelse(
    is.na(mapping$Parent_pathway) | mapping$Parent_pathway == "",
    NA_character_,
    pretty_pathway(mapping$Parent_pathway)
  )
  mapping$Direction_definition <- direction_notes[[nm]]

  # Main pathways retain original fgsea statistics and leading-edge genes.
  main <- mapping[mapping$Collapse_status == "Main", , drop=FALSE]
  main <- main[order(main$padj, -abs(main$NES)), , drop=FALSE]

  write.csv(mapping,
            file.path(outdir, paste0(nm, "_collapse_mapping_ALL_significant.csv")),
            row.names=FALSE)
  write.csv(main,
            file.path(outdir, paste0(nm, "_collapsed_main_pathways.csv")),
            row.names=FALSE)

  cat(" -> ", nrow(main), " main pathways retained; ",
      nrow(sig)-nrow(main), " redundant pathways collapsed.\n", sep="")

  summary_parts[[nm]] <- data.frame(
    Comparison=nm,
    Significant_before_collapse=nrow(sig),
    Main_pathways_after_collapse=nrow(main),
    Redundant_pathways_collapsed=nrow(sig)-nrow(main),
    Positive_main=sum(main$NES > 0, na.rm=TRUE),
    Negative_main=sum(main$NES < 0, na.rm=TRUE),
    Direction_definition=direction_notes[[nm]],
    stringsAsFactors=FALSE
  )

  if (nrow(main) > 0) {
    tmp <- main
    combined_main[[nm]] <- tmp
  }
  combined_mapping[[nm]] <- mapping
}

summary_out <- do.call(rbind, summary_parts)
rownames(summary_out) <- NULL
write.csv(summary_out, file.path(outdir, "NP_fgsea_collapse_summary.csv"), row.names=FALSE)

if (length(combined_main) > 0) {
  main_all <- do.call(rbind, combined_main)
  rownames(main_all) <- NULL
  # Place comparison first for easier figure preparation.
  main_all <- main_all[, c("Comparison", setdiff(names(main_all), "Comparison")), drop=FALSE]
  write.csv(main_all, file.path(outdir, "ALL_comparisons_collapsed_main_pathways.csv"), row.names=FALSE)
}

if (length(combined_mapping) > 0) {
  map_all <- do.call(rbind, combined_mapping)
  rownames(map_all) <- NULL
  write.csv(map_all, file.path(outdir, "ALL_comparisons_collapse_mapping.csv"), row.names=FALSE)
}

cat("\nNONREDUNDANT PATHWAY COLLAPSE COMPLETE\n")
cat("=====================================\n")
print(summary_out, row.names=FALSE)
cat("\nMethod: fgsea::collapsePathways; conditional dependence threshold p = 0.05; nperm = 2000.\n")
cat("This replaces the earlier Jaccard-similarity reduction for the nonparametric preranked pathway analysis.\n")
cat("No limma statistics were used in pathway collapsing.\n")
cat("Results folder:", normalizePath(outdir), "\n")
