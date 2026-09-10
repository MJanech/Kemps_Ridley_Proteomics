# Reviewer validation: tie-robust NP-ranked Reactome fgsea
# Uses ONLY original nonparametric framework: Mann-Whitney U (S4/S5)
# and paired Wilcoxon signed-rank (S6/S7). No limma statistics are used.
#
# Purpose: remove arbitrary ordering of tied preranked statistics in fgsea.
# Primary rank score is unchanged:
#   rank-biserial correlation * [-log10(original NP p-value)]
# Exact ties in that primary score are broken deterministically using a
# Hodges-Lehmann nonparametric location-shift estimate:
#   S4/S5: median of all pairwise (B - A) differences
#   S6/S7: Wilcoxon paired Hodges-Lehmann pseudomedian (Walsh averages)
# A tiny numerical offset is then added ONLY within tied primary-score groups,
# small enough that it cannot change the ordering of distinct primary scores.
# fgsea is run serially to avoid Windows/RStudio parallel serialization warnings.

required <- c("readxl", "fgsea", "msigdbr", "BiocParallel")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Please install required package(s): ", paste(missing_pkgs, collapse=", "))
}

library(readxl)
library(fgsea)
library(msigdbr)
library(BiocParallel)

args <- commandArgs(trailingOnly = TRUE)
xlsx <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  "Supplemental_Tables_1_through_8.xlsx"
}
if (!file.exists(xlsx)) {
  stop("Cannot find ", xlsx, " in working directory: ", getwd())
}

outdir <- "NP_ranked_fgsea_tie_robust"
dir.create(outdir, showWarnings=FALSE)

trim <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}

hl_independent <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  median(as.vector(outer(b, a, "-")), na.rm=TRUE)
}

hl_paired <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b)
  d <- b[ok] - a[ok]
  if (length(d) == 0) return(NA_real_)
  w <- outer(d, d, "+") / 2
  median(w[upper.tri(w, diag=TRUE)], na.rm=TRUE)
}

read_np_comp <- function(sheet, skip, nA, nB, labelA, labelB, test_name, paired=FALSE) {
  d <- read_excel(xlsx, sheet=sheet, skip=skip, .name_repair="unique")
  names(d)[1:5] <- c("Row_Number","Protein","Gene","Accession","MW")
  d$Gene <- trim(d$Gene)
  d$Accession <- trim(d$Accession)
  d$Protein <- trim(d$Protein)

  Acols <- 6:(5+nA)
  Bcols <- (6+nA):(5+nA+nB)
  A <- as.data.frame(lapply(d[Acols], as.numeric))
  B <- as.data.frame(lapply(d[Bcols], as.numeric))

  pcol <- grep("^p\\.value$", names(d))
  qcol <- grep("^padj_BH$", names(d))
  rcol <- grep("Rank-biserial Correlation", names(d), fixed=TRUE)
  if (length(pcol) != 1) stop("Could not uniquely identify p.value in ", sheet)
  if (length(qcol) != 1) stop("Could not uniquely identify padj_BH in ", sheet)
  if (length(rcol) != 1) stop("Could not uniquely identify rank-biserial correlation in ", sheet)

  d$p_raw <- as.numeric(d[[pcol]])
  d$padj_BH_original <- as.numeric(d[[qcol]])
  d$rank_biserial <- as.numeric(d[[rcol]])

  safe_p <- pmax(d$p_raw, .Machine$double.xmin, na.rm=FALSE)
  d$NP_primary_score <- d$rank_biserial * (-log10(safe_p))

  # Nonparametric tie-break effect estimate computed directly from the same data.
  if (paired) {
    d$HL_shift_B_minus_A <- vapply(seq_len(nrow(d)), function(i) {
      hl_paired(unlist(A[i, ], use.names=FALSE), unlist(B[i, ], use.names=FALSE))
    }, numeric(1))
  } else {
    d$HL_shift_B_minus_A <- vapply(seq_len(nrow(d)), function(i) {
      hl_independent(unlist(A[i, ], use.names=FALSE), unlist(B[i, ], use.names=FALSE))
    }, numeric(1))
  }

  d$Test <- test_name
  d$Comparison <- sheet
  d$Label_A <- labelA
  d$Label_B <- labelB
  d
}

S4 <- read_np_comp("Table S4", skip=2, nA=9,  nB=15,
                   labelA="OC2_0", labelB="OC4_0", test_name="Mann-Whitney U", paired=FALSE)
S5 <- read_np_comp("Table S5", skip=2, nA=9,  nB=15,
                   labelA="OC2_1", labelB="OC4_1", test_name="Mann-Whitney U", paired=FALSE)
S6 <- read_np_comp("Table S6", skip=3, nA=9,  nB=9,
                   labelA="OC2_0", labelB="OC2_1", test_name="paired Wilcoxon signed-rank", paired=TRUE)
S7 <- read_np_comp("Table S7", skip=3, nA=15, nB=15,
                   labelA="OC4_0", labelB="OC4_1", test_name="paired Wilcoxon signed-rank", paired=TRUE)

comps <- list(S4=S4, S5=S5, S6=S6, S7=S7)

direction_notes <- c(
  S4="Positive NES = higher in OC4_0; negative NES = higher in OC2_0",
  S5="Positive NES = higher in OC4_1; negative NES = higher in OC2_1",
  S6="Positive NES = higher at OC2_1 release; negative NES = higher at OC2_0 admission",
  S7="Positive NES = higher at OC4_1 release; negative NES = higher at OC4_0 admission"
)

# Collapse multiple accessions to one gene symbol using the strongest absolute
# primary NP score. If tied, use smaller p, then larger absolute HL shift.
collapse_to_gene <- function(nm, d) {
  z <- d[!is.na(d$Gene) & is.finite(d$NP_primary_score) & is.finite(d$p_raw) & is.finite(d$rank_biserial),
         c("Gene","Protein","Accession","p_raw","padj_BH_original",
           "rank_biserial","NP_primary_score","HL_shift_B_minus_A")]
  z$Abs_primary <- abs(z$NP_primary_score)
  z$Abs_HL <- abs(z$HL_shift_B_minus_A)
  z$Abs_HL[!is.finite(z$Abs_HL)] <- -Inf
  z <- z[order(-z$Abs_primary, z$p_raw, -z$Abs_HL, z$Accession), , drop=FALSE]
  z <- z[!duplicated(z$Gene), , drop=FALSE]
  z
}

# Break ONLY exact ties in primary score. Distinct primary scores never cross.
make_tie_robust_scores <- function(z) {
  primary <- z$NP_primary_score
  hl <- z$HL_shift_B_minus_A
  genes <- z$Gene

  uv <- sort(unique(primary[is.finite(primary)]))
  gaps <- diff(uv)
  gaps <- gaps[is.finite(gaps) & gaps > 0]
  eps <- if (length(gaps) > 0) min(gaps) / 1000 else 1e-10

  final <- primary
  tied_vals <- unique(primary[duplicated(primary) | duplicated(primary, fromLast=TRUE)])

  for (v in tied_vals) {
    idx <- which(primary == v)
    if (length(idx) < 2) next
    # Larger HL shift (B-A) ranks higher within a tied NP score.
    # Gene symbol is a deterministic last-resort tiebreaker only if HL is also tied.
    ord <- order(-hl[idx], genes[idx], na.last=TRUE)
    idx_ord <- idx[ord]
    offsets <- seq(0.4, -0.4, length.out=length(idx_ord)) * eps
    final[idx_ord] <- primary[idx_ord] + offsets
  }

  # Defensive final uniqueness check: if exact ties somehow remain, add a far smaller
  # deterministic offset by gene name without altering any meaningful magnitude.
  if (anyDuplicated(final)) {
    eps2 <- eps / 1e6
    vals2 <- unique(final[duplicated(final) | duplicated(final, fromLast=TRUE)])
    for (v in vals2) {
      idx <- which(final == v)
      ord <- order(genes[idx])
      idx_ord <- idx[ord]
      offsets <- seq(0.4, -0.4, length.out=length(idx_ord)) * eps2
      final[idx_ord] <- final[idx_ord] + offsets
    }
  }

  list(score=final, epsilon=eps)
}

cat("Loading MSigDB Reactome gene sets...\n")
reactome_df <- tryCatch(
  msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME"),
  error=function(e) msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
)
reactome_df <- reactome_df[!is.na(reactome_df$gene_symbol) & reactome_df$gene_symbol != "", ]
pathways <- split(reactome_df$gene_symbol, reactome_df$gs_name)
pathways <- lapply(pathways, unique)
cat("Reactome gene sets loaded:", length(pathways), "\n")

set.seed(20260822)
summary_parts <- list()

for (nm in names(comps)) {
  z <- collapse_to_gene(nm, comps[[nm]])
  tieinfo <- make_tie_robust_scores(z)
  z$NP_tie_robust_score <- tieinfo$score
  z <- z[order(z$NP_tie_robust_score, decreasing=TRUE), , drop=FALSE]

  before_ties <- sum(duplicated(z$NP_primary_score[z$NP_primary_score != 0]))
  before_pct <- 100 * before_ties / nrow(z)
  after_ties <- sum(duplicated(z$NP_tie_robust_score[z$NP_tie_robust_score != 0]))

  write.csv(z, file.path(outdir, paste0(nm, "_tie_robust_ranked_genes.csv")), row.names=FALSE)

  ranks <- z$NP_tie_robust_score
  names(ranks) <- z$Gene
  ranks <- sort(ranks, decreasing=TRUE)

  cat("Running tie-robust NP-ranked fgsea for", nm,
      "with", length(ranks), "genes; original tied nonzero stats =",
      sprintf("%.2f%%", before_pct), "...\n")

  res <- fgsea(
    pathways=pathways,
    stats=ranks,
    minSize=5,
    maxSize=500,
    eps=0,
    BPPARAM=BiocParallel::SerialParam()
  )
  res <- as.data.frame(res)
  if ("leadingEdge" %in% names(res)) {
    res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse=";"), character(1))
  }
  if (nrow(res) > 0) {
    res <- res[order(res$padj, res$pval), , drop=FALSE]
    res$Direction <- ifelse(res$NES > 0, "Positive", "Negative")
    res$Direction_definition <- direction_notes[[nm]]
  }

  write.csv(res, file.path(outdir, paste0(nm, "_Reactome_tie_robust_fgsea_ALL.csv")), row.names=FALSE)
  sig <- res[!is.na(res$padj) & res$padj < 0.05, , drop=FALSE]
  write.csv(sig, file.path(outdir, paste0(nm, "_Reactome_tie_robust_fgsea_FDR05.csv")), row.names=FALSE)
  write.csv(head(res, 30), file.path(outdir, paste0(nm, "_Reactome_tie_robust_fgsea_TOP30.csv")), row.names=FALSE)

  summary_parts[[nm]] <- data.frame(
    Comparison=nm,
    NP_test=unique(comps[[nm]]$Test)[1],
    Ranked_unique_genes=length(ranks),
    Original_tied_nonzero_stats=before_ties,
    Original_tied_nonzero_percent=round(before_pct, 2),
    Ties_after_HL_tiebreak=after_ties,
    Tiebreak_epsilon=tieinfo$epsilon,
    Pathways_tested=nrow(res),
    FDR_lt_0_05=nrow(sig),
    FDR_positive=sum(sig$NES > 0, na.rm=TRUE),
    FDR_negative=sum(sig$NES < 0, na.rm=TRUE),
    Min_FDR=if (nrow(res)>0 && any(is.finite(res$padj))) min(res$padj, na.rm=TRUE) else NA_real_,
    Direction_definition=direction_notes[[nm]],
    stringsAsFactors=FALSE
  )
}

summary_out <- do.call(rbind, summary_parts)
rownames(summary_out) <- NULL
write.csv(summary_out, file.path(outdir, "NP_tie_robust_fgsea_summary.csv"), row.names=FALSE)

cat("\nTIE-ROBUST NONPARAMETRIC FGSEA COMPLETE\n")
cat("======================================\n")
print(summary_out, row.names=FALSE)
cat("\nPrimary ranking metric remains:\n")
cat("  rank-biserial correlation * [-log10(original Mann-Whitney/Wilcoxon p-value)]\n")
cat("Exact ties were broken ONLY with a Hodges-Lehmann nonparametric effect estimate.\n")
cat("No limma statistics were used. fgsea was run serially.\n")
cat("Results folder:", normalizePath(outdir), "\n")
