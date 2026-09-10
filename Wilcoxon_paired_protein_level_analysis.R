# ==============================================================================
# Wilcoxon_paired_protein_level_analysis.R
#
# Purpose
#   Reproduce the paired protein-level Wilcoxon signed-rank analyses reported in:
#     Table S6: OC2 admission vs release
#     Table S7: OC4 admission vs release
#
# Input
#   Supplemental_Tables_1_through_8.xlsx
#   Sheet "Table S3" contains the study-wide weighted spectral counts (wSC).
#
# Statistical framework
#   - Two-sided paired Wilcoxon signed-rank test
#   - Exact conditional p-values are calculated with exactRankTests::wilcox.exact,
#     matching the original analysis used for Tables S6 and S7.
#   - Zero wSC values are retained as non-detections; no abundance values are
#     imputed.
#   - Comparison-specific filtering follows Table S1:
#       OC2: proteins with >12 zeros across the 18 paired values are excluded.
#       OC4: proteins with >24 zeros across the 30 paired values are excluded.
#   - Benjamini-Hochberg correction is applied separately within each comparison.
#   - Paired differences and effect-size direction are defined as:
#       release - admission
#     Therefore, positive values indicate higher abundance at release and
#     negative values indicate higher abundance at admission.
#
# Expected validation
#   Table S6: 350 proteins retained;  87 proteins at BH FDR < 0.05
#   Table S7: 347 proteins retained; 153 proteins at BH FDR < 0.05
#
# Usage
#   Place this script and the Supplemental Tables workbook in the working
#   directory and run:
#       source("Wilcoxon_paired_protein_level_analysis.R")
#
#   An alternative workbook path can be supplied from the command line:
#       Rscript Wilcoxon_paired_protein_level_analysis.R path/to/workbook.xlsx
# ==============================================================================

required_packages <- c("readxl", "exactRankTests")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required. Install it with:\n",
      "install.packages('", pkg, "')"
    )
  }
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "Supplemental_Tables_1_through_8.xlsx"
output_dir <- "Wilcoxon_paired_results"

if (!file.exists(input_file)) {
  stop("Input workbook not found: ", input_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Paired sample definitions used in the manuscript
# ------------------------------------------------------------------------------

OC2_IDS <- c("114", "115", "116", "133", "137", "162", "192", "249", "278")

OC4_IDS <- c(
  "081", "087", "088", "098", "100", "107", "108", "110",
  "125", "193", "203", "204", "254", "265", "269"
)

OC2_ADMISSION <- paste("Sample", paste0(OC2_IDS, "A"))
OC2_RELEASE   <- paste("Sample", paste0(OC2_IDS, "B"))

OC4_ADMISSION <- paste("Sample", paste0(OC4_IDS, "A"))
OC4_RELEASE   <- paste("Sample", paste0(OC4_IDS, "B"))

MAX_ZEROS_OC2 <- 12L  # Table S1: 12/18 zeros permitted
MAX_ZEROS_OC4 <- 24L  # Table S1: 24/30 zeros permitted

# ------------------------------------------------------------------------------
# Read the study-wide wSC matrix from Table S3
# ------------------------------------------------------------------------------

# Row 4 of the Excel worksheet contains the variable names, so skip the title
# and the two group-label rows above it.
wsc <- readxl::read_excel(input_file, sheet = "Table S3", skip = 3)
names(wsc) <- trimws(names(wsc))

required_columns <- c(
  "Row", "Identified Proteins", "Gene Name", "Accession Number",
  "Molecular Weight (kDa)",
  OC2_ADMISSION, OC2_RELEASE, OC4_ADMISSION, OC4_RELEASE
)

missing_columns <- setdiff(required_columns, names(wsc))

if (length(missing_columns) > 0) {
  stop(
    "The following required columns were not found in Table S3:\n",
    paste(missing_columns, collapse = "\n")
  )
}

if (nrow(wsc) != 488L) {
  stop("Expected 488 proteins in Table S3, but found ", nrow(wsc), ".")
}

metadata_columns <- c(
  "Row", "Identified Proteins", "Gene Name",
  "Accession Number", "Molecular Weight (kDa)"
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

numeric_matrix <- function(data, columns) {
  x <- as.matrix(data[, columns, drop = FALSE])
  storage.mode(x) <- "numeric"

  if (anyNA(x)) {
    stop("Unexpected missing/non-numeric abundance value found in sample columns.")
  }

  x
}

signed_fold_change <- function(numerator, denominator) {
  # Preserve direction while expressing decreases as negative reciprocal FC.
  # Example: ratio 0.5 -> -2; ratio 2 -> +2.
  if (is.na(numerator) || is.na(denominator)) return(NA_real_)
  if (numerator == 0 && denominator == 0) return(NA_real_)
  if (denominator == 0) return(Inf)
  if (numerator == 0) return(-Inf)

  ratio <- numerator / denominator
  if (ratio >= 1) ratio else -1 / ratio
}

paired_rank_biserial <- function(difference) {
  # Rank-biserial correlation for paired data.
  # Zero differences are omitted from the signed-rank calculation.
  difference <- difference[difference != 0]

  if (length(difference) == 0) {
    return(NA_real_)
  }

  ranks <- rank(abs(difference), ties.method = "average")
  w_positive <- sum(ranks[difference > 0])
  w_negative <- sum(ranks[difference < 0])

  (w_positive - w_negative) / (w_positive + w_negative)
}

run_paired_wilcoxon <- function(
  data,
  admission_columns,
  release_columns,
  max_zeros,
  comparison_label,
  group_label
) {

  if (length(admission_columns) != length(release_columns)) {
    stop("Admission and release sample counts do not match for ", comparison_label, ".")
  }

  admission_all <- numeric_matrix(data, admission_columns)
  release_all <- numeric_matrix(data, release_columns)

  # The paired columns are deliberately ordered by the same turtle IDs.
  # This check prevents accidental reordering of A/B samples.
  admission_ids <- sub("^Sample ([0-9]+)A$", "\\1", admission_columns)
  release_ids   <- sub("^Sample ([0-9]+)B$", "\\1", release_columns)

  if (!identical(admission_ids, release_ids)) {
    stop("Admission/release pairing is not aligned for ", comparison_label, ".")
  }

  total_zero_count <- rowSums(cbind(admission_all, release_all) == 0)

  # Filtering follows Table S1. Zero values remain in the analysis; the filter
  # removes only proteins exceeding the comparison-specific zero threshold.
  keep <- total_zero_count <= max_zeros

  data_f <- data[keep, , drop = FALSE]
  admission <- admission_all[keep, , drop = FALSE]
  release <- release_all[keep, , drop = FALSE]
  zero_count <- total_zero_count[keep]

  n_proteins <- nrow(data_f)

  p_value <- numeric(n_proteins)
  signed_rank_v <- numeric(n_proteins)

  for (i in seq_len(n_proteins)) {

    # Put release first so the test statistic follows the same direction as the
    # descriptive statistics and effect size. The two-sided p-value is unchanged
    # if the paired vectors are reversed.
    test <- exactRankTests::wilcox.exact(
      x = release[i, ],
      y = admission[i, ],
      paired = TRUE,
      alternative = "two.sided"
    )

    p_value[i] <- test$p.value
    signed_rank_v[i] <- unname(test$statistic)
  }

  padj_bh <- p.adjust(p_value, method = "BH")

  median_admission <- apply(admission, 1, median)
  iqr_admission <- apply(admission, 1, IQR)
  median_release <- apply(release, 1, median)
  iqr_release <- apply(release, 1, IQR)

  mean_admission <- rowMeans(admission)
  sd_admission <- apply(admission, 1, sd)
  mean_release <- rowMeans(release)
  sd_release <- apply(release, 1, sd)

  median_fc <- mapply(
    signed_fold_change,
    numerator = median_release,
    denominator = median_admission
  )

  mean_fc <- mapply(
    signed_fold_change,
    numerator = mean_release,
    denominator = mean_admission
  )

  paired_difference <- release - admission

  median_paired_difference <- apply(paired_difference, 1, median)
  iqr_paired_difference <- apply(paired_difference, 1, IQR)

  number_increased <- rowSums(paired_difference > 0)
  number_decreased <- rowSums(paired_difference < 0)
  number_unchanged <- rowSums(paired_difference == 0)
  nonzero_pairs <- rowSums(paired_difference != 0)

  rank_biserial <- apply(
    paired_difference,
    1,
    paired_rank_biserial
  )

  results <- data.frame(
    data_f[, metadata_columns, drop = FALSE],
    data_f[, c(admission_columns, release_columns), drop = FALSE],
    Zero_Count_Total = zero_count,
    p.value = p_value,
    padj_BH = padj_bh,
    Median_Admission = median_admission,
    IQR_Admission = iqr_admission,
    Median_Release = median_release,
    IQR_Release = iqr_release,
    Signed_Median_Fold_Change_Release_vs_Admission = median_fc,
    Median_Paired_Difference_Release_minus_Admission = median_paired_difference,
    IQR_Paired_Difference = iqr_paired_difference,
    Number_Increased = number_increased,
    Number_Decreased = number_decreased,
    Number_Unchanged = number_unchanged,
    Nonzero_Pairs = nonzero_pairs,
    Wilcoxon_Signed_Rank_V_Release_vs_Admission = signed_rank_v,
    Rank_Biserial_Correlation_Release_vs_Admission = rank_biserial,
    Mean_Admission = mean_admission,
    SD_Admission = sd_admission,
    Mean_Release = mean_release,
    SD_Release = sd_release,
    Signed_Mean_Fold_Change_Release_vs_Admission = mean_fc,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Match the Supporting Information presentation: lowest BH-adjusted p-values
  # first while preserving the original Table S3 protein order within ties.
  original_order <- seq_len(nrow(results))
  results <- results[order(results$padj_BH, original_order), , drop = FALSE]
  rownames(results) <- NULL

  message(
    comparison_label, ": ",
    nrow(results), " proteins retained; ",
    sum(results$padj_BH < 0.05, na.rm = TRUE),
    " proteins at BH FDR < 0.05."
  )

  results
}

# ------------------------------------------------------------------------------
# Table S6: paired OC2 admission vs release
# ------------------------------------------------------------------------------

S6 <- run_paired_wilcoxon(
  data = wsc,
  admission_columns = OC2_ADMISSION,
  release_columns = OC2_RELEASE,
  max_zeros = MAX_ZEROS_OC2,
  comparison_label = "Table S6 - OC2 admission vs release",
  group_label = "OC2"
)

# ------------------------------------------------------------------------------
# Table S7: paired OC4 admission vs release
# ------------------------------------------------------------------------------

S7 <- run_paired_wilcoxon(
  data = wsc,
  admission_columns = OC4_ADMISSION,
  release_columns = OC4_RELEASE,
  max_zeros = MAX_ZEROS_OC4,
  comparison_label = "Table S7 - OC4 admission vs release",
  group_label = "OC4"
)

# ------------------------------------------------------------------------------
# Reproducibility checks against the finalized manuscript results
# ------------------------------------------------------------------------------

stopifnot(
  nrow(S6) == 350L,
  sum(S6$padj_BH < 0.05, na.rm = TRUE) == 87L,
  nrow(S7) == 347L,
  sum(S7$padj_BH < 0.05, na.rm = TRUE) == 153L
)

# The smallest exact two-sided p-values are fixed by the paired sample sizes:
# 9 OC2 pairs -> 2 / 2^9  = 0.00390625
# 15 OC4 pairs -> 2 / 2^15 = 0.00006103515625
stopifnot(
  isTRUE(all.equal(min(S6$p.value), 0.00390625, tolerance = 1e-12)),
  isTRUE(all.equal(min(S7$p.value), 0.00006103515625, tolerance = 1e-12))
)

# ------------------------------------------------------------------------------
# Save machine-readable results
# ------------------------------------------------------------------------------

write.csv(
  S6,
  file.path(output_dir, "TableS6_Wilcoxon_OC2_Admission_vs_Release.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  S7,
  file.path(output_dir, "TableS7_Wilcoxon_OC4_Admission_vs_Release.csv"),
  row.names = FALSE,
  na = ""
)

summary_table <- data.frame(
  Comparison = c(
    "OC2 admission vs release",
    "OC4 admission vs release"
  ),
  Test = c(
    "Paired Wilcoxon signed-rank",
    "Paired Wilcoxon signed-rank"
  ),
  Number_of_Paired_Turtles = c(
    length(OC2_IDS),
    length(OC4_IDS)
  ),
  Total_Proteins = c(488L, 488L),
  Maximum_Zeros_Allowed = c(
    MAX_ZEROS_OC2,
    MAX_ZEROS_OC4
  ),
  Proteins_Retained = c(
    nrow(S6),
    nrow(S7)
  ),
  BH_FDR_lt_0.05 = c(
    sum(S6$padj_BH < 0.05, na.rm = TRUE),
    sum(S7$padj_BH < 0.05, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(output_dir, "Wilcoxon_paired_analysis_summary.csv"),
  row.names = FALSE
)

message("Paired Wilcoxon analysis completed successfully.")
message("Results written to: ", normalizePath(output_dir))
