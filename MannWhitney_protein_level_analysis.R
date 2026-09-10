# ==============================================================================
# MannWhitney_protein_level_analysis.R
#
# Purpose
#   Reproduce the protein-level Mann-Whitney U analyses reported in:
#     Table S4: OC2 vs OC4 at admission
#     Table S5: OC2 vs OC4 after rehabilitation
#
# Input
#   Supplemental_Tables_1_through_8.xlsx
#   Sheet "Table S3" contains the study-wide weighted spectral counts (wSC).
#
# Statistical framework
#   - Two-sided Mann-Whitney U / Wilcoxon rank-sum test
#   - Zero values are retained as non-detections; no abundance values are imputed.
#   - Proteins with >21 zero values across the 24 samples in a comparison are
#     excluded, matching Table S1 (maximum permitted zeros = 21/24).
#   - Benjamini-Hochberg correction is applied separately within each comparison.
#   - R's wilcox.test() default exact behavior is retained: exact p-values are
#     used when possible; tied values are handled with the normal approximation
#     and continuity correction.
#   - The reported U statistic is oriented as OC4 vs OC2 so that the
#     rank-biserial correlation is positive when OC4 is higher.
#
# Expected validation
#   Table S4: 424 proteins retained; 17 proteins at BH FDR < 0.05
#   Table S5: 328 proteins retained;  0 proteins at BH FDR < 0.05
#
# Usage
#   Place this script and the Supplemental Tables workbook in the working
#   directory and run:
#       source("MannWhitney_protein_level_analysis.R")
#
#   An alternative workbook path can be supplied from the command line:
#       Rscript MannWhitney_protein_level_analysis.R path/to/workbook.xlsx
# ==============================================================================

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. Install it with install.packages('readxl').")
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "Supplemental_Tables_1_through_8.xlsx"
output_dir <- "MannWhitney_results"

if (!file.exists(input_file)) {
  stop("Input workbook not found: ", input_file)
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Sample definitions used in the manuscript
# ------------------------------------------------------------------------------

OC2_ADMISSION <- paste(
  "Sample",
  c("114A", "115A", "116A", "133A", "137A", "162A", "192A", "249A", "278A")
)

OC2_RELEASE <- paste(
  "Sample",
  c("114B", "115B", "116B", "133B", "137B", "162B", "192B", "249B", "278B")
)

OC4_ADMISSION <- paste(
  "Sample",
  c("081A", "087A", "088A", "098A", "100A", "107A", "108A", "110A",
    "125A", "193A", "203A", "204A", "254A", "265A", "269A")
)

OC4_RELEASE <- paste(
  "Sample",
  c("081B", "087B", "088B", "098B", "100B", "107B", "108B", "110B",
    "125B", "193B", "203B", "204B", "254B", "265B", "269B")
)

MAX_ZEROS <- 21L  # Table S1: 21/24 zeros permitted for S4 and S5

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

run_mann_whitney <- function(data, oc2_columns, oc4_columns,
                             comparison_label, oc2_label, oc4_label) {

  oc2_all <- numeric_matrix(data, oc2_columns)
  oc4_all <- numeric_matrix(data, oc4_columns)

  total_zero_count <- rowSums(cbind(oc2_all, oc4_all) == 0)

  # Filtering follows Table S1. Zero values are retained in the statistical
  # test; the filter only removes proteins with too many zero values to provide
  # sufficient information for the comparison.
  keep <- total_zero_count <= MAX_ZEROS

  data_f <- data[keep, , drop = FALSE]
  oc2 <- oc2_all[keep, , drop = FALSE]
  oc4 <- oc4_all[keep, , drop = FALSE]
  zero_count <- total_zero_count[keep]

  n_proteins <- nrow(data_f)

  p_value <- numeric(n_proteins)
  u_statistic <- numeric(n_proteins)

  for (i in seq_len(n_proteins)) {
    # Put OC4 first so the reported statistic is U(OC4 vs OC2).
    # The two-sided p-value is unchanged by group order.
    test <- suppressWarnings(
      wilcox.test(
        x = oc4[i, ],
        y = oc2[i, ],
        paired = FALSE,
        alternative = "two.sided",
        exact = NULL,
        correct = TRUE
      )
    )

    p_value[i] <- test$p.value
    u_statistic[i] <- unname(test$statistic)
  }

  padj_bh <- p.adjust(p_value, method = "BH")

  median_oc2 <- apply(oc2, 1, median)
  iqr_oc2 <- apply(oc2, 1, IQR)
  median_oc4 <- apply(oc4, 1, median)
  iqr_oc4 <- apply(oc4, 1, IQR)

  mean_oc2 <- rowMeans(oc2)
  sd_oc2 <- apply(oc2, 1, sd)
  mean_oc4 <- rowMeans(oc4)
  sd_oc4 <- apply(oc4, 1, sd)

  median_fc <- mapply(signed_fold_change, median_oc4, median_oc2)
  mean_fc <- mapply(signed_fold_change, mean_oc4, mean_oc2)

  # Positive values mean greater abundance/ranks in OC4.
  rank_biserial <- (2 * u_statistic / (ncol(oc4) * ncol(oc2))) - 1

  results <- data.frame(
    data_f[, metadata_columns, drop = FALSE],
    data_f[, c(oc2_columns, oc4_columns), drop = FALSE],
    Zero_Count_Total = zero_count,
    p.value = p_value,
    padj_BH = padj_bh,
    Median_OC2 = median_oc2,
    IQR_OC2 = iqr_oc2,
    Median_OC4 = median_oc4,
    IQR_OC4 = iqr_oc4,
    Signed_Median_Fold_Change_OC4_vs_OC2 = median_fc,
    Mean_OC2 = mean_oc2,
    SD_OC2 = sd_oc2,
    Mean_OC4 = mean_oc4,
    SD_OC4 = sd_oc4,
    Signed_Mean_Fold_Change_OC4_vs_OC2 = mean_fc,
    Mann_Whitney_U_OC4_vs_OC2 = u_statistic,
    Rank_Biserial_Correlation_OC4_vs_OC2 = rank_biserial,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Match the SI presentation: lowest BH-adjusted p-values first while
  # preserving the original protein order when adjusted p-values are tied.
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
# Table S4: OC2 vs OC4 at admission
# ------------------------------------------------------------------------------

S4 <- run_mann_whitney(
  data = wsc,
  oc2_columns = OC2_ADMISSION,
  oc4_columns = OC4_ADMISSION,
  comparison_label = "Table S4 - OC2 vs OC4 at admission",
  oc2_label = "OC2_0",
  oc4_label = "OC4_0"
)

# ------------------------------------------------------------------------------
# Table S5: OC2 vs OC4 after rehabilitation
# ------------------------------------------------------------------------------

S5 <- run_mann_whitney(
  data = wsc,
  oc2_columns = OC2_RELEASE,
  oc4_columns = OC4_RELEASE,
  comparison_label = "Table S5 - OC2 vs OC4 after rehabilitation",
  oc2_label = "OC2_1",
  oc4_label = "OC4_1"
)

# ------------------------------------------------------------------------------
# Reproducibility checks against the finalized manuscript results
# ------------------------------------------------------------------------------

stopifnot(
  nrow(S4) == 424L,
  sum(S4$padj_BH < 0.05, na.rm = TRUE) == 17L,
  nrow(S5) == 328L,
  sum(S5$padj_BH < 0.05, na.rm = TRUE) == 0L
)

# ------------------------------------------------------------------------------
# Save machine-readable results
# ------------------------------------------------------------------------------

write.csv(
  S4,
  file.path(output_dir, "TableS4_MannWhitney_OC2_vs_OC4_Admission.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  S5,
  file.path(output_dir, "TableS5_MannWhitney_OC2_vs_OC4_Release.csv"),
  row.names = FALSE,
  na = ""
)

summary_table <- data.frame(
  Comparison = c("OC2 vs OC4 at admission", "OC2 vs OC4 after rehabilitation"),
  Test = c("Mann-Whitney U", "Mann-Whitney U"),
  Total_Proteins = c(488L, 488L),
  Maximum_Zeros_Allowed = c(MAX_ZEROS, MAX_ZEROS),
  Proteins_Retained = c(nrow(S4), nrow(S5)),
  BH_FDR_lt_0.05 = c(
    sum(S4$padj_BH < 0.05, na.rm = TRUE),
    sum(S5$padj_BH < 0.05, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(output_dir, "MannWhitney_analysis_summary.csv"),
  row.names = FALSE
)

message("Mann-Whitney analysis completed successfully.")
message("Results written to: ", normalizePath(output_dir))
