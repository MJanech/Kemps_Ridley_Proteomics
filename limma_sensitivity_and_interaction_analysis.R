# ==============================================================================
# limma_sensitivity_and_interaction_analysis.R
#
# Purpose
#   Reproduce the empirical-Bayes limma analyses added during peer review of the
#   Kemp's ridley sea turtle plasma proteomics study.
#
#   This script performs:
#     Table S4: OC2 vs OC4 at admission
#     Table S5: OC2 vs OC4 after rehabilitation/release
#     Table S6: paired OC2 admission vs release
#     Table S7: paired OC4 admission vs release
#     Group x sampling-occasion interaction:
#       (OC4 release - OC4 admission) - (OC2 release - OC2 admission)
#
# Relationship to the primary analysis
#   Mann-Whitney U and paired Wilcoxon signed-rank tests remain the primary
#   protein-level analyses. limma is used here as a sensitivity analysis and
#   to directly evaluate whether the admission-to-release change differs
#   between OC2 and OC4.
#
# Input
#   Supplemental_Tables_1_through_8.xlsx
#
# Statistical framework
#   - wSC values are transformed as log2(wSC + 1).
#   - The +1 is a mathematical offset for the logarithm, not imputation.
#   - Numerical zero values exported from Scaffold are retained.
#   - No abundance values are imputed.
#   - Standard limma and limma-trend are both reported.
#   - limma-trend is the primary limma sensitivity analysis.
#   - Benjamini-Hochberg correction is applied separately to each comparison.
#   - FDR < 0.05 is considered significant.
#
# Protein universes
#   S4-S7 use the proteins retained in the corresponding finalized
#   protein-level comparison after the original zero-value filtering.
#
#   The interaction uses a conservative universe consisting of proteins retained
#   in all four primary comparisons (S4-S7). This gives 306 proteins.
#
# Interaction model
#   Repeated measurements are handled with turtle fixed effects. Two separate
#   within-turtle change coefficients are estimated:
#       OC2_change = OC2 release - OC2 admission
#       OC4_change = OC4 release - OC4 admission
#   The interaction contrast is:
#       OC4_change - OC2_change
#
# Expected results with the manuscript analysis environment
#   R 4.3.2
#   limma 3.58.1
#
#   Comparison                         NP     Standard limma   limma-trend
#   S4 OC2 vs OC4 admission            17          28             34
#   S5 OC2 vs OC4 release               0           0              0
#   S6 OC2 paired                       87         132            145
#   S7 OC4 paired                      153         158            162
#   Group x time interaction             0           0              0
#
#   Primary-test findings reproduced by limma-trend:
#       S4: 17 / 17
#       S6: 84 / 87
#       S7: 141 / 153
#
# Usage
#   Place this script and the Supplemental Tables workbook in the working
#   directory and run:
#       source("limma_sensitivity_and_interaction_analysis.R")
#
#   Or provide the workbook path from the command line:
#       Rscript limma_sensitivity_and_interaction_analysis.R path/to/workbook.xlsx
# ==============================================================================

required_packages <- c("readxl", "limma")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required. Install it before running this script."
    )
  }
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "Supplemental_Tables_1_through_8.xlsx"
output_dir <- "limma_sensitivity_results"

if (!file.exists(input_file)) {
  stop("Input workbook not found: ", input_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("R version: ", R.version.string, "\n", sep = "")
cat("limma version: ", as.character(utils::packageVersion("limma")), "\n", sep = "")

if (as.character(utils::packageVersion("limma")) != "3.58.1") {
  warning(
    "The manuscript analysis used limma 3.58.1. ",
    "Small numerical differences may occur with another version."
  )
}

# ------------------------------------------------------------------------------
# Sample definitions
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

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

check_columns <- function(data, required, label) {
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0) {
    stop(
      label, " is missing required columns:\n",
      paste(missing_columns, collapse = "\n")
    )
  }
}

to_numeric_matrix <- function(data, columns, label) {
  check_columns(data, columns, label)

  x <- as.matrix(data[, columns, drop = FALSE])
  storage.mode(x) <- "numeric"

  if (anyNA(x)) {
    stop(
      "Unexpected missing/non-numeric abundance values were found in ",
      label, "."
    )
  }

  x
}

agreement_class <- function(np_sig, trend_sig) {
  ifelse(
    np_sig & trend_sig, "Both",
    ifelse(
      np_sig & !trend_sig, "Primary test only",
      ifelse(!np_sig & trend_sig, "limma-trend only", "Neither")
    )
  )
}

# ------------------------------------------------------------------------------
# Unpaired limma sensitivity analysis: Tables S4 and S5
# ------------------------------------------------------------------------------

run_unpaired_limma <- function(
  sheet_name,
  oc2_columns,
  oc4_columns,
  comparison_label
) {

  dat <- readxl::read_excel(input_file, sheet = sheet_name, skip = 2)

  metadata_required <- c(
    "Row Number", "Identified.Proteins", "Gene.Name",
    "Accession.Number", "p.value", "padj_BH"
  )
  check_columns(dat, metadata_required, sheet_name)

  abundance <- to_numeric_matrix(
    dat,
    c(oc2_columns, oc4_columns),
    sheet_name
  )

  group <- factor(
    c(
      rep("OC2", length(oc2_columns)),
      rep("OC4", length(oc4_columns))
    ),
    levels = c("OC2", "OC4")
  )

  design <- stats::model.matrix(~ group)

  # Numerical zeros are retained. +1 is only the offset required before log2.
  y <- log2(abundance + 1)

  fit <- limma::lmFit(y, design)

  fit_standard <- limma::eBayes(
    fit,
    trend = FALSE,
    robust = FALSE
  )

  fit_trend <- limma::eBayes(
    fit,
    trend = TRUE,
    robust = FALSE
  )

  coefficient <- "groupOC4"

  standard <- limma::topTable(
    fit_standard,
    coef = coefficient,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  trend <- limma::topTable(
    fit_trend,
    coef = coefficient,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  np_sig <- dat$padj_BH < 0.05
  trend_sig <- trend$adj.P.Val < 0.05

  out <- data.frame(
    Row = dat[["Row Number"]],
    Protein = dat[["Identified.Proteins"]],
    Gene = dat[["Gene.Name"]],
    Accession = trimws(dat[["Accession.Number"]]),
    Primary_p = dat[["p.value"]],
    Primary_BH_FDR = dat[["padj_BH"]],
    limma_log2FC_OC4_minus_OC2 = trend$logFC,
    limma_AveExpr = trend$AveExpr,
    standard_t = standard$t,
    standard_p = standard$P.Value,
    standard_BH_FDR = standard$adj.P.Val,
    trend_t = trend$t,
    trend_p = trend$P.Value,
    trend_BH_FDR = trend$adj.P.Val,
    Primary_FDR_lt_0.05 = np_sig,
    Trend_FDR_lt_0.05 = trend_sig,
    Primary_reproduced_by_trend = np_sig & trend_sig,
    Agreement_class = agreement_class(np_sig, trend_sig),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  message(
    comparison_label, ": ",
    nrow(out), " proteins; primary = ",
    sum(np_sig, na.rm = TRUE),
    ", standard limma = ",
    sum(standard$adj.P.Val < 0.05, na.rm = TRUE),
    ", limma-trend = ",
    sum(trend_sig, na.rm = TRUE),
    "."
  )

  list(
    results = out,
    primary_n = sum(np_sig, na.rm = TRUE),
    standard_n = sum(standard$adj.P.Val < 0.05, na.rm = TRUE),
    trend_n = sum(trend_sig, na.rm = TRUE),
    primary_reproduced = sum(np_sig & trend_sig, na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# Paired limma sensitivity analysis: Tables S6 and S7
# ------------------------------------------------------------------------------

run_paired_limma <- function(
  sheet_name,
  admission_columns,
  release_columns,
  comparison_label
) {

  dat <- readxl::read_excel(input_file, sheet = sheet_name, skip = 3)

  metadata_required <- c(
    "Row Number", "Identified.Proteins", "Gene.Name",
    "Accession.Number", "p.value", "padj_BH",
    "Rank-biserial Correlation"
  )
  check_columns(dat, metadata_required, sheet_name)

  if (length(admission_columns) != length(release_columns)) {
    stop("Paired sample counts do not match for ", comparison_label, ".")
  }

  admission_ids <- sub("^Sample ([0-9]+)A$", "\\1", admission_columns)
  release_ids   <- sub("^Sample ([0-9]+)B$", "\\1", release_columns)

  if (!identical(admission_ids, release_ids)) {
    stop("Admission and release samples are not correctly paired for ", comparison_label, ".")
  }

  abundance <- to_numeric_matrix(
    dat,
    c(admission_columns, release_columns),
    sheet_name
  )

  n_pairs <- length(admission_columns)

  subject <- factor(rep(seq_len(n_pairs), 2))

  time <- factor(
    rep(c("Admission", "Release"), each = n_pairs),
    levels = c("Admission", "Release")
  )

  design <- stats::model.matrix(~ subject + time)

  y <- log2(abundance + 1)

  fit <- limma::lmFit(y, design)

  fit_standard <- limma::eBayes(
    fit,
    trend = FALSE,
    robust = FALSE
  )

  fit_trend <- limma::eBayes(
    fit,
    trend = TRUE,
    robust = FALSE
  )

  coefficient <- "timeRelease"

  standard <- limma::topTable(
    fit_standard,
    coef = coefficient,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  trend <- limma::topTable(
    fit_trend,
    coef = coefficient,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  np_sig <- dat$padj_BH < 0.05
  trend_sig <- trend$adj.P.Val < 0.05

  out <- data.frame(
    Row = dat[["Row Number"]],
    Protein = dat[["Identified.Proteins"]],
    Gene = dat[["Gene.Name"]],
    Accession = trimws(dat[["Accession.Number"]]),
    Primary_p = dat[["p.value"]],
    Primary_BH_FDR = dat[["padj_BH"]],
    Primary_rank_biserial = dat[["Rank-biserial Correlation"]],
    limma_log2FC_Release_minus_Admission = trend$logFC,
    limma_AveExpr = trend$AveExpr,
    standard_t = standard$t,
    standard_p = standard$P.Value,
    standard_BH_FDR = standard$adj.P.Val,
    trend_t = trend$t,
    trend_p = trend$P.Value,
    trend_BH_FDR = trend$adj.P.Val,
    Primary_FDR_lt_0.05 = np_sig,
    Trend_FDR_lt_0.05 = trend_sig,
    Primary_reproduced_by_trend = np_sig & trend_sig,
    Agreement_class = agreement_class(np_sig, trend_sig),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  message(
    comparison_label, ": ",
    nrow(out), " proteins; primary = ",
    sum(np_sig, na.rm = TRUE),
    ", standard limma = ",
    sum(standard$adj.P.Val < 0.05, na.rm = TRUE),
    ", limma-trend = ",
    sum(trend_sig, na.rm = TRUE),
    "."
  )

  list(
    results = out,
    primary_n = sum(np_sig, na.rm = TRUE),
    standard_n = sum(standard$adj.P.Val < 0.05, na.rm = TRUE),
    trend_n = sum(trend_sig, na.rm = TRUE),
    primary_reproduced = sum(np_sig & trend_sig, na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# Run S4-S7 sensitivity analyses
# ------------------------------------------------------------------------------

S4 <- run_unpaired_limma(
  sheet_name = "Table S4",
  oc2_columns = OC2_ADMISSION,
  oc4_columns = OC4_ADMISSION,
  comparison_label = "S4 OC2 vs OC4 at admission"
)

S5 <- run_unpaired_limma(
  sheet_name = "Table S5",
  oc2_columns = OC2_RELEASE,
  oc4_columns = OC4_RELEASE,
  comparison_label = "S5 OC2 vs OC4 after rehabilitation"
)

S6 <- run_paired_limma(
  sheet_name = "Table S6",
  admission_columns = OC2_ADMISSION,
  release_columns = OC2_RELEASE,
  comparison_label = "S6 OC2 paired admission vs release"
)

S7 <- run_paired_limma(
  sheet_name = "Table S7",
  admission_columns = OC4_ADMISSION,
  release_columns = OC4_RELEASE,
  comparison_label = "S7 OC4 paired admission vs release"
)

# ------------------------------------------------------------------------------
# Group x sampling-occasion interaction
# ------------------------------------------------------------------------------

s3 <- readxl::read_excel(input_file, sheet = "Table S3", skip = 3)
s4_raw <- readxl::read_excel(input_file, sheet = "Table S4", skip = 2)
s5_raw <- readxl::read_excel(input_file, sheet = "Table S5", skip = 2)
s6_raw <- readxl::read_excel(input_file, sheet = "Table S6", skip = 3)
s7_raw <- readxl::read_excel(input_file, sheet = "Table S7", skip = 3)

check_columns(
  s3,
  c("Identified Proteins", "Gene Name", "Accession Number"),
  "Table S3"
)

for (obj_name in c("s4_raw", "s5_raw", "s6_raw", "s7_raw")) {
  obj <- get(obj_name)
  check_columns(obj, "Accession.Number", obj_name)
}

# Conservative interaction universe: retained in all four primary comparisons.
common_accessions <- Reduce(
  intersect,
  list(
    trimws(s4_raw[["Accession.Number"]]),
    trimws(s5_raw[["Accession.Number"]]),
    trimws(s6_raw[["Accession.Number"]]),
    trimws(s7_raw[["Accession.Number"]])
  )
)

s3_accession <- trimws(s3[["Accession Number"]])
interaction_data <- s3[s3_accession %in% common_accessions, , drop = FALSE]

interaction_samples <- c(
  OC2_ADMISSION, OC2_RELEASE,
  OC4_ADMISSION, OC4_RELEASE
)

interaction_y <- to_numeric_matrix(
  interaction_data,
  interaction_samples,
  "Table S3 interaction data"
)

interaction_y <- log2(interaction_y + 1)

# Sample metadata are in the same order as interaction_samples.
subject <- factor(c(
  paste0("OC2_", OC2_IDS),
  paste0("OC2_", OC2_IDS),
  paste0("OC4_", OC4_IDS),
  paste0("OC4_", OC4_IDS)
))

group <- factor(
  c(
    rep("OC2", 2 * length(OC2_IDS)),
    rep("OC4", 2 * length(OC4_IDS))
  ),
  levels = c("OC2", "OC4")
)

time <- factor(
  c(
    rep("Admission", length(OC2_IDS)),
    rep("Release", length(OC2_IDS)),
    rep("Admission", length(OC4_IDS)),
    rep("Release", length(OC4_IDS))
  ),
  levels = c("Admission", "Release")
)

# One fixed-effect column per turtle absorbs the turtle-specific baseline.
subject_design <- stats::model.matrix(~ 0 + subject)
colnames(subject_design) <- make.names(colnames(subject_design))

OC2_change <- as.numeric(group == "OC2" & time == "Release")
OC4_change <- as.numeric(group == "OC4" & time == "Release")

interaction_design <- cbind(
  subject_design,
  OC2_change = OC2_change,
  OC4_change = OC4_change
)

if (qr(interaction_design)$rank != ncol(interaction_design)) {
  stop("The interaction design is rank deficient.")
}

interaction_fit <- limma::lmFit(
  interaction_y,
  interaction_design
)

interaction_contrast <- limma::makeContrasts(
  Interaction = OC4_change - OC2_change,
  levels = interaction_design
)

interaction_fit2 <- limma::contrasts.fit(
  interaction_fit,
  interaction_contrast
)

interaction_standard_fit <- limma::eBayes(
  interaction_fit2,
  trend = FALSE,
  robust = FALSE
)

interaction_trend_fit <- limma::eBayes(
  interaction_fit2,
  trend = TRUE,
  robust = FALSE
)

interaction_standard <- limma::topTable(
  interaction_standard_fit,
  coef = "Interaction",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

interaction_trend <- limma::topTable(
  interaction_trend_fit,
  coef = "Interaction",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

interaction_out <- data.frame(
  Protein = interaction_data[["Identified Proteins"]],
  Gene = interaction_data[["Gene Name"]],
  Accession = trimws(interaction_data[["Accession Number"]]),
  Interaction_log2_effect =
    interaction_trend$logFC,
  standard_t =
    interaction_standard$t,
  standard_p =
    interaction_standard$P.Value,
  standard_BH_FDR =
    interaction_standard$adj.P.Val,
  trend_t =
    interaction_trend$t,
  trend_p =
    interaction_trend$P.Value,
  trend_BH_FDR =
    interaction_trend$adj.P.Val,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

interaction_standard_n <- sum(
  interaction_out$standard_BH_FDR < 0.05,
  na.rm = TRUE
)

interaction_trend_n <- sum(
  interaction_out$trend_BH_FDR < 0.05,
  na.rm = TRUE
)

message(
  "Group x time interaction: ",
  nrow(interaction_out), " proteins; standard limma = ",
  interaction_standard_n,
  ", limma-trend = ",
  interaction_trend_n,
  "."
)

# ------------------------------------------------------------------------------
# Reproducibility checks against the finalized reviewer-validation results
# ------------------------------------------------------------------------------

stopifnot(
  nrow(S4$results) == 424L,
  S4$primary_n == 17L,
  S4$standard_n == 28L,
  S4$trend_n == 34L,
  S4$primary_reproduced == 17L,

  nrow(S5$results) == 328L,
  S5$primary_n == 0L,
  S5$standard_n == 0L,
  S5$trend_n == 0L,

  nrow(S6$results) == 350L,
  S6$primary_n == 87L,
  S6$standard_n == 132L,
  S6$trend_n == 145L,
  S6$primary_reproduced == 84L,

  nrow(S7$results) == 347L,
  S7$primary_n == 153L,
  S7$standard_n == 158L,
  S7$trend_n == 162L,
  S7$primary_reproduced == 141L,

  nrow(interaction_out) == 306L,
  nrow(interaction_design) == 48L,
  ncol(interaction_design) == 26L,
  qr(interaction_design)$rank == 26L,
  interaction_fit$df.residual[1] == 22L,
  interaction_standard_n == 0L,
  interaction_trend_n == 0L
)

# ------------------------------------------------------------------------------
# Save machine-readable outputs
# ------------------------------------------------------------------------------

write.csv(
  S4$results,
  file.path(output_dir, "S4_limma_sensitivity.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  S5$results,
  file.path(output_dir, "S5_limma_sensitivity.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  S6$results,
  file.path(output_dir, "S6_limma_sensitivity.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  S7$results,
  file.path(output_dir, "S7_limma_sensitivity.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  interaction_out,
  file.path(output_dir, "group_time_interaction_limma.csv"),
  row.names = FALSE,
  na = ""
)

summary_table <- data.frame(
  Comparison = c(
    "S4 OC2 vs OC4 at admission",
    "S5 OC2 vs OC4 after rehabilitation",
    "S6 OC2 paired admission vs release",
    "S7 OC4 paired admission vs release",
    "Group x sampling-occasion interaction"
  ),
  Proteins_tested = c(
    nrow(S4$results),
    nrow(S5$results),
    nrow(S6$results),
    nrow(S7$results),
    nrow(interaction_out)
  ),
  Primary_test_FDR_lt_0.05 = c(
    S4$primary_n,
    S5$primary_n,
    S6$primary_n,
    S7$primary_n,
    NA_integer_
  ),
  Standard_limma_FDR_lt_0.05 = c(
    S4$standard_n,
    S5$standard_n,
    S6$standard_n,
    S7$standard_n,
    interaction_standard_n
  ),
  Limma_trend_FDR_lt_0.05 = c(
    S4$trend_n,
    S5$trend_n,
    S6$trend_n,
    S7$trend_n,
    interaction_trend_n
  ),
  Primary_hits_reproduced_by_limma_trend = c(
    S4$primary_reproduced,
    NA_integer_,
    S6$primary_reproduced,
    S7$primary_reproduced,
    NA_integer_
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(output_dir, "limma_sensitivity_summary.csv"),
  row.names = FALSE,
  na = ""
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message("limma sensitivity and interaction analyses completed successfully.")
message("Results written to: ", normalizePath(output_dir))
