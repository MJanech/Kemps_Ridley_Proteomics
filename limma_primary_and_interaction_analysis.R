# ==============================================================================
# limma_primary_and_interaction_analysis.R
#
# Purpose
#   Reproduce the primary empirical-Bayes limma-trend protein-level analyses
#   for the Kemp's ridley sea turtle plasma proteomics study and the
#   group x sampling-occasion interaction analysis.
#
#   This script performs:
#     Table S4: OC2 vs OC4 at admission
#     Table S5: OC2 vs OC4 after rehabilitation/release
#     Table S6: paired OC2 admission vs release
#     Table S7: paired OC4 admission vs release
#     Group x sampling-occasion interaction:
#       (OC4 release - OC4 admission) - (OC2 release - OC2 admission)
#
# Statistical hierarchy used in the revised manuscript
#   - limma-trend is the PRIMARY protein-level differential-abundance analysis.
#   - Mann-Whitney U and paired Wilcoxon signed-rank tests are retained as
#     NONPARAMETRIC SENSITIVITY ANALYSES.
#   - Standard limma is also reported for comparison.
#
# Input
#   Supplemental_Tables_1_through_8_Limma_Primary.xlsx
#
# Statistical framework
#   - wSC values are transformed as log2(wSC + 1).
#   - The +1 is a mathematical offset for the logarithm, not imputation.
#   - Numerical zero values exported from Scaffold are retained.
#   - No abundance values are imputed.
#   - limma-trend uses empirical-Bayes moderation with trend = TRUE.
#   - Benjamini-Hochberg correction is applied separately to each comparison.
#   - FDR < 0.05 is considered significant.
#
# Protein universes
#   S4-S7 use the proteins retained in the corresponding finalized
#   protein-level comparison after the comparison-specific detection filtering.
#   The same retained proteins are used for limma-trend and the nonparametric
#   sensitivity analyses.
#
#   The interaction uses a conservative universe consisting of proteins retained
#   in all four protein-level comparisons (S4-S7). This gives 306 proteins.
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
#   Comparison                         limma-trend   Standard limma   NP sensitivity
#   S4 OC2 vs OC4 admission                 34             28              17
#   S5 OC2 vs OC4 release                    0              0               0
#   S6 OC2 paired                           145            132              87
#   S7 OC4 paired                           162            158             153
#   Group x time interaction                  0              0              NA
#
#   Nonparametric sensitivity findings reproduced by primary limma-trend:
#       S4: 17 / 17
#       S6: 84 / 87
#       S7: 141 / 153
#
# Usage
#   Place this script and the revised Supplemental Tables workbook in the
#   working directory and run:
#       source("limma_primary_and_interaction_analysis.R")
#
#   Or provide the workbook path from the command line:
#       Rscript limma_primary_and_interaction_analysis.R path/to/workbook.xlsx
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
input_file <- if (length(args) >= 1) {
  args[1]
} else {
  "Supplemental_Tables_1_through_8_Limma_Primary.xlsx"
}

output_dir <- "limma_primary_results"

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

normalize_name <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

find_column <- function(data, candidates) {
  normalized <- normalize_name(names(data))
  targets <- normalize_name(candidates)
  hit <- which(normalized %in% targets)

  if (length(hit) == 0) {
    return(NA_character_)
  }

  names(data)[hit[1]]
}

require_column <- function(data, candidates, label) {
  hit <- find_column(data, candidates)

  if (is.na(hit)) {
    stop(
      label, " is missing a required column. Expected one of:\n",
      paste(candidates, collapse = "\n")
    )
  }

  hit
}

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

agreement_class <- function(sensitivity_sig, primary_sig) {
  ifelse(
    sensitivity_sig & primary_sig, "Both",
    ifelse(
      sensitivity_sig & !primary_sig, "Nonparametric sensitivity only",
      ifelse(!sensitivity_sig & primary_sig, "limma-trend primary only", "Neither")
    )
  )
}

# ------------------------------------------------------------------------------
# Unpaired primary limma-trend analyses: Tables S4 and S5
# ------------------------------------------------------------------------------

run_unpaired_limma <- function(
  sheet_name,
  oc2_columns,
  oc4_columns,
  comparison_label
) {

  dat <- readxl::read_excel(
    input_file,
    sheet = sheet_name,
    skip = 2,
    .name_repair = "unique"
  )

  row_col <- require_column(dat, c("Row Number", "Row"), sheet_name)
  protein_col <- require_column(
    dat,
    c("Identified.Proteins", "Identified Proteins", "Protein"),
    sheet_name
  )
  gene_col <- require_column(
    dat,
    c("Gene.Name", "Gene Name", "Gene Symbol", "Gene"),
    sheet_name
  )
  accession_col <- require_column(
    dat,
    c("Accession.Number", "Accession Number", "Accession"),
    sheet_name
  )
  np_p_col <- require_column(
    dat,
    c("Mann-Whitney p", "Mann–Whitney p", "p.value"),
    sheet_name
  )
  np_fdr_col <- require_column(
    dat,
    c("Mann-Whitney BH FDR", "Mann–Whitney BH FDR", "padj_BH"),
    sheet_name
  )
  rb_col <- require_column(
    dat,
    c(
      "Rank-biserial Correlation (positive = OC4 higher)",
      "Rank-biserial Correlation"
    ),
    sheet_name
  )

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

  sensitivity_sig <- dat[[np_fdr_col]] < 0.05
  primary_sig <- trend$adj.P.Val < 0.05

  out <- data.frame(
    Row = dat[[row_col]],
    Protein = dat[[protein_col]],
    Gene = dat[[gene_col]],
    Accession = trimws(dat[[accession_col]]),
    limma_log2FC_OC4_minus_OC2 = trend$logFC,
    limma_AveExpr = trend$AveExpr,
    limma_trend_t = trend$t,
    limma_trend_p = trend$P.Value,
    limma_trend_BH_FDR = trend$adj.P.Val,
    Primary_limma_trend_FDR_lt_0.05 = primary_sig,
    standard_limma_t = standard$t,
    standard_limma_p = standard$P.Value,
    standard_limma_BH_FDR = standard$adj.P.Val,
    Nonparametric_sensitivity_p = dat[[np_p_col]],
    Nonparametric_sensitivity_BH_FDR = dat[[np_fdr_col]],
    Nonparametric_rank_biserial = dat[[rb_col]],
    Sensitivity_FDR_lt_0.05 = sensitivity_sig,
    Sensitivity_reproduced_by_primary = sensitivity_sig & primary_sig,
    Agreement_class = agreement_class(sensitivity_sig, primary_sig),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  message(
    comparison_label, ": ",
    nrow(out), " proteins; primary limma-trend = ",
    sum(primary_sig, na.rm = TRUE),
    ", standard limma = ",
    sum(standard$adj.P.Val < 0.05, na.rm = TRUE),
    ", nonparametric sensitivity = ",
    sum(sensitivity_sig, na.rm = TRUE),
    "."
  )

  list(
    results = out,
    primary_n = sum(primary_sig, na.rm = TRUE),
    standard_n = sum(standard$adj.P.Val < 0.05, na.rm = TRUE),
    sensitivity_n = sum(sensitivity_sig, na.rm = TRUE),
    sensitivity_reproduced = sum(sensitivity_sig & primary_sig, na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# Paired primary limma-trend analyses: Tables S6 and S7
# ------------------------------------------------------------------------------

run_paired_limma <- function(
  sheet_name,
  admission_columns,
  release_columns,
  comparison_label
) {

  dat <- readxl::read_excel(
    input_file,
    sheet = sheet_name,
    skip = 3,
    .name_repair = "unique"
  )

  row_col <- require_column(dat, c("Row Number", "Row"), sheet_name)
  protein_col <- require_column(
    dat,
    c("Identified.Proteins", "Identified Proteins", "Protein"),
    sheet_name
  )
  gene_col <- require_column(
    dat,
    c("Gene.Name", "Gene Name", "Gene Symbol", "Gene"),
    sheet_name
  )
  accession_col <- require_column(
    dat,
    c("Accession.Number", "Accession Number", "Accession"),
    sheet_name
  )
  np_p_col <- require_column(
    dat,
    c("Wilcoxon p", "p.value"),
    sheet_name
  )
  np_fdr_col <- require_column(
    dat,
    c("Wilcoxon BH FDR", "padj_BH"),
    sheet_name
  )
  rb_col <- require_column(
    dat,
    c(
      "Rank-biserial Correlation (positive = Rehabilitation higher)",
      "Rank-biserial Correlation"
    ),
    sheet_name
  )

  if (length(admission_columns) != length(release_columns)) {
    stop("Paired sample counts do not match for ", comparison_label, ".")
  }

  admission_ids <- sub("^Sample ([0-9]+)A$", "\\1", admission_columns)
  release_ids   <- sub("^Sample ([0-9]+)B$", "\\1", release_columns)

  if (!identical(admission_ids, release_ids)) {
    stop(
      "Admission and release samples are not correctly paired for ",
      comparison_label, "."
    )
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

  sensitivity_sig <- dat[[np_fdr_col]] < 0.05
  primary_sig <- trend$adj.P.Val < 0.05

  out <- data.frame(
    Row = dat[[row_col]],
    Protein = dat[[protein_col]],
    Gene = dat[[gene_col]],
    Accession = trimws(dat[[accession_col]]),
    limma_log2FC_Release_minus_Admission = trend$logFC,
    limma_AveExpr = trend$AveExpr,
    limma_trend_t = trend$t,
    limma_trend_p = trend$P.Value,
    limma_trend_BH_FDR = trend$adj.P.Val,
    Primary_limma_trend_FDR_lt_0.05 = primary_sig,
    standard_limma_t = standard$t,
    standard_limma_p = standard$P.Value,
    standard_limma_BH_FDR = standard$adj.P.Val,
    Nonparametric_sensitivity_p = dat[[np_p_col]],
    Nonparametric_sensitivity_BH_FDR = dat[[np_fdr_col]],
    Nonparametric_rank_biserial = dat[[rb_col]],
    Sensitivity_FDR_lt_0.05 = sensitivity_sig,
    Sensitivity_reproduced_by_primary = sensitivity_sig & primary_sig,
    Agreement_class = agreement_class(sensitivity_sig, primary_sig),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  message(
    comparison_label, ": ",
    nrow(out), " proteins; primary limma-trend = ",
    sum(primary_sig, na.rm = TRUE),
    ", standard limma = ",
    sum(standard$adj.P.Val < 0.05, na.rm = TRUE),
    ", nonparametric sensitivity = ",
    sum(sensitivity_sig, na.rm = TRUE),
    "."
  )

  list(
    results = out,
    primary_n = sum(primary_sig, na.rm = TRUE),
    standard_n = sum(standard$adj.P.Val < 0.05, na.rm = TRUE),
    sensitivity_n = sum(sensitivity_sig, na.rm = TRUE),
    sensitivity_reproduced = sum(sensitivity_sig & primary_sig, na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# Run S4-S7 primary limma-trend analyses
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

s3 <- readxl::read_excel(
  input_file,
  sheet = "Table S3",
  skip = 3,
  .name_repair = "unique"
)
s4_raw <- readxl::read_excel(
  input_file,
  sheet = "Table S4",
  skip = 2,
  .name_repair = "unique"
)
s5_raw <- readxl::read_excel(
  input_file,
  sheet = "Table S5",
  skip = 2,
  .name_repair = "unique"
)
s6_raw <- readxl::read_excel(
  input_file,
  sheet = "Table S6",
  skip = 3,
  .name_repair = "unique"
)
s7_raw <- readxl::read_excel(
  input_file,
  sheet = "Table S7",
  skip = 3,
  .name_repair = "unique"
)

s3_protein_col <- require_column(
  s3,
  c("Identified Proteins", "Identified.Proteins", "Protein"),
  "Table S3"
)
s3_gene_col <- require_column(
  s3,
  c("Gene Name", "Gene.Name", "Gene Symbol", "Gene"),
  "Table S3"
)
s3_accession_col <- require_column(
  s3,
  c("Accession Number", "Accession.Number", "Accession"),
  "Table S3"
)

s4_accession_col <- require_column(
  s4_raw,
  c("Accession.Number", "Accession Number", "Accession"),
  "Table S4"
)
s5_accession_col <- require_column(
  s5_raw,
  c("Accession.Number", "Accession Number", "Accession"),
  "Table S5"
)
s6_accession_col <- require_column(
  s6_raw,
  c("Accession.Number", "Accession Number", "Accession"),
  "Table S6"
)
s7_accession_col <- require_column(
  s7_raw,
  c("Accession.Number", "Accession Number", "Accession"),
  "Table S7"
)

# Conservative interaction universe: retained in all four protein-level comparisons.
common_accessions <- Reduce(
  intersect,
  list(
    trimws(s4_raw[[s4_accession_col]]),
    trimws(s5_raw[[s5_accession_col]]),
    trimws(s6_raw[[s6_accession_col]]),
    trimws(s7_raw[[s7_accession_col]])
  )
)

s3_accession <- trimws(s3[[s3_accession_col]])
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
  Protein = interaction_data[[s3_protein_col]],
  Gene = interaction_data[[s3_gene_col]],
  Accession = trimws(interaction_data[[s3_accession_col]]),
  Interaction_log2_effect = interaction_trend$logFC,
  limma_trend_t = interaction_trend$t,
  limma_trend_p = interaction_trend$P.Value,
  limma_trend_BH_FDR = interaction_trend$adj.P.Val,
  standard_limma_t = interaction_standard$t,
  standard_limma_p = interaction_standard$P.Value,
  standard_limma_BH_FDR = interaction_standard$adj.P.Val,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

interaction_primary_n <- sum(
  interaction_out$limma_trend_BH_FDR < 0.05,
  na.rm = TRUE
)

interaction_standard_n <- sum(
  interaction_out$standard_limma_BH_FDR < 0.05,
  na.rm = TRUE
)

message(
  "Group x time interaction: ",
  nrow(interaction_out), " proteins; primary limma-trend = ",
  interaction_primary_n,
  ", standard limma = ",
  interaction_standard_n,
  "."
)

# ------------------------------------------------------------------------------
# Reproducibility checks against the finalized results
# ------------------------------------------------------------------------------

stopifnot(
  nrow(S4$results) == 424L,
  S4$primary_n == 34L,
  S4$standard_n == 28L,
  S4$sensitivity_n == 17L,
  S4$sensitivity_reproduced == 17L,

  nrow(S5$results) == 328L,
  S5$primary_n == 0L,
  S5$standard_n == 0L,
  S5$sensitivity_n == 0L,

  nrow(S6$results) == 350L,
  S6$primary_n == 145L,
  S6$standard_n == 132L,
  S6$sensitivity_n == 87L,
  S6$sensitivity_reproduced == 84L,

  nrow(S7$results) == 347L,
  S7$primary_n == 162L,
  S7$standard_n == 158L,
  S7$sensitivity_n == 153L,
  S7$sensitivity_reproduced == 141L,

  nrow(interaction_out) == 306L,
  nrow(interaction_design) == 48L,
  ncol(interaction_design) == 26L,
  qr(interaction_design)$rank == 26L,
  interaction_fit$df.residual[1] == 22L,
  interaction_primary_n == 0L,
  interaction_standard_n == 0L
)

# ------------------------------------------------------------------------------
# Save machine-readable outputs
# ------------------------------------------------------------------------------

write.csv(
  S4$results,
  file.path(output_dir, "S4_limma_primary.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  S5$results,
  file.path(output_dir, "S5_limma_primary.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  S6$results,
  file.path(output_dir, "S6_limma_primary.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  S7$results,
  file.path(output_dir, "S7_limma_primary.csv"),
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
  Primary_limma_trend_FDR_lt_0.05 = c(
    S4$primary_n,
    S5$primary_n,
    S6$primary_n,
    S7$primary_n,
    interaction_primary_n
  ),
  Standard_limma_FDR_lt_0.05 = c(
    S4$standard_n,
    S5$standard_n,
    S6$standard_n,
    S7$standard_n,
    interaction_standard_n
  ),
  Nonparametric_sensitivity_FDR_lt_0.05 = c(
    S4$sensitivity_n,
    S5$sensitivity_n,
    S6$sensitivity_n,
    S7$sensitivity_n,
    NA_integer_
  ),
  Sensitivity_hits_reproduced_by_primary = c(
    S4$sensitivity_reproduced,
    NA_integer_,
    S6$sensitivity_reproduced,
    S7$sensitivity_reproduced,
    NA_integer_
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(output_dir, "limma_primary_summary.csv"),
  row.names = FALSE,
  na = ""
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message("Primary limma-trend and interaction analyses completed successfully.")
message("Results written to: ", normalizePath(output_dir))
