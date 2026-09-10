# ==============================================================================
# clinical_mixed_models_S2_S10.R
#
# Purpose
#   Reproduce the clinical/morphometric mixed-effects analyses reported in
#   Supplemental Table S2 and Supplemental Table S10B.
#
# Model for each outcome
#   Response ~ Time * Oiling + (1 | Turtle)
#
# Final analysis specification
#   - Time: Admission vs Rehabilitation (within-turtle repeated measure)
#   - Oiling: OC2 vs OC4 (between-turtle factor)
#   - Turtle: random intercept
#   - Weight is analyzed on the original scale.
#   - All eight clinical chemistry variables are log2-transformed.
#   - Missing measurements are not imputed.
#   - Models are fit by restricted maximum likelihood (REML).
#   - Sum-to-zero contrasts are used for Time and Oiling.
#   - Main fixed-effect tests use Type III F tests with Satterthwaite
#     denominator degrees of freedom (lmerTest).
#   - Prespecified simple-effect contrasts use Satterthwaite denominator
#     degrees of freedom (emmeans).
#   - Benjamini-Hochberg correction is applied across all 27 main fixed-effect
#     tests (9 outcomes x 3 effects), and separately across all 36 simple-effect
#     contrasts (9 outcomes x 4 contrasts).
#
# Input
#   Supplemental_Table_S10_Clinical_Data_and_Mixed_Models.xlsx
#   worksheet: S10A Individual Data
#
# Run from a terminal / command prompt
#   Rscript clinical_mixed_models_S2_S10.R \
#     Supplemental_Table_S10_Clinical_Data_and_Mixed_Models.xlsx
#
# Or in RStudio
#   1. Put this script and the workbook in the same folder.
#   2. Set that folder as the working directory.
#   3. Source this script.
#
# Outputs (written to Clinical_mixed_model_results/ next to the input workbook)
#   S10B_main_effects.csv
#   S10B_simple_effects.csv
#   TableS2_BH_pvalues.csv
#   mixed_model_diagnostics.csv
#   package_versions.csv
#   sessionInfo.txt
# ==============================================================================

options(digits = 17)

# ------------------------------------------------------------------------------
# Package checks
# ------------------------------------------------------------------------------
required_packages <- c("readxl", "lme4", "lmerTest", "emmeans")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    "\nInstall them once with:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

# ------------------------------------------------------------------------------
# Input/output paths
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (!interactive() && length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  "Supplemental_Table_S10_Clinical_Data_and_Mixed_Models.xlsx"
}

if (!file.exists(input_file)) {
  stop(
    "Input workbook not found: ", input_file,
    "\nEither place the workbook in the current working directory or pass its path ",
    "as the first command-line argument."
  )
}

input_file <- normalizePath(input_file, mustWork = TRUE)
input_sheet <- "S10A Individual Data"
output_dir <- file.path(dirname(input_file), "Clinical_mixed_model_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Read and validate the individual-level data
# ------------------------------------------------------------------------------
# In S10A, rows 1-3 are title/note/blank rows and row 4 contains column names.
dat <- readxl::read_excel(input_file, sheet = input_sheet, skip = 3)
names(dat) <- trimws(names(dat))

required_columns <- c(
  "Sample", "Time Point", "OC", "Group", "Weight, kg", "Glc, mg/dL",
  "Ca, mg/dL", "P, mg/dL", "UA, mg/dL", "AST, U/L", "CPK, U/L",
  "ALP, U/L", "GGT, U/L", "Days Since Admission"
)

missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0) {
  stop(
    "The input worksheet is missing required column(s):\n",
    paste(missing_columns, collapse = "\n")
  )
}

# Guard against accidental duplicate turtle/time records.
duplicate_key <- duplicated(dat[c("Sample", "Time Point")])
if (any(duplicate_key)) {
  stop(
    "Duplicate Sample + Time Point rows were found in S10A. ",
    "The analysis was stopped to prevent silent duplication."
  )
}

# Confirm only expected factor labels occur.
expected_times <- c("Admission", "Rehabilitation")
expected_oiling <- c("OC2", "OC4")

unexpected_times <- setdiff(unique(na.omit(dat[["Time Point"]])), expected_times)
unexpected_oiling <- setdiff(unique(na.omit(dat[["OC"]])), expected_oiling)

if (length(unexpected_times) > 0) {
  stop("Unexpected Time Point value(s): ", paste(unexpected_times, collapse = ", "))
}
if (length(unexpected_oiling) > 0) {
  stop("Unexpected OC value(s): ", paste(unexpected_oiling, collapse = ", "))
}

dat$Turtle <- factor(as.character(dat$Sample))
dat$Time <- factor(dat[["Time Point"]], levels = expected_times)
dat$Oiling <- factor(dat$OC, levels = expected_oiling)

# Sum-to-zero contrasts are required for the Type III marginal main-effect tests.
contrasts(dat$Time) <- contr.sum(2)
contrasts(dat$Oiling) <- contr.sum(2)

analytes <- data.frame(
  Analyte = c(
    "Weight, kg", "Glc, mg/dL", "Ca, mg/dL", "P, mg/dL", "UA, mg/dL",
    "AST, U/L", "CPK, U/L", "ALP, U/L", "GGT, U/L"
  ),
  Analysis_scale = c("Raw", rep("log2", 8)),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || length(x) > 5000) return(NA_real_)
  unname(stats::shapiro.test(x)$p.value)
}

fit_one_outcome <- function(analyte, scale_label) {
  d <- dat[!is.na(dat[[analyte]]), c("Turtle", "Time", "Oiling", analyte)]
  names(d)[4] <- "Observed"

  d$Observed <- as.numeric(d$Observed)

  if (scale_label == "Raw") {
    d$Response <- d$Observed
  } else {
    if (any(d$Observed <= 0, na.rm = TRUE)) {
      stop(
        "Cannot log2-transform ", analyte,
        " because at least one non-missing value is <= 0."
      )
    }
    d$Response <- log2(d$Observed)
  }

  # Preserve the intended two-level factors and contrast coding after subsetting.
  d$Time <- factor(d$Time, levels = expected_times)
  d$Oiling <- factor(d$Oiling, levels = expected_oiling)
  contrasts(d$Time) <- contr.sum(2)
  contrasts(d$Oiling) <- contr.sum(2)

  fit <- lmerTest::lmer(
    Response ~ Time * Oiling + (1 | Turtle),
    data = d,
    REML = TRUE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000)
    )
  )

  # Type III marginal tests with Satterthwaite denominator degrees of freedom.
  a <- stats::anova(fit, type = 3, ddf = "Satterthwaite")
  effect_rows <- c("Time", "Oiling", "Time:Oiling")
  effect_labels <- c("Time", "Oiling", "Time x Oiling")

  if (!all(effect_rows %in% rownames(a))) {
    stop("Expected fixed-effect rows were not returned for ", analyte)
  }

  main <- data.frame(
    Analyte = analyte,
    Analysis_scale = scale_label,
    n_observations = nrow(d),
    n_turtles = length(unique(d$Turtle)),
    Effect = effect_labels,
    NumDF = as.numeric(a[effect_rows, "NumDF"]),
    DenDF = as.numeric(a[effect_rows, "DenDF"]),
    F = as.numeric(a[effect_rows, "F value"]),
    Raw_p = as.numeric(a[effect_rows, "Pr(>F)"]),
    stringsAsFactors = FALSE
  )

  # Prespecified simple effects.
  # revpairwise yields Rehabilitation - Admission and OC4 - OC2 because of the
  # factor level order defined above.
  emm_time <- emmeans::emmeans(
    fit,
    specs = ~ Time | Oiling,
    lmer.df = "satterthwaite"
  )
  time_ctr <- as.data.frame(
    summary(
      emmeans::contrast(emm_time, method = "revpairwise", adjust = "none"),
      infer = c(FALSE, TRUE),
      adjust = "none"
    )
  )

  required_time_cols <- c("Oiling", "estimate", "SE", "df", "t.ratio", "p.value")
  missing_time_cols <- setdiff(required_time_cols, names(time_ctr))
  if (length(missing_time_cols) > 0) {
    stop(
      "emmeans time-contrast output for ", analyte,
      " is missing expected column(s): ",
      paste(missing_time_cols, collapse = ", ")
    )
  }
  if (nrow(time_ctr) != 2) {
    stop("Expected 2 time-within-oiling contrasts for ", analyte,
         "; received ", nrow(time_ctr), ".")
  }

  time_rows <- data.frame(
    Analyte = analyte,
    Contrast = paste0(time_ctr$Oiling, ": Rehabilitation - Admission"),
    Analysis_scale = scale_label,
    Estimate = as.numeric(time_ctr$estimate),
    SE = as.numeric(time_ctr$SE),
    DenDF = as.numeric(time_ctr$df),
    t = as.numeric(time_ctr$t.ratio),
    Raw_p = as.numeric(time_ctr$p.value),
    stringsAsFactors = FALSE
  )

  emm_oiling <- emmeans::emmeans(
    fit,
    specs = ~ Oiling | Time,
    lmer.df = "satterthwaite"
  )
  oil_ctr <- as.data.frame(
    summary(
      emmeans::contrast(emm_oiling, method = "revpairwise", adjust = "none"),
      infer = c(FALSE, TRUE),
      adjust = "none"
    )
  )

  required_oil_cols <- c("Time", "estimate", "SE", "df", "t.ratio", "p.value")
  missing_oil_cols <- setdiff(required_oil_cols, names(oil_ctr))
  if (length(missing_oil_cols) > 0) {
    stop(
      "emmeans oiling-contrast output for ", analyte,
      " is missing expected column(s): ",
      paste(missing_oil_cols, collapse = ", ")
    )
  }
  if (nrow(oil_ctr) != 2) {
    stop("Expected 2 oiling-within-time contrasts for ", analyte,
         "; received ", nrow(oil_ctr), ".")
  }

  oil_rows <- data.frame(
    Analyte = analyte,
    Contrast = paste0(oil_ctr$Time, ": OC4 - OC2"),
    Analysis_scale = scale_label,
    Estimate = as.numeric(oil_ctr$estimate),
    SE = as.numeric(oil_ctr$SE),
    DenDF = as.numeric(oil_ctr$df),
    t = as.numeric(oil_ctr$t.ratio),
    Raw_p = as.numeric(oil_ctr$p.value),
    stringsAsFactors = FALSE
  )

  vc <- as.data.frame(lme4::VarCorr(fit))
  random_sd <- vc$sdcor[vc$grp == "Turtle"][1]

  diagnostics <- data.frame(
    Analyte = analyte,
    Analysis_scale = scale_label,
    n_observations = nrow(d),
    n_turtles = length(unique(d$Turtle)),
    Singular_fit = lme4::isSingular(fit, tol = 1e-5),
    Random_intercept_SD = as.numeric(random_sd),
    Residual_SD = as.numeric(stats::sigma(fit)),
    Residual_Shapiro_p = safe_shapiro(stats::residuals(fit)),
    stringsAsFactors = FALSE
  )

  list(
    main = main,
    simple = rbind(time_rows, oil_rows),
    diagnostics = diagnostics
  )
}

# ------------------------------------------------------------------------------
# Fit all nine models
# ------------------------------------------------------------------------------
main_list <- vector("list", nrow(analytes))
simple_list <- vector("list", nrow(analytes))
diag_list <- vector("list", nrow(analytes))

for (i in seq_len(nrow(analytes))) {
  fit_i <- fit_one_outcome(
    analyte = analytes$Analyte[i],
    scale_label = analytes$Analysis_scale[i]
  )
  main_list[[i]] <- fit_i$main
  simple_list[[i]] <- fit_i$simple
  diag_list[[i]] <- fit_i$diagnostics
}

main <- do.call(rbind, main_list)
simple <- do.call(rbind, simple_list)
diagnostics <- do.call(rbind, diag_list)

# ------------------------------------------------------------------------------
# Multiple-testing correction
# ------------------------------------------------------------------------------
# One BH family for the 27 primary fixed-effect tests.
main$BH_p <- stats::p.adjust(main$Raw_p, method = "BH")

# A separate BH family for the 36 prespecified simple-effect contrasts.
simple$BH_p <- stats::p.adjust(simple$Raw_p, method = "BH")

# Stable publication order.
main$Analyte <- factor(main$Analyte, levels = analytes$Analyte)
main$Effect <- factor(main$Effect, levels = c("Time", "Oiling", "Time x Oiling"))
main <- main[order(main$Analyte, main$Effect), ]
main$Analyte <- as.character(main$Analyte)
main$Effect <- as.character(main$Effect)
rownames(main) <- NULL

contrast_order <- c(
  "OC2: Rehabilitation - Admission",
  "OC4: Rehabilitation - Admission",
  "Admission: OC4 - OC2",
  "Rehabilitation: OC4 - OC2"
)
simple$Analyte <- factor(simple$Analyte, levels = analytes$Analyte)
simple$Contrast <- factor(simple$Contrast, levels = contrast_order)
simple <- simple[order(simple$Analyte, simple$Contrast), ]
simple$Analyte <- as.character(simple$Analyte)
simple$Contrast <- as.character(simple$Contrast)
rownames(simple) <- NULL

# Table S2 summary: BH-adjusted p-values from the same 27-test family.
summary_s2 <- data.frame(
  Analyte = analytes$Analyte,
  Analysis_scale = analytes$Analysis_scale,
  Time_BH_p = NA_real_,
  Oiling_BH_p = NA_real_,
  Time_x_Oiling_BH_p = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(summary_s2))) {
  z <- main[main$Analyte == summary_s2$Analyte[i], ]
  summary_s2$Time_BH_p[i] <- z$BH_p[z$Effect == "Time"]
  summary_s2$Oiling_BH_p[i] <- z$BH_p[z$Effect == "Oiling"]
  summary_s2$Time_x_Oiling_BH_p[i] <- z$BH_p[z$Effect == "Time x Oiling"]
}

# ------------------------------------------------------------------------------
# Write reproducible outputs
# ------------------------------------------------------------------------------
utils::write.csv(
  main,
  file.path(output_dir, "S10B_main_effects.csv"),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  simple,
  file.path(output_dir, "S10B_simple_effects.csv"),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  summary_s2,
  file.path(output_dir, "TableS2_BH_pvalues.csv"),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  diagnostics,
  file.path(output_dir, "mixed_model_diagnostics.csv"),
  row.names = FALSE,
  na = ""
)

package_versions <- data.frame(
  Component = c("R", required_packages),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    vapply(required_packages, function(x) as.character(utils::packageVersion(x)), character(1))
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  package_versions,
  file.path(output_dir, "package_versions.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

cat("\nClinical mixed-effects analysis completed successfully.\n")
cat("Input: ", input_file, "\n", sep = "")
cat("Output directory: ", normalizePath(output_dir), "\n\n", sep = "")
cat("Table S2 BH-adjusted p-values:\n")
print(summary_s2, row.names = FALSE)
cat("\nAny singular fits are reported in mixed_model_diagnostics.csv.\n")
