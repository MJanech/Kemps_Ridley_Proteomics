# PCA / PERMANOVA / PERMDISP figures for publication
# Kemp's ridley plasma proteomics
#
# INPUT:
#   Supplemental_Tables_1_through_8.xlsx
#   Sheet: Table S3
#
# OUTPUT:
#   Results are written in the same folder as the input workbook.
#
# NOTES:
# - Numerical zeros are retained as Scaffold-exported non-detections.
# - No missing values are imputed.
# - PCA is centered and unit-variance scaled.
# - PERMANOVA/PERMDISP use Bray-Curtis dissimilarity.
# - Admission-to-release permutations are blocked within turtle.
# - The four-group PERMANOVA is retained as a descriptive/global overview.

required <- c("readxl", "ggplot2", "vegan", "permute", "ggrepel", "patchwork")

missing_pkgs <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Please install required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nFor example:\ninstall.packages(c(",
    paste0('"', missing_pkgs, '"', collapse = ", "),
    "))"
  )
}

library(readxl)
library(ggplot2)
library(vegan)
library(permute)

set.seed(20260824)

# -------------------------------------------------------------------------
# 1. Locate the workbook
# -------------------------------------------------------------------------
input_name <- "Supplemental_Tables_1_through_8.xlsx"
args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path(getwd(), input_name)
}

if (!file.exists(input_file)) {
  stop(
    "Input workbook not found: ", input_file, "\n",
    "Place ", input_name, " in the working directory or provide its path as the first command-line argument."
  )
}

if (!grepl("\\.xls[x]?$", input_file, ignore.case = TRUE)) {
  stop("The selected input does not look like an Excel workbook: ", input_file)
}

cat("Using input workbook:\n", input_file, "\n\n")

out_dir <- dirname(normalizePath(input_file, winslash = "/", mustWork = TRUE))
cat("Results will be written to:\n", out_dir, "\n\n")

# -------------------------------------------------------------------------
# 2. Read Table S3
# -------------------------------------------------------------------------
dat <- read_excel(input_file, sheet = "Table S3", skip = 3)

sample_cols <- 9:56

if (ncol(dat) < max(sample_cols)) {
  stop("Table S3 does not contain the expected 48 sample columns.")
}

abund <- as.matrix(dat[, sample_cols])
storage.mode(abund) <- "numeric"

if (anyNA(abund)) {
  stop(
    "Unexpected NA values were found in the 48-sample wSC matrix.\n",
    "Review Table S3 before PCA."
  )
}

sample_names <- names(dat)[sample_cols]

protein_label <- if (ncol(dat) >= 2) as.character(dat[[2]]) else paste0("Protein_", seq_len(nrow(dat)))
gene_label    <- if (ncol(dat) >= 3) as.character(dat[[3]]) else NA_character_
accession     <- if ("Accession Number" %in% names(dat)) as.character(dat[["Accession Number"]]) else NA_character_

X <- t(abund)

# -------------------------------------------------------------------------
# 3. Sample metadata
# -------------------------------------------------------------------------
group4 <- factor(
  c(
    rep("OC2 Admission", 9),
    rep("OC2 Release", 9),
    rep("OC4 Admission", 15),
    rep("OC4 Release", 15)
  ),
  levels = c("OC2 Admission", "OC2 Release", "OC4 Admission", "OC4 Release")
)

oiling <- factor(
  c(rep("OC2", 18), rep("OC4", 30)),
  levels = c("OC2", "OC4")
)

occasion <- factor(
  c(
    rep("Admission", 9),
    rep("Release", 9),
    rep("Admission", 15),
    rep("Release", 15)
  ),
  levels = c("Admission", "Release")
)

sample_id <- sub("^Sample\\s+", "", sample_names)
turtle_id <- sub("[AB]$", "", sample_id)

if (length(unique(turtle_id)) != 24) {
  warning(
    "Turtle IDs could not be recovered cleanly from sample names. ",
    "Using known paired ordering instead."
  )
  turtle_id <- c(
    paste0("OC2_", seq_len(9)),
    paste0("OC2_", seq_len(9)),
    paste0("OC4_", seq_len(15)),
    paste0("OC4_", seq_len(15))
  )
}

meta <- data.frame(
  Sample = sample_names,
  Sample_ID = sample_id,
  Turtle = factor(turtle_id),
  Oiling = oiling,
  Occasion = occasion,
  Group = group4,
  stringsAsFactors = TRUE
)

# -------------------------------------------------------------------------
# 4. PCA helper
# -------------------------------------------------------------------------
run_pca <- function(M, label, protein_index) {
  sds <- apply(M, 2, sd)
  keep_var <- is.finite(sds) & sds > 0

  M2 <- M[, keep_var, drop = FALSE]
  kept_index <- protein_index[keep_var]

  pc <- prcomp(M2, center = TRUE, scale. = TRUE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)

  scores <- data.frame(
    Analysis = label,
    Sample = sample_names,
    Sample_ID = sample_id,
    Turtle = turtle_id,
    Group = group4,
    Oiling = oiling,
    Occasion = occasion,
    PC1 = pc$x[, 1],
    PC2 = pc$x[, 2],
    stringsAsFactors = FALSE
  )

  loadings <- data.frame(
    Analysis = label,
    Protein_Row = kept_index,
    Protein = protein_label[kept_index],
    Gene = gene_label[kept_index],
    Accession = accession[kept_index],
    PC1_Loading = pc$rotation[, 1],
    PC2_Loading = pc$rotation[, 2],
    Absolute_PC1_Loading = abs(pc$rotation[, 1]),
    Absolute_PC2_Loading = abs(pc$rotation[, 2]),
    stringsAsFactors = FALSE
  )

  loadings$Combined_PC1_PC2_Strength <-
    sqrt(loadings$PC1_Loading^2 + loadings$PC2_Loading^2)

  summary <- data.frame(
    Analysis = label,
    Proteins_Analyzed = ncol(M2),
    PC1_Percent = ve[1],
    PC2_Percent = ve[2],
    PC1_PC2_Combined_Percent = ve[1] + ve[2],
    stringsAsFactors = FALSE
  )

  list(
    pc = pc,
    scores = scores,
    loadings = loadings,
    summary = summary,
    variance = ve,
    matrix = M2,
    protein_index = kept_index
  )
}

all_protein_index <- seq_len(nrow(abund))

pca_original <- run_pca(X, "Original: raw wSC, all proteins", all_protein_index)
X_log <- log2(X + 1)
pca_log_all <- run_pca(X_log, "Sensitivity: log2(wSC+1), all proteins", all_protein_index)

group_indices <- split(seq_len(nrow(X)), group4)

detected50_every_group <- vapply(
  seq_len(ncol(X)),
  function(j) {
    all(vapply(
      group_indices,
      function(ii) mean(X[ii, j] > 0) >= 0.50,
      logical(1)
    ))
  },
  logical(1)
)

X_log_50 <- X_log[, detected50_every_group, drop = FALSE]
protein_index_50 <- all_protein_index[detected50_every_group]
pca_log_50 <- run_pca(X_log_50, "Sensitivity: log2(wSC+1), >=50% detected in every group", protein_index_50)

pca_results <- list(
  Raw_wSC = pca_original,
  Log2_all = pca_log_all,
  Log2_detected50 = pca_log_50
)

# -------------------------------------------------------------------------
# 5. PERMANOVA / PERMDISP helper
# -------------------------------------------------------------------------
make_block_permutations <- function(blocks, nperm = 9999) {
  ctrl <- how(nperm = nperm)
  setBlocks(ctrl) <- blocks
  ctrl
}

run_multivariate_test <- function(M, analysis_name, comparison_name, rows, group, blocks = NULL, nperm = 9999) {
  mat <- M[rows, , drop = FALSE]
  grp <- droplevels(factor(group[rows]))

  if (any(mat < 0, na.rm = TRUE)) {
    stop("Negative values detected before Bray-Curtis calculation.")
  }

  d <- vegdist(mat, method = "bray")

  if (is.null(blocks)) {
    perms <- nperm
    perm_label <- "Unrestricted"
  } else {
    block_factor <- droplevels(factor(blocks[rows]))
    perms <- make_block_permutations(block_factor, nperm = nperm)
    perm_label <- "Blocked within turtle"
  }

  ad <- adonis2(d ~ grp, permutations = perms, by = NULL)
  bd <- betadisper(d, group = grp, type = "median", bias.adjust = TRUE)
  dp <- permutest(bd, permutations = perms)

  data.frame(
    Analysis = analysis_name,
    Comparison = comparison_name,
    N_Samples = nrow(mat),
    N_Proteins = ncol(mat),
    Distance = "Bray-Curtis",
    Permutation_Scheme = perm_label,
    PERMANOVA_F = unname(ad$F[1]),
    PERMANOVA_R2 = unname(ad$R2[1]),
    PERMANOVA_p = unname(ad$`Pr(>F)`[1]),
    PERMDISP_F = unname(dp$tab[1, "F"]),
    PERMDISP_p = unname(dp$tab[1, "Pr(>F)"]),
    stringsAsFactors = FALSE
  )
}

run_all_multivariate <- function(M, analysis_name) {
  idx_all <- seq_len(nrow(M))
  idx_adm <- which(meta$Occasion == "Admission")
  idx_rel <- which(meta$Occasion == "Release")
  idx_oc2 <- which(meta$Oiling == "OC2")
  idx_oc4 <- which(meta$Oiling == "OC4")

  rbind(
    run_multivariate_test(M, analysis_name, "All four groups (global overview)", idx_all, meta$Group, blocks = NULL),
    run_multivariate_test(M, analysis_name, "OC2 vs OC4 at admission", idx_adm, meta$Oiling, blocks = NULL),
    run_multivariate_test(M, analysis_name, "OC2 vs OC4 at release", idx_rel, meta$Oiling, blocks = NULL),
    run_multivariate_test(M, analysis_name, "OC2 admission vs release (paired)", idx_oc2, meta$Occasion, blocks = meta$Turtle),
    run_multivariate_test(M, analysis_name, "OC4 admission vs release (paired)", idx_oc4, meta$Occasion, blocks = meta$Turtle)
  )
}

multivar_results <- rbind(
  run_all_multivariate(X[, pca_original$protein_index, drop = FALSE], "Original: raw wSC, all proteins"),
  run_all_multivariate(X_log[, pca_log_all$protein_index, drop = FALSE], "Sensitivity: log2(wSC+1), all proteins"),
  run_all_multivariate(X_log_50[, match(pca_log_50$protein_index, protein_index_50), drop = FALSE], "Sensitivity: log2(wSC+1), >=50% detected in every group")
)

# -------------------------------------------------------------------------
# 6. Save numerical PCA results
# -------------------------------------------------------------------------
all_scores <- do.call(rbind, lapply(pca_results, `[[`, "scores"))
all_loadings <- do.call(rbind, lapply(pca_results, `[[`, "loadings"))
all_summary <- do.call(rbind, lapply(pca_results, `[[`, "summary"))

all_centroids <- aggregate(cbind(PC1, PC2) ~ Analysis + Group, data = all_scores, FUN = mean)

write.csv(all_summary, file.path(out_dir, "PCA_sensitivity_summary.csv"), row.names = FALSE)
write.csv(all_scores, file.path(out_dir, "PCA_sensitivity_scores.csv"), row.names = FALSE)
write.csv(all_centroids, file.path(out_dir, "PCA_sensitivity_group_centroids.csv"), row.names = FALSE)
write.csv(all_loadings, file.path(out_dir, "PCA_sensitivity_loadings.csv"), row.names = FALSE)
write.csv(multivar_results, file.path(out_dir, "PCA_PERMANOVA_PERMDISP_results.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# 7. Sensitivity PCA plots
# -------------------------------------------------------------------------
global_stats_for <- function(label) {
  z <- multivar_results[
    multivar_results$Analysis == label &
      multivar_results$Comparison == "All four groups (global overview)",
    ,
    drop = FALSE
  ]

  if (nrow(z) != 1) return("")

  paste0(
    "Global 4-group overview: PERMANOVA R² = ", sprintf("%.3f", z$PERMANOVA_R2),
    ", p = ", format.pval(z$PERMANOVA_p, digits = 2, eps = 0.0001),
    "; PERMDISP p = ", format.pval(z$PERMDISP_p, digits = 2, eps = 0.0001)
  )
}

plot_pca <- function(result, filename, title_text) {
  d <- result$scores

  p <- ggplot(d, aes(x = PC1, y = PC2, colour = Group, shape = Group)) +
    stat_ellipse(aes(group = Group, colour = Group), type = "t", level = 0.95, linewidth = 0.7, show.legend = FALSE) +
    geom_point(size = 2.8) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.35) +
    labs(
      title = title_text,
      subtitle = global_stats_for(unique(d$Analysis)),
      x = sprintf("PC1 (%.1f%%)", result$variance[1]),
      y = sprintf("PC2 (%.1f%%)", result$variance[2]),
      colour = "Group",
      shape = "Group"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9),
      legend.position = "right"
    )

  ggsave(file.path(out_dir, filename), p, width = 8.0, height = 6.5, units = "in", dpi = 600, bg = "white")
}

plot_pca(pca_original, "PCA_original_raw_wSC_95CI.png", "Original PCA: raw wSC")
plot_pca(pca_log_all, "PCA_log2_all_proteins_95CI.png", "PCA sensitivity: log2(wSC + 1)")
plot_pca(pca_log_50, "PCA_log2_detected50_each_group_95CI.png", "PCA sensitivity: log2(wSC + 1), >=50% detected in every group")

# -------------------------------------------------------------------------
# 8. Publication-quality Figures 4 and 5
# -------------------------------------------------------------------------
publication_colors <- c(
  "OC2 Admission" = "#EF7594",
  "OC2 Release"   = "#7FBDF7",
  "OC4 Admission" = "#FEDB6F",
  "OC4 Release"   = "#41A99C"
)

publication_group_labels <- c(
  "OC2 Admission" = "OC2_0",
  "OC2 Release"   = "OC2_1",
  "OC4 Admission" = "OC4_0",
  "OC4 Release"   = "OC4_1"
)

publication_family <- "Arial"

if (.Platform$OS.type == "windows") {
  try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
}

format_axis_percent <- function(x) {
  if (abs(x - round(x)) < 0.05) sprintf("%.0f%%", x) else sprintf("%.1f%%", x)
}

format_publication_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  paste0("= ", sprintf("%.3f", p))
}

get_publication_stats <- function(comparison_name) {
  z <- multivar_results[
    multivar_results$Analysis == "Original: raw wSC, all proteins" &
      multivar_results$Comparison == comparison_name,
    ,
    drop = FALSE
  ]

  if (nrow(z) != 1) {
    stop("Could not uniquely locate raw-wSC multivariate statistics for: ", comparison_name)
  }
  z
}

make_stats_label <- function(z) {
  paste0(
    "R² = ", sprintf("%.3f", z$PERMANOVA_R2),
    "\n",
    "p ", format_publication_p(z$PERMANOVA_p)
  )
}

run_publication_pca <- function(rows, analysis_label, orientation_reference = NULL, target_pc1_sign = NULL, target_pc2_sign = NULL) {
  M <- X[rows, , drop = FALSE]
  sds <- apply(M, 2, sd)
  keep_var <- is.finite(sds) & sds > 0
  M2 <- M[, keep_var, drop = FALSE]

  ord <- vegan::rda(M2, scale = TRUE)
  eig <- vegan::eigenvals(ord, model = "all")
  variance_percent <- 100 * eig / sum(eig)

  site_scores <- as.data.frame(vegan::scores(ord, display = "sites", choices = 1:2, scaling = 2))
  names(site_scores)[1:2] <- c("PC1", "PC2")

  d <- data.frame(
    Analysis = analysis_label,
    Sample = sample_names[rows],
    Sample_ID = sample_id[rows],
    Turtle = turtle_id[rows],
    Group = droplevels(group4[rows]),
    Oiling = droplevels(oiling[rows]),
    Occasion = droplevels(occasion[rows]),
    PC1 = site_scores$PC1,
    PC2 = site_scores$PC2,
    stringsAsFactors = FALSE
  )

  if (!is.null(orientation_reference)) {
    ref_row <- which(d$Sample_ID == orientation_reference)
    if (length(ref_row) == 1) {
      if (!is.null(target_pc1_sign) && is.finite(d$PC1[ref_row])) {
        if (sign(d$PC1[ref_row]) != sign(target_pc1_sign)) d$PC1 <- -d$PC1
      }
      if (!is.null(target_pc2_sign) && is.finite(d$PC2[ref_row])) {
        if (sign(d$PC2[ref_row]) != sign(target_pc2_sign)) d$PC2 <- -d$PC2
      }
    }
  }

  list(scores = d, variance = variance_percent, ordination = ord, proteins_analyzed = ncol(M2))
}

idx_all <- seq_len(nrow(X))
idx_adm <- which(meta$Occasion == "Admission")
idx_rel <- which(meta$Occasion == "Release")
idx_oc2 <- which(meta$Oiling == "OC2")
idx_oc4 <- which(meta$Oiling == "OC4")

pub_fig4 <- run_publication_pca(idx_all, "Figure 4: all four groups", "116A", -1, 1)
pub_fig5A <- run_publication_pca(idx_adm, "Figure 5A: OC2 vs OC4 at admission", "116A", -1, 1)
pub_fig5B <- run_publication_pca(idx_rel, "Figure 5B: OC2 vs OC4 at release", "249B", -1, 1)
pub_fig5C <- run_publication_pca(idx_oc2, "Figure 5C: OC2 admission vs release", "116A", -1, 1)
pub_fig5D <- run_publication_pca(idx_oc4, "Figure 5D: OC4 admission vs release", "269B", 1, 1)

# -------------------------------------------------------------------------
# Label rule for publication figures
# -------------------------------------------------------------------------
# A sample is labeled ONLY when it falls outside the exact 95% ellipse polygon
# drawn by ggplot2::stat_ellipse() for its own group. Points inside or on the
# visible ellipse are not labeled.
#
# The ellipse polygon used for classification is generated by ggplot2 itself,
# using the same type, level, and number of segments used in the figure.

ellipse_segments <- 200

point_on_polygon_boundary <- function(px, py, vx, vy, tol = 1e-7) {
  n <- length(vx)
  if (n < 2) return(FALSE)

  scale_xy <- max(diff(range(vx)), diff(range(vy)), 1)
  eps <- tol * scale_xy

  for (i in seq_len(n)) {
    j <- if (i == n) 1 else i + 1

    ax <- vx[i]
    ay <- vy[i]
    bx <- vx[j]
    by <- vy[j]

    abx <- bx - ax
    aby <- by - ay
    apx <- px - ax
    apy <- py - ay

    ab2 <- abx^2 + aby^2
    if (ab2 == 0) next

    t <- (apx * abx + apy * aby) / ab2
    t <- max(0, min(1, t))

    cx <- ax + t * abx
    cy <- ay + t * aby

    if (sqrt((px - cx)^2 + (py - cy)^2) <= eps) {
      return(TRUE)
    }
  }

  FALSE
}

point_inside_polygon <- function(px, py, vx, vy) {
  # Treat points on the polygon boundary as inside.
  if (point_on_polygon_boundary(px, py, vx, vy)) {
    return(TRUE)
  }

  n <- length(vx)
  inside <- FALSE
  j <- n

  for (i in seq_len(n)) {
    crosses <- ((vy[i] > py) != (vy[j] > py))

    if (crosses) {
      x_intersect <- (vx[j] - vx[i]) * (py - vy[i]) /
        (vy[j] - vy[i]) + vx[i]

      if (px < x_intersect) {
        inside <- !inside
      }
    }

    j <- i
  }

  inside
}

get_visible_ellipse_outliers <- function(result) {
  d <- result$scores
  d$Group <- droplevels(d$Group)

  outside_ids <- character(0)

  # Build and test each treatment group independently. This makes the
  # classification use the exact same 95% ellipse geometry that is drawn
  # for that group in the publication figure.
  for (g in levels(d$Group)) {
    dg <- d[d$Group == g, , drop = FALSE]

    if (nrow(dg) < 3) next

    ellipse_plot <- ggplot(
      dg,
      aes(x = PC1, y = PC2)
    ) +
      stat_ellipse(
        type = "t",
        level = 0.95,
        geom = "polygon",
        segments = ellipse_segments
      )

    built <- ggplot_build(ellipse_plot)
    poly <- built$data[[1]]

    poly <- poly[
      is.finite(poly$x) & is.finite(poly$y),
      ,
      drop = FALSE
    ]

    if (nrow(poly) < 3) next

    for (i in seq_len(nrow(dg))) {
      inside <- point_inside_polygon(
        dg$PC1[i],
        dg$PC2[i],
        poly$x,
        poly$y
      )

      if (!inside) {
        outside_ids <- c(outside_ids, dg$Sample_ID[i])
      }
    }
  }

  unique(outside_ids)
}

fig4_label_ids  <- get_visible_ellipse_outliers(pub_fig4)
fig5A_label_ids <- get_visible_ellipse_outliers(pub_fig5A)
fig5B_label_ids <- get_visible_ellipse_outliers(pub_fig5B)
fig5C_label_ids <- get_visible_ellipse_outliers(pub_fig5C)
fig5D_label_ids <- get_visible_ellipse_outliers(pub_fig5D)

cat("\nSamples labeled because they fall outside their own visible 95% ellipse:\n")
cat("Figure 4 :", paste(fig4_label_ids, collapse = ", "), "\n")
cat("Figure 5A:", paste(fig5A_label_ids, collapse = ", "), "\n")
cat("Figure 5B:", paste(fig5B_label_ids, collapse = ", "), "\n")
cat("Figure 5C:", paste(fig5C_label_ids, collapse = ", "), "\n")
cat("Figure 5D:", paste(fig5D_label_ids, collapse = ", "), "\n\n")

# Save an audit table so the exact samples labeled in each panel are recorded.
ellipse_label_audit <- rbind(
  data.frame(Panel = "Figure 4",  Sample_ID = fig4_label_ids),
  data.frame(Panel = "Figure 5A", Sample_ID = fig5A_label_ids),
  data.frame(Panel = "Figure 5B", Sample_ID = fig5B_label_ids),
  data.frame(Panel = "Figure 5C", Sample_ID = fig5C_label_ids),
  data.frame(Panel = "Figure 5D", Sample_ID = fig5D_label_ids)
)

write.csv(
  ellipse_label_audit,
  file.path(out_dir, "PCA_ellipse_outside_labels.csv"),
  row.names = FALSE
)

publication_theme <- function(base_size = 11) {
  theme_classic(base_family = publication_family, base_size = base_size) +
    theme(
      axis.title = element_text(family = publication_family, face = "bold", colour = "black", size = base_size + 4),
      axis.text = element_text(family = publication_family, colour = "black", size = base_size + 2),
      axis.line = element_line(colour = "black", linewidth = 0.70),
      axis.ticks = element_line(colour = "black", linewidth = 0.60),
      legend.title = element_blank(),
      legend.text = element_text(family = publication_family, face = "bold", colour = "black", size = base_size + 1),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.background = element_blank(),
      legend.position = "right",
      plot.tag = element_text(family = publication_family, face = "bold", colour = "black", size = base_size + 10),
      plot.tag.position = c(0.02, 0.98),
      plot.margin = margin(8, 10, 8, 8)
    )
}

make_publication_pca_plot <- function(
    result,
    comparison_name,
    label_ids = character(0),
    panel_tag = NULL,
    point_size = 2.8,
    label_size = 3.2,
    base_size = 11,
    repel_force = 0.70,
    repel_box_padding = 0.18,
    repel_point_padding = 0.10,
    repel_segment_color = "grey70",
    repel_max_time = 2.0) {

  d <- result$scores
  stat_row <- get_publication_stats(comparison_name)
  stat_text <- make_stats_label(stat_row)
  label_data <- d[d$Sample_ID %in% label_ids, , drop = FALSE]

  group_breaks <- intersect(names(publication_colors), levels(droplevels(d$Group)))

  p <- ggplot(d, aes(x = PC1, y = PC2, colour = Group, fill = Group)) +
    stat_ellipse(
      aes(group = Group),
      type = "t",
      level = 0.95,
      geom = "polygon",
      segments = ellipse_segments,
      alpha = 0.22,
      linewidth = 0.85,
      show.legend = FALSE
    ) +
    geom_point(shape = 21, size = point_size, stroke = 0.35, alpha = 0.82) +
    scale_colour_manual(values = publication_colors, breaks = group_breaks, labels = publication_group_labels[group_breaks], drop = TRUE, name = NULL, guide = "none") +
    scale_fill_manual(
      values = publication_colors,
      breaks = group_breaks,
      labels = publication_group_labels[group_breaks],
      drop = TRUE,
      name = NULL,
      guide = guide_legend(override.aes = list(shape = 21, size = 5, alpha = 0.65))
    ) +
    labs(
      x = paste0("PC1 (", format_axis_percent(result$variance[1]), ")"),
      y = paste0("PC2 (", format_axis_percent(result$variance[2]), ")"),
      tag = panel_tag
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = stat_text,
      family = publication_family,
      fontface = "bold",
      colour = "black",
      size = if (is.null(panel_tag)) 5.7 else 4.6,
      lineheight = 1.45,
      hjust = 1.08,
      vjust = 1.18
    ) +
    publication_theme(base_size = base_size) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.15))) +
    scale_y_continuous(expand = expansion(mult = c(0.07, 0.12))) +
    coord_cartesian(clip = "off")

  if (nrow(label_data) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = label_data,
        aes(x = PC1, y = PC2, label = Sample_ID),
        inherit.aes = FALSE,
        family = publication_family,
        colour = "black",
        size = label_size,
        seed = 20260824,
        box.padding = repel_box_padding,
        point.padding = repel_point_padding,
        force = repel_force,
        max.overlaps = Inf,
        max.time = repel_max_time,
        segment.color = repel_segment_color,
        segment.size = 0.25,
        min.segment.length = 0.08,
        show.legend = FALSE
      )
  }

  p
}

figure4_publication <- make_publication_pca_plot(
  pub_fig4,
  "All four groups (global overview)",
  label_ids = fig4_label_ids,
  panel_tag = NULL,
  point_size = 3.3,
  label_size = 4.8,
  base_size = 12.5,
  repel_force = 0.90,
  repel_box_padding = 0.22,
  repel_point_padding = 0.12,
  repel_max_time = 2.5
)

figure5A_publication <- make_publication_pca_plot(
  pub_fig5A,
  "OC2 vs OC4 at admission",
  label_ids = fig5A_label_ids,
  panel_tag = "A",
  point_size = 2.7,
  label_size = 3.7,
  base_size = 11,
  repel_force = 0.85,
  repel_max_time = 2.0
)

figure5B_publication <- make_publication_pca_plot(
  pub_fig5B,
  "OC2 vs OC4 at release",
  label_ids = fig5B_label_ids,
  panel_tag = "B",
  point_size = 2.7,
  label_size = 3.6,
  base_size = 11,
  repel_force = 0.95,
  repel_box_padding = 0.20,
  repel_point_padding = 0.12,
  repel_max_time = 2.5
)

figure5C_publication <- make_publication_pca_plot(
  pub_fig5C,
  "OC2 admission vs release (paired)",
  label_ids = fig5C_label_ids,
  panel_tag = "C",
  point_size = 2.6,
  label_size = 3.55,
  base_size = 11,
  repel_force = 1.00,
  repel_box_padding = 0.22,
  repel_point_padding = 0.13,
  repel_max_time = 2.8
)

figure5D_publication <- make_publication_pca_plot(
  pub_fig5D,
  "OC4 admission vs release (paired)",
  label_ids = fig5D_label_ids,
  panel_tag = "D",
  point_size = 2.6,
  label_size = 3.50,
  base_size = 11,
  repel_force = 1.05,
  repel_box_padding = 0.22,
  repel_point_padding = 0.13,
  repel_max_time = 3.0
)

figure5_publication <-
  (figure5A_publication + figure5B_publication) /
  (figure5C_publication + figure5D_publication) +
  patchwork::plot_layout(guides = "keep")

ggsave(file.path(out_dir, "Figure4_PCA_publication_600dpi.png"), figure4_publication, width = 7.4, height = 6.0, units = "in", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure4_PCA_publication.pdf"), figure4_publication, width = 7.4, height = 6.0, units = "in", device = grDevices::cairo_pdf, bg = "white")

ggsave(file.path(out_dir, "Figure5_PCA_publication_600dpi.png"), figure5_publication, width = 11.0, height = 8.5, units = "in", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure5_PCA_publication.pdf"), figure5_publication, width = 11.0, height = 8.5, units = "in", device = grDevices::cairo_pdf, bg = "white")

for (nm in c("A", "B", "C", "D")) {
  panel_plot <- get(paste0("figure5", nm, "_publication"))
  ggsave(
    file.path(out_dir, paste0("Figure5", nm, "_PCA_publication_600dpi.png")),
    panel_plot,
    width = 5.5,
    height = 4.25,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

publication_scores <- rbind(
  transform(pub_fig4$scores, Panel = "Figure 4"),
  transform(pub_fig5A$scores, Panel = "Figure 5A"),
  transform(pub_fig5B$scores, Panel = "Figure 5B"),
  transform(pub_fig5C$scores, Panel = "Figure 5C"),
  transform(pub_fig5D$scores, Panel = "Figure 5D")
)

write.csv(publication_scores, file.path(out_dir, "Publication_Figures4_5_PCA_scores.csv"), row.names = FALSE)

publication_variance <- data.frame(
  Figure_Panel = c("Figure 4", "Figure 5A", "Figure 5B", "Figure 5C", "Figure 5D"),
  Proteins_Analyzed = c(pub_fig4$proteins_analyzed, pub_fig5A$proteins_analyzed, pub_fig5B$proteins_analyzed, pub_fig5C$proteins_analyzed, pub_fig5D$proteins_analyzed),
  PC1_Percent = c(pub_fig4$variance[1], pub_fig5A$variance[1], pub_fig5B$variance[1], pub_fig5C$variance[1], pub_fig5D$variance[1]),
  PC2_Percent = c(pub_fig4$variance[2], pub_fig5A$variance[2], pub_fig5B$variance[2], pub_fig5C$variance[2], pub_fig5D$variance[2]),
  stringsAsFactors = FALSE
)

write.csv(publication_variance, file.path(out_dir, "Publication_Figures4_5_PCA_variance.csv"), row.names = FALSE)

publication_readme <- c(
  "Kemp's ridley plasma proteomics - publication Figures 4 and 5",
  "",
  "Source data:",
  "- Supplemental_Tables_1_through_8.xlsx, Table S3",
  "- Original Scaffold weighted spectral counts (wSC)",
  "- Numerical zeros retained as Scaffold-exported non-detections",
  "- No missing values imputed",
  "",
  "PCA:",
  "- vegan::rda(scale=TRUE)",
  "- site scores displayed with scaling=2",
  "- zero-variance proteins removed separately within each PCA",
  "- 95% t-distribution ellipses",
  "- sample labels are shown only when a point falls outside the exact 95% ellipse polygon drawn for its own group",
  "- PCA axis signs may be flipped only to match the orientation of the original manuscript figures; sign flips do not alter results",
  "",
  "Figure 4:",
  "- all 48 samples",
  "- all four oiling-category/sampling-occasion groups",
  "- global four-group PERMANOVA is descriptive because turtles are sampled repeatedly",
  "",
  "Figure 5:",
  "- A: OC2 vs OC4 at admission",
  "- B: OC2 vs OC4 at release",
  "- C: OC2 admission vs release; permutations blocked within turtle",
  "- D: OC4 admission vs release; permutations blocked within turtle",
  "",
  "Statistics:",
  "- PERMANOVA/PERMDISP use Bray-Curtis dissimilarity",
  "- 9,999 permutations",
  "- PERMDISP uses spatial median with bias.adjust=TRUE",
  "- statistics displayed on Figures 4 and 5 are read directly from the raw-wSC multivariate results generated earlier in this script",
  "",
  "Files:",
  "- Figure4_PCA_publication_600dpi.png",
  "- Figure4_PCA_publication.pdf",
  "- Figure5_PCA_publication_600dpi.png",
  "- Figure5_PCA_publication.pdf",
  "- Figure5A-D_PCA_publication_600dpi.png",
  "- Publication_Figures4_5_PCA_scores.csv",
  "- Publication_Figures4_5_PCA_variance.csv"
)

writeLines(publication_readme, file.path(out_dir, "Publication_Figures4_5_README.txt"))

cat("\nPublication Figures 4 and 5 were generated.\n")
cat("Figure 4: Figure4_PCA_publication_600dpi.png / .pdf\n")
cat("Figure 5: Figure5_PCA_publication_600dpi.png / .pdf\n\n")

# -------------------------------------------------------------------------
# 9. README for sensitivity outputs
# -------------------------------------------------------------------------
readme <- c(
  "Kemp's ridley plasma proteomics - PCA / PERMANOVA / PERMDISP sensitivity analysis",
  "",
  paste("Input workbook:", input_file),
  "",
  "PCA analyses:",
  "1. Original raw Scaffold weighted spectral counts, all proteins.",
  "2. log2(wSC + 1), all proteins.",
  "3. log2(wSC + 1), proteins detected in >=50% of samples within every group.",
  "",
  "PCA plotting:",
  "- PCA is centered and unit-variance scaled.",
  "- Ellipses are drawn at the 95% level using ggplot2::stat_ellipse(type='t').",
  "",
  "Multivariate testing:",
  "- PERMANOVA and PERMDISP use Bray-Curtis dissimilarity.",
  "- PERMDISP uses the spatial median with bias.adjust=TRUE.",
  "- OC2 and OC4 admission-to-release permutations are blocked within turtle.",
  "- The four-group test is a descriptive/global overview because the dataset contains repeated samples.",
  "- Comparison-specific tests should be emphasized for inference.",
  "",
  "Zero values:",
  "- Numerical zeros are values exported by Scaffold for non-detections.",
  "- No missing values are imputed.",
  "",
  "Output files:",
  "- PCA_sensitivity_summary.csv",
  "- PCA_sensitivity_scores.csv",
  "- PCA_sensitivity_group_centroids.csv",
  "- PCA_sensitivity_loadings.csv",
  "- PCA_PERMANOVA_PERMDISP_results.csv",
  "- PCA_original_raw_wSC_95CI.png",
  "- PCA_log2_all_proteins_95CI.png",
  "- PCA_log2_detected50_each_group_95CI.png"
)

writeLines(readme, file.path(out_dir, "PCA_PERMANOVA_PERMDISP_README.txt"))

# -------------------------------------------------------------------------
# 10. Console report
# -------------------------------------------------------------------------
cat("\n====================================================\n")
cat("PCA / PERMANOVA / PERMDISP ANALYSIS COMPLETE\n")
cat("====================================================\n\n")

cat("PCA variance explained:\n")
print(all_summary, row.names = FALSE)

cat("\nPERMANOVA / PERMDISP results:\n")
print(multivar_results, row.names = FALSE)

cat(
  "\nZero wSC cells in original matrix: ",
  sum(X == 0), " / ", length(X),
  " (", sprintf("%.2f", 100 * mean(X == 0)), "%)\n",
  sep = ""
)

cat(
  "\nProteins passing >=50% detection in every group: ",
  sum(detected50_every_group), " / ", ncol(X), "\n",
  sep = ""
)

cat("\nAll result files were saved to:\n")
cat(out_dir, "\n\n")
